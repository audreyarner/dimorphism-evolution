#!/usr/bin/env Rscript

#===================================================================================
#--trait-SDS (tSDS) analysis for SexDiff-SNPs compared to phenotype-associated SNPs
#===================================================================================

library(dplyr)

#reads in the table with LD blocks
LDBlocks <- read.table("LDBlocks.txt", header = T)

#function to prune SNPs in the Pdiff4 category for each phenotype
getPrunedPdiff = function(SNPs){
  Phenotype <- read.table(SNPs, header = T)
  PhenotypeName <- gsub("Pdiff4.tsv", "", SNPs)
    
  Phenotype$CHR <- as.numeric(Phenotype$CHR)
  LDBlocks$CHR <- as.numeric(LDBlocks$CHR)

  Phenotype$POS <- as.numeric(Phenotype$POS)
  LDBlocks$Block_Start <- as.numeric(LDBlocks$Block_Start)
  LDBlocks$Block_End <- as.numeric(LDBlocks$Block_End)

  PrunedSNPsOutput <- rep('a', nrow(LDBlocks))

  for(i in 1:nrow(LDBlocks)){
    BlockSNPs <- Phenotype %>% filter((LDBlocks$CHR[i] == Phenotype$CHR) & (LDBlocks$Block_Start[i] <= Phenotype$POS) & (LDBlocks$Block_End[i] >= Phenotype$POS))
    OrderbyP <- BlockSNPs[order(BlockSNPs$pdiff.FDR),]
    OrderbyP$rsid <- as.character(OrderbyP$rsid)
    PrunedSNPsOutput[i] <- OrderbyP$rsid[1]
  }

  PrunedSigSNPs <- Phenotype[Phenotype$rsid %in% PrunedSNPsOutput,]
  
  PrunedSigSNPs[,"Index"] <- NA
  
  for (i in 1:nrow(PrunedSigSNPs)){
    #which() will give the position of the gene that overlaps, which we can then use to find the gene name
    SNP.Overlap = which((LDBlocks$CHR == PrunedSigSNPs$CHR[i]) & (LDBlocks$Block_Start <= PrunedSigSNPs$POS[i]) &(LDBlocks$Block_End>=PrunedSigSNPs$POS[i]))
    if(length(SNP.Overlap>=1)){
      PrunedSigSNPs$Index[i] = LDBlocks$Index[SNP.Overlap]
    } else {SNPs$GeneOverlap[i] = 0}
  } 
  
  write.table(PrunedSigSNPs, file=paste(PhenotypeName, "PdiffPrunedSNPs.txt", sep=""), sep="\t", row.names=F)
}

getPrunedPdiff("HeightPdiff4.tsv")
getPrunedPdiff("HipCircPdiff4.tsv")
getPrunedPdiff("BodyFatPercentPdiff4.tsv")
getPrunedPdiff("WaistCircPdiff4.tsv")
getPrunedPdiff("BodyMassPdiff4.tsv")

gettSDSPerm = function(Pdiff, SigPheno){
  PrunedPdiff <- read.table(Pdiff, header = T)
  PhenotypeName <- gsub("PdiffPrunedSNPs.txt", "", Pdiff)
  
  #get average trait-SDS for SNPs disproportionately impacting Female height
  PrunedPdiffFemale <- PrunedPdiff %>% filter(Female.pval <= Male.pval)
  AveragetSDSFemale <- mean(PrunedPdiffFemale$Female.tSDS)
  print(paste("Female Average tSDS is", AveragetSDSFemale))
  TotalNumberSNPsF <- nrow(PrunedPdiffFemale)
  print(paste("Female number SNPs is", TotalNumberSNPsF))
  
  #get average tSDS for SNPs disproportionately impacting male height
  PrunedPdiffMale <- PrunedPdiff %>% filter(Male.pval <= Female.pval)
  AveragetSDSMale <- mean(PrunedPdiffMale$Male.tSDS)
  print(paste("Male Average tSDS is", AveragetSDSMale))
  TotalNumberSNPsM <- nrow(PrunedPdiffMale)
  print(paste("Male number SNPs is", TotalNumberSNPsM))
  
  #average tSDS permutation Female
  PhenoSNPs <- read.table(SigPheno, header = T)
  
  PhenoSNPs$traitSDS <- NA
  
  for (i in 1:nrow(PhenoSNPs)){
    #which() will give the position of the gene that overlaps, which we can then use to find the gene name
    if(PhenoSNPs$Female.pval[i]<PhenoSNPs$Male.pval[i]){
      PhenoSNPs$traitSDS[i] = PhenoSNPs$Female.tSDS[i]
    } else {PhenoSNPs$traitSDS[i] = PhenoSNPs$Male.tSDS[i]}
  }
  
  t1 <- rep(NA, 10000)
  
  for(j in 1:10000){
    table1 <- PhenoSNPs[sample(nrow(PhenoSNPs), TotalNumberSNPsF),]
    t1[j] <- mean(table1$traitSDS)
  }
  
  write.table(t1, file=paste(PhenotypeName, "FemaletSDSPerm.txt", sep=""), sep="/t", row.names=F)
  
  df <- as.data.frame(t1)
  names(df)[1] <- "tSDS"
  
  CalcPercent <- df %>% filter(tSDS>=AveragetSDSFemale)
  HighPercent <-nrow(CalcPercent)/10000
  
  print(paste("Female high Percentile is", HighPercent))
  
  CalcPercent <- df %>% filter(tSDS<=AveragetSDSFemale)
  LowPercent <-nrow(CalcPercent)/10000
  
  print(paste("Female low Percentile is", LowPercent, sep = ' '))
  
  #average tSDS permutation male
  
  t1 <- rep(NA, 10000)
  
  for(j in 1:10000)
  {
    table1 <- PhenoSNPs[sample(nrow(PhenoSNPs), TotalNumberSNPsM),]
    t1[j] <- mean(table1$traitSDS)
  }
  
  write.table(t1, file=paste(PhenotypeName, "MaletSDSPerm.txt", sep=""), sep="\t", row.names=F)
  
  df <- as.data.frame(t1)
  names(df)[1] <- "tSDS"
  
  CalcPercent <- df %>% filter(tSDS>=AveragetSDSMale)
  HighPercent <-nrow(CalcPercent)/10000
  
  print(paste("Male high Percentile is", HighPercent, sep = ' '))
  
  CalcPercent <- df %>% filter(tSDS<=AveragetSDSMale)
  LowPercent <-nrow(CalcPercent)/10000
  
  print(paste("Male low Percentile is", LowPercent, sep = ' '))
  
}

gettSDSPerm("HeightPdiffPrunedSNPs.txt", "HeightPrunedSNPs.txt")
gettSDSPerm("HipCircPdiffPrunedSNPs.txt", "HipCircPrunedSNPs.txt")
gettSDSPerm("WaistCircPdiffPrunedSNPs.txt", "WaistCircPrunedSNPs.txt")
gettSDSPerm("BodyFatPercentPdiffPrunedSNPs.txt", "BodyFatPercentPrunedSNPs.txt")
gettSDSPerm("BodyMassPdiffPrunedSNPs.txt", "BodyMassPrunedSNPs.txt")


