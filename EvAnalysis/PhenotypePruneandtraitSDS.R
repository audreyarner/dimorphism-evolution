#!/usr/bin/env Rscript

#------------------------------------------------------------------------------------------------------
#this script is used for our trait-SDS analysis of phenotype wide scores
#------------------------------------------------------------------------------------------------------

#Sexually dimorphic SNPs are pruned

library(dplyr)

#linkage disequilibrium blocks
LDBlocks <- read.table("NEWall_fourier_ls.bed")

#function to prune SNPs in the Pdiff4 category for each phenotype
getPrunedPdiff = function(SNPs){
  Phenotype <- read.table(SNPs, header = T)
  PhenotypeName <- gsub("withPDIFF.FDR.tsv", "", SNPs)
  Phenotype$CHR <- as.numeric(Phenotype$CHR)
  LDBlocks$V1 <- as.numeric(LDBlocks$V1)
  
  Phenotype$POS <- as.numeric(Phenotype$POS)
  LDBlocks$V2 <- as.numeric(LDBlocks$V2)
  LDBlocks$V3 <- as.numeric(LDBlocks$V3)
  
  Phenotype[,"pval"] <- NA

  Phenotype <- Phenotype %>% filter(WOMEN.pval < 5e-8)
  
  for (i in 1:nrow(Phenotype)){
    if(Phenotype$MEN.pval[i]<Phenotype$WOMEN.pval[i]){
      Phenotype$pval[i]=.9
    } else {Phenotype$pval[i]=Phenotype$WOMEN.pval[i]}
  }
  
  PrunedSNPsOutput <- rep('a', nrow(LDBlocks))
  
  for(i in 1:nrow(LDBlocks)){
    BlockSNPs <- Phenotype %>% filter((LDBlocks$V1[i] == Phenotype$CHR) & (LDBlocks$V2[i] <= Phenotype$POS) & (LDBlocks$V3[i] >= Phenotype$POS))
    OrderbyP <- BlockSNPs[order(BlockSNPs$pval),]
    OrderbyP$rsid <- as.character(OrderbyP$rsid)
    PrunedSNPsOutput[i] <- OrderbyP$rsid[1]
  }
  
  PrunedSigSNPs <- Phenotype[Phenotype$rsid %in% PrunedSNPsOutput,]
  write.table(PrunedSigSNPs, file=paste(PhenotypeName, "PrunedSNPs.tsv", sep=""), row.names=F)
}

getPrunedPdiff("HeightwithPDIFF.FDR.tsv")
getPrunedPdiff("BodyMasswithPDIFF.FDR.tsv")
getPrunedPdiff("HipCircwithPDIFF.FDR.tsv")
getPrunedPdiff("BodyFatPercentwithPDIFF.FDR.tsv")
getPrunedPdiff("WaistCircwithPDIFF.FDR.tsv")

gentraitSDS = function(PrunedPheno){
  PrunedSNPs <- read.table(PrunedPheno, header = T)
  PhenotypeName <- gsub("PrunedSNPs.tsv", "", PrunedPheno)

  #Permutation of Female p-diff t-SDS score to genome-wide tSDS scores
  
  PrunedSNPs[,"tSDS"] <- NA
  
  for (i in 1:nrow(PrunedSNPs)){
    if(PrunedSNPs$MEN.pval[i]<PrunedSNPs$WOMEN.pval[i]){
      PrunedSNPs$tSDS[i]=PrunedSNPs$MEN.tSDS[i]
    } else {PrunedSNPs$tSDS[i]=PrunedSNPs$WOMEN.tSDS[i]}
  }
  
  AveragetSDS <- mean(PrunedSNPs$tSDS)
  
  print(paste(PhenotypeName, "average tSDS is", AveragetSDS))
}

gentraitSDS("HeightPrunedSNPs.tsv")
gentraitSDS("BodyMassPrunedSNPs.tsv")
gentraitSDS("HipCircPrunedSNPs.tsv")
gentraitSDS("BodyFatPercentPrunedSNPs.tsv")
gentraitSDS("WaistCircPrunedSNPs.tsv")

