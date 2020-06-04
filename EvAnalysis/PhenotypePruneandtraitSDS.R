#!/usr/bin/env Rscript

#------------------------------------------------------------------------------------------------------
#this script is used for our trait-SDS analysis of phenotype wide scores
#------------------------------------------------------------------------------------------------------

#Sexually dimorphic SNPs are pruned

library(dplyr)

#linkage disequilibrium blocks
LDBlocks <- read.table("LDBlocks.txt", header =T)

#function to prune SNPs in the Pdiff4 category for each phenotype
getPrunedPdiff = function(SNPs){
  Phenotype <- read.table(SNPs, header = T)
  PhenotypeName <- gsub("withPDIFF.FDR.tsv", "", SNPs)
  
  Phenotype$CHR <- as.numeric(Phenotype$CHR)
  LDBlocks$CHR <- as.numeric(LDBlocks$CHR)
  
  Phenotype$POS <- as.numeric(Phenotype$POS)
  LDBlocks$Block_Start <- as.numeric(LDBlocks$Block_Start)
  LDBlocks$Block_End <- as.numeric(LDBlocks$Block_End)
  
  Phenotype[,"pval"] <- NA
  
  for (i in 1:nrow(Phenotype)){
    if(Phenotype$Male.pval[i]<Phenotype$Female.pval[i]){
      Phenotype$pval[i]=Phenotype$Male.pval[i]
    } else {Phenotype$pval[i]=Phenotype$Female.pval[i]}
  }

  PrunedSNPsOutput <- rep('a', nrow(LDBlocks))
  
  for(i in 1:nrow(LDBlocks)){
    BlockSNPs <- Phenotype %>% filter((LDBlocks$CHR[i] == Phenotype$CHR) & (LDBlocks$Block_Start[i] <= Phenotype$POS) & (LDBlocks$Block_End[i] >= Phenotype$POS))
    OrderbyP <- BlockSNPs[order(BlockSNPs$pval),]
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
  
  write.table(PrunedSigSNPs, file=paste(PhenotypeName, "PrunedSNPs.txt", sep=""), sep="\t", row.names=F, quote=F)
}

getPrunedPdiff("HeightwithPDIFF.FDR.tsv")
getPrunedPdiff("BodyMasswithPDIFF.FDR.tsv")
getPrunedPdiff("HipCircwithPDIFF.FDR.tsv")
getPrunedPdiff("BodyFatPercentwithPDIFF.FDR.tsv")
getPrunedPdiff("WaistCircwithPDIFF.FDR.tsv")

gentraitSDS = function(PrunedPheno){
  PrunedSNPs <- read.table(PrunedPheno, header = T)
  PhenotypeName <- gsub("PrunedSNPs.txt", "", PrunedPheno)

  #Permutation of Female p-diff t-SDS score to genome-wide tSDS scores
  
  PrunedSNPs[,"tSDS"] <- NA
  
  for (i in 1:nrow(PrunedSNPs)){
    if(PrunedSNPs$Male.pval[i]<PrunedSNPs$Female.pval[i]){
      PrunedSNPs$tSDS[i]=PrunedSNPs$Male.tSDS[i]
    } else {PrunedSNPs$tSDS[i]=PrunedSNPs$Female.tSDS[i]}
  }
  
  AveragetSDS <- mean(PrunedSNPs$tSDS)
  
  print(paste(PhenotypeName, "average tSDS is", AveragetSDS))
}

gentraitSDS("HeightPrunedSNPs.txt")
gentraitSDS("BodyMassPrunedSNPs.txt")
gentraitSDS("HipCircPrunedSNPs.txt")
gentraitSDS("BodyFatPercentPrunedSNPs.txt")
gentraitSDS("WaistCircPrunedSNPs.txt")
