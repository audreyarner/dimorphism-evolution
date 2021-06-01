#!/usr/bin/env Rscript

#===================================================================================
#--Analysis to determine if SNPs are in genic or intergenic regions
#===================================================================================

library(dplyr)

#Genes from the GO database
AllGOGenes <- read.table("AllGenesGRCh37.txt", header = TRUE)

#function to run a permutation on SexDiff and phenotype-associated SNPs in genic/intergenic regions
getGenicPerm = function(PdiffPruned, PhenotypePruned){
  Pdiff<-read.table(PdiffPruned, header =T)
  All<-read.table(PhenotypePruned, header =T)
  PhenotypeName <- gsub("PdiffPrunedSNPs.txt", "", PdiffPruned)
  
  #Identify if Phenotype SNPs are in intergenic or genic regions
  All[,"GeneGroup"] <- NA
  
  for(i in 1:nrow(All)){ 
    Gene.Overlap = which((All$CHR[i] == AllGOGenes$CHR) & (AllGOGenes$TenThousandStart <= All$POS[i]) &(AllGOGenes$TenThousandEnd>=All$POS[i]))
    if(length(Gene.Overlap >=1)){
      All$GeneGroup[i]="G"
    } else {All$GeneGroup[i]="IG"}
  }
  
  #identify if SexDiff-associated SNPs are in intergenic or genic regions
  Pdiff[,"GeneGroup"] <- NA
  
  for(i in 1:nrow(Pdiff)){ 
    Gene.Overlap = which((Pdiff$CHR[i] == AllGOGenes$CHR) & (AllGOGenes$TenThousandStart <= Pdiff$POS[i]) &(AllGOGenes$TenThousandEnd>=Pdiff$POS[i]))
    if(length(Gene.Overlap >=1)){
      Pdiff$GeneGroup[i]="G"
    } else {Pdiff$GeneGroup[i]="IG"}
  }
  
  #SNPs disproportionately impacting female phenotype
  PrunedPdiffFemale <- Pdiff %>% filter(Female.pval <= Male.pval)
  
  #Find number of SNPs in intergenic regions
  PrunedPdiffFemale2 <- PrunedPdiffFemale %>% filter(GeneGroup=="IG")
  print(nrow(PrunedPdiffFemale2))
  
  NumberFemale <- nrow(PrunedPdiffFemale)
  
  FemaleGenic <- PrunedPdiffFemale %>% filter(GeneGroup == "G")
  NumberFemaleGenic <- nrow(FemaleGenic)
  
  #permutation to identify if there are more SNPs in genic regions than expected by chance
  t1 <- rep(NA, 10000)
  
  for(j in 1:10000){
    table1 <- All[sample(nrow(All), NumberFemale),]
    table2 <- table1 %>% filter(GeneGroup == "G")
    t1[j] <- nrow(table2)
  }
  
  write.table(t1, file=paste(PhenotypeName, "FemaleGenicPerm.txt", sep=""), sep="/t", row.names=F)
  
  df <- as.data.frame(t1)
  names(df)[1] <- "Genic"
  
  #Calculate percent of permutation
  CalcPercent <- df %>% filter(Genic>=NumberFemaleGenic)
  HighPercent <-nrow(CalcPercent)/10000
  LowPercent <-1-HighPercent
  
  print(paste("Female genic percentile is", HighPercent))
  print(paste("Female genic lower percentile is", LowPercent))
  
  #SNPs disproportionately impacting male phenotype
  PrunedPdiffMale <- Pdiff %>% filter(Male.pval <= Female.pval)
 
  #Find number of SNPs in intergenic regions 
  PrunedPdiffMale2 <- PrunedPdiffMale %>% filter(GeneGroup=="IG")
  print(nrow(PrunedPdiffMale2))
  
  NumberMale <- nrow(PrunedPdiffMale)
  MaleGenic <- PrunedPdiffMale %>% filter(GeneGroup == "G")
  NumberMaleGenic <- nrow(MaleGenic)
  
  t1 <- rep(NA, 10000)
  
  for(j in 1:10000){
    table1 <- All[sample(nrow(All), NumberMale),]
    table2 <- table1 %>% filter(GeneGroup == "G")
    t1[j] <- nrow(table2)
  }
  
  write.table(t1, file=paste(PhenotypeName, "MaleGenicPerm.txt", sep=""), sep="/t", row.names=F)
  
  df <- as.data.frame(t1)
  names(df)[1] <- "Genic"
  
  CalcPercent <- df %>% filter(Genic>=NumberMaleGenic)
  HighPercent <-nrow(CalcPercent)/10000
  LowPercent <-1-HighPercent
  
  print(paste(PhenotypeName, "Male genic percentile is", HighPercent))
  print(paste("Male genic lower percentile is", LowPercent))
}

getGenicPerm("WaistCircPdiffPrunedSNPs.txt", "WaistCircPrunedSNPs.txt")
getGenicPerm("BodyMassPdiffPrunedSNPs.txt", "BodyMassPrunedSNPs.txt")
getGenicPerm("BodyFatPercentPdiffPrunedSNPs.txt", "BodyFatPercentPrunedSNPs.txt")
getGenicPerm("HeightPdiffPrunedSNPs.txt", "HeightPrunedSNPs.txt")
getGenicPerm("HipCircPdiffPrunedSNPs.txt", "HipCircPrunedSNPs.txt")

#list of permutation p-values
list <- c(.0584,.2489,.0735, .0024,.1556,.1093,.0749,.0499,.3102,.0089)

#multiply p-values by two to account for two-sided test
list2 <- list *2

#FDR adjustment of p-values to correct for multiple testing
p.adjust(list2, method="fdr")
