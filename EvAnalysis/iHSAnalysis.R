#!/usr/bin/env Rscript

#===================================================================================
#--trait-SDS (iHS) analysis for SexDiff-SNPs compared to phenotype-associated SNPs
#===================================================================================

library(dplyr)

getiHSPerm = function(Pdiff, SigPheno){
  PrunedPdiff <- read.table(Pdiff, header = T)
  PhenotypeName <- gsub("iHSScores.txt", "", Pdiff)
  
  PrunedPdiff <- mutate(PrunedPdiff, AbsStiHS=abs(PrunedPdiff$StiHS))
  
  #get average trait-SDS for SNPs disproportionately impacting Female height
  PrunedPdiffFemale1 <- PrunedPdiff %>% filter(Female.pval <= Male.pval)
  PrunedPdiffFemale <- na.omit(PrunedPdiffFemale1)
  AverageiHSFemale <- mean(PrunedPdiffFemale$AbsStiHS)
  print(paste("Female Average iHS is", AverageiHSFemale))
  TotalNumberSNPsF <- nrow(PrunedPdiffFemale)
  print(paste("Female number SNPs is", TotalNumberSNPsF))
  
  #get average iHS for SNPs disproportionately impacting male height
  PrunedPdiffMale1 <- PrunedPdiff %>% filter(Male.pval <= Female.pval)
  PrunedPdiffMale <- na.omit(PrunedPdiffMale1)
  AverageiHSMale <- mean(PrunedPdiffMale$AbsStiHS)
  print(paste("Male Average iHS is", AverageiHSMale))
  TotalNumberSNPsM <- nrow(PrunedPdiffMale)
  print(paste("Male number SNPs is", TotalNumberSNPsM))
  
  #average iHS permutation Female
  PhenoSNPs <- read.table(SigPheno, header = T)
  
  PhenoSNPs <- mutate(PhenoSNPs, AbsStiHS=abs(PhenoSNPs$StiHS))
  PhenoSNPs <- na.omit(PhenoSNPs)
  
  AverageiHSPheno <- mean(PhenoSNPs$AbsStiHS)
  print(paste("Pheno Average iHS is", AverageiHSPheno))
  TotalNumberSNPs <- nrow(PhenoSNPs)
  print(paste("Pheno number SNPs is", TotalNumberSNPs))

  t1 <- rep(NA, 10000)
  
  for(j in 1:10000){
    table1 <- PhenoSNPs[sample(nrow(PhenoSNPs), TotalNumberSNPsF),]
    t1[j] <- mean(table1$AbsStiHS)
  }
  
  write.table(t1, file=paste(PhenotypeName, "FemaleiHSPerm.txt", sep=""), sep="/t", row.names=F)
  
  df <- as.data.frame(t1)
  names(df)[1] <- "iHS"
  
  CalcPercent <- df %>% filter(iHS>=AverageiHSFemale)
  HighPercent <-nrow(CalcPercent)/10000
  
  print(paste("Female high Percentile is", HighPercent))
  
  #average iHS permutation male
  
  t1 <- rep(NA, 10000)
  
  for(j in 1:10000)
  {
    table1 <- PhenoSNPs[sample(nrow(PhenoSNPs), TotalNumberSNPsM),]
    t1[j] <- mean(table1$AbsStiHS)
  }
  
  write.table(t1, file=paste(PhenotypeName, "MaleiHSPerm.txt", sep=""), sep="\t", row.names=F)
  
  df <- as.data.frame(t1)
  names(df)[1] <- "iHS"
  
  CalcPercent <- df %>% filter(iHS>=AverageiHSMale)
  HighPercent <-nrow(CalcPercent)/10000
  
  print(paste("Male high Percentile is", HighPercent, sep = ' '))
  
}

getiHSPerm("HeightPdiffiHSScores.txt", "HeightiHSScores.txt")
getiHSPerm("HipCircPdiffiHSScores.txt", "HipCirciHSScores.txt")
getiHSPerm("WaistCircPdiffiHSScores.txt", "WaistCirciHSScores.txt")
getiHSPerm("BodyFatPercentPdiffiHSScores.txt", "BodyFatPercentiHSScores.txt")
getiHSPerm("BodyMassPdiffiHSScores.txt", "BodyMassiHSScores.txt")
