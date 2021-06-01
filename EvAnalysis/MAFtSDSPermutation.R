#!/usr/bin/env Rscript

#===================================================================================
#--trait-SDS (tSDS) analysis selecting for close MAF
#===================================================================================

library(dplyr)  

getMAFtSDS = function(MAFPdiff,MAFPheno){
  #read in both SexDiff and Phenotype-associated SNPs
  Pdiff <- read.table(MAFPdiff, header =T)
  Pheno <- read.table(MAFPheno, header = T)
  
  PhenotypeName <- gsub("MAF.txt","",MAFPheno)
  
  Pheno$rsid <- as.character(Pheno$rsid)
  
  #filter for SexDiff SNPs disproportionately impacting females
  PdiffF <- Pdiff %>% filter(Female.pval <= Male.pval)
  
  AveragetSDSFemale <- mean(PdiffF$Female.tSDS)
  
  #permutation analysis 
  SDS <- rep(NA, 10000)
  
  for(j in 1:10000){
    t1<-rep(NA,nrow(PdiffF))
    t2<-rep(NA,nrow(PdiffF))
    rows <- sample(nrow(PdiffF))
    #shuffle SexDiffSNPs so that they aren't in the same order each time to be matched to other SNPs with a close MAF
    PdiffFShuffle <- Pdiff[rows,]
    Pheno2 <- Pheno
    for(i in 1:nrow(PdiffFShuffle)){
      Pheno2 <- mutate(Pheno2, Absmaf=abs((PdiffFShuffle$Female.maf[i]-Pheno2$Female.maf)))
      #identify SNPs with MAF within 0.05
      Pheno1 <- Pheno2 %>% filter(Absmaf <= 0.05)
      #randomly select one of these SNPs 
      Selection <- Pheno1[sample(nrow(Pheno1), 1),]
      t1[i]<-Selection$Female.tSDS[1]
      t2[i]<-Selection$rsid[1]
      #take out SNP that was randomly selected within 0.05 MAF so that it isn't selected again
      Pheno2<-Pheno2[!(Pheno2$rsid %in% t2),]
    }
    SDS[j]<-mean(t1)
  }
  
  write.table(SDS, file=paste(PhenotypeName,"FemaleMAFtSDSPerm.txt", sep=""), sep="\t", row.names=F)
  
  df <- as.data.frame(SDS)
  names(df)[1] <- "tSDS"
  
  #identify percent of SNPs with trait-SDS above the observed value
  CalcPercent <- df %>% filter(tSDS>=AveragetSDSFemale)
  HighPercent <-nrow(CalcPercent)/10000
  
  print(paste(PhenotypeName, "Female high Percentile is", HighPercent))
  
  #repeat, but for SNPs disproportionately impacting males
  PdiffM <- Pdiff %>% filter(Male.pval <= Female.pval)
  
  AveragetSDSMale <- mean(PdiffM$Male.tSDS)
  
  SDS <- rep(NA, 10000)
  
  for(j in 1:10000){
    t1<-rep(NA,nrow(PdiffM))
    t2<-rep(NA,nrow(PdiffM))
    rows <- sample(nrow(PdiffM))
    PdiffMShuffle <- Pdiff[rows,]
    Pheno2 <- Pheno
    for(i in 1:nrow(PdiffMShuffle)){
      Pheno2 <- mutate(Pheno2, Absmaf=abs((PdiffMShuffle$Male.maf[i]-Pheno2$Male.maf)))
      Pheno1 <- Pheno2 %>% filter(Absmaf <= 0.05)
      Selection <- Pheno1[sample(nrow(Pheno1), 1),]
      t1[i]<-Selection$Male.tSDS[1]
      t2[i]<-Selection$rsid[1]
      Pheno2<-Pheno2[!(Pheno2$rsid %in% t2),]
    }
    SDS[j]<-mean(t1)
  }
  
  write.table(SDS, file=paste(PhenotypeName,"MaleMAFtSDSPerm.txt", sep=""), sep="\t", row.names=F)
  
  df <- as.data.frame(SDS)
  names(df)[1] <- "tSDS"
  
  CalcPercent <- df %>% filter(tSDS>=AveragetSDSMale)
  HighPercent <-nrow(CalcPercent)/10000
  
  print(paste(PhenotypeName, "Male high Percentile is", HighPercent))
}

getMAFtSDS("HeightPdiffMAF.txt", "HeightMAF.txt")
getMAFtSDS("WaistCircPdiffMAF.txt", "WaistCircMAF.txt")
getMAFtSDS("HipCircPdiffMAF.txt", "HipCircMAF.txt")
getMAFtSDS("BodyMassPdiffMAF.txt", "BodyMassMAF.txt")
getMAFtSDS("BodyFatPercentPdiffMAF.txt", "BodyFatPercentMAF.txt")



