#!/usr/bin/env Rscript

#===================================================================================
#--Get iHS scores for each of the SNPs from 
#===================================================================================

library(dplyr)
library(tidyr)

chr1 <- read.table("chr1.1kg.p3.allPops.iHS.txt", header =T)
chr2 <- read.table("chr2.1kg.p3.allPops.iHS.txt", header =T)
chr3 <- read.table("chr3.1kg.p3.allPops.iHS.txt", header =T)
chr4 <- read.table("chr4.1kg.p3.allPops.iHS.txt", header =T)
chr5 <- read.table("chr5.1kg.p3.allPops.iHS.txt", header =T)
chr6 <- read.table("chr6.1kg.p3.allPops.iHS.txt", header =T)
chr7 <- read.table("chr7.1kg.p3.allPops.iHS.txt", header =T)
chr8 <- read.table("chr8.1kg.p3.allPops.iHS.txt", header =T)
chr9 <- read.table("chr9.1kg.p3.allPops.iHS.txt", header =T)
chr10 <- read.table("chr10.1kg.p3.allPops.iHS.txt", header =T)
chr11 <- read.table("chr11.1kg.p3.allPops.iHS.txt", header =T)
chr12 <- read.table("chr12.1kg.p3.allPops.iHS.txt", header =T)
chr13 <- read.table("chr13.1kg.p3.allPops.iHS.txt", header =T)
chr14 <- read.table("chr14.1kg.p3.allPops.iHS.txt", header =T)
chr15 <- read.table("chr15.1kg.p3.allPops.iHS.txt", header =T)
chr16 <- read.table("chr16.1kg.p3.allPops.iHS.txt", header =T)
chr17 <- read.table("chr17.1kg.p3.allPops.iHS.txt", header =T)
chr18 <- read.table("chr18.1kg.p3.allPops.iHS.txt", header =T)
chr19 <- read.table("chr19.1kg.p3.allPops.iHS.txt", header =T)
chr20 <- read.table("chr20.1kg.p3.allPops.iHS.txt", header =T)
chr21 <- read.table("chr21.1kg.p3.allPops.iHS.txt", header =T)
chr22 <- read.table("chr22.1kg.p3.allPops.iHS.txt", header =T)

#concatenate iHS values for all chromosomes into one data frame
AllChrom <- rbind(chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22)

write.table(AllChrom, "AllChromiHS.txt", row.names=F, quote=F)

iHS <- AllChrom

#function to identify the iHS value of each SNP 
getiHS = function(Pheno){
  
  Trait <- read.table(Pheno, header=T)
  
  PhenotypeName <- gsub("PrunedSNPs.txt", "", Pheno)
  
  iHS$stdIHS <- as.character(iHS$stdIHS)
  iHS$RSNUM <- as.character(iHS$RSNUM)
  Trait$rsid <- as.character(Trait$rsid)
  
  Trait[,"StiHS"] <- NA
  
  for (i in 1:nrow(Trait)){
    SNP.Overlap = which(iHS$RSNUM == Trait$rsid[i])
    if(length(SNP.Overlap)>=1){
      Trait$StiHS[i] = iHS$stdIHS[SNP.Overlap]
    }
  }
  
  Trait2 <- Trait %>% filter(StiHS != 0)
  print(paste(nrow(Trait2), sep=" "))

  if (nrow(Trait2)>=1){
    Trait3 <- separate_rows(Trait2, StiHS, sep = "\\|")
  
    Trait4 = Trait3[seq(19, nrow(Trait3), 26),]

    write.table(Trait4, file=paste(PhenotypeName, "iHSScores.txt", sep=""), sep="\t", row.names=F) 
    }
}

getiHS("WaistCircPrunedSNPs.txt")
getiHS("HipCircPrunedSNPs.txt")
getiHS("HeightPrunedSNPs.txt")
getiHS("BodyFatPercentPrunedSNPs.txt")
getiHS("BodyMassPrunedSNPs.txt")

getiHS("WaistCircPdiffPrunedSNPs.txt")
getiHS("HipCircPdiffPrunedSNPs.txt")
getiHS("HeightPdiffPrunedSNPs.txt")
getiHS("BodyFatPercentPdiffPrunedSNPs.txt")
getiHS("BodyMassPdiffPrunedSNPs.txt")



