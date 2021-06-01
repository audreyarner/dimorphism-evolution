#!/usr/bin/env Rscript

#========================================================================================
#--Identify differences in beta values for male & female SexDiff-associated SNPs & plots
#========================================================================================

library(dplyr)
library(ggplot2)

logBeta = function(PrunedSNPs, PhenoSNPs){
  PrunedPdiff <- read.table(PrunedSNPs, header = T)
  PhenotypeName <- gsub("PdiffPrunedSNPs.txt", "", PrunedSNPs)

  #get average trait-SDS for SNPs disproportionately impacting female height
  PrunedPdiffWomen <- PrunedPdiff %>% filter(Female.pval <= Male.pval)
  PrunedPdiffWomenBeta <- PrunedPdiffWomen %>% select(rsid, Female.beta, Male.beta, Female.tSDS, Male.tSDS)

  PrunedPdiffWomenBeta <- mutate(PrunedPdiffWomenBeta, Abs.Women = abs(PrunedPdiffWomenBeta$Female.beta))
  PrunedPdiffWomenBeta <- mutate(PrunedPdiffWomenBeta, Abs.Men = abs(PrunedPdiffWomenBeta$Male.beta))

  PrunedPdiffWomenBeta <- mutate(PrunedPdiffWomenBeta, Diff = (PrunedPdiffWomenBeta$Abs.Women-PrunedPdiffWomenBeta$Abs.Men))
  PrunedPdiffWomenBeta <- mutate(PrunedPdiffWomenBeta, Ratio = log2(PrunedPdiffWomenBeta$Abs.Women/PrunedPdiffWomenBeta$Abs.Men))

  PrunedPdiffWomenBeta$Pheno <- PhenotypeName
  PrunedPdiffWomenBeta$Sex <- "Female"
  
  print(mean(PrunedPdiffWomenBeta$Ratio))

  write.table(PrunedPdiffWomenBeta, file=paste(PhenotypeName, "FemaleLogBeta.txt", sep=""), sep="\t", row.names=F)

  #get average tSDS for SNPs disproportionately impacting male height
  PrunedPdiffMen <- PrunedPdiff %>% filter(Male.pval <= Female.pval)
  PrunedPdiffMenBeta <- PrunedPdiffMen %>% select(rsid, Female.beta, Male.beta, Female.tSDS, Male.tSDS)

  PrunedPdiffMenBeta <- mutate(PrunedPdiffMenBeta, Abs.Women = abs(PrunedPdiffMenBeta$Female.beta))
  PrunedPdiffMenBeta <- mutate(PrunedPdiffMenBeta, Abs.Men = abs(PrunedPdiffMenBeta$Male.beta))

  PrunedPdiffMenBeta <- mutate(PrunedPdiffMenBeta, Diff = (PrunedPdiffMenBeta$Abs.Men-PrunedPdiffMenBeta$Abs.Women))
  PrunedPdiffMenBeta <- mutate(PrunedPdiffMenBeta, Ratio = log2(PrunedPdiffMenBeta$Abs.Women/PrunedPdiffMenBeta$Abs.Men))

  PrunedPdiffMenBeta$Pheno <- PhenotypeName
  PrunedPdiffMenBeta$Sex <- "Male"

  print(mean(PrunedPdiffMenBeta$Ratio))
  
  write.table(PrunedPdiffMenBeta, file=paste(PhenotypeName, "MaleLogBeta.txt", sep=""), sep="\t", row.names=F)

  #get average tSDS for SNPs disproportionately impacting male height
  PrunedPheno <- read.table(PhenoSNPs, header = T)
  PrunedPheno <- PrunedPheno %>% select(rsid,Female.beta, Male.beta, Female.tSDS, Male.tSDS)
  
  PrunedPheno <- mutate(PrunedPheno, Abs.Women = abs(PrunedPheno$Female.beta))
  PrunedPheno <- mutate(PrunedPheno, Abs.Men = abs(PrunedPheno$Male.beta))
  
  PrunedPheno <- mutate(PrunedPheno, Diff = (PrunedPheno$Abs.Men-PrunedPheno$Abs.Women))
  PrunedPheno <- mutate(PrunedPheno, Ratio = log2(PrunedPheno$Abs.Women/PrunedPheno$Abs.Men))
  
  PrunedPheno$Pheno <- PhenotypeName
  PrunedPheno$Sex <- "All"
  
  write.table(PrunedPheno, file=paste(PhenotypeName, "AllLogBeta.txt", sep=""), sep="\t", row.names=F)
  
}

logBeta("HeightPdiffPrunedSNPs.txt", "HeightPrunedSNPs.txt")
logBeta("HipCircPdiffPrunedSNPs.txt", "HipCircPrunedSNPs.txt")
logBeta("WaistCircPdiffPrunedSNPs.txt", "WaistCircPrunedSNPs.txt")
logBeta("BodyFatPercentPdiffPrunedSNPs.txt", "BodyFatPercentPrunedSNPs.txt")
logBeta("BodyMassPdiffPrunedSNPs.txt", "BodyMassPrunedSNPs.txt")

HeightF<-read.table("HeightFemaleLogBeta.txt", header =T)
HeightM<-read.table("HeightMaleLogBeta.txt", header =T)
HeightAll<-read.table("HeightAllLogBeta.txt", header =T)
HipCircF<-read.table("HipCircFemaleLogBeta.txt", header =T)
HipCircM<-read.table("HipCircMaleLogBeta.txt", header =T)
HipCircAll<-read.table("HipCircAllLogBeta.txt", header =T)
WaistCircF<-read.table("WaistCircFemaleLogBeta.txt", header =T)
WaistCircM<-read.table("WaistCircMaleLogBeta.txt", header =T)
WaistCircAll<-read.table("WaistCircAllLogBeta.txt", header =T)
BodyFatPercentF<-read.table("BodyFatPercentFemaleLogBeta.txt", header =T)
BodyFatPercentM<-read.table("BodyFatPercentMaleLogBeta.txt", header =T)
BodyFatPercentAll<-read.table("BodyFatPercentAllLogBeta.txt", header =T)
BodyMassF<-read.table("BodyMassFemaleLogBeta.txt", header =T)
BodyMassM<-read.table("BodyMassMaleLogBeta.txt", header =T)
BodyMassAll<-read.table("BodyMassAllLogBeta.txt", header =T)

AllPhenoBeta <- rbind(HeightAll, HeightF, HeightM, HipCircAll, HipCircF, HipCircM, WaistCircAll, WaistCircF, WaistCircM, BodyFatPercentAll, BodyFatPercentF, BodyFatPercentM, BodyMassAll, BodyMassF, BodyMassM)

AllPhenoBeta$Pheno <- factor(AllPhenoBeta$Pheno, levels=c("Height", "BodyMass", "HipCirc", "BodyFatPercent", "WaistCirc"))

Plot <- ggplot(AllPhenoBeta, aes(x=Pheno,y=Ratio, fill=Sex, color = Sex)) + 
  geom_hline(yintercept=0, color = "#efefef")+
  geom_violin(lwd=.75, width = .7) +
  geom_point(position = position_jitterdodge(seed = 1, dodge.width = 0.7)) +
  scale_fill_manual(values=c("#e4f1e1", "#f3cbd3", "#f8ea73")) + 
  scale_color_manual(values=c("#0d585f", "#c94f7c", "#DDCC77")) +
  theme_classic() +
  ylab("Log2(Ratio)") +
  stat_summary(fun=median, geom="point", size = 3) #, color = c("#0d585f", "#c94f7c", "#DDCC77","#0d585f", "#c94f7c", "#DDCC77","#0d585f", "#c94f7c", "#DDCC77","#0d585f", "#c94f7c", "#DDCC77", "#DDCC77","#0d585f", "#c94f7c"))

Plot

ggsave(Plot, file="2BetaViolin8-24.pdf", height=6, width =12, useDingbats=FALSE)
