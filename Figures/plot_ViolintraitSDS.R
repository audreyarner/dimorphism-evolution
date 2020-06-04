#!/usr/bin/env Rscript

#===================================================================================
#--Violin plot of trait-SDS scores
#===================================================================================

library(ggplot2)
library(dplyr)

#data frames of pruned phenotype-associated SNPs
BodyMassAll <- read.table("BodyMassPrunedSNPs.txt", header = T)
BodyFatPercentAll <- read.table("BodyFatPercentPrunedSNPs.txt", header = T)
HeightAll <- read.table("HeightPrunedSNPs.txt", header = T)
HipCircAll <- read.table("HipCircPrunedSNPs.txt", header = T)
WaistCircAll <- read.table("WaistCircPrunedSNPs.txt", header = T)


#data frames of pruned SexDiff-associated SNPs
BodyMassPruned <- read.table("BodyMassPdiffPrunedSNPs.txt", header = T)
BodyFatPercentPruned <- read.table("BodyFatPercentPdiffPrunedSNPs.txt", header = T)
HeightPruned <- read.table("HeightPdiffPrunedSNPs.txt", header = T)
HipCircPruned <- read.table("HipCircPdiffPrunedSNPs.txt", header = T)
WaistCircPruned <- read.table("WaistCircPdiffPrunedSNPs.txt", header = T)

#add column for phenotype and group, which will determine color/position in plot
BodyMassAll$Pheno <- "BodyMass"
BodyFatPercentAll$Pheno <- "BodyFatPercent"
HeightAll$Pheno <- "Height"
HipCircAll$Pheno <- "HipCirc"
WaistCircAll$Pheno <- "WaistCirc"

BodyMassPruned$Pheno <- "BodyMass"
BodyFatPercentPruned$Pheno <- "BodyFatPercent"
HeightPruned$Pheno <- "Height"
HipCircPruned$Pheno <- "HipCirc"
WaistCircPruned$Pheno <- "WaistCirc"

BodyMassAll$Group <- "All"
BodyFatPercentAll$Group <- "All"
HeightAll$Group <- "All"
HipCircAll$Group <- "All"
WaistCircAll$Group <- "All"

#Assigns SexDiff-associated SNP to group based on if the SNP is more significant in males or females
BodyMassPruned$Group <- NA

for (i in 1:nrow(BodyMassPruned)){
  if(BodyMassPruned$Female.pval[i]<BodyMassPruned$Male.pval[i]){
    BodyMassPruned$Group[i] = "Female"
  } else {BodyMassPruned$Group[i] = "Male"}
}

BodyFatPercentPruned$Group <- NA

for (i in 1:nrow(BodyFatPercentPruned)){
  if(BodyFatPercentPruned$Female.pval[i]<BodyFatPercentPruned$Male.pval[i]){
    BodyFatPercentPruned$Group[i] = "Female"
  } else {BodyFatPercentPruned$Group[i] = "Male"}
}

HeightPruned$Group <- NA

for (i in 1:nrow(HeightPruned)){
  if(HeightPruned$Female.pval[i]<HeightPruned$Male.pval[i]){
    HeightPruned$Group[i] = "Female"
  } else {HeightPruned$Group[i] = "Male"}
}

HipCircPruned$Group <- NA

for (i in 1:nrow(HipCircPruned)){
  if(HipCircPruned$Female.pval[i]<HipCircPruned$Male.pval[i]){
    HipCircPruned$Group[i] = "Female"
  } else {HipCircPruned$Group[i] = "Male"}
}


WaistCircPruned$Group <- NA

for (i in 1:nrow(WaistCircPruned)){
  if(WaistCircPruned$Female.pval[i]<WaistCircPruned$Male.pval[i]){
    WaistCircPruned$Group[i] = "Female"
  } else {WaistCircPruned$Group[i] = "Male"}
}

#Assigns tSDS score of SexDiff-associated SNP based on whether the female or male association value is more significant
BodyMassPruned$tSDS <- NA

for (i in 1:nrow(BodyMassPruned)){
  if(BodyMassPruned$Female.pval[i]<BodyMassPruned$Male.pval[i]){
    BodyMassPruned$tSDS[i] = BodyMassPruned$Female.tSDS[i]
  } else {BodyMassPruned$tSDS[i] = BodyMassPruned$Male.tSDS[i]}
}

BodyFatPercentPruned$tSDS <- NA

for (i in 1:nrow(BodyFatPercentPruned)){
  if(BodyFatPercentPruned$Female.pval[i]<BodyFatPercentPruned$Male.pval[i]){
    BodyFatPercentPruned$tSDS[i] = BodyFatPercentPruned$Female.tSDS[i]
  } else {BodyFatPercentPruned$tSDS[i] = BodyFatPercentPruned$Male.tSDS[i]}
}

HeightPruned$tSDS <- NA

for (i in 1:nrow(HeightPruned)){
  if(HeightPruned$Female.pval[i]<HeightPruned$Male.pval[i]){
    HeightPruned$tSDS[i] = HeightPruned$Female.tSDS[i]
  } else {HeightPruned$tSDS[i] = HeightPruned$Male.tSDS[i]}
}

HipCircPruned$tSDS <- NA

for (i in 1:nrow(HipCircPruned)){
  if(HipCircPruned$Female.pval[i]<HipCircPruned$Male.pval[i]){
    HipCircPruned$tSDS[i] = HipCircPruned$Female.tSDS[i]
  } else {HipCircPruned$tSDS[i] = HipCircPruned$Male.tSDS[i]}
}

WaistCircPruned$tSDS <- NA

for (i in 1:nrow(WaistCircPruned)){
  if(WaistCircPruned$Female.pval[i]<WaistCircPruned$Male.pval[i]){
    WaistCircPruned$tSDS[i] = WaistCircPruned$Female.tSDS[i]
  } else {WaistCircPruned$tSDS[i] = WaistCircPruned$Male.tSDS[i]}
}

#Assigns tSDS of phenotype-associated SNPs based on whether the male or female association value is more significant
BodyMassAll$tSDS <- NA

for (i in 1:nrow(BodyMassAll)){
  if(BodyMassAll$Female.pval[i]<BodyMassAll$Male.pval[i]){
    BodyMassAll$tSDS[i] = BodyMassAll$Female.tSDS[i]
  } else {BodyMassAll$tSDS[i] = BodyMassAll$Male.tSDS[i]}
}

BodyFatPercentAll$tSDS <- NA

for (i in 1:nrow(BodyFatPercentAll)){
  if(BodyFatPercentAll$Female.pval[i]<BodyFatPercentAll$Male.pval[i]){
    BodyFatPercentAll$tSDS[i] = BodyFatPercentAll$Female.tSDS[i]
  } else {BodyFatPercentAll$tSDS[i] = BodyFatPercentAll$Male.tSDS[i]}
}

HeightAll$tSDS <- NA

for (i in 1:nrow(HeightAll)){
  if(HeightAll$Female.pval[i]<HeightAll$Male.pval[i]){
    HeightAll$tSDS[i] = HeightAll$Female.tSDS[i]
  } else {HeightAll$tSDS[i] = HeightAll$Male.tSDS[i]}
}

HipCircAll$tSDS <- NA

for (i in 1:nrow(HipCircAll)){
  if(HipCircAll$Female.pval[i]<HipCircAll$Male.pval[i]){
    HipCircAll$tSDS[i] = HipCircAll$Female.tSDS[i]
  } else {HipCircAll$tSDS[i] = HipCircAll$Male.tSDS[i]}
}

WaistCircAll$tSDS <- NA

for (i in 1:nrow(WaistCircAll)){
  if(WaistCircAll$Female.pval[i]<WaistCircAll$Male.pval[i]){
    WaistCircAll$tSDS[i] = WaistCircAll$Female.tSDS[i]
  } else {WaistCircAll$tSDS[i] = WaistCircAll$Male.tSDS[i]}
}

#select only relevant columns needed for Violin plot (traitSDS, Group, Phenotype)
HeightAll <- HeightAll %>% select(tSDS,Group,Pheno)
HipCircAll <- HipCircAll %>% select(tSDS,Group,Pheno)
BodyMassAll <- BodyMassAll %>% select(tSDS,Group,Pheno)
BodyFatPercentAll <- BodyFatPercentAll %>% select(tSDS,Group,Pheno)
WaistCircAll <- WaistCircAll %>% select(tSDS,Group,Pheno)

HeightPruned <- HeightPruned %>% select(tSDS,Group,Pheno)
HipCircPruned <- HipCircPruned %>% select(tSDS,Group,Pheno)
BodyMassPruned <- BodyMassPruned %>% select(tSDS,Group,Pheno)
BodyFatPercentPruned <- BodyFatPercentPruned %>% select(tSDS,Group,Pheno)
WaistCircPruned <- WaistCircPruned %>% select(tSDS,Group,Pheno)

#combine into one dataframe
AlltSDS <- rbind(HeightPruned, HipCircPruned, BodyMassPruned, BodyFatPercentPruned, WaistCircPruned, HeightAll, HipCircAll, BodyMassAll, BodyFatPercentAll, WaistCircAll)

#turn phenotype into factors so they can be ordered according to order in manuscript
AlltSDS$Pheno <- factor(AlltSDS$Pheno, levels=c("Height", "BodyMass", "HipCirc", "BodyFatPercent", "WaistCirc"))

Plot <- ggplot(AlltSDS, aes(x=Pheno,y=tSDS, fill=Group, color = Group)) + 
  geom_hline(yintercept=0, color = "#efefef")+
  geom_violin(trim =FALSE, lwd=.75, width = .7) +
  ylim(-4,6) +
  scale_fill_manual(values=c("#e4f1e1", "#f3cbd3", "#f8ea73")) + 
  scale_color_manual(values=c("#0d585f", "#c94f7c", "#DDCC77")) +
  theme_classic() + 
  stat_summary(fun.y=median, geom="point", size = 3) #, color = c("#0d585f", "#c94f7c", "#DDCC77","#0d585f", "#c94f7c", "#DDCC77","#0d585f", "#c94f7c", "#DDCC77","#0d585f", "#c94f7c", "#DDCC77", "#DDCC77","#0d585f", "#c94f7c"))

Plot

ggsave(Plot, file="ViolinPlot.pdf", height=7, width =12, useDingbats = FALSE)

