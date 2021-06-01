#!/usr/bin/env Rscript

#===================================================================================
#--Violin plot of trait-SDS scores
#===================================================================================

library(ggplot2)
library(dplyr)

#data frames of pruned phenotype-associated SNPs
BodyMassAll <- read.table("BodyMassiHSScores.txt", header = T)
BodyFatPercentAll <- read.table("BodyFatPercentiHSScores.txt", header = T)
HeightAll <- read.table("HeightiHSScores.txt", header = T)
HipCircAll <- read.table("HipCirciHSScores.txt", header = T)
WaistCircAll <- read.table("WaistCirciHSScores.txt", header = T)


#data frames of pruned SexDiff-associated SNPs
BodyMassPruned <- read.table("BodyMassPdiffiHSScores.txt", header = T)
BodyFatPercentPruned <- read.table("BodyFatPercentPdiffiHSScores.txt", header = T)
HeightPruned <- read.table("HeightPdiffiHSScores.txt", header = T)
HipCircPruned <- read.table("HipCircPdiffiHSScores.txt", header = T)
WaistCircPruned <- read.table("WaistCircPdiffiHSScores.txt", header = T)

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

#select only relevant columns needed for Violin plot (traitSDS, Group, Phenotype)
HeightAll <- HeightAll %>% select(StiHS,Group,Pheno)
HipCircAll <- HipCircAll %>% select(StiHS,Group,Pheno)
BodyMassAll <- BodyMassAll %>% select(StiHS,Group,Pheno)
BodyFatPercentAll <- BodyFatPercentAll %>% select(StiHS,Group,Pheno)
WaistCircAll <- WaistCircAll %>% select(StiHS,Group,Pheno)

HeightPruned <- HeightPruned %>% select(StiHS,Group,Pheno)
HipCircPruned <- HipCircPruned %>% select(StiHS,Group,Pheno)
BodyMassPruned <- BodyMassPruned %>% select(StiHS,Group,Pheno)
BodyFatPercentPruned <- BodyFatPercentPruned %>% select(StiHS,Group,Pheno)
WaistCircPruned <- WaistCircPruned %>% select(StiHS,Group,Pheno)

#combine into one dataframe
AlltSDS <- rbind(HeightPruned, HipCircPruned, BodyMassPruned, BodyFatPercentPruned, WaistCircPruned, HeightAll, HipCircAll, BodyMassAll, BodyFatPercentAll, WaistCircAll)

AlltSDS <- mutate(AlltSDS, AbsStiHS=abs(AlltSDS$StiHS))

#turn phenotype into factors so they can be ordered according to order in manuscript
AlltSDS$Pheno <- factor(AlltSDS$Pheno, levels=c("Height", "BodyMass", "HipCirc", "BodyFatPercent", "WaistCirc"))

Plot <- ggplot(AlltSDS, aes(x=Pheno,y=AbsStiHS, fill=Group, color = Group)) + 
  geom_violin(lwd=.75, width = .7) +
  geom_point(position = position_jitterdodge(seed = 1, dodge.width = 0.7, jitter.width=0.15)) +
  scale_fill_manual(values=c("#d0d0d0", "#F6D79C", "#C0F3FF")) + 
  scale_color_manual(values=c("#3d3d3d", "#FF9F00", "#00B0EB")) +
  theme_classic() +
  stat_summary(fun=median, geom="point", size = 3) #, color = c("#0d585f", "#c94f7c", "#DDCC77","#0d585f", "#c94f7c", "#DDCC77","#0d585f", "#c94f7c", "#DDCC77","#0d585f", "#c94f7c", "#DDCC77", "#DDCC77","#0d585f", "#c94f7c"))

Plot

ggsave(Plot, file="ViolinPlotiHSBFP2-28.pdf", height=6, width =12, useDingbats = FALSE)

