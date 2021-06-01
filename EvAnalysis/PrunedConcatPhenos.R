#!/usr/bin/env Rscript

#------------------------------------------------------------------------------------------------------
#this script is used to get concatenated SexDiff- and phenotype-associated SNPs
#------------------------------------------------------------------------------------------------------

#read in tables of SexDiff-associated SNPs
WaistCirc <- read.table("WaistCircPdiffPrunedSNPs.txt", header =T)

#include phenotype name as a column
WaistCirc <- mutate(WaistCirc, Pheno="WaistCirc")

#Identify if SNPs are more significantly associated with females or males
for (i in 1:nrow(WaistCirc)){
  if(WaistCirc$Female.pval[i]<WaistCirc$Male.pval[i]){
    WaistCirc$Sex[i] = "Female"
  } else {WaistCirc$Sex[i] = "Male"}
}

Height <- read.table("HeightPdiffPrunedSNPs.txt", header =T)

Height <- mutate(Height, Pheno="Height")

for (i in 1:nrow(Height)){
  if(Height$Female.pval[i]<Height$Male.pval[i]){
    Height$Sex[i] = "Female"
  } else {Height$Sex[i] = "Male"}
}

HipCirc <- read.table("HipCircPdiffPrunedSNPs.txt", header =T)

HipCirc <- mutate(HipCirc, Pheno="HipCirc")

for (i in 1:nrow(HipCirc)){
  if(HipCirc$Female.pval[i]<HipCirc$Male.pval[i]){
    HipCirc$Sex[i] = "Female"
  } else {HipCirc$Sex[i] = "Male"}
}

BodyMass <- read.table("BodyMassPdiffPrunedSNPs.txt", header =T)

BodyMass <- mutate(BodyMass, Pheno="BodyMass")

for (i in 1:nrow(BodyMass)){
  if(BodyMass$Female.pval[i]<BodyMass$Male.pval[i]){
    BodyMass$Sex[i] = "Female"
  } else {BodyMass$Sex[i] = "Male"}
}

BodyFatPercent <- read.table("BodyFatPercentPdiffPrunedSNPs.txt", header =T)

BodyFatPercent <- mutate(BodyFatPercent, Pheno="BodyFatPercent")

for (i in 1:nrow(BodyFatPercent)){
  if(BodyFatPercent$Female.pval[i]<BodyFatPercent$Male.pval[i]){
    BodyFatPercent$Sex[i] = "Female"
  } else {BodyFatPercent$Sex[i] = "Male"}
}

All <- rbind(WaistCirc,Height,HipCirc,BodyMass,BodyFatPercent)
All <- All %>% select(rsid,CHR,POS,Pheno,Sex)

write.table(All, "SexDiffSNPs.txt", quote=F,row.names=F,sep="\t")

WaistCirc <- read.table("WaistCircPrunedSNPs.txt", header =T)

#include phenotype name as a column
WaistCirc <- mutate(WaistCirc, Pheno="WaistCirc")

BodyFatPercent <- read.table("BodyFatPercentPrunedSNPs.txt", header =T)

BodyFatPercent <- mutate(BodyFatPercent, Pheno="BodyFatPercent")

BodyMass <- read.table("BodyMassPrunedSNPs.txt", header =T)

BodyMass <- mutate(BodyMass, Pheno="BodyMass")

Height <- read.table("HeightPrunedSNPs.txt", header =T)

Height <- mutate(Height, Pheno="Height")

HipCirc <- read.table("HipCircPrunedSNPs.txt", header =T)

HipCirc <- mutate(HipCirc, Pheno="HipCirc")

All <- rbind(WaistCirc,Height,HipCirc,BodyMass,BodyFatPercent)

All <- All %>% select(rsid,CHR,POS,Pheno)

write.table(All,"AllPrunedPheno.txt", quote=F, row.names=F)
