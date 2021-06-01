#!/usr/bin/env Rscript

#===================================================================================
#--Enrichment permutation analysis of all anthropometric phenotype-associated SNPs
#===================================================================================

library(dplyr)

#each of the dataframes of SexDiff SNPs
Height <- read.table("HeightwithPDIFF.FDR.tsv", header =T)
BodyMass <- read.table("BodyMasswithPDIFF.FDR.tsv", header =T)
HipCirc <- read.table("HipCircwithPDIFF.FDR.tsv", header =T)
WaistCirc <- read.table("WaistCircwithPDIFF.FDR.tsv", header =T)
BodyFatPercent <- read.table("BodyFatPercentwithPDIFF.FDR.tsv", header =T)

#All SexDiff-associated SNPs associated with at least one of our phenotypes
SNPs <- rbind(Height, BodyMass, HipCirc, WaistCirc, BodyFatPercent)

write.table(SNPs, "AllSigAnthroSNPs.txt", sep='\t', row.names = F)

#All genes associated with a GO term (n=17946)
AllGenesTable <- read.table("AllGenesGRCh37.txt", header = TRUE)

#Make sure all variables are numeric
AllGenesTable$TenThousandStart <- as.numeric(AllGenesTable$TenThousandStart)
AllGenesTable$TenThousandEnd <- as.numeric(AllGenesTable$TenThousandEnd)
AllGenesTable$CHR<-as.numeric(AllGenesTable$CHR)
AllGenesTable$Gene<-as.character(AllGenesTable$Gene)

#Genes associated with sexual differentiation (n=259)
SDGenesTable <- read.table("GOSexDifferenceGRCh37.txt", header = TRUE)

#Make sure all variables are numeric
SDGenesTable$TenThousandStart <- as.numeric(SDGenesTable$TenThousandStart)
SDGenesTable$TenThousandEnd <- as.numeric(SDGenesTable$TenThousandEnd)
SDGenesTable$CHR<-as.numeric(SDGenesTable$CHR)
SDGenesTable$Gene<-as.character(SDGenesTable$Gene)

#Make a column that will tell us if a gene is involved in sexual differentiation (SD) or not (O)
AllGenesTable[,"GeneGroup"] <- NA

for(i in 1:nrow(AllGenesTable)){ 
  Gene.Overlap = which(AllGenesTable$Gene[i] == SDGenesTable$Gene)
  if(length(Gene.Overlap) >=1){
    AllGenesTable$GeneGroup[i]="SD"
  } else {AllGenesTable$GeneGroup[i]="O"}
}

#Make sure all variables are numeric
SNPs$POS <- as.numeric(SNPs$POS)
SNPs$CHR<-as.numeric(SNPs$CHR)

#Finding the non sexual differentiation genes
Intervals <- AllGenesTable 

#makes a new column in our dataframe, which we will put if the SNP is in a gene, and if so, which one
Intervals[,"HaveSNP"] <- NA

for (i in 1:nrow(Intervals)){
  SNP.Overlap = which((Intervals$CHR[i] == SNPs$CHR) & (Intervals$TenThousandStart[i] <= SNPs$POS) &(Intervals$TenThousandEnd[i]>=SNPs$POS))
  if(length(SNP.Overlap)>=1){
    Intervals$HaveSNP[i] = 1
  } else {Intervals$HaveSNP[i] = 0}
}

WithSNP <- Intervals %>% filter(HaveSNP ==1)

############################################################
#Identify SNPs at our first FDR level
AllPdiff4 <- SNPs %>% filter(pdiff.FDR <= 0.001)
AllUniqueSigSNPs <- subset(AllPdiff4, !duplicated(rsid))

WithSNP[,"HaveSNP"] <- NA

for (i in 1:nrow(WithSNP)){
  SNP.Overlap = which((WithSNP$CHR[i] == AllPdiff4$CHR) & (WithSNP$TenThousandStart[i] <= AllPdiff4$POS) &(WithSNP$TenThousandEnd[i]>=AllPdiff4$POS))
  if(length(SNP.Overlap)>=1){
    WithSNP$HaveSNP[i] = 1
  } else {WithSNP$HaveSNP[i] = 0}
}

WithPdiffSNP <- WithSNP %>% filter(HaveSNP ==1)
SDGenes <- WithPdiffSNP %>% filter(GeneGroup=="SD")
NumberSD <- nrow(SDGenes)

OGenes <- WithPdiffSNP %>% filter(GeneGroup=="O")
NumberO <- nrow(OGenes)

TotalGenes <- NumberSD + NumberO

#Identify the genes that SNPs that are associated with a phenotype, but not at an FDR threshold, are in
#Remove SNPs in AllPdiff4
SNPsNOPdiff4 <-SNPs[!(SNPs$rsid %in% AllPdiff4$rsid),]

WithSNP[,"HaveSNP"] <- NA

for (i in 1:nrow(WithSNP)){
  SNP.Overlap = which((WithSNP$CHR[i] == SNPsNOPdiff4$CHR) & (WithSNP$TenThousandStart[i] <= SNPsNOPdiff4$POS) &(WithSNP$TenThousandEnd[i]>=SNPsNOPdiff4$POS))
  if(length(SNP.Overlap)>=1){
    WithSNP$HaveSNP[i] = 1
  } else {WithSNP$HaveSNP[i] = 0}
}

WithPdiffSNP <- WithSNP %>% filter(HaveSNP ==1)
SDGenes <- WithPdiffSNP %>% filter(GeneGroup=="SD")
nrow(SDGenes)

OGenes <- WithPdiffSNP %>% filter(GeneGroup=="O")
nrow(OGenes)

#Permutation to see if there is an overrepresentation of SNPs in sexual differentiation genes at this FDR level than we would see by chance
PermutationTable <- rep(NA, 10000)

for(j in 1:10000){
  table1 <- WithSNP[sample(nrow(WithSNP), TotalGenes),]
  table2 <- table1 %>% filter(GeneGroup =="SD")
  PermutationTable[j] <- nrow(table2)
}

PermutationTable2 <- data.frame(PermutationTable)
names(PermutationTable2)[1] <- "NumberGenes"

#give the P-value that our number of observed sexual differentiation genes (n=9) at this FDR was due to chance
CalcPercent <- PermutationTable2 %>% filter(NumberGenes>=NumberSD)
nrow(CalcPercent)/10000

write.table(PermutationTable2, "Pdiff4AnthroPermutationTable.txt", sep="\t", row.names = F)

############################################################
#same analyses as above, repeated for FDR = 0.005
AllPdiff3 <- SNPs %>% filter(pdiff.FDR <= 0.005)
AllUniqueSigSNPs <- subset(AllPdiff3, !duplicated(rsid))

WithSNP[,"HaveSNP"] <- NA

for (i in 1:nrow(WithSNP)){
  SNP.Overlap = which((WithSNP$CHR[i] == AllPdiff3$CHR) & (WithSNP$TenThousandStart[i] <= AllPdiff3$POS) &(WithSNP$TenThousandEnd[i]>=AllPdiff3$POS))
  if(length(SNP.Overlap)>=1){
    WithSNP$HaveSNP[i] = 1
  } else {WithSNP$HaveSNP[i] = 0}
}

WithPdiffSNP <- WithSNP %>% filter(HaveSNP ==1)
SDGenes <- WithPdiffSNP %>% filter(GeneGroup=="SD")
NumberSD <- nrow(SDGenes)

OGenes <- WithPdiffSNP %>% filter(GeneGroup=="O")
NumberO <- nrow(OGenes)

TotalGenes <- NumberSD + NumberO

#Identify the genes that SNPs that are associated with a phenotype, but not at an FDR threshold, are in
#Remove SNPs in AllPdiff3
SNPsNOPdiff3 <-SNPs[!(SNPs$rsid %in% AllPdiff3$rsid),]

WithSNP[,"HaveSNP"] <- NA

for (i in 1:nrow(WithSNP)){
  SNP.Overlap = which((WithSNP$CHR[i] == SNPsNOPdiff3$CHR) & (WithSNP$TenThousandStart[i] <= SNPsNOPdiff3$POS) &(WithSNP$TenThousandEnd[i]>=SNPsNOPdiff3$POS))
  if(length(SNP.Overlap)>=1){
    WithSNP$HaveSNP[i] = 1
  } else {WithSNP$HaveSNP[i] = 0}
}

WithPdiffSNP <- WithSNP %>% filter(HaveSNP ==1)
SDGenes <- WithPdiffSNP %>% filter(GeneGroup=="SD")
nrow(SDGenes)

OGenes <- WithPdiffSNP %>% filter(GeneGroup=="O")
nrow(OGenes)

#Permutation to see if there is an overrepresentation of SNPs in sexual differentiation genes at this FDR level than we would see by chance
PermutationTable <- rep(NA, 10000)

for(j in 1:10000){
  table1 <- WithSNP[sample(nrow(WithSNP), TotalGenes),]
  table2 <- table1 %>% filter(GeneGroup =="SD")
  PermutationTable[j] <- nrow(table2)
}

PermutationTable2 <- data.frame(PermutationTable)
names(PermutationTable2)[1] <- "NumberGenes"

#give the P-value that our number of observed sexual differentiation genes (n=9) at this FDR was due to chance
CalcPercent <- PermutationTable2 %>% filter(NumberGenes>=NumberSD)
nrow(CalcPercent)/10000

write.table(PermutationTable2, "Pdiff3AnthroPermutationTable.txt", sep="\t", row.names = F)
############################################################

AllPdiff2 <- SNPs %>% filter(pdiff.FDR <= 0.01)
AllUniqueSigSNPs <- subset(AllPdiff2, !duplicated(rsid))
WithSNP[,"HaveSNP"] <- NA

for (i in 1:nrow(WithSNP)){
  SNP.Overlap = which((WithSNP$CHR[i] == AllPdiff2$CHR) & (WithSNP$TenThousandStart[i] <= AllPdiff2$POS) &(WithSNP$TenThousandEnd[i]>=AllPdiff2$POS))
  if(length(SNP.Overlap)>=1){
    WithSNP$HaveSNP[i] = 1
  } else {WithSNP$HaveSNP[i] = 0}
}

WithPdiffSNP <- WithSNP %>% filter(HaveSNP ==1)
SDGenes <- WithPdiffSNP %>% filter(GeneGroup=="SD")
NumberSD <- nrow(SDGenes)

OGenes <- WithPdiffSNP %>% filter(GeneGroup=="O")
NumberO <- nrow(OGenes)

TotalGenes <- NumberSD + NumberO

#Identify the genes that SNPs that are associated with a phenotype, but not at an FDR threshold, are in
#Remove SNPs in AllPdiff2
SNPsNOPdiff2 <-SNPs[!(SNPs$rsid %in% AllPdiff2$rsid),]

WithSNP[,"HaveSNP"] <- NA

for (i in 1:nrow(WithSNP)){
  SNP.Overlap = which((WithSNP$CHR[i] == SNPsNOPdiff2$CHR) & (WithSNP$TenThousandStart[i] <= SNPsNOPdiff2$POS) &(WithSNP$TenThousandEnd[i]>=SNPsNOPdiff2$POS))
  if(length(SNP.Overlap)>=1){
    WithSNP$HaveSNP[i] = 1
  } else {WithSNP$HaveSNP[i] = 0}
}

WithPdiffSNP <- WithSNP %>% filter(HaveSNP ==1)
SDGenes <- WithPdiffSNP %>% filter(GeneGroup=="SD")
nrow(SDGenes)

OGenes <- WithPdiffSNP %>% filter(GeneGroup=="O")
nrow(OGenes)

#Permutation to see if there is an overrepresentation of SNPs in sexual differentiation genes at this FDR level than we would see by chance
PermutationTable <- rep(NA, 10000)

for(j in 1:10000){
  table1 <- WithSNP[sample(nrow(WithSNP), TotalGenes),]
  table2 <- table1 %>% filter(GeneGroup =="SD")
  PermutationTable[j] <- nrow(table2)
}

PermutationTable2 <- data.frame(PermutationTable)
names(PermutationTable2)[1] <- "NumberGenes"

#give the P-value that our number of observed sexual differentiation genes (n=9) at this FDR was due to chance
CalcPercent <- PermutationTable2 %>% filter(NumberGenes>=NumberSD)
nrow(CalcPercent)/10000

write.table(PermutationTable2, "Pdiff2AnthroPermutationTable.txt", sep="\t", row.names = F)
############################################################

AllPdiff1 <- SNPs %>% filter(pdiff.FDR <= 0.05)
AllUniqueSigSNPs <- subset(AllPdiff1, !duplicated(rsid))
WithSNP[,"HaveSNP"] <- NA

for (i in 1:nrow(WithSNP)){
  SNP.Overlap = which((WithSNP$CHR[i] == AllPdiff1$CHR) & (WithSNP$TenThousandStart[i] <= AllPdiff1$POS) &(WithSNP$TenThousandEnd[i]>=AllPdiff1$POS))
  if(length(SNP.Overlap)>=1){
    WithSNP$HaveSNP[i] = 1
  } else {WithSNP$HaveSNP[i] = 0}
}

WithPdiffSNP <- WithSNP %>% filter(HaveSNP ==1)
SDGenes <- WithPdiffSNP %>% filter(GeneGroup=="SD")
NumberSD <- nrow(SDGenes)

OGenes <- WithPdiffSNP %>% filter(GeneGroup=="O")
NumberO <- nrow(OGenes)

TotalGenes <- NumberSD + NumberO

#Identify the genes that SNPs that are associated with a phenotype, but not at an FDR threshold, are in
#Remove SNPs in AllPdiff1
SNPsNOPdiff1 <-SNPs[!(SNPs$rsid %in% AllPdiff1$rsid),]

WithSNP[,"HaveSNP"] <- NA

for (i in 1:nrow(WithSNP)){
  SNP.Overlap = which((WithSNP$CHR[i] == SNPsNOPdiff1$CHR) & (WithSNP$TenThousandStart[i] <= SNPsNOPdiff1$POS) &(WithSNP$TenThousandEnd[i]>=SNPsNOPdiff1$POS))
  if(length(SNP.Overlap)>=1){
    WithSNP$HaveSNP[i] = 1
  } else {WithSNP$HaveSNP[i] = 0}
}

WithPdiffSNP <- WithSNP %>% filter(HaveSNP ==1)
SDGenes <- WithPdiffSNP %>% filter(GeneGroup=="SD")
nrow(SDGenes)

OGenes <- WithPdiffSNP %>% filter(GeneGroup=="O")
nrow(OGenes)

#Permutation to see if there is an overrepresentation of SNPs in sexual differentiation genes at this FDR level than we would see by chance
PermutationTable <- rep(NA, 10000)

for(j in 1:10000){
  table1 <- WithSNP[sample(nrow(WithSNP), TotalGenes),]
  table2 <- table1 %>% filter(GeneGroup =="SD")
  PermutationTable[j] <- nrow(table2)
}

PermutationTable2 <- data.frame(PermutationTable)
names(PermutationTable2)[1] <- "NumberGenes"

#give the P-value that our number of observed sexual differentiation genes (n=9) at this FDR was due to chance
CalcPercent <- PermutationTable2 %>% filter(NumberGenes>=NumberSD)
nrow(CalcPercent)/10000

write.table(PermutationTable2, "Pdiff1AnthroPermutationTable.txt", sep="\t", row.names = F)