#!/usr/bin/env Rscript

#===================================================================================
#--Enrichment permutation analysis of all anthropometric phenotype-associated SNPs
#===================================================================================

library(dplyr)
library(ggplot2)

#each of the dataframes of SexDiff SNPs
Height <- read.table("HeightwithPDIFF.FDR.tsv", header =T)
BodyMass <- read.table("BodyMasswithPDIFF.FDR.tsv", header =T)
HipCirc <- read.table("HipCircwithPDIFF.FDR.tsv", header =T)
WaistCirc <- read.table("WaistCircwithPDIFF.FDR.tsv", header =T)
BodyFatPercent <- read.table("BodyFatPercentwithPDIFF.FDR.tsv", header =T)

#All SNPs associated with at least one of our phenotypes
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
  if(length(Gene.Overlap >=1)){
    AllGenesTable$GeneGroup[i]="SD"
  } else {AllGenesTable$GeneGroup[i]="O"}
}

#Make sure all variables are numeric
SNPs$POS <- as.numeric(SNPs$POS)
SNPs$CHR<-as.numeric(SNPs$CHR)

#Finding the non sexual differentiation genes
Intervals <- AllGenesTable %>% filter(GeneGroup == "O")

#makes a new column in our dataframe, which we will put if the SNP is in a gene, and if so, which one
SNPs[,"GeneOverlap"] <- NA

for (i in 1:nrow(SNPs)){
  #which() will give the position of the gene that overlaps, which we can then use to find the gene name
  SNP.Overlap = which((Intervals$CHR == SNPs$CHR[i]) & (Intervals$TenThousandStart <= SNPs$POS[i]) &(Intervals$TenThousandEnd>=SNPs$POS[i]))
  if(length(SNP.Overlap>=1)){
    SNPs$GeneOverlap[i] = Intervals$Gene[SNP.Overlap]
  } else {SNPs$GeneOverlap[i] = 0}
}

#Number of unique non-sexual differentiation genes represented by all SNPs (n=2138)
#1 is subtracted for 0
length(unique(SNPs$GeneOverlap)) - 1

#Make dataframe of those unique gene names
NonSDGenes <- data.frame(unique(SNPs$GeneOverlap))
names(NonSDGenes)[1] <- "GeneName"
NonSDGenes[,"GeneType"] <- "O"
NonSDGenes <- NonSDGenes %>% filter(GeneName != 0)

#Finding the sexual differentiation genes
Intervals <- AllGenesTable %>% filter(GeneGroup == "SD")

#makes a new column in our dataframe, which we will put if the SNP is in a gene, and if so, which one
SNPs[,"GeneOverlap"] <- NA

for (i in 1:nrow(SNPs)){
  #which() will give the position of the gene that overlaps, which we can then use to find the gene name
  SNP.Overlap = which((Intervals$CHR == SNPs$CHR[i]) & (Intervals$TenThousandStart <= SNPs$POS[i]) &(Intervals$TenThousandEnd>=SNPs$POS[i]))
  if(length(SNP.Overlap>=1)){
    SNPs$GeneOverlap[i] = Intervals$Gene[SNP.Overlap]
  } else {SNPs$GeneOverlap[i] = 0}
}

#Number of unique non-sexual differentiation genes represented by all SNPs (n=53)
#1 is subtracted for 0
length(unique(SNPs$GeneOverlap)) - 1

#Make dataframe of those unique gene names
SDGenes <- data.frame(unique(SNPs$GeneOverlap))
names(SDGenes)[1] <- "GeneName"
SDGenes[,"GeneType"] <- "SD"
SDGenes <- SDGenes %>% filter(GeneName != 0)

AllGenes <- rbind(SDGenes, NonSDGenes)

############################################################
#Identify SNPs at our first FDR level
AllPdiff4 <- SNPs %>% filter(pdiff.FDR <= 0.001)
AllUniqueSigSNPs <- subset(AllPdiff4, !duplicated(rsid))

Intervals <- AllGenesTable %>% filter(GeneGroup == "O")

#makes a new column in our dataframe, which we will put if the SNP is in a gene, and if so, which one
AllPdiff4[,"GeneOverlap"] <- NA

for (i in 1:nrow(AllPdiff4)){
  #which() will give the position of the gene that overlaps, which we can then use to find the gene name
  SNP.Overlap = which((Intervals$CHR == AllPdiff4$CHR[i]) & (Intervals$TenThousandStart <= AllPdiff4$POS[i]) &(Intervals$TenThousandEnd>=AllPdiff4$POS[i]))
  if(length(SNP.Overlap>=1)){
    AllPdiff4$GeneOverlap[i] = Intervals$Gene[SNP.Overlap]
  } else {AllPdiff4$GeneOverlap[i] = 0}
}

#n=126
NumberO <- length(unique(AllPdiff4$GeneOverlap)) - 1

Intervals <- AllGenesTable %>% filter(GeneGroup == "SD")

#makes a new column in our dataframe, which we will put if the SNP is in a gene, and if so, which one
AllPdiff4[,"GeneOverlap"] <- NA

for (i in 1:nrow(AllPdiff4)){
  #which() will give the position of the gene that overlaps, which we can then use to find the gene name
  SNP.Overlap = which((Intervals$CHR == AllPdiff4$CHR[i]) & (Intervals$TenThousandStart <= AllPdiff4$POS[i]) &(Intervals$TenThousandEnd>=AllPdiff4$POS[i]))
  if(length(SNP.Overlap>=1)){
    AllPdiff4$GeneOverlap[i] = Intervals$Gene[SNP.Overlap]
  } else {AllPdiff4$GeneOverlap[i] = 0}
}

#n=9
NumberSD <- length(unique(AllPdiff4$GeneOverlap)) - 1

TotalGenes <- NumberSD + NumberO

#Identify the genes that SNPs that are associated with a phenotype, but not at an FDR threshold, are in
#Remove SNPs in AllPdiff4
SNPsNOPdiff4 <-SNPs[!(SNPs$rsid %in% AllPdiff4$rsid),]

Intervals <- AllGenesTable %>% filter(GeneGroup == "O")

#makes a new column in our dataframe, which we will put if the SNP is in a gene, and if so, which one
SNPsNOPdiff4[,"GeneOverlap"] <- NA

for (i in 1:nrow(SNPsNOPdiff4)){
  #which() will give the position of the gene that overlaps, which we can then use to find the gene name
  SNP.Overlap = which((Intervals$CHR == SNPsNOPdiff4$CHR[i]) & (Intervals$TenThousandStart <= SNPsNOPdiff4$POS[i]) &(Intervals$TenThousandEnd>=SNPsNOPdiff4$POS[i]))
  if(length(SNP.Overlap>=1)){
    SNPsNOPdiff4$GeneOverlap[i] = Intervals$Gene[SNP.Overlap]
  } else {SNPsNOPdiff4$GeneOverlap[i] = 0}
}

#n=2115
length(unique(SNPsNOPdiff4$GeneOverlap)) - 1

Intervals <- AllGenesTable %>% filter(GeneGroup == "SD")

SNPsNOPdiff4[,"GeneOverlap"] <- NA

for (i in 1:nrow(SNPsNOPdiff4)){
  #which() will give the position of the gene that overlaps, which we can then use to find the gene name
  SNP.Overlap = which((Intervals$CHR == SNPsNOPdiff4$CHR[i]) & (Intervals$TenThousandStart <= SNPsNOPdiff4$POS[i]) &(Intervals$TenThousandEnd>=SNPsNOPdiff4$POS[i]))
  if(length(SNP.Overlap>=1)){
    SNPsNOPdiff4$GeneOverlap[i] = Intervals$Gene[SNP.Overlap]
  } else {SNPsNOPdiff4$GeneOverlap[i] = 0}
}

#n=52
length(unique(SNPsNOPdiff4$GeneOverlap)) - 1

#Permutation to see if there is an overrepresentation of SNPs in sexual differentiation genes at this FDR level than we would see by chance
PermutationTable <- rep(NA, 10000)

for(j in 1:10000){
  #take 135 genes because that was the total number of genes represented by SNPs from Pdiff 4
  table1 <- AllGenes[sample(nrow(AllGenes), TotalGenes),]
  table2 <- table1 %>% filter(GeneType =="SD")
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
nrow(AllUniqueSigSNPs)

Intervals <- AllGenesTable %>% filter(GeneGroup == "O")

#makes a new column in our dataframe, which we will put if the SNP is in a gene, and if so, which one
AllPdiff3[,"GeneOverlap"] <- NA

for (i in 1:nrow(AllPdiff3)){
  #which() will give the position of the gene that overlaps, which we can then use to find the gene name
  SNP.Overlap = which((Intervals$CHR == AllPdiff3$CHR[i]) & (Intervals$TenThousandStart <= AllPdiff3$POS[i]) &(Intervals$TenThousandEnd>=AllPdiff3$POS[i]))
  if(length(SNP.Overlap>=1)){
    AllPdiff3$GeneOverlap[i] = Intervals$Gene[SNP.Overlap]
  } else {AllPdiff3$GeneOverlap[i] = 0}
}

NumberO <- length(unique(AllPdiff3$GeneOverlap)) - 1

Intervals <- AllGenesTable %>% filter(GeneGroup == "SD")

#makes a new column in our dataframe, which we will put if the SNP is in a gene, and if so, which one
AllPdiff3[,"GeneOverlap"] <- NA

for (i in 1:nrow(AllPdiff3)){
  #which() will give the position of the gene that overlaps, which we can then use to find the gene name
  SNP.Overlap = which((Intervals$CHR == AllPdiff3$CHR[i]) & (Intervals$TenThousandStart <= AllPdiff3$POS[i]) &(Intervals$TenThousandEnd>=AllPdiff3$POS[i]))
  if(length(SNP.Overlap>=1)){
    AllPdiff3$GeneOverlap[i] = Intervals$Gene[SNP.Overlap]
  } else {AllPdiff3$GeneOverlap[i] = 0}
}

NumberSD <- length(unique(AllPdiff3$GeneOverlap)) - 1

TotalGenes <- NumberSD + NumberO

SNPsNOPdiff3 <-SNPs[!(SNPs$rsid %in% AllPdiff3$rsid),]

Intervals <- AllGenesTable %>% filter(GeneGroup == "O")

#makes a new column in our dataframe, which we will put if the SNP is in a gene, and if so, which one
SNPsNOPdiff3[,"GeneOverlap"] <- NA

for (i in 1:nrow(SNPsNOPdiff3)){
  #which() will give the position of the gene that overlaps, which we can then use to find the gene name
  SNP.Overlap = which((Intervals$CHR == SNPsNOPdiff3$CHR[i]) & (Intervals$TenThousandStart <= SNPsNOPdiff3$POS[i]) &(Intervals$TenThousandEnd>=SNPsNOPdiff3$POS[i]))
  if(length(SNP.Overlap>=1)){
    SNPsNOPdiff3$GeneOverlap[i] = Intervals$Gene[SNP.Overlap]
  } else {SNPsNOPdiff3$GeneOverlap[i] = 0}
}

#n=2075
length(unique(SNPsNOPdiff3$GeneOverlap)) - 1

Intervals <- AllGenesTable %>% filter(GeneGroup == "SD")

#makes a new column in our dataframe, which we will put if the SNP is in a gene, and if so, which one
SNPsNOPdiff3[,"GeneOverlap"] <- NA

for (i in 1:nrow(SNPsNOPdiff3)){
  #which() will give the position of the gene that overlaps, which we can then use to find the gene name
  SNP.Overlap = which((Intervals$CHR == SNPsNOPdiff3$CHR[i]) & (Intervals$TenThousandStart <= SNPsNOPdiff3$POS[i]) &(Intervals$TenThousandEnd>=SNPsNOPdiff3$POS[i]))
  if(length(SNP.Overlap>=1)){
    SNPsNOPdiff3$GeneOverlap[i] = Intervals$Gene[SNP.Overlap]
  } else {SNPsNOPdiff3$GeneOverlap[i] = 0}
}

#n=51
length(unique(SNPsNOPdiff3$GeneOverlap)) - 1

PermutationTable <- rep(NA, 10000)

for(j in 1:10000){
  table1 <- AllGenes[sample(nrow(AllGenes), TotalGenes),]
  table2 <- table1 %>% filter(GeneType =="SD")
  PermutationTable[j] <- nrow(table2)
}

PermutationTable2 <- data.frame(PermutationTable)
names(PermutationTable2)[1] <- "NumberGenes"

CalcPercent <- PermutationTable2 %>% filter(NumberGenes>=NumberSD)
nrow(CalcPercent)/10000

write.table(PermutationTable2, "Pdiff3AnthroPermutationTable.txt", sep="\t", row.names = F)

############################################################

AllPdiff2 <- SNPs %>% filter(pdiff.FDR <= 0.01)
AllUniqueSigSNPs <- subset(AllPdiff2, !duplicated(rsid))
nrow(AllUniqueSigSNPs)

Intervals <- AllGenesTable %>% filter(GeneGroup == "O")

#makes a new column in our dataframe, which we will put if the SNP is in a gene, and if so, which one
AllPdiff2[,"GeneOverlap"] <- NA

for (i in 1:nrow(AllPdiff2)){
  #which() will give the position of the gene that overlaps, which we can then use to find the gene name
  SNP.Overlap = which((Intervals$CHR == AllPdiff2$CHR[i]) & (Intervals$TenThousandStart <= AllPdiff2$POS[i]) &(Intervals$TenThousandEnd>=AllPdiff2$POS[i]))
  if(length(SNP.Overlap>=1)){
    AllPdiff2$GeneOverlap[i] = Intervals$Gene[SNP.Overlap]
  } else {AllPdiff2$GeneOverlap[i] = 0}
}

NumberO <- length(unique(AllPdiff2$GeneOverlap)) - 1

Intervals <- AllGenesTable %>% filter(GeneGroup == "SD")

#makes a new column in our dataframe, which we will put if the SNP is in a gene, and if so, which one
AllPdiff2[,"GeneOverlap"] <- NA

for (i in 1:nrow(AllPdiff2)){
  #which() will give the position of the gene that overlaps, which we can then use to find the gene name
  SNP.Overlap = which((Intervals$CHR == AllPdiff2$CHR[i]) & (Intervals$TenThousandStart <= AllPdiff2$POS[i]) &(Intervals$TenThousandEnd>=AllPdiff2$POS[i]))
  if(length(SNP.Overlap>=1)){
    AllPdiff2$GeneOverlap[i] = Intervals$Gene[SNP.Overlap]
  } else {AllPdiff2$GeneOverlap[i] = 0}
}

NumberSD <- length(unique(AllPdiff2$GeneOverlap)) - 1

TotalGenes <- NumberSD + NumberO

SNPsNOPdiff2 <-SNPs[!(SNPs$rsid %in% AllPdiff2$rsid),]

Intervals <- AllGenesTable %>% filter(GeneGroup == "O")

#makes a new column in our dataframe, which we will put if the SNP is in a gene, and if so, which one
SNPsNOPdiff2[,"GeneOverlap"] <- NA

for (i in 1:nrow(SNPsNOPdiff2)){
  #which() will give the position of the gene that overlaps, which we can then use to find the gene name
  SNP.Overlap = which((Intervals$CHR == SNPsNOPdiff2$CHR[i]) & (Intervals$TenThousandStart <= SNPsNOPdiff2$POS[i]) &(Intervals$TenThousandEnd>=SNPsNOPdiff2$POS[i]))
  if(length(SNP.Overlap>=1)){
    SNPsNOPdiff2$GeneOverlap[i] = Intervals$Gene[SNP.Overlap]
  } else {SNPsNOPdiff2$GeneOverlap[i] = 0}
}

length(unique(SNPsNOPdiff2$GeneOverlap)) - 1

Intervals <- AllGenesTable %>% filter(GeneGroup == "SD")

#makes a new column in our dataframe, which we will put if the SNP is in a gene, and if so, which one
SNPsNOPdiff2[,"GeneOverlap"] <- NA

for (i in 1:nrow(SNPsNOPdiff2)){
  #which() will give the position of the gene that overlaps, which we can then use to find the gene name
  SNP.Overlap = which((Intervals$CHR == SNPsNOPdiff2$CHR[i]) & (Intervals$TenThousandStart <= SNPsNOPdiff2$POS[i]) &(Intervals$TenThousandEnd>=SNPsNOPdiff2$POS[i]))
  if(length(SNP.Overlap>=1)){
    SNPsNOPdiff2$GeneOverlap[i] = Intervals$Gene[SNP.Overlap]
  } else {SNPsNOPdiff2$GeneOverlap[i] = 0}
}

length(unique(SNPsNOPdiff2$GeneOverlap)) - 1


PermutationTable <- rep(NA, 10000)

for(j in 1:10000){
  table1 <- AllGenes[sample(nrow(AllGenes), TotalGenes),]
  table2 <- table1 %>% filter(GeneType =="SD")
  PermutationTable[j] <- nrow(table2)
}

PermutationTable2 <- data.frame(PermutationTable)
names(PermutationTable2)[1] <- "NumberGenes"

CalcPercent <- PermutationTable2 %>% filter(NumberGenes>=NumberSD)
nrow(CalcPercent)/10000

write.table(PermutationTable2, "Pdiff2AnthroPermutationTable.txt", sep="\t", row.names = F)

############################################################

AllPdiff1 <- SNPs %>% filter(pdiff.FDR <= 0.05)
AllUniqueSigSNPs <- subset(AllPdiff1, !duplicated(rsid))
nrow(AllUniqueSigSNPs)

Intervals <- AllGenesTable %>% filter(GeneGroup == "O")

#makes a new column in our dataframe, which we will put if the SNP is in a gene, and if so, which one
AllPdiff1[,"GeneOverlap"] <- NA

for (i in 1:nrow(AllPdiff1)){
  #which() will give the position of the gene that overlaps, which we can then use to find the gene name
  SNP.Overlap = which((Intervals$CHR == AllPdiff1$CHR[i]) & (Intervals$TenThousandStart <= AllPdiff1$POS[i]) &(Intervals$TenThousandEnd>=AllPdiff1$POS[i]))
  if(length(SNP.Overlap>=1)){
    AllPdiff1$GeneOverlap[i] = Intervals$Gene[SNP.Overlap]
  } else {AllPdiff1$GeneOverlap[i] = 0}
}

NumberO <- length(unique(AllPdiff1$GeneOverlap)) - 1

Intervals <- AllGenesTable %>% filter(GeneGroup == "SD")

#makes a new column in our dataframe, which we will put if the SNP is in a gene, and if so, which one
AllPdiff1[,"GeneOverlap"] <- NA

for (i in 1:nrow(AllPdiff1)){
  #which() will give the position of the gene that overlaps, which we can then use to find the gene name
  SNP.Overlap = which((Intervals$CHR == AllPdiff1$CHR[i]) & (Intervals$TenThousandStart <= AllPdiff1$POS[i]) &(Intervals$TenThousandEnd>=AllPdiff1$POS[i]))
  if(length(SNP.Overlap>=1)){
    AllPdiff1$GeneOverlap[i] = Intervals$Gene[SNP.Overlap]
  } else {AllPdiff1$GeneOverlap[i] = 0}
}

NumberSD <- length(unique(AllPdiff1$GeneOverlap)) - 1

TotalGenes <- NumberSD + NumberO

SNPsNOPdiff1 <-SNPs[!(SNPs$rsid %in% AllPdiff1$rsid),]

Intervals <- AllGenesTable %>% filter(GeneGroup == "O")

#makes a new column in our dataframe, which we will put if the SNP is in a gene, and if so, which one
SNPsNOPdiff1[,"GeneOverlap"] <- NA

for (i in 1:nrow(SNPsNOPdiff1)){
  #which() will give the position of the gene that overlaps, which we can then use to find the gene name
  SNP.Overlap = which((Intervals$CHR == SNPsNOPdiff1$CHR[i]) & (Intervals$TenThousandStart <= SNPsNOPdiff1$POS[i]) &(Intervals$TenThousandEnd>=SNPsNOPdiff1$POS[i]))
  if(length(SNP.Overlap>=1)){
    SNPsNOPdiff1$GeneOverlap[i] = Intervals$Gene[SNP.Overlap]
  } else {SNPsNOPdiff1$GeneOverlap[i] = 0}
}

length(unique(SNPsNOPdiff1$GeneOverlap)) - 1

Intervals <- AllGenesTable %>% filter(GeneGroup == "SD")

#makes a new column in our dataframe, which we will put if the SNP is in a gene, and if so, which one
SNPsNOPdiff1[,"GeneOverlap"] <- NA

for (i in 1:nrow(SNPsNOPdiff1)){
  #which() will give the position of the gene that overlaps, which we can then use to find the gene name
  SNP.Overlap = which((Intervals$CHR == SNPsNOPdiff1$CHR[i]) & (Intervals$TenThousandStart <= SNPsNOPdiff1$POS[i]) &(Intervals$TenThousandEnd>=SNPsNOPdiff1$POS[i]))
  if(length(SNP.Overlap>=1)){
    SNPsNOPdiff1$GeneOverlap[i] = Intervals$Gene[SNP.Overlap]
  } else {SNPsNOPdiff1$GeneOverlap[i] = 0}
}

length(unique(SNPsNOPdiff1$GeneOverlap)) - 1


PermutationTable <- rep(NA, 10000)

for(j in 1:10000){
  table1 <- AllGenes[sample(nrow(AllGenes), TotalGenes),]
  table2 <- table1 %>% filter(GeneType =="SD")
  PermutationTable[j] <- nrow(table2)
}

PermutationTable2 <- data.frame(PermutationTable)
names(PermutationTable2)[1] <- "NumberGenes"

CalcPercent <- PermutationTable2 %>% filter(NumberGenes>=NumberSD)
nrow(CalcPercent)/10000

write.table(PermutationTable2, "Pdiff1AnthroPermutationTable.txt", sep="\t" row.names = F)