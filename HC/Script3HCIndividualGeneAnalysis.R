#!/usr/bin/env Rscript

#-----------------------------------------------------------------------------------------
#this script is used to find SNPs within +/- 10000 base pairs (bp) of genes in GO:0007548
#-----------------------------------------------------------------------------------------

#packages needed for this script
library(dplyr)
library(ggplot2)

#reads in the table with our genes of interest
Intervals <- read.table("/storage/home/ama6560/scratch/GOwGeneNamesNoRepeats6-13.txt", header = TRUE)
#All SNPs that hit the genome wide significance threshold (5x10^-8)
EntireSample <- read.table("HipCircwithPDIFF.FDR.tsv", header = TRUE)

#makes sure formatting is all correct
Intervals$TenThousandStart <- as.numeric(Intervals$TenThousandStart)
Intervals$TenThousandEnd <- as.numeric(Intervals$TenThousandEnd)
Intervals$CHR<-as.numeric(Intervals$CHR)
Intervals$Gene<-as.character(Intervals$Gene)

EntireSample$POS <- as.numeric(EntireSample$POS)
EntireSample$CHR<-as.numeric(EntireSample$CHR)

#makes a new column in our dataframe, which we will put if the SNP is in a gene, and if so, which one
EntireSample[,"GeneOverlap"] <- NA

#puts the gene name next to SNPs that are within intervals for the entire dataset.
for (i in 1:nrow(EntireSample)){
  # print(i)
  #which() will give the position of the gene that overlaps, which we can then use to find the gene name
  SNP.Overlap = which((Intervals$CHR == EntireSample$CHR[i]) & (Intervals$TenThousandStart <= EntireSample$POS[i]) &(Intervals$TenThousandEnd>=EntireSample$POS[i]))
  if(length(SNP.Overlap>=1)){
    EntireSample$GeneOverlap[i] = Intervals$Gene[SNP.Overlap]
  } else {EntireSample$GeneOverlap[i] = 0}
}

PercentTable <- rep(NA, 4)
NumberGenes <- rep(NA, 4)

#------------------------------------------------------------------------------
#reads table of SNPs that passed a pdiff FDR of P=0.001
PdiffCutoff <- read.table("HipCircPdiff4.tsv", header = TRUE)

#makes sure the data is how we want it
PdiffCutoff$POS <- as.numeric(PdiffCutoff$POS)
PdiffCutoff$CHR<-as.numeric(PdiffCutoff$CHR)

#makes a new column in our dataframe, which we will put if the SNP is in a gene, and if so, which one
PdiffCutoff[,"GeneOverlap"] <- NA

#puts the gene name next to SNPs that are within intervals for the cutoff.
for (i in 1:nrow(PdiffCutoff)){
  SNP.Overlap = which((Intervals$CHR == PdiffCutoff$CHR[i]) & (Intervals$TenThousandStart <= PdiffCutoff$POS[i]) &(Intervals$TenThousandEnd>=PdiffCutoff$POS[i]))
  if(length(SNP.Overlap>=1)){
    PdiffCutoff$GeneOverlap[i] = Intervals$Gene[SNP.Overlap]
  } else {PdiffCutoff$GeneOverlap[i] = 0}
}

#this gives us the number of unique genes our SNPs of interest are in
signumber <- length(unique(PdiffCutoff$GeneOverlap)) - 1
NumberGenes[1] <- signumber

#permutation 10,000 times
GeneTable <- rep(NA, 10000) #list that we will put our number of genes in

for(i in 1:10000){
  table1 <- EntireSample[sample(nrow(EntireSample), nrow(PdiffCutoff)), ] #takes a random sample from the set of SNPs that passes the genome wide cutoff with the number of SNPs that is significant at out threshold
  GeneTable[i] <- (length(unique(table1$GeneOverlap)) - 1) #you have to subtract one because of the zero that's hanging around
}

write.table(GeneTable, "HipCircPdiff4Permutation.tsv")

FinishedTable <- read.table("HipCircPdiff4Permutation.tsv")

CalcPercent <- FinishedTable %>% filter(x>signumber)

PercentTable[1] <- (nrow(CalcPercent)/10000)

#------------------------------------------------------------------------------
#reads table of SNPs that passed a pdiff FDR of P=0.005
PdiffCutoff <- read.table("HipCircPdiff3.tsv", header = TRUE)

#makes sure the data is how we want it
PdiffCutoff$POS <- as.numeric(PdiffCutoff$POS)
PdiffCutoff$CHR<-as.numeric(PdiffCutoff$CHR)

#makes a new column in our dataframe, which we will put if the SNP is in a gene, and if so, which one
PdiffCutoff[,"GeneOverlap"] <- NA

#puts the gene name next to SNPs that are within intervals for the cutoff.
for (i in 1:nrow(PdiffCutoff)){
  SNP.Overlap = which((Intervals$CHR == PdiffCutoff$CHR[i]) & (Intervals$TenThousandStart <= PdiffCutoff$POS[i]) &(Intervals$TenThousandEnd>=PdiffCutoff$POS[i]))
  if(length(SNP.Overlap>=1)){
    PdiffCutoff$GeneOverlap[i] = Intervals$Gene[SNP.Overlap]
  } else {PdiffCutoff$GeneOverlap[i] = 0}
}

#this gives us the number of unique genes our SNPs of interest are in
signumber <- length(unique(PdiffCutoff$GeneOverlap)) - 1
NumberGenes[2] <- signumber

#permutation 10,000 times
GeneTable <- rep(NA, 10000) #list that we will put our number of genes in

for(i in 1:10000){
  table1 <- EntireSample[sample(nrow(EntireSample), nrow(PdiffCutoff)), ] #takes a random sample from the set of SNPs that passes the genome wide cutoff with the number of SNPs that is significant at out threshold
  GeneTable[i] <- (length(unique(table1$GeneOverlap)) - 1) #you have to subtract one because of the zero that's hanging around
}

write.table(GeneTable, "HipCircPdiff3Permutation.tsv")

FinishedTable <- read.table("HipCircPdiff3Permutation.tsv")

CalcPercent <- FinishedTable %>% filter(x>signumber)

PercentTable[2] <- (nrow(CalcPercent)/10000)

#------------------------------------------------------------------------------
#reads table of SNPs that passed a pdiff FDR of P=0.01
PdiffCutoff <- read.table("HipCircPdiff2.tsv", header = TRUE)

#makes sure the data is how we want it
PdiffCutoff$POS <- as.numeric(PdiffCutoff$POS)
PdiffCutoff$CHR<-as.numeric(PdiffCutoff$CHR)

#makes a new column in our dataframe, which we will put if the SNP is in a gene, and if so, which one
PdiffCutoff[,"GeneOverlap"] <- NA

#puts the gene name next to SNPs that are within intervals for the cutoff.
for (i in 1:nrow(PdiffCutoff)){
  SNP.Overlap = which((Intervals$CHR == PdiffCutoff$CHR[i]) & (Intervals$TenThousandStart <= PdiffCutoff$POS[i]) &(Intervals$TenThousandEnd>=PdiffCutoff$POS[i]))
  if(length(SNP.Overlap>=1)){
    PdiffCutoff$GeneOverlap[i] = Intervals$Gene[SNP.Overlap]
  } else {PdiffCutoff$GeneOverlap[i] = 0}
}

#this gives us the number of unique genes our SNPs of interest are in
signumber <- length(unique(PdiffCutoff$GeneOverlap)) - 1
NumberGenes[3] <- signumber

#permutation 10,000 times
GeneTable <- rep(NA, 10000) #list that we will put our number of genes in

for(i in 1:10000){
  table1 <- EntireSample[sample(nrow(EntireSample), nrow(PdiffCutoff)), ] #takes a random sample from the set of SNPs that passes the genome wide cutoff with the number of SNPs that is significant at out threshold
  GeneTable[i] <- (length(unique(table1$GeneOverlap)) - 1) #you have to subtract one because of the zero that's hanging around
}

write.table(GeneTable, "HipCircPdiff2Permutation.tsv")

FinishedTable <- read.table("HipCircPdiff2Permutation.tsv")

CalcPercent <- FinishedTable %>% filter(x>signumber)

PercentTable[3] <- (nrow(CalcPercent)/10000)

#------------------------------------------------------------------------------
#reads table of SNPs that passed a pdiff FDR of P=0.05
PdiffCutoff <- read.table("HipCircPdiff1.tsv", header = TRUE)

#makes sure the data is how we want it
PdiffCutoff$POS <- as.numeric(PdiffCutoff$POS)
PdiffCutoff$CHR<-as.numeric(PdiffCutoff$CHR)

#makes a new column in our dataframe, which we will put if the SNP is in a gene, and if so, which one
PdiffCutoff[,"GeneOverlap"] <- NA

#puts the gene name next to SNPs that are within intervals for the cutoff.
for (i in 1:nrow(PdiffCutoff)){
  SNP.Overlap = which((Intervals$CHR == PdiffCutoff$CHR[i]) & (Intervals$TenThousandStart <= PdiffCutoff$POS[i]) &(Intervals$TenThousandEnd>=PdiffCutoff$POS[i]))
  if(length(SNP.Overlap>=1)){
    PdiffCutoff$GeneOverlap[i] = Intervals$Gene[SNP.Overlap]
  } else {PdiffCutoff$GeneOverlap[i] = 0}
}

#this gives us the number of unique genes our SNPs of interest are in
signumber <- length(unique(PdiffCutoff$GeneOverlap)) - 1
NumberGenes[4] <- signumber

#permutation 10,000 times
GeneTable <- rep(NA, 10000) #list that we will put our number of genes in

for(i in 1:10000){
  table1 <- EntireSample[sample(nrow(EntireSample), nrow(PdiffCutoff)), ] #takes a random sample from the set of SNPs that passes the genome wide cutoff with the number of SNPs that is significant at out threshold
  GeneTable[i] <- (length(unique(table1$GeneOverlap)) - 1) #you have to subtract one because of the zero that's hanging around
}

write.table(GeneTable, "HipCircPdiff1Permutation.tsv")

FinishedTable <- read.table("HipCircPdiff1Permutation.tsv")

CalcPercent <- FinishedTable %>% filter(x>signumber)

PercentTable[4] <- (nrow(CalcPercent)/10000)

write.table(PercentTable, "HipCircPermutationPercentile.tsv")
write.table(NumberGenes, "HipCircNumberGenesPdiff.tsv")
