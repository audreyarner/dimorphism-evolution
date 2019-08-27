#!/usr/bin/env Rscript

#-----------------------------------------------------------------------------------------
#this script is used for our permutation analysis
#-----------------------------------------------------------------------------------------

#packages needed for this script
library(dplyr)
library(ggplot2)

#pulls the file from each phenotype that gives the SNPs crossing the genome-wide threshold
files = list.files(path="/storage/home/ama6560/scratch/Enrichment", pattern = "withPDIFF.FDR.tsv")

#list of genes and their chromosomal position (hg19) from GO: 0007548
Intervals <- read.table("GOwGeneNamesNoRepeats6-13.txt", header = TRUE)

#makes sure formatting is all correct
Intervals$TenThousandStart <- as.numeric(Intervals$TenThousandStart)
Intervals$TenThousandEnd <- as.numeric(Intervals$TenThousandEnd)
Intervals$CHR<-as.numeric(Intervals$CHR)
Intervals$Gene<-as.character(Intervals$Gene)

#the input to the function is a phenotype file with the list of SNPs that cross the genome-wide threshold
genEnrichment = function(phenotype){
  #read the file
  SNPs <- read.table(phenotype, header = TRUE)
  #Grabs the name of the phenotype from the file
  PhenotypeName <- gsub("withPDIFF.FDR.tsv", "", phenotype)

  SNPs$POS <- as.numeric(SNPs$POS)
  SNPs$CHR<-as.numeric(SNPs$CHR)

  #makes a new column in our dataframe, which we will put if the SNP is in a gene, and if so, which one
  SNPs[,"GeneOverlap"] <- NA

  for (i in 1:nrow(SNPs)){
    #which() will give the position of the gene that overlaps, which we can then use to find the gene name
    SNP.Overlap = which((Intervals$CHR == SNPs$CHR[i]) & (Intervals$TenThousandStart <= SNPs$POS[i]) &(Intervals$TenThousandEnd>=SNPs$POS[i]))
    if(length(SNP.Overlap>=1)){
      SNPs$GeneOverlap[i] = Intervals$Gene[SNP.Overlap]
    } else {SNPs$GeneOverlap[i] = 0}
  }
  
  #initializes tables where we can put our permutation p-value and the number of genes with significant SNPs
  PercentTable <- rep(NA, 4)
  NumberGenes <- rep(NA, 4)

  #orders the height SNPs first by chromosome, then position
  SNPs2 <- SNPs[order(SNPs[,2], SNPs[,3]),]

  #makes an additional row that tells what each SNP's chromosomal ordering is
  SNPs2 <- SNPs2 %>% mutate(Index = row_number())

  #selects only the rsid and index value (to then add to the sexually dimorphic SNPs)
  SNPs3 <- SNPs2

  ###############################
  #filters out only the SNPs that have a pdiff FDR of .001
  SexDiSNPs <- SNPs3 %>% filter(pdiff.FDR <= .001)

  if (nrow(SexDiSNPs)<1){ #the phenotype may not have SNPs with the pdiff value, which this if statement checks for 
    print(paste(PhenotypeName, "no SNPs Pdiff4"))
  } else { #if there are SNPs, we can continue with the analysis

    #Gives us the number of unique genes our SNPs of interest are in
    signumber <- length(unique(SexDiSNPs$GeneOverlap)) - 1 #have to subtract 1 because of the NAs that are in the table
    NumberGenes[1] <- signumber

    GeneTable <- rep(NA, 10000) #list that we will put our number of genes in

    #here begins the permutation analysis
    for(j in 1:10000){
      #takes a random sample, from which we'll have our starting position
      SampleSNPs <- sample(1:nrow(SNPs3), 1)
      #changes the index of all the SNPs so that the first will start with the index + 1
      SexDiSNPs3=mutate(SexDiSNPs, Index = SampleSNPs + SexDiSNPs$Index) 
      #loops through each row to see if it passes the max index value from the entire height dataset
        for (i in 1:nrow(SexDiSNPs3)){
          if(SexDiSNPs3$Index[i] > max(SNPs3$Index)){ #changes the index value to loop back to the front of the dataset if it is past the max index value 
            SexDiSNPs3$Index[i] = (SexDiSNPs3$Index[i]-max(SNPs3$Index))
          }
        }
        #selects only the new index values that started with the random number
        Indexes <- SexDiSNPs3 %>% select(Index)
        PermDataset <- merge(Indexes, SNPs2)
        #Gets the number of unique gene names
        GeneTable[j] <- (length(unique(PermDataset$GeneOverlap)) - 1) #you have to subtract one because of the zero that's hanging around
      }
    write.table(GeneTable, file = paste(PhenotypeName, "Pdiff4Permutation.tsv", sep=""))

    df <- as.data.frame(GeneTable)

    CalcPercent <- df %>% filter(GeneTable>=signumber)
    PercentTable[1] <-(nrow(CalcPercent)/10000)
  }
  ###############################
  SexDiSNPs <- SNPs3 %>% filter(pdiff.FDR <= .005)
  if (nrow(SexDiSNPs)<1){
    print(paste(PhenotypeName, "no SNPs Pdiff3"))
  } else {

    #this gives us the number of unique genes our SNPs of interest are in
    signumber <- length(unique(SexDiSNPs$GeneOverlap)) - 1
    NumberGenes[2] <- signumber

    GeneTable <- rep(NA, 10000) #list that we will put our number of genes in

    #here begins the permutation analysis
    for(j in 1:10000){
      #takes a random sample, from which we'll have our starting position
      SampleSNPs <- sample(1:nrow(SNPs3), 1)
      #changes the index of all the SNPs so that the first will start with the index + 1
      SexDiSNPs3=mutate(SexDiSNPs, Index = SampleSNPs + SexDiSNPs$Index) 
      #loops through each row to see if it passes the max index value from the entire height dataset
        for (i in 1:nrow(SexDiSNPs3)){
          if(SexDiSNPs3$Index[i] > max(SNPs3$Index)){
            SexDiSNPs3$Index[i] = (SexDiSNPs3$Index[i]-max(SNPs3$Index))
          }
        }
        #selects only the new index values that started with the random number
        Indexes <- SexDiSNPs3 %>% select(Index)
        PermDataset <- merge(Indexes, SNPs2)
        #Gets the number of unique gene names
        GeneTable[j] <- (length(unique(PermDataset$GeneOverlap)) - 1) #you have to subtract one because of the zero that's hanging around
      }
    write.table(GeneTable, file = paste(PhenotypeName, "Pdiff3Permutation.tsv", sep=""))

    df <- as.data.frame(GeneTable)

    CalcPercent <- df %>% filter(GeneTable>=signumber)
    PercentTable[2] <-(nrow(CalcPercent)/10000)
  }
    ###############################
  SexDiSNPs <- SNPs3 %>% filter(pdiff.FDR <= .01)
  if (nrow(SexDiSNPs)<1){
    print(paste(PhenotypeName, "no SNPs Pdiff2"))
  } else {

    #this gives us the number of unique genes our SNPs of interest are in
    signumber <- length(unique(SexDiSNPs$GeneOverlap)) - 1
    NumberGenes[3] <- signumber

    GeneTable <- rep(NA, 10000) #list that we will put our number of genes in

    #here begins the permutation analysis
    for(j in 1:10000){
      #takes a random sample, from which we'll have our starting position
      SampleSNPs <- sample(1:nrow(SNPs3), 1)
      #changes the index of all the SNPs so that the first will start with the index + 1
      SexDiSNPs3=mutate(SexDiSNPs, Index = SampleSNPs + SexDiSNPs$Index) 
      #loops through each row to see if it passes the max index value from the entire height dataset
        for (i in 1:nrow(SexDiSNPs3)){
          if(SexDiSNPs3$Index[i] > max(SNPs3$Index)){
            SexDiSNPs3$Index[i] = (SexDiSNPs3$Index[i]-max(SNPs3$Index))
          }
        }
        #selects only the new index values that started with the random number
        Indexes <- SexDiSNPs3 %>% select(Index)
        PermDataset <- merge(Indexes, SNPs2)
        #Gets the number of unique gene names
        GeneTable[j] <- (length(unique(PermDataset$GeneOverlap)) - 1) #you have to subtract one because of the zero that's hanging around
      }
    write.table(GeneTable, file = paste(PhenotypeName, "Pdiff2Permutation.tsv", sep=""))

    df <- as.data.frame(GeneTable)

    CalcPercent <- df %>% filter(GeneTable>=signumber)
    PercentTable[3] <-(nrow(CalcPercent)/10000)
  }
    ###############################
  SexDiSNPs <- SNPs3 %>% filter(pdiff.FDR <= .05)
  if (nrow(SexDiSNPs)<1){
    print(paste(PhenotypeName, "no SNPs Pdiff1"))
  } else {
    #this gives us the number of unique genes our SNPs of interest are in
    signumber <- length(unique(SexDiSNPs$GeneOverlap)) - 1
    NumberGenes[4] <- signumber

    GeneTable <- rep(NA, 10000) #list that we will put our number of genes in

    #here begins the permutation analysis
    for(j in 1:10000){
      #takes a random sample, from which we'll have our starting position
      SampleSNPs <- sample(1:nrow(SNPs3), 1)
      #changes the index of all the SNPs so that the first will start with the index + 1
      SexDiSNPs3=mutate(SexDiSNPs, Index = SampleSNPs + SexDiSNPs$Index) 
      #loops through each row to see if it passes the max index value from the entire height dataset
        for (i in 1:nrow(SexDiSNPs3)){
          if(SexDiSNPs3$Index[i] > max(SNPs3$Index)){
            SexDiSNPs3$Index[i] = (SexDiSNPs3$Index[i]-max(SNPs3$Index))
          }
        }
        #selects only the new index values that started with the random number
        Indexes <- SexDiSNPs3 %>% select(Index)
        PermDataset <- merge(Indexes, SNPs2)
        #Gets the number of unique gene names
        GeneTable[j] <- (length(unique(PermDataset$GeneOverlap)) - 1) #you have to subtract one because of the zero that's hanging around
      }
    write.table(GeneTable, file = paste(PhenotypeName, "Pdiff1Permutation.tsv", sep=""))

    df <- as.data.frame(GeneTable)

    CalcPercent <- df %>% filter(GeneTable>=signumber)
    PercentTable[4] <-(nrow(CalcPercent)/10000)
  }

  write.table(NumberGenes, file = paste(PhenotypeName, "NumberGenes.tsv", sep=""))
  write.table(PercentTable, file = paste(PhenotypeName, "PermutationPercentile.tsv", sep=""))
}


for (i in 1:length(files)){
  NumRows <- read.table(files[i], header = TRUE)
  if(nrow(NumRows)<1){ #makes sure that the phenotype has significant SNPs
    print(paste(files[i], "No significant SNPs"))
  } else {genEnrichment(files[i])} #if it does have SNPs, apply the function to it
}