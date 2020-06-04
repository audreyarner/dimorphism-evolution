#!/usr/bin/env Rscript

#=========================================================================
#--Identifies the tSDS score for each SNP
#=========================================================================
#change lines 31, 75, 93, 95, 98, 100 with the phenotype

#required for %>% and the rbind functions
library(dplyr)

#function to return an SDS value of the opposite sign if the derived allele is not the same as the tested allele in the GWAS
#it is also used to return an SDS value of the opposite sign if the beta value is negative
#--this is because a positive tSDS score means there is an increased frequency in the TRAIT INCREASING allele
traitSDS = function(SDS){
	return(SDS*(-1))
}

#the chromosome, position, and two alleles are all concatenated together in the original file and separated by ":"
#--this function returns just the final value, which is the effect allele 
GiveAllele = function(information){
	return(gsub(".*:", "", information))
}

#The UK10K data set has all the SDS values as computed in Field et al. from the Pritchard lab
#this file can be retrieved from https://datadryad.org/resource/doi:10.5061/dryad.kd58f/1
#tSDS file should have 4451436 rows
UK10KDataSet <- read.table("/storage/home/ama6560/work/SDS_UK10K_n3195_release_Sep_19_2016.txt", header = TRUE)
#the phenotype data you are interested in. Should have ~13 million rows
Phenotype <- read.table("49_irnt.FemalewithRSID.tsv", header = TRUE)
nrow(Phenotype)

#Extracts the effect allele from the variant string
Phenotype.Clean <- mutate(Phenotype, EffectAllele=GiveAllele(Phenotype$variant))

#makes the data easier to combine by making the common column name the same
names(UK10KDataSet)[3] <- "rsid"

#merges the phenotype and SDS values vertically based on the MarkerName
#--this gets rid of all SNPs that don't have a corresponding tSDS score
PhenotypeandSDS <- merge(UK10KDataSet, Phenotype.Clean, by = "rsid")

#ensures the derived allele from the SDS analysis and Tested_Allele inputs are both characters 
PhenotypeandSDS$DA<-as.character(PhenotypeandSDS$DA)
PhenotypeandSDS$Tested_Allele<-as.character(PhenotypeandSDS$EffectAllele)

#filters based on if the SDS sign corresponds to the correct allele
CorrectAllele <- PhenotypeandSDS %>% filter(DA == Tested_Allele)
IncorrectAllele <- PhenotypeandSDS %>% filter(DA != Tested_Allele)

#makes a new column with the correct SDS score (tSDS)
#--for the dataframe where the tested allele was not the derived allele, this flips the sign.
ChangedSDS = mutate(IncorrectAllele, FIRSTtSDS = traitSDS(IncorrectAllele$SDS))
CorrectAllele$FIRSTtSDS = CorrectAllele$SDS

#re-merges the previous two datatables horizontally
PhenotypeandtSDS <- rbind(CorrectAllele, ChangedSDS)

#renames dataframe so that it can be checked on if running locally
PhenotypeandtSDS2 <- PhenotypeandtSDS

#filters based on whether the effect is positive or negative
CorrectEffect <- PhenotypeandtSDS2 %>% filter(beta > 0)
IncorrectEffect <- PhenotypeandtSDS2 %>% filter(beta <= 0)

#flips the sign if the beta value is negative so that only trait-increasing alleles are positive
IncorrectEffect2 <- mutate(IncorrectEffect, tSDS = traitSDS(IncorrectEffect$FIRSTtSDS))
CorrectEffect$tSDS = CorrectEffect$FIRSTtSDS

#re-merges the previous two datatables horizontally
Final.Phenotype.tSDS <- rbind(CorrectEffect, IncorrectEffect2)

#a complete table with every SNP and corresponding tSDS score
write.csv(Final.Phenotype.tSDS, "FemaleHipCircwithtSDS.csv", row.names = FALSE)

#====================================================================
#--The following creates a graph as in Field et al., 2016 Figure 4A
#====================================================================

#orders the table by increasing P value
OrderedbyP <- Final.Phenotype.tSDS[order(Final.Phenotype.tSDS$pval),]

#creates a new matrix with bins of 1000 SNPs, which has a mean taken
#may not be perfectly divisible by 1000
MeantSDS <- colMeans(matrix(OrderedbyP$tSDS, nrow=1000))
MeanP <- colMeans(matrix(OrderedbyP$pval, nrow=1000))

#merges the two binned means vertically
FinalBinSDSWOOD <- cbind(MeantSDS, MeanP)

#writes a .csv file for use in the genome-wide plot
write.csv(FinalBinSDSWOOD, "FemaleHipCircSDSBins.csv", row.names = FALSE)

#====================================================================
#--The following repeats the process for males
#====================================================================

Phenotype <- read.table("49_irnt.MalewithRSID.tsv", header = TRUE)
nrow(Phenotype)

#Extracts the effect allele from the variant string
Phenotype.Clean <- mutate(Phenotype, EffectAllele=GiveAllele(Phenotype$variant))

#makes the data easier to combine by making the common column name the same
names(UK10KDataSet)[3] <- "rsid"

#merges the phenotype and SDS values vertically based on the MarkerName
#--this gets rid of all SNPs that don't have a corresponding tSDS score
PhenotypeandSDS <- merge(UK10KDataSet, Phenotype.Clean, by = "rsid")

#ensures the derived allele from the SDS analysis and Tested_Allele inputs are both characters 
PhenotypeandSDS$DA<-as.character(PhenotypeandSDS$DA)
PhenotypeandSDS$Tested_Allele<-as.character(PhenotypeandSDS$EffectAllele)

#filters based on if the SDS sign corresponds to the correct allele
CorrectAllele <- PhenotypeandSDS %>% filter(DA == Tested_Allele)
IncorrectAllele <- PhenotypeandSDS %>% filter(DA != Tested_Allele)

#makes a new column with the correct SDS score (tSDS)
#--for the dataframe where the tested allele was not the derived allele, this flips the sign.
ChangedSDS = mutate(IncorrectAllele, FIRSTtSDS = traitSDS(IncorrectAllele$SDS))
CorrectAllele$FIRSTtSDS = CorrectAllele$SDS

#re-merges the previous two datatables horizontally
PhenotypeandtSDS <- rbind(CorrectAllele, ChangedSDS)

#renames dataframe so that it can be checked on if running locally
PhenotypeandtSDS2 <- PhenotypeandtSDS

#filters based on whether the effect is positive or negative
CorrectEffect <- PhenotypeandtSDS2 %>% filter(beta > 0)
IncorrectEffect <- PhenotypeandtSDS2 %>% filter(beta <= 0)

#flips the sign if the beta value is negative so that only trait-increasing alleles are positive
IncorrectEffect2 <- mutate(IncorrectEffect, tSDS = traitSDS(IncorrectEffect$FIRSTtSDS))
CorrectEffect$tSDS = CorrectEffect$FIRSTtSDS

#re-merges the previous two datatables horizontally
Final.Phenotype.tSDS <- rbind(CorrectEffect, IncorrectEffect2)

#a complete table with every SNP and corresponding tSDS score
write.csv(Final.Phenotype.tSDS, "MaleHipCircwithtSDS.csv", row.names = FALSE)

#====================================================================
#--The following creates a graph as in Field et al., 2016 Figure 4A
#====================================================================

#orders the table by increasing P value
OrderedbyP <- Final.Phenotype.tSDS[order(Final.Phenotype.tSDS$pval),]

#creates a new matrix with bins of 1000 SNPs, which has a mean taken
#may not be perfectly divisible by 1000
MeantSDS <- colMeans(matrix(OrderedbyP$tSDS, nrow=1000))
MeanP <- colMeans(matrix(OrderedbyP$pval, nrow=1000))

#merges the two binned means vertically
FinalBinSDSWOOD <- cbind(MeantSDS, MeanP)

#writes a .csv file for use in the genome-wide plot
write.csv(FinalBinSDSWOOD, "MaleHipCircSDSBins.csv", row.names = FALSE)

FemalePheno <- read.csv("FemaleHipCircwithtSDS.csv", header = TRUE)
MalePheno <- read.csv("MaleHipCircwithtSDS.csv", header = TRUE)

#takes out the remaining SNPs that are low confidence variants or have a minor allele frequency of less than 0.05
Phenotype1Clean <- MalePheno %>% filter(low_confidence_variant == "false")
nrow(Phenotype1Clean)
Phenotype1Clean <- Phenotype1Clean %>% filter(minor_AF >= .05)
nrow(Phenotype1Clean)
Phenotype2Clean <- FemalePheno %>% filter(low_confidence_variant == "false")
nrow(Phenotype2Clean)
Phenotype2Clean <- Phenotype2Clean %>% filter(minor_AF >= .05)
nrow(Phenotype2Clean)

#total table with tSDS, FDR, and Bonferroni
write.table(Phenotype1Clean, "MaleHipCirc.tSDS.FDR.tsv", row.names = FALSE)
write.table(Phenotype2Clean, "FemaleHipCirc.tSDS.FDR.tsv", row.names = FALSE)

MaleSmall <- Phenotype1Clean
FemaleSmall <- Phenotype2Clean

#selects only the columns needed to create the Manhattan plot
MaleSmall <- MaleSmall %>% select(rsid, CHR, POS, pval)
FemaleSmall <- FemaleSmall %>% select(rsid, CHR, POS, pval)

#filters out high pvalues so they 
MaleSmall <- MaleSmall %>% filter(-log10(pval) >1 )
FemaleSmall <- FemaleSmall %>% filter(-log10(pval) >1 )

write.table(MaleSmall, "HipCircMaleManhattan.tsv", row.names = FALSE)
write.table(FemaleSmall, "HipCircFemaleManhattan.tsv", row.names = FALSE)

#---------------------------------------------------------------------
#this part of the script is used to get the Spearman correlation coefficient (rho)
#---------------------------------------------------------------------

HipCircMale2 <- read.table("MaleHipCirc.tSDS.FDR.tsv", header = TRUE)
HipCircFemale2 <- read.table ("FemaleHipCirc.tSDS.FDR.tsv", header = TRUE)

HipCircMale2 <- HipCircMale2 %>% select(rsid, CHR, POS, AA, DA, variant, n_complete_samples, beta, se, tstat, pval, tSDS)
HipCircFemale2 <- HipCircFemale2 %>% select(rsid, CHR, POS, AA, DA, variant, n_complete_samples, beta, se, tstat, pval, tSDS)

#changes the name of the column from Beta to differentiate for later
names(HipCircMale2)[8] <- "BETAM"
names(HipCircFemale2)[8] <- "BETAW"

#makes a table with only the marker name and beta values
BETAHipCircMale <- HipCircMale2 %>% select(rsid, BETAM) #need to put in the MarkerName
BETAHipCircFemale <- HipCircFemale2 %>% select(rsid, BETAW) #need to put in the MarkerName

#joins the tables so that there is a column of markernames with the beta values for
##Male and Female right next to each other
BETAHipCirc <- merge(BETAHipCircMale, BETAHipCircFemale, by = "rsid") 

#Spearman rank correlation test
HipCircCorr1 <- cor.test(~ BETAM + BETAW, data = BETAHipCirc, method = 'spearman', exact = FALSE)
HipCircCorr1
#grabs the rho value, needed later for the pdiff equation 
HipCircCorr1$estimate

#adds a column to both the Male and Female dataframes of the rho value
HipCircMale2=mutate(HipCircMale2, rho=HipCircCorr1$estimate)
HipCircFemale2=mutate(HipCircFemale2, rho=HipCircCorr1$estimate)

#makes a table for these that can be uploaded
write.table(HipCircMale2, "HipCircMaleLessCol.tsv", row.names= FALSE, quote = FALSE)
write.table(HipCircFemale2, "HipCircFemaleLessCol.tsv", row.names= FALSE, quote = FALSE)

#-----------------------------------------
#this part of the script is used to get the pdiff value
#-----------------------------------------

HipCircMale2 <- read.table("HipCircMaleLessCol.tsv", header = TRUE)
HipCircFemale2 <- read.table("HipCircFemaleLessCol.tsv", header = TRUE)

tvalueF=function(bMale,bFemale,SEMale,SEFemale,rho){
  return ((bMale - bFemale) / (sqrt((SEMale)^2 + (SEFemale)^2 - (2*rho*SEMale*SEFemale))))
}

#creates a function to give the pvalue from the tvalue
#make sure to include not only the values needed to calculate pvalue, but also those for tvalue
pvalueF=function(bMale,bFemale,SEMale,SEFemale,rho,DegreesFreedom){
  return(2*pt(-abs(tvalueF(bMale,bFemale,SEMale,SEFemale,rho)), (DegreesFreedom-1)))
}

names(HipCircMale2)[7] <- "Male.N"
names(HipCircMale2)[8] <- "Male.beta"
names(HipCircMale2)[9] <- "Male.se"
names(HipCircMale2)[10] <- "Male.tstat"
names(HipCircMale2)[11] <- "Male.pval"
names(HipCircMale2)[12] <- "Male.tSDS"

names(HipCircFemale2)[7] <- "Female.N"
names(HipCircFemale2)[8] <- "Female.beta"
names(HipCircFemale2)[9] <- "Female.se"
names(HipCircFemale2)[10] <- "Female.tstat"
names(HipCircFemale2)[11] <- "Female.pval"
names(HipCircFemale2)[12] <- "Female.tSDS"

#concatenates the Male and Female values into the same datatable
UKBBHipCirc <- merge(HipCircMale2, HipCircFemale2)

#cleaning data
UKBBHipCirc$Male.beta<-as.numeric(as.character(UKBBHipCirc$Male.beta)) #beta
UKBBHipCirc$Female.beta<-as.numeric(as.character(UKBBHipCirc$Female.beta)) #beta
UKBBHipCirc$Male.se <-as.numeric(as.character(UKBBHipCirc$Male.se)) #standarderror
UKBBHipCirc$Female.se <-as.numeric(as.character(UKBBHipCirc$Female.se)) #standarderror
UKBBHipCirc$Male.N<-as.numeric(as.character(UKBBHipCirc$Male.N)) #number of individuals
UKBBHipCirc$Female.N<-as.numeric(as.character(UKBBHipCirc$Female.N)) #number of individuals
UKBBHipCirc$Male.tSDS<-as.numeric(as.character(UKBBHipCirc$Male.tSDS)) #tSDS value
UKBBHipCirc$Female.tSDS<-as.numeric(as.character(UKBBHipCirc$Female.tSDS)) #tSDS value
UKBBHipCirc$rho<-as.numeric(as.character(UKBBHipCirc$rho)) #rho value from previous script


#mutate a new variable "pval" which holds the result of the function pvalueF
#make sure that the arguMalets passed into pvalueF are in the right order
#pvalueF(bMale,bFemale,SEMale,SEFemale,rho,DegreesFreedom) <-- these are the arguMalets in order
#the arguMalet for DegreesFreedom is (BMIadjRandallRaw$Male.N+BMIadjRandallRaw$Female.N)
UKBBHipCirc=mutate(UKBBHipCirc,pdiff=pvalueF(UKBBHipCirc$Male.beta,UKBBHipCirc$Female.beta,UKBBHipCirc$Male.se,UKBBHipCirc$Female.se,UKBBHipCirc$rho,(UKBBHipCirc$Male.N+UKBBHipCirc$Female.N)))
UKBBHipCirc=mutate(UKBBHipCirc,tdiff=tvalueF(UKBBHipCirc$Male.beta,UKBBHipCirc$Female.beta,UKBBHipCirc$Male.se,UKBBHipCirc$Female.se,UKBBHipCirc$rho))

write.table(UKBBHipCirc, "UKBBHipCircwithPDIFF.tsv", row.names = FALSE, quote = FALSE)

#------------------------------------------------------
#this part of the script is used to choose SNPs at various cutoffs
#------------------------------------------------------

library(dplyr) #needed for %>%

#reads in the dataframe
ConcatenatedSexHipCirc <- read.table("UKBBHipCircwithPDIFF.tsv", header = TRUE)

#assigns the genome-wide significance value, P=5x10-8, as our cutoff line
cutoff <- 5e-8

#applies the cutoff to the Female values and makes a new dataframe
FemaleSig<- ConcatenatedSexHipCirc %>% filter(Female.pval <= cutoff)

#applies the cutoff to the Male values and makes a new dataframe
MaleSig<- ConcatenatedSexHipCirc %>% filter(Male.pval <= cutoff)

AllSig <- merge(MaleSig, FemaleSig, all=TRUE)

#calculates the FDR for the pdiff values
HipCircPdiffFDR <- unname(unlist(AllSig[,"pdiff"]))
fdr <- p.adjust(HipCircPdiffFDR, method = "fdr")
AllSig$pdiff.FDR <- fdr

write.table(AllSig, "HipCircwithPDIFF.FDR.tsv", row.names = FALSE)

#applies different FDR cutoffs for pdiff and writes a table
Ex1 <- AllSig %>% filter(pdiff.FDR <= .05)
write.table(Ex1, "HipCircPdiff1.tsv", row.names = FALSE)

Ex2 <- AllSig %>% filter(pdiff.FDR <= .01)
write.table(Ex2, "HipCircPdiff2.tsv", row.names = FALSE)

Ex3 <- AllSig %>% filter(pdiff.FDR <= .005)
write.table(Ex3, "HipCircPdiff3.tsv", row.names = FALSE)

Ex4 <- AllSig %>% filter(pdiff.FDR <= .001)
write.table(Ex4, "HipCircPdiff4.tsv", row.names = FALSE)


#------------------------------------------------------
#this part of the script calculates SDS
#------------------------------------------------------

tSDSTableMale <- rep(NA, 4)
SNPsTableMale <- rep(NA, 4)

tSDSTableFemale <- rep(NA, 4)
SNPsTableFemale <- rep(NA, 4)

FDRtable <- read.table("HipCircPdiff4.tsv", header = TRUE)

FDRtableFemale <- FDRtable %>% filter(Female.pval <= Male.pval)

tSDSTableFemale[1] <- mean(FDRtableFemale$Female.tSDS)
SNPsTableFemale[1] <- nrow(FDRtableFemale)

FDRtableMale <- FDRtable %>% filter(Male.pval < Female.pval)

tSDSTableMale[1] <- mean(FDRtableMale$Male.tSDS)
SNPsTableMale[1] <- nrow(FDRtableMale)

##########################
FDRtable <- read.table("HipCircPdiff3.tsv", header = TRUE)

FDRtableFemale <- FDRtable %>% filter(Female.pval <= Male.pval)

tSDSTableFemale[2] <- mean(FDRtableFemale$Female.tSDS)
SNPsTableFemale[2] <- nrow(FDRtableFemale)

FDRtableMale <- FDRtable %>% filter(Male.pval < Female.pval)

tSDSTableMale[2] <- mean(FDRtableMale$Male.tSDS)
SNPsTableMale[2] <- nrow(FDRtableMale)
##########################
FDRtable <- read.table("HipCircPdiff2.tsv", header = TRUE)

FDRtableFemale <- FDRtable %>% filter(Female.pval <= Male.pval)

tSDSTableFemale[3] <- mean(FDRtableFemale$Female.tSDS)
SNPsTableFemale[3] <- nrow(FDRtableFemale)

FDRtableMale <- FDRtable %>% filter(Male.pval < Female.pval)

tSDSTableMale[3] <- mean(FDRtableMale$Male.tSDS)
SNPsTableMale[3] <- nrow(FDRtableMale)


##########################
FDRtable <- read.table("HipCircPdiff1.tsv", header = TRUE)

FDRtableFemale <- FDRtable %>% filter(Female.pval <= Male.pval)

tSDSTableFemale[4] <- mean(FDRtableFemale$Female.tSDS)
SNPsTableFemale[4] <- nrow(FDRtableFemale)

FDRtableMale <- FDRtable %>% filter(Male.pval < Female.pval)

tSDSTableMale[4] <- mean(FDRtableMale$Male.tSDS)
SNPsTableMale[4] <- nrow(FDRtableMale)

write.table(tSDSTableMale, "HipCirctSDSTableMale.txt")
write.table(SNPsTableMale, "HipCircSNPsTableMale.txt")

write.table(tSDSTableFemale, "HipCirctSDSTableFemale.txt")
write.table(SNPsTableFemale, "HipCircSNPsTableFemale.txt")

