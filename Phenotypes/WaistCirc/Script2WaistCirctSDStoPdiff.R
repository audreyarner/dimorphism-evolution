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
Phenotype <- read.table("48_irnt.FemalewithRSID.tsv", header = TRUE)
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
write.csv(Final.Phenotype.tSDS, "FemaleWaistCircwithtSDS.csv", row.names = FALSE)

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
write.csv(FinalBinSDSWOOD, "FemaleWaistCircSDSBins.csv", row.names = FALSE)

#====================================================================
#--The following repeats the process for males
#====================================================================

Phenotype <- read.table("48_irnt.MalewithRSID.tsv", header = TRUE)
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
write.csv(Final.Phenotype.tSDS, "MaleWaistCircwithtSDS.csv", row.names = FALSE)

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
write.csv(FinalBinSDSWOOD, "MaleWaistCircSDSBins.csv", row.names = FALSE)

FemalePheno <- read.csv("FemaleWaistCircwithtSDS.csv", header = TRUE)
MalePheno <- read.csv("MaleWaistCircwithtSDS.csv", header = TRUE)

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
write.table(Phenotype1Clean, "MaleWaistCirc.tSDS.FDR.tsv", row.names = FALSE)
write.table(Phenotype2Clean, "FemaleWaistCirc.tSDS.FDR.tsv", row.names = FALSE)

MaleSmall <- Phenotype1Clean
FemaleSmall <- Phenotype2Clean

#selects only the columns needed to create the Manhattan plot
MaleSmall <- MaleSmall %>% select(rsid, CHR, POS, pval)
FemaleSmall <- FemaleSmall %>% select(rsid, CHR, POS, pval)

#filters out high pvalues so they 
MaleSmall <- MaleSmall %>% filter(-log10(pval) >1 )
FemaleSmall <- FemaleSmall %>% filter(-log10(pval) >1 )

write.table(MaleSmall, "WaistCircMaleManhattan.tsv", row.names = FALSE)
write.table(FemaleSmall, "WaistCircFemaleManhattan.tsv", row.names = FALSE)

#---------------------------------------------------------------------
#this part of the script is used to get the Spearman correlation coefficient (rho)
#---------------------------------------------------------------------

WaistCircMale2 <- read.table("MaleWaistCirc.tSDS.FDR.tsv", header = TRUE)
WaistCircFemale2 <- read.table ("FemaleWaistCirc.tSDS.FDR.tsv", header = TRUE)

WaistCircMale2 <- WaistCircMale2 %>% select(rsid, CHR, POS, AA, DA, variant, n_complete_samples, beta, se, tstat, pval, tSDS)
WaistCircFemale2 <- WaistCircFemale2 %>% select(rsid, CHR, POS, AA, DA, variant, n_complete_samples, beta, se, tstat, pval, tSDS)

#changes the name of the column from Beta to differentiate for later
names(WaistCircMale2)[8] <- "BETAM"
names(WaistCircFemale2)[8] <- "BETAW"

#makes a table with only the marker name and beta values
BETAWaistCircMale <- WaistCircMale2 %>% select(rsid, BETAM) #need to put in the MarkerName
BETAWaistCircFemale <- WaistCircFemale2 %>% select(rsid, BETAW) #need to put in the MarkerName

#joins the tables so that there is a column of markernames with the beta values for
##Male and Female right next to each other
BETAWaistCirc <- merge(BETAWaistCircMale, BETAWaistCircFemale, by = "rsid") 

#Spearman rank correlation test
WaistCircCorr1 <- cor.test(~ BETAM + BETAW, data = BETAWaistCirc, method = 'spearman', exact = FALSE)
WaistCircCorr1
#grabs the rho value, needed later for the pdiff equation 
WaistCircCorr1$estimate

#adds a column to both the Male and Female dataframes of the rho value
WaistCircMale2=mutate(WaistCircMale2, rho=WaistCircCorr1$estimate)
WaistCircFemale2=mutate(WaistCircFemale2, rho=WaistCircCorr1$estimate)

#makes a table for these that can be uploaded
write.table(WaistCircMale2, "WaistCircMaleLessCol.tsv", row.names= FALSE, quote = FALSE)
write.table(WaistCircFemale2, "WaistCircFemaleLessCol.tsv", row.names= FALSE, quote = FALSE)

#-----------------------------------------
#this part of the script is used to get the pdiff value
#-----------------------------------------

WaistCircMale2 <- read.table("WaistCircMaleLessCol.tsv", header = TRUE)
WaistCircFemale2 <- read.table("WaistCircFemaleLessCol.tsv", header = TRUE)

tvalueF=function(bMale,bFemale,SEMale,SEFemale,rho){
  return ((bMale - bFemale) / (sqrt((SEMale)^2 + (SEFemale)^2 - (2*rho*SEMale*SEFemale))))
}

#creates a function to give the pvalue from the tvalue
#make sure to include not only the values needed to calculate pvalue, but also those for tvalue
pvalueF=function(bMale,bFemale,SEMale,SEFemale,rho,DegreesFreedom){
  return(2*pt(-abs(tvalueF(bMale,bFemale,SEMale,SEFemale,rho)), (DegreesFreedom-1)))
}

names(WaistCircMale2)[7] <- "Male.N"
names(WaistCircMale2)[8] <- "Male.beta"
names(WaistCircMale2)[9] <- "Male.se"
names(WaistCircMale2)[10] <- "Male.tstat"
names(WaistCircMale2)[11] <- "Male.pval"
names(WaistCircMale2)[12] <- "Male.tSDS"

names(WaistCircFemale2)[7] <- "Female.N"
names(WaistCircFemale2)[8] <- "Female.beta"
names(WaistCircFemale2)[9] <- "Female.se"
names(WaistCircFemale2)[10] <- "Female.tstat"
names(WaistCircFemale2)[11] <- "Female.pval"
names(WaistCircFemale2)[12] <- "Female.tSDS"

#concatenates the Male and Female values into the same datatable
UKBBWaistCirc <- merge(WaistCircMale2, WaistCircFemale2)

#cleaning data
UKBBWaistCirc$Male.beta<-as.numeric(as.character(UKBBWaistCirc$Male.beta)) #beta
UKBBWaistCirc$Female.beta<-as.numeric(as.character(UKBBWaistCirc$Female.beta)) #beta
UKBBWaistCirc$Male.se <-as.numeric(as.character(UKBBWaistCirc$Male.se)) #standarderror
UKBBWaistCirc$Female.se <-as.numeric(as.character(UKBBWaistCirc$Female.se)) #standarderror
UKBBWaistCirc$Male.N<-as.numeric(as.character(UKBBWaistCirc$Male.N)) #number of individuals
UKBBWaistCirc$Female.N<-as.numeric(as.character(UKBBWaistCirc$Female.N)) #number of individuals
UKBBWaistCirc$Male.tSDS<-as.numeric(as.character(UKBBWaistCirc$Male.tSDS)) #tSDS value
UKBBWaistCirc$Female.tSDS<-as.numeric(as.character(UKBBWaistCirc$Female.tSDS)) #tSDS value
UKBBWaistCirc$rho<-as.numeric(as.character(UKBBWaistCirc$rho)) #rho value from previous script


#mutate a new variable "pval" which holds the result of the function pvalueF
#make sure that the arguMalets passed into pvalueF are in the right order
#pvalueF(bMale,bFemale,SEMale,SEFemale,rho,DegreesFreedom) <-- these are the arguMalets in order
#the arguMalet for DegreesFreedom is (BMIadjRandallRaw$Male.N+BMIadjRandallRaw$Female.N)
UKBBWaistCirc=mutate(UKBBWaistCirc,pdiff=pvalueF(UKBBWaistCirc$Male.beta,UKBBWaistCirc$Female.beta,UKBBWaistCirc$Male.se,UKBBWaistCirc$Female.se,UKBBWaistCirc$rho,(UKBBWaistCirc$Male.N+UKBBWaistCirc$Female.N)))
UKBBWaistCirc=mutate(UKBBWaistCirc,tdiff=tvalueF(UKBBWaistCirc$Male.beta,UKBBWaistCirc$Female.beta,UKBBWaistCirc$Male.se,UKBBWaistCirc$Female.se,UKBBWaistCirc$rho))

write.table(UKBBWaistCirc, "UKBBWaistCircwithPDIFF.tsv", row.names = FALSE, quote = FALSE)

#------------------------------------------------------
#this part of the script is used to choose SNPs at various cutoffs
#------------------------------------------------------

library(dplyr) #needed for %>%

#reads in the dataframe
ConcatenatedSexWaistCirc <- read.table("UKBBWaistCircwithPDIFF.tsv", header = TRUE)

#assigns the genome-wide significance value, P=5x10-8, as our cutoff line
cutoff <- 5e-8

#applies the cutoff to the Female values and makes a new dataframe
FemaleSig<- ConcatenatedSexWaistCirc %>% filter(Female.pval <= cutoff)

#applies the cutoff to the Male values and makes a new dataframe
MaleSig<- ConcatenatedSexWaistCirc %>% filter(Male.pval <= cutoff)

AllSig <- merge(MaleSig, FemaleSig, all=TRUE)

#calculates the FDR for the pdiff values
WaistCircPdiffFDR <- unname(unlist(AllSig[,"pdiff"]))
fdr <- p.adjust(WaistCircPdiffFDR, method = "fdr")
AllSig$pdiff.FDR <- fdr

write.table(AllSig, "WaistCircwithPDIFF.FDR.tsv", row.names = FALSE)

#applies different FDR cutoffs for pdiff and writes a table
Ex1 <- AllSig %>% filter(pdiff.FDR <= .05)
write.table(Ex1, "WaistCircPdiff1.tsv", row.names = FALSE)

Ex2 <- AllSig %>% filter(pdiff.FDR <= .01)
write.table(Ex2, "WaistCircPdiff2.tsv", row.names = FALSE)

Ex3 <- AllSig %>% filter(pdiff.FDR <= .005)
write.table(Ex3, "WaistCircPdiff3.tsv", row.names = FALSE)

Ex4 <- AllSig %>% filter(pdiff.FDR <= .001)
write.table(Ex4, "WaistCircPdiff4.tsv", row.names = FALSE)


#------------------------------------------------------
#this part of the script calculates SDS
#------------------------------------------------------

tSDSTableMale <- rep(NA, 4)
SNPsTableMale <- rep(NA, 4)

tSDSTableFemale <- rep(NA, 4)
SNPsTableFemale <- rep(NA, 4)

FDRtable <- read.table("WaistCircPdiff4.tsv", header = TRUE)

FDRtableFemale <- FDRtable %>% filter(Female.pval <= Male.pval)

tSDSTableFemale[1] <- mean(FDRtableFemale$Female.tSDS)
SNPsTableFemale[1] <- nrow(FDRtableFemale)

FDRtableMale <- FDRtable %>% filter(Male.pval < Female.pval)

tSDSTableMale[1] <- mean(FDRtableMale$Male.tSDS)
SNPsTableMale[1] <- nrow(FDRtableMale)

##########################
FDRtable <- read.table("WaistCircPdiff3.tsv", header = TRUE)

FDRtableFemale <- FDRtable %>% filter(Female.pval <= Male.pval)

tSDSTableFemale[2] <- mean(FDRtableFemale$Female.tSDS)
SNPsTableFemale[2] <- nrow(FDRtableFemale)

FDRtableMale <- FDRtable %>% filter(Male.pval < Female.pval)

tSDSTableMale[2] <- mean(FDRtableMale$Male.tSDS)
SNPsTableMale[2] <- nrow(FDRtableMale)
##########################
FDRtable <- read.table("WaistCircPdiff2.tsv", header = TRUE)

FDRtableFemale <- FDRtable %>% filter(Female.pval <= Male.pval)

tSDSTableFemale[3] <- mean(FDRtableFemale$Female.tSDS)
SNPsTableFemale[3] <- nrow(FDRtableFemale)

FDRtableMale <- FDRtable %>% filter(Male.pval < Female.pval)

tSDSTableMale[3] <- mean(FDRtableMale$Male.tSDS)
SNPsTableMale[3] <- nrow(FDRtableMale)


##########################
FDRtable <- read.table("WaistCircPdiff1.tsv", header = TRUE)

FDRtableFemale <- FDRtable %>% filter(Female.pval <= Male.pval)

tSDSTableFemale[4] <- mean(FDRtableFemale$Female.tSDS)
SNPsTableFemale[4] <- nrow(FDRtableFemale)

FDRtableMale <- FDRtable %>% filter(Male.pval < Female.pval)

tSDSTableMale[4] <- mean(FDRtableMale$Male.tSDS)
SNPsTableMale[4] <- nrow(FDRtableMale)

write.table(tSDSTableMale, "WaistCirctSDSTableMale.txt")
write.table(SNPsTableMale, "WaistCircSNPsTableMale.txt")

write.table(tSDSTableFemale, "WaistCirctSDSTableFemale.txt")
write.table(SNPsTableFemale, "WaistCircSNPsTableFemale.txt")

