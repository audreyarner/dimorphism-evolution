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
Phenotype <- read.table("21002_irnt.FemalewithRSID.tsv", header = TRUE)
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
write.csv(Final.Phenotype.tSDS, "FemaleBodyMasswithtSDS.csv", row.names = FALSE)

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
write.csv(FinalBinSDSWOOD, "FemaleBodyMassSDSBins.csv", row.names = FALSE)

#====================================================================
#--The following repeats the process for males
#====================================================================

Phenotype <- read.table("21002_irnt.MalewithRSID.tsv", header = TRUE)
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
write.csv(Final.Phenotype.tSDS, "MaleBodyMasswithtSDS.csv", row.names = FALSE)

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
write.csv(FinalBinSDSWOOD, "MaleBodyMassSDSBins.csv", row.names = FALSE)

FemalePheno <- read.csv("FemaleBodyMasswithtSDS.csv", header = TRUE)
MalePheno <- read.csv("MaleBodyMasswithtSDS.csv", header = TRUE)

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
write.table(Phenotype1Clean, "MaleBodyMass.tSDS.FDR.tsv", row.names = FALSE)
write.table(Phenotype2Clean, "FemaleBodyMass.tSDS.FDR.tsv", row.names = FALSE)

MaleSmall <- Phenotype1Clean
FemaleSmall <- Phenotype2Clean

#selects only the columns needed to create the Manhattan plot
MaleSmall <- MaleSmall %>% select(rsid, CHR, POS, pval)
FemaleSmall <- FemaleSmall %>% select(rsid, CHR, POS, pval)

#filters out high pvalues so they 
MaleSmall <- MaleSmall %>% filter(-log10(pval) >1 )
FemaleSmall <- FemaleSmall %>% filter(-log10(pval) >1 )

write.table(MaleSmall, "BodyMassMaleManhattan.tsv", row.names = FALSE)
write.table(FemaleSmall, "BodyMassFemaleManhattan.tsv", row.names = FALSE)

#---------------------------------------------------------------------
#this part of the script is used to get the Spearman correlation coefficient (rho)
#---------------------------------------------------------------------

BodyMassMale2 <- read.table("MaleBodyMass.tSDS.FDR.tsv", header = TRUE)
BodyMassFemale2 <- read.table ("FemaleBodyMass.tSDS.FDR.tsv", header = TRUE)

BodyMassMale2 <- BodyMassMale2 %>% select(rsid, CHR, POS, AA, DA, variant, n_complete_samples, beta, se, tstat, pval, tSDS)
BodyMassFemale2 <- BodyMassFemale2 %>% select(rsid, CHR, POS, AA, DA, variant, n_complete_samples, beta, se, tstat, pval, tSDS)

#changes the name of the column from Beta to differentiate for later
names(BodyMassMale2)[8] <- "BETAM"
names(BodyMassFemale2)[8] <- "BETAW"

#makes a table with only the marker name and beta values
BETABodyMassMale <- BodyMassMale2 %>% select(rsid, BETAM) #need to put in the MarkerName
BETABodyMassFemale <- BodyMassFemale2 %>% select(rsid, BETAW) #need to put in the MarkerName

#joins the tables so that there is a column of markernames with the beta values for
##Male and Female right next to each other
BETABodyMass <- merge(BETABodyMassMale, BETABodyMassFemale, by = "rsid") 

#Spearman rank correlation test
BodyMassCorr1 <- cor.test(~ BETAM + BETAW, data = BETABodyMass, method = 'spearman', exact = FALSE)
BodyMassCorr1
#grabs the rho value, needed later for the pdiff equation 
BodyMassCorr1$estimate

#adds a column to both the Male and Female dataframes of the rho value
BodyMassMale2=mutate(BodyMassMale2, rho=BodyMassCorr1$estimate)
BodyMassFemale2=mutate(BodyMassFemale2, rho=BodyMassCorr1$estimate)

#makes a table for these that can be uploaded
write.table(BodyMassMale2, "BodyMassMaleLessCol.tsv", row.names= FALSE, quote = FALSE)
write.table(BodyMassFemale2, "BodyMassFemaleLessCol.tsv", row.names= FALSE, quote = FALSE)

#-----------------------------------------
#this part of the script is used to get the pdiff value
#-----------------------------------------

BodyMassMale2 <- read.table("BodyMassMaleLessCol.tsv", header = TRUE)
BodyMassFemale2 <- read.table("BodyMassFemaleLessCol.tsv", header = TRUE)

tvalueF=function(bMale,bFemale,SEMale,SEFemale,rho){
  return ((bMale - bFemale) / (sqrt((SEMale)^2 + (SEFemale)^2 - (2*rho*SEMale*SEFemale))))
}

#creates a function to give the pvalue from the tvalue
#make sure to include not only the values needed to calculate pvalue, but also those for tvalue
pvalueF=function(bMale,bFemale,SEMale,SEFemale,rho,DegreesFreedom){
  return(2*pt(-abs(tvalueF(bMale,bFemale,SEMale,SEFemale,rho)), (DegreesFreedom-1)))
}

names(BodyMassMale2)[7] <- "Male.N"
names(BodyMassMale2)[8] <- "Male.beta"
names(BodyMassMale2)[9] <- "Male.se"
names(BodyMassMale2)[10] <- "Male.tstat"
names(BodyMassMale2)[11] <- "Male.pval"
names(BodyMassMale2)[12] <- "Male.tSDS"

names(BodyMassFemale2)[7] <- "Female.N"
names(BodyMassFemale2)[8] <- "Female.beta"
names(BodyMassFemale2)[9] <- "Female.se"
names(BodyMassFemale2)[10] <- "Female.tstat"
names(BodyMassFemale2)[11] <- "Female.pval"
names(BodyMassFemale2)[12] <- "Female.tSDS"

#concatenates the Male and Female values into the same datatable
UKBBBodyMass <- merge(BodyMassMale2, BodyMassFemale2)

#cleaning data
UKBBBodyMass$Male.beta<-as.numeric(as.character(UKBBBodyMass$Male.beta)) #beta
UKBBBodyMass$Female.beta<-as.numeric(as.character(UKBBBodyMass$Female.beta)) #beta
UKBBBodyMass$Male.se <-as.numeric(as.character(UKBBBodyMass$Male.se)) #standarderror
UKBBBodyMass$Female.se <-as.numeric(as.character(UKBBBodyMass$Female.se)) #standarderror
UKBBBodyMass$Male.N<-as.numeric(as.character(UKBBBodyMass$Male.N)) #number of individuals
UKBBBodyMass$Female.N<-as.numeric(as.character(UKBBBodyMass$Female.N)) #number of individuals
UKBBBodyMass$Male.tSDS<-as.numeric(as.character(UKBBBodyMass$Male.tSDS)) #tSDS value
UKBBBodyMass$Female.tSDS<-as.numeric(as.character(UKBBBodyMass$Female.tSDS)) #tSDS value
UKBBBodyMass$rho<-as.numeric(as.character(UKBBBodyMass$rho)) #rho value from previous script


#mutate a new variable "pval" which holds the result of the function pvalueF
#make sure that the arguMalets passed into pvalueF are in the right order
#pvalueF(bMale,bFemale,SEMale,SEFemale,rho,DegreesFreedom) <-- these are the arguMalets in order
#the arguMalet for DegreesFreedom is (BMIadjRandallRaw$Male.N+BMIadjRandallRaw$Female.N)
UKBBBodyMass=mutate(UKBBBodyMass,pdiff=pvalueF(UKBBBodyMass$Male.beta,UKBBBodyMass$Female.beta,UKBBBodyMass$Male.se,UKBBBodyMass$Female.se,UKBBBodyMass$rho,(UKBBBodyMass$Male.N+UKBBBodyMass$Female.N)))
UKBBBodyMass=mutate(UKBBBodyMass,tdiff=tvalueF(UKBBBodyMass$Male.beta,UKBBBodyMass$Female.beta,UKBBBodyMass$Male.se,UKBBBodyMass$Female.se,UKBBBodyMass$rho))

write.table(UKBBBodyMass, "UKBBBodyMasswithPDIFF.tsv", row.names = FALSE, quote = FALSE)

#------------------------------------------------------
#this part of the script is used to choose SNPs at various cutoffs
#------------------------------------------------------

library(dplyr) #needed for %>%

#reads in the dataframe
ConcatenatedSexBodyMass <- read.table("UKBBBodyMasswithPDIFF.tsv", header = TRUE)

#assigns the genome-wide significance value, P=5x10-8, as our cutoff line
cutoff <- 5e-8

#applies the cutoff to the Female values and makes a new dataframe
FemaleSig<- ConcatenatedSexBodyMass %>% filter(Female.pval <= cutoff)

#applies the cutoff to the Male values and makes a new dataframe
MaleSig<- ConcatenatedSexBodyMass %>% filter(Male.pval <= cutoff)

AllSig <- merge(MaleSig, FemaleSig, all=TRUE)

#calculates the FDR for the pdiff values
BodyMassPdiffFDR <- unname(unlist(AllSig[,"pdiff"]))
fdr <- p.adjust(BodyMassPdiffFDR, method = "fdr")
AllSig$pdiff.FDR <- fdr

write.table(AllSig, "BodyMasswithPDIFF.FDR.tsv", row.names = FALSE)

#applies different FDR cutoffs for pdiff and writes a table
Ex1 <- AllSig %>% filter(pdiff.FDR <= .05)
write.table(Ex1, "BodyMassPdiff1.tsv", row.names = FALSE)

Ex2 <- AllSig %>% filter(pdiff.FDR <= .01)
write.table(Ex2, "BodyMassPdiff2.tsv", row.names = FALSE)

Ex3 <- AllSig %>% filter(pdiff.FDR <= .005)
write.table(Ex3, "BodyMassPdiff3.tsv", row.names = FALSE)

Ex4 <- AllSig %>% filter(pdiff.FDR <= .001)
write.table(Ex4, "BodyMassPdiff4.tsv", row.names = FALSE)


#------------------------------------------------------
#this part of the script calculates SDS
#------------------------------------------------------

tSDSTableMale <- rep(NA, 4)
SNPsTableMale <- rep(NA, 4)

tSDSTableFemale <- rep(NA, 4)
SNPsTableFemale <- rep(NA, 4)

FDRtable <- read.table("BodyMassPdiff4.tsv", header = TRUE)

FDRtableFemale <- FDRtable %>% filter(Female.pval <= Male.pval)

tSDSTableFemale[1] <- mean(FDRtableFemale$Female.tSDS)
SNPsTableFemale[1] <- nrow(FDRtableFemale)

FDRtableMale <- FDRtable %>% filter(Male.pval < Female.pval)

tSDSTableMale[1] <- mean(FDRtableMale$Male.tSDS)
SNPsTableMale[1] <- nrow(FDRtableMale)

##########################
FDRtable <- read.table("BodyMassPdiff3.tsv", header = TRUE)

FDRtableFemale <- FDRtable %>% filter(Female.pval <= Male.pval)

tSDSTableFemale[2] <- mean(FDRtableFemale$Female.tSDS)
SNPsTableFemale[2] <- nrow(FDRtableFemale)

FDRtableMale <- FDRtable %>% filter(Male.pval < Female.pval)

tSDSTableMale[2] <- mean(FDRtableMale$Male.tSDS)
SNPsTableMale[2] <- nrow(FDRtableMale)
##########################
FDRtable <- read.table("BodyMassPdiff2.tsv", header = TRUE)

FDRtableFemale <- FDRtable %>% filter(Female.pval <= Male.pval)

tSDSTableFemale[3] <- mean(FDRtableFemale$Female.tSDS)
SNPsTableFemale[3] <- nrow(FDRtableFemale)

FDRtableMale <- FDRtable %>% filter(Male.pval < Female.pval)

tSDSTableMale[3] <- mean(FDRtableMale$Male.tSDS)
SNPsTableMale[3] <- nrow(FDRtableMale)


##########################
FDRtable <- read.table("BodyMassPdiff1.tsv", header = TRUE)

FDRtableFemale <- FDRtable %>% filter(Female.pval <= Male.pval)

tSDSTableFemale[4] <- mean(FDRtableFemale$Female.tSDS)
SNPsTableFemale[4] <- nrow(FDRtableFemale)

FDRtableMale <- FDRtable %>% filter(Male.pval < Female.pval)

tSDSTableMale[4] <- mean(FDRtableMale$Male.tSDS)
SNPsTableMale[4] <- nrow(FDRtableMale)

write.table(tSDSTableMale, "BodyMasstSDSTableMale.txt")
write.table(SNPsTableMale, "BodyMassSNPsTableMale.txt")

write.table(tSDSTableFemale, "BodyMasstSDSTableFemale.txt")
write.table(SNPsTableFemale, "BodyMassSNPsTableFemale.txt")

