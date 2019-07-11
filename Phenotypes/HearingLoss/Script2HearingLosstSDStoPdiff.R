#!/usr/bin/env Rscript

#=========================================================================
#--Identifies the tSDS score for each SNP
#=========================================================================
#CHANGE PHENOTYPE INPUT FILE AND PHENOTYPE NAME IF DIFFERENT PHENOTYPE!!!

#required for %>% and the rbind functions
library(dplyr)
#required for graphs at the end
library(ggplot2)

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
UK10KDataSet <- read.table("/storage/home/ama6560/scratch/HeightDimorphism/SDS_UK10K_n3195_release_Sep_19_2016.txt", header = TRUE)
#the phenotype data you are interested in. Should have ~13 million rows
Phenotype <- read.table("H90.gwas.FemalewithRSID.tsv", header = TRUE)
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
write.csv(Final.Phenotype.tSDS, "WomenHearingLosswithtSDS.csv", row.names = FALSE)

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
write.csv(FinalBinSDSWOOD, "WomenHearingLossSDSBins.csv", row.names = FALSE)

ReplicatedFieldGraph <- read.csv("WomenHearingLossSDSBins.csv")

#makes a graph of the genome wide tSDS score of bins of 1000 SNPs, as seen in Field et al Figure 4A
ReplicatedFieldGraph %>% ggplot(aes(x=MeanP, y=MeantSDS)) + geom_point() + ggtitle("Phenotype tSDS score observed genome-wide")+ labs(x="Phenotype signifcance (rank)", y = "Phenotype-increasing tSDS (bin average)") + scale_x_reverse() + ylim(-.3,.75)

ggsave("WomenHearingLoss.tSDSGraph.png")

#====================================================================
#--The following repeats the process for males
#====================================================================

Phenotype <- read.table("H90.gwas.MalewithRSID.tsv", header = TRUE)
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
write.csv(Final.Phenotype.tSDS, "MenHearingLosswithtSDS.csv", row.names = FALSE)

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
write.csv(FinalBinSDSWOOD, "MenHearingLossSDSBins.csv", row.names = FALSE)

ReplicatedFieldGraph <- read.csv("MenHearingLossSDSBins.csv")

#makes a graph of the genome wide tSDS score of bins of 1000 SNPs, as seen in Field et al Figure 4A
ReplicatedFieldGraph %>% ggplot(aes(x=MeanP, y=MeantSDS)) + geom_point() + ggtitle("Phenotype tSDS score observed genome-wide")+ labs(x="Phenotype signifcance (rank)", y = "Phenotype-increasing tSDS (bin average)") + scale_x_reverse() + ylim(-.3,.75)

ggsave("MenHearingLoss.tSDSGraph.png")

WomenPheno <- read.csv("WomenHearingLosswithtSDS.csv", header = TRUE)
MenPheno <- read.csv("MenHearingLosswithtSDS.csv", header = TRUE)

#takes out the remaining SNPs that are low confidence variants or have a minor allele frequency of less than 0.05
Phenotype1Clean <- MenPheno %>% filter(low_confidence_variant == "false")
nrow(Phenotype1Clean)
Phenotype1Clean <- Phenotype1Clean %>% filter(minor_AF >= .05)
nrow(Phenotype1Clean)
Phenotype2Clean <- WomenPheno %>% filter(low_confidence_variant == "false")
nrow(Phenotype2Clean)
Phenotype2Clean <- Phenotype2Clean %>% filter(minor_AF >= .05)
nrow(Phenotype2Clean)

#calculates the FDR for men
HearingLossPdiffvalues <- unname(unlist(Phenotype1Clean[,"pval"]))
fdr <- p.adjust(HearingLossPdiffvalues, method = "fdr")
Phenotype1Clean$MEN.pval.FDR <- fdr

#calculates the FDR for women
HearingLossPdiffvalues <- unname(unlist(Phenotype2Clean[,"pval"]))
fdr <- p.adjust(HearingLossPdiffvalues, method = "fdr")
Phenotype2Clean$WOMEN.pval.FDR <- fdr

#total table with tSDS, FDR, and Bonferroni
write.table(Phenotype1Clean, "MenHearingLoss.tSDS.FDR.tsv", row.names = FALSE)
write.table(Phenotype2Clean, "WomenHearingLoss.tSDS.FDR.tsv", row.names = FALSE)

MenSmall <- Phenotype1Clean
WomenSmall <- Phenotype2Clean

#selects only the columns needed to create the Manhattan plot
MenSmall <- MenSmall %>% select(rsid, CHR, POS, pval, MEN.pval.FDR)
WomenSmall <- WomenSmall %>% select(rsid, CHR, POS, pval, WOMEN.pval.FDR)

#filters out high pvalues so they 
MenSmall <- MenSmall %>% filter(-log10(pval) >1 )
WomenSmall <- WomenSmall %>% filter(-log10(pval) >1 )

write.table(MenSmall, "HearingLossMenManhattan.tsv", row.names = FALSE)
write.table(WomenSmall, "HearingLossWomenManhattan.tsv", row.names = FALSE)

#---------------------------------------------------------------------
#this part of the script is used to get the Spearman correlation coefficient (rho)
#---------------------------------------------------------------------

HearingLossMen <- read.table("MenHearingLoss.tSDS.FDR.tsv", header = TRUE)
HearingLossWomen <- read.table ("WomenHearingLoss.tSDS.FDR.tsv", header = TRUE)

HearingLossMen2 <- HearingLossMen %>% select(rsid, CHR, POS, AA, DA, variant, n_complete_samples, beta, se, tstat, pval, tSDS, MEN.pval.FDR)
HearingLossWomen2 <- HearingLossWomen %>% select(rsid, CHR, POS, AA, DA, variant, n_complete_samples, beta, se, tstat, pval, tSDS, WOMEN.pval.FDR)

#changes the name of the column from Beta to differentiate for later
names(HearingLossMen2)[8] <- "BETAM"
names(HearingLossWomen2)[8] <- "BETAW"

#makes a table with only the marker name and beta values
BETAHearingLossMen <- HearingLossMen2 %>% select(rsid, BETAM) #need to put in the MarkerName
BETAHearingLossWomen <- HearingLossWomen2 %>% select(rsid, BETAW) #need to put in the MarkerName

#joins the tables so that there is a column of markernames with the beta values for
##men and women right next to each other
BETAHearingLoss <- merge(BETAHearingLossMen, BETAHearingLossWomen, by = "rsid") 

#Spearman rank correlation test
HearingLossCorr1 <- cor.test(~ BETAM + BETAW, data = BETAHearingLoss, method = 'spearman', exact = FALSE)
HearingLossCorr1
#grabs the rho value, needed later for the pdiff equation 
HearingLossCorr1$estimate

#adds a column to both the men and women dataframes of the rho value
HearingLossMen2=mutate(HearingLossMen2, rho=HearingLossCorr1$estimate)
HearingLossWomen2=mutate(HearingLossWomen2, rho=HearingLossCorr1$estimate)

#makes a table for these that can be uploaded
write.table(HearingLossMen2, "HearingLossMenLessCol.tsv", row.names= FALSE, quote = FALSE)
write.table(HearingLossWomen2, "HearingLossWomenLessCol.tsv", row.names= FALSE, quote = FALSE)

#-----------------------------------------
#this part of the script is used to get the pdiff value
#-----------------------------------------

HearingLossMen2 <- read.table("HearingLossMenLessCol.tsv", header = TRUE)
HearingLossWomen2 <- read.table("HearingLossWomenLessCol.tsv", header = TRUE)

tvalueF=function(bmen,bwomen,SEMen,SEWomen,rho){
  return ((bmen - bwomen) / (sqrt((SEMen)^2 + (SEWomen)^2 - (2*rho*SEMen*SEWomen))))
}

#creates a function to give the pvalue from the tvalue
#make sure to include not only the values needed to calculate pvalue, but also those for tvalue
pvalueF=function(bmen,bwomen,SEMen,SEWomen,rho,DegreesFreedom){
  return(2*pt(-abs(tvalueF(bmen,bwomen,SEMen,SEWomen,rho)), (DegreesFreedom-1)))
}

names(HearingLossMen2)[7] <- "MEN.N"
names(HearingLossMen2)[8] <- "MEN.beta"
names(HearingLossMen2)[9] <- "MEN.se"
names(HearingLossMen2)[10] <- "MEN.tstat"
names(HearingLossMen2)[11] <- "MEN.pval"
names(HearingLossMen2)[12] <- "MEN.tSDS"

names(HearingLossWomen2)[7] <- "WOMEN.N"
names(HearingLossWomen2)[8] <- "WOMEN.beta"
names(HearingLossWomen2)[9] <- "WOMEN.se"
names(HearingLossWomen2)[10] <- "WOMEN.tstat"
names(HearingLossWomen2)[11] <- "WOMEN.pval"
names(HearingLossWomen2)[12] <- "WOMEN.tSDS"

#concatenates the men and women values into the same datatable
UKBBHearingLoss <- merge(HearingLossMen2, HearingLossWomen2)

#cleaning data
UKBBHearingLoss$MEN.beta<-as.numeric(as.character(UKBBHearingLoss$MEN.beta)) #beta
UKBBHearingLoss$WOMEN.beta<-as.numeric(as.character(UKBBHearingLoss$WOMEN.beta)) #beta
UKBBHearingLoss$MEN.se <-as.numeric(as.character(UKBBHearingLoss$MEN.se)) #standarderror
UKBBHearingLoss$WOMEN.se <-as.numeric(as.character(UKBBHearingLoss$WOMEN.se)) #standarderror
UKBBHearingLoss$MEN.N<-as.numeric(as.character(UKBBHearingLoss$MEN.N)) #number of individuals
UKBBHearingLoss$WOMEN.N<-as.numeric(as.character(UKBBHearingLoss$WOMEN.N)) #number of individuals
UKBBHearingLoss$MEN.tSDS<-as.numeric(as.character(UKBBHearingLoss$MEN.tSDS)) #tSDS value
UKBBHearingLoss$WOMEN.tSDS<-as.numeric(as.character(UKBBHearingLoss$WOMEN.tSDS)) #tSDS value
UKBBHearingLoss$rho<-as.numeric(as.character(UKBBHearingLoss$rho)) #rho value from previous script


#mutate a new variable "pval" which holds the result of the function pvalueF
#make sure that the arguments passed into pvalueF are in the right order
#pvalueF(bmen,bwomen,SEMen,SEWomen,rho,DegreesFreedom) <-- these are the arguments in order
#the argument for DegreesFreedom is (BMIadjRandallRaw$MEN.N+BMIadjRandallRaw$WOMEN.N)
UKBBHearingLoss=mutate(UKBBHearingLoss,pdiff=pvalueF(UKBBHearingLoss$MEN.beta,UKBBHearingLoss$WOMEN.beta,UKBBHearingLoss$MEN.se,UKBBHearingLoss$WOMEN.se,UKBBHearingLoss$rho,(UKBBHearingLoss$MEN.N+UKBBHearingLoss$WOMEN.N)))

write.table(UKBBHearingLoss, "UKBBHearingLosswithPDIFF.tsv", row.names = FALSE, quote = FALSE)

#------------------------------------------------------
#this part of the script is used to choose SNPs at various cutoffs
#------------------------------------------------------

library(dplyr) #needed for %>%

#reads in the dataframe
ConcatenatedSexHearingLoss <- read.table("UKBBHearingLosswithPDIFF.tsv", header = TRUE)

#assigns the genome-wide significance value, P=5x10-8, as our cutoff line
cutoff <- 5e-8

#applies the cutoff to the women values and makes a new dataframe
WomenSig<- ConcatenatedSexHearingLoss %>% filter(WOMEN.pval <= cutoff)

#applies the cutoff to the men values and makes a new dataframe
MenSig<- ConcatenatedSexHearingLoss %>% filter(MEN.pval <= cutoff)

AllSig <- merge(MenSig, WomenSig, all=TRUE)

#calculates the FDR for the pdiff values
HearingLossPdiffFDR <- unname(unlist(AllSig[,"pdiff"]))
fdr <- p.adjust(HearingLossPdiffFDR, method = "fdr")
AllSig$pdiff.FDR <- fdr

write.table(AllSig, "HearingLosswithPDIFF.FDR.tsv", row.names = FALSE)

#applies different FDR cutoffs for pdiff and writes a table
Ex1 <- AllSig %>% filter(pdiff.FDR <= .05)
write.table(Ex1, "HearingLossPdiff1.tsv", row.names = FALSE)

Ex2 <- AllSig %>% filter(pdiff.FDR <= .01)
write.table(Ex2, "HearingLossPdiff2.tsv", row.names = FALSE)

Ex3 <- AllSig %>% filter(pdiff.FDR <= .005)
write.table(Ex3, "HearingLossPdiff3.tsv", row.names = FALSE)

Ex4 <- AllSig %>% filter(pdiff.FDR <= .001)
write.table(Ex4, "HearingLossPdiff4.tsv", row.names = FALSE)


#------------------------------------------------------
#this part of the script calculates SDS
#------------------------------------------------------

tSDSTableMen <- rep(NA, 4)
SNPsTableMen <- rep(NA, 4)

tSDSTableWomen <- rep(NA, 4)
SNPsTableWomen <- rep(NA, 4)

FDRtable <- read.table("HearingLossPdiff4.tsv", header = TRUE)

FDRtableWomen <- FDRtable %>% filter(WOMEN.pval <= MEN.pval)

tSDSTableWomen[1] <- mean(FDRtableWomen$WOMEN.tSDS)
SNPsTableWomen[1] <- nrow(FDRtableWomen)

FDRtableMen <- FDRtable %>% filter(MEN.pval < WOMEN.pval)

tSDSTableMen[1] <- mean(FDRtableMen$MEN.tSDS)
SNPsTableMen[1] <- nrow(FDRtableMen)

##########################
FDRtable <- read.table("HearingLossPdiff3.tsv", header = TRUE)

FDRtableWomen <- FDRtable %>% filter(WOMEN.pval <= MEN.pval)

tSDSTableWomen[2] <- mean(FDRtableWomen$WOMEN.tSDS)
SNPsTableWomen[2] <- nrow(FDRtableWomen)

FDRtableMen <- FDRtable %>% filter(MEN.pval < WOMEN.pval)

tSDSTableMen[2] <- mean(FDRtableMen$MEN.tSDS)
SNPsTableMen[2] <- nrow(FDRtableMen)
##########################
FDRtable <- read.table("HearingLossPdiff2.tsv", header = TRUE)

FDRtableWomen <- FDRtable %>% filter(WOMEN.pval <= MEN.pval)

tSDSTableWomen[3] <- mean(FDRtableWomen$WOMEN.tSDS)
SNPsTableWomen[3] <- nrow(FDRtableWomen)

FDRtableMen <- FDRtable %>% filter(MEN.pval < WOMEN.pval)

tSDSTableMen[3] <- mean(FDRtableMen$MEN.tSDS)
SNPsTableMen[3] <- nrow(FDRtableMen)


##########################
FDRtable <- read.table("HearingLossPdiff1.tsv", header = TRUE)

FDRtableWomen <- FDRtable %>% filter(WOMEN.pval <= MEN.pval)

tSDSTableWomen[4] <- mean(FDRtableWomen$WOMEN.tSDS)
SNPsTableWomen[4] <- nrow(FDRtableWomen)

FDRtableMen <- FDRtable %>% filter(MEN.pval < WOMEN.pval)

tSDSTableMen[4] <- mean(FDRtableMen$MEN.tSDS)
SNPsTableMen[4] <- nrow(FDRtableMen)

write.table(tSDSTableMen, "HearingLosstSDSTableMen.txt")
write.table(SNPsTableMen, "HearingLossSNPsTableMen.txt")

write.table(tSDSTableWomen, "HearingLosstSDSTableWomen.txt")
write.table(SNPsTableWomen, "HearingLossSNPsTableWomen.txt")

