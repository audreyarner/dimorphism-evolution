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
UK10KDataSet <- read.table("/storage/home/ama6560/work/SDS_UK10K_n3195_release_Sep_19_2016.txt", header = TRUE)
#the phenotype data you are interested in. Should have ~13 million rows
Phenotype <- read.table("1309_irnt.gwas.FemalewithRSID.tsv", header = TRUE)
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
write.csv(Final.Phenotype.tSDS, "WomenOsteoporosiswithtSDS.csv", row.names = FALSE)

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
write.csv(FinalBinSDSWOOD, "WomenOsteoporosisSDSBins.csv", row.names = FALSE)

ReplicatedFieldGraph <- read.csv("WomenOsteoporosisSDSBins.csv")

#makes a graph of the genome wide tSDS score of bins of 1000 SNPs, as seen in Field et al Figure 4A
ReplicatedFieldGraph %>% ggplot(aes(x=MeanP, y=MeantSDS)) + geom_point() + ggtitle("Phenotype tSDS score observed genome-wide")+ labs(x="Phenotype signifcance (rank)", y = "Phenotype-increasing tSDS (bin average)") + scale_x_reverse() + ylim(-.3,.75)

ggsave("WomenOsteoporosis.tSDSGraph.png")

#====================================================================
#--The following repeats the process for males
#====================================================================

Phenotype <- read.table("1309_irnt.gwas.MalewithRSID.tsv", header = TRUE)
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
write.csv(Final.Phenotype.tSDS, "MenOsteoporosiswithtSDS.csv", row.names = FALSE)

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
write.csv(FinalBinSDSWOOD, "MenOsteoporosisSDSBins.csv", row.names = FALSE)

ReplicatedFieldGraph <- read.csv("MenOsteoporosisSDSBins.csv")

#makes a graph of the genome wide tSDS score of bins of 1000 SNPs, as seen in Field et al Figure 4A
ReplicatedFieldGraph %>% ggplot(aes(x=MeanP, y=MeantSDS)) + geom_point() + ggtitle("Phenotype tSDS score observed genome-wide")+ labs(x="Phenotype signifcance (rank)", y = "Phenotype-increasing tSDS (bin average)") + scale_x_reverse() + ylim(-.3,.75)

ggsave("MenOsteoporosis.tSDSGraph.png")

WomenPheno <- read.csv("WomenOsteoporosiswithtSDS.csv", header = TRUE)
MenPheno <- read.csv("MenOsteoporosiswithtSDS.csv", header = TRUE)

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
OsteoporosisPdiffvalues <- unname(unlist(Phenotype1Clean[,"pval"]))
fdr <- p.adjust(OsteoporosisPdiffvalues, method = "fdr")
Phenotype1Clean$MEN.pval.FDR <- fdr

#calculates the FDR for women
OsteoporosisPdiffvalues <- unname(unlist(Phenotype2Clean[,"pval"]))
fdr <- p.adjust(OsteoporosisPdiffvalues, method = "fdr")
Phenotype2Clean$WOMEN.pval.FDR <- fdr

#total table with tSDS, FDR, and Bonferroni
write.table(Phenotype1Clean, "MenOsteoporosis.tSDS.FDR.tsv", row.names = FALSE)
write.table(Phenotype2Clean, "WomenOsteoporosis.tSDS.FDR.tsv", row.names = FALSE)

MenSmall <- Phenotype1Clean
WomenSmall <- Phenotype2Clean

#selects only the columns needed to create the Manhattan plot
MenSmall <- MenSmall %>% select(rsid, CHR, POS, pval, MEN.pval.FDR)
WomenSmall <- WomenSmall %>% select(rsid, CHR, POS, pval, WOMEN.pval.FDR)

#filters out high pvalues so they 
MenSmall <- MenSmall %>% filter(-log10(pval) >1 )
WomenSmall <- WomenSmall %>% filter(-log10(pval) >1 )

write.table(MenSmall, "OsteoporosisMenManhattan.tsv", row.names = FALSE)
write.table(WomenSmall, "OsteoporosisWomenManhattan.tsv", row.names = FALSE)

#---------------------------------------------------------------------
#this part of the script is used to get the Spearman correlation coefficient (rho)
#---------------------------------------------------------------------

OsteoporosisMen <- read.table("MenOsteoporosis.tSDS.FDR.tsv", header = TRUE)
OsteoporosisWomen <- read.table ("WomenOsteoporosis.tSDS.FDR.tsv", header = TRUE)

OsteoporosisMen2 <- OsteoporosisMen %>% select(rsid, CHR, POS, AA, DA, variant, n_complete_samples, beta, se, tstat, pval, tSDS, MEN.pval.FDR)
OsteoporosisWomen2 <- OsteoporosisWomen %>% select(rsid, CHR, POS, AA, DA, variant, n_complete_samples, beta, se, tstat, pval, tSDS, WOMEN.pval.FDR)

#changes the name of the column from Beta to differentiate for later
names(OsteoporosisMen2)[8] <- "BETAM"
names(OsteoporosisWomen2)[8] <- "BETAW"

#makes a table with only the marker name and beta values
BETAOsteoporosisMen <- OsteoporosisMen2 %>% select(rsid, BETAM) #need to put in the MarkerName
BETAOsteoporosisWomen <- OsteoporosisWomen2 %>% select(rsid, BETAW) #need to put in the MarkerName

#joins the tables so that there is a column of markernames with the beta values for
##men and women right next to each other
BETAOsteoporosis <- merge(BETAOsteoporosisMen, BETAOsteoporosisWomen, by = "rsid") 

#Spearman rank correlation test
OsteoporosisCorr1 <- cor.test(~ BETAM + BETAW, data = BETAOsteoporosis, method = 'spearman', exact = FALSE)
OsteoporosisCorr1
#grabs the rho value, needed later for the pdiff equation 
OsteoporosisCorr1$estimate

#adds a column to both the men and women dataframes of the rho value
OsteoporosisMen2=mutate(OsteoporosisMen2, rho=OsteoporosisCorr1$estimate)
OsteoporosisWomen2=mutate(OsteoporosisWomen2, rho=OsteoporosisCorr1$estimate)

#makes a table for these that can be uploaded
write.table(OsteoporosisMen2, "OsteoporosisMenLessCol.tsv", row.names= FALSE, quote = FALSE)
write.table(OsteoporosisWomen2, "OsteoporosisWomenLessCol.tsv", row.names= FALSE, quote = FALSE)

#-----------------------------------------
#this part of the script is used to get the pdiff value
#-----------------------------------------

OsteoporosisMen2 <- read.table("OsteoporosisMenLessCol.tsv", header = TRUE)
OsteoporosisWomen2 <- read.table("OsteoporosisWomenLessCol.tsv", header = TRUE)

tvalueF=function(bmen,bwomen,SEMen,SEWomen,rho){
  return ((bmen - bwomen) / (sqrt((SEMen)^2 + (SEWomen)^2 - (2*rho*SEMen*SEWomen))))
}

#creates a function to give the pvalue from the tvalue
#make sure to include not only the values needed to calculate pvalue, but also those for tvalue
pvalueF=function(bmen,bwomen,SEMen,SEWomen,rho,DegreesFreedom){
  return(2*pt(-abs(tvalueF(bmen,bwomen,SEMen,SEWomen,rho)), (DegreesFreedom-1)))
}

names(OsteoporosisMen2)[7] <- "MEN.N"
names(OsteoporosisMen2)[8] <- "MEN.beta"
names(OsteoporosisMen2)[9] <- "MEN.se"
names(OsteoporosisMen2)[10] <- "MEN.tstat"
names(OsteoporosisMen2)[11] <- "MEN.pval"
names(OsteoporosisMen2)[12] <- "MEN.tSDS"

names(OsteoporosisWomen2)[7] <- "WOMEN.N"
names(OsteoporosisWomen2)[8] <- "WOMEN.beta"
names(OsteoporosisWomen2)[9] <- "WOMEN.se"
names(OsteoporosisWomen2)[10] <- "WOMEN.tstat"
names(OsteoporosisWomen2)[11] <- "WOMEN.pval"
names(OsteoporosisWomen2)[12] <- "WOMEN.tSDS"

#concatenates the men and women values into the same datatable
UKBBOsteoporosis <- merge(OsteoporosisMen2, OsteoporosisWomen2)

#cleaning data
UKBBOsteoporosis$MEN.beta<-as.numeric(as.character(UKBBOsteoporosis$MEN.beta)) #beta
UKBBOsteoporosis$WOMEN.beta<-as.numeric(as.character(UKBBOsteoporosis$WOMEN.beta)) #beta
UKBBOsteoporosis$MEN.se <-as.numeric(as.character(UKBBOsteoporosis$MEN.se)) #standarderror
UKBBOsteoporosis$WOMEN.se <-as.numeric(as.character(UKBBOsteoporosis$WOMEN.se)) #standarderror
UKBBOsteoporosis$MEN.N<-as.numeric(as.character(UKBBOsteoporosis$MEN.N)) #number of individuals
UKBBOsteoporosis$WOMEN.N<-as.numeric(as.character(UKBBOsteoporosis$WOMEN.N)) #number of individuals
UKBBOsteoporosis$MEN.tSDS<-as.numeric(as.character(UKBBOsteoporosis$MEN.tSDS)) #tSDS value
UKBBOsteoporosis$WOMEN.tSDS<-as.numeric(as.character(UKBBOsteoporosis$WOMEN.tSDS)) #tSDS value
UKBBOsteoporosis$rho<-as.numeric(as.character(UKBBOsteoporosis$rho)) #rho value from previous script


#mutate a new variable "pval" which holds the result of the function pvalueF
#make sure that the arguments passed into pvalueF are in the right order
#pvalueF(bmen,bwomen,SEMen,SEWomen,rho,DegreesFreedom) <-- these are the arguments in order
#the argument for DegreesFreedom is (BMIadjRandallRaw$MEN.N+BMIadjRandallRaw$WOMEN.N)
UKBBOsteoporosis=mutate(UKBBOsteoporosis,pdiff=pvalueF(UKBBOsteoporosis$MEN.beta,UKBBOsteoporosis$WOMEN.beta,UKBBOsteoporosis$MEN.se,UKBBOsteoporosis$WOMEN.se,UKBBOsteoporosis$rho,(UKBBOsteoporosis$MEN.N+UKBBOsteoporosis$WOMEN.N)))

write.table(UKBBOsteoporosis, "UKBBOsteoporosiswithPDIFF.tsv", row.names = FALSE, quote = FALSE)

#------------------------------------------------------
#this part of the script is used to choose SNPs at various cutoffs
#------------------------------------------------------

library(dplyr) #needed for %>%

#reads in the dataframe
ConcatenatedSexOsteoporosis <- read.table("UKBBOsteoporosiswithPDIFF.tsv", header = TRUE)

#assigns the genome-wide significance value, P=5x10-8, as our cutoff line
cutoff <- 5e-8

#applies the cutoff to the women values and makes a new dataframe
WomenSig<- ConcatenatedSexOsteoporosis %>% filter(WOMEN.pval <= cutoff)

#applies the cutoff to the men values and makes a new dataframe
MenSig<- ConcatenatedSexOsteoporosis %>% filter(MEN.pval <= cutoff)

AllSig <- merge(MenSig, WomenSig, all=TRUE)

#calculates the FDR for the pdiff values
OsteoporosisPdiffFDR <- unname(unlist(AllSig[,"pdiff"]))
fdr <- p.adjust(OsteoporosisPdiffFDR, method = "fdr")
AllSig$pdiff.FDR <- fdr

write.table(AllSig, "OsteoporosiswithPDIFF.FDR.tsv", row.names = FALSE)

#applies different FDR cutoffs for pdiff and writes a table
Ex1 <- AllSig %>% filter(pdiff.FDR <= .05)
write.table(Ex1, "OsteoporosisPdiff1.tsv", row.names = FALSE)

Ex2 <- AllSig %>% filter(pdiff.FDR <= .01)
write.table(Ex2, "OsteoporosisPdiff2.tsv", row.names = FALSE)

Ex3 <- AllSig %>% filter(pdiff.FDR <= .005)
write.table(Ex3, "OsteoporosisPdiff3.tsv", row.names = FALSE)

Ex4 <- AllSig %>% filter(pdiff.FDR <= .001)
write.table(Ex4, "OsteoporosisPdiff4.tsv", row.names = FALSE)


#------------------------------------------------------
#this part of the script calculates SDS
#------------------------------------------------------

tSDSTableMen <- rep(NA, 4)
SNPsTableMen <- rep(NA, 4)

tSDSTableWomen <- rep(NA, 4)
SNPsTableWomen <- rep(NA, 4)

FDRtable <- read.table("OsteoporosisPdiff4.tsv", header = TRUE)

FDRtableWomen <- FDRtable %>% filter(WOMEN.pval <= MEN.pval)

tSDSTableWomen[1] <- mean(FDRtableWomen$WOMEN.tSDS)
SNPsTableWomen[1] <- nrow(FDRtableWomen)

FDRtableMen <- FDRtable %>% filter(MEN.pval < WOMEN.pval)

tSDSTableMen[1] <- mean(FDRtableMen$MEN.tSDS)
SNPsTableMen[1] <- nrow(FDRtableMen)

##########################
FDRtable <- read.table("OsteoporosisPdiff3.tsv", header = TRUE)

FDRtableWomen <- FDRtable %>% filter(WOMEN.pval <= MEN.pval)

tSDSTableWomen[2] <- mean(FDRtableWomen$WOMEN.tSDS)
SNPsTableWomen[2] <- nrow(FDRtableWomen)

FDRtableMen <- FDRtable %>% filter(MEN.pval < WOMEN.pval)

tSDSTableMen[2] <- mean(FDRtableMen$MEN.tSDS)
SNPsTableMen[2] <- nrow(FDRtableMen)
##########################
FDRtable <- read.table("OsteoporosisPdiff2.tsv", header = TRUE)

FDRtableWomen <- FDRtable %>% filter(WOMEN.pval <= MEN.pval)

tSDSTableWomen[3] <- mean(FDRtableWomen$WOMEN.tSDS)
SNPsTableWomen[3] <- nrow(FDRtableWomen)

FDRtableMen <- FDRtable %>% filter(MEN.pval < WOMEN.pval)

tSDSTableMen[3] <- mean(FDRtableMen$MEN.tSDS)
SNPsTableMen[3] <- nrow(FDRtableMen)


##########################
FDRtable <- read.table("OsteoporosisPdiff1.tsv", header = TRUE)

FDRtableWomen <- FDRtable %>% filter(WOMEN.pval <= MEN.pval)

tSDSTableWomen[4] <- mean(FDRtableWomen$WOMEN.tSDS)
SNPsTableWomen[4] <- nrow(FDRtableWomen)

FDRtableMen <- FDRtable %>% filter(MEN.pval < WOMEN.pval)

tSDSTableMen[4] <- mean(FDRtableMen$MEN.tSDS)
SNPsTableMen[4] <- nrow(FDRtableMen)

write.table(tSDSTableMen, "OsteoporosistSDSTableMen.txt")
write.table(SNPsTableMen, "OsteoporosisSNPsTableMen.txt")

write.table(tSDSTableWomen, "OsteoporosistSDSTableWomen.txt")
write.table(SNPsTableWomen, "OsteoporosisSNPsTableWomen.txt")

