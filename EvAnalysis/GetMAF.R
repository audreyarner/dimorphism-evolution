#!/usr/bin/env Rscript

#=========================================================================
#--Get minor allele frequency needed for evolutionary analysis
#=========================================================================
library(dplyr)

getMAF=function(PrunedSNPs, MAFTableFemale, MAFTableMale){
	Pruned <- read.table(PrunedSNPs, header =T)
	PhenotypeName <- gsub("PrunedSNPs.txt", "", PrunedSNPs)

	MAF.Female <- read.csv(MAFTableFemale, header =T)
	MAF.Male <- read.csv(MAFTableMale, header =T)

	MAF.Female2 <- MAF.Female %>% select(rsid,minor_AF)
	names(MAF.Female2)[2] <- "Female.maf"

	MAF.Male2 <- MAF.Male %>% select(rsid,minor_AF)
	names(MAF.Male2)[2] <- "Male.maf"

	BothMAF <- merge(MAF.Female2, MAF.Male2)

	SNPsandMAF <- merge(BothMAF, Pruned, by = "rsid")

	write.table(SNPsandMAF, file=paste(PhenotypeName, "MAF.txt", sep=""),sep="\t", row.names =F)

}

getMAF("HeightPdiffPrunedSNPs.txt", "FemaleHeightwithtSDS.csv", "MaleHeightwithtSDS.csv")
getMAF("HeightPrunedSNPs.txt", "FemaleHeightwithtSDS.csv", "MaleHeightwithtSDS.csv")
getMAF("HipCircPdiffPrunedSNPs.txt", "FemaleHipCircwithtSDS.csv", "MaleHipCircwithtSDS.csv")
getMAF("HipCircPrunedSNPs.txt", "FemaleHipCircwithtSDS.csv", "MaleHipCircwithtSDS.csv")
getMAF("WaistCircPdiffPrunedSNPs.txt", "FemaleWaistCircwithtSDS.csv", "MaleWaistCircwithtSDS.csv")
getMAF("WaistCircPrunedSNPs.txt", "FemaleWaistCircwithtSDS.csv", "MaleWaistCircwithtSDS.csv")
getMAF("BodyMassPdiffPrunedSNPs.txt", "FemaleBodyMasswithtSDS.csv", "MaleBodyMasswithtSDS.csv")
getMAF("BodyMassPrunedSNPs.txt", "FemaleBodyMasswithtSDS.csv", "MaleBodyMasswithtSDS.csv")
getMAF("BodyFatPercentPdiffPrunedSNPs.txt", "FemaleBodyFatPercentwithtSDS.csv", "MaleBodyFatPercentwithtSDS.csv")
getMAF("BodyFatPercentPrunedSNPs.txt", "FemaleBodyFatPercentwithtSDS.csv", "MaleBodyFatPercentwithtSDS.csv")
