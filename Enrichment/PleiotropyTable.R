#!/usr/bin/env Rscript

#------------------------------------------------------------------------------------------------------
#this script is used to create a table with the phenotypes each SNP is associated with
#------------------------------------------------------------------------------------------------------

#GWAS Catalog downloaded from webiste
GWASCatalog <- read.delim("/Users/audreyarner/Box\ Sync/Lab/HeightDimorphism/PaperWriting/PLOSRevisions/GWASCatalog.tsv", header =T, sep="\t", quote = "")
#make a column with the indexes
GWASCatalog$names <- rownames(GWASCatalog)

GWASCatalog$DISEASE.TRAIT <- as.character(GWASCatalog$DISEASE.TRAIT)

#concatenated list of SexDiff-associated SNPs
SexDiff <- read.table("SexDiffSNPs.txt", header =T)

#find phenotypes associated with each SNP
SexDiff[,"Phenotypes"] <- NA

for(i in 1:nrow(SexDiff)){ 
  SNP.Overlap = which(SexDiff$rsid[i] == GWASCatalog$SNPS) 
  if(length(SNP.Overlap >=1)){
    SNPs <-GWASCatalog[(GWASCatalog$names %in% SNP.Overlap),]
    Pheno <- unique(SNPs$DISEASE.TRAIT)
    OneString <- toString(Pheno)
    SexDiff$Phenotypes[i]=OneString
  }
}

write.table(SexDiff, "SexDiffSNPsPleiotropy.txt", sep="\t", row.names = F, quote =F)

#Concatenated phenotype-associated SNPs
Pheno1 <- read.table("AllPrunedPheno.txt", header =T)

#find phenotypes associated with each SNP
Pheno1[,"Phenotypes"] <- NA

for(i in 1:nrow(Pheno1)){ 
  SNP.Overlap = which(Pheno1$rsid[i] == GWASCatalog$SNPS) 
  if(length(SNP.Overlap >=1)){
    SNPs <-GWASCatalog[(GWASCatalog$names %in% SNP.Overlap),]
    Pheno <- unique(SNPs$DISEASE.TRAIT)
    OneString <- toString(Pheno)
    Pheno1$Phenotypes[i]=OneString
  }
}

write.table(Pheno1, "PhenotypeSNPsPleiotropy.txt", sep="\t", row.names = F, quote =F)
