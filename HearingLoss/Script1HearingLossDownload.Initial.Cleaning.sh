#!/bin/sh

#===============================================
#--FILE DOWNLOAD AND EARLY DATA CLEANING
#===============================================

#downloads the GWAS data from the Neale lab
wget https://www.dropbox.com/s/js788fn3jyrry4r/H90.gwas.imputed_v3.female.tsv.bgz?dl=0 -O H90.gwas.imputed_v3.female.tsv.bgz
wget https://www.dropbox.com/s/hfbbdkk7g1oqp8v/H90.gwas.imputed_v3.male.tsv.bgz?dl=0 -O H90.gwas.imputed_v3.male.tsv.bgz
#downloads all the variants, which are of the same number and order in all GWAS files downloaded from the Neale lab
wget https://www.dropbox.com/s/puxks683vb0omeg/variants.tsv.bgz?dl=0 -O variants.tsv.bgz

#renames the files so that they can be decompressed
mv H90.gwas.imputed_v3.female.tsv.bgz H90.gwas.imputed_v3.female.tsv.gz
mv H90.gwas.imputed_v3.male.tsv.bgz H90.gwas.imputed_v3.male.tsv.gz
mv variants.tsv.bgz variants.tsv.gz

#decompresses files
gunzip H90.gwas.imputed_v3.female.tsv.gz
gunzip H90.gwas.imputed_v3.male.tsv.gz
gunzip variants.tsv.gz

#isolates the rsid column from the variant file
cut -f6 variants.tsv > variantRSID.tsv

#adds the RSID to the phenotype files
paste variantRSID.tsv H90.gwas.imputed_v3.female.tsv > H90.gwas.FemalewithRSID.tsv
paste variantRSID.tsv H90.gwas.imputed_v3.male.tsv > H90.gwas.MalewithRSID.tsv

exit;