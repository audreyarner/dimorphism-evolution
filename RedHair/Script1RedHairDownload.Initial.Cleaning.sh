#!/bin/sh

#===============================================
#--FILE DOWNLOAD AND EARLY DATA CLEANING
#===============================================

#downloads the GWAS data from the Neale lab
wget https://www.dropbox.com/s/od69ykxn35o5047/1747_2.gwas.imputed_v3.female.tsv.bgz?dl=0 -O 1747_2.gwas.imputed_v3.female.tsv.bgz
wget https://www.dropbox.com/s/456jrd5yl0ljy0g/1747_2.gwas.imputed_v3.male.tsv.bgz?dl=0 -O 1747_2.gwas.imputed_v3.male.tsv.bgz
#downloads all the variants, which are of the same number and order in all GWAS files downloaded from the Neale lab
wget https://www.dropbox.com/s/puxks683vb0omeg/variants.tsv.bgz?dl=0 -O variants.tsv.bgz

#renames the files so that they can be decompressed
mv 1747_2.gwas.imputed_v3.female.tsv.bgz 1747_2.gwas.imputed_v3.female.tsv.gz
mv 1747_2.gwas.imputed_v3.male.tsv.bgz 1747_2.gwas.imputed_v3.male.tsv.gz
mv variants.tsv.bgz variants.tsv.gz

#decompresses files
gunzip 1747_2.gwas.imputed_v3.female.tsv.gz
gunzip 1747_2.gwas.imputed_v3.male.tsv.gz
gunzip variants.tsv.gz

#isolates the rsid column from the variant file
cut -f6 variants.tsv > variantRSID.tsv

#adds the RSID to the phenotype files
paste variantRSID.tsv 1747_2.gwas.imputed_v3.female.tsv > 1747_2.gwas.FemalewithRSID.tsv
paste variantRSID.tsv 1747_2.gwas.imputed_v3.male.tsv > 1747_2.gwas.MalewithRSID.tsv

exit;