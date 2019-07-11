#!/bin/sh

#===============================================
#--FILE DOWNLOAD AND EARLY DATA CLEANING
#===============================================

#downloads the GWAS data from the Neale lab
wget https://www.dropbox.com/s/q10shqmyx9llwr8/1777.gwas.imputed_v3.female.tsv.bgz?dl=0 -O 1777.gwas.imputed_v3.female.tsv.bgz
wget https://www.dropbox.com/s/riqpno8r6aygkak/1777.gwas.imputed_v3.male.tsv.bgz?dl=0 -O 1777.gwas.imputed_v3.male.tsv.bgz
#downloads all the variants, which are of the same number and order in all GWAS files downloaded from the Neale lab
wget https://www.dropbox.com/s/puxks683vb0omeg/variants.tsv.bgz?dl=0 -O variants.tsv.bgz

#renames the files so that they can be decompressed
mv 1777.gwas.imputed_v3.female.tsv.bgz 1777.gwas.imputed_v3.female.tsv.gz
mv 1777.gwas.imputed_v3.male.tsv.bgz 1777.gwas.imputed_v3.male.tsv.gz
mv variants.tsv.bgz variants.tsv.gz

#decompresses files
gunzip 1777.gwas.imputed_v3.female.tsv.gz
gunzip 1777.gwas.imputed_v3.male.tsv.gz
gunzip variants.tsv.gz

#isolates the rsid column from the variant file
cut -f6 variants.tsv > variantRSID.tsv

#adds the RSID to the phenotype files
paste variantRSID.tsv 1777.gwas.imputed_v3.female.tsv > 1777.gwas.FemalewithRSID.tsv
paste variantRSID.tsv 1777.gwas.imputed_v3.male.tsv > 1777.gwas.MalewithRSID.tsv

exit;