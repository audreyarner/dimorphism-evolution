#!/bin/sh

#===============================================
#--FILE DOWNLOAD AND EARLY DATA CLEANING
#===============================================

#downloads the GWAS data from the Neale lab
wget https://www.dropbox.com/s/s2n334txw4ruz6j/21002_irnt.gwas.imputed_v3.female.tsv.bgz?dl=0 -O 21002_irnt.gwas.imputed_v3.female.tsv.bgz
wget https://www.dropbox.com/s/afgenopl375dhy4/21002_irnt.gwas.imputed_v3.male.tsv.bgz?dl=0 -O 21002_irnt.gwas.imputed_v3.male.tsv.bgz
#downloads all the variants, which are of the same number and order in all GWAS files downloaded from the Neale lab
wget https://www.dropbox.com/s/puxks683vb0omeg/variants.tsv.bgz?dl=0 -O variants.tsv.bgz

#renames the files so that they can be decompressed
mv 21002_irnt.gwas.imputed_v3.male.tsv.bgz 21002_irnt.gwas.imputed_v3.male.tsv.gz
mv 21002_irnt.gwas.imputed_v3.female.tsv.bgz 21002_irnt.gwas.imputed_v3.female.tsv.gz
mv variants.tsv.bgz variants.tsv.gz

#decompresses files
gunzip 21002_irnt.gwas.imputed_v3.male.tsv.gz
gunzip 21002_irnt.gwas.imputed_v3.female.tsv.gz
gunzip variants.tsv.gz

#isolates the rsid column from the variant file
cut -f6 variants.tsv > variantRSID.tsv

#adds the RSID to the phenotype files
paste variantRSID.tsv 21002_irnt.gwas.imputed_v3.male.tsv > 21002_irnt.MalewithRSID.tsv
paste variantRSID.tsv 21002_irnt.gwas.imputed_v3.female.tsv > 21002_irnt.FemalewithRSID.tsv

exit;