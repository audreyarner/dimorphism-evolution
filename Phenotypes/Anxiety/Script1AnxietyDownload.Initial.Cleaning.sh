#!/bin/sh

#===============================================
#--FILE DOWNLOAD AND EARLY DATA CLEANING
#===============================================

#downloads the GWAS data from the Neale lab
wget https://www.dropbox.com/s/00x1bdkjqa8begx/20002_1287.gwas.imputed_v3.female.tsv.bgz?dl=0 -O 20002_1287.gwas.imputed_v3.female.tsv.bgz
wget https://www.dropbox.com/s/22c8baie6irqz5v/20002_1287.gwas.imputed_v3.male.tsv.bgz?dl=0 -O 20002_1287.gwas.imputed_v3.male.tsv.bgz
#downloads all the variants, which are of the same number and order in all GWAS files downloaded from the Neale lab
wget https://www.dropbox.com/s/puxks683vb0omeg/variants.tsv.bgz?dl=0 -O variants.tsv.bgz

#renames the files so that they can be decompressed
mv 20002_1287.gwas.imputed_v3.female.tsv.bgz 20002_1287.gwas.imputed_v3.female.tsv.gz
mv 20002_1287.gwas.imputed_v3.male.tsv.bgz 20002_1287.gwas.imputed_v3.male.tsv.gz
mv variants.tsv.bgz variants.tsv.gz

#decompresses files
gunzip 20002_1287.gwas.imputed_v3.female.tsv.gz
gunzip 20002_1287.gwas.imputed_v3.male.tsv.gz
gunzip variants.tsv.gz

#isolates the rsid column from the variant file
cut -f6 variants.tsv > variantRSID.tsv

#adds the RSID to the phenotype files
paste variantRSID.tsv 20002_1287.gwas.imputed_v3.female.tsv > 1287_irnt.gwas.FemalewithRSID.tsv
paste variantRSID.tsv 20002_1287.gwas.imputed_v3.male.tsv > 1287_irnt.gwas.MalewithRSID.tsv

exit;