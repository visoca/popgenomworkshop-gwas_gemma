#!/bin/bash
#$ -l h_rt=1:00:00
#$ -j y
#$ -o gemma_relmatrix.log

GEMMA='gemma'
DIR="/data/$USER/gwas_gemma"
GENOTYPES='fha.bbgeno.gz'
PHENOTYPES='data/fha.pheno'
# centered matrix preferred in general, accounts better for population structure
# standardized matrix preferred if SNPs with lower MAF have larger effects 
MATRIXTYPE=1 # 1=centered matrix, 2=standardized matrix
OUTBASE='relmatrix'

hostname
uname -a
date
echo "----------------------------------------------------------"
echo
cd $DIR

$GEMMA \
-g $GENOTYPES \
-p $PHENOTYPES \
-gk $MATRIXTYPE \
-o $OUTBASE

echo
echo "----------------------------------------------------------"
date
