#!/bin/bash
#$ -l h_rt=4:00:00
#$ -l rmem=8g
#$ -j y
#$ -o bbgeno.log


# Estimate genotype probabilities from genotype likelihoods in bcf/vcf
# and convert to mean genotype format (required for gemma)

BCF2BBGENO="/data/$USER/gwas_gemma/scripts/bcf2bbgeno.pl"

INDIR="/data/$USER/fst_hmm"

# input file
INPUT='data/fha.vcf.gz'
# output file
OUTPUT='fha.bbgeno'


hostname
uname -a
date
echo "----------------------------------------------------------" 
echo

$BCF2BBGENO -i $INPUT -o $OUTPUT -p H-W -s -r

echo
echo "----------------------------------------------------------"
date

