#!/bin/bash
#$ -l h_rt=07:00:00
#$ -l rmem=2g
#$ -j y
#$ -o gemma_bslmm.log

GEMMA='gemma'

DIR="/data/$USER/gwas_gemma"
GENOTYPES='fha.bbgeno.gz'
PHENOTYPES='data/fha.pheno'
RELMATRIX='output/relmatrix.cXX.txt'
BSLMM=1 # 1=BSLMM, 2=standard ridge regression/GBLUP (no mcmc), 3=probit BSLMM (requires 0/1 phenotypes)
OUTBASE='bslmm'

# Number of variants that will be used (required for setting some priors and proposals below)
NVARS=$(grep "number of analyzed SNPs" $DIR/output/relmatrix.log.txt | perl -pe 's/.*\= //g')

# priors
# ----------------
# h -> approximation to PVE: proportion of phenotypic variance explained by loci
HMIN=0
HMAX=1

# rho -> approximation to PGE: proportion of genetic variance explained by sparse effect terms (~major effect loci)
# rho=0 -> pure LMM, highly polygenic; rho=1 => pure BVSR, few loci
RHOMIN=0
RHOMAX=1

# pi -> proportion of variants with non-zero effects (random + sparse effects)
PIMIN=$(Rscript -e 'cat(log10(1/'$NVARS'))') # log10(1/no_variants)
PIMAX=0 # log10(1)

# gamma -> Number of variants with sparse effects (~ number of major effect loci)
GAMMAMIN=0
GAMMAMAX=300
# ----------------

# proposals
# ----------------
# don't need to tweak them unless you have convergence problems
GEOMMEAN=2000
HSTEP=$(Rscript -e 'cat(min(c(10/sqrt('$NVARS'),1)))') # 0-1, default: min(10/sqrt(no_variants),1)
RHOSTEP=$(Rscript -e 'cat(min(c(10/sqrt('$NVARS'),1)))') # 0-1, default: min(10/sqrt(no_variants),1)
PISTEP=$(Rscript -e 'cat(min(c(5/sqrt('$NVARS'),1)))') # 0-1, default: min(5/sqrt(n),1)
# ----------------

# chain parameters
# ----------------
BURNIN=250000 # No MCMC initial steps to be discarded (suggested: 10-25% MCCMC length)
MCMCLEN=1000000 # No MCMC steps after burnin
RECORDPACE=100 # Record states every X steps
WRITEPACE=1000 # Write to file every X steps (suggested: >=MCMCLEN/1000)
# ----------------

# QC filters
# ----------------
MAF='0.01' # exclude very rare variants
# ----------------

hostname
uname -a
date
echo "----------------------------------------------------------"
echo
cd $DIR

$GEMMA \
-g $GENOTYPES \
-p $PHENOTYPES \
-k $RELMATRIX \
-bslmm $BSLMM \
-w $BURNIN \
-s $MCMCLEN \
-rpace $RECORDPACE \
-wpace $WRITEPACE \
-maf $MAF \
-o $OUTBASE

# This is a commented example in case you would like
# to specify specify priors and proposals
# $GEMMA \
# -g $GENOTYPES \
# -p $PHENOTYPES \
# -k $RELMATRIX \
# -bslmm $BSLMM \
# -hmin $HMIN \
# -hmax $HMAX \
# -rmin $RHOMIN \
# -rmax $RHOMAX \
# -pmin $PIMIN \
# -pmax $PIMAX \
# -gmean $GEOMMEAN \
# -hscale $HSTEP \
# -rscale $RHOSTEP \
# -pscale $PISTEP \
# -w $BURNIN \
# -s $MCMCLEN \
# -rpace $RECORDPACE \
# -wpace $WRITEPACE \
# -maf $MAF \
# -o $OUTBASE

echo
echo "----------------------------------------------------------"
date
