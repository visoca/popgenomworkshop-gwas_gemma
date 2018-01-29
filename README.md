*Population Genomics Workshop 2018, University of Sheffield*

# Genetic architecture of traits using Genome Wide Association (GWA) mapping with GEMMA
#### Victor Soria-Carrasco

The aim of this practical is use a multi-locus method to study the genetic architecture of colour pattern phenotype (i.e. white stripe) in *Timema cristinae stick* insects. In particular, we will use the Bayesian Sparse Linear Mixed Model (BSLMM) model implemented in the program ```GEMMA```, a piece of software that also allow taking into account any potential population structure. This model is a hybrid approach that allows inferring whether the genetic basis of a trait is highly polygenic or rather determined by a few loci (i.e. oligogenic), as well as identifying the loci that show the stronger association with the trait.  We will use the dataset of [Comeault et al. 2015](http://www.cell.com/current-biology/fulltext/S0960-9822(15)00661-2), which consists of 602 individuals sampled from a phenotypically variable population called FHA. We will infer the proportion of phenotypic variance explained by the genotypes as well as what fraction of such variance can be explained by a reduced number of major effect loci.

Data and scripts used during this practical will be available in a shared folder in Iceberg, the University of Sheffield HPC cluster, which will allow faster transfers. However, in case attendees would like to use it for practice after the workshop, scripts are also available in [this repository](https://github.com/visoca/popgenomworkshop-gwas_gemma) and data can be downloaded [here](https://drive.google.com/file/d/12m3cHEXGUj9AK2YUJBO8a8VCsdKF7GrM/view?usp=sharing) (and results [here]()).

### Resources
* [Zhou et al. 2013](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003264) - Paper where the Bayesian Sparse Linear Mixed Model approach implemented in ```GEMMA``` is described.
* [Zhou Lab website](http://www.xzlab.org/software.html) - You can obtain GEMMA and a few other interesting programs for GWAS here.
* [GEMMA manual](http://www.xzlab.org/software/GEMMAmanual.pdf)

## 1. Initial set up
We are going to create a working directory in a dedicated space in the HPC cluster (/data/$USER) and copy the necessary scripts and data files to run this practical.

Connect to Iceberg HPC cluster (change *myuser* by your username):
```bash
ssh myuser@iceberg.sheffield.ac.uk
```
Request an interactive session:
```bash
qrsh
```
#### Important note
***
This tutorial relies on having access to a number of programs. The easiest way is to have your account configured to use the Genomics Software Repository. If that is the case you should see the following message when you get an interactive session with ```qrsh```:
```
  Your account is set up to use the Genomics Software Repository
     More info: http://soria-carrasco.staff.shef.ac.uk/softrepo
```
If you don't get that message, follow the instructions [here](http://soria-carrasco.staff.shef.ac.uk/softrepo/) to set up your account.

In addition, if you want to configure the ```nano``` text editor to have syntax highlighting and line numbering, you can configure it this way:
```bash
cat /usr/local/extras/Genomics/workshops/January2018/.nanorc >> /home/$USER/.nanorc
```
***

Change to your data directory:
```bash
cd /data/$USER/
```
Create the working directory for this practical:
```bash
mkdir gwas_gemma
```
Copy scripts:
```bash
cp -r /usr/local/extras/Genomics/workshops/January2018/gwas_gemma/scripts ./gwas_gemma
```
Copy data:
```bash
cp -r /usr/local/extras/Genomics/workshops/January2018/gwas_gemma/data ./gwas_gemma/
```
To check you are getting the expected results, you can also copy a directory containing all the files that should be produced after running this practical:
```bash
cp -r /usr/local/extras/Genomics/workshops/January2018/gwas_gemma/results ./gwas_gemma/
```
Change to the working directory we are going to use for this practical:
```bash
cd gwas_gemma
```

## 2. Data formatting
Now let's have a look at the data.
You should have the following input files:
```bash
ls -lh data
```
>total 1.2G<br>
>-rw-r--r-- 1 myuser cs 347M Jan 28 18:41 fha.bbgeno.gz<br>
>-rw-r--r-- 1 myuser cs 6.5K Jan 28 18:40 fha.bbgeno.ids.txt<br>
>-rw-r--r-- 1 myuser cs 9.9K Jan 28 18:40 fha.pheno<br>
>-rw-r--r-- 1 myuser cs 1.3K Jan 28 18:40 fha.pheno2<br>
>-rw-r--r-- 1 myuser cs 876M Jan 28 18:41 fha.vcf.gz<br>

There is a vcf file containing single nucleotide polymorphisms (SNPs) from RAD data of 602 individuals from a single polymorphic population of *Timema cristinae* (population code FHA). vcf is a very popular format for genetic variants, you can find more info [here](http://www.internationalgenome.org/wiki/Analysis/vcf4.0/). You can have a look at the file content with the following commands:
```bash
gzip -dc data/fha.vcf.gz | less -S
# or with bcftools
bcftools view data/fha.vcf.gz | less -S
# excluding long header
bcftools view -H data/fha.vcf.gz | less -S
```
There are two files containing the phenotypes of the same inviduals in the same order than in the vcf file (NA when the phenotype is missing). This phenotypes encode the dorsal white stripe either as a continuous trait (standardised white area; ```fha.pheno```) or as a discerte binary trait (presence/absence; ```fha.pheno2```). You can have a look at the files content:
```bash
head data/fha.pheno
```
>0.866078916198945<br>
>-1.17516992488642<br>
>NA<br>
>NA<br>
>NA<br>
>-3.11627693813939<br>
>-0.348977465694598<br>
>-0.663007099264362<br>
>NA<br>
>NA<br>

```bash
head data/fha.pheno2
```
>1<br>
>1<br>
>NA<br>
>NA<br>
>NA<br>
>0<br>
>1<br>
>1<br>
>NA<br>
>NA<br>

```gemma```, the program we will be using to carry out multi-variant GWA, needs two input files: one with the genotypes and another one with the phenotypes. Most GWA programs rely on called genotypes, but ```gemma``` can work with genotype probabilities, which allows incorporating genotype uncertainty in the analyses. ```gemma``` accepts a number of formats for genotypes, we are going to use the mean genotype format (based on the BIMBAM format), where genotypes are encoded as mean genotypes. A mean genotype is a value between 0 and 2 that can be interpreted as the minor allele dosage: 0 is homozygous for the major allele, 1 is a heterozygote, and 2 is a homozygote for the minor allele. Thus, intermediate values reflect the uncertainty in genotypes. 

We are going to use the custom Perl script ```bcf2bbgeno.pl``` to get such mean genotypes. First, the script removes all the SNPs that have more than two alleles. This is common practice, because most of the models used for inference are based in biallelic variants. Then, the script calculates empirical posterior genotype probabilities from the genotype likelihoods in the vcf file under the assumption that the population is in Hardy-Weinberg equilibrium (HWE). Specifically, the script uses inferred allele frequencies to set HWE priors:

*p*(AA) = *p*<sup>2</sup>; *p*(aa) = (1-*p*)<sup>2</sup>; *p*(Aa) = 2*p*(1-*p*)

being *p* the allele frequency of major/reference allele A. Genotype likelihoods are multiplied by these priors to obtain genotype posterior probabilities that are then encoded as mean genotypes and saved to a .bbgeno file.

You can get some info about how to run the Perl script:
```bash
# show help
scripts/bcf2bbgeno.pl -h
```

Now let's calculate the genotype posterior probabilites and save them in mean genotype format. This process may take quite a while (over 15 minutes), so you may want to either put it a bash script and submit a batch job to the cluster (see ```bcf2bbgeno.sh```) or skip this step altogether (see below).

```bash
# Execute script
scripts/bcf2bbgeno.pl -i data/fha.vcf.gz -o fha.bbgeno -p H-W -s -r
```
When the script finishes, given that ```gemma``` can handle gzipped files, it is a good idea to compress the output file to save some space:
```bash
# compress the file to save some space
gzip fha.bbgeno
```

As stated before, this step may take a while and it may be worth to simply skip it and use the .bbgeno.gz and .bbgeno.ids.txt files from the results folder:
```bash
cp results/*.bbgeno.* ./
```

In any case, you should have now two new files: fha.bbgeno.gz, which contains the mean genotypes, and fha.bbgeno.ids.txt, which contains the ids of the samples (i.e. individuals) in the same order as shown in the genotypes file. You can have a look at them:

```bash
gzip -dc fha.bbgeno.gz | less -S
```
>lg13_ord45_scaf428-158031 C T 0.01108 0.00279 0.00140 0.00140 0.01108 0.00018 0.01108 0.00555<br>
>lg13_ord45_scaf428-48027 T C 0.00299 0.00597 0.00150 0.00019 0.00597 0.00038 0.00009 0.00075<br>
>lg13_ord45_scaf428-80879 G A 0.00163 0.00163 0.00325 0.00041 0.00325 0.00041 0.00082 0.00163<br>
>lg13_ord45_scaf428-94107 G A 0.00358 0.00358 0.02771 0.01408 0.05530 0.00045 0.01408 0.00358<br>
>lg13_ord45_scaf428-158069 T A 0.16220 0.04246 0.02174 0.99552 0.16220 0.00279 0.16220 0.08128<br>
>lg13_ord45_scaf428-325672 A T 0.00784 0.00393 0.00197 0.00197 0.00784 0.00099 0.00784 0.00025<br>

The first column is the SNP id, second one is major/reference allele and the third one is the minor/alternate allele. Having major and minor or reference and alternate alleles depends on how you called variants. This does not matter for downstream analyses with ```gemma```, but should be born in mind when interpreting the results. The rest of the columns are the mean genotypes for all the individuals (one for each).

```bash
less -S fha.ids.txt
```
>2013FHA001<br>
>2013FHA002<br>
>2013FHA004<br>
>2013FHA005<br>
>2013FHA006<br>
>2013FHA007<br>

## 3. Running ```gemma``` to fit a Bayesian sparse linear mixed model (BSLMM)

```gemma```, standing for Genome-wide Efficient Mixed Model Association, is a complex piece of software with many options. It is highly recommended that you read carefully the [manual](http://www.xzlab.org/software/GEMMAmanual.pdf). Let's have a look at the options, with emphasis on the specific ones we will be using.
```bash
# general help
gemma -h
# quick guide
gemma -h 1
# input files
gemma -h 2
```
We will be using the Bayesian sparse linear mixed model (BSLMM):
```bash
gemma -h 9
```
but first we will have to calculate the relatedness matrix to take into account population structure:
```bash
gemma -h 8
```
and it would also be interesting do some filtering such as excluding rare variants, which are difficult to tell apart from sequencing errors:
```bash
gemma -h 3
```

```gemma``` uses a relatedness matrix to account for any population structure in the data that may affect results. The easiest way to obtain such matrix is using ```gemma``` itself. We are going to use the submission script ```gemma_relmatrix.sh```. Let's have a look with the nano text editor:
```bash
nano scripts/gemma_relmatrix.sh
```bash
#!/bin/bash
#$ -l h_rt=1:00:00
#$ -j y
#$ -o gemma_relmatrix.log

GEMMA='gemma'
DIR="/data/$USER/gwas_gemma"
GENOTYPES='data/fha.bbgeno.gz'
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
```

And then let's submit the job (it should take just a few minutes):
```bash
qsub scripts/gemma_relmatrix.sh
```
You can check the status of your jobs with ```Qstat```. When finished, it should have two files in the newly created output directory: relmatrix.log.txt, a log file, and relmatrix.cXX.txt, the actual relatedness matrix. If you have a look at them, they should look like these:
```bash
cat output/relmatrix.log.txt
```
>##<br>
>## GEMMA Version = 0.94<br>
>##<br>
>## Command Line Input = -g data/fha.bbgeno.gz -p data/fha.pheno -gk 1 -o relmatrix <br>
>##<br>
>## Summary Statistics:<br>
>## number of total individuals = 602<br>
>## number of analyzed individuals = 546<br>
>## number of covariates = 1<br>
>## number of phenotypes = 1<br>
>## number of total SNPs = 518232<br>
>## number of analyzed SNPs = 346660<br>
>##<br>
>## Computation Time:<br>
>## total computation time = 3.91983 min <br>
>## computation time break down: <br>
>##      time on calculating relatedness matrix = 2.46767 min <br>
>##<br>
><br>

```bash
less -S output/relmatrix.cXX.txt
```
>0.08650617109   0.000333571409  0.001472054408  0.001656484245  0.001901346128  0.001199397691  0.001921586938  0.00<br>
>0.000333571409  0.08827073595   0.0008395834056 -0.0008735407094        0.0004308437331 0.004303127573  0.0014308843<br>
>0.001472054408  0.0008395834056 0.09971276805   0.003820448145  0.001814713502  -0.0001891657497        0.0022152800<br>
>0.001656484245  -0.0008735407094        0.003820448145  0.10577918      0.002068707695  -7.400802932e-05        0.00<br>
>0.001901346128  0.0004308437331 0.001814713502  0.002068707695  0.05202678922   -8.732093775e-05        0.0021426359<br>
>0.001199397691  0.004303127573  -0.0001891657497        -7.400802932e-05        -8.732093775e-05        0.1186797893<br>

Now we are going to run ```GEMMA``` to fit a BSLMM model. Since this will take quite some time (~15 minutes), we are going to submit it as batch job to the cluster queue. We will be using the script ```gemma_bslmm.sh```. Let's have a look with nano:
```bash
nano scripts/gemma_bslmm.sh
```
```bash
#!/bin/bash
#$ -l h_rt=07:00:00
#$ -l rmem=2g
#$ -l mem=4g
#$ -j y
#$ -o gemma_bslmm.log

GEMMA='gemma'

DIR="/data/$USER/gwas_gemma"
GENOTYPES='data/fha.bbgeno.gz'
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
```

## 4. Analysing ```gemma``` BSLMM output
