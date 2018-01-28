*Population Genomics Workshop 2018, University of Sheffield*

# Genetic architecture of traits using Genome Wide Association (GWA) mapping with GEMMA
#### Victor Soria-Carrasco

The aim of this practical is ....

Data and scripts used during this practical will be available in a shared folder in Iceberg, the University of Sheffield HPC cluster, which will allow faster transfers. However, in case attendees would like to use it for practice after the workshop, scripts are also available in [this repository](https://github.com/visoca/popgenomworkshop-gwas_gemma) and data can be downloaded [here]() (and results [here]()).

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
gzip -dc fha.vcf.gz | less -S
# or with bcftools
bcftools view fha.vcf.gz | less -S
# excluding long header
bcftools view -H fha.vcf.gz | less -S
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

```gemma```, the program we will be using to carry out multi-variant GWA, needs two input files: one with the genotypes and another one with the phenotypes. Most GWA programs rely on called genotypes, but ```gemma``` can work with genotype probabilities, which allows incorporating genotype uncertainty in the analyses. ```gemma``` accepts a number of formats for genotypes, we are going to use the mean genotype format (based on the BIMBAM format), where genotypes are encoded as posterior mean genotype probabilities. A posterior mean genotype is a value between 0 to 2 that can be interpreted as the minor allele dosage: 0 is homozygous for the major allele, 1 is a heterozygote, and 2 is a homozygote for the minor allele.
We are going to use the custom Perl script ```bcf2bbgeno.pl``` to calculate empirical mean genotype posterior probabilities from the genotype likelihoods in the VCF. We will use inferred allele frequencies to set Hardy-Weinberg Equilibrium priors (i.e. *p*(AA) = *p*<sup>2</sup>; *p*(aa) = (1-*p*)<sup>2</sup>; *p*(Aa) = 2*p*(1-*p*); being *p* the allele frequency of major/reference allele A).

You can some info about how to run the Perl script:
```bash
# show help
scripts/bcf2bbgeno.pl -h
```

Now let's calculate the mean genotypes. This may take a while.
```bash
# Execute script
scripts/bcf2bbgeno.pl -i data/fha.vcf.gz -o fha.bbgeno -p H-W -s -r
```
And then compress the output file to save some space
```bash
# then compress the file to save some space
$ gzip fha.bbgeno
```
