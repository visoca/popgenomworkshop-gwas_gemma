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

>total 2.5G<br>
>-rw-r--r-- 1 myuser cs  22K Jan 26 12:31 lg_ord_sca_length.dsv<br>
>-rw-r--r-- 1 myuser cs 768M Jan 26 12:31 timemaHVA.gl<br>
>-rw-r--r-- 1 myuser cs 467M Jan 26 12:31 timemaHVA.vcf.gz<br>
>-rw-r--r-- 1 myuser cs 781M Jan 26 12:31 timemaHVC.gl<br>
>-rw-r--r-- 1 myuser cs 476M Jan 26 12:31 timemaHVC.vcf.gz<br>

There is a vcf file containing single nucleotide polymorphisms (SNPs) from RAD data of 605 individuals from a single polymorphic population of *Timema cristinae* (population code FHA). vcf is a very popular format for genetic variants, you can find more info [here](http://www.internationalgenome.org/wiki/Analysis/vcf4.0/). You can have a look at the file content with the following commands:
```bash
gzip -dc fha.vcf.gz | less -S
# or with bcftools
bcftools view fha.vcf.gz | less -S
# excluding long header
bcftools view -H fha.vcf.gz | less -S
```
This file has to be converted to...
