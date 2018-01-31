#!/usr/bin/env perl

# (c) Victor Soria-Carrasco
# Last modified: 31/01/2018 01:01:15

# Description:
# This script estimates genotype probabilities and 
# converts a vcf/bcf file to BIMBAM mean genotype format 
# (i.e. alternate allele dosage)
# 
# Mean genotype probabilities can be estimated using a uniform or
# a Hardy-Weinberg prior
# Be aware that SNPs with >2 alleles will be excluded

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use File::Path qw(make_path);

my $version='1.3-2018.01.31';

# Changelog
# 03/03/2016
#	- Fixed automatic output filename
#	- Added check for existence of input file
#	- Added support for bcftools 1.2.x
#	- Added support for .vcf.gz and .vcf.bz2 files
#	- If AF is not available, it now calculates  allele frequencies from 
#	  allele counts (AC/AN) - but throws warnings

&author;

my $bcftools='bcftools';
# my $bcftools='/usr/local/extras/Genomics/apps/bcftools/1.2/bin/bcftools';

my $bcftools_version=`$bcftools |& grep Version | awk '{print \$2}'`;


# Read arguments and set up output directories/files
# =============================================================================
my $prior='H-W';
print "\n";
my $samples=0;
my $inv=0;
my ($infile, $outfile);
GetOptions(
    'i|I=s'  => \$infile,
    'o|O=s'  => \$outfile,
    'p|P=s'  => \$prior,  
    's|S'    => \$samples,
	'r|R'    => \$inv,
    'h|help' => \&usage
)
or (print "\nERROR: Argument invalid\n" and &usage);

&usage if (!defined($infile) || ($prior ne 'uniform' && $prior ne 'H-W'));

$infile=File::Spec->rel2abs($infile);

# Check existence input file
die ("\nFile $infile does not exist\n\n")
	if (!-e $infile);

# Determine if format is vcf or bcf
my $input=$infile;
if ($input =~ m/\.bcf$/i){
	if ($bcftools_version=~ m/^0\.1/){ # bcftools 0.1.x
		$input="$bcftools view -N $infile |";
	}
	else{ # bcftools v 1.2
		$input="$bcftools view -e 'REF=\"N\"' $infile |";
	}
}
elsif ($input =~ m/\.vcf\.gz$/i){
	$input="gzip -cd $infile |";
}
elsif ($input =~ m/\.vcf\.bz2$/i){
	$input="lbzip2 -cdk $infile |";
}
elsif ($input =~ m/\.vcf$/i){
	$input="cat $infile |";
}
else{
	die ("\nFormat of input file $infile not recognised\n\n");
}

# Define output file, 
# same name that input file, but bb extension if not provided
if (!defined($outfile)){
	$outfile=$infile;
	$outfile=~ s/\.[bv]cf(\.gz|\.bz2)?$/\.bbgeno/;
}
else{
	$outfile=File::Spec->rel2abs($outfile);
}

my $warnfile="$outfile.warnings";

# =============================================================================

# Read vcf/bcf file
# =============================================================================
my @warnings;
my %genotypes;
my @ids=();
my $nsnps=0;
open (FILE, $input) or
	die ("\nCan't open input file $infile\n\n");
while (<FILE>){
	next if (/^\#\#/);
	
	chomp;
	my @aux=split(/\s/,$_);
	# Header
	if (/^\#CHROM/){
		# print STDERR "Header found\n";
		foreach my $i (9..$#aux){
			$aux[$i]=~ s/(\.sorted)?(\.[s|b]am)?//g;
			push (@ids, $aux[$i]);
		}
	}
	else{
		print "SNPs processed: $nsnps\r" if ($nsnps%100==0);
		$|++; # flush buffer

		if ($aux[4]!~ m/\,/ && # exclude multiallelic SNPs (> 2 alleles) 
			$aux[3] ne 'N'){   # exclude SNPs where ref = N

			# print "LINE: $_\n";
			my @line=();
			push(@line, $aux[0]."_".$aux[1]); # Add SNP id as "scaffold_position"
			# Get allele frequencies (there could be multiple estimates if bcf/vcf file was filtered after calling SNPs)
			my @afs=();
			while ($aux[7]=~ m/AF(1)?\=([0-1]\.?[0-9]*e?\-?[0-9]*)\;/g){
				push (@afs, $2);
			}
			
			# If no allele frequency is found, estimate it from allele counts
			if (!defined($afs[0])){
				push (@warnings, "Warning: No allele frequency found for $aux[0]:$aux[1]. It will be estimated from allele counts (AC/AN)");
				$aux[7]=~ m/AC\=([0-9]+)/;
				my $ac=$1;
				$aux[7]=~ m/AN\=([0-9]+)/;
				my $an=$1;
				# print "$aux[0]:$aux[1] INFO $aux[7] - AC $ac - AN $an - AF ".($ac/$an)."\n";
				push (@afs, $ac/$an);
			}
			
			# exclude private SNPs: all individuals have the same allele, but it is different from reference
			next if ($inv && $afs[$#afs] == 1);	

			# Use last estimated allele frequency
			my $af=$afs[$#afs];
			# Get alleles states
			my $alleles=" $aux[4] $aux[3]";						
			if ($af > 0.5){ # reference allele is actually minor allele
				$alleles=" $aux[3] $aux[4]";
			}
			$genotypes{$aux[0]}{$aux[1]}.=$alleles;

			# get column index for genotype likelihoods from FORMAT field	
			# print "CHECK $aux[8]\n";
			my @aux2=split(/\:/,$aux[8]);
			my $plcol=0;
			while ($aux2[$plcol] ne "PL"){
				$plcol++;
			}
			
			# Get mean genotype probabilities
			foreach my $i (9..$#aux){

				# Convert genotype Phred-scaled likelihoods to composite genotype probability
				my $cgp;
				if ($prior eq 'uniform'){
					$cgp=pl2cgp_uni($aux[$i],$af,$plcol); 
				}
				elsif ($prior eq 'H-W'){
					$cgp=pl2cgp_hw($aux[$i],$af,$plcol);
				}
				$genotypes{$aux[0]}{$aux[1]}.=" ".sprintf("%.5f",$cgp);
			}
			$nsnps++;
		}
	}
	# last if ($nsnps > 1000); # testing
}
close (FILE);
print "\n\n";
# =============================================================================

# Write output to file
# =============================================================================
# Otuput genotypes
my $size;
open (OUT, ">$outfile")
	or die ("\nCan't write to output file: $outfile\n\n");
foreach my $sca (keys %genotypes){
	foreach my $pos (keys %{$genotypes{$sca}}){
		my $thissize=()=$genotypes{$sca}{$pos}=~ / /gi;
		die ("\nSomething went wrong. Different number of samples for different SNPs: $size != ".$thissize."\n\n")
			if (defined($size) && $thissize != $size);
		$size=$thissize;
		
		print OUT "$sca-$pos".$genotypes{$sca}{$pos}."\n";
	}
}
close (OUT);

print "\nFormat conversion finished. File saved as: $outfile\n\n";

# Output sample order
if ($samples){
	my $outsamples=$outfile.'.ids.txt';
	open (OUT, ">$outsamples")
		or die ("\nCan't write to output file: $outsamples\n\n");
		foreach my $i (@ids){
			print OUT "$i\n";
		}
	close (OUT);
	print "\nSamples (ids) order saved as: $outsamples\n\n";
}

# Output warnings
if (scalar(@warnings)>0){
	open (OUT, ">$warnfile")
		or die ("\nCan't write to warnings file: $warnfile\n\n");
		print OUT join ("\n",@warnings);
	close (OUT);
	print "\nBE AWARE: Warnings saved as: $warnfile\n\n";
}
# =============================================================================


# ==============================================================================
# ==============================================================================
# ============================== SUBROUTINES ===================================
# ==============================================================================
# ==============================================================================



# Convert Phred-scaled genotype likelihoods (pl) to
# composite genotype probabilities (cgp)
# ==============================================================================

# using uniform priors 
# -----------------------------------------------------------------------------
sub pl2cgp_uni{
	my $gt=shift;
	my $af=shift;
	my $plcol=shift;

	my @aux=split(/\:/,$gt);
	$gt=$aux[0];
	my $pl=$aux[$plcol];

	my @pls=(1/3,1/3,1/3);;
	if ($gt ne './.'){
		@pls=split(/\,/,$pl);
	}	
	
	# if ref allele is not the major allele
	# likelihoods must be swapped
	# we don't need to change AF because it
	# is not used for the computation of probabilities
	if ($af>0.5){
		my @aux=@pls;
		$pls[0]=$aux[2];
		$pls[2]=$aux[0];
	}

	my $sum=0;
	my @gps=();
	foreach my $pl (@pls){
		my $gl=10**(-$pl/10);
		push (@gps, $gl);
		$sum+=$gl;
	}
	foreach my $pr (@gps){
		$pr=$pr/$sum;
	}
	my $cgp=$gps[2]*2+$gps[1]; # minor allele dosage

	return($cgp);
}
# -----------------------------------------------------------------------------

# using Hardy-Weinberg priors
# -----------------------------------------------------------------------------
sub pl2cgp_hw{
	my $gt=shift;
	my $af=shift;
	my $plcol=shift;
	
	my @aux=split(/\:/,$gt);
	$gt=$aux[0];
	my $pl=$aux[$plcol];

	my @pls=(1/3,1/3,1/3);;
	if ($gt ne './.'){
		@pls=split(/\,/,$pl);
	}

	# if ref allele is not the major allele
	# likelihoods must be swapped
	# minor allele frequency is then 1-af
	if ($af>0.5){
		my @aux=@pls;
		$pls[0]=$aux[2];
		$pls[2]=$aux[0];
		$af=1-$af;
	}

	my @gps=();
	# 1st term => Phred-descalation
	# 2nd term => Hardy-Weinberg prior
	$gps[0]=(10**(-$pls[0]/10)) * ((1-$af)**2);
	$gps[1]=(10**(-$pls[1]/10)) * (2*$af*(1-$af));
	$gps[2]=(10**(-$pls[2]/10)) * ($af**2);

	my $sum=0;
	foreach my $pr (@gps){
		$sum+=$pr
	}
	if ($sum==0){
		print "PL: $pl - AF: $af - @gps\n";	
	}
	foreach my $pr (@gps){
		$pr=$pr/$sum;
	}
	my $cgp=$gps[2]*2+$gps[1]; # minor allele dosage 

	return($cgp);
}
# -----------------------------------------------------------------------------
# ==============================================================================

# Show copyright
# ==============================================================================
sub author{
    print "\n";
	print "#########################################\n";
	print "  ".basename($0)." version $version     \n";
    print "  (c) Victor Soria-Carrasco             \n";
    print "  victor.soria.carrasco\@gmail.com      \n";
	print "#########################################\n";
}
# ==============================================================================

# Show usage
# ==============================================================================
sub usage{
	print "\n";
    print "  Usage:\n";
    print "    ".basename($0)."\n";
    print "      -i <SNPs input file (bcf/vcf format)>\n";
    print "      -o <output file> (optional)\n";
	print "      -p <prior for genotype probabilities> (uniform | H-W ; optional, default: H-W)\n";
	print "      -s (output sample order to a file)\n";
	print "      -r (remove private SNPs for this dataset (i.e. same allele in all individuals, but different from reference))\n";
    print "\n";
	print "  Example:\n";
	print "      ".basename($0)." -i variants.bcf\n";
    print "\n\n";
    exit;
}

