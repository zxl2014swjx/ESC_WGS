#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my %opts;

GetOptions(\%opts,"c=s","n=s","s=s","o=s","f=s","h");

my $usage =<< "USAGE";
Program: $0 -c cancer.mpileup -n normal.mpileup -s Gender -o OutputDir -f config_file
            -c mpileup file for cancer sample
            -n mpileup file for normal sample
            -s Gender, 'XY' for man and 'XX' for woman
            -o output directory
            -f config file
            -h help
USAGE

die $usage if(!$opts{'c'} || !$opts{'n'} || !$opts{'s'} || !$opts{'o'} || !$opts{'f'} || $opts{'h'});

my $cancer_mpileup = $opts{'c'};
my $normal_mpileup = $opts{'n'};
my $output_dir = $opts{'o'};
my $config_file = $opts{'f'};
my $gender = $opts{'s'};

open(OUT, ">$config_file");
print OUT "[general]\n\n";
print OUT "chrLenFile = /mnt/X500/farmers/suyao/NL200/data/reference/hs37d5.len\n";
print OUT "window = 500000\n";
print OUT "step = 100000\n";
print OUT "ploidy = 2\n";
print OUT "outputDir = $output_dir\n";
print OUT "chrFiles = /mnt/X500/farmers/suyao/NL200/data/reference/chromosome\n";
print OUT "sex = $gender\n";
print OUT "breakPointType = 2\n";
print OUT "maxThreads = 10\n";
print OUT "breakPointThreshold = 0.8\n";
print OUT "noisyData = FALSE\n";
print OUT "printNA = FALSE\n";
print OUT "readCountThreshold = 10\n";
print OUT "[sample]\n\n";
print OUT "mateFile = $cancer_mpileup\n";
print OUT "inputFormat = pileup\n";
print OUT "mateOrientation = FR\n";
print OUT "[control]\n\n";
print OUT "mateFile = $normal_mpileup\n";
print OUT "inputFormat = pileup\n";
print OUT "mateOrientation = FR\n";
print OUT "[BAF]\n\n";
print OUT "SNPfile = /mnt/X500/farmers/suyao/NL200/tools/BNC/BNC/program/NoahCare/db/alignment/gatk_bundle/dbsnp_138.b37.del100.vcf.gz\n";
print OUT "minimalCoveragePerPosition = 5\n";

close(OUT);












