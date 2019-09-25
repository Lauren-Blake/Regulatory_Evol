#!/usr/bin/perl
use strict;
use warnings;

## Julien Roux
## Sep 22, 2014
# - Run on spudhead
# - Launches overlap_features_DMRs_0.2.pl for all DMRs, and all the random DMRs files of same length / CpG density

scalar @ARGV >= 1 or die<<USAGE;
Usage:
perl launch_all_overlap_features_DMRs_0.2_control_CpG_density.pl <DMR type>
USAGE
my $DMRtype = $ARGV[0];
print "Launching overlap of genomic features for $DMRtype DMRs\n";

# first create BED files of DMR positions
my @files = glob("../bsseq/DMRs/$DMRtype/*_DMRs.txt");

foreach my $file (@files){
  print $file, "\n";
  $file =~ m/\.\.\/bsseq\/DMRs\/$DMRtype\/(.+\_DMRs)\.txt/;
  my $DMRfile = $1;

#  ## only relaunch some of the jobs
#  next unless ($DMRtype eq 'individuals');
#  $DMRfile =~ m/(.+)\_\d\_\d\_DMRs/;
#  next unless ($1 eq 'Rhesus');

  # launch script overlap_features_DMRs_0.2.pl
  # already done by launch_all_overlap_features_DMRs_0.2.pl for species and tissues
  # system("echo \"perl overlap_features_DMRs_0.2.pl ./DMRs/$DMRtype/$DMRfile.bed\" | qsub -l h_vmem=4g -v PATH -cwd -N $DMRfile -o ./DMRs/$DMRtype/$DMRfile.out -e ./DMRs/$DMRtype/$DMRfile.err");

  # launch script overlap_features_DMRs_0.2.pl for 100 sets of random DMRs of same length/CpG density
  my $randomDMRfile = $DMRfile;
  $randomDMRfile=~ s/DMRs/randomDMRs_/;

  for (my $i=1; $i <= 100; $i++) {
    system("echo \"perl overlap_features_DMRs_0.2.pl ./DMRs/$DMRtype/randomized_control_CpG_density/$randomDMRfile$i.bed\" | qsub -l h_vmem=4g -v PATH -cwd -N $randomDMRfile$i -o ./DMRs/$DMRtype/randomized_control_CpG_density/$randomDMRfile$i.out -e ./DMRs/$DMRtype/randomized_control_CpG_density/$randomDMRfile$i.err");
  }
  #last; # testing
}


