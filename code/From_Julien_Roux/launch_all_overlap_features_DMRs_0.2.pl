#!/usr/bin/perl
use strict;
use warnings;

## Julien Roux
## Oct 10, 2014
# - Run on spudhead
# - Launches overlap_features_DMRs_0.2.pl for all DMRs, and all the random DMRs files of same length (no control for CpG density)

# first create BED files of DMR positions
my @files = glob("../bsseq/DMRs/*/*_DMRs.txt");

foreach my $file (@files){
  print $file, "\n";
  $file =~ m/\.\.\/bsseq\/DMRs\/(.+)\/(.+\_DMRs)\.txt/;
  my $DMRtype = $1;
  my $DMRfile = $2;

  # ## only relaunch some of the jobs, e.g Rhesus individual DMRs
  # next unless ($DMRtype eq 'individuals');
  # $DMRfile =~ m/(.+)\_\d\_\d\_DMRs/;
  # next unless ($1 eq 'Rhesus');

  # launch script overlap_features_DMRs_0.2.pl: 
  # for individualDMRs: already done by launch_all_overlap_features_DMRs_0.2_control_CpG_density.pl
  unless ($DMRtype eq 'individuals'){
    system("echo \"perl overlap_features_DMRs_0.2.pl ./DMRs/$DMRtype/$DMRfile.bed\" | qsub -l h_vmem=4g -v PATH -cwd -N $DMRfile -o ./DMRs/$DMRtype/$DMRfile.out -e ./DMRs/$DMRtype/$DMRfile.err");
  }

  # launch script overlap_features_DMRs_0.2.pl for 100 sets of random DMRs of same length
  my $randomDMRfile = $DMRfile;
  $randomDMRfile=~ s/DMRs/randomDMRs_/;

  for (my $i=1; $i <= 100; $i++) {
    system("echo \"perl overlap_features_DMRs_0.2.pl ./DMRs/$DMRtype/randomized/$randomDMRfile$i.bed\" | qsub -l h_vmem=4g -v PATH -cwd -N $randomDMRfile$i -o ./DMRs/$DMRtype/randomized/$randomDMRfile$i.out -e ./DMRs/$DMRtype/randomized/$randomDMRfile$i.err");
  }
  #last; # testing
}


