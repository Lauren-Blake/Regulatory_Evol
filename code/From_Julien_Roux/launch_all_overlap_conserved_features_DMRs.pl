#!/usr/bin/perl
use strict;
use warnings;

## Julien Roux
## Mar 11, 2016
# - Run on spudhead
# - Launches overlap_conserved_features_DMRs.pl for all DMRs

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

  # launch script overlap_conserved_features_DMRs.pl (species and tissues): 
  unless ($DMRtype eq 'individuals'){
    system("echo \"perl overlap_conserved_features_DMRs.pl ./DMRs/$DMRtype/$DMRfile.bed\" | qsub -l h_vmem=4g -v PATH -cwd -N $DMRfile -o ./DMRs/$DMRtype/$DMRfile.out -e ./DMRs/$DMRtype/$DMRfile.err");
  }

  #   # launch script overlap_conserved_features_DMRs.pl for 100 sets of random DMRs of same length
  #   my $randomDMRfile = $DMRfile;
  #   $randomDMRfile=~ s/DMRs/randomDMRs_/;
  
  #   for (my $i=1; $i <= 100; $i++) {
  #     system("echo \"perl overlap_conserved_features_DMRs.pl ./DMRs/$DMRtype/randomized/$randomDMRfile$i.bed\" | qsub -l h_vmem=4g -v PATH -cwd -N $randomDMRfile$i -o ./DMRs/$DMRtype/randomized/$randomDMRfile$i.out -e ./DMRs/$DMRtype/randomized/$randomDMRfile$i.err");
  #   }
  #last; # testing
}


