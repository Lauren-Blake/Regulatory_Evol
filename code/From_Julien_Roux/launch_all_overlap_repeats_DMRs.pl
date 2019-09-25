#!/usr/bin/perl
use strict;
use warnings;

## Julien Roux
## Aug 20, 2015
# - Run on spudhead
# - Launches overlap_repeats_DMRs.pl for all DMRs, and all the random DMRs files of same length (no control for CpG density)

# first create BED files of DMR positions
my @files = glob("../bsseq/DMRs/*/*_DMRs.txt");

foreach my $file (@files){
  print $file, "\n";
  $file =~ m/\.\.\/bsseq\/DMRs\/(.+)\/(.+\_DMRs)\.txt/;
  my $DMRtype = $1;
  my $DMRfile = $2;

  # Individual DMRs not launched (not in paper)
  next if ($DMRtype eq 'individuals');

  # Some DMR files in ../bsseq should not be taken into account
  next if ($DMRfile =~ m/\_overlap\_/);

  # ## only relaunch some of the jobs, e.g Rhesus individual DMRs
  # next unless ($DMRtype eq 'individuals');
  # $DMRfile =~ m/(.+)\_\d\_\d\_DMRs/;
  # next unless ($1 eq 'Rhesus');

  # launch script overlap_repeats_DMRs.pl: 
  system("echo \"perl overlap_repeats_DMRs.pl ./DMRs/$DMRtype/$DMRfile.bed\" | qsub -l h_vmem=4g -v PATH -cwd -N $DMRfile -o ./DMRs/$DMRtype/$DMRfile.out -e ./DMRs/$DMRtype/$DMRfile.err");

  # launch script for 100 sets of random DMRs of same length
  my $randomDMRfile = $DMRfile;
  $randomDMRfile=~ s/DMRs/randomDMRs_/;

  for (my $i=1; $i <= 100; $i++) {
    system("echo \"perl overlap_repeats_DMRs.pl ./DMRs/$DMRtype/randomized/$randomDMRfile$i.bed\" | qsub -l h_vmem=4g -v PATH -cwd -N $randomDMRfile$i -o ./DMRs/$DMRtype/randomized/$randomDMRfile$i.out -e ./DMRs/$DMRtype/randomized/$randomDMRfile$i.err");
  }
  #last; # testing
}


