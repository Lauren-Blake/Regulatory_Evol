#!/usr/bin/perl
use strict;
use warnings;

## Julien Roux
## Oct 8, 2015
# - Run on spudhead
# - Launches bedtools to calculate CpG density of DMRs for all DMRs, and all the random DMRs files of same length (no control for CpG density!)

# first get files of DMR positions
my @files = glob("./DMRs/*/*_DMRs.txt");

foreach my $file (@files){
  print $file, "\n";
  $file =~ m/\.\/DMRs\/(.+)\/(.+\_DMRs)\.txt/;
  my $DMRtype = $1;
  my $DMRfile = $2;

  # Individual DMRs not launched (not in paper)
  next if ($DMRtype eq 'individuals');

  # Only species DMRs launched this time
  next if ($DMRtype eq 'tissues');

  # Some DMR files in ./bsseq should not be taken into account
  next if ($DMRfile =~ m/\_overlap\_/);

  # ## used to relaunch only some of the jobs, e.g Rhesus individual DMRs
  # $DMRfile =~ m/(.+)\_\d\_\d\_DMRs/;
  # next unless ($1 eq 'Rhesus');
  # next unless ($DMRfile =~ m/Specific/);

  ## File with all the sets of CpG used for each contrast
  my $CpG_file = $DMRfile;
  $CpG_file=~ s/\_DMRs/\_tstat_allSites\.bed/;

  my $outfile = $DMRfile.'_number_CpG_sites.bed';

  # launch bedtools for DMR file
  system("echo \"coverageBed -a DMRs/$DMRtype/$CpG_file -b ../annotation/DMRs/$DMRtype/$DMRfile.bed -counts > DMRs/$DMRtype/number_CpG_sites/$outfile\" | qsub -l h_vmem=1g -v PATH -cwd -N $DMRfile -o ./DMRs/$DMRtype/number_CpG_sites/$DMRfile.out -e ./DMRs/$DMRtype/number_CpG_sites/$DMRfile.err");


  # launch bedtools for 100 sets of random DMRs of same length
  my $randomDMRfile = $DMRfile;
  $randomDMRfile=~ s/DMRs/randomDMRs_/;

  for (my $i=1; $i <= 100; $i++) {
    $outfile = $randomDMRfile.$i.'_number_CpG_sites.bed';
    system("echo \"coverageBed -a DMRs/$DMRtype/$CpG_file -b ../annotation/DMRs/$DMRtype/randomized/$randomDMRfile$i.bed -counts > DMRs/$DMRtype/number_CpG_sites/$outfile\" | qsub -l h_vmem=1g -v PATH -cwd -N $randomDMRfile$i -o ./DMRs/$DMRtype/number_CpG_sites/$randomDMRfile$i.out -e ./DMRs/$DMRtype/number_CpG_sites/$randomDMRfile$i.err");
  }
  # last; # testing
}


