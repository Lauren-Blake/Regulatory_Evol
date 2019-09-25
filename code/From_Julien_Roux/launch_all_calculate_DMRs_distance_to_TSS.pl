#!/usr/bin/perl
use strict;
use warnings;

## Julien Roux
## Oct 15, 2015
# - Run on spudhead
# - Launches closestBed to calculate distance to TSS for all DMRs, and all the random DMRs files of same length + control for CpG density

# first get all files of DMR positions
my @files = glob("../bsseq/DMRs/*/*_DMRs.txt");

foreach my $file (@files){
  print $file, "\n";
  $file =~ m/DMRs\/(.+)\/(.+\_DMRs)\.txt/;
  my $DMRtype = $1;
  my $DMRfile = $2;

  # Individual DMRs not launched (not in paper)
  next if ($DMRtype eq 'individuals');

  # Already done
  next if ($DMRtype eq 'tissues');

  # Some DMR files in ./bsseq should not be taken into account
  next if ($DMRfile =~ m/\_overlap\_/);

  # ## used to relaunch only some of the jobs, e.g Rhesus individual DMRs
  # $DMRfile =~ m/(.+)\_\d\_\d\_DMRs/;
  # next unless ($1 eq 'Rhesus');
  # next unless ($DMRfile =~ m/Specific/);
  # next if ($DMRfile =~ m/Human_heart_liver/); # already done before

  # in the ~/Methylation/annotation/DMRs/ folder
  my $outfile = $DMRfile.'_distanceToTSS.bed';



  # launch bedtools for DMR file
  system("echo \"closestBed -a ./DMRs/$DMRtype/$DMRfile.bed -b TSS_gencode_v19.bed -d -t first > ./DMRs/$DMRtype/$outfile\" | qsub -l h_vmem=1g -v PATH -cwd -N $DMRfile -o ./DMRs/$DMRtype/$DMRfile.out -e ./DMRs/$DMRtype/$DMRfile.err");

  # launch bedtools for 100 sets of random DMRs of same length
  my $randomDMRfile = $DMRfile;
  $randomDMRfile=~ s/DMRs/randomDMRs_/;

  for (my $i=1; $i <= 100; $i++) {
    $outfile = $randomDMRfile.$i.'_distanceToTSS.bed';
    system("echo \"closestBed -a ./DMRs/$DMRtype/randomized/$randomDMRfile$i.bed -b TSS_gencode_v19.bed -d -t first > ./DMRs/$DMRtype/randomized/$outfile\" | qsub -l h_vmem=1g -v PATH -cwd -N $DMRfile -o ./DMRs/$DMRtype/$DMRfile.out -e ./DMRs/$DMRtype/$DMRfile.err");
  }

  # launch bedtools for 100 sets of random DMRs of same length and same density
  my $randomDMRfile = $DMRfile;
  $randomDMRfile=~ s/DMRs/randomDMRs_/;

  for (my $i=1; $i <= 100; $i++) {
    $outfile = $randomDMRfile.$i.'_distanceToTSS.bed';
    system("echo \"closestBed -a ./DMRs/$DMRtype/randomized_control_CpG_density/$randomDMRfile$i.bed -b TSS_gencode_v19.bed -d -t first > ./DMRs/$DMRtype/randomized_control_CpG_density/$outfile\" | qsub -l h_vmem=1g -v PATH -cwd -N $DMRfile -o ./DMRs/$DMRtype/$DMRfile.out -e ./DMRs/$DMRtype/$DMRfile.err");
  }
  # last #testing;
}


