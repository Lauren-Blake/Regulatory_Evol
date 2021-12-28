#!/usr/bin/perl
use strict;
use warnings;

## Julien Roux
## Oct 14, 2015
# - Run on spudhead
# - Launches bigWigAverageOverBed to calculate mean RS score (GERP) for all DMRs, and all the random DMRs files of same length + control for CpG density

# first get all files of DMR positions
my @files = glob("../bsseq/DMRs/*/*_DMRs.txt");

foreach my $file (@files){
  print $file, "\n";
  $file =~ m/DMRs\/(.+)\/(.+\_DMRs)\.txt/;
  my $DMRtype = $1;
  my $DMRfile = $2;

  # Individual DMRs not launched (not in paper)
  next if ($DMRtype eq 'individuals');

  # Species DMRs not launched (for now)
  # next if ($DMRtype eq 'species');

  # Tissue DMRs not launched (for now)
  next if ($DMRtype eq 'tissues');

  # Some DMR files in ./bsseq should not be taken into account
  next if ($DMRfile =~ m/\_overlap\_/);

  # ## used to relaunch only some of the jobs, e.g Rhesus individual DMRs
  # $DMRfile =~ m/(.+)\_\d\_\d\_DMRs/;
  # next unless ($1 eq 'Rhesus');
  # next unless ($DMRfile =~ m/Specific/);
  # next if ($DMRfile =~ m/Human_heart_liver/); # already done before

  # in the ~/Methylation/annotation/DMRs/ folder
  my $outfile1 = $DMRfile.'_meanSeqConservation.bed';
  # in the ~/Methylation/annotation/hg19.GERP_scores folder
  my $outfile2 = $DMRfile.'_meanSeqConservation.out';

  # launch bedtools for DMR file
  system("echo \"awk \'{print \\\$0 \\\"\\t\\\" NR}\' ./DMRs/$DMRtype/$DMRfile.bed | bigWigAverageOverBed hg19.GERP_scores/All_hg19_RS.bw stdin hg19.GERP_scores/results/$outfile2 -bedOut=DMRs/$DMRtype/$outfile1\" | qsub -l h_vmem=4g -v PATH -cwd -N $DMRfile -o ./DMRs/$DMRtype/$DMRfile.out -e ./DMRs/$DMRtype/$DMRfile.err");

  # launch bedtools for 100 sets of random DMRs of same length
  my $randomDMRfile = $DMRfile;
  $randomDMRfile=~ s/DMRs/randomDMRs_/;

  for (my $i=1; $i <= 100; $i++) {
    $outfile1 = $randomDMRfile.$i.'_meanSeqConservation.bed';
    $outfile2 = $randomDMRfile.$i.'_meanSeqConservation.out';
    system("echo \"awk \'{print\ \\\$0 \\\"\\t\\\" NR}\' ./DMRs/$DMRtype/randomized/$randomDMRfile$i.bed | bigWigAverageOverBed hg19.GERP_scores/All_hg19_RS.bw stdin hg19.GERP_scores/results/$outfile2 -bedOut=DMRs/$DMRtype/randomized/$outfile1\" | qsub -l h_vmem=4g -v PATH -cwd -N $DMRfile -o ./DMRs/$DMRtype/$DMRfile.out -e ./DMRs/$DMRtype/$DMRfile.err");
  }

  # launch bedtools for 100 sets of random DMRs of same length and same density
  for (my $i=1; $i <= 100; $i++) {
    $outfile1 = $randomDMRfile.$i.'_meanSeqConservation.bed';
    $outfile2 = $randomDMRfile.$i.'_meanSeqConservation.out';
    system("echo \"awk \'{print \\\$0 \\\"\\t\\\" NR}\' ./DMRs/$DMRtype/randomized_control_CpG_density/$randomDMRfile$i.bed | bigWigAverageOverBed hg19.GERP_scores/All_hg19_RS.bw stdin hg19.GERP_scores/results/$outfile2 -bedOut=DMRs/$DMRtype/randomized_control_CpG_density/$outfile1\" | qsub -l h_vmem=4g -v PATH -cwd -N $DMRfile -o ./DMRs/$DMRtype/$DMRfile.out -e ./DMRs/$DMRtype/$DMRfile.err");
  }
}

