#!/usr/bin/perl
use strict;
use warnings;

## Julien Roux
## Sep 22, 2014
# - Run on spudhead
# - Collect all DMRs files, export them in BED format and launches the overlap to genomic features (script overlap_features_DMRs.pl)
# - Also generates 100 files of random DMRs (random segments of same length on same chromosome)

# first create BED files of DMR positions
my @files = glob("../bsseq/DMRs/*/*_DMRs.txt");

foreach my $file (@files){

  print $file, "\n";
  $file =~ m/\.\.\/bsseq\/DMRs\/(.+)\/(.+\_DMRs)\.txt/;
  # DMR type -> file
  my $DMRtype = $1;
  my $DMRfile = $2;


#  ## only relaunch some of the jobs
#  next unless ($DMRtype eq 'individuals');
#  $DMRfile =~ m/(.+)\_\d\_\d\_DMRs/;
#  next unless ($1 eq 'Rhesus');


  # Export BED file
  # BED end coordinate are 0-based exclusive (== 1-based)
  system("sed 1d $file | awk 'BEGIN{OFS=\"\t\";}{ print \$1,\$2,\$3+1 }' > ./DMRs/$DMRtype/$DMRfile.bed");

  # launch script overlap_features_DMRs.pl
  system("echo \"perl overlap_features_DMRs.pl ./DMRs/$DMRtype/$DMRfile.bed\" | qsub -l h_vmem=4g -v PATH -cwd -N $DMRfile -o ./DMRs/$DMRtype/$DMRfile.out -e ./DMRs/$DMRtype/$DMRfile.err");

  # Create 100 BED files of random DMRs of same length
  # First, read hg19.genome to get the length of each chromosome
  my %chromosomes;
  open(GENOME, "hg19.genome" ) or die ("Cannot open chromosome size file\n");
  my $header = <GENOME>;
  while (defined (my $line = <GENOME>)) {
    chomp $line;
    my @tmp = split("\t", $line);
    # chr -> position
    $chromosomes{$tmp[0]} = $tmp[1];
  }
  close GENOME;

  my @DMRs;
  open(IN, $file) or die ("Cannot open DMRs file\n");
  $header = <IN>;
  while (defined (my $line = <IN>)) {
    chomp $line;
    my @tmp = split("\t", $line);
    my $chr = $tmp[0];
    my $length = $tmp[2] - $tmp[1] + 1;
    push(@DMRs, [$chr, $length]); ## build an array of array
  }
  close IN;

  # base name for randomDMRs files
  my $randomDMRfile = $DMRfile;
  $randomDMRfile=~ s/DMRs/randomDMRs_/;

  # 100 randomizations
  for (my $i=1; $i <= 100; $i++) {
    open(OUT, '>', "./DMRs/$DMRtype/randomized/$randomDMRfile".$i.".bed") or die ("Cannot open OUT random DMRs file\n");

    # For each DMR, draw random segments of same length on same chromosome
    for (my $j=0; $j < scalar(@DMRs); $j++) {
      my $chr = $DMRs[$j][0];
      my $DMRlength = $DMRs[$j][1];
      my $range = $chromosomes{$chr} - $DMRlength;
      my $randomStart = int(rand($range)); # draw an integer between 0 and size of chromosome - length of DMR
      my $randomEnd = $randomStart + $DMRlength;
      print OUT "$chr\t$randomStart\t$randomEnd\n";
    }
    close OUT;

    # launch script overlap_features_DMRs.pl on each of these random DMRs files
    system("echo \"perl overlap_features_DMRs.pl ./DMRs/$DMRtype/randomized/$randomDMRfile$i.bed\" | qsub -l h_vmem=4g -v PATH -cwd -N $randomDMRfile$i -o ./DMRs/$DMRtype/randomized/$randomDMRfile$i.out -e ./DMRs/$DMRtype/randomized/$randomDMRfile$i.err");
  }

  #last; # testing
}

# TO DO? Better random set: list of random DMRs of same length and same number of CpG (same density)!

## TO DO: load in R
## Boxplot of proportion of random DMRs for each feature. Add point for real DMRs. This gives an idea of significance.


