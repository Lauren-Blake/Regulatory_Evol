#!/usr/bin/perl
use strict;
use warnings;

## Julien Roux
## June 20, 2013
## Generate a map of CpGs from chimp/rhesus to human positions
## identify if these positions are also CpG in human
## This script makes use of cytosine report files in folder ~/Methylation/CpG_maps_to_hg19/cytosine_reports

## Generate BED files of CpG positions in a chromosome

scalar @ARGV >= 2 or die<<USAGE;
Usage:
perl generate_bed_by_chr.pl <chimp/rhesus> <chr>
USAGE
my $species = $ARGV[0];
my $chr = $ARGV[1];

my $dir;
my $chr_file;
if ($species eq 'chimp'){
  $dir = './cytosine_report/chimp/C4K1_cytosine_report/';
  $chr_file = "C4K1_trimmed.fq.gz_bismark.deduplicated.CpG_report.chr$chr.txt.gz";
}
if ($species eq 'rhesus'){
  $dir = './cytosine_report/rhesus/R1K1_cytosine_report/';
  $chr_file = "R1K1_trimmed.fq.gz_bismark.deduplicated.CpG_report.chr$chr.txt.gz";
}

print "Generating BED file of CpG positons for chromosome $chr...\n";
my $index = 1;
## read the cytosine file (0-based) and filter to keep only CpGs
## Append filtered lines to a compressed file with index as added column
open(IN, "zcat $dir$chr_file |") or die "Can't read input file for chromosome $chr\n";
my $outfile = "by_chr_qsub/$species\_all_CpGs_$chr.bed.gz";
open(OUT, "| gzip -9 - > $outfile") or die "Can't open CpG output file";
while (defined (my $line = <IN>)) {
  chomp $line;
  my @tmp = split("\t", $line);
  next if ($tmp[5] ne 'CG'); ## only keeping CpGs
  print OUT $tmp[0], "\t", $tmp[1], "\t", $tmp[1]+1, "\t", $tmp[2], "\t", $tmp[6], "\t", $index, "\n";
  $index++;
}
close IN;
close OUT;

exit;


