#!/usr/bin/perl
use strict;
use warnings;
$| = 1;

## Julien Roux
## Aug 28, 2014
## Takes a chimp of rhesus bedGraph (collapsed) file and convert it to human coordinates using map file in ../CpG_maps_to_hg19
## This script is adapted for the output of identify_strand_problems.pl

scalar @ARGV >= 3 or die<<USAGE;
Usage:
perl map_to_hg19.pl <compressed bedGraph> <path to compressed map files> <species>
USAGE
## e.g. perl map_to_hg19.pl ./Methylation_3_C4K1/C4K1.problematic_positions.bedGraph.gz ../CpG_maps_to_hg19/by_chr_qsub/ chimp

my $bedgraph = $ARGV[0];
my $path_to_maps = $ARGV[1];
my $species = $ARGV[2];

my %to_map;
my $last_chr;
open(BED, "zcat $bedgraph |") or die "Can't read bedgraph file\n";
my $outfile1 = $bedgraph;
$outfile1 =~ s/\.bedGraph\.gz/_mapped\.bedGraph/; ## this files need to be sorted so it shoubn't be compressed now
my $outfile2 = $bedgraph;
$outfile2 =~ s/\.bedGraph\.gz/_mapped_CpG\.bedGraph/; ## this files need to be sorted so it shoubn't be compressed now
my $outfile3 = $bedgraph;
$outfile3 =~ s/\.bedGraph\.gz/_unmapped.bedGraph\.gz/;
open(OUT1, '>', "$outfile1") or die "Can't open output file 1";
open(OUT2, '>', "$outfile2") or die "Can't open output file 2";
open(OUT3, "| gzip -9 - > $outfile3") or die "Can't open output file 3";
print "Reading bedGraph file...\n";
while (defined (my $line = <BED>)) {
  chomp $line;

  ## this file is sorted, so we can read it linearly and keep things in memory for each chromosome. When we reach a new chromosome, we map positions by reading the appropriate map file, and then we can empty the memory and start over with another chromosome.
  my @tmp = split("\t", $line);
  ## we don't want to treat lambda phage chromosome
  next if ($tmp[0] =~ /gi\|215104\|gb\|J02459\.1/);

  ## if the chromosome is different than on previous line
  if ((defined $last_chr) and ($last_chr ne $tmp[0])){
    map_chr($last_chr, %to_map);

    ## erase data to map
    %to_map = ();
  }

  ## store the positions to map in hash %to_map
  $to_map{$tmp[1]}++;
  $last_chr = $tmp[0];
}
close BED;

## Don't forget to map the last chromosome
map_chr($last_chr, %to_map);

sub map_chr {
  my ($chromosome, %coords_to_map) = @_;
  print "Mapping positions for chromosome $chromosome...\n";

  ## read map file
  open(IN, "zcat $path_to_maps/$species\_ToHg19_$chromosome.map.gz |") or die "Can't read map file for chromosome $chromosome\n";
  while (defined (my $line = <IN>)) {
    chomp $line;
    my @tmp = split("\t", $line);
    if (exists $coords_to_map{$tmp[2]}) {
      ## print new bedGraph entry
      ## strand info is taken from the map, but the original strand info is also printed (4th column)
      print OUT1 $tmp[6], "\t", $tmp[7], "\t", $tmp[9], "\n";

      ## print another file with only CpG positions
      if ($tmp[10] =~ m/^CG/) {
        print OUT2 $tmp[6], "\t", $tmp[7], "\t", $tmp[9], "\n";
      }
      delete $coords_to_map{$tmp[2]};
    }
  }
  close IN;

  ## print unmapped data to another file
  print "Printing unmapped positions for chromosome $chromosome...\n";
  foreach my $start (sort keys %coords_to_map) {
    print OUT3 $chromosome, "\t", $start, "\t-\n";
  }
}

close OUT1;
close OUT2;
close OUT3;

exit;


