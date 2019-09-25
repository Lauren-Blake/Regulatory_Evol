#!/usr/bin/perl
use strict;
use warnings;

## Julien Roux
## Aug 19 2015
## Read bed of annotated repeats in human, and intersect with file of DMRs positions (or file of random DMRs positions)

# This script differs from overlap_features_DMRs_0.2.pl on these points:
# - No need for of minimum 20% of DMR length, because we want to include as much overlaps as possible, and many repeats are short. Human-specific repeats are even not included in the reference genome often
# - different classes and families tested

scalar @ARGV >= 1 or die<<USAGE;
Usage:
perl overlap_repeats_DMRs.pl <DMR BED file with full path>
USAGE
my $DMRfile = $ARGV[0];

# record all DMRs from BED file
print "Loading list of DMR positions to overlap (file $DMRfile)...\n";
my %DMRs;
open(DMR, $DMRfile ) or die ("Cannot open DMR file\n");
# my $header = <DMR>;  ## We missed one DMR because of this line!
while (defined (my $line = <DMR>)) {
  chomp $line;
  my @tmp = split("\t", $line);
  # chr -> start -> end
  $DMRs{$tmp[0]}->{$tmp[1]}->{$tmp[2]} = ();
}
close DMR;


## repeats_full.txt.gz
## UCSC RepeatMasker
print "Querying Full repeat file...\n";
open (BED, "zcat repeats_full.txt.gz | tail -n +2 | cut -f1,2,3 | intersectBed -u -a $DMRfile -b - |") or die ("Bedtools query failed\n");
while (defined (my $line = <BED>)) {
  chomp $line;
  my @tmp = split("\t", $line);

  # tag each position
  $DMRs{$tmp[0]}->{$tmp[1]}->{$tmp[2]}->{'All_repeats'} = 1;
}
close BED;

## Grep different subtypes of repeats:

## name of repeat -> expression to grep
my %types = (
  'Low complexity' => 'Low_complexity', # different subtypes
  'Simple repeat' => 'Simple_repeat', # different subtypes
  'DNA/TcMar' => "\"DNA\tTcMar\"", # TcMar, TcMar-Mariner, TcMar-Tc2, TcMar-Tigger
  'DNA/hAT' => "\"DNA\thAT\"", #hAT, hAT-Blackjack, hAT-Charlie, hAT-Tip100
  'LINE/L1' => "\"LINE\tL1\"",
  'LINE/L2' => "\"LINE\tL2\"",
  'LTR/ERV1' => "\"LTR\tERV1\"",
  'LTR/ERVL' => "\"LTR\tERVL\"", ##ERVL and ERVL-MaLR
  'SINE/Alu' => "\"SINE\tAlu\"",
  'SINE/MIR' => "\"SINE\tMIR\"",
);
# Note some families are followed by a "?". We include them
# Note some classes are followed by a "?". We ignore them

foreach my $type (keys %types){
  print "Querying $type repeats...\n";
  open (BED, "zcat repeats_full.txt.gz | tail -n +2 | grep $types{$type} | cut -f1,2,3 | intersectBed -u -a $DMRfile -b - |") or die ("Bedtools query failed\n");
  while (defined (my $line = <BED>)) {
    chomp $line;
    my @tmp = split("\t", $line);
    
    # tag each position
    $DMRs{$tmp[0]}->{$tmp[1]}->{$tmp[2]}->{$type} = 1;
  }
  close BED;
}

print "Querying human-specific Alu repeats...\n";
open (BED, "grep 'Alu' repeats_human_specific_L1+Alu_hg19.bed | intersectBed -u -a $DMRfile -b - |") or die ("Bedtools query failed\n");
while (defined (my $line = <BED>)) {
  chomp $line;
  my @tmp = split("\t", $line);
  
  # tag each position
  $DMRs{$tmp[0]}->{$tmp[1]}->{$tmp[2]}->{'Human-specific SINE/Alu'} = 1;
}
close BED;

print "Querying human-specific L1 repeats...\n";
open (BED, "grep 'L1' repeats_human_specific_L1+Alu_hg19.bed | intersectBed -u -a $DMRfile -b - |") or die ("Bedtools query failed\n");
while (defined (my $line = <BED>)) {
  chomp $line;
  my @tmp = split("\t", $line);
  
  # tag each position
  $DMRs{$tmp[0]}->{$tmp[1]}->{$tmp[2]}->{'Human-specific SINE/L1'} = 1;
}
close BED;


# export matrix of features for all positions ######################################################
print "  Exporting matrix of features...\n";

my $outfile = $DMRfile;
$outfile =~ s/\.bed//;
open(OUT, "| gzip - > ".$outfile."_repeats.gz" ) or die ("Cannot open OUT file\n");
## headers
my @features = ('All_repeats', keys %types, 'Human-specific SINE/Alu', 'Human-specific LINE/L1');
print OUT "chromosome\tstart\tend\t", join("\t", @features), "\n";

foreach my $chromosome (sort keys %DMRs) {
  foreach my $start (sort {$a <=> $b} keys %{$DMRs{$chromosome}}) {
    foreach my $end (sort {$a <=> $b} keys %{$DMRs{$chromosome}->{$start}}) {
      print OUT $chromosome, "\t", $start, "\t", $end;
      foreach my $feature (@features){
        if (exists $DMRs{$chromosome}->{$start}->{$end}->{$feature}){
          print OUT "\t1";
        }
        else {
          print OUT "\t0";
        }
      }
      print OUT "\n";
    }
  }
}
close OUT;

exit;
