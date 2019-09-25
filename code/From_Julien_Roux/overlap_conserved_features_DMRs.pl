#!/usr/bin/perl
use strict;
use warnings;

## Julien Roux
## Mar 11, 2016
## Read bed of conserved features between human/chimp/macaque genomes, and intersect with file of DMRs positions (or file of random DMRs positions)
## Conserved features are the following
## - promoters
## - CpG islands
## - Repeats

# This script differs from overlap_features_DMRs.pl on these points:
# - Instead of outputing 0 or 1, we output the name or coordinate of the orthologous gene / CGI repeat

scalar @ARGV >= 1 or die<<USAGE;
Usage:
perl overlap_conserved_features_DMRs.pl <DMR BED file with full path>
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

## Using Bedtools intersect with conserved promoters
## - The overlap is done only with the conserved portion of promoter, not the full promoter
## - 1b minimum overlap (20% overlap would be too stringent since we overlap with conserved portion of promoter)
## - DMRs mapping the same promoter are identified. Instead of 0/1, we output NA / Ensembl Gene ID associated with promoter

## Human promoters (2kb upstream)
print "Querying human promoter positions...\n";
open (BED, "intersectBed -wa -wb -a $DMRfile -b orthologous_promoters/promoters_2kb_gencode_v19.bed |");
while (defined (my $line = <BED>)) {
  chomp $line;
  my @tmp = split("\t", $line);

  # tag each position. The gene corresponding to the overlapping promoter is recorded
  $DMRs{$tmp[0]}->{$tmp[1]}->{$tmp[2]}->{'promoter'}->{$tmp[6]} = ();
}
close BED;

## Human promoters conserved in chimp
## Only use the intersection to do the overlap
print "Querying human promoter positions conserved in chimp...\n";
open (BED, 'awk \'{ print $1 "\t" $6 "\t" $7 "\t" $4 "\t\t" $5 "\t" $2 "\t" $3 "\t" $8}\' orthologous_promoters/promoters_2kb_gencode_v19_conserved_panTro3.bed | intersectBed -wa -wb -a '.$DMRfile.' -b - |');
while (defined (my $line = <BED>)) {
  chomp $line;
  my @tmp = split("\t", $line);

  # tag each position. The overlapping promoter coordinates is recorded
  $DMRs{$tmp[0]}->{$tmp[1]}->{$tmp[2]}->{'promoter_conserved_panTro3'}->{$tmp[6]} = ();
}
close BED;

## Human promoters conserved in rhesus
## Only use the intersection to do the overlap
print "Querying human promoter positions conserved in rhesus...\n";
open (BED, 'awk \'{ print $1 "\t" $6 "\t" $7 "\t" $4 "\t\t" $5 "\t" $2 "\t" $3 "\t" $8}\' orthologous_promoters/promoters_2kb_gencode_v19_conserved_rheMac2.bed | intersectBed -wa -wb -a '.$DMRfile.' -b - |');
while (defined (my $line = <BED>)) {
  chomp $line;
  my @tmp = split("\t", $line);

  # tag each position. The overlapping promoter coordinates is recorded
  $DMRs{$tmp[0]}->{$tmp[1]}->{$tmp[2]}->{'promoter_conserved_rheMac2'}->{$tmp[6]} = ();
}
close BED;

## Human promoters conserved in chimp+rhesus
## Only use the intersection to do the overlap
print "Querying human promoter positions conserved in chimp+rhesus...\n";
open (BED, 'awk \'{ print $1 "\t" $6 "\t" $7 "\t" $4 "\t\t" $5 "\t" $2 "\t" $3 "\t" $8}\' orthologous_promoters/promoters_2kb_gencode_v19_conserved_panTro3+rheMac2.bed | intersectBed -wa -wb -a '.$DMRfile.' -b - |');
while (defined (my $line = <BED>)) {
  chomp $line;
  my @tmp = split("\t", $line);

  # tag each position. The overlapping promoter coordinates is recorded
  $DMRs{$tmp[0]}->{$tmp[1]}->{$tmp[2]}->{'promoter_conserved_panTro3+rheMac2'}->{$tmp[6]} = ();
}
close BED;

## Human CGI 
print "Querying human CGI positions...\n";
open (BED, "intersectBed -wa -wb -a $DMRfile -b CpG_islands_UCSC_hg19.bed |");
while (defined (my $line = <BED>)) {
  chomp $line;
  my @tmp = split("\t", $line);

  # tag each position. The overlapping CGI coordinates is recorded
  $DMRs{$tmp[0]}->{$tmp[1]}->{$tmp[2]}->{'CGI'}->{"$tmp[3]:$tmp[4]-$tmp[5]"} = ();
}
close BED;

## Human CGI conserved in chimp
## Only use the intersection to do the overlap
print "Querying human CGI positions conserved in chimp...\n";
open (BED, 'awk \'{ print $1 "\t" $4 "\t" $5 "\t" $2 "\t" $3 }\' orthologous_CGI/CpG_islands_UCSC_hg19_conserved_panTro3.bed | intersectBed -wa -wb -a '.$DMRfile.' -b - |');
while (defined (my $line = <BED>)) {
  chomp $line;
  my @tmp = split("\t", $line);

  # tag each position. The overlapping CGI coordinates is recorded
  $DMRs{$tmp[0]}->{$tmp[1]}->{$tmp[2]}->{'CGI_conserved_panTro3'}->{"$tmp[3]:$tmp[6]-$tmp[7]"} = ();
}
close BED;

## Human CGI conserved in rhesus
## Only use the intersection to do the overlap
print "Querying human CGI positions conserved in rhesus...\n";
open (BED, 'awk \'{ print $1 "\t" $4 "\t" $5 "\t" $2 "\t" $3}\' orthologous_CGI/CpG_islands_UCSC_hg19_conserved_rheMac2.bed | intersectBed -wa -wb -a '.$DMRfile.' -b - |');
while (defined (my $line = <BED>)) {
  chomp $line;
  my @tmp = split("\t", $line);

  # tag each position. The overlapping CGI coordinates is recorded
  $DMRs{$tmp[0]}->{$tmp[1]}->{$tmp[2]}->{'CGI_conserved_rheMac2'}->{"$tmp[3]:$tmp[6]-$tmp[7]"} = ();
}
close BED;

## Human CGI conserved in chimp+rhesus
## Only use the intersection to do the overlap
print "Querying human CGI positions conserved in chimp+rhesus...\n";
open (BED, 'awk \'{ print $1 "\t" $4 "\t" $5 "\t" $2 "\t" $3}\' orthologous_CGI/CpG_islands_UCSC_hg19_conserved_panTro3+rheMac2.bed | intersectBed -wa -wb -a '.$DMRfile.' -b - |');
while (defined (my $line = <BED>)) {
  chomp $line;
  my @tmp = split("\t", $line);

  # tag each position. The overlapping CGI coordinates is recorded
  $DMRs{$tmp[0]}->{$tmp[1]}->{$tmp[2]}->{'CGI_conserved_panTro3+rheMac2'}->{"$tmp[3]:$tmp[6]-$tmp[7]"} = ();
}
close BED;

## Human repeats 
print "Querying human repeats positions...\n";
open (BED, "intersectBed -wa -wb -a $DMRfile -b repeats_full.txt.gz |");
while (defined (my $line = <BED>)) {
  chomp $line;
  my @tmp = split("\t", $line);

  # tag each position. The overlapping repeats coordinates is recorded
  $DMRs{$tmp[0]}->{$tmp[1]}->{$tmp[2]}->{'repeats'}->{"$tmp[3]:$tmp[4]-$tmp[5]"} = ();
}
close BED;

## Human repeats conserved in chimp
## Only use the intersection to do the overlap
print "Querying human repeats positions conserved in chimp...\n";
open (BED, 'awk \'{ print $1 "\t" $8 "\t" $9 "\t" $2 "\t" $3 }\' orthologous_repeats/repeats_hg19_conserved_panTro3.bed | intersectBed -wa -wb -a '.$DMRfile.' -b - |');
while (defined (my $line = <BED>)) {
  chomp $line;
  my @tmp = split("\t", $line);

  # tag each position. The overlapping repeats coordinates is recorded
  $DMRs{$tmp[0]}->{$tmp[1]}->{$tmp[2]}->{'repeats_conserved_panTro3'}->{"$tmp[3]:$tmp[6]-$tmp[7]"} = ();
}
close BED;

## Human repeats conserved in rhesus
## Only use the intersection to do the overlap
print "Querying human repeats positions conserved in rhesus...\n";
open (BED, 'awk \'{ print $1 "\t" $8 "\t" $9 "\t" $2 "\t" $3}\' orthologous_repeats/repeats_hg19_conserved_rheMac2.bed | intersectBed -wa -wb -a '.$DMRfile.' -b - |');
while (defined (my $line = <BED>)) {
  chomp $line;
  my @tmp = split("\t", $line);

  # tag each position. The overlapping repeats coordinates is recorded
  $DMRs{$tmp[0]}->{$tmp[1]}->{$tmp[2]}->{'repeats_conserved_rheMac2'}->{"$tmp[3]:$tmp[6]-$tmp[7]"} = ();
}
close BED;

## Human repeats conserved in chimp+rhesus
## Only use the intersection to do the overlap
print "Querying human repeats positions conserved in chimp+rhesus...\n";
open (BED, 'awk \'{ print $1 "\t" $8 "\t" $9 "\t" $2 "\t" $3}\' orthologous_repeats/repeats_hg19_conserved_panTro3+rheMac2.bed | intersectBed -wa -wb -a '.$DMRfile.' -b - |');
while (defined (my $line = <BED>)) {
  chomp $line;
  my @tmp = split("\t", $line);

  # tag each position. The overlapping repeats coordinates is recorded
  $DMRs{$tmp[0]}->{$tmp[1]}->{$tmp[2]}->{'repeats_conserved_panTro3+rheMac2'}->{"$tmp[3]:$tmp[6]-$tmp[7]"} = ();
}
close BED;

# export matrix of features for all positions ######################################################
print "Exporting matrix of features...\n";

my $outfile = $DMRfile;
$outfile =~ s/\.bed//;
open(OUT, "| gzip - > ".$outfile."_conserved_features.gz" ) or die ("Cannot open OUT file\n");
## headers
my @features = ('promoter', 'promoter_conserved_panTro3', 'promoter_conserved_rheMac2', 'promoter_conserved_panTro3+rheMac2', 'CGI', 'CGI_conserved_panTro3', 'CGI_conserved_rheMac2', 'CGI_conserved_panTro3+rheMac2', 'repeats', 'repeats_conserved_panTro3', 'repeats_conserved_rheMac2', 'repeats_conserved_panTro3+rheMac2');
print OUT "chromosome\tstart\tend\t", join("\t", @features), "\n";

foreach my $chromosome (sort keys %DMRs) {
  foreach my $start (sort {$a <=> $b} keys %{$DMRs{$chromosome}}) {
    foreach my $end (sort {$a <=> $b} keys %{$DMRs{$chromosome}->{$start}}) {
      print OUT $chromosome, "\t", $start, "\t", $end;
      foreach my $feature (@features){
        if (exists $DMRs{$chromosome}->{$start}->{$end}->{$feature}){
          print OUT "\t", join(',', keys %{$DMRs{$chromosome}->{$start}->{$end}->{$feature}});
        }
        else {
          print OUT "\tNA";
        }
      }
      print OUT "\n";
    }
  }
}
close OUT;
exit;

