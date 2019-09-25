#!/usr/bin/perl
use strict;
use warnings;

## Julien Roux
## Sep 19, 2014
## Read bed of annotated features of human genomes, and intersect with file of DMRs positions

## Notes:
## - Tested on a set of Human-chimp DMRs in liver
## - We don't care about strand because our data are projected on + strand and methylation is usually symetrical

scalar @ARGV >= 1 or die<<USAGE;
Usage:
perl overlap_features_DMRs.pl <DMR BED file with full path>
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

## Using Bedtools intersect + filtering GTF with grep for each feature

## GTF format definition
# 1-based format (BED is 0-based for start position, 1-based for end-position)
# http://www.gencodegenes.org/gencodeformat.html
# E.g.: features delineated in GTF file:
# CDS
# exon
# gene
# Selenocysteine
# start_codon
# stop_codon
# transcript
# UTR

## We're tagging
# From protein-coding genes
#   exonic positions
#   exonic positions from first exon
#   CDS positions
#   UTR positions
#   intronic positions
#   promoter

## TO DO? Tag positions from pseudogenes? non-coding RNA?
# This seems a bit of a mess:
# see http://www.gencodegenes.org/gencode_biotypes.html


## Bedtools v2.17.0
# query giving the number of overlaping features for each position: 
# zcat gencode.v19.annotation.gtf.gz | grep "exon_number 1" | grep "gene_type \"protein_coding\"" | intersectBed -c -a $DMRfile -b stdin | less

# only return the original BED lines with a match on 2nd file
# zcat gencode.v19.annotation.gtf.gz | grep "exon_number 1" | grep "gene_type \"protein_coding\"" | intersectBed -u -a $DMRfile -b stdin | wc

print "Querying GENCODE GTF file for exons positions...\n";
open (GTF, "zcat gencode.v19.annotation.gtf.gz | grep \"exon_id\" | grep \"gene_type \\\"protein_coding\\\"\" | intersectBed -u -a $DMRfile -b stdin |") or die ("Bedtools query failed\n");
while (defined (my $line = <GTF>)) {
  chomp $line;
  my @tmp = split("\t", $line);

  # tag each position
  $DMRs{$tmp[0]}->{$tmp[1]}->{$tmp[2]}->{'exon'} = 1;
}
close GTF;

print "Querying GENCODE GTF file for FIRST exons positions...\n";
open (GTF, "zcat gencode.v19.annotation.gtf.gz | grep \"exon_number 1\" | grep \"gene_type \\\"protein_coding\\\"\" | intersectBed -u -a $DMRfile -b stdin |") or die ("Bedtools query failed\n");
while (defined (my $line = <GTF>)) {
  chomp $line;
  my @tmp = split("\t", $line);

  # tag each position
  $DMRs{$tmp[0]}->{$tmp[1]}->{$tmp[2]}->{'first_exon'} = 1;
}
close GTF;
# TO DO? last exon?

print "Querying GENCODE GTF file for CDS positions...\n";
open (GTF, "zcat gencode.v19.annotation.gtf.gz | grep \"\tCDS\t\" | grep \"gene_type \\\"protein_coding\\\"\" | intersectBed -u -a $DMRfile -b stdin |") or die ("Bedtools query failed\n");
while (defined (my $line = <GTF>)) {
  chomp $line;
  my @tmp = split("\t", $line);

  # tag each position
  $DMRs{$tmp[0]}->{$tmp[1]}->{$tmp[2]}->{'CDS'} = 1;
}
close GTF;

print "Querying GENCODE GTF file for UTR positions...\n";
# open (GTF, "zcat gencode.v19.annotation.gtf.gz | grep \"\tUTR\t\" | grep \"gene_type \\\"protein_coding\\\"\" | intersectBed -u -a $DMRfile -b stdin |") or die ("Bedtools query failed\n");
# while (defined (my $line = <GTF>)) {
#   chomp $line;
#   my @tmp = split("\t", $line);

#   # tag each position
#   $DMRs{$tmp[0]}->{$tmp[1]}->{$tmp[2]}->{'UTR'} = 1;
# }
# close GTF;

# Considering 3' UTR and 5' UTR separately. See http://davetang.org/muse/2013/01/18/defining-genomic-regions/
# perl UTR.pl gencode.v19.annotation.gtf.gz > UTRs.bed
open (GTF, "grep \"3_UTR\" UTRs.bed | intersectBed -u -a $DMRfile -b stdin |") or die ("Bedtools query failed\n");
while (defined (my $line = <GTF>)) {
  chomp $line;
  my @tmp = split("\t", $line);

  # tag each position
  $DMRs{$tmp[0]}->{$tmp[1]}->{$tmp[2]}->{'3_UTR'} = 1;
}
close GTF;

open (GTF, "grep \"5_UTR\" UTRs.bed | intersectBed -u -a $DMRfile -b stdin |") or die ("Bedtools query failed\n");
while (defined (my $line = <GTF>)) {
  chomp $line;
  my @tmp = split("\t", $line);

  # tag each position
  $DMRs{$tmp[0]}->{$tmp[1]}->{$tmp[2]}->{'5_UTR'} = 1;
}
close GTF;


print "Querying BED for intronic positions...\n";
## Not so trivial:
## See http://davetang.org/muse/2013/01/18/defining-genomic-regions/
# zcat gencode.v19.annotation.gtf.gz | grep "gene_type \"protein_coding\"" | awk 'BEGIN{OFS="\t";} $3=="exon" {print $1,$4-1,$5}' | sortBed | mergeBed -i - | gzip -9 > gencode_v19_exon_merged.bed.gz
# zcat gencode.v19.annotation.gtf.gz | grep "gene_type \"protein_coding\"" | awk 'BEGIN{OFS="\t";} $3=="gene" {print $1,$4-1,$5}' | sortBed | subtractBed -a stdin -b gencode_v19_exon_merged.bed.gz | gzip -9 > gencode_v19_intron.bed.gz
# Intersecting the two files shouldn't produce any output: 
# intersectBed -a gencode_v19_exon_merged.bed.gz -b gencode_v19_intron.bed.gz

open (BED, "intersectBed -u -a $DMRfile -b gencode_v19_intron.bed.gz |") or die ("Bedtools query failed\n");
while (defined (my $line = <BED>)) {
  chomp $line;
  my @tmp = split("\t", $line);

  # tag each position
  $DMRs{$tmp[0]}->{$tmp[1]}->{$tmp[2]}->{'intron'} = 1;
}
close BED;
# TO DO? first intron? last intron?


print "Querying BED file for intergenic positions...\n";
## Agin from: http://davetang.org/muse/2013/01/18/defining-genomic-regions/
# mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
#         "select chrom, size from hg19.chromInfo"  > hg19.genome
# Here we want all genes (not only protein-coding genes)
# zcat gencode.v19.annotation.gtf.gz | awk 'BEGIN{OFS="\t";} $3=="gene" {print $1,$4-1,$5}' | sortBed | complementBed -i stdin -g hg19.genome | gzip -9 > gencode_v19_intergenic.bed.gz

open (BED, "intersectBed -u -a $DMRfile -b gencode_v19_intergenic.bed.gz |") or die ("Bedtools query failed\n");
while (defined (my $line = <BED>)) {
  chomp $line;
  my @tmp = split("\t", $line);

  # tag each position
  $DMRs{$tmp[0]}->{$tmp[1]}->{$tmp[2]}->{'intergenic'} = 1;
}
close BED;


## conserved_elements.bed.gz (UCSC)
#     Table phastConsElements46wayPlacenta
#       Log-odds score equal to its log probability under the conserved model minus its log probability under the non-conserved model.
#       The "score" field associated with this track contains transformed log-odds scores, taking values between 0 and 1000. (The scores are transformed using a monotonic function of the form a * log(x) + b.)

print "Querying BED file for conserved positions...\n";
open (BED, "intersectBed -u -a $DMRfile -b conserved_elements.bed.gz |") or die ("Bedtools query failed\n");
while (defined (my $line = <BED>)) {
  chomp $line;
  my @tmp = split("\t", $line);

  # tag each position
  $DMRs{$tmp[0]}->{$tmp[1]}->{$tmp[2]}->{'conserved'} = 1;
}
close BED;

##   CpG_islands_UCSC_hg19_full.txt / CpG_islands_UCSC_hg19.bed (UCSC)
print "Querying BED file for CpG islands...\n";
open (BED, "intersectBed -u -a $DMRfile -b CpG_islands_UCSC_hg19.bed |") or die ("Bedtools query failed\n");
while (defined (my $line = <BED>)) {
  chomp $line;
  my @tmp = split("\t", $line);

  # tag each position
  $DMRs{$tmp[0]}->{$tmp[1]}->{$tmp[2]}->{'CpG_island'} = 1;
}
close BED;

## CGI shores
print "Querying BED file for CpG islands shores...\n";
open (BED, "bedtools flank -i CpG_islands_UCSC_hg19.bed -g hg19.genome -b 2000 | intersectBed -u -a $DMRfile -b stdin |") or die ("Bedtools query failed\n");
while (defined (my $line = <BED>)) {
  chomp $line;
  my @tmp = split("\t", $line);

  # tag each position
  $DMRs{$tmp[0]}->{$tmp[1]}->{$tmp[2]}->{'CpG_island_shore'} = 1;
}
close BED;

## CGI shelves
print "Querying BED file for CpG islands shelves...\n";
## we first merge the CGI and the CGI shores files, then we take 2kb of flanking regions
open (BED, "bedtools flank -i CpG_islands_UCSC_hg19.bed -g hg19.genome -b 2000 | cat - CpG_islands_UCSC_hg19.bed | sortBed | mergeBed -i - | bedtools flank -i stdin -g hg19.genome -b 2000 | intersectBed -u -a $DMRfile -b stdin |") or die ("Bedtools query failed\n");
while (defined (my $line = <BED>)) {
  chomp $line;
  my @tmp = split("\t", $line);

  # tag each position
  $DMRs{$tmp[0]}->{$tmp[1]}->{$tmp[2]}->{'CpG_island_shelf'} = 1;
}
close BED;



## repeats.bed.gz
## UCSC RepeatMasker track
print "Querying BED file for repeats positions...\n";
open (BED, "intersectBed -u -a $DMRfile -b repeats.bed.gz |") or die ("Bedtools query failed\n");
while (defined (my $line = <BED>)) {
  chomp $line;
  my @tmp = split("\t", $line);

  # tag each position
  $DMRs{$tmp[0]}->{$tmp[1]}->{$tmp[2]}->{'repeat'} = 1;
}
close BED;
# TO DO? This is encompassing all types of repeats & TEs. Should we filter this?


## promoters_2kb_gencode_v19_UCSC.bed
#     -> promoter: 2kb from TSS 
#     http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=388678057_1GDXJVaocliVHtxeLYPjxwHQYFCl&boolshad.hgta_printCustomTrackHeaders=0&hgta_ctName=tb_wgEncodeGencodeCompV19&hgta_ctDesc=table+browser+query+on+wgEncodeGencodeCompV19&hgta_ctVis=pack&hgta_ctUrl=&fbQual=upstreamAll&fbUpBases=2000&fbExonBases=0&fbIntronBases=0&fbDownBases=200&hgta_doGetBed=get+BED
# TO DO: should we filter for promoters of protein coding genes?

## Or do this directly from GTF file: see http://davetang.org/muse/2013/01/18/defining-genomic-regions/
# perl promoter.pl gencode.v19.annotation.gtf.gz 2000 > promoters_2kb_gencode_v19.bed
print "Querying BED file for promoter positions...\n";
open (BED, "intersectBed -u -a $DMRfile -b promoters_2kb_gencode_v19.bed |") or die ("Bedtools query failed\n");
while (defined (my $line = <BED>)) {
  chomp $line;
  my @tmp = split("\t", $line);

  # tag each position
  $DMRs{$tmp[0]}->{$tmp[1]}->{$tmp[2]}->{'promoter'} = 1;
}
close BED;

## considering only promoters of protein_coding genes
# perl promoter_coding.pl gencode.v19.annotation.gtf.gz 2000 > promoters_coding_2kb_gencode_v19.bed
print "Querying BED file for coding promoter positions...\n";
open (BED, "intersectBed -u -a $DMRfile -b promoters_coding_2kb_gencode_v19.bed |") or die ("Bedtools query failed\n");
while (defined (my $line = <BED>)) {
  chomp $line;
  my @tmp = split("\t", $line);

  # tag each position
  $DMRs{$tmp[0]}->{$tmp[1]}->{$tmp[2]}->{'promoter_coding'} = 1;
}
close BED;

# Proximal promoters
# perl promoter.pl gencode.v19.annotation.gtf.gz 250 > promoters_250bp_gencode_v19.bed
print "Querying BED file for proximal promoter positions...\n";
open (BED, "intersectBed -u -a $DMRfile -b promoters_250bp_gencode_v19.bed |") or die ("Bedtools query failed\n");
while (defined (my $line = <BED>)) {
  chomp $line;
  my @tmp = split("\t", $line);

  # tag each position
  $DMRs{$tmp[0]}->{$tmp[1]}->{$tmp[2]}->{'promoter_proximal'} = 1;
}
close BED;


# export matrix of features for all positions ######################################################
print "  Exporting matrix of features...\n";

my $outfile = $DMRfile;
$outfile =~ s/\.bed//;
open(OUT, "| gzip - > ".$outfile."_15_features.gz" ) or die ("Cannot open OUT file\n");
## headers
my @features = ('exon', 'first_exon', 'CDS', '3_UTR', '5_UTR', 'intron', 'intergenic', 'conserved', 'CpG_island', 'CpG_island_shore', 'CpG_island_shelf', 'repeat', 'promoter', 'promoter_coding', 'promoter_proximal');
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
