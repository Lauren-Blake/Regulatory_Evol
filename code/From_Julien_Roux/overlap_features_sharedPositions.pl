#!/usr/bin/perl
use strict;
use warnings;

## Julien Roux
## Sep 12, 2014
## Read bed of annotated features of human genomes, and intersect with file of interesting positions

## Notes:
## - First we use test set of ~5M shared CpG between species: this can be changed and a different set of positions loaded!
## - We don't care about strand because our data are projected on + strand and methylation is usually symetrical


# record all CpG sites
print "Loading list of positions to overlap...\n";
my %CpGs;

## open(BED, "zcat ../bsseq/combine_data_perl/evidence.gz | cut -f1,2 | " ) or die ("Cannot open bedGraph file\n");
# Directly create a file of studied positions:
# zcat ../bsseq/combine_data_perl/evidence.gz | awk '{ print $1 "\t" $2 "\t" $2+1}' > sharedPositions.bed
# remove header
# gzip -9 sharedPositions.bed.gz

open(BED, "zcat sharedPositions.bed.gz | " ) or die ("Cannot open bedGraph file\n"); ## same thing
# my $header = <BED>;  ## We missed one position because of this line!
while (defined (my $line = <BED>)) {
  chomp $line;
  my @tmp = split("\t", $line);
  # chr -> position
  $CpGs{$tmp[0]}->{$tmp[1]} = ();
}
close BED;

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
# zcat gencode.v19.annotation.gtf.gz | grep "exon_number 1" | grep "gene_type \"protein_coding\"" | intersectBed -c -a sharedPositions.bed.gz -b stdin | less

# only get the matching positions
# zcat gencode.v19.annotation.gtf.gz | grep "exon_number 1" | grep "gene_type \"protein_coding\"" | intersectBed -u -a sharedPositions.bed.gz -b stdin | wc

print "Querying GENCODE GTF file for exons positions...\n";
open (GTF, "zcat gencode.v19.annotation.gtf.gz | grep \"exon_id\" | grep \"gene_type \\\"protein_coding\\\"\" | intersectBed -u -a sharedPositions.bed.gz -b stdin |");
while (defined (my $line = <GTF>)) {
  chomp $line;
  my @tmp = split("\t", $line);

  # tag each position
  $CpGs{$tmp[0]}->{$tmp[1]}->{'exon'} = 1;
}
close GTF;

print "Querying GENCODE GTF file for FIRST exons positions...\n";
open (GTF, "zcat gencode.v19.annotation.gtf.gz | grep \"exon_number 1\" | grep \"gene_type \\\"protein_coding\\\"\" | intersectBed -u -a sharedPositions.bed.gz -b stdin |");
while (defined (my $line = <GTF>)) {
  chomp $line;
  my @tmp = split("\t", $line);

  # tag each position
  $CpGs{$tmp[0]}->{$tmp[1]}->{'first_exon'} = 1;
}
close GTF;
# TO DO? last exon?

print "Querying GENCODE GTF file for CDS positions...\n";
open (GTF, "zcat gencode.v19.annotation.gtf.gz | grep \"\tCDS\t\" | grep \"gene_type \\\"protein_coding\\\"\" | intersectBed -u -a sharedPositions.bed.gz -b stdin |");
while (defined (my $line = <GTF>)) {
  chomp $line;
  my @tmp = split("\t", $line);

  # tag each position
  $CpGs{$tmp[0]}->{$tmp[1]}->{'CDS'} = 1;
}
close GTF;

print "Querying GENCODE GTF file for UTR positions...\n";
# open (GTF, "zcat gencode.v19.annotation.gtf.gz | grep \"\tUTR\t\" | grep \"gene_type \\\"protein_coding\\\"\" | intersectBed -u -a sharedPositions.bed.gz -b stdin |");
# while (defined (my $line = <GTF>)) {
#   chomp $line;
#   my @tmp = split("\t", $line);

#   # tag each position
#   $CpGs{$tmp[0]}->{$tmp[1]}->{'UTR'} = 1;
# }
# close GTF;

# Considering 3' UTR and 5' UTR separately. See http://davetang.org/muse/2013/01/18/defining-genomic-regions/
# perl UTR.pl gencode.v19.annotation.gtf.gz > UTRs.bed
open (GTF, "grep \"3_UTR\" UTRs.bed | intersectBed -u -a sharedPositions.bed.gz -b stdin |");
while (defined (my $line = <GTF>)) {
  chomp $line;
  my @tmp = split("\t", $line);

  # tag each position
  $CpGs{$tmp[0]}->{$tmp[1]}->{'3_UTR'} = 1;
}
close GTF;

open (GTF, "grep \"5_UTR\" UTRs.bed | intersectBed -u -a sharedPositions.bed.gz -b stdin |");
while (defined (my $line = <GTF>)) {
  chomp $line;
  my @tmp = split("\t", $line);

  # tag each position
  $CpGs{$tmp[0]}->{$tmp[1]}->{'5_UTR'} = 1;
}
close GTF;


print "Querying BED for intronic positions...\n";
## Not so trivial:
## See http://davetang.org/muse/2013/01/18/defining-genomic-regions/
# zcat gencode.v19.annotation.gtf.gz | grep "gene_type \"protein_coding\"" | awk 'BEGIN{OFS="\t";} $3=="exon" {print $1,$4-1,$5}' | sortBed | mergeBed -i - | gzip -9 > gencode_v19_exon_merged.bed.gz
# zcat gencode.v19.annotation.gtf.gz | grep "gene_type \"protein_coding\"" | awk 'BEGIN{OFS="\t";} $3=="gene" {print $1,$4-1,$5}' | sortBed | subtractBed -a stdin -b gencode_v19_exon_merged.bed.gz | gzip -9 > gencode_v19_intron.bed.gz
# Intersecting the two files shouldn't produce any output: 
# intersectBed -a gencode_v19_exon_merged.bed.gz -b gencode_v19_intron.bed.gz

open (BED, "intersectBed -u -a sharedPositions.bed.gz -b gencode_v19_intron.bed.gz |");
while (defined (my $line = <BED>)) {
  chomp $line;
  my @tmp = split("\t", $line);

  # tag each position
  $CpGs{$tmp[0]}->{$tmp[1]}->{'intron'} = 1;
}
close BED;
# TO DO? first intron? last intron?


print "Querying BED file for intergenic positions...\n";
## Agin from: http://davetang.org/muse/2013/01/18/defining-genomic-regions/
# mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
#         "select chrom, size from hg19.chromInfo"  > hg19.genome
# Here we want all genes (not only protein-coding genes)
# zcat gencode.v19.annotation.gtf.gz | awk 'BEGIN{OFS="\t";} $3=="gene" {print $1,$4-1,$5}' | sortBed | complementBed -i stdin -g hg19.genome | gzip -9 > gencode_v19_intergenic.bed.gz

open (BED, "intersectBed -u -a sharedPositions.bed.gz -b gencode_v19_intergenic.bed.gz |");
while (defined (my $line = <BED>)) {
  chomp $line;
  my @tmp = split("\t", $line);

  # tag each position
  $CpGs{$tmp[0]}->{$tmp[1]}->{'intergenic'} = 1;
}
close BED;


## conserved_elements.bed.gz (UCSC)
#     Table phastConsElements46wayPlacenta
#       Log-odds score equal to its log probability under the conserved model minus its log probability under the non-conserved model.
#       The "score" field associated with this track contains transformed log-odds scores, taking values between 0 and 1000. (The scores are transformed using a monotonic function of the form a * log(x) + b.)

print "Querying BED file for conserved positions...\n";
open (BED, "intersectBed -u -a sharedPositions.bed.gz -b conserved_elements.bed.gz |");
while (defined (my $line = <BED>)) {
  chomp $line;
  my @tmp = split("\t", $line);

  # tag each position
  $CpGs{$tmp[0]}->{$tmp[1]}->{'conserved'} = 1;
}
close BED;

##   CpG_islands_UCSC_hg19_full.txt / CpG_islands_UCSC_hg19.bed (UCSC)
print "Querying BED file for human CpG islands...\n";
open (BED, "intersectBed -u -a sharedPositions.bed.gz -b CpG_islands_UCSC_hg19.bed |");
while (defined (my $line = <BED>)) {
  chomp $line;
  my @tmp = split("\t", $line);

  # tag each position
  $CpGs{$tmp[0]}->{$tmp[1]}->{'CpG_island'} = 1;
}
close BED;

## CGI shores
print "Querying BED file for CpG islands shores...\n";
open (BED, "bedtools flank -i CpG_islands_UCSC_hg19.bed -g hg19.genome -b 2000 | intersectBed -u -a sharedPositions.bed.gz -b stdin |");
while (defined (my $line = <BED>)) {
  chomp $line;
  my @tmp = split("\t", $line);

  # tag each position
  $CpGs{$tmp[0]}->{$tmp[1]}->{'CpG_island_shore'} = 1;
}
close BED;

## CGI shelves
print "Querying BED file for CpG islands shelves...\n";
## we first merge the CGI and the CGI shores files, then we take 2kb of flanking regions
open (BED, "bedtools flank -i CpG_islands_UCSC_hg19.bed -g hg19.genome -b 2000 | cat - CpG_islands_UCSC_hg19.bed | sortBed | mergeBed -i - | bedtools flank -i stdin -g hg19.genome -b 2000 | intersectBed -u -a sharedPositions.bed.gz -b stdin |");
while (defined (my $line = <BED>)) {
  chomp $line;
  my @tmp = split("\t", $line);

  # tag each position
  $CpGs{$tmp[0]}->{$tmp[1]}->{'CpG_island_shelf'} = 1;
}
close BED;

## TO DO: CGI / CGI shores / CGI shelves in chimp / rhesus
## We need to use the CGI data in human coordinates
## See orthologous_CGI/ folder

# ## CpG islands chimp
# ## CpG_islands_UCSC_panTro3.bed (UCSC)
# print "Querying BED file for chimp CpG islands...\n";
# open (BED, "intersectBed -u -a sharedPositions.bed.gz -b CpG_islands_UCSC_panTro3.bed |");
# while (defined (my $line = <BED>)) {
#   chomp $line;
#   my @tmp = split("\t", $line);

#   # tag each position
#   $CpGs{$tmp[0]}->{$tmp[1]}->{'CpG_island_chimp'} = 1;
# }
# close BED;

# ## CGI shores chimp
# ## We need to filter Y, random and unmapped chromosomes
# print "Querying BED file for CpG islands shores...\n";
# open (BED, "bedtools flank -i CpG_islands_UCSC_panTro3_filtered.bed -g panTro3_chrs_length.genome -b 2000 | intersectBed -u -a sharedPositions.bed.gz -b stdin |");
# while (defined (my $line = <BED>)) {
#   chomp $line;
#   my @tmp = split("\t", $line);

#   # tag each position
#   $CpGs{$tmp[0]}->{$tmp[1]}->{'CpG_island_shore_chimp'} = 1;
# }
# close BED;

# ## CGI shelves chimp
# print "Querying BED file for CpG islands shelves...\n";
# ## we first merge the CGI and the CGI shores files, then we take 2kb of flanking regions
# open (BED, "bedtools flank -i CpG_islands_UCSC_panTro3_filtered.bed -g panTro3_chrs_length.genome -b 2000 | cat - CpG_islands_UCSC_panTro3_filtered.bed | sortBed | mergeBed -i - | bedtools flank -i stdin -g panTro3_chrs_length.genome -b 2000 | intersectBed -u -a sharedPositions.bed.gz -b stdin |");
# while (defined (my $line = <BED>)) {
#   chomp $line;
#   my @tmp = split("\t", $line);

#   # tag each position
#   $CpGs{$tmp[0]}->{$tmp[1]}->{'CpG_island_shelf_chimp'} = 1;
# }
# close BED;


# ## CpG islands rhesus
# ## CpG_islands_UCSC_rheMac2.bed (UCSC)
# print "Querying BED file for rhesus CpG islands...\n";
# open (BED, "intersectBed -u -a sharedPositions.bed.gz -b CpG_islands_UCSC_rheMac2.bed |");
# while (defined (my $line = <BED>)) {
#   chomp $line;
#   my @tmp = split("\t", $line);

#   # tag each position
#   $CpGs{$tmp[0]}->{$tmp[1]}->{'CpG_island_rhesus'} = 1;
# }
# close BED;

# ## CGI shores rhesus
# ## We need to filter Y, random and unmapped chromosomes
# print "Querying BED file for CpG islands shores...\n";
# open (BED, "bedtools flank -i CpG_islands_UCSC_rheMac2_filtered.bed -g rheMac2_chrs_length.genome -b 2000 | intersectBed -u -a sharedPositions.bed.gz -b stdin |");
# while (defined (my $line = <BED>)) {
#   chomp $line;
#   my @tmp = split("\t", $line);

#   # tag each position
#   $CpGs{$tmp[0]}->{$tmp[1]}->{'CpG_island_shore_rhesus'} = 1;
# }
# close BED;

# ## CGI shelves rhesus
# print "Querying BED file for CpG islands shelves...\n";
# ## we first merge the CGI and the CGI shores files, then we take 2kb of flanking regions
# open (BED, "bedtools flank -i CpG_islands_UCSC_rheMac2_filtered.bed -g rheMac2_chrs_length.genome -b 2000 | cat - CpG_islands_UCSC_rheMac2_filtered.bed | sortBed | mergeBed -i - | bedtools flank -i stdin -g rheMac2_chrs_length.genome -b 2000 | intersectBed -u -a sharedPositions.bed.gz -b stdin |");
# while (defined (my $line = <BED>)) {
#   chomp $line;
#   my @tmp = split("\t", $line);

#   # tag each position
#   $CpGs{$tmp[0]}->{$tmp[1]}->{'CpG_island_shelf_rhesus'} = 1;
# }
# close BED;


## repeats.bed.gz
## UCSC RepeatMasker track
print "Querying BED file for repeats positions...\n";
open (BED, "intersectBed -u -a sharedPositions.bed.gz -b repeats.bed.gz |");
while (defined (my $line = <BED>)) {
  chomp $line;
  my @tmp = split("\t", $line);

  # tag each position
  $CpGs{$tmp[0]}->{$tmp[1]}->{'repeat'} = 1;
}
close BED;
# TO DO? This is encompassing all types of repeats & TEs. Should we filter this?

## Promoters (2kb upstream)
## Or do this directly from GTF file: see http://davetang.org/muse/2013/01/18/defining-genomic-regions/
# perl promoter.pl gencode.v19.annotation.gtf.gz 2000 > promoters_2kb_gencode_v19.bed
print "Querying BED file for promoter positions...\n";
open (BED, "intersectBed -u -a sharedPositions.bed.gz -b promoters_2kb_gencode_v19.bed |");
while (defined (my $line = <BED>)) {
  chomp $line;
  my @tmp = split("\t", $line);

  # tag each position
  $CpGs{$tmp[0]}->{$tmp[1]}->{'promoter'} = 1;
}
close BED;

## considering only promoters of protein_coding genes
# perl promoter_coding.pl gencode.v19.annotation.gtf.gz 2000 > promoters_coding_2kb_gencode_v19.bed
print "Querying BED file for coding promoter positions...\n";
open (BED, "intersectBed -u -a sharedPositions.bed.gz -b promoters_coding_2kb_gencode_v19.bed |");
while (defined (my $line = <BED>)) {
  chomp $line;
  my @tmp = split("\t", $line);

  # tag each position
  $CpGs{$tmp[0]}->{$tmp[1]}->{'promoter_coding'} = 1;
}
close BED;

# Proximal promoters
# perl promoter.pl gencode.v19.annotation.gtf.gz 250 > promoters_250bp_gencode_v19.bed
print "Querying BED file for proximal promoter positions...\n";
open (BED, "intersectBed -u -a sharedPositions.bed.gz -b promoters_250bp_gencode_v19.bed |");
while (defined (my $line = <BED>)) {
  chomp $line;
  my @tmp = split("\t", $line);

  # tag each position
  $CpGs{$tmp[0]}->{$tmp[1]}->{'promoter_proximal'} = 1;
}
close BED;


# export matrix of features for all positions ######################################################
print "  Exporting matrix of features...\n";
open(OUT, "| gzip - > 15_features_sharedPositions.gz" ) or die ("Cannot open OUT file\n");
## headers
my @features = ('exon', 'first_exon', 'CDS', '3_UTR', '5_UTR', 'intron', 'intergenic', 'conserved', 'CpG_island', 'CpG_island_shore', 'CpG_island_shelf', 'repeat', 'promoter', 'promoter_coding', 'promoter_proximal');
print OUT "chromosome\tposition\t", join("\t", @features), "\n";

foreach my $chromosome (sort keys %CpGs) {
  foreach my $position (sort {$a <=> $b} keys %{$CpGs{$chromosome}}) {
    print OUT $chromosome, "\t", $position;
    foreach my $feature (@features){
      if (exists $CpGs{$chromosome}->{$position}->{$feature}){
        print OUT "\t1";
      }
      else {
        print OUT "\t0";
      }
    }
    print OUT "\n";
  }
}
close OUT;


## not used (too slow) #####################################################################################
# print "Reading GENCODE GTF file...\n";
# open(GTF, "zcat gencode.v19.annotation.gtf.gz | ") or die "GTF file can't be opened\n";
# # Be careful, GTF is 1-based!

# my %codingExons;
# my %firstExons;
# while (defined (my $line = <GTF>)) {
#   next if ($line =~ m/^#/);
#   chomp $line;
#   my @tmp = split("\t", $line);

#   $tmp[8] =~ m/gene\_type\s\"(.+?)\"\;/; # we want a non-greedy regex, thus the ?
#   my $type = $1;

#   if (($tmp[2] eq "exon") and ($type eq "protein_coding")){
#     # chromosome -> start = end
#     $codingExons{$tmp[0]}->{$tmp[3]-1} = $tmp[4]-1;

#     $tmp[8] =~ m/exon\_number\s(\d+)/;
#     my $exon_number = $1;
#     if ($exon_number eq 1){
#       $firstExons{$tmp[0]}->{$tmp[3]-1} = $tmp[4]-1;
#     }
#   }
# }
# close GTF;

# # intersect
# foreach my $chr (sort keys %CpGs){
#  POS:
#   foreach my $position (sort {$a <=> $b} keys %{$CpGs{$chr}}){
#     my $flag = 0;
#     foreach my $firstExonStart (keys %{$firstExons{$chr}}){
#       if (($firstExonStart <= $position) and ($firstExons{$chr}->{$firstExonStart} >= $position)){
#         $flag = 1;
#       }
#       next POS if ($flag eq 1);
#     }
#     print "$chr\t$position\t$flag\n";
#   }
# }
# ## same length as original bed file? 5228469
# ## problem: very slow!
exit;
