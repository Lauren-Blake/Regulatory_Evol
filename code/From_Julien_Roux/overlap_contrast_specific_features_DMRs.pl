#!/usr/bin/perl
use strict;
use warnings;

## Julien Roux
## Aug 21 2015
## intersect file of DMRs (or file of random DMRs) to features specific to the same contrast:
#           DE exons in same contrasts
#           Differentially spliced exons
#           Enhancers from FANTOM (human) in specific tissues
#           tDMRs (for sDMRs) and vice-versa
## Not launched on tissue-specific or species-specific contrasts


scalar @ARGV >= 1 or die<<USAGE;
Usage:
perl overlap_contrast_specific_features_DMRs.pl <DMR BED file with full path>
USAGE
my $DMRfile = $ARGV[0];

## TO DO: should work for randomized files
## less DMRs/species/randomized/ChimpRhesus_heart_randomDMRs_98.

$DMRfile =~ m/DMRs\/([a-z]+)\/?(.+)?\/(.+)\_(random)?DMRs(\_.+)?\.bed/;
# This regex is made to accomodate randomDMR files
my $DMRtype = $1;
my $DMRcontrast = $3;
# If random DMR file, we need to edit $DMRcontrast
$DMRcontrast =~ s/\_randomDMRs/\_DMRs/;
$DMRcontrast =~ s/DMRs\_.+/DMRs/;

my $species;
my $tissue1;
my $tissue2;
if ($DMRtype eq 'tissues'){
  $DMRcontrast =~ m/(.+)\_(.+)\_(.+)/;
  $species = $1;
  $tissue1 = $2;
  $tissue2 = $3;
  # print "$species $tissue1 $tissue2\n";
}

my $tissue;
my $species1;
my $species2;
if ($DMRtype eq 'species'){
  ## Need CamelCase matching regex
  # http://stackoverflow.com/questions/815787/what-perl-regex-can-match-camelcase-words
  $DMRcontrast =~ m/([A-Z]{1}[a-z]+)([A-Z]{1}[a-z]+)\_(.+)/;
  $species1 = $1;
  $species2 = $2;
  $tissue = $3;
  #print "$DMRcontrast $species1 $species2 $tissue\n";
}

# record all DMRs from BED file
print "Loading list of DMR positions to overlap (file $DMRfile)...\n";
my %DMRs;
open(DMR, $DMRfile ) or die ("Cannot open DMR file\n");
# my $header = <DMR>; ## We missed one DMR because of this line!
while (defined (my $line = <DMR>)) {
  chomp $line;
  my @tmp = split("\t", $line);
  # chr -> start -> end
  $DMRs{$tmp[0]}->{$tmp[1]}->{$tmp[2]} = ();
}
close DMR;

## DE exons
## We need a correspondance to the file names (in ../limma/gene_lists/DMRs_comparison/)
my %correspondance = (
  # tissues
  heart => {
    kidney => 'KidneyvsHeart',
    liver => 'HeartvsLiver',
    lung => 'HeartvsLung',
  },
  liver => {
    kidney => 'KidneyvsLiver',
    lung => 'LungvsLiver',
  },
  lung =>  {
    kidney => 'KidneyvsLung',
  },

  # species
  Human => {
    Chimp => 'ChimpvsHuman',
    Rhesus => 'HumanvsRhesus',
  },
  Chimp => {
    Rhesus => 'ChimpvsRhesus',
  },
);

## Check if BED file was already generated
my $filename = 'DE_exons/'.$DMRtype.'/'.$DMRcontrast.'.bed';
unless (-s $filename) {
  print "Generating BED files of significanly DE exons...\n";
  my $file;
  if ($DMRtype eq 'tissues') {
    $file = '../limma/gene_lists/DMRs_comparison/'.$species.$correspondance{$tissue1}->{$tissue2}.'.txt';
  }
  if ($DMRtype eq 'species') {
    $file = '../limma/gene_lists/DMRs_comparison/'.$correspondance{$species1}->{$species2}.'_'.ucfirst($tissue).'.txt';
  }
  open(DE, $file ) or die ("Cannot open DE file\n");
  # record significant DE exons
  my %DEsignif;
  while (defined (my $line = <DE>)) {
    chomp $line;
    my @tmp = split("\t", $line);
    # genes   logFC   AveExpr t       P.Value adj.P.Val       B
    if ($tmp[5] < 0.01) {       # FDR 1%
      my $exon = $tmp[0];
      $exon =~ s/\:/\-/;
      $DEsignif{$exon} = ();
    }
  }
  close DE;
  # print scalar keys %DEsignif, " DE exons\n";

  # generate a BED file of signficant DE orthoExons
  open(IN, '../orthoExon/orthoExons_hg19.bed' ) or die ("Cannot open file\n");
  open(OUT, '>', $filename ) or die ("Cannot open OUT file\n");
  while (defined (my $line = <IN>)) {
    chomp $line;
    my @tmp = split("\t", $line);
    if (exists $DEsignif{$tmp[3]}) {
      print OUT $line, "\n";
    }
  }
  close IN;
  close OUT;
}

# now we can launch bedtools
print "Querying DE exons file...\n";
open (BED, "intersectBed -u -a $DMRfile -b DE_exons/$DMRtype/$DMRcontrast.bed -f 0.2 |") or die ("Bedtools query failed\n");
while (defined (my $line = <BED>)) {
  chomp $line;
  my @tmp = split("\t", $line);

  # tag each position
  $DMRs{$tmp[0]}->{$tmp[1]}->{$tmp[2]}->{'Diff. expr. exons'} = 1;
}
close BED;

## Differentially splice exons
## Check if BED file was already generated
$filename = 'DS_exons/'.$DMRtype.'/'.$DMRcontrast.'.bed';
unless (-s $filename) {
  print "Generating BED files of significanly differentially spliced exons...\n";
  ## We use the same correspondance hash as DE exons for the file names);
  my $file;
  if ($DMRtype eq 'tissues') {
    $file = '../limma/differential_splicing/DMRs_comparison/'.$species.$correspondance{$tissue1}->{$tissue2}.'_exons.txt';
  }
  if ($DMRtype eq 'species') {
    $file = '../limma/differential_splicing/DMRs_comparison/'.$correspondance{$species1}->{$species2}.'_'.ucfirst($tissue).'_exons.txt';
  }
  open(DS, $file ) or die ("Cannot open DS file\n");
  # record significant DS exons
  my %DSsignif;
  while (defined (my $line = <DS>)) {
    chomp $line;
    my @tmp = split("\t", $line);
    # ExonID / gene / exon number  / log fold-change / t statistic / p-value / FDR
    if ($tmp[6] < 0.01) {       # FDR 1%
      my $exon = $tmp[0];
      $exon =~ s/\:/\-/;
      $DSsignif{$exon} = ();
    }
  }
  close DS;
  # print scalar keys %DSsignif, " DS exons\n";

  # generate a BED file of signficant DS orthoExons
  open(IN, '../orthoExon/orthoExons_hg19.bed' ) or die ("Cannot open file\n");
  open(OUT, '>', 'DS_exons/'.$DMRtype.'/'.$DMRcontrast.'.bed' ) or die ("Cannot open OUT file\n");
  while (defined (my $line = <IN>)) {
    chomp $line;
    my @tmp = split("\t", $line);
    if (exists $DSsignif{$tmp[3]}) {
      print OUT $line, "\n";
    }
  }
  close IN;
  close OUT;
}

# now we can launch bedtools
print "Querying differentially spliced exons file...\n";
open (BED, "intersectBed -u -a $DMRfile -b DS_exons/$DMRtype/$DMRcontrast.bed -f 0.2 |") or die ("Bedtools query failed\n");
while (defined (my $line = <BED>)) {
  chomp $line;
  my @tmp = split("\t", $line);

  # tag each position
  $DMRs{$tmp[0]}->{$tmp[1]}->{$tmp[2]}->{'Diff. used exons'} = 1;
}
close BED;


## Overlap with enhancers from FANTOM (human)
# Ubiquitous, Tissue-expressed and tissue-specific enhancers for the tissues of interest

my %enhancerFiles = (
  'Heart-specific enhancers' => 'UBERON_0000948_heart_differentially_expressed_enhancers.bed',
  'Heart enhancers' => 'UBERON:0000948_heart_expressed_enhancers.bed',
  'Lung-specific enhancers' => 'UBERON_0002048_lung_differentially_expressed_enhancers.bed',
  'Lung enhancers' => 'UBERON:0002048_lung_expressed_enhancers.bed',
  'Liver-specific enhancers' => 'UBERON_0002107_liver_differentially_expressed_enhancers.bed',
  'Liver enhancers' => 'UBERON:0002107_liver_expressed_enhancers.bed',
  'Kidney-specific enhancers' => 'UBERON_0002113_kidney_differentially_expressed_enhancers.bed',
  'Kidney enhancers' => 'UBERON:0002113_kidney_expressed_enhancers.bed',
  'Ubiquitous enhancers' => 'Ubiquitous_enhancers_S10.bed',
);

foreach my $enhancerFile (sort keys %enhancerFiles){
  # now we can launch bedtools
  print "Querying $enhancerFile file...\n";
  # Here overlap of 20% not required (because enhancers are short)
  open (BED, "tail -n +2 FANTOM_enhancers/". $enhancerFiles{$enhancerFile} ."| intersectBed -u -a $DMRfile -b - |") or die ("Bedtools query failed\n");
  while (defined (my $line = <BED>)) {
    chomp $line;
    my @tmp = split("\t", $line);

    # tag each position
    $DMRs{$tmp[0]}->{$tmp[1]}->{$tmp[2]}->{$enhancerFile} = 1;
  }
  close BED;
}


## Overlap of tDMRs with sDMRs
my %DMR2contrasts; # record the contrast for output file
if ($DMRtype eq 'tissues'){
  my @DMR2files = glob("../bsseq/DMRs/species/*_DMRs.txt");
  foreach my $DMR2file (@DMR2files){
    $DMR2file =~ m/\.\.\/bsseq\/DMRs\/species\/(.+)\_DMRs\.txt/;
    my $DMR2contrast = $1;
    next if ($DMR2contrast =~ m/\_overlap\_/); ## we do not want the txt files generated from overlap of different DMR file
    next unless ($DMR2contrast =~ m/$species/); ## we want only sDMRs involving the species of the tDMR
    next unless (($DMR2contrast =~ m/$tissue1/) or ($DMR2contrast =~ m/$tissue2/)); ## we want only sDMRs in the tissues involved in the tDMR

    # reformatting the contrats names for output
    my $contrast = $DMR2contrast;
    $contrast =~ s/(.+)Chimp/$1 vs\. Chimp/; # we do not want ot match the "Chimp" at beginning of string
    $contrast =~ s/Rhesus/ vs\. Rhesus/;
    $contrast =~ s/Specific/\-specific/;
    $contrast =~ s/\_/ in /;
    $contrast .= ' species-DMRs';
    $DMR2contrasts{$contrast} = ();
    print "Querying $contrast file...\n";

    # Here overlap of 20% is required
    open (BED, "intersectBed -u -a $DMRfile -b DMRs/species/".$DMR2contrast."_DMRs.bed -f 0.2 |") or die ("Bedtools query failed\n");
    while (defined (my $line = <BED>)) {
      chomp $line;
      my @tmp = split("\t", $line);

      # tag each position
      $DMRs{$tmp[0]}->{$tmp[1]}->{$tmp[2]}->{$contrast} = 1;
    }
    close BED;
  }
}


## Overlap of sDMRs with tDMRs
if ($DMRtype eq 'species'){
  my @DMR2files = glob("../bsseq/DMRs/tissues/*_DMRs.txt");
  foreach my $DMR2file (@DMR2files){
    $DMR2file =~ m/\.\.\/bsseq\/DMRs\/tissues\/(.+)\_DMRs\.txt/;
    my $DMR2contrast = $1;
    next if ($DMR2contrast =~ m/\_overlap\_/); ## we do not want the txt files generated from overlap of different DMR file
    next unless ($DMR2contrast =~ m/$tissue/); ## we want only tDMRs involving the tissue of the sDMR
    next unless (($DMR2contrast =~ m/$species1/) or ($DMR2contrast =~ m/$species2/)); ## we want only tDMRs in the species involved in the sDMR

    # reformatting the contrats names for output
    my $contrast = $DMR2contrast;
    $contrast =~ s/(heart|liver|lung)\_kidney/$1 vs\. kidney/;
    $contrast =~ s/(heart|liver)\_lung/$1 vs\. lung/;
    $contrast =~ s/(heart)\_liver/$1 vs\. liver/;
    $contrast =~ s/Specific/\-specific/;
    $contrast =~ s/\_/ /;
    $contrast .= ' tissue-DMRs';
    $DMR2contrasts{$contrast} = ();
    print "Querying $contrast file...\n";

    # Here overlap of 20% is required
    open (BED, "intersectBed -u -a $DMRfile -b DMRs/tissues/".$DMR2contrast."_DMRs.bed -f 0.2 |") or die ("Bedtools query failed\n");
    while (defined (my $line = <BED>)) {
      chomp $line;
      my @tmp = split("\t", $line);

      # tag each position
      $DMRs{$tmp[0]}->{$tmp[1]}->{$tmp[2]}->{$contrast} = 1;
    }
    close BED;
  }
}

# export matrix of features for all positions ######################################################
print "Exporting matrix of features...\n";
my $outfile = $DMRfile;
$outfile =~ s/\.bed//;
open(OUT, "| gzip - > ".$outfile."_contrast_specific_features.gz" ) or die ("Cannot open OUT file\n");
## headers
my @features = ('Diff. expr. exons', 'Diff. used exons',  keys %enhancerFiles, keys %DMR2contrasts);
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

