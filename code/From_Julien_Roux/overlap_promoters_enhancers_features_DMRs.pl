#!/usr/bin/perl
use strict;
use warnings;

## Julien Roux
## Sep 1 2015
## intersect file of DMRs (or file of random DMRs) to these features
#          Enhancers from FANTOM (human) in specific tissues (already done above but here again to be able to compare to promoters)
#          Promoters from FANTOM (human) in specific tissues
#          Enhancers from VISTA (human) in heart / anywhere during dvlpmt
#          ENCODE promoters and enhancers from cell lines: Hepg2 (liver), Hsmm (muscle), Nhlf (lung fibroblast)
#          NIH Roadmap epigenomics chromatin states


scalar @ARGV >= 1 or die<<USAGE;
Usage:
perl overlap_promoters_enhancers_features_DMRs.pl <DMR BED file with full path>
USAGE
my $DMRfile = $ARGV[0];

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

## Overlap with enhancers from FANTOM (human)
# Ubiquitous, Tissue-expressed and tissue-specific enhancers for the tissues of interest
my %FANTOMenhancerFiles = (
  'FANTOM Heart-specific enhancers' => 'UBERON_0000948_heart_differentially_expressed_enhancers.bed',
  'FANTOM Heart enhancers' => 'UBERON:0000948_heart_expressed_enhancers.bed',
  'FANTOM Lung-specific enhancers' => 'UBERON_0002048_lung_differentially_expressed_enhancers.bed',
  'FANTOM Lung enhancers' => 'UBERON:0002048_lung_expressed_enhancers.bed',
  'FANTOM Liver-specific enhancers' => 'UBERON_0002107_liver_differentially_expressed_enhancers.bed',
  'FANTOM Liver enhancers' => 'UBERON:0002107_liver_expressed_enhancers.bed',
  'FANTOM Kidney-specific enhancers' => 'UBERON_0002113_kidney_differentially_expressed_enhancers.bed',
  'FANTOM Kidney enhancers' => 'UBERON:0002113_kidney_expressed_enhancers.bed',
  'FANTOM Ubiquitous enhancers' => 'Ubiquitous_enhancers_S10.bed',
);

foreach my $enhancerFile (sort keys %FANTOMenhancerFiles){
  # now we can launch bedtools
  print "Querying $enhancerFile file...\n";
  # Here overlap of 20% not required (because enhancers are short)
  open (BED, "tail -n +2 FANTOM_enhancers/". $FANTOMenhancerFiles{$enhancerFile} ."| intersectBed -u -a $DMRfile -b - |") or die ("Bedtools query failed\n");
  while (defined (my $line = <BED>)) {
    chomp $line;
    my @tmp = split("\t", $line);

    # tag each position
    $DMRs{$tmp[0]}->{$tmp[1]}->{$tmp[2]}->{$enhancerFile} = 1;
  }
  close BED;
}

## Overlap with promoters from FANTOM (human)
# Top 1000 ubiquitous and diff expr. in tissue of interest of now (asked expressed in tissue of interest)
my %FANTOMpromoterFiles = (
  'FANTOM Heart-specific promoters' => 'UBERON:0000948_heart_differentially_expressed_promoters.txt',
  #'FANTOM Heart promoters' => '',
  'FANTOM Lung-specific promoters' => 'UBERON:0002048_lung_differentially_expressed_promoters.txt',
  #'FANTOM Lung promoters' => '',
  'FANTOM Liver-specific promoters' => 'UBERON:0002107_liver_differentially_expressed_promoters.txt',
  #'FANTOM Liver promoters' => '',
  'FANTOM Kidney-specific promoters' => 'UBERON:0002113_kidney_differentially_expressed_promoters.txt',
  #'FANTOM Kidney promoters' => '',
  'FANTOM Ubiquitous promoters' => 'top1000_ubiquitous_promoters.txt',
);

foreach my $promoterFile (sort keys %FANTOMpromoterFiles){
  # now we can launch bedtools
  print "Querying $promoterFile file...\n";
  # Here overlap of 20% not required (because promoters are short)
  open (BED, "cat FANTOM_promoters/". $FANTOMpromoterFiles{$promoterFile} ." | cut -f1 | cut -f1,2 -d\':\' --output-delimiter=\'\t\' | cut -f1,3 -d\'.\' --output-delimiter=\'\t\' | cut -f1,2 -d\',\' --output-delimiter=\'\t\' | intersectBed -u -a $DMRfile -b - |") or die ("Bedtools query failed\n");
  while (defined (my $line = <BED>)) {
    chomp $line;
    my @tmp = split("\t", $line);
    # tag each position
    $DMRs{$tmp[0]}->{$tmp[1]}->{$tmp[2]}->{$promoterFile} = 1;
  }
  close BED;
}


## Overlap with enhancers from VISTA (human)
# all detected enhancers, heart enhancers
my %VISTAenhancerFiles = (
  'VISTA enhancers' => 'all.bed',
  'VISTA Heart enhancers' => 'heart.bed',
);

foreach my $enhancerFile (sort keys %VISTAenhancerFiles){
  # now we can launch bedtools
  print "Querying $enhancerFile file...\n";
  # Here overlap of 20% not required (because enhancers are short)
  open (BED, "intersectBed -u -a $DMRfile -b VISTA_enhancers/". $VISTAenhancerFiles{$enhancerFile}. "|") or die ("Bedtools query failed\n");
  while (defined (my $line = <BED>)) {
    chomp $line;
    my @tmp = split("\t", $line);

    # tag each position
    $DMRs{$tmp[0]}->{$tmp[1]}->{$tmp[2]}->{$enhancerFile} = 1;
  }
  close BED;
}

## ENCODE enhancers
## Problem: no end position for these files. For now we just look at the intersection with start coordinates
my %ENCODEenhancerFiles = (
  'ENCODE Hepg2 enhancers' => 'predicted_enhancer_human/Hepg2.enhancer.final.bed',
  'ENCODE Hsmm enhancers' => 'predicted_enhancer_human/Hsmm.enhancer.final.bed',
  'ENCODE Nhlf enhancers' => 'predicted_enhancer_human/Nhlf.enhancer.final.bed',
);

foreach my $enhancerFile (sort keys %ENCODEenhancerFiles){
  # now we can launch bedtools
  print "Querying $enhancerFile file...\n";
  # Here overlap of 20% not required (because enhancers are short)
  open (BED, "intersectBed -u -a $DMRfile -b ENCODE_enhancers/". $ENCODEenhancerFiles{$enhancerFile} ." |") or die ("Bedtools query failed\n");
  while (defined (my $line = <BED>)) {
    chomp $line;
    my @tmp = split("\t", $line);

    # tag each position
    $DMRs{$tmp[0]}->{$tmp[1]}->{$tmp[2]}->{$enhancerFile} = 1;
  }
  close BED;
}

## ENCODE promoters
my %ENCODEpromoterFiles = (
  'ENCODE Hepg2 promoters' => 'predicted_promoter_human/Hepg2.H3k4me3_peaks.xls.config.en.2fold.1rpkm.bed',
  'ENCODE Hsmm promoters' => 'predicted_promoter_human/Hsmm.H3k4me3_peaks.xls.config.en.2fold.1rpkm.bed',
  'ENCODE Nhlf promoters' => 'predicted_promoter_human/Nhlf.H3k4me3_peaks.xls.config.en.2fold.1rpkm.bed',
);

foreach my $promoterFile (sort keys %ENCODEpromoterFiles){
  # now we can launch bedtools
  print "Querying $promoterFile file...\n";
  # Here overlap of 20% not required (because promoters are short)
  open (BED, "intersectBed -u -a $DMRfile -b ENCODE_enhancers/". $ENCODEpromoterFiles{$promoterFile} ." |") or die ("Bedtools query failed\n");
  while (defined (my $line = <BED>)) {
    chomp $line;
    my @tmp = split("\t", $line);

    # tag each position
    $DMRs{$tmp[0]}->{$tmp[1]}->{$tmp[2]}->{$promoterFile} = 1;
  }
  close BED;
}


## NIH Roadmap epigenomics chromatin state
my %RoadmapFiles = (
  'Liver' => 'E066',
  'Fetal Heart' => 'E083',
  'Heart Right Atrium' => 'E104',
  'Heart Right Ventricle' => 'E105',
  'Heart Left Ventricle' => 'E095',
  'Fetal Kidney' => 'E086',
  'Fetal Lung' => 'E088',
  'Lung' => 'E096',
);

my %chromatinStates = (
  '1_TssA' => 'Active TSS',
  '2_TssAFlnk' => 'Flanking Active TSS',
  '3_TxFlnk' => 'Transcr. at gene 5\' and 3\'',
  '4_Tx' => 'Strong transcription',
  '5_TxWk' => 'Weak transcription',
  '6_EnhG' => 'Genic enhancers',
  '7_Enh' => 'Enhancers',
  '8_ZNF/Rpts' => 'ZNF genes & repeats',
  '9_Het' => 'Heterochromatin',
  '10_TssBiv' => 'Bivalent/Poised TSS',
  '11_BivFlnk' => 'Flanking Bivalent TSS/Enh',
  '12_EnhBiv' => 'Bivalent Enhancer',
  '13_ReprPC' => 'Repressed PolyComb',
  '14_ReprPCWk' => 'Weak Repressed PolyComb',
  '15_Quies' => 'Quiescent/Low',
);

foreach my $file (sort keys %RoadmapFiles){
  foreach my $state (sort keys %chromatinStates){
  # now we can launch bedtools
    print "Querying $state state in $file file...\n";
    open (BED, "zcat RoadmapEpigenomics_chromatin_states/". $RoadmapFiles{$file} ."_15_coreMarks_mnemonics.bed.gz | grep $state | intersectBed -u -a $DMRfile -b - |") or die ("Bedtools query failed\n");
    while (defined (my $line = <BED>)) {
      chomp $line;
      my @tmp = split("\t", $line);

      # tag each position
      $DMRs{$tmp[0]}->{$tmp[1]}->{$tmp[2]}->{$file.' / '.$chromatinStates{$state}} = 1;
    }
    close BED;
  }
}



# export matrix of features for all positions ######################################################
print "Exporting matrix of features...\n";
my $outfile = $DMRfile;
$outfile =~ s/\.bed//;
open(OUT, "| gzip - > ".$outfile."_promoters_enhancers_features.gz" ) or die ("Cannot open OUT file\n");
## headers
my @features = (keys %FANTOMenhancerFiles,  keys %FANTOMpromoterFiles, keys %VISTAenhancerFiles,  keys %ENCODEenhancerFiles, keys %ENCODEpromoterFiles);
foreach my $file (sort keys %RoadmapFiles){
  foreach my $state (sort keys %chromatinStates){
    push(@features, $file.' / '.$chromatinStates{$state});
  }
}

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

