#!/usr/bin/perl
use strict;
use warnings;
use Scalar::Util qw(looks_like_number);
$| = 1;

## Julien Roux
## Aug 27, 2014
# Create new bedGraph methylation summary file, with + and - strand collapsed (CpG methylation is usually symetrical)
# see collapse_strand.pl from which this script was modified

# The aim of this script is to keep track of the number of unmethylated reads on + and - strands
# Because we suspect that CpG -> TpG SNPs could be mistaken to unmethylated positions. This bias if present should be observed only on 1 strand, so we can filter the sites with evidence for unmethylated reads on only 1 strand and see (with PCA) if the problem persists. Notably this could be responsible for global %methylation differences between species.

## cartoon to clarify problem:
## original sequence
# top strand    CpG (ref genome)
# bottom strand GpC

## converted sequence
# top strand    TpG -> mapped to TpG (T to C change)
# bottom strand GpT -> mapped to CpA (G to A change)

## CpG -> TpG SNP on top strand
# top strand    TpG -> mapped to TpG (looks like converted, i.e. unmethylated read)
# bottom strand ApC -> mapped to GpC (with mismatch, looks like methylated) or not mapped

## CpG -> TpG SNP on bottom strand
# top strand    CpA -> mapped to CpG (with mismatch, looks like methylated) or not mapped
# bottom strand GpT -> mapped to GpT (looks like converted, i.e. unmethylated read)

## Output: similar to ~/Methylation/bismakr/collapse_strands.pl, but the number of unmethylated reads should be split into + and - strands
## Note: to be consistent with the mistake on collapse_strands.pl (some Cs on - strand were reported), we report the same mistake. The number of lines in output files should be the same


scalar @ARGV >= 3 or die<<USAGE;
Usage:
perl collapse_evidence_strand.pl <folder/sample> <compressed input bedGraph file> <compressed input CpG_OT file>
USAGE
# e.g. perl collapse_evidence_strand ./H1H1 H1H1_trimmed.fq.gz_bismark.deduplicated.bedGraph.gz CpG_OT_H1H1_trimmed.fq.gz_bismark.deduplicated.txt.gz

my $sample = $ARGV[0];
my $bedgraph = $ARGV[1];
my $cpg = $ARGV[2];

## Read CpG_OT file to know which CpG were encountered on top strand (keep in hash in memory)
print "Reading $cpg to know which CpGs are on + strand...\n";
open(IN, "zcat $cpg |") or die "Can't read CpG file";
my $header = <IN>;
my %top_cpgs;
my $i = 0;
while (defined (my $line = <IN>)) {
  $i++;
  ## print the script progress
  if ($i % 1000000 eq 0) {
    print "\t".$i." lines parsed\n";
  }

  chomp $line;
  my @tmp = split("\t", $line);

  ## remove this once problem identified?
  if (!looks_like_number($tmp[3])){
    print "The following line seems to have a problem in file $cpg: $line\n";
  }
  ## chromosome -> coordinate of C of CpGs (1 based)
  $top_cpgs{$tmp[2]}->{$tmp[3]-1}++
}
close IN;

print "Reading $bedgraph file and collapsing into $sample.unmethylated_evidence.bedGraph.gz...\n";
open(IN, "zcat $bedgraph |") or die "Can't read bedgraph file";
open(OUT, "| gzip -9 - > $sample.unmethylated_evidence.bedGraph.gz") or die "Can't open output file";
my $previous_line = <IN>;
chomp($previous_line);
$i = 0;
my $j = 0;
while (defined (my $line = <IN>)) {
  $i++; ## the counters may be wrong by 1, use them just as indicator
  $j++;
  ## print the script progress
  if ($i % 1000000 eq 0) {
    print "\t".$i." CpGs parsed\n";
  }

  chomp $line;
  my @previous = split("\t", $previous_line);
  my @tmp = split("\t", $line);

  ## if coordinate of line follows the one of previous line
  ## if previous chr = chr and previous start+1 = start
  if (($previous[0] eq $tmp[0]) and ($previous[1]+1 eq $tmp[1])){
    ## if the previous line was on top strand (here the coordinates are 0-based)
    if (exists $top_cpgs{$previous[0]}->{$previous[1]}){
      ## we can collapse reads and print with @previous coordinates
      my $meth = $previous[4] + $tmp[4];
      my $unmeth_top = $previous[5];
      my $unmeth_bottom = $tmp[5];
      my $percent = $meth / ($meth + $unmeth_top + $unmeth_bottom) * 100;
      print OUT "$previous[0]\t$previous[1]\t$previous[2]\t$percent\t$meth\t$unmeth_top\t$unmeth_bottom\t+/-\tmerged\n";

      ## Now we don't want the current line to be used again because it was collapsed: we read a new line
      $line = <IN>;
      chomp($line) if (defined $line);
      $i++;
    }
    else {
      ## the @previous read was on - strand: don't merge
      ## number of unmethylated reads on top strand = 0

      ## Note: in the collpase_strands.pl there was an error, we forgot to subtract 1 to these coordinates
      ## see identify_strand_problems.pl
      ## The correct code should be:
      # print OUT "$previous[0]\t", $previous[1]-1, "\t", $previous[2]-1, "\t$previous[3]\t$previous[4]\t0\t$previous[5]\t-\tnot_merged\n";
      ## But it's preferable to be consistent with previous files (the probematic positions are flaggd ith other script)
      print OUT "$previous[0]\t$previous[1]\t$previous[2]\t$previous[3]\t$previous[4]\t0\t$previous[5]\t-\tnot_merged\n";
    }
  }
  else {
    if (exists $top_cpgs{$previous[0]}->{$previous[1]}){
      ## the @previous read was on + strand: don't merge
      ## number of unmethylated reads on bottom strand = 0
      print OUT "$previous[0]\t$previous[1]\t$previous[2]\t$previous[3]\t$previous[4]\t$previous[5]\t0\t+\tnot_merged\n";
    }
    else {
      ## the @previous read was on - strand: we substract 1 to the coordinates to print the coordinates on the CpG on the + strand
      ## number of unmethylated reads on top strand = 0
      print OUT "$previous[0]\t", $previous[1]-1, "\t", $previous[2]-1, "\t$previous[3]\t$previous[4]\t0\t$previous[5]\t-\tnot_merged\n";
    }
  }
  $previous_line = $line;
}
## in most cases, the last line is not written down (except if it is collapsed with the previous line. In this case, $line and $previous_line will be undef
if (defined $previous_line){
  my @previous = split("\t", $previous_line);
  if (exists $top_cpgs{$previous[0]}->{$previous[1]}){
    ## the @previous read was on + strand
    print OUT "$previous[0]\t$previous[1]\t$previous[2]\t$previous[3]\t$previous[4]\t$previous[5]\t0\t+\tnot_merged\n";
  }
  else {
    ## the @previous read was on - strand: we substract 1 to the coordinates to print the coordinates on the CpG on the + strand
    print OUT "$previous[0]\t", $previous[1]-1, "\t", $previous[2]-1, "\t$previous[3]\t$previous[4]\t0\t$previous[5]\t-\tnot_merged\n";
  }
}
close IN;
close OUT;

print "$i CpGs parsed, $j CpGs in final collapsed file.\n";
exit;
