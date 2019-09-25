#!/usr/bin/perl
use strict;
use warnings;
use Scalar::Util qw(looks_like_number);
$| = 1;

## Julien Roux
## Aug 29, 2014
# In the collapse_strands.pl script there was a small error:
# When 2 CpGs are following each other, and the first one has only data on the - strand, we should have subtracted 1 to print the coordinates on the top strand, but we didn't. This should be pretty rare (there should be 2 consecutive CpG sites, with data only on - strand for the first one, and some data on the + strand for the second one).
# This script (modified from collapse_strands.pl) is here to list the problematic positions
# All problematic positions are on the - strand. When technical samples are gathered, there could be some missing coverage on positions on + strand: see ~/methylation/bsseq/combine_data_perl/check_strand_evidences.pl 

# (Another way to detect (some of) these problematic positions with R: check if CpG positions have consecutive start positions (should not happen). This identifies 1629 positions out of 5228469 in positions conserved between the 2 species. In comparison, a difference of 2 (consecutive CpG) identifies 131812 sites.)

scalar @ARGV >= 3 or die<<USAGE;
Usage:
perl identify_strand_problems.pl <folder/sample> <compressed input bedGraph file> <compressed input CpG_OT file>
USAGE
# e.g. perl identify_strand_problems.pl ./H1H1 H1H1_trimmed.fq.gz_bismark.deduplicated.bedGraph.gz CpG_OT_H1H1_trimmed.fq.gz_bismark.deduplicated.txt.gz

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
  ## chromosome -> coordinate of C of CpGs (they are 1 based so we subtract 1)
  $top_cpgs{$tmp[2]}->{$tmp[3]-1}++
}
close IN;

print "Reading $bedgraph file and reporting erroneous positions into $sample.problematic_positions.bedGraph.gz...\n";
open(IN, "zcat $bedgraph |") or die "Can't read bedgraph file";
open(OUT, "| gzip -9 - > $sample.problematic_positions.bedGraph.gz") or die "Can't open output file";
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
    ## if the previous line was on top strand: no problem!
    if (exists $top_cpgs{$previous[0]}->{$previous[1]}){
      ## Now we don't want the current line to be used again because it was collapsed: we read a new line
      $line = <IN>;
      chomp($line) if (defined $line);
      $i++;
    }
    else {
      # Here was the problem!
      # report this position in the OUT file:
      # chr / position
      print OUT "$previous[0]\t", $previous[1], "\t-\n";
    }
  }
  # else: no problem

  $previous_line = $line;
}
close IN;
close OUT;

print "$i CpGs parsed, $j CpGs in final collapsed file.\n";
exit;
