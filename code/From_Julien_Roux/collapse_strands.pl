#!/usr/bin/perl
use strict;
use warnings;
use Scalar::Util qw(looks_like_number);
$| = 1;

## Julien Roux
## June 13, 2013
## Create new bedGraph methylation summary file, with + and - strand collapsed (CpG methylation is usually symetrical)
## The 16g asked with qsub should be enough memory for this script...
## Nov 29, 2013: Even if there is nothing to collapse, when there is data only on the - strand, coordinates of the + strand should be the output. Otherwise, when we load multiple samples with bsseq, there will be samples with data on the + and samples with data on the - strand: the data will not be collapsed.
## /!\ so even if column 7 indicates -, the coordinates are on the + strand! Column 7 indicates where the original data come from
## Aug 27: There was an additional (small) mistake, see below
die("Beware. This script was modified to correct a small mistake that was reporting Cs on the - strand. The output will differ from the files used so far. Maybe the old files should be backed-up before!\n")

scalar @ARGV >= 3 or die<<USAGE;
Usage:
perl collapse_strand.pl <folder/sample> <compressed input bedGraph file> <compressed input CpG_OT file>
USAGE
# e.g. perl collapse_strand ./H1H1 H1H1_trimmed.fq.gz_bismark.deduplicated.bedGraph.gz CpG_OT_H1H1_trimmed.fq.gz_bismark.deduplicated.txt.gz

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

print "Reading $bedgraph file and collapsing into $sample.bedGraph.gz...\n";
open(IN, "zcat $bedgraph |") or die "Can't read bedgraph file";
open(OUT, "| gzip -9 - > $sample.bedGraph.gz") or die "Can't open output file";
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
      my $unmeth = $previous[5] + $tmp[5];
      my $percent = $meth / ($meth + $unmeth) * 100;
      print OUT "$previous[0]\t$previous[1]\t$previous[2]\t$percent\t$meth\t$unmeth\t+/-\tmerged\n";

      ## Now we don't want the current line to be used again because it was collapsed: we read a new line
      $line = <IN>;
      chomp($line) if (defined $line);
      $i++;
    }
    else {
      ## the @previous read was on - strand: don't merge
      # print OUT $previous_line, "\t-\tnot_merged\n";

      ## This is wrong! as below, we should have subtracted 1 to print the coordinates on the top strand!
      ## This happens when 2 CpGs are follwoing each other, and the first one has only data on the - strand: should be pretty rare
      ## It is not an error strictly speaking, but some data could have been gathered and were not maybe. Also it can be confusing because these positions will be in teh final bedgraph files although they are not Cs on + strand... I'll make a script to identify these problematic positions and possibly remove them from final analyses...
      ## corrected output:
      print OUT "$previous[0]\t", $previous[1]-1, "\t", $previous[2]-1, "\t$previous[3]\t$previous[4]\t$previous[5]\t-\tnot_merged\n";
    }
  }
  else {
    if (exists $top_cpgs{$previous[0]}->{$previous[1]}){
      ## the @previous read was on + strand: don't touch
      print OUT $previous_line, "\t+\tnot_merged\n";
    }
    else {
      ## the @previous read was on - strand: we substract 1 to the coordinates to print the coordinates on the CpG on the + strand
      print OUT "$previous[0]\t", $previous[1]-1, "\t", $previous[2]-1, "\t$previous[3]\t$previous[4]\t$previous[5]\t-\tnot_merged\n";
    }
  }
  $previous_line = $line;
}
## in most cases, the last line is not written down (except if it is collapsed with the previous line. In this case, $line and $previous_line will be undef
if (defined $previous_line){
  my @previous = split("\t", $previous_line);
  if (exists $top_cpgs{$previous[0]}->{$previous[1]}){
    ## the @previous read was on + strand
    print OUT $previous_line, "\t+\tnot_merged\n";
  }
  else {
    ## the @previous read was on - strand: we substract 1 to the coordinates to print the coordinates on the CpG on the + strand
    print OUT "$previous[0]\t", $previous[1]-1, "\t", $previous[2]-1, "\t$previous[3]\t$previous[4]\t$previous[5]\t-\tnot_merged\n";
  }
}
close IN;
close OUT;

print "$i CpGs parsed, $j CpGs in final collapsed file.\n";
exit;
