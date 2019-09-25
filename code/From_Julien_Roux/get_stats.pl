#!/usr/bin/perl
use strict;
use warnings;

## Julien Roux
## June 7th, 2012
# Extract basic stats from finished bismark run: did everything go fine?

my $outdir = $ARGV[0];
my $sample = $ARGV[1];

print "Extracting basic statistics for the run...\n";

## TO DO? Extract stats from cutadapt?

# total number of reads processed (prep_reads.info)
open(IN, '<', "$outdir/$sample\_trimmed.fq.gz_bismark_SE_report.txt") or die ("Cannot open mapping report\n");
my $total;
my $efficiency;
my $top;
my $bottom;
my $methylated;
while (defined (my $line = <IN>)) {
  chomp $line;
  if ($line =~ m/Sequences analysed in total\:\t(\d+)/){
    $total = $1;
  }
  if ($line =~ m/Mapping efficiency\:\t(\d+\.\d+\%)/){
    $efficiency = $1;
  }
  if ($line =~ m/CT\/CT\:\t(\d+)\t\(\(converted\) top strand\)/){
    $top = $1;
  }
  if ($line =~ m/CT\/GA\:\t(\d+)\t\(\(converted\) bottom strand\)/){
    $bottom = $1;
  }
  if ($line =~ m/C methylated in CpG context\:\t(\d+\.\d+\%)/){
    $methylated = $1;
  }
}
close IN;

# time of execution
# my $time;
# open(IN, '<', ) or die ("Cannot open file\n");
# while (defined (my $line = <IN>)) {
#   chomp $line;
#   if ($line =~ m//){
#     $time = $1;
#   }
# }
# close IN;

$total = commify($total);
$top = commify($top);
$bottom = commify($bottom);

# output
open(OUT, '>>', $outdir.$sample.'.report') or die 'Cannot open OUT file';
print OUT "\nBasic statistics from bismark mapping:\n  Total number of reads processed: $total\n  Mapping efficiency: $efficiency\n  Number of reads mapped to top strand: $top\n  Number of reads mapped to bottom strand: $bottom\n  Pourcentage of Cytosines methylated in CpG context: $methylated\n";
close OUT;

sub commify {
  my $text = reverse $_[0];
  $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
  return scalar reverse $text
}

exit;

