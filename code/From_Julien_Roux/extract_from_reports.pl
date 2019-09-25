#!/usr/bin/perl
use strict;
use warnings;

## Julien Roux
## Mar 28, 2014
## Extract all stats from .reports files (generated after each tophat run) and put them in a tab-delimited file reports.txt

## record all .report files to consider
my %samples;
foreach my $file ( glob("./*/*.report") ){
  print "$file\n";
  $file =~ m/\.\/(.+)\/.+\.report/;
  my $sample = $1;

  # open the .report file
  open(IN, $file) or die ("Cannot open report file\n");
  my $flag = 0;
  while (my $line = <IN>){
    chomp($line);
    if (($line =~ m/Total\snumber\sof\sreads\sprocessed\:\s(\d+)/) and ($flag eq 0)) {
      $samples{$sample}->{'total_sequenced'} = $1;
    }
    if ($line =~ m/Percentage\sof\sbps\strimmed\:\s(.+)\%/) {
      $samples{$sample}->{'percent_trimmed'} = $1;
    }
    if ($line =~ m/Number\sof\sreads\sremoved\sbecause\sthey\swere\sshorter\sthan\s20bp\:\s(\d+)/) {
      $samples{$sample}->{'number_reads_removed'} = $1;
    }
    if ($line =~ m/Minimum\sread\slength\:\s(\d+)/) {
      $samples{$sample}->{'min_length'} = $1;
    }
    if ($line =~ m/Maximum\sread\slength\:\s(\d+)/) {
      $samples{$sample}->{'max_length'} = $1;
    }
    if ($line =~ m/Basic statistics from tophat mapping/) {
      $flag = 1;
    }
    if (($line =~ m/Total\snumber\sof\sreads\sprocessed\:\s(\d+)/) and ($flag eq 1)) {
      $samples{$sample}->{'total_reads_tophat'} = $1;
    }
    if ($line =~ m/Total\snumber\sof\sreads\smapped\:\s(\d+)/) {
      $samples{$sample}->{'total_mapped'} = $1;
    }
    if ($line =~ m/Percentage\sof\sreads\smapped\:\s(.+)\%/){
      $samples{$sample}->{'percent_mapped'} = $1;
    }
    if ($line =~ m/Percentage\sof\smapped\sreads\soverlapping\sa\sjunction\:\s(.+)\%/) {
      $samples{$sample}->{'total_overlap_1junc'} = $1;
    }
    if ($line =~ m/Percentage\sof\smapped\sreads\soverlapping\smore\sthan\sone\sjunction\:\s(.+)\%/) {
      $samples{$sample}->{'total_overlap_>1junc'} = $1;
    }
    if ($line =~ m/Number\sof\sjunctions\:\s(\d+)/) {
      $samples{$sample}->{'num_junc'} = $1;
    }
    if ($line =~ m/Number\sof\sinsertions\:\s(\d+)/) {
      $samples{$sample}->{'num_ins'} = $1;
    }
    if ($line =~ m/Number\sof\sdeletions\:\s(\d+)/) {
      $samples{$sample}->{'num_del'} = $1;
    }
    if ($line =~ m/Time\sfor\smapping\:\s(.+)/) {
      $samples{$sample}->{'time'} = $1;
    }
  }
  close IN;
}

#write in reports.txt file
open(OUT, ">reports.txt") or die ("Cannot open reports.txt file\n");
## header
print OUT "SampleID\tTotal number of reads sequenced\tPercentage of bps trimmed (quality and adapters)\tNumber of reads shorter than 20bp removed\tMinimum read length after trimming\tMaximum read length after trimming\tTotal number of reads processed in tophat\tTotal number of reads mapped\tPercentage of reads mapped\tPercentage of mapped reads overlapping a junction\tPercentage of mapped reads overlapping more than one junction\tNumber of junctions\tNumber of insertions\tNumber of deletions\tTime for mapping\n";
foreach my $sample (sort keys %samples){
  print OUT "$sample\t", $samples{$sample}->{'total_sequenced'}, "\t", $samples{$sample}->{'percent_trimmed'}*100, "\t", $samples{$sample}->{'number_reads_removed'}, "\t", $samples{$sample}->{'min_length'}, "\t", $samples{$sample}->{'max_length'}, "\t", $samples{$sample}->{'total_reads_tophat'}, "\t", $samples{$sample}->{'total_mapped'}, "\t", $samples{$sample}->{'percent_mapped'}, "\t", $samples{$sample}->{'total_overlap_1junc'}, "\t", $samples{$sample}->{'total_overlap_>1junc'}, "\t", $samples{$sample}->{'num_junc'}, "\t", $samples{$sample}->{'num_ins'}, "\t", $samples{$sample}->{'num_del'}, "\t", $samples{$sample}->{'time'}, "\n";
}
close OUT;
exit;


