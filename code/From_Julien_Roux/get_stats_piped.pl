#!/usr/bin/perl
use strict;
use warnings;

## Julien Roux
## April 4, 2012
# Extract basic stats from finished tophat run

# to improve memory efficiency, the BAM file is read by direct pipe into this script
# command line to add to launch_tophat.pl:
# samtools view $outdir/accepted_hits.bam | perl get_stats_piped.pl $outdir $sample

my $path = $ARGV[0];
if ($path !~ m/\/$/){
  $path .= '/';
}
my $sample = $ARGV[1];

print "Extracting basic statistics for the run...\n";

# total number of mapped reads
my $mapped = 0;
my %count_overlap;

## input from samtools view (line by line)
while (<STDIN>) {
  my $line = $_;
  chomp($line);
  $mapped++; ## one more read mapped

  # reads mapping to a junction: these have a N in their cigar line (e.g., 87M1000N20M)
  # we want to output the number of reads which are mapped to 0 junctions, at least 1 junc, and more than 1 junc
  # We first extract the string of number of Ns per line
  my @tmp = split("\t", $line);
  my $cigar = $tmp[5];
  my $overlap = $cigar =~ tr/N//;
  $count_overlap{$overlap}++;
}

my $no_overlap = $count_overlap{"0"};
my $at_least_one_overlap = 0;
foreach my $overlap (keys %count_overlap){
  #print "$overlap = $count_overlap{$overlap}\n";
  next if ($overlap eq "0");
  $at_least_one_overlap += $count_overlap{$overlap};
}
my $more_than_one_overlap = $at_least_one_overlap - $count_overlap{"1"};
#print "0 $no_overlap\n";
#print "1+ $at_least_one_overlap\n";
#print "2+ $more_than_one_overlap\n";

# number of junctions/insertions/deletions found
# wc without header
my $junctions = `tail -n+2 $path/junctions.bed | wc -l`;
my $insertions = `tail -n+2 $path/insertions.bed | wc -l`;
my $deletions = `tail -n+2 $path/deletions.bed | wc -l`;
chomp($junctions);
chomp($insertions);
chomp($deletions);

# total number of reads trimmed
my $sequenced_reads;
my $percent_bps_trimmed;
my $bps_trimmed = 0;
my $bps_total;
my $removed_reads;
open(IN, '<', $path.$sample.'.fastq.gz_trimming_report.txt') or die ("Cannot open trimming report file\n");
while (defined (my $line = <IN>)) {
  chomp $line;
  if ($line =~ m/(\d+) sequences processed in total/){
    $sequenced_reads = $1;
  }
  if ($line =~ m/Sequences removed because they became shorter than the length cutoff of \d+ bp:\s+(\d+)\s\(\d+\.\d+\%\)/){
    $removed_reads = $1;
  }
  if ($line =~ m/Quality-trimmed bases:\s+(\d+)\sbp/){
    $bps_trimmed += $1;
  }
  if ($line =~ m/Trimmed bases:\s+(\d+)\sbp/){
    $bps_trimmed += $1;
  }
  if ($line =~ m/Processed bases:\s+(\d+)\s/){
    $bps_total += $1;
  }
}
close IN;

# total number of reads processed (prep_reads.info)
my $min_length;
my $max_length;
my $total_reads;
open(IN, '<', $path.'prep_reads.info') or die ("Cannot open prep_reads.info file\n");
while (defined (my $line = <IN>)) {
  chomp $line;
  if ($line =~ m/min\_read\_len\=(\d+)/){
    $min_length = $1;
  }
  if ($line =~ m/max\_read\_len\=(\d+)/){
    $max_length = $1;
  }
  if ($line =~ m/reads\_out\=(\d+)/){
    $total_reads = $1;
  }
}
close IN;

# time of execution (std.err)
my $time;
open(IN, '<', $path.$sample.'.err') or die ("Cannot open std.err file\n");
while (defined (my $line = <IN>)) {
  chomp $line;
  if ($line =~ m/Run\scomplete\:\s(\d+\:\d+\:\d+)\selapsed/){
    $time = $1;
  }
}
close IN;

# output
open(OUT, '>>', $path.$sample.'.report') or die 'Cannot open OUT file';
print OUT "Basic statistics from trimGalore adapter and quality trimming:
  Total number of reads processed: $sequenced_reads
  Percentage of bps trimmed: ", $bps_trimmed/$bps_total*100, "%
  Number of reads removed because they were shorter than 20bp: $removed_reads
Basic statistics from tophat mapping:
  Minimum read length: $min_length
  Maximum read length: $max_length
  Total number of reads processed: $total_reads
  Total number of reads mapped: $mapped
  Percentage of reads mapped: ", $mapped/$total_reads*100, "\%
  Percentage of mapped reads overlapping a junction: ", $at_least_one_overlap/$mapped*100, "\%
  Percentage of mapped reads overlapping more than one junction: ", $more_than_one_overlap/$mapped*100, "\%
  Number of junctions: $junctions
  Number of insertions: $insertions
  Number of deletions: $deletions
  Time for mapping: $time\n";
close OUT;
exit;

