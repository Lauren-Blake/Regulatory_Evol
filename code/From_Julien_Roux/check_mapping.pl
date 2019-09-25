#!/usr/bin/perl
use strict;
use warnings;

## Julien Roux
## Feb 11, 2014
## Launch a series of checks to be sure that samples were processed correctly
## Sometines corrupted or incomplete BAM files are generated for example

scalar @ARGV >= 1 or die<<USAGE;
Usage:
perl check_mapping.pl <Flow cell_sample>
USAGE

my $directory = $ARGV[0];
$directory = substr($directory, 0, -1) if ($directory =~ m/\/$/);
$directory =~ m/(Methylation\_\d+)\_(.+)/;
my $flowcell = $1;
my $sample = $2;
my $condition;
my $warning = 0;

print "Checking if sample was already treated...\n";
open(IN, '<', "checked_samples.txt") or die ("Cannot open report file\n");
my $first = 0; # record if this is the first sample checked, so that headers are printed in the checked_samples.txt file
while (defined (my $line = <IN>)) {
  $first++;
  chomp $line;
  my @tmp = split("\t", $line);
  if (($tmp[0] eq $flowcell) and ($tmp[1] eq $sample)){
    die '/!\ This sample was previously checked!', "\n";
  }
}

## check that the index number is not in disagreement with other technical replicates
print "Reading information file...\n";
open(IN, '<', "../raw_data/SamplesDirectories.txt") or die ("Cannot open info file\n");
my $header = <IN>;
my %index;
while (defined (my $line = <IN>)) {
  chomp $line;
  my @tmp = split("\t", $line);
  next if ($tmp[1] ne 'BS-seq'); # we don't want to treat the RNA-seq data
  $index{$tmp[8]}++ if ($tmp[3] eq $sample);
  $condition = $tmp[4] if ($tmp[3] eq $sample);
}
close IN;
if (scalar keys %index > 1){
  print '/!\ The index indicated for $sample samples differ!', "\n";
  $warning = 1;
}

######### TrimGalore ##########
my $input_sequenced;
open(IN, '<', '../raw_data/SamplesDirectories+rescued+counts.txt') or die ("Cannot open info file\n");
while (defined (my $line = <IN>)) {
  chomp $line;
  my @tmp = split("\t", $line);
  next if ($tmp[1] ne 'BS-seq'); # we don't want to treat the RNA-seq data
  if (($tmp[0] eq $flowcell) and ($tmp[3] eq $sample)){
    $input_sequenced = $tmp[10];
  }
}
close IN;
print "Number of reads sequenced: ", commify($input_sequenced), "\n";

my $input_trimGalore;
my $output_trimGalore;
my $number_removed;
my $percent_trimmed;
if (-s "$directory/$sample.fastq.gz_trimming_report.txt"){
  open(IN, '<', "$directory/$sample.fastq.gz_trimming_report.txt") or die ("Cannot open trimming report file\n");
  while (defined (my $line = <IN>)) {
    chomp $line;
    if ($line =~ m/(\d+) sequences processed in total/){
      $input_trimGalore = $1;
    }
    if ($line =~ m/Sequences removed because they became shorter than the length cutoff of 15 bp:\s+(\d+)\s\((\d+\.\d+)\%\)/){
      $number_removed = $1;
      $percent_trimmed = $2;
    }
  }
  close IN;
}
else {
  print '/!\ No report file found for TrimGalore. Using .err file', "\n";
  open(IN, '<', "$directory/$sample.err") or die ("Cannot open .err file\n");
  while (defined (my $line = <IN>)) {
    chomp $line;
    if ($line =~ m/(\d+) sequences processed in total/){
      $input_trimGalore = $1;
    }
    if ($line =~ m/Sequences removed because they became shorter than the length cutoff of 15 bp:\s+(\d+)\s\((\d+\.\d+)\%\)/){
      $number_removed = $1;
      $percent_trimmed = $2;
    }
  }
  close IN;
}
$output_trimGalore = $input_trimGalore - $number_removed;
# print "Number of reads analyzed by TrimGalore: ", commify($input_trimGalore), "\n";
if ($input_trimGalore ne $input_sequenced){
  print '/!\ TrimGalore didn\'t process all sequenced reads?', "\n";
  $warning = 1;
}
print "Number of reads after trimming: ", commify($output_trimGalore), " (",$percent_trimmed, "%)\n";

######### bismark ##########
my $input_bismark;
my $output_bismark;
my $percent_mapped;
my $top;
my $bottom;
if (-s "$directory/$sample\_trimmed.fq.gz_bismark_SE_report.txt"){
  open(IN, '<', "$directory/$sample\_trimmed.fq.gz_bismark_SE_report.txt") or die ("Cannot open bismark report\n");
  while (defined (my $line = <IN>)) {
    chomp $line;
    if ($line =~ m/Sequences analysed in total\:\t(\d+)/) {
      $input_bismark = $1;
    }
    if ($line =~ m/Number of alignments with a unique best hit from the different alignments:\t(\d+)/) {
      $output_bismark = $1;
    }
    if ($line =~ m/Mapping efficiency\:\t(\d+\.\d+)\%/) {
      $percent_mapped = $1;
    }
    if ($line =~ m/CT\/CT\:\t(\d+)\t\(\(converted\) top strand\)/) {
      $top = $1;
    }
    if ($line =~ m/CT\/GA\:\t(\d+)\t\(\(converted\) bottom strand\)/) {
      $bottom = $1;
    }
  }
  close IN;
}
else {
  print '/!\ No report file found for bismark mapping', "\n";
  $warning = 1;
}
#print "Number of reads analyzed by bismark: ", commify($input_bismark), "\n";
if ($output_trimGalore ne $input_bismark){
  print '/!\ Bismark didn\'t process all trimmed reads?', "\n";
  $warning = 1;
}
print "Number of reads mapped by bismark: ", commify($output_bismark), " (", $percent_mapped, "%)\n";

# report a problem if 5% difference between 2 numbers
#print "Number of reads mapped on top strand: ", commify($top), "\n";
#print "Number of reads mapped on bottom strand: ", commify($bottom), "\n";
if (($top > 1.05*$bottom) or ($top < 0.95*$bottom)){
  print '/!\ Problem with mapping on top vs. bottom strand with bismark?', " top: $top, bottom: $bottom\n";
  $warning = 1;
}

######### deduplication ##########
my $input_dedup;
my $removed_dedup;
my $output_dedup;
my $percent_dedup;
if (-s "$directory/$sample\_trimmed.fq.gz_bismark.deduplication_report.txt"){
  open(IN, '<', "$directory/$sample\_trimmed.fq.gz_bismark.deduplication_report.txt") or die ("Cannot open dedup report file\n");
  while (defined (my $line = <IN>)) {
    chomp $line;
    if ($line =~ m/Total number of alignments analysed in .+:\t(\d+)/){
      $input_dedup = $1;
    }
    if ($line =~ m/Total number duplicated alignments removed:\t(\d+)\s\((\d+\.\d+)\%\)/){
      $removed_dedup = $1;
      $percent_dedup = $2;
    }
  }
  close IN;
}
else {
  print '/!\ No report file found for deduplication', "\n";
  $warning = 1;
}
$output_dedup = $input_dedup - $removed_dedup;
#print "Number of reads analyzed by deduplication script: ", commify($input_dedup), "\n";
#The number is always slightly lower than numbe rof reported mapped reads. Report a problem if more than 1/1000 difference
if ($input_dedup < 0.999*$output_bismark){
  print '/!\ Deduplication script didn\'t process all mapped reads?', "\n";
  $warning = 1;
}
print "Number of reads after deduplication: ", commify($output_dedup), " (", $percent_dedup, "%)\n";

######### methylation extractor ##########
# No report file here: read in .err and .out file
my $input_meth_extractor;
open(IN, '<', "$directory/$sample.err") or die ("Cannot open .err file\n");
while (defined (my $line = <IN>)) {
  chomp $line;
  if ($line =~ m/Processed (\d+) lines from $directory\/$sample\_trimmed\.fq\.gz\_bismark\.deduplicated\.bam in total/){
    $input_meth_extractor = $1;
  }
}
close IN;
#print "Number of reads analyzed by methylation extractor: ", commify($input_meth_extractor), "\n";
if ($output_dedup ne $input_meth_extractor){
  print '/!\ Methylation extractor didn\'t process all deduplicated reads?', "\n";
  $warning = 1;
}

my $number_conversions;
my $percent_meth_CpG;
my $percent_meth_non_CpG;
open(IN, '<', "$directory/$sample.out") or die ("Cannot open .out file\n");
while (defined (my $line = <IN>)) {
  chomp $line;
  if ($line =~ m/Total C to T conversions in CpG context:\t(\d+)/){
    $number_conversions = $1;
  }
  if ($line =~ m/C methylated in CpG context:\t(\d+\.\d+)\%/){
    $percent_meth_CpG = $1;
  }
  if ($line =~ m/C methylated in non-CpG context:\t(\d+\.\d+)\%/){
    $percent_meth_non_CpG = $1;
  }
}
close IN;
#print "Number of unmethylated CpGs (methylation extractor): ", commify($number_conversions), "\n";
print "Percent methylation in CpG context: $percent_meth_CpG\%\n";
print "Percent methylation in non-CpG context: $percent_meth_non_CpG\%\n";

## launch UNIX command to get the number of unmethylated CpG in each sample
my $number_unmethylated = `zcat $directory/$sample.bedGraph.gz | cut -f6 | paste -s -d + - | bc`;
chomp($number_unmethylated);
#print "Number of unmethylated CpGs (collapsed bedGraph file): ", commify($number_unmethylated), "\n";
if ($number_conversions ne $number_unmethylated){
  print '/!\ Problem in number of CpGs reported in begGraph file?', "\n";
  $warning = 1;
}

######### write report file ##########
open(OUT, '>>', 'checked_samples.txt') or die 'Cannot open OUT file';
# print headers:
if ($first eq 0){
  print OUT "Flow-cell\tSample\tCondition\tNumber of sequenced reads\tNumber of reads after trimming\tPercentage of reads removed after trimming\tNumber of mapped reads\tMapping efficiency\tNumber of reads after deduplication\tPercentage of duplication\tPercentage of methylation in CpG context\tPercentage of methylation in non-CpG context\tWarnings\n";
}

# first check that all variables are defined
if ((!defined $input_sequenced) or (!defined $output_trimGalore) or (!defined $percent_trimmed) or (!defined $output_bismark) or (!defined $percent_mapped) or (!defined $output_dedup) or (!defined $percent_dedup) or (!defined $percent_meth_CpG) or (!defined $percent_meth_non_CpG)){
  print OUT "$flowcell\t$sample\t$condition\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t1\n";
}
else {
  print OUT "$flowcell\t$sample\t$condition\t", commify($input_sequenced), "\t", commify($output_trimGalore), "\t", $percent_trimmed, "\t", commify($output_bismark), "\t", $percent_mapped, "\t", commify($output_dedup), "\t", $percent_dedup, "\t", $percent_meth_CpG, "\t", $percent_meth_non_CpG, "\t", $warning, "\n";
}
close OUT;
exit;

sub commify {
  my $text = reverse $_[0];
  $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
  return scalar reverse $text
}
