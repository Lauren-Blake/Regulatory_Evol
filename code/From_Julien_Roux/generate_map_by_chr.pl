#!/usr/bin/perl
use strict;
use warnings;

## Julien Roux
## June 20, 2013
## Generate a map of CpGs from chimp/rhesus to human positions
## identify if these positions are also CpG in human
## This script is called after the genration of BED files and the liftOver steps

scalar @ARGV >= 1 or die<<USAGE;
Usage:
perl generate_map.pl_by_chr.pl <chimp/rhesus> <chromosome>
USAGE
my $species = $ARGV[0];
my $chr = $ARGV[1];

print "Reading original BED file for chromosome $chr...\n";
my $original_file = "by_chr_qsub/$species\_all_CpGs_$chr.bed.gz";
my %original;
open(IN, "zcat $original_file |") or die "Can't open BED file";
while (defined (my $line = <IN>)) {
  chomp $line;
  my @tmp = split("\t", $line);

  ## $tmp[5] is the index
  $original{$tmp[5]}->{'chr'} = $tmp[0];
  $original{$tmp[5]}->{'start'} = $tmp[1];
  $original{$tmp[5]}->{'strand'} = $tmp[3];
  $original{$tmp[5]}->{'context'} = $tmp[4];
}
close IN;

## reading liftOver result files
## filter out positions that were not lifted Over back to original genome (or lifted over but not at the same position)
print "Reading liftOver result file for chromosome $chr (back)...\n";
my %back;
open(IN, "zcat by_chr_qsub/$species\_ToHg19_$chr\_back.bed.gz |") or die "Can't read liftOver (back) results\n";
while (defined (my $line = <IN>)) {
  chomp $line;
  my @tmp = split("\t", $line);
  if (($original{$tmp[5]}->{'chr'} eq $tmp[0]) and ($original{$tmp[5]}->{'start'} eq $tmp[1])){
    $back{$tmp[5]}++;
  }
}

## keep in memory all positions that will be printed
my %lift;
print "Reading liftOver result file for chromosome $chr...\n";
open(IN, "zcat by_chr_qsub/$species\_ToHg19_$chr.bed.gz |") or die "Can't read liftOver results\n";
while (defined (my $line = <IN>)) {
  chomp $line;
  my @tmp = split("\t", $line);
  ## if this position was not lifted over back or not to the same position: no 1-to-1 orthology
  next if (!exists $back{$tmp[5]});

  ## chr -> start = strand / index
  $lift{$tmp[0]}->{$tmp[1]}->{'index'} = $tmp[5];
  $lift{$tmp[0]}->{$tmp[1]}->{'strand'} = $tmp[3];
  ## Note: strand doesn't seem to be correct in liftOver output files. It seems to be taken as a external column and forwarded to the output
}
close IN;

print "Reading human cytosine context file and printing out map for chromosome $chr...\n";
## open output file
## open(OUT, "| gzip -9 - > by_chr_qsub/$species\_ToHg19_$chr.map.gz") or die "Can't open CpG output file"; ## we will sort and zip the file later
open(OUT, '>', "by_chr_qsub/$species\_ToHg19_$chr.map") or die "Can't open CpG output file";
print OUT "index\torigChr\torigStart\torigEnd\torigStrand\torigContext\thg19Chr\thg19Start\thg19End\thg19Strand\thg19Context\n";

## read human cytosine reports to report all cytosines contexts
opendir(DIR, "./cytosine_report/human/") or die "can't opendir: $!";
while (defined(my $chr_file = readdir(DIR))) {
  next if ($chr_file !~ m/.+\.chr(chr.+)\.txt\.gz/);
  my $chr_human = $1;
  next if (!exists $lift{$chr_human}); ## no need to read to chromosome if no position was lifted over on it
  print "  Reading cytosine reports for human chromosome $chr_human...\n";

  open(IN, "zcat ./cytosine_report/human/$chr_file |") or die "Can't read input file for chromosome $chr_human\n";
  while (defined (my $line = <IN>)) {
    chomp $line;
    my @tmp = split("\t", $line);
    next if (!exists $lift{$tmp[0]}->{$tmp[1]}); # no CpG was lifted over at this position

    my $index = $lift{$tmp[0]}->{$tmp[1]}->{'index'};
    my $human_context = $tmp[6];
    ## print map
    print OUT $index, "\t", $original{$index}->{'chr'}, "\t", $original{$index}->{'start'}, "\t", $original{$index}->{'start'}+1, "\t", $original{$index}->{'strand'}, "\t", $original{$index}->{'context'}, "\t", $tmp[0], "\t", $tmp[1], "\t", $tmp[1]+1, "\t", $lift{$tmp[0]}->{$tmp[1]}->{'strand'}, "\t", $human_context, "\n";
    delete $lift{$tmp[0]}->{$tmp[1]};
  }
  close IN;
}
closedir(DIR);

print "Reading remaining positions to map (human other context)...\n";
## print all remaining positions (other context)
foreach my $left_chr (keys %lift){
  foreach my $left_start (keys %{$lift{$left_chr}}){
    my $index = $lift{$left_chr}->{$left_start}->{'index'};
    print OUT $index, "\t", $original{$index}->{'chr'}, "\t", $original{$index}->{'start'}, "\t", $original{$index}->{'start'}+1, "\t", $original{$index}->{'strand'}, "\t", $original{$index}->{'context'}, "\t", $left_chr, "\t", $left_start, "\t", $left_start, "\t", $lift{$left_chr}->{$left_start}->{'strand'}, "\tother\n";
  }
}
## Note: see comment about strand above. The hg19Strand is wrong here...
## TO DO: end coordinate should be $left_start + 1
close OUT;
exit;
