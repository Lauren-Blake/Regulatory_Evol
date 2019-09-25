#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;

## Julien Roux
## May 30, 2015
## Input 2 DMRs files, and report for each DMR in first file if it overlaps a DMR in the second file (vector of 0 and 1)
## Option: report 1 only if same direction only

scalar @ARGV == 3 or die<<USAGE;
Usage:
perl overlap_DMRs_files.pl <full path to file1> <full path to file2> <direction: same|opposing|NA>
USAGE

my $file1 = $ARGV[0];
my $file2 = $ARGV[1];
my $direction = $ARGV[2];

print "Reading DMR file 2...\n";
my %DMRs;
open(IN, '<', $file2) or die ("Cannot open DMRs file 2\n");
my $header = <IN>;
while (defined (my $line = <IN>)) {
  chomp $line;
  my @tmp = split("\t", $line);
  # chr -> start (used as ID) -> start/end
  $DMRs{$tmp[0]}->{$tmp[1]}->{'start'} = $tmp[1];
  $DMRs{$tmp[0]}->{$tmp[1]}->{'end'} = $tmp[2];
  $DMRs{$tmp[0]}->{$tmp[1]}->{'direction'} = $tmp[15];
}
close IN;

## Open OUT file
my $outfile = $file1;
$outfile =~ s/.txt/_overlap_/;
$outfile .= basename($file2);
open(OUT, '>', $outfile) or die ("Cannot open OUT file\n");

print "Comparing to DMR file 1...\n";
open(IN, '<', $file1) or die ("Cannot open DMRs file 1\n");
$header = <IN>;
while (defined (my $line = <IN>)) {
  chomp $line;
  my @tmp = split("\t", $line);
  my $DMR1_chr = $tmp[0];
  my $DMR1_start = $tmp[1];
  my $DMR1_end = $tmp[2];
  my $DMR1_direction = $tmp[15];

  # check if this DMR overlap one in file 2
  my $overlap = 0;
 DMR2:
  foreach my $dmr2 (keys %{$DMRs{$DMR1_chr}}) {
    # check coordinates overlap
    if (($DMR1_start <= $DMRs{$DMR1_chr}->{$dmr2}->{'end'}) and ($DMR1_end >= $DMRs{$DMR1_chr}->{$dmr2}->{'start'})){
      if (($direction eq 'same') and ($DMR1_direction eq $DMRs{$DMR1_chr}->{$dmr2}->{'direction'})){
        $overlap = 1;
        last DMR2;
      }
      elsif (($direction eq 'opposing') and ($DMR1_direction ne $DMRs{$DMR1_chr}->{$dmr2}->{'direction'})){
        $overlap = 1;
        last DMR2;
      }
      elsif ($direction eq 'NA'){
        $overlap = 1;
        last DMR2;
      }
    }
  }
  print OUT "$DMR1_chr\t$DMR1_start\t$DMR1_end\t$overlap\n";
}
close IN;
close OUT;
exit;
