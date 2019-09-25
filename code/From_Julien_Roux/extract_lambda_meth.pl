#!/usr/bin/perl
use strict;
use warnings;

## Julien Roux
## Feb 11, 2014
## Read bedgraph files of technical replicates for each sample and extract the percentage of methylation of lambda phage DNA (indicates conversion efficiency) and the mean coverage on lambda phage DNA
## lambda phage chromosome: gi|215104|gb|J02459.1|e

## read annotation file
my %samples;
open(IN, '<', '../raw_data/SamplesDirectories.txt') or die ("Cannot open info file\n");
while (defined (my $line = <IN>)) {
  chomp $line;
  my @tmp = split("\t", $line);
  next if ($tmp[1] ne 'BS-seq');
  # flow-cell -> technical replicate = sample
  $samples{$tmp[0]}->{$tmp[3]} = $tmp[4];
}
close IN;

open(OUT, '>', 'conversion_efficiency.txt') or die ("Cannot open OUT file\n");
foreach my $flow_cell (keys %samples) {
  foreach my $tech_rep (keys %{$samples{$flow_cell}}) {
    print "  Technical replicate $tech_rep on $flow_cell.\n";

    ## open gzipped bedgraph files
    open(BED, "zcat ../bismark/".$flow_cell."_".$tech_rep."/".$tech_rep.".bedGraph.gz |" ) or die ("Cannot open bedGraph file\n");

    my $sum = 0;
    my $cov = 0;
    my $count = 0;
    while (defined (my $line = <BED>)) {
      chomp $line;
      my @tmp = split("\t", $line);
      next if ($tmp[0] ne 'gi|215104|gb|J02459.1|e');

      $count++;
      $sum += $tmp[3];
      $cov += $tmp[4]+$tmp[5];
    }
    close BED;

    print OUT "$flow_cell\t$tech_rep\t", $samples{$flow_cell}->{$tech_rep}, "\t", 100-($sum/$count), "\t", $cov/$count, "\n";
  }
}
close OUT;


exit;
