#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw( min max );
$| = 1;

## Julien Roux
## Apr 20, 2016
# Use orthology_map_allCpGs.txt.gz to calculate the CpG density in 1kb windows, and their orthlogy / evolution

print "Recording human/chimp/rhesus CpGs and their orthology/conservation...\n";
my %orthology;
open(IN, "zcat orthology_map_allCpGs.txt.gz |") or die "Can't read input file\n";
## testing:
## open(IN, "zcat orthology_map_allCpGs.txt.gz | grep chr20 | ") or die "Can't read input file\n";
my $count = 0;
my $header = <IN>;
while (defined (my $line = <IN>)) {
  $count++;
  ## print the script progress
  if ($count % 1000000 eq 0) {
    print "  ".$count." CpGs parsed\n";
  }
  chomp $line;
  my @tmp = split("\t", $line);

  # chr -> position
  $orthology{$tmp[0]}->{$tmp[1]}->{'human'} = $tmp[2];
  $orthology{$tmp[0]}->{$tmp[1]}->{'chimp'} = $tmp[3];
  $orthology{$tmp[0]}->{$tmp[1]}->{'rhesus'} = $tmp[5];
}
close IN;

print "Going through 1kb tiles to calculate CpG density...\n";
open(OUT, "| gzip - > density_allCpGs_1kb.txt.gz" ) or die ("Cannot open OUT file\n");
print OUT "chromosome\tstart\tend\tnumHuman\tnumChimp\tnumRhesus\tnumHumanRefRhesus\tnumChimpRefRhesus\tnum0Human1Chimp\tnum0Human1Rhesus\n";
foreach my $chr (sort keys %orthology){
  # For each chromosome, from min to max position, split in 1kb windows
  my $chrStart = min keys %{$orthology{$chr}};
  my $chrEnd   = max keys %{$orthology{$chr}};
  for (my $i = $chrStart; $i <= $chrEnd; $i += 1000) {
    print OUT "$chr\t$i\t", $i+1000;
    # for each CpG position within window look if a hash position exists
    my $numHuman           = 0;
    my $numChimp           = 0;
    my $numRhesus          = 0;
    my $numHumanRefRhesus  = 0;
    my $numChimpRefRhesus  = 0;
    my $num0Human1Chimp    = 0;
    my $num0Human1Rhesus   = 0;

    for my $position ($i .. $i+999){
      if (exists $orthology{$chr}->{$position}) {
        # record number CpGs in human, chimp and rhesus
        $numHuman  += $orthology{$chr}->{$position}->{'human'};
        $numChimp  += $orthology{$chr}->{$position}->{'chimp'};
        $numRhesus += $orthology{$chr}->{$position}->{'rhesus'};

        if ($orthology{$chr}->{$position}->{'rhesus'} eq 1){
          $numHumanRefRhesus += $orthology{$chr}->{$position}->{'human'};
          $numChimpRefRhesus += $orthology{$chr}->{$position}->{'chimp'};
          if ($orthology{$chr}->{$position}->{'human'} eq 0){
            $num0Human1Rhesus++;
          }
        }
        if (($orthology{$chr}->{$position}->{'chimp'} eq 1) and ($orthology{$chr}->{$position}->{'human'} eq 0)){
          $num0Human1Chimp++;
        }
      }
    }
    # We display:
    # - Number of CpG in human, chimp and rhesus
    #   Note: number of CpGs in chimp and rhesus is not really meaning absence, since there could be no orthology, missing genome part, etc
    # - Number of CpGs in human and chimp, only considering the CpGs present in rhesus
    #   Note: again number of CpGs in chimp coul dbe underestimated
    # - Number of CpGs not seen in human but seen in chimp or macaque: we know they were gained
    print OUT "\t", $numHuman, "\t", $numChimp, "\t", $numRhesus, "\t", $numHumanRefRhesus, "\t", $numChimpRefRhesus, "\t", $num0Human1Chimp, "\t", $num0Human1Rhesus, "\n";
  }
}
close OUT;
## Then use R script to calculate and use:
# - densities
# - absolute number of CpGs lost compared to rhesus
# - number of CpGs gained compared to human
# - only use windows with a least 1 human, 1 chimp and 1 rhesus CpG detected?


