#!/usr/bin/perl
use strict;
use warnings;

## Julien Roux
## Aug 20, 2014
## Extract the intersection of parwise species (small-scale) DMRs in each tissue to see if they correspond to species specific DMRs
## DMRs intrsect if they share at least CpG site

my @tissues = ("heart", "kidney", "liver", "lung");
my @species = ("Human", "Chimp", "Rhesus");
my @speciesPairs = ("HumanChimp", "ChimpRhesus", "HumanRhesus");
my %DMRs;
print "Reading DMRs...\n";

foreach my $tissue (@tissues){
  foreach my $pair (@speciesPairs){
    my $file = $pair.'_'.$tissue.'_DMRs.txt';
    print "  $file\n";
    open(IN, '<', './DMRs/species/'.$file) or die ("Cannot open DMRs file\n");
    my $header = <IN>;
    while (defined (my $line = <IN>)) {
      chomp $line;
      my @tmp = split("\t", $line);
      # tissue -> pair -> chr -> start (used as ID) -> start/end
      $DMRs{$tissue}->{$pair}->{$tmp[0]}->{$tmp[1]}->{'start'} = $tmp[1];
      $DMRs{$tissue}->{$pair}->{$tmp[0]}->{$tmp[1]}->{'end'} = $tmp[2];
      $DMRs{$tissue}->{$pair}->{$tmp[0]}->{$tmp[1]}->{'direction'} = $tmp[15];
    }
    close IN;
  }
  foreach my $species (@species[0..1]){
    my $file = $species.'Specific_'.$tissue.'_DMRs.txt';
    print "  $file\n";
    open(IN, '<', './DMRs/species/'.$file) or die ("Cannot open DMRs file\n");
    my $header = <IN>;
    while (defined (my $line = <IN>)) {
      chomp $line;
      my @tmp = split("\t", $line);
      # tissue -> species -> chr -> start (used as ID) -> start/end
      $DMRs{$tissue}->{$species}->{$tmp[0]}->{$tmp[1]}->{'start'} = $tmp[1];
      $DMRs{$tissue}->{$species}->{$tmp[0]}->{$tmp[1]}->{'end'} = $tmp[2];
      $DMRs{$tissue}->{$species}->{$tmp[0]}->{$tmp[1]}->{'direction'} = $tmp[15];
    }
    close IN;
  }
}

print "Human-specific DMRs...\n";
## First intersect Human-chimp and Human-rhesus
foreach my $tissue (@tissues){
  foreach my $chr (sort keys %{$DMRs{$tissue}->{'HumanChimp'}}){
    print "$chr:\n";
    foreach my $dmr1 (keys %{$DMRs{$tissue}->{'HumanChimp'}->{$chr}}){
       # look at coordinates of all DMRs in Human-Rhesus (same chromosome)
      foreach my $dmr2 (keys %{$DMRs{$tissue}->{'HumanRhesus'}->{$chr}}){
        # check overlap and same direction
        if (($DMRs{$tissue}->{'HumanChimp'}->{$chr}->{$dmr1}->{'start'} <= $DMRs{$tissue}->{'HumanRhesus'}->{$chr}->{$dmr2}->{'end'}) and ($DMRs{$tissue}->{'HumanChimp'}->{$chr}->{$dmr1}->{'end'} >= $DMRs{$tissue}->{'HumanRhesus'}->{$chr}->{$dmr2}->{'start'}) and ($DMRs{$tissue}->{'HumanChimp'}->{$chr}->{$dmr1}->{'direction'} eq $DMRs{$tissue}->{'HumanRhesus'}->{$chr}->{$dmr2}->{'direction'})){
          # if there is an overlap, add it to HumanOverlap DMRs in %DMRs
          my $new_start = (sort { $a <=> $b } ($DMRs{$tissue}->{'HumanChimp'}->{$chr}->{$dmr1}->{'start'}, $DMRs{$tissue}->{'HumanRhesus'}->{$chr}->{$dmr2}->{'start'}))[-1]; # max start
          my $new_end = (sort { $a <=> $b } ($DMRs{$tissue}->{'HumanChimp'}->{$chr}->{$dmr1}->{'end'}, $DMRs{$tissue}->{'HumanRhesus'}->{$chr}->{$dmr2}->{'end'}))[0]; # min end
          $DMRs{$tissue}->{'HumanOverlap'}->{$chr}->{$new_start}->{'start'} = $new_start;
          $DMRs{$tissue}->{'HumanOverlap'}->{$chr}->{$new_start}->{'end'} = $new_end;
          $DMRs{$tissue}->{'HumanOverlap'}->{$chr}->{$new_start}->{'direction'} = $DMRs{$tissue}->{'HumanChimp'}->{$chr}->{$dmr1}->{'direction'};
        }
      }
    }

    # compare the above DMRs to Chimp-Rhesus contrast: they should not be significant there
  DMR1:
    foreach my $dmr1 (keys %{$DMRs{$tissue}->{'HumanOverlap'}->{$chr}}){
       # look at coordinates of all DMRs in Chimp-Rhesus (same chromosome)
      foreach my $dmr2 (keys %{$DMRs{$tissue}->{'ChimpRhesus'}->{$chr}}){
        # check overlap and same direction
        if (($DMRs{$tissue}->{'HumanOverlap'}->{$chr}->{$dmr1}->{'start'} <= $DMRs{$tissue}->{'ChimpRhesus'}->{$chr}->{$dmr2}->{'end'}) and ($DMRs{$tissue}->{'HumanOverlap'}->{$chr}->{$dmr1}->{'end'} >= $DMRs{$tissue}->{'ChimpRhesus'}->{$chr}->{$dmr2}->{'start'})){
          # if there is an overlap, delete this DMR
          delete $DMRs{$tissue}->{'HumanOverlap'}->{$chr}->{$dmr1};
          next DMR1;
        }
      }
      # print $chr, "\t", $DMRs{$tissue}->{'HumanOverlap'}->{$chr}->{$dmr1}->{'start'}, "\t", $DMRs{$tissue}->{'HumanOverlap'}->{$chr}->{$dmr1}->{'end'}, "\n";
    }

    print "  Pairwise overlap: ", $chr, " ", scalar keys %{$DMRs{$tissue}->{'HumanOverlap'}->{$chr}}, " DMRs\n";
    print "  Human specific: ", $chr, " ", scalar keys %{$DMRs{$tissue}->{'Human'}->{$chr}}, " DMRs\n";
    # Sensitivity is higher when human specific DMRs are gotten from BSseq directly

    # Compare DMRs obtained by overlap to human specific DMRs (directly with BSseq):
    my $count = 0;
  DMR2:
    foreach my $dmr1 (keys %{$DMRs{$tissue}->{'Human'}->{$chr}}){
      foreach my $dmr2 (keys %{$DMRs{$tissue}->{'HumanOverlap'}->{$chr}}){
        # check overlap and same direction
        if (($DMRs{$tissue}->{'Human'}->{$chr}->{$dmr1}->{'start'} <= $DMRs{$tissue}->{'HumanOverlap'}->{$chr}->{$dmr2}->{'end'}) and ($DMRs{$tissue}->{'Human'}->{$chr}->{$dmr1}->{'end'} >= $DMRs{$tissue}->{'HumanOverlap'}->{$chr}->{$dmr2}->{'start'})){
          # if there is an overlap, implement
          $count++;
          next DMR2;
        }
      }
    }
    print "  ", $count, " human specific are overlapping DMRs from the overlap method\n";

    ## in the other direction
    $count = 0;
  DMR3:
    foreach my $dmr1 (keys %{$DMRs{$tissue}->{'HumanOverlap'}->{$chr}}){
      foreach my $dmr2 (keys %{$DMRs{$tissue}->{'Human'}->{$chr}}){
        # check overlap and same direction
        if (($DMRs{$tissue}->{'HumanOverlap'}->{$chr}->{$dmr1}->{'start'} <= $DMRs{$tissue}->{'Human'}->{$chr}->{$dmr2}->{'end'}) and ($DMRs{$tissue}->{'HumanOverlap'}->{$chr}->{$dmr1}->{'end'} >= $DMRs{$tissue}->{'Human'}->{$chr}->{$dmr2}->{'start'})){
          # if there is an overlap, implement
          $count++;
          next DMR3;
        }
      }
    }
    print "  ", $count, " DMRs from the overlap method are overlapping human specific DMRs\n";
    # problem: this shoudl give the same number as the loop above... why is it sometimes differing by 1 or 2?

  }
}
# TO DOs:
# - fix issue above
# - export to text file the overlap DMRs
# - do the same for chimp-specific DMRs

exit;


