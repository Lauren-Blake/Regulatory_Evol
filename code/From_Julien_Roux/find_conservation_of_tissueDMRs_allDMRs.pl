#!/usr/bin/perl
use strict;
use warnings;

## Julien Roux
## Nov 27, 2014
## find the conservation of tissue and tissue-specific DMRs between species
## Modified from script find_conservation_of_tissueDMRs.pl to output the conservaton status of each DMR

my @contrastsToCheck = ("heartSpecific", "kidneySpecific", "liverSpecific", "lungSpecific", "heart_kidney", "heart_liver", "heart_lung", "liver_kidney", "liver_lung", "lung_kidney");

foreach my $contrast (@contrastsToCheck){
  print "Checking $contrast DMRs...\n";

  my @species = ("Human", "Chimp", "Rhesus");

  # reading DMRs in each species
  my %DMRs;
  foreach my $species (@species){
    print "  Reading DMRs in $species...\n";

    my $file = $species.'_'.$contrast.'_DMRs.txt';
    #print "  $file\n";
    open(IN, '<', './DMRs/tissues/'.$file) or die ("Cannot open DMRs file\n");
    my $header = <IN>;
    while (defined (my $line = <IN>)) {
      chomp $line;
      my @tmp = split("\t", $line);
      # species -> chr -> start (used as ID) -> start/end
      $DMRs{$species}->{$tmp[0]}->{$tmp[1]}->{'start'} = $tmp[1];
      $DMRs{$species}->{$tmp[0]}->{$tmp[1]}->{'end'} = $tmp[2];
      $DMRs{$species}->{$tmp[0]}->{$tmp[1]}->{'direction'} = $tmp[15];
      $DMRs{$species}->{$tmp[0]}->{$tmp[1]}->{'overlap'} = 0;
    }
    close IN;
  }

  # comparing DMRs:
  # we compare all species, but only some comparisons make sense
  # How we code them:
  # for human: unique (0), human-chimp (1), human-rhesus (2) and human-chimp-rhesus (4)
  # for chimp: unique (0), human-chimp (1), chimp-rhesus (3) and human-chimp-rhesus (4)
  # for rhesus: unique (0), human-chimp-rhesus (4)

  # All possible transitions:
  # 0 -> 1 (HC) -> 4 (HR)
  #             -> 4 (CR)
  # 0 -> 2 (HR) -> 4 (HC)
  #             -> 4 (CR)
  # 0 -> 3 (CR) -> 4 (HC)
  #             -> 4 (HR)

  foreach my $species1 (@species){
  SPECIES2:
    foreach my $species2 (@species){
      next SPECIES2 if ($species1 eq $species2);
      print "  Comparing DMRs in $species1 vs $species2...\n";

      foreach my $chr (sort keys %{$DMRs{$species1}}) {
        #print "    $chr\n";
      DMR1:
        foreach my $dmr1 (keys %{$DMRs{$species1}->{$chr}}) {

          # no need to screen this DMR if it is shared between all species
          next if ($DMRs{$species1}->{$chr}->{$dmr1}->{'overlap'} eq 4);

          foreach my $dmr2 (keys %{$DMRs{$species2}->{$chr}}) {
            # check coordinates overlap and same direction
            if (($DMRs{$species1}->{$chr}->{$dmr1}->{'start'} <= $DMRs{$species2}->{$chr}->{$dmr2}->{'end'}) and ($DMRs{$species1}->{$chr}->{$dmr1}->{'end'} >= $DMRs{$species2}->{$chr}->{$dmr2}->{'start'}) and ($DMRs{$species1}->{$chr}->{$dmr1}->{'direction'} eq $DMRs{$species2}->{$chr}->{$dmr2}->{'direction'})) {

              # if no overlap seen/tested before
              if (((($species1 eq 'Human') and ($species2 eq 'Chimp')) or (($species1 eq 'Chimp') and ($species2 eq 'Human'))) and ($DMRs{$species1}->{$chr}->{$dmr1}->{'overlap'} eq 0)){
                $DMRs{$species1}->{$chr}->{$dmr1}->{'overlap'} = 1;
              }
              if (((($species1 eq 'Human') and ($species2 eq 'Rhesus')) or (($species1 eq 'Rhesus') and ($species2 eq 'Human'))) and ($DMRs{$species1}->{$chr}->{$dmr1}->{'overlap'} eq 0)){
                $DMRs{$species1}->{$chr}->{$dmr1}->{'overlap'} = 2;
              }
              if (((($species1 eq 'Rhesus') and ($species2 eq 'Chimp')) or (($species1 eq 'Chimp') and ($species2 eq 'Rhesus'))) and ($DMRs{$species1}->{$chr}->{$dmr1}->{'overlap'} eq 0)){
                $DMRs{$species1}->{$chr}->{$dmr1}->{'overlap'} = 3;
              }

              # if some overlap was already seen for this DMR
              # human-chimp was seen before
              if (((($species1 eq 'Rhesus') and ($species2 eq 'Chimp')) or (($species1 eq 'Chimp') and ($species2 eq 'Rhesus')) or (($species1 eq 'Rhesus') and ($species2 eq 'Human')) or (($species1 eq 'Human') and ($species2 eq 'Rhesus'))) and ($DMRs{$species1}->{$chr}->{$dmr1}->{'overlap'} eq 1)){
                $DMRs{$species1}->{$chr}->{$dmr1}->{'overlap'} = 4;
              }
              # human-rhesus was seen before
              if (((($species1 eq 'Rhesus') and ($species2 eq 'Chimp')) or (($species1 eq 'Chimp') and ($species2 eq 'Rhesus')) or (($species1 eq 'Chimp') and ($species2 eq 'Human')) or (($species1 eq 'Human') and ($species2 eq 'Chimp'))) and ($DMRs{$species1}->{$chr}->{$dmr1}->{'overlap'} eq 2)){
                $DMRs{$species1}->{$chr}->{$dmr1}->{'overlap'} = 4;
              }
              # chimp-rhesus was seen before
              if (((($species1 eq 'Human') and ($species2 eq 'Chimp')) or (($species1 eq 'Chimp') and ($species2 eq 'Human')) or (($species1 eq 'Rhesus') and ($species2 eq 'Human')) or (($species1 eq 'Human') and ($species2 eq 'Rhesus'))) and ($DMRs{$species1}->{$chr}->{$dmr1}->{'overlap'} eq 3)){
                $DMRs{$species1}->{$chr}->{$dmr1}->{'overlap'} = 4;
              }
              next DMR1;
            }
          }
        }
      }
    }
  }

  # printing out conservation status of all DMRs
  foreach my $species (@species){
    my $outfile = $species.'_'.$contrast.'_DMRs_conservation.txt';
    open(OUT, '>', './DMRs/tissues/'.$outfile) or die ("Cannot open OUT file\n");
    print OUT "chr\tstart\tend\tconservation\n";

    # Rereading DMRs to print them in same order
    my $file = $species.'_'.$contrast.'_DMRs.txt';
    open(IN, '<', './DMRs/tissues/'.$file) or die ("Cannot open DMRs file\n");
    my $header = <IN>;
    while (defined (my $line = <IN>)) {
      chomp $line;
      my @tmp = split("\t", $line);

      print OUT "$tmp[0]\t$tmp[1]\t$tmp[2]\t", $DMRs{$species}->{$tmp[0]}->{$tmp[1]}->{'overlap'}, "\n";
    }
    close IN;
    close OUT;
  }
}
exit;
