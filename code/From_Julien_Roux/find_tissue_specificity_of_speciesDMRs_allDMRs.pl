#!/usr/bin/perl
use strict;
use warnings;

## Julien Roux
## Nov 27, 2014
## find the overlap of species and species-specific DMRs between tissues: do species DMRs tend to affect only subset of tissues? or many tissues?
## Similar to find_tissue_specificity_of_speciesDMRs.pl, but outputs the tissue-specificty status for each DMR

my @contrastsToCheck = ("HumanSpecific", "ChimpSpecific", "HumanChimp", "ChimpRhesus", "HumanRhesus");

foreach my $contrast (@contrastsToCheck){
  print "Checking $contrast DMRs...\n";

  my @tissues = ("heart", "kidney", "liver", "lung");

  # reading DMRs in each tissue
  my %DMRs;
  foreach my $tissue (@tissues){
    print "  Reading DMRs in $tissue...\n";

    my $file = $contrast.'_'.$tissue.'_DMRs.txt';
    #print "  $file\n";
    open(IN, '<', './DMRs/species/'.$file) or die ("Cannot open DMRs file\n");
    my $header = <IN>;
    while (defined (my $line = <IN>)) {
      chomp $line;
      my @tmp = split("\t", $line);
      # tissue -> chr -> start (used as ID) -> start/end
      $DMRs{$tissue}->{$tmp[0]}->{$tmp[1]}->{'start'} = $tmp[1];
      $DMRs{$tissue}->{$tmp[0]}->{$tmp[1]}->{'end'} = $tmp[2];
      $DMRs{$tissue}->{$tmp[0]}->{$tmp[1]}->{'direction'} = $tmp[15];
      $DMRs{$tissue}->{$tmp[0]}->{$tmp[1]}->{'overlap'} = 0;
    }
    close IN;
  }

  # comparing DMRs
  foreach my $tissue1 (@tissues){
  TISSUE2:
    foreach my $tissue2 (@tissues){
      next TISSUE2 if ($tissue1 eq $tissue2);
      print "  Comparing DMRs in $tissue1 vs $tissue2...\n";

      foreach my $chr (sort keys %{$DMRs{$tissue1}}) {
        #print "    $chr\n";
      DMR1:
        foreach my $dmr1 (keys %{$DMRs{$tissue1}->{$chr}}) {
          foreach my $dmr2 (keys %{$DMRs{$tissue2}->{$chr}}) {
            # check coordinates overlap and same direction
            if (($DMRs{$tissue1}->{$chr}->{$dmr1}->{'start'} <= $DMRs{$tissue2}->{$chr}->{$dmr2}->{'end'}) and ($DMRs{$tissue1}->{$chr}->{$dmr1}->{'end'} >= $DMRs{$tissue2}->{$chr}->{$dmr2}->{'start'}) and ($DMRs{$tissue1}->{$chr}->{$dmr1}->{'direction'} eq $DMRs{$tissue2}->{$chr}->{$dmr2}->{'direction'})) {

              # only increment DMR in tissue1 since DMRs in tissue2 will be screened again
              $DMRs{$tissue1}->{$chr}->{$dmr1}->{'overlap'}++;
              next DMR1;
            }
          }
        }
      }
    }
  }

  # output tissue-specificity for each DMR
  foreach my $tissue (@tissues){

    my $outfile = $contrast.'_'.$tissue.'_DMRs_tissueSpecificity.txt';
    open(OUT, '>', './DMRs/species/'.$outfile) or die ("Cannot open OUT file\n");
    print OUT "chr\tstart\tend\tnumberSharedTissues\n"; 

    # Rereading DMRs to print them in same order
    my $file = $contrast.'_'.$tissue.'_DMRs.txt';
    #print "  $file\n";
    open(IN, '<', './DMRs/species/'.$file) or die ("Cannot open DMRs file\n");
    my $header = <IN>;
    while (defined (my $line = <IN>)) {
      chomp $line;
      my @tmp = split("\t", $line);

      print OUT "$tmp[0]\t$tmp[1]\t$tmp[2]\t", $DMRs{$tissue}->{$tmp[0]}->{$tmp[1]}->{'overlap'}, "\n";
    }
    close IN;
    close OUT;
  }
}
exit;


