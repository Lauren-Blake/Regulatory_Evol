#!/usr/bin/perl
use strict;
use warnings;
$| = 1;

## Julien Roux
## Apr 12, 2016
# Reads all CpGs in human
# Reads all CpGs from chimp and macaque, mapped to hg19
# Makes a file with orthology status for all CpGs (1 if found in species, 0 otherwise). Beware, a 0 in chimp or macaque does mean that this CpG was not used in out dataset in these species. It doesn't mean that it is absent from these species. It is possible that the genome is incomplete, or that we could not assess the orthology with high confidence.

# first record the position of all CpGs on + strand in human genome. ~28M CpGs in total
print "Recording human CpGs...\n";
my %CpG_map;
open(IN, "zcat ../annotation/CpG_map/CpG_map.txt.gz |") or die "Can't read input file\n";
my $count = 0;
while (defined (my $line = <IN>)) {
  $count++;
  ## print the script progress
  if ($count % 1000000 eq 0) {
    print "  ".$count." CpGs parsed\n";
  }
  chomp $line;
  my @tmp = split("\t", $line);
  # Remove lambda phage
  if ($tmp[0] ne 'gi|215104|gb|J02459.1|e'){
    # chr -> position
    $CpG_map{$tmp[0]}->{$tmp[1]}->{'human'} = ();
  }
}
close IN;

## record all orthologous positions in chimp and macaque
## These only include cytosines in CpG context in chimp and macaque. The context can be of course different at human orthologous position
foreach my $species ('chimp', 'rhesus'){
  print "Recording ", $species, " CpGs...\n";
  my @files = glob( "by_chr_qsub/".$species."_ToHg19_chr*.map.gz" );
  foreach my $map (@files){
    print "  Parsing $map\n";
    open(IN, "zcat $map |") or die "Can't read input file\n";
    my $header = <IN>;
    while (defined (my $line = <IN>)) {
      chomp $line;
      my @tmp = split("\t", $line);
      ## All CpGs not among human positions is either on - strand and in CpG context, or in other contexts (other, CAA, ). We record hg19Context to check. But see README.txt file: the hg19Strand is erroneous in this file, it is in fact the strand of the original position (so meaningless here)
      ## Here we record only 1 of the 2 strands. Always + strand is taken, or - if there is no orthologous mapping for + strand
      if ($tmp[4] eq '+'){
        $CpG_map{$tmp[6]}->{$tmp[7]}->{$species} = $tmp[10];
      }
      else {
        if (!exists $CpG_map{$tmp[6]}->{$tmp[7]-1}->{$species}){
          $CpG_map{$tmp[6]}->{$tmp[7]}->{$species} = $tmp[10];
        }
      }
    }
    close IN;
  }
}

print "Going through orthology map to reconcile +/- CpG positions...\n";
# Go through CpG_map hash to reconcile positions on - strand. There are some chimp or macaque + positions, that are mapped to human - strand in CpG context (hg19). These should be moved to + strand to see the orthology with human + position. For positions mapped to - strand but not on a hg19 CpG, there is no need to do this, since there is no orthology.
foreach my $chr (sort keys %CpG_map){
  print "  $chr\n";
  # for each CpG position
  foreach my $start (sort {$a <=> $b} keys %{$CpG_map{$chr}}) {
    foreach my $species ('chimp', 'rhesus'){
      if ((exists $CpG_map{$chr}->{$start}->{$species}) and ($CpG_map{$chr}->{$start}->{$species} =~ m/^CG/) and (exists $CpG_map{$chr}->{$start-1}->{'human'}) and (!exists $CpG_map{$chr}->{$start}->{'human'})){
        $CpG_map{$chr}->{$start-1}->{$species} = '-'; # We record that it was a position mapepd to - strand
        delete $CpG_map{$chr}->{$start}->{$species};
      }
    }
  }
}

print "Printing orthology map...\n";
open(OUT, "| gzip - > orthology_map_allCpGs.txt.gz" ) or die ("Cannot open OUT file\n");
print OUT "chromosome\tposition\thuman_CpG\tchimp_orthologous_CpG\trhesus_orthologous_CpG\n";
foreach my $chr (sort keys %CpG_map){
  print "  $chr\n";
  # for each CpG position
  foreach my $start (sort {$a <=> $b} keys %{$CpG_map{$chr}}) {
    # remove "empty" positions (probably created when the existence of the hash key is queried)
    next if ((!exists $CpG_map{$chr}->{$start}->{'human'}) and (!exists $CpG_map{$chr}->{$start}->{'chimp'}) and (!exists $CpG_map{$chr}->{$start}->{'rhesus'}));
    # else:
    print OUT "$chr\t$start";
    if (exists $CpG_map{$chr}->{$start}->{'human'}){
      print OUT "\t1";
    }
    else {
      print OUT "\t0";
      # Check: CpGs not seen in human should all be in non-CpG context, or not Cs
      if ((exists $CpG_map{$chr}->{$start}->{'chimp'}) and ( $CpG_map{$chr}->{$start}->{'chimp'} =~ m/^CG/)){
        print "Chimp not human: $chr $start / hg19 context: ", $CpG_map{$chr}->{$start}->{'chimp'}, ". Please check!\n";
      }
      if ((exists $CpG_map{$chr}->{$start}->{'rhesus'}) and ( $CpG_map{$chr}->{$start}->{'rhesus'} =~ m/^CG/)){
        print "Macaque not human: $chr $start / hg19 context: ", $CpG_map{$chr}->{$start}->{'rhesus'}, ". Please check!\n";
      }
    }
    foreach my $species ('chimp', 'rhesus'){
      if (exists $CpG_map{$chr}->{$start}->{$species}){
        print OUT "\t1";
        # show positions on - strand
        if ($CpG_map{$chr}->{$start}->{$species} eq '-'){
          print OUT "\t-";
        }
        else {
          print OUT "\t";
        }
      }
      else {
        print OUT "\t0\t";
      }
    }
    print OUT "\n";
  }
}
close OUT;

