#!/usr/bin/perl
use strict;
use warnings;

## Julien Roux
## Aug 25, 2014
## Extract the number of DMRs found by bsseq on permuted sets of samples

## Tissue
opendir(DIR, './DMRs/tissues/permutations/') or die "can't opendir: $!";
my @files = sort readdir(DIR);
my %permutations;
while ( my $file = shift @files ) {
  if ($file =~ m/(\w+)\_(\w+)\_(\w+)\_permut\d+\_DMRs\.txt/){
    #print "$file\n";
    my $species = $1;
    my $tissue1 = $2;
    my $tissue2 = $3;
    $permutations{$species}->{$tissue1}->{$tissue2}->{'num'}++;
    ## each permutation was tested twice (see README.txt), so we'll divide this number by 2 at the end. The min, mean and max number of DMRs is unchanged

    open(IN, '<', './DMRs/tissues/permutations/'.$file) or die ("Cannot open DMRs file\n");
    my $header = <IN>;
    my $count_DMRs = 0;
    while (defined (my $line = <IN>)) {
      #chomp $line;
      #my @tmp = split("\t", $line);
      $count_DMRs++;
    }
    close IN;
    push(@{$permutations{$species}->{$tissue1}->{$tissue2}->{'count'}}, $count_DMRs);
  }
}
closedir(DIR);

open(OUT, '>', 'number_permuted_tissueDMRs.txt');
print OUT "Species\tTissue 1\tTissue 2\tnumber permutations\tmin number DMRs\tmean number DMRs\tmax number DMRs\n";
foreach my $species (keys %permutations){
  foreach my $tissue1 (keys %{$permutations{$species}}){
    foreach my $tissue2 (keys %{$permutations{$species}->{$tissue1}}){
      print OUT "$species\t$tissue1\t$tissue2\t", $permutations{$species}->{$tissue1}->{$tissue2}->{'num'}/2, "\t", min(@{$permutations{$species}->{$tissue1}->{$tissue2}->{'count'}}), "\t", mean(@{$permutations{$species}->{$tissue1}->{$tissue2}->{'count'}}), "\t", max(@{$permutations{$species}->{$tissue1}->{$tissue2}->{'count'}}), "\n";
    }
  }
}
close OUT;

## Species
opendir(DIR, './DMRs/species/permutations/') or die "can't opendir: $!";
@files = sort readdir(DIR);
%permutations = ();
while ( my $file = shift @files ) {
  if ($file =~ m/([H,C,R]\w+)([H,C,R]\w+)\_(\w+)\_permut\d+\_DMRs\.txt/){
    #print "$file\n";
    my $species1 = $1;
    my $species2 = $2;
    my $tissue = $3;
    $permutations{$tissue}->{$species1}->{$species2}->{'num'}++;

    open(IN, '<', './DMRs/species/permutations/'.$file) or die ("Cannot open DMRs file\n");
    my $header = <IN>;
    my $count_DMRs = 0;
    while (defined (my $line = <IN>)) {
      #chomp $line;
      #my @tmp = split("\t", $line);
      $count_DMRs++;
    }
    close IN;
    push(@{    $permutations{$tissue}->{$species1}->{$species2}->{'count'}}, $count_DMRs);
  }
}
closedir(DIR);

open(OUT, '>', 'number_permuted_speciesDMRs.txt');
print OUT "Tissue\tSpecies 1\tSpecies 2\tnumber permutations\tmin number DMRs\tmean number DMRs\tmax number DMRs\n";
foreach my $tissue (keys %permutations){
  foreach my $species1 (keys %{$permutations{$tissue}}){
    foreach my $species2 (keys %{$permutations{$tissue}->{$species1}}){
      print OUT "$tissue\t$species1\t$species2\t", $permutations{$tissue}->{$species1}->{$species2}->{'num'}/2, "\t", min(@{$permutations{$tissue}->{$species1}->{$species2}->{'count'}}), "\t", mean(@{$permutations{$tissue}->{$species1}->{$species2}->{'count'}}), "\t", max(@{$permutations{$tissue}->{$species1}->{$species2}->{'count'}}), "\n";
    }
  }
}
close OUT;


## Individual
opendir(DIR, './DMRs/individuals/permutations/') or die "can't opendir: $!";
@files = sort readdir(DIR);
%permutations = ();
while ( my $file = shift @files ) {
  if ($file =~ m/(\w+)\_(\w+)\_(\w+)\_permut\d+\_DMRs\.txt/){
    #print "$file\n";
    my $species = $1;
    my $individual1 = $2;
    my $individual2 = $3;
    $permutations{$species}->{$individual1}->{$individual2}->{'num'}++;

    open(IN, '<', './DMRs/individuals/permutations/'.$file) or die ("Cannot open DMRs file\n");
    my $header = <IN>;
    my $count_DMRs = 0;
    while (defined (my $line = <IN>)) {
      #chomp $line;
      #my @tmp = split("\t", $line);
      $count_DMRs++;
    }
    close IN;
    push(@{$permutations{$species}->{$individual1}->{$individual2}->{'count'}}, $count_DMRs);
  }
}
closedir(DIR);

open(OUT, '>', 'number_permuted_individualDMRs.txt');
print OUT "Species\tIndividual 1\tIndividual 2\tnumber permutations\tmin number DMRs\tmean number DMRs\tmax number DMRs\n";
foreach my $species (keys %permutations){
  foreach my $individual1 (keys %{$permutations{$species}}){
    foreach my $individual2 (keys %{$permutations{$species}->{$individual1}}){
      print OUT "$species\t$individual1\t$individual2\t", $permutations{$species}->{$individual1}->{$individual2}->{'num'}/2, "\t", min(@{$permutations{$species}->{$individual1}->{$individual2}->{'count'}}), "\t", mean(@{$permutations{$species}->{$individual1}->{$individual2}->{'count'}}), "\t", max(@{$permutations{$species}->{$individual1}->{$individual2}->{'count'}}), "\n";
    }
  }
}
close OUT;
exit;

sub mean {
  my $sum = eval join '+', @_;
  return $sum/scalar(@_);
}
sub min {
  my @array = sort { $a <=> $b } @_;
  return $array[0];
}
sub max {
  my @array = sort { $a <=> $b } @_;
  return $array[-1];
}



