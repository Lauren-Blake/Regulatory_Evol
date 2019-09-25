#!/usr/bin/perl
use strict;
use warnings;

## Julien Roux
## Mar 19, 2014
## Extract the number of DMRs found by bsseq
## Tissue DMRs, Species DMRs, Individual DMRs etc

## Tissue
open(OUT, '>', 'number_tissueDMRs.txt');
print OUT "Species\tTissue 1\tTissue 2\tLocal correction\thyperDMRs\thypoDMRs\n";
opendir(DIR, './DMRs/tissues') or die "can't opendir: $!";
my @files = sort readdir(DIR);
while ( my $file = shift @files ) {
  if ($file =~ m/(\w+)\_(\w+)\_(\w+)\_DMRs(\.uncorrected)?\.txt/){
    #print "$file\n";
    my $species = $1;
    my $tissue1 = $2;
    my $tissue2 = $3;
    my $correction = 'yes';
    $correction = 'no' if (defined $4);

    open(IN, '<', './DMRs/tissues/'.$file) or die ("Cannot open DMRs file\n");
    my $header = <IN>;
    my $count_hyper = 0;
    my $count_hypo = 0;
    while (defined (my $line = <IN>)) {
      chomp $line;
      my @tmp = split("\t", $line);
      if ($tmp[15] eq 'hyper'){ $count_hyper++ }
      if ($tmp[15] eq 'hypo'){ $count_hypo++ }
    }
    close IN;
    print OUT "$species\t$tissue1\t$tissue2\t$correction\t$count_hyper\t$count_hypo\n";
  }
  if ($file =~ m/(\w+)\_(\w+)Specific\_DMRs\.txt/){
    #print "$file\n";
    my $species = $1;
    my $tissue = $2;

    open(IN, '<', './DMRs/tissues/'.$file) or die ("Cannot open DMRs file\n");
    my $header = <IN>;
    my $count_hyper = 0;
    my $count_hypo = 0;
    while (defined (my $line = <IN>)) {
      chomp $line;
      my @tmp = split("\t", $line);
      if ($tmp[15] eq 'hyper'){ $count_hyper++ }
      if ($tmp[15] eq 'hypo'){ $count_hypo++ }
    }
    close IN;
    print OUT "$species\t$tissue\tNA\tyes\t$count_hyper\t$count_hypo\n";
  }
}
closedir(DIR);
close OUT;

## Species
open(OUT, '>', 'number_speciesDMRs.txt');
print OUT "Tissue\tSpecies 1\tSpecies 2\tLocal correction\thyperDMRs\thypoDMRs\n";
opendir(DIR, './DMRs/species') or die "can't opendir: $!";
@files = sort readdir(DIR);
while ( my $file = shift @files ) {
  if ($file =~ m/([H,C,R]\w+)([H,C,R]\w+)\_(\w+)\_DMRs(\.uncorrected)?\.txt/){
    #print "$file\n";
    my $tissue = $3;
    my $species1 = $1;
    my $species2 = $2;
    my $correction = 'yes';
    $correction = 'no' if (defined $4);

    open(IN, '<', './DMRs/species/'.$file) or die ("Cannot open DMRs file\n");
    my $header = <IN>;
    my $count_hyper = 0;
    my $count_hypo = 0;
    while (defined (my $line = <IN>)) {
      chomp $line;
      my @tmp = split("\t", $line);
      if ($tmp[15] eq 'hyper'){ $count_hyper++ }
      if ($tmp[15] eq 'hypo'){ $count_hypo++ }
    }
    close IN;
    print OUT "$tissue\t$species1\t$species2\t$correction\t$count_hyper\t$count_hypo\n";
  }
  if ($file =~ m/(\w+)Specific\_(\w+)\_DMRs\.txt/){
    #print "$file\n";
    my $species = $1;
    my $tissue = $2;

    open(IN, '<', './DMRs/species/'.$file) or die ("Cannot open DMRs file\n");
    my $header = <IN>;
    my $count_hyper = 0;
    my $count_hypo = 0;
    while (defined (my $line = <IN>)) {
      chomp $line;
      my @tmp = split("\t", $line);
      if ($tmp[15] eq 'hyper'){ $count_hyper++ }
      if ($tmp[15] eq 'hypo'){ $count_hypo++ }
    }
    close IN;
    print OUT "$tissue\t$species\tNA\tyes\t$count_hyper\t$count_hypo\n";
  }
}
closedir(DIR);
close OUT;

## Individual
open(OUT, '>', 'number_individualDMRs.txt');
opendir(DIR, './DMRs/individuals') or die "can't opendir: $!";
print OUT "Species\tIndividual 1\tIndividual 2\tLocal correction\thyperDMRs\thypoDMRs\n";
@files = sort readdir(DIR);
while ( my $file = shift @files ) {
  if ($file =~ m/(\w+)\_(\w+)\_(\w+)\_DMRs(\.uncorrected)?\.txt/){
    #print "$file\n";
    my $species = $1;
    my $individual1 = $2;
    my $individual2 = $3;
    my $correction = 'yes';
    $correction = 'no' if (defined $4);

    open(IN, '<', './DMRs/individuals/'.$file) or die ("Cannot open DMRs file\n");
    my $header = <IN>;
    my $count_hyper = 0;
    my $count_hypo = 0;
    while (defined (my $line = <IN>)) {
      chomp $line;
      my @tmp = split("\t", $line);
      if ($tmp[15] eq 'hyper'){ $count_hyper++ }
      if ($tmp[15] eq 'hypo'){ $count_hypo++ }
    }
    close IN;
    print OUT "$species\t$individual1\t$individual2\t$correction\t$count_hyper\t$count_hypo\n";
  }
}
exit;


