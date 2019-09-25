#!/usr/bin/perl
use strict;
use warnings;

## Julien Roux
## June 20, 2013
## Generate a map of CpGs from chimp/rhesus to human positions
## identify if these positions are also CpG in human
## This script makes use of cytosine report files in folder ~/Methylation/CpG_maps_to_hg19/cytosine_reports

scalar @ARGV >= 1 or die<<USAGE;
Usage:
perl generate_map.pl <chimp/rhesus>
USAGE
my $species = $ARGV[0];

my $dir;
my $chain;
my $chain_back;
if ($species eq 'chimp'){
  $dir = './cytosine_report/chimp/C4K1_cytosine_report/';
  ## $dir = './cytosine_report/chimp/test/'; ## for testing purposes
  $chain = '/data/share/tridnase/LiftOverChain/panTro3ToHg19.over.chain.gz';
  $chain_back = '/data/share/tridnase/LiftOverChain/hg19ToPanTro3.over.chain.gz';
}
if ($species eq 'rhesus'){
  $dir = './cytosine_report/rhesus/R1K1_cytosine_report/';
  $chain = '/data/share/tridnase/LiftOverChain/rheMac2ToHg19.over.chain.gz';
  $chain_back = '/data/share/tridnase/LiftOverChain/hg19ToRheMac2.over.chain.gz';
}

my $index = 1;
my $outfile = "$species\_all_CpGs.bed.gz";
my %original;
open(OUT, "| gzip -9 - > $outfile") or die "Can't open CpG output file";
opendir(DIR, $dir) or die "can't opendir: $!";
while (defined(my $chr_file = readdir(DIR))) {
  next if ($chr_file !~ m/.+\.chr(chr.+)\.txt\.gz/);
  my $chr = $1;
  print "Reading cytosine reports for chromosome $chr...\n";

  ## read the cytosine file (0-based) and filter to keep only CpGs
  ## Append filtered lines to a compressed file with index as added column
  open(IN, "zcat $dir$chr_file |") or die "Can't read input file for chromosome $chr\n";
  while (defined (my $line = <IN>)) {
    chomp $line;
    my @tmp = split("\t", $line);
    next if ($tmp[5] ne 'CG'); ## only keeping CpGs
    print OUT $tmp[0], "\t", $tmp[1], "\t", $tmp[1]+1, "\t", $tmp[2], "\t", $tmp[6], "\t", $index, "\n";
    ## index -> chr
    ## index -> position
    $original{$index}->{'chr'} = $tmp[0];
    $original{$index}->{'start'} = $tmp[1];
    $original{$index}->{'strand'} = $tmp[2];
    $original{$index}->{'context'} = $tmp[3];
    $index++;
  }
  close IN;
}
closedir(DIR);
close OUT;

## launch liftOver in both directions
print "Lifting over from $species to human...\n";

my $command = "/home/jroux/bin/liftOver -bedPlus=4 -minMatch=0.1 -tab $outfile $chain $species\_ToHg19.bed $species\_ToHg19_unmapped.bed; gzip -9 $species\_ToHg19.bed; gzip -9 $species\_ToHg19_unmapped.bed";
#Launch the command
system($command)==0
    or warn "Failed to launch liftOver\n";

print "Lifting over from human back to $species...\n";
$command = "/home/jroux/bin/liftOver -bedPlus=4 -minMatch=0.1 -tab $species\_ToHg19.bed.gz $chain_back $species\_ToHg19_back.bed $species\_ToHg19_back_unmapped.bed; gzip -9 $species\_ToHg19_back.bed; gzip -9 $species\_ToHg19_back_unmapped.bed";
#Launch the command
system($command)==0
    or warn "Failed to launch liftOver (back to $species genome)\n";


## read both files
## filter out positions thta were not lifted Over back to original genome (or lifted over but not at the same position)
my %back;
open(IN, "zcat $species\_ToHg19_back.bed.gz |") or die "Can't read liftOver (back) results\n";
while (defined (my $line = <IN>)) {
  chomp $line;
  my @tmp = split("\t", $line);
  if (($original{$tmp[5]}->{'chr'} eq $tmp[0]) and ($original{$tmp[5]}->{'start'} eq $tmp[1])){
    $back{$tmp[5]}++;
  }
}


## read human cytosine reports to keep in meory all cytosines contexts
my %human_context;
opendir(DIR, "./cytosine_report/human/") or die "can't opendir: $!";
while (defined(my $chr_file = readdir(DIR))) {
  next if ($chr_file !~ m/.+\.chr(chr.+)\.txt\.gz/);
  my $chr = $1;
  print "Reading cytosine reports for human chromosome $chr...\n";

  open(IN, "zcat ./cytosine_report/human/$chr_file |") or die "Can't read input file for chromosome $chr\n";
  while (defined (my $line = <IN>)) {
    chomp $line;
    my @tmp = split("\t", $line);
    ## chr -> start -> context
    $human_context{$tmp[0]}->{$tmp[1]} = $tmp[6];
  }
  close IN;
  # last; ## TO DO: remove this
}
closedir(DIR);

open(OUT, "| gzip -9 - > $species\_ToHg19.map.gz") or die "Can't open CpG output file";
open(IN, "zcat $species\_ToHg19.bed.gz |") or die "Can't read liftOver results\n";
print OUT "index\torigChr\torigStart\torigEnd\torigStrand\torigContext\thg19Chr\thg19Start\thg19End\thg19Strand\thg19Context\n";
while (defined (my $line = <IN>)) {
  chomp $line;
  my @tmp = split("\t", $line);
  ## if this position was not lifted over back or not to the same position: no 1-to-1 orthology
  next if (!exists $back{$tmp[5]});

  my $human_context;
  if (exists $human_context{$tmp[0]}->{$tmp[1]}){
    $human_context = $human_context{$tmp[0]}->{$tmp[1]};
  }
  else {
    $human_context = 'other';
  }

  ## print map
  print OUT $tmp[5], "\t", $original{$tmp[5]}->{'chr'}, "\t", $original{$tmp[5]}->{'start'}, "\t", $original{$tmp[5]}->{'start'}+1, "\t", $original{$tmp[5]}->{'strand'}, "\t", $original{$tmp[5]}->{'context'}, "\t", $tmp[0], "\t", $tmp[1], "\t", $tmp[2], "\t", $tmp[3], "\t", $human_context, "\n";
}
close IN;
close OUT;

exit;


