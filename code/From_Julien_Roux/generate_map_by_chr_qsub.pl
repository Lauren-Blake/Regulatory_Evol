#!/usr/bin/perl
use strict;
use warnings;

## Julien Roux
## June 20, 2013
## Generate a map of CpGs from chimp/rhesus to human positions
## identify if these positions are also CpG in human
## This script launches perl script generate_bed_by_chr.pl, then lifOver, liftOver back, and then perl script generate_map_by_chr.pl

scalar @ARGV >= 1 or die<<USAGE;
Usage:
perl generate_map_by_chr.pl <chimp/rhesus>
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
my %original;
opendir(DIR, $dir) or die "can't opendir: $!";
while (defined(my $chr_file = readdir(DIR))) {
  next if ($chr_file !~ m/.+\.chr(chr.+)\.txt\.gz/);
  my $chr = $1;
  my $bedfile = "$species\_all_CpGs_$chr.bed.gz";

  ## launch liftOver in both directions
  print "Generating BED file, then lifting over from $species to human and from human back to $species ($chr)...\n";
  ## full command
  my $command = "echo \"perl generate_bed_by_chr.pl $species $chr; /home/jroux/bin/liftOver -bedPlus=4 -minMatch=0.1 -tab by_chr_qsub/$bedfile $chain by_chr_qsub/$species\_ToHg19_$chr.bed by_chr_qsub/$species\_ToHg19_$chr\_unmapped.bed; gzip -9 by_chr_qsub/$species\_ToHg19_$chr.bed; gzip -9 by_chr_qsub/$species\_ToHg19_$chr\_unmapped.bed; /home/jroux/bin/liftOver -bedPlus=4 -minMatch=0.1 -tab by_chr_qsub/$species\_ToHg19_$chr.bed.gz $chain_back by_chr_qsub/$species\_ToHg19_$chr\_back.bed by_chr_qsub/$species\_ToHg19_$chr\_back_unmapped.bed; gzip -9 by_chr_qsub/$species\_ToHg19_$chr\_back.bed; gzip -9 by_chr_qsub/$species\_ToHg19_$chr\_back_unmapped.bed; perl generate_map_by_chr.pl $species $chr; sort -k1,1n by_chr_qsub/$species\_ToHg19_$chr.map | gzip -9 > by_chr_qsub/$species\_ToHg19_$chr.map.gz; rm by_chr_qsub/$species\_ToHg19_$chr.map\" | qsub -l h_vmem=16g -v PATH -cwd -N liftOver_$chr -o by_chr_qsub/map_$species\_$chr.out -e by_chr_qsub/map_$species\_$chr.err\n";

  ## only beginning
  ## my $command = "echo \"perl generate_bed_by_chr.pl $species $chr; /home/jroux/bin/liftOver -bedPlus=4 -minMatch=0.1 -tab by_chr_qsub/$bedfile $chain by_chr_qsub/$species\_ToHg19_$chr.bed by_chr_qsub/$species\_ToHg19_$chr\_unmapped.bed; gzip -9 by_chr_qsub/$species\_ToHg19_$chr.bed; gzip -9 by_chr_qsub/$species\_ToHg19_$chr\_unmapped.bed; /home/jroux/bin/liftOver -bedPlus=4 -minMatch=0.1 -tab by_chr_qsub/$species\_ToHg19_$chr.bed.gz $chain_back by_chr_qsub/$species\_ToHg19_$chr\_back.bed by_chr_qsub/$species\_ToHg19_$chr\_back_unmapped.bed; gzip -9 by_chr_qsub/$species\_ToHg19_$chr\_back.bed; gzip -9 by_chr_qsub/$species\_ToHg19_$chr\_back_unmapped.bed\" | qsub -l h_vmem=8g -v PATH -cwd -N liftOver_$chr -o by_chr_qsub/map_$species\_$chr.out -e by_chr_qsub/map_$species\_$chr.err\n";
      
  ## only end
  ##my $command = "echo \"perl generate_map_by_chr.pl $species $chr; sort -k1,1n by_chr_qsub/$species\_ToHg19_$chr.map | gzip -9 > by_chr_qsub/$species\_ToHg19_$chr.map.gz; rm by_chr_qsub/$species\_ToHg19_$chr.map\" | qsub -l h_vmem=16g -v PATH -cwd -N liftOver_$chr -o by_chr_qsub/map_$species\_$chr.out -e by_chr_qsub/map_$species\_$chr.err\n";

  #Launch the command
  system($command)==0
      or warn "Failed to launch command\n";
}
closedir(DIR);
exit;


