#!perl -w
use strict;
# Julien Roux, Apr 30 2015
# generate sub-files regarding only a region of interest
# This is working species by species. For comparative plots, things are more complex. For now: plot orthologous regions separately and assemble ad hoc

# BEDgraph file of coverage
# BED file of junctions discovered by tophat
# subpart of GTF file of RefSeq models

# Output in folder prepared_regions/

scalar @ARGV == 2 or die<<USAGE;
Usage:
perl prepare_region.pl <comma-separated samples IDs> <region>
USAGE
# e.g., perl prepare_region.pl H1H,H2H,H3H,H4H chr1:500000-1000000

my @samples = split(",", $ARGV[0]);
my $region = $ARGV[1];
$region =~ m/(.+):(\d+)-(\d+)/;
my $chr = $1;
my $start = $2;
my $end = $3;

## make output directory
my $outdir = './prepared_regions/'.$region.'_'.$ARGV[0].'/';
unless(-d $outdir){
  system("mkdir -p $outdir");
}

foreach my $sample (@samples){
  print "Generating files for $sample...\n";
  ## open genome_coverage.bed.gz
  ## generate new output file
  open(IN, 'zcat ./'.$sample.'/genome_coverage.bed.gz |') or die ("Cannot open bedgraph file\n");
  open(OUT, '| gzip - > '.$outdir.'genome_coverage_'.$sample.'.bedgraph.gz');
  my $header = <IN>;
  print OUT $header;
  while (defined (my $line = <IN>)) {
    chomp $line;
    my @tmp = split("\t", $line);
    next if ($tmp[0] ne $chr);
    # we only keep segments overlapping with the region
    next if ($tmp[2] < $start);
    last if ($tmp[1] > $end); ## files are sorted, so once we passed the region we can stop scanning

    # if the line is in the region:
    print OUT $line, "\n";
  }
  close IN;
  close OUT;

  ## open junctions.bed
  ## generate new output file
  open(IN, './'.$sample.'/junctions.bed') or die ("Cannot open junctions bed file\n");
  open(OUT, '| gzip - > '.$outdir.'junctions_'.$sample.'.bed.gz');
  $header = <IN>;
  chomp $header;
  $header =~ s/TopHat junctions/$sample junctions/;
  $header =~ s/name=junctions/name=\"$sample junctions\"/;
  $header .= " graphType=junctions\n";
  print OUT $header;
  while (defined (my $line = <IN>)) {
    chomp $line;
    my @tmp = split("\t", $line);
    next if ($tmp[0] ne $chr);
    # we only keep segments overlapping with the region
    next if ($tmp[2] < $start);
    last if ($tmp[1] > $end); ## files are sorted, so once we passed the region we can stop scanning

    # if the line is in the region:
    print OUT $line, "\n";
  }
  close IN;
  close OUT;
}


## TO DO : - change this and generate annotation depending on the species
##         - plot Refseq annotation for each species
##         - plot orthologous exons coordinates

my $refseq;
my $exons;
if ($samples[0] =~ /^H/){
  $refseq = 'refSeq_hg19.gtf';
  $exons = 'orthoExons_hg19.bed';
}
elsif ($samples[0] =~ /^C/){
  $refseq = 'refSeq_panTro3.gtf';
  $exons = 'orthoExons_panTro3.bed';
}
elsif ($samples[0] =~ /^R/){
  $refseq = 'refSeq_rheMac2.gtf';
  $exons = 'orthoExons_rheMac2.bed';
}

## open GTF
## generate new output file
print "Generating Refseq gene models GTF...\n";
open(IN, '/mnt/lustre/home/jroux/data/GTF/'.$refseq) or die ("Cannot open GTF file\n");
open(OUT, '| gzip - > '.$outdir.$refseq.'.gz');
while (defined (my $line = <IN>)) {
  chomp $line;
  my @tmp = split("\t", $line);
  next if ($tmp[0] ne $chr);
  # we only keep genes overlapping with the region
  next if ($tmp[4] < $start);
  next if ($tmp[3] > $end);

  # if the line is in the region:
  print OUT $line, "\n";
}
close IN;
close OUT;

open(IN, '../orthoExon/'.$exons) or die ("Cannot open bed file\n");
open(OUT, '| gzip - > '.$outdir.$exons.'.gz');
while (defined (my $line = <IN>)) {
  chomp $line;
  my @tmp = split("\t", $line);
  next if ($tmp[0] ne $chr);
  # we only keep segments overlapping with the region
  next if ($tmp[2] < $start);
  next if ($tmp[1] > $end);

  # if the line is in the region:
  print OUT $line, "\n";
}
close IN;
close OUT;
exit;
