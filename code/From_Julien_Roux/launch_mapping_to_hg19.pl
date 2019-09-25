#!/usr/bin/perl
use strict;
use warnings;

## Julien Roux
## June 21, 2013 :-)
## Launch mapping of bedgraph CpG positions from rhesus and chimp to human positions
## We start from the bedgraph where + and - strand have been merged and use mapping files in ~/Methylation/CpG_maps_to_hg19/by_chr_qsub/

print "Reading information file...\n";
open(IN, '<', "../raw_data/SamplesDirectories.txt") or die ("Cannot open info file\n");
my %to_map;
my $header = <IN>;
while (defined (my $line = <IN>)) {
  chomp $line;
  my @tmp = split("\t", $line);
  next if ($tmp[1] ne 'BS-seq'); # we don't want to treat the RNA-seq data
  my $sample = $tmp[3];
  my $flowcell = $tmp[0];
  next if ($sample =~ m/^H/); # we don't want to map human

  ## check that the bedgraph.gz file exists and is not empty
  if (-s "$flowcell\_$sample/$sample.bedGraph.gz"){
    if (-s "$flowcell\_$sample/$sample\_mapped.bedGraph.gz"){
     print "Sample $sample in flow-cell seems to be already mapped to human.\n";
    }
    else {
      print "Sample $sample in flow-cell $flowcell seems OK. The mapping to human will be launched.\n";
      $to_map{"$flowcell\_$sample"} = $sample;
    }
  }
  else {
    print "Sample $sample in flow-cell $flowcell is not ready: not bedGraph file found.\n";
  }
}
close IN;

foreach my $folder (sort keys %to_map){
  my $sample = $to_map{$folder};
  my $species;
  if ($sample =~ m/^C/){ $species = 'chimp'; }
  if ($sample =~ m/^R/){ $species = 'rhesus'; }

  ## launch script, sort files and compress them
  my $command = "echo \"perl map_to_hg19.pl $folder/$sample\.bedGraph.gz ../CpG_maps_to_hg19/by_chr_qsub/ $species; sort -k1,1 -k2,2n $folder/$sample\_mapped\.bedGraph | gzip -9 > $folder/$sample\_mapped\.bedGraph.gz; rm $folder/$sample\_mapped\.bedGraph; sort -k1,1 -k2,2n $folder/$sample\_mapped_CpG\.bedGraph | gzip -9 > $folder/$sample\_mapped_CpG\.bedGraph.gz; rm $folder/$sample\_mapped_CpG\.bedGraph;\" | qsub -l h_vmem=8g -v PATH -cwd -N map_$sample -o $folder/$sample.out -e $folder/$sample.err";
  #print "$command\n";

  ## Launch the command
  system($command)==0
     or warn "Failed to submit the job to the cluster\n";
}
exit;

