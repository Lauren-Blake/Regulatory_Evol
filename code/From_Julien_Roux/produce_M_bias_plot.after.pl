#!/usr/bin/perl
use strict;
use warnings;

## Julien Roux
## July 26th, 2013
## This is an annex script to launch just the R script if it previously bugged
## This is a subpart of the script launch_meth_extractor.pl

print "Reading information file...\n";
open(IN, '<', "../raw_data/SamplesDirectories.txt") or die ("Cannot open info file\n");
my %samples;
my $header = <IN>;
while (defined (my $line = <IN>)) {
  chomp $line;
  my @tmp = split("\t", $line);
  next if ($tmp[1] ne 'BS-seq'); # we don't want to treat the RNA-seq data
  my $sample = $tmp[3];
  my $flowcell = $tmp[0];
  my $ignore = $tmp[10];

  if ((-s "$flowcell\_$sample/$sample\_trimmed.fq.gz_bismark.deduplicated.M-bias.after.txt") and !(-s "$flowcell\_$sample/$sample\_trimmed.fq.gz_bismark.deduplicated.M-bias.after.pdf")){
    print "Plotting will be launched for sample $sample in flow-cell $flowcell.\n";
    $samples{$sample}->{$flowcell} = $ignore;
  }
}
close IN;

foreach my $sample (sort keys %samples){
  foreach my $flowcell (sort keys %{$samples{$sample}}){
    # launch the job
    my $command = "echo \"R CMD BATCH --no-restore --no-save '--args flowcell=\\\"$flowcell\\\" sample=\\\"$sample\\\"' produce_M_bias_plot.after.R $flowcell\_$sample/$sample\_trimmed.fq.gz_bismark.deduplicated.M-bias.Rout\" | qsub -l h_vmem=512m -v PATH -cwd -N meth_$sample -o $flowcell\_$sample/$sample.out -e $flowcell\_$sample/$sample.err";
    # print "$command\n";

    #Launch the command
    system($command)==0
        or warn "Failed to submit the job to the cluster\n";
  }
}
exit;

