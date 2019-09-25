#!/usr/bin/perl
use strict;
use warnings;

## Julien Roux
## July 19th, 2013
## Launch methylation extractor to produce (only) a M-bias plot
## This is useful to know how many bps should be ignored before actually extracting the methylation info with launch_dedup_meth_extractor.pl script

## if we want ot launch only one sample
my $single_sample;
if (defined $ARGV[0]){
  $single_sample = $ARGV[0];
}

my %mapped;
opendir(DIR, '.') or die "can't opendir: $!";
while (defined(my $sample_dir = readdir(DIR))) {
  ## if we treat only 1 sample
  next if ((defined $single_sample) and ($single_sample ne $sample_dir));

  if ($sample_dir =~ m/(Methylation\_\d+)\_(.+)/){
    my $flowcell = $1;
    my $sample = $2;

    unless (-s "$flowcell\_$sample/$sample.report"){
      print "No report file for $flowcell\_$sample/$sample.report ???\n";
      next;
    } 
    ## check that the sample was fully mapped (report produced get_stats.pl)
    open(IN, '<', "$flowcell\_$sample/$sample.report") or die ("Cannot open report file\n");
    my $flag = 0;
    while (defined (my $line = <IN>)) {
      chomp $line;
      if ($line =~ m/Basic statistics from bismark mapping:/){
        $flag = 1;
      }
      if ($line =~ m/Number of reads mapped to top strand: 0/){
        warn("$flowcell\_$sample: there seem to be no reads mapped to top strand. Deduplication and methylation extractor not launched on this sample.\n");
        $flag = 0;
      }
      if ($line =~ m/Number of reads mapped to bottom strand: 0/){
        warn("$flowcell\_$sample: there seem to be no reads mapped to bottom strand. Deduplication and methylation extractor not launched on this sample.\n");
        $flag = 0;
      }
    }
    if ($flag eq 1){
      $mapped{$sample}->{$flowcell} = ();
    }
  }
}
closedir(DIR);

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

  if (!exists $mapped{$sample}->{$flowcell}){
    print "Sample $sample in flow-cell $flowcell is not mapped, not done mapping, or not correctly mapped: use launch_bismark.pl first!\n";
  }
  elsif (-s "$flowcell\_$sample/$sample\_trimmed.fq.gz_bismark.deduplicated.bam_splitting_report.txt"){
    print "Sample $sample in flow-cell $flowcell was already passed through methylation extractor. It will not be relaunched here.\n";
  }
  elsif (-s "$flowcell\_$sample/$sample\_trimmed.fq.gz_bismark.deduplicated.M-bias.before.txt"){
    print "Sample $sample in flow-cell $flowcell was already passed through dedup and produce M-bias plot. It will not be relaunched here.\n";
  }
  else {
    print "Sample $sample in flow-cell $flowcell seems to be mapped correctly. Methylation extractor will be launched.\n";
    ## sample -> flow cell = nb of bps to ignore
    $samples{$sample}->{$flowcell} = $ignore;
  }
}
close IN;

foreach my $sample (sort keys %samples){
  foreach my $flowcell (sort keys %{$samples{$sample}}){
    # launch the job
    my $command = "echo \"perl /home/jroux/bin/bismark_v0.8.1/deduplicate_bismark_alignment_output.pl -s --bam $flowcell\_$sample/$sample\_trimmed.fq.gz_bismark.bam; /home/jroux/bin/bismark_v0.8.1/bismark_methylation_extractor_0.8.1a_CpG_only --single-end -o $flowcell\_$sample/ --report --mbias_only $flowcell\_$sample/$sample\_trimmed.fq.gz_bismark.deduplicated.bam; mv $flowcell\_$sample/$sample\_trimmed.fq.gz_bismark.deduplicated.M-bias.txt $flowcell\_$sample/$sample\_trimmed.fq.gz_bismark.deduplicated.M-bias.before.txt; rm $flowcell\_$sample/$sample\_trimmed.fq.gz_bismark.deduplicated.bam_splitting_report.txt; R CMD BATCH --no-restore --no-save '--args flowcell=\\\"$flowcell\\\" sample=\\\"$sample\\\"' produce_M_bias_plot.before.R $flowcell\_$sample/$sample\_trimmed.fq.gz_bismark.deduplicated.M-bias.Rout\" | qsub -l h_vmem=32g -v PATH -cwd -q batch64.q -N meth_$sample -o $flowcell\_$sample/$sample.out -e $flowcell\_$sample/$sample.err";
    # print "$command\n";

    # add to report file
    open(OUT, '>>', "$flowcell\_$sample/$sample.report") or die 'Cannot open OUT file';
    print OUT "\nCommand line submitted to the cluster (deduplication and methylation extractor):\n$command\n\n";
    close OUT;

    #Launch the command
    system($command)==0
        or warn "Failed to submit the job to the cluster\n";
  }
}
exit;

