#!/usr/bin/perl
use strict;
use warnings;

## Julien Roux
## May 31, 2013
## Deduplication of mapped reads should already have been done (check)
## Launch methylation extractor and compression of output files
## Note (June 17th, 2013): we use a modified bismark_methylation_extractor_CpG_only script that doesn't write huge Non_CpG files (that we remove anyway)... These were problematic because disk space is limited on the cluster

my %mapped;
opendir(DIR, '.') or die "can't opendir: $!";
while (defined(my $sample_dir = readdir(DIR))) {
  if ($sample_dir =~ m/(Methylation\_\d+)\_(.+)/){
    my $flowcell = $1;
    my $sample = $2;

    ## check that the sample was deduplicated 
    my $flag = 0;
    if (-s "$flowcell\_$sample/$sample\_trimmed.fq.gz_bismark.deduplication_report.txt"){
      open(IN, '<', "$flowcell\_$sample/$sample\_trimmed.fq.gz_bismark.deduplication_report.txt") or die ("Cannot open dedup report file\n");
      while (defined (my $line = <IN>)) {
        chomp $line;
        if ($line =~ m/Total number of alignments analysed in/){
          $flag = 1;
        }
      }
    }
    if ($flag eq 1){
      $mapped{$sample}->{$flowcell} = ();
    }
    else {
      print "Sample $sample in flow-cell $flowcell was apparently not deduplicated. Run dedup_and_produce_M_bias_plot.pl before launch_meth_extractor.pl\n";
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
    print "Sample $sample in flow-cell $flowcell is not deduplicated: use launch_bismark.pl / dedup_and_produce_M_bias_plot.pl first!\n";
  }
  elsif((!defined $ignore) or ($ignore eq '')){
    print "No indication of 5'end bps to ignore for sample $sample in flow-cell $flowcell. Was there a problem with dedup_and_produce_M_bias_plot.pl?\n";
  }
  #elsif (-s "$flowcell\_$sample/$sample\_trimmed.fq.gz_bismark.deduplicated.bam_splitting_report.txt"){
  elsif (-e "$flowcell\_$sample/$sample\_trimmed.fq.gz_bismark.deduplicated.bam_splitting_report.txt"){
    print "Sample $sample in flow-cell $flowcell was already passed through methylation extractor. It will not be relaunched here.\n";
  }
  else {
    print "Methylation extractor will be launched on sample $sample in flow-cell $flowcell\n";
    ## sample -> flow cell = nb of bps to ignore
    $samples{$sample}->{$flowcell} = $ignore;
  }
}
close IN;

foreach my $sample (sort keys %samples){
  foreach my $flowcell (sort keys %{$samples{$sample}}){
    # launch the job
    my $genome;
    if ($sample =~ m/^H/) {     ## human
      $genome = "/home/jroux/data/indexes/Bowtie/hg19_norandom_noM_noY+Lambda_prepared_bismark/";
    }
    elsif ($sample =~ m/^C/) {  ## chimp
      $genome = "/home/jroux/data/indexes/Bowtie/panTro3_nonrandom+Lambda_prepared_bismark/";
    }
    elsif ($sample =~ m/^R/) {  ## rhesus macaque
      $genome = "/home/jroux/data/indexes/Bowtie/rheMac2_norandom_noUr+Lambda_prepared_bismark/";
    }

    my $command = "echo \"/home/jroux/bin/bismark_v0.8.1/bismark_methylation_extractor_0.8.1a_CpG_only --single-end --merge_non_CpG -o $flowcell\_$sample/ --report --bedGraph --counts --ignore $samples{$sample}->{$flowcell} --genome_folder $genome $flowcell\_$sample/$sample\_trimmed.fq.gz_bismark.deduplicated.bam; gzip -f -9 $flowcell\_$sample/$sample\_trimmed.fq.gz_bismark.deduplicated.bedGraph; rm $flowcell\_$sample/Non_CpG_CTOB_$sample\_trimmed.fq.gz_bismark.deduplicated.txt; rm $flowcell\_$sample/Non_CpG_CTOT_$sample\_trimmed.fq.gz_bismark.deduplicated.txt; rm $flowcell\_$sample/Non_CpG_OB_$sample\_trimmed.fq.gz_bismark.deduplicated.bam.txt; rm $flowcell\_$sample/Non_CpG_OT_$sample\_trimmed.fq.gz_bismark.deduplicated.txt; rm $flowcell\_$sample/CpG_CTOB_$sample\_trimmed.fq.gz_bismark.deduplicated.txt; rm $flowcell\_$sample/CpG_CTOT_$sample\_trimmed.fq.gz_bismark.deduplicated.txt; gzip -f -9 $flowcell\_$sample/CpG_OB_$sample\_trimmed.fq.gz_bismark.deduplicated.txt; gzip -f -9 $flowcell\_$sample/CpG_OT_$sample\_trimmed.fq.gz_bismark.deduplicated.txt; mv $flowcell\_$sample/$sample\_trimmed.fq.gz_bismark.deduplicated.M-bias.txt $flowcell\_$sample/$sample\_trimmed.fq.gz_bismark.deduplicated.M-bias.after.txt; R CMD BATCH --no-restore --no-save '--args flowcell=\\\"$flowcell\\\" sample=\\\"$sample\\\"' produce_M_bias_plot.after.R $flowcell\_$sample/$sample\_trimmed.fq.gz_bismark.deduplicated.M-bias.Rout; perl collapse_strands.pl $flowcell\_$sample/$sample $flowcell\_$sample/$sample\_trimmed.fq.gz_bismark.deduplicated.bedGraph.gz $flowcell\_$sample/CpG_OT_$sample\_trimmed.fq.gz_bismark.deduplicated.txt.gz\" | qsub -l h_vmem=16g -v PATH -cwd -q batch64.q -N meth_$sample -o $flowcell\_$sample/$sample.out -e $flowcell\_$sample/$sample.err";
    # print "$command\n";

    # add to report file
    open(OUT, '>>', "$flowcell\_$sample/$sample.report") or die 'Cannot open OUT file';
    print OUT "\nCommand line submitted to the cluster (methylation extractor):\n$command\n\n";
    close OUT;

    #Launch the command
    system($command)==0
        or warn "Failed to submit the job to the cluster\n";
  }
}
exit;

