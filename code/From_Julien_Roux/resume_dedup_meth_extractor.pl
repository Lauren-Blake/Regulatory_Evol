#!/usr/bin/perl
use strict;
use warnings;

## Julien Roux
## June 17, 2013
## - script similar to launch_dedup_meth_extractor.pl, but only 1 sample is launched (taken as argument), and only part of the whole command is launched (in case the command was interupted for any reason)
## - If we put KG as argument: we want to write down results in /KG/jroux/Methylation/bismark/ sup disk space
## - we use a modified bismark_methylation_extractor_CpG_only script that doesn't write huge Non_CpG files (that we remove anyway...)

die "before re-using this script, chnage the version number for bismark, and check that the KG harddrive doesn't have any problems (some drives crashed last time)\n";

scalar @ARGV >= 1 or die<<USAGE;
Usage:
perl resume_dedup_meth_extractor.pl <sample> <KG>
USAGE
my $sample = $ARGV[0];
my $kg = 0;
if ($ARGV[1] eq 'KG'){
  $kg = 1;
  print "Output directory is changed to /KG/jroux/Methylation/bismark/$sample/\n";
}

my %mapped;
opendir(DIR, '.') or die "can't opendir: $!";
while (defined(my $sample_dir = readdir(DIR))) {
  if ($sample_dir =~ m/$sample/){
    $sample =~ m/(Methylation\_\d+)\_(.+)/;
    my $flowcell = $1;
    my $sampleID = $2;

    ## check that the sample was fully mapped (report produced get_stats.pl)
    open(IN, '<', "$sample/$sampleID.report") or die ("Cannot open report file\n");
    my $flag = 0;
    while (defined (my $line = <IN>)) {
      chomp $line;
      if ($line =~ m/Basic statistics from bismark mapping:/){
        $flag = 1;
      }
      if ($line =~ m/Number of reads mapped to top strand: 0/){
        warn("$sample: there seem to be no reads mapped to top strand. Deduplication and methylation extractor not launched on this sample.\n");
        $flag = 0;
      }
      if ($line =~ m/Number of reads mapped to bottom strand: 0/){
        warn("$sample: there seem to be no reads mapped to bottom strand. Deduplication and methylation extractor not launched on this sample.\n");
        $flag = 0;
      }
    }
    if ($flag eq 1){
      $mapped{$sampleID}->{$flowcell} = ();
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
  my $sampleID = $tmp[3];
  my $flowcell = $tmp[0];
  my $ignore = $tmp[10];

  if (exists $mapped{$sampleID}->{$flowcell}){
    print "Sample $sampleID in flow-cell $flowcell seems to be mapped correctly. Methylation extractor will be resumed.\n";
    ## sampleID -> flow cell = nb of bps to ignore
    $samples{$sampleID}->{$flowcell} = $ignore;
  }
}
close IN;

foreach my $sampleID (sort keys %samples){
  foreach my $flowcell (sort keys %{$samples{$sampleID}}){
    # launch the job
    my $genome;
    if ($sampleID =~ m/^H/) {     ## human
      $genome = "/home/jroux/data/indexes/Bowtie/hg19_norandom_noM_noY+Lambda_prepared_bismark/";
    }
    elsif ($sampleID =~ m/^C/) {  ## chimp
      $genome = "/home/jroux/data/indexes/Bowtie/panTro3_nonrandom+Lambda_prepared_bismark/";
    }
    elsif ($sampleID =~ m/^R/) {  ## rhesus macaque
      $genome = "/home/jroux/data/indexes/Bowtie/rheMac2_norandom_noUr+Lambda_prepared_bismark/";
    }

    my $command;
    ## if the deduplication looks fine:
    if ((-s "$flowcell\_$sampleID/$sampleID\_trimmed.fq.gz_bismark.deduplication_report.txt") and !(-s "$flowcell\_$sampleID/$sampleID\_trimmed.fq.gz_bismark.deduplicated.bam_splitting_report.txt")){
      open(IN, '<', "$flowcell\_$sampleID/$sampleID\_trimmed.fq.gz_bismark.deduplication_report.txt") or die ("Cannot open dedup report file\n");
      while (defined (my $line = <IN>)) {
        chomp $line;
        if ($line =~ m/Total number of alignments analysed in/){
          print "Sample $sampleID in flow-cell $flowcell seems to be correctly deduplicated, but not methylation extracted.\n";
          if ($kg eq 0){
            $command = "echo \"/home/jroux/bin/bismark_v0.7.12/bismark_methylation_extractor_CpG_only --single-end --merge_non_CpG -o $flowcell\_$sampleID/ --report --bedGraph --counts --ignore $samples{$sampleID}->{$flowcell} --genome_folder $genome $flowcell\_$sampleID/$sampleID\_trimmed.fq.gz_bismark.deduplicated.bam; gzip -9 $flowcell\_$sampleID/$sampleID\_trimmed.fq.gz_bismark.deduplicated.bedGraph; rm $flowcell\_$sampleID/Non_CpG_CTOB_$sampleID\_trimmed.fq.gz_bismark.deduplicated.txt; rm $flowcell\_$sampleID/Non_CpG_CTOT_$sampleID\_trimmed.fq.gz_bismark.deduplicated.txt; rm $flowcell\_$sampleID/Non_CpG_OB_$sampleID\_trimmed.fq.gz_bismark.deduplicated.bam.txt; rm $flowcell\_$sampleID/Non_CpG_OT_$sampleID\_trimmed.fq.gz_bismark.deduplicated.txt; rm $flowcell\_$sampleID/CpG_CTOB_$sampleID\_trimmed.fq.gz_bismark.deduplicated.txt; rm $flowcell\_$sampleID/CpG_CTOT_$sampleID\_trimmed.fq.gz_bismark.deduplicated.txt; gzip -9 $flowcell\_$sampleID/CpG_OB_$sampleID\_trimmed.fq.gz_bismark.deduplicated.txt; gzip -9 $flowcell\_$sampleID/CpG_OT_$sampleID\_trimmed.fq.gz_bismark.deduplicated.txt; perl collapse_strands.pl $flowcell\_$sampleID/$sampleID $flowcell\_$sampleID/$sampleID\_trimmed.fq.gz_bismark.deduplicated.bedGraph.gz $flowcell\_$sampleID/CpG_OT_$sampleID\_trimmed.fq.gz_bismark.deduplicated.txt.gz\" | qsub -l h_vmem=16g -v PATH -cwd -N meth_$sampleID -o $flowcell\_$sampleID/$sampleID.out -e $flowcell\_$sampleID/$sampleID.err";
          }
          else {
            $command = "echo \"mkdir -p /KG/jroux/Methylation/bismark/$sample; /home/jroux/bin/bismark_v0.7.12/bismark_methylation_extractor_CpG_only --single-end --merge_non_CpG -o /KG/jroux/Methylation/bismark/$flowcell\_$sampleID/ --report --bedGraph --counts --ignore $samples{$sampleID}->{$flowcell} --genome_folder $genome $flowcell\_$sampleID/$sampleID\_trimmed.fq.gz_bismark.deduplicated.bam; gzip -9 /KG/jroux/Methylation/bismark/$flowcell\_$sampleID/$sampleID\_trimmed.fq.gz_bismark.deduplicated.bedGraph; rm /KG/jroux/Methylation/bismark/$flowcell\_$sampleID/Non_CpG_CTOB_$sampleID\_trimmed.fq.gz_bismark.deduplicated.txt; rm /KG/jroux/Methylation/bismark/$flowcell\_$sampleID/Non_CpG_CTOT_$sampleID\_trimmed.fq.gz_bismark.deduplicated.txt; rm /KG/jroux/Methylation/bismark/$flowcell\_$sampleID/Non_CpG_OB_$sampleID\_trimmed.fq.gz_bismark.deduplicated.bam.txt; rm /KG/jroux/Methylation/bismark/$flowcell\_$sampleID/Non_CpG_OT_$sampleID\_trimmed.fq.gz_bismark.deduplicated.txt; rm /KG/jroux/Methylation/bismark/$flowcell\_$sampleID/CpG_CTOB_$sampleID\_trimmed.fq.gz_bismark.deduplicated.txt; rm /KG/jroux/Methylation/bismark/$flowcell\_$sampleID/CpG_CTOT_$sampleID\_trimmed.fq.gz_bismark.deduplicated.txt; gzip -9 /KG/jroux/Methylation/bismark/$flowcell\_$sampleID/CpG_OB_$sampleID\_trimmed.fq.gz_bismark.deduplicated.txt; gzip -9 /KG/jroux/Methylation/bismark/$flowcell\_$sampleID/CpG_OT_$sampleID\_trimmed.fq.gz_bismark.deduplicated.txt; perl collapse_strands.pl /KG/jroux/Methylation/bismark/$flowcell\_$sampleID/$sampleID /KG/jroux/Methylation/bismark/$flowcell\_$sampleID/$sampleID\_trimmed.fq.gz_bismark.deduplicated.bedGraph.gz /KG/jroux/Methylation/bismark/$flowcell\_$sampleID/CpG_OT_$sampleID\_trimmed.fq.gz_bismark.deduplicated.txt.gz; mv /KG/jroux/Methylation/bismark/$flowcell\_$sampleID/* $flowcell\_$sampleID/; rmdir /KG/jroux/Methylation/bismark/$flowcell\_$sampleID/\" | qsub -l h_vmem=16g -v PATH -cwd -N meth_$sampleID -o $flowcell\_$sampleID/$sampleID\_KG.out -e $flowcell\_$sampleID/$sampleID\_KG.err";
          }
        }
      }
      close IN;
    }
    elsif ((-s "$flowcell\_$sampleID/$sampleID\_trimmed.fq.gz_bismark.deduplicated.bam_splitting_report.txt") and !(-s "$flowcell\_$sampleID/$sampleID\_trimmed.fq.gz_bismark.deduplicated.bedGraph.gz")){
      print "There was a problem with methylation extractor for sample $sampleID in flow-cell $flowcell. It will be relaunched.\n";
      if ($kg eq 0){
        $command = "echo \"rm $flowcell\_$sampleID/CpG_*; rm $flowcell\_$sampleID/Non_CpG_*; rm $flowcell\_$sampleID/$sampleID\_trimmed.fq.gz_bismark.deduplicated.bam_splitting_report.txt; rm $flowcell\_$sampleID/*.bedGraph.gz; rm $flowcell\_$sampleID/*.bedGraph; /home/jroux/bin/bismark_v0.7.12/bismark_methylation_extractor_CpG_only --single-end --merge_non_CpG -o $flowcell\_$sampleID/ --report --bedGraph --counts --ignore $samples{$sampleID}->{$flowcell} --genome_folder $genome $flowcell\_$sampleID/$sampleID\_trimmed.fq.gz_bismark.deduplicated.bam; gzip -9 $flowcell\_$sampleID/$sampleID\_trimmed.fq.gz_bismark.deduplicated.bedGraph; rm $flowcell\_$sampleID/Non_CpG_CTOB_$sampleID\_trimmed.fq.gz_bismark.deduplicated.txt; rm $flowcell\_$sampleID/Non_CpG_CTOT_$sampleID\_trimmed.fq.gz_bismark.deduplicated.txt; rm $flowcell\_$sampleID/Non_CpG_OB_$sampleID\_trimmed.fq.gz_bismark.deduplicated.bam.txt; rm $flowcell\_$sampleID/Non_CpG_OT_$sampleID\_trimmed.fq.gz_bismark.deduplicated.txt; rm $flowcell\_$sampleID/CpG_CTOB_$sampleID\_trimmed.fq.gz_bismark.deduplicated.txt; rm $flowcell\_$sampleID/CpG_CTOT_$sampleID\_trimmed.fq.gz_bismark.deduplicated.txt; gzip -9 $flowcell\_$sampleID/CpG_OB_$sampleID\_trimmed.fq.gz_bismark.deduplicated.txt; gzip -9 $flowcell\_$sampleID/CpG_OT_$sampleID\_trimmed.fq.gz_bismark.deduplicated.txt; perl collapse_strands.pl $flowcell\_$sampleID/$sampleID $flowcell\_$sampleID/$sampleID\_trimmed.fq.gz_bismark.deduplicated.bedGraph.gz $flowcell\_$sampleID/CpG_OT_$sampleID\_trimmed.fq.gz_bismark.deduplicated.txt.gz\" | qsub -l h_vmem=16g -v PATH -cwd -N meth_$sampleID -o $flowcell\_$sampleID/$sampleID.out -e $flowcell\_$sampleID/$sampleID.err";
      }
      else {
        $command = "echo \"mkdir -p /KG/jroux/Methylation/bismark/$sample; rm $flowcell\_$sampleID/CpG_*; rm $flowcell\_$sampleID/Non_CpG_*; rm $flowcell\_$sampleID/$sampleID\_trimmed.fq.gz_bismark.deduplicated.bam_splitting_report.txt; rm $flowcell\_$sampleID/*.bedGraph.gz; rm $flowcell\_$sampleID/*.bedGraph; rm /KG/jroux/Methylation/bismark/$flowcell\_$sampleID/CpG_*; rm /KG/jroux/Methylation/bismark/$flowcell\_$sampleID/Non_CpG_*; rm /KG/jroux/Methylation/bismark/$flowcell\_$sampleID/$sampleID\_trimmed.fq.gz_bismark.deduplicated.bam_splitting_report.txt; rm /KG/jroux/Methylation/bismark/$flowcell\_$sampleID/*.bedGraph.gz; rm /KG/jroux/Methylation/bismark/$flowcell\_$sampleID/*.bedGraph; /home/jroux/bin/bismark_v0.7.12/bismark_methylation_extractor_CpG_only --single-end --merge_non_CpG -o /KG/jroux/Methylation/bismark/$flowcell\_$sampleID/ --report --bedGraph --counts --ignore $samples{$sampleID}->{$flowcell} --genome_folder $genome $flowcell\_$sampleID/$sampleID\_trimmed.fq.gz_bismark.deduplicated.bam; gzip -9 /KG/jroux/Methylation/bismark/$flowcell\_$sampleID/$sampleID\_trimmed.fq.gz_bismark.deduplicated.bedGraph; rm /KG/jroux/Methylation/bismark/$flowcell\_$sampleID/Non_CpG_CTOB_$sampleID\_trimmed.fq.gz_bismark.deduplicated.txt; rm /KG/jroux/Methylation/bismark/$flowcell\_$sampleID/Non_CpG_CTOT_$sampleID\_trimmed.fq.gz_bismark.deduplicated.txt; rm /KG/jroux/Methylation/bismark/$flowcell\_$sampleID/Non_CpG_OB_$sampleID\_trimmed.fq.gz_bismark.deduplicated.bam.txt; rm /KG/jroux/Methylation/bismark/$flowcell\_$sampleID/Non_CpG_OT_$sampleID\_trimmed.fq.gz_bismark.deduplicated.txt; rm /KG/jroux/Methylation/bismark/$flowcell\_$sampleID/CpG_CTOB_$sampleID\_trimmed.fq.gz_bismark.deduplicated.txt; rm /KG/jroux/Methylation/bismark/$flowcell\_$sampleID/CpG_CTOT_$sampleID\_trimmed.fq.gz_bismark.deduplicated.txt; gzip -9 /KG/jroux/Methylation/bismark/$flowcell\_$sampleID/CpG_OB_$sampleID\_trimmed.fq.gz_bismark.deduplicated.txt; gzip -9 /KG/jroux/Methylation/bismark/$flowcell\_$sampleID/CpG_OT_$sampleID\_trimmed.fq.gz_bismark.deduplicated.txt; perl collapse_strands.pl /KG/jroux/Methylation/bismark/$flowcell\_$sampleID/$sampleID /KG/jroux/Methylation/bismark/$flowcell\_$sampleID/$sampleID\_trimmed.fq.gz_bismark.deduplicated.bedGraph.gz /KG/jroux/Methylation/bismark/$flowcell\_$sampleID/CpG_OT_$sampleID\_trimmed.fq.gz_bismark.deduplicated.txt.gz; mv /KG/jroux/Methylation/bismark/$flowcell\_$sampleID/* $flowcell\_$sampleID/; rmdir /KG/jroux/Methylation/bismark/$flowcell\_$sampleID/\" | qsub -l h_vmem=16g -v PATH -cwd -N meth_$sampleID -o $flowcell\_$sampleID/$sampleID\_KG.out -e $flowcell\_$sampleID/$sampleID\_KG.err";
      }
    }
    ## this is never going to happen since a small collapsed bedgraph file seem to be created...
    elsif ((-s "$flowcell\_$sampleID/$sampleID\_trimmed.fq.gz_bismark.deduplicated.bedGraph.gz") and !(-s "$flowcell\_$sampleID/$sampleID.bedGraph.gz")){
      print "There was a problem with script to collapse strands. It will be relaunched.\n";
      if ($kg eq 0){
        $command = "echo \"perl collapse_strands.pl $flowcell\_$sampleID/$sampleID $flowcell\_$sampleID/$sampleID\_trimmed.fq.gz_bismark.deduplicated.bedGraph.gz $flowcell\_$sampleID/CpG_OT_$sampleID\_trimmed.fq.gz_bismark.deduplicated.txt.gz\" | qsub -l h_vmem=16g -v PATH -cwd -N meth_$sampleID -o $flowcell\_$sampleID/$sampleID.out -e $flowcell\_$sampleID/$sampleID.err";
      }
      else {
        $command = "echo \"mkdir -p /KG/jroux/Methylation/bismark/$sample; perl collapse_strands.pl /KG/jroux/Methylation/bismark/$flowcell\_$sampleID/$sampleID $flowcell\_$sampleID/$sampleID\_trimmed.fq.gz_bismark.deduplicated.bedGraph.gz $flowcell\_$sampleID/CpG_OT_$sampleID\_trimmed.fq.gz_bismark.deduplicated.txt.gz; mv /KG/jroux/Methylation/bismark/$flowcell\_$sampleID/* $flowcell\_$sampleID/; rmdir /KG/jroux/Methylation/bismark/$flowcell\_$sampleID/\" | qsub -l h_vmem=16g -v PATH -cwd -N meth_$sampleID -o $flowcell\_$sampleID/$sampleID\_KG.out -e $flowcell\_$sampleID/$sampleID\_KG.err";
      }
    }
    else {
      print "Couldn't identify where the problem was... Can you check manually?\n";
    }
    print "Command submitted to the cluster:\n$command\n";

    # add to report file
    open(OUT, '>>', "$flowcell\_$sampleID/$sampleID.report") or die 'Cannot open OUT file';
    print OUT "\nCommand line submitted to the cluster (resumed deduplication and methylation extractor):\n$command\n\n";
    close OUT;

    #Launch the command
    system($command)==0
        or warn "Failed to submit the job to the cluster\n";
  }
}
exit;
