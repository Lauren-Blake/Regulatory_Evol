#!/usr/bin/perl
use strict;
use warnings;

## Julien Roux
## Aug 27, 2014
## Launch collapse_evidence_strands.pl and identify_strand_problems.pl

print "Reading information file...\n";
open(IN, '<', "../raw_data/SamplesDirectories.txt") or die ("Cannot open info file\n");
my %samples;
my $header = <IN>;
while (defined (my $line = <IN>)) {
  chomp $line;
  my @tmp = split("\t", $line);
  next if ($tmp[1] ne 'BS-seq'); # we don't want to treat the RNA-seq data

  my $sample = $tmp[3];
  ## next if ($sample ne 'C1H2');  # testing:

  my $flowcell = $tmp[0];
  print "Scripts will be launched on sample $sample in flow-cell $flowcell\n";
  ## sample -> flow cell++
  $samples{$sample}->{$flowcell}++;
}
close IN;

my $command;
foreach my $sample (sort keys %samples){
  foreach my $flowcell (sort keys %{$samples{$sample}}){
    # launch the job
    if ($sample =~ m/^H/){ 
      # we do not launch mapping after bedGraph files generation
      $command = "echo \"perl collapse_evidence_strands.pl $flowcell\_$sample/$sample $flowcell\_$sample/$sample\_trimmed.fq.gz_bismark.deduplicated.bedGraph.gz $flowcell\_$sample/CpG_OT_$sample\_trimmed.fq.gz_bismark.deduplicated.txt.gz; perl identify_strand_problems.pl $flowcell\_$sample/$sample $flowcell\_$sample/$sample\_trimmed.fq.gz_bismark.deduplicated.bedGraph.gz $flowcell\_$sample/CpG_OT_$sample\_trimmed.fq.gz_bismark.deduplicated.txt.gz; \" | qsub -l h_vmem=16g -v PATH -cwd -N collapse_$sample -o $flowcell\_$sample/$sample.out -e $flowcell\_$sample/$sample.err";
    }
    else {
      my $species;
      if ($sample =~ m/^C/){ $species = 'chimp'; }
      if ($sample =~ m/^R/){ $species = 'rhesus'; }
      $command = "echo \"perl collapse_evidence_strands.pl $flowcell\_$sample/$sample $flowcell\_$sample/$sample\_trimmed.fq.gz_bismark.deduplicated.bedGraph.gz $flowcell\_$sample/CpG_OT_$sample\_trimmed.fq.gz_bismark.deduplicated.txt.gz; perl identify_strand_problems.pl $flowcell\_$sample/$sample $flowcell\_$sample/$sample\_trimmed.fq.gz_bismark.deduplicated.bedGraph.gz $flowcell\_$sample/CpG_OT_$sample\_trimmed.fq.gz_bismark.deduplicated.txt.gz; perl map_to_hg19.pl $flowcell\_$sample/$sample\.unmethylated_evidence.bedGraph.gz ../CpG_maps_to_hg19/by_chr_qsub/ $species; sort -k1,1 -k2,2n $flowcell\_$sample/$sample\.unmethylated_evidence_mapped.bedGraph | gzip -9 > $flowcell\_$sample/$sample\.unmethylated_evidence_mapped.bedGraph.gz; rm $flowcell\_$sample/$sample\.unmethylated_evidence_mapped.bedGraph; sort -k1,1 -k2,2n $flowcell\_$sample/$sample\.unmethylated_evidence_mapped_CpG.bedGraph | gzip -9 > $flowcell\_$sample/$sample\.unmethylated_evidence_mapped_CpG.bedGraph.gz; rm $flowcell\_$sample/$sample\.unmethylated_evidence_mapped_CpG.bedGraph; perl map_to_hg19_problematic_positions.pl $flowcell\_$sample/$sample\.problematic_positions.bedGraph.gz ../CpG_maps_to_hg19/by_chr_qsub/ $species; sort -k1,1 -k2,2n $flowcell\_$sample/$sample\.problematic_positions_mapped.bedGraph | gzip -9 > $flowcell\_$sample/$sample\.problematic_positions_mapped.bedGraph.gz; rm $flowcell\_$sample/$sample\.problematic_positions_mapped.bedGraph; sort -k1,1 -k2,2n $flowcell\_$sample/$sample\.problematic_positions_mapped_CpG.bedGraph | gzip -9 > $flowcell\_$sample/$sample\.problematic_positions_mapped_CpG.bedGraph.gz; rm $flowcell\_$sample/$sample\.problematic_positions_mapped_CpG.bedGraph;\" | qsub -l h_vmem=16g -v PATH -cwd -N collapse_$sample -o $flowcell\_$sample/$sample.out -e $flowcell\_$sample/$sample.err";

      # mapping alone
      #$command = "echo \"perl map_to_hg19.pl $flowcell\_$sample/$sample\.unmethylated_evidence.bedGraph.gz ../CpG_maps_to_hg19/by_chr_qsub/ $species; sort -k1,1 -k2,2n $flowcell\_$sample/$sample\.unmethylated_evidence_mapped.bedGraph | gzip -9 > $flowcell\_$sample/$sample\.unmethylated_evidence_mapped.bedGraph.gz; rm $flowcell\_$sample/$sample\.unmethylated_evidence_mapped.bedGraph; sort -k1,1 -k2,2n $flowcell\_$sample/$sample\.unmethylated_evidence_mapped_CpG.bedGraph | gzip -9 > $flowcell\_$sample/$sample\.unmethylated_evidence_mapped_CpG.bedGraph.gz; rm $flowcell\_$sample/$sample\.unmethylated_evidence_mapped_CpG.bedGraph; perl map_to_hg19_problematic_positions.pl $flowcell\_$sample/$sample\.problematic_positions.bedGraph.gz ../CpG_maps_to_hg19/by_chr_qsub/ $species; sort -k1,1 -k2,2n $flowcell\_$sample/$sample\.problematic_positions_mapped.bedGraph | gzip -9 > $flowcell\_$sample/$sample\.problematic_positions_mapped.bedGraph.gz; rm $flowcell\_$sample/$sample\.problematic_positions_mapped.bedGraph; sort -k1,1 -k2,2n $flowcell\_$sample/$sample\.problematic_positions_mapped_CpG.bedGraph | gzip -9 > $flowcell\_$sample/$sample\.problematic_positions_mapped_CpG.bedGraph.gz; rm $flowcell\_$sample/$sample\.problematic_positions_mapped_CpG.bedGraph;\" | qsub -l h_vmem=8g -v PATH -cwd -N collapse_$sample -o $flowcell\_$sample/$sample.out -e $flowcell\_$sample/$sample.err";
    #system($command)==0
    #    or warn "Failed to submit the job to the cluster\n";
    }

    #Launch the command
    # print "$command\n";
    system($command)==0
        or warn "Failed to submit the job to the cluster\n";
  }
}
exit;
