#!/usr/bin/perl
use strict;
use warnings;

## Julien Roux
## Feb 7, 2014
## Relaunch bismark only (starting from trimmed file) for a single sample
## Store temp on KG (when disk is quite full on main cluster)
## /KG/jroux/Methylation/bismark/

## TO DO: write even final files to KG if this still doesn't work

scalar @ARGV >= 2 or die<<USAGE;
Usage:
perl relaunch_bismark.pl <Flow_cell_Sample> <Genome>
USAGE

my $sample_flow_cell = $ARGV[0];
$sample_flow_cell =~ m/(Methylation\_\d+)\_(.+)/;
my $sample = $2;
my $flowcell = $1;
my $file_info = '/home/jroux/Methylation/raw_data/SamplesDirectories.txt';
my $genome = $ARGV[1];

my $outdir = $sample_flow_cell.'/';
my $KGdir = '/KG/jroux/Methylation/bismark/'.$sample_flow_cell.'/';
unless(-d $KGdir){
  system("mkdir -p $KGdir");
}

my $command = "echo \"/home/jroux/bin/bismark_v0.8.1/bismark -n 1 --path_to_bowtie /home/jroux/bin/bowtie-1.0.0/ --solexa1.3-quals --unmapped --ambiguous --bam --temp_dir $KGdir -o $outdir $genome $outdir$sample\_trimmed.fq.gz; perl get_stats.pl $outdir $sample; gzip -f -9 $outdir$sample\_trimmed.fq.gz_ambiguous_reads.txt; gzip -f -9 $outdir$sample\_trimmed.fq.gz_unmapped_reads.txt\" | qsub -l h_vmem=48g -v PATH -q batch64.q -cwd -N rebis_$sample -o $outdir$sample.out -e $outdir$sample.err";
#print "$command\n";

# print small report
open(OUT, '>>', $outdir.$sample.'.report') or die 'Cannot open OUT file';
print OUT "Relaunching mapping (after trimming and fastqc)...\nSample: $sample\nFlow-cell: $flowcell\nInformation file: $file_info\nMapped to genome: $genome\n";
print OUT "Command line submitted to the cluster:\n$command\n\n";
close OUT;

#Launch the command
system($command)==0
    or warn "Failed to submit the job to the cluster\n";
exit;

