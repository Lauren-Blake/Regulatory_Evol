#!/usr/bin/perl
use strict;
use warnings;

## Julien Roux
## May 24, 2013
## Map BS-seq data to converted genomes using Bismark
## take name of flow-cell + technical replicate to map (e.g. Methylation_1 C1Li1) as arguments

scalar @ARGV >= 4 or die<<USAGE;
Usage:
perl launch_bismark.pl <Sample> <Flow cell> <Info_file> <Genome>
USAGE
# e.g. perl launch_bismark.pl H1Lu1 Methylation_4 ../raw_data/SamplesDirectories.txt ../../data/indexes/Bowtie/hg19_norandom_noM_noY+Lambda_prepared_bismark/

my $sample = $ARGV[0];
my $flowcell = $ARGV[1];
my $file_info = $ARGV[2];
my $genome = $ARGV[3];

# test if the sample was already mapped (to prevent erasing the previous files)
my $flag = 0;
opendir(DIR, '.') or die "can't opendir: $!";
while (defined(my $sample_dir = readdir(DIR))) {
  if ($sample_dir eq $flowcell."_".$sample){
    $flag = 1;
  }
}
closedir(DIR);
if ($flag eq 1){
  die "Sample $sample was already mapped previously. Please remove or rename directory.\n";
}

# path to the raw data
$file_info =~ m/(.+)\/\w+\.txt/;
my $raw_data_path = $1;


print "Reading information file...\n";
my $file;
my $fastqc;
my $adapter;
open(IN, '<', $file_info) or die ("Cannot open info file\n");
my $header = <IN>;
while (defined (my $line = <IN>)) {
  chomp $line;
  my @tmp = split("\t", $line);
  next if ($tmp[1] ne 'BS-seq'); # we want to read only the BS-seq rows

  if (($tmp[3] eq $sample) and ($tmp[0] eq $flowcell)){
    ## flow cell 5 was sequenced elsewhere (UChicago core facility) so the fastq files are splitted. We read the folder content to get them all.
    if ($flowcell eq 'Methylation_5'){ 
      opendir(DIR, $raw_data_path.'/'.$tmp[9]) or die "can't opendir: $!";
      while (defined(my $sub_fastq = readdir(DIR))) {
        if ($sub_fastq =~ m/fastq\.gz/){
          $file .= $raw_data_path.'/'.$tmp[9].'/'.$sub_fastq.' ';
          $sub_fastq =~ s/\.fastq.gz//;
          $fastqc .= '~/clusterhome/Methylation/FASTQC/'.$tmp[9].'/'.$sub_fastq.'_fastqc/fastqc_report.html; ';
        }
      }
      closedir(DIR);
    }
    else {
      $file = $raw_data_path.'/'.$flowcell.'/s_'.$tmp[2].'_sequence.txt.gz';
      $fastqc = '~/clusterhome/Methylation/FASTQC/'.$flowcell.'/s_'.$tmp[2].'_sequence_fastqc/fastqc_report.html';
    }
    $adapter = $tmp[8];
    # print "adapter: $adapter\n";
  }
}
close IN;

if (!defined $file){ die "This sample was not found on this flow cell (or something else wrong in the information file)\n"; }
if (!defined $adapter){ die "The adapter for this sample was not found on the information file\n"; }

my $adapter_seq;
open(IN, '<', "../raw_data/adapters.txt") or die ("Cannot open adapter file\n");
while (defined (my $line = <IN>)) {
  chomp $line;
  my @tmp = split("\t", $line);
  # search for sequence of adapter (we know its number)
  if ($tmp[0] eq $adapter){
    $adapter_seq = $tmp[1];
  }
}
close IN;
if (!defined $adapter_seq ){ die "Adapter undefined!\n"; }


# create the output directory if it doesn't exist
my $outdir = "./".$flowcell."_".$sample.'/';
unless(-d $outdir){
  system("mkdir -p $outdir");
}

## quality trim + cutadapt: trim_galore
##   cutadapt: -b to filter adapter sequences on 3' and 5' ends
##   remove end sequence if 3 letters similar to adapter (default is agressive: 1)
##   remove all read if remaining length < 15bp
##   FASTQC is run after trimming and cutadapt

## bismark
# /home/jroux/bin/bismark_v0.8.1/bismark -n 1 --path_to_bowtie /home/jroux/bin/bowtie-1.0.0/ --solexa1.3-quals --unmapped --ambiguous --bam -o ./ /home/jroux/data/indexes/Bowtie/hg19_norandom_noM_noY+Lambda_prepared_bismark/ /home/jroux/Methylation/raw_data/Methylation_1/s_5_sequence.txt.gz
## The tests indicated that mapping with bowtie1 is much faster
##  2 or 3 mismatches were the best mapping, but larger number of methylated CHG and CHH
##  1 mismatch was almost as good, less CHG and CHH methylated, and much faster
##  Bismakr with Bowtie2 was much slower and less reads mapped

my $command;
if ($flowcell eq 'Methylation_5'){ 
  # We create a temporary file and then remove it. This can potentially be problematic with the disk space on the cluster...
  $command = "echo \"cat $file > $outdir$sample.fastq.gz; /home/jroux/bin/trim_galore_v0.2.8/trim_galore_mod -q 20 --fastqc --phred33 -a $adapter_seq --gzip --length 15 --output_dir $outdir --stringency 3 $outdir$sample.fastq.gz; rm $outdir$sample.fastq.gz; touch $outdir$sample.fastq.gz; /home/jroux/bin/bismark_v0.8.1/bismark -n 1 --path_to_bowtie /home/jroux/bin/bowtie-1.0.0/ --phred33-quals --unmapped --ambiguous --bam --temp_dir $outdir -o $outdir $genome $outdir$sample\_trimmed.fq.gz; perl get_stats.pl $outdir $sample; gzip -f -9 $outdir$sample\_trimmed.fq.gz_ambiguous_reads.txt; gzip -f -9 $outdir$sample\_trimmed.fq.gz_unmapped_reads.txt\" | qsub -l h_vmem=48g -v PATH -cwd -N bis_$sample -o $outdir$sample.out -e $outdir$sample.err";
}
else {
  $command = "echo \"ln -s $file $outdir$sample.fastq.gz; /home/jroux/bin/trim_galore_v0.2.8/trim_galore_mod -q 20 --fastqc --phred64 -a $adapter_seq --gzip --length 15 --output_dir $outdir --stringency 3 $outdir$sample.fastq.gz; /home/jroux/bin/bismark_v0.8.1/bismark -n 1 --path_to_bowtie /home/jroux/bin/bowtie-1.0.0/ --solexa1.3-quals --unmapped --ambiguous --bam --temp_dir $outdir -o $outdir $genome $outdir$sample\_trimmed.fq.gz; perl get_stats.pl $outdir $sample; gzip -f -9 $outdir$sample\_trimmed.fq.gz_ambiguous_reads.txt; gzip -f -9 $outdir$sample\_trimmed.fq.gz_unmapped_reads.txt\" | qsub -l h_vmem=48g -v PATH -q batch64.q -cwd -N bis_$sample -o $outdir$sample.out -e $outdir$sample.err";
}
#print "$command\n";

# print small report
open(OUT, '>', $outdir.$sample.'.report') or die 'Cannot open OUT file';
print OUT "Sample: $sample\nFlow-cell: $flowcell\nInformation file: $file_info\nMapped to genome: $genome\n";
print OUT "Fastq file(s) used for the mapping of $sample: $file\nThis fastq file was linked to $outdir$sample.fastq.gz\n\n";
print OUT "Command line submitted to the cluster:\n$command\n\n";
print OUT "Fastqc file before quality and adapter trimming: $fastqc\n";
print OUT "Fastqc file after quality and adapter trimming: ~/clusterhome/Methylation/bismark/$outdir$sample\_trimmed.fq_fastqc/fastqc_report.html\n";
close OUT;

#Launch the command
system($command)==0
    or warn "Failed to submit the job to the cluster\n";
exit;

