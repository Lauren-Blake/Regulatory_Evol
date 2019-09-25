#!/usr/bin/perl
use strict;
use warnings;

## Julien Roux
## Mar 25, 2014
## Map samples with tophat after adapter removal with TrimGalore
## take name of sample to map (e.g. C1Li) as first argument
## retrieve original fastq files + rescued reads for all technical replicates of this sample using info in file in second argument
## Usually this file is ../raw_data/SamplesDirectories.txt
## map them to the trasncriptome and then genome of the 3 species using tophat (and qsub)

scalar @ARGV >= 4 or die<<USAGE;
Usage:
perl launch_tophat.pl <Sample> <Info_file> <Genome> <Transcriptome> <Mismatches_number_to_genome>
USAGE
# e.g. perl launch_tophat.pl H1Li ../raw_data/SamplesDirectories.txt /mnt/lustre/home/jroux/data/indexes/Bowtie2/Homo_sapiens.GRCh37.71.dna.toplevel ~/data/indexes/Bowtie2/Homo_sapiens.GRCh37.71.gtf 3

my $sample = $ARGV[0];
my $file_info = $ARGV[1];
my $genome = $ARGV[2];
my $transcriptome = $ARGV[3];
my $mismatch = 3; #parameter --initial-read-mismatches
if (defined $ARGV[4]){ $mismatch = $ARGV[4] };

# test if the sample was already mapped (to prevents erasing the previous files)
my $flag = 0;
opendir(DIR, '.') or die "can't opendir: $!";
while (defined(my $sample_dir = readdir(DIR))) {
  if ($sample_dir eq $sample){
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
my %sample;
my %index; ##just to check it is the same index used in all these technical replicates
open(IN, '<', $file_info) or die ("Cannot open info file\n");
my $header = <IN>;
while (defined (my $line = <IN>)) {
  chomp $line;
  my @tmp = split("\t", $line);
  next if ($tmp[1] ne 'RNA-seq'); # we don't want to treat the bisulfite data

  if ($tmp[4] eq $sample){
    # path to the data
    my $path = $tmp[9];
    if ($tmp[0] eq 'RNA-seq_4'){ ## this flow cell was sequenced by the core so the file path and file names are a bit special
      $sample{$raw_data_path.'/'.$path.$tmp[3].'_'.$tmp[7].'_L00'.$tmp[2].'_R1_001.fastq.gz'}->{'flow_cell'} = $tmp[0];
      $sample{$raw_data_path.'/'.$path.$tmp[3].'_'.$tmp[7].'_L00'.$tmp[2].'_R1_001.fastq.gz'}->{'lane'} = $tmp[2];
      $sample{$raw_data_path.'/'.$path.$tmp[3].'_'.$tmp[7].'_L00'.$tmp[2].'_R1_001.fastq.gz'}->{'index'} = $tmp[7];
      # path to the rescued dataa
      $sample{$raw_data_path.'/Rescued/'.$path.$tmp[3].'_'.$tmp[7].'_L00'.$tmp[2].'_R1_001.fastq.gz'}->{'flow_cell'} = $tmp[0];
      $sample{$raw_data_path.'/Rescued/'.$path.$tmp[3].'_'.$tmp[7].'_L00'.$tmp[2].'_R1_001.fastq.gz'}->{'lane'} = $tmp[2];
      $sample{$raw_data_path.'/Rescued/'.$path.$tmp[3].'_'.$tmp[7].'_L00'.$tmp[2].'_R1_001.fastq.gz'}->{'index'} = $tmp[7];
    }
    else { ## flow cells 1, 2 and 3
      $sample{$raw_data_path.'/'.$path.'s_'.$tmp[2].'_sequence.txt.gz'}->{'flow_cell'} = $tmp[0];
      $sample{$raw_data_path.'/'.$path.'s_'.$tmp[2].'_sequence.txt.gz'}->{'lane'} = $tmp[2];
      $sample{$raw_data_path.'/'.$path.'s_'.$tmp[2].'_sequence.txt.gz'}->{'index'} = $tmp[7];
      # path to the rescued data
      $sample{$raw_data_path.'/Rescued/'.$path.'s_'.$tmp[2].'_sequence.txt.gz'}->{'flow_cell'} = $tmp[0];
      $sample{$raw_data_path.'/Rescued/'.$path.'s_'.$tmp[2].'_sequence.txt.gz'}->{'lane'} = $tmp[2];
      $sample{$raw_data_path.'/Rescued/'.$path.'s_'.$tmp[2].'_sequence.txt.gz'}->{'index'} = $tmp[7];
    }
    # record the index sequence
    $index{$tmp[7]}++;
  }
}
close IN;

my $index;
if (scalar keys %index gt 1){
  die "More than one index was used for different technical replicates!\n";
}
else {
  $index = $index{(keys %index)[0]};
}

#retrieve adapter sequence for cutadapt
my $adapter;
open(IN, '<', "../raw_data/adapters.txt") or die ("Cannot open adapter file\n");
while (defined (my $line = <IN>)) {
  chomp $line;
  my @tmp = split("\t", $line);
  # check correspondence of adapter number and sequence
  if ($tmp[0] eq $index){
    $adapter = $tmp[1];
  }
}
close IN;
if (!defined $adapter){
  die "Adapter undefined!\n";
}

# create the output directory if it doesn't exist
my $outdir = "./".$sample.'/';
unless(-d $outdir){
  system("mkdir -p $outdir");
}

#create the string of all fastq files to map
my $all_files;
foreach my $fastq (sort keys %sample){
   $all_files .= $fastq." ";
}
$all_files = substr($all_files, 0, -1); #remove last character " "
# print "$all_files\n";

## create the string used as command
# * quality trim + cutadapt: trim_galore
#   trim_galore modified to use cutadapt with option -b to filter adapter sequences on 3' and 5' ends
#   remove end sequence if 3 letters similar to adapter (default is agressive: 1)
#   remove all read if remaining lenght < 20bp
#   Not possible to pipe input into trimgalore: need to create a temporary fast.gz file
#   FASTQC is run after trimming and cutadapt: removed base composition biases

# * Launch get_stats_piped.pl script to extract basic stats from finished tophat run
# * Use "qsub -pe simple_pe 4" to request multi-threading (tophat -p 4): this is actually not great since this submission queue is saturated: removed this option
# * Use -q batch64.q to prevent to be run on bigmem machines (full disk apparently)

# Adjust tophat parameters:
# * Most important: "--read-mismatches" ("-N") and --read-edit-dist
#   I was using 3 mismatches because I found this gave the best mapping. Tomas and Irene think it's too much: we'll use 2 mismatches
# * -g 1: report uniquely mapped reads
# * swith to --bowtie-n mode (mismatches counted in seed region only)
# * -p 4 multithreading. Set to 8 doesn't seem to go faster
# * -m 1: 1 mismatch allowed in splice anchor region: is 2 better? No!
# * --library-type fr-unstranded to specify library type
# * tophat version 2:
#       uses Bowtie 2
#       maps to transcriptome first (-n 3 max mismatches. Is it still used?): Add option -G GTF_file for the first run using this transcriptome file to build the index
# * -x 1: do not use this! (unique mapping to transcriptome). We must allow multiple hits to different transcripts of a same gene
# * -I 50000: max intron size. See Hong et al. (2006) Mol Biol Evol, 23, 2392-2404.
# * --b2-very-sensitive: removed because did not change anything
# * --b2-very-fast: It is indeed faster and the mapping is still as good. A few less junctions are found.
# * --b2-N 1: didn't improve alignement in my previous tests!
# * --coverage-search: Enables the coverage based search for junctions. Use when coverage search is disabled by default (such as for reads 75bp or longer), for maximum sensitivity. The mapping does not seem to be better, but ~1000 more junctions are found. The runtime is much longer (X2)!
# * -a  anchor size: Default is 8. Increase to 15? http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjourneal.pone.0013875. Maybe use in conjucntion with -m 2? no this doesn't change anything!
# * --no-novel-indels: doesn't seem to work
# * --read-mismatches 6??
# * quality: the lab sequencer is Phred-64 encoded (option --solexa1.3-quals), while the core facility is Phred-33 (no need for any option): we have to convert these to the same format
# SUMMARY: important parameters are mapping to transcriptome, using tophat2, very fast, and coverage (for more junctions but slower!)

# /!\ because the reads are 50bp long, the coverage search is enabled, so we disable it with --no-coverage-search (it's taking too long!)

## building of the command line
my $command = "echo \"cat $all_files > $outdir$sample.fastq.gz; /home/jroux/bin/trim_galore_v0.2.8/trim_galore_mod -q 20 --fastqc --phred64 -a $adapter --gzip --length 20 --output_dir $outdir --stringency 3 $outdir$sample.fastq.gz; /home/jroux/bin/tophat-2.0.8b.Linux_x86_64/tophat2 --solexa1.3-quals --library-type fr-unstranded --bowtie-n -g 1 -m 1 -N $mismatch --read-edit-dist $mismatch -n 3 -i 70 -I 50000 --b2-very-fast --transcriptome-index $transcriptome --no-coverage-search -o $outdir $genome $outdir$sample\_trimmed.fq.gz; samtools view $outdir/accepted_hits.bam | perl get_stats_piped.pl $outdir $sample; rm $outdir$sample.fastq.gz\" | qsub -l h_vmem=8g -q batch64.q -v PATH -cwd -N th_$sample -o $outdir$sample.out -e $outdir$sample.err\n";
#print "$command\n";

# print small report
open(OUT, '>', $outdir.$sample.'.report') or die 'Cannot open OUT file';
print OUT "Sample: $sample\nInformation file: $file_info\nMapped to transcriptome (3 mismatches): $transcriptome\nMapped to genome: $genome\nNumber of mismatches allowed to genome: $mismatch\n";
print OUT "Fastq files used for the mapping of $sample:\n";
foreach my $fastq (sort keys %sample){
  print OUT "  ", $fastq, "\n";
  print OUT "    With Index ", $sample{$fastq}->{'index'}, " on lane ", $sample{$fastq}->{'lane'}, " of flow cell ", $sample{$fastq}->{'flow_cell'}, "\n";
}
print OUT "Command line submitted to the cluster:\n$command\n";
close OUT;

## Launch the command
system($command)==0
    or warn "Failed to submit the job to the cluster\n";

exit;
