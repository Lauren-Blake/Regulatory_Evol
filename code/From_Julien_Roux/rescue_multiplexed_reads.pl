#!/usr/bin/perl
use strict;
use warnings;
$| = 1; # no buffering

## Julien Roux
## March 21, 2013
## Script inspired from Russ's script Internal_garbagecollector.pl v2 (Gilad lab wiki)

## for a given flow cell:
## find all the non demultiplexed fastq files in the raw data directory
## foreach lane, attribute non indexed reads to the closest index with a limit of difference fixed by the user in argument
## The rescued reads are written to files in Rescued directory (same structure as Demultiplexed directory)
## (This does not require to copy the original files to the user account)

scalar @ARGV == 2 or die<<USAGE;
Usage:
perl -w rescue_multiplexed_reads.pl <flow cell> <limit>
USAGE

my $flow_cell = $ARGV[0];
# maximum number of mismatches allowed for a successful assignment (1 usually)
my $distlimit = $ARGV[1];
print "Rescue reads with a mismatch distance limit of $distlimit\n";

## record all fastq files to use
print "Recording all reference indexes and file paths...\n";
my $raw_data = "/mnt/lustre/home/jroux/Methylation/raw_data/";

my %samples;
open(IN, '<', $raw_data.'SamplesDirectories.txt') or die "couldn't open IN file\n";
my $header = <IN>;
while (defined (my $line = <IN>)) {
  chomp $line;
  my @tmp = split("\t", $line);
  ## keep only selected flow cell
  next if ((defined $flow_cell) and ($tmp[0] ne $flow_cell));

  ## lane -> name = folder / index sequence
  if (defined $tmp[9]){
    $samples{$tmp[2]}->{$tmp[3]}->{'folder'} = $tmp[9];
    $samples{$tmp[2]}->{$tmp[3]}->{'index'} = $tmp[7];
  }
}
close IN;

## rescue the non demultiplexed reads lane by lane
foreach my $lane (keys %samples){
  ## used for stats
  my $count_reads = 0;
  my $count_rescued = 0;

  print "Rescuing: samples in lane $lane\n";

  ## we need to collect the directory where the non demultiplexed files are
  my $indir = ();

  ## open all out files (compressed)
  foreach my $sample (keys %{$samples{$lane}}){
    my $dir = $samples{$lane}->{$sample}->{'folder'};
    unless (defined $indir){
      $indir = $dir;
      if ($flow_cell eq 'RNA-seq_4'){ ## this flow cell was sequenced by the core facility, the path and files names are special
        $indir =~ s{$flow_cell\/Sample_YG-JR-S\d+\/}{$flow_cell\/Undetermines_Indices\/};
      }
      else {
        $indir =~ s{$flow_cell\/\d{3}\/(.+)}{$flow_cell\/unknown\/$1};
      }
    }
    $dir = 'Rescued/'.$dir;
    unless(-d $dir){
      system("mkdir -p $dir") or print "mkdir failed\n";
    }

    my $filename;
    if ($flow_cell eq 'RNA-seq_4'){ ## this flow cell was sequenced by the core facility, the path and files names are special
      $filename = $dir.$sample.'\_'.$samples{$lane}->{$sample}->{'index'}.'\_L00'.$lane.'\_R1_001.fastq.gz';
    }
    else {
      $filename = $dir.'s_'.$lane.'_sequence.txt.gz';
    }
    open_fh($filename, $samples{$lane}->{$sample}->{'index'});
  }

  ## open file of reads to rescue
  if ($flow_cell eq 'RNA-seq_4'){ ## this flow cell was sequenced by the core facility, the path and files names are special
    open(IN, 'zcat '.$indir.'lane'.$lane.'_Undetermined_L00'.$lane.'_R1_001.fastq.gz |') or die ("Cannot open fastq file\n");
  }
  else {
    open(IN, 'zcat '.$indir.'s_'.$lane.'_sequence.txt.gz |') or die ("Cannot open fastq file\n");
  }
  while (defined (my $line = <IN>)) {
    chomp $line;
    if ($line =~ m/^@/) {
      $count_reads++;
      ## print the script progress
      if ($count_reads % 1000000 eq 0) {
        print "\t".$count_reads." read parsed...\n";
      }

      ## record the index for this read
      my $indexseq;
      if ($flow_cell eq 'RNA-seq_4'){ ## this flow cell was sequenced by the core facility, the path and files names are special
        $line =~ m/\s1\:N\:0\:([ACGTN]{6})/;
        $indexseq = $1;
      }
      else {
        $line =~ m/\#([ACGTN]{6})\/1/;
        $indexseq = $1;
      }
      next if (!defined $indexseq);

      # find the closest index in the indexes used in this lane
      my $closest = checkseqs($indexseq, \%{$samples{$lane}}, $distlimit);
      next if ($closest eq 'NA'); #go to the next line for ambiguous sequences
      $count_rescued++;

      ## printing the rest of the fastq entry for this read:
      ## print the read name and add the closest real index to it
      my $read_name = $line."#takenAsIndex:$closest";
      print_fh($closest, $read_name);
      ## read and print the read sequence
      my $nextline = <IN>;
      chomp($nextline);
      print_fh($closest, $nextline);
      ## print the read name again starting with a +
      $nextline = <IN>;
      chomp($nextline);
      $read_name = $nextline."#takenAsIndex:$closest";
      if ($flow_cell eq 'RNA-seq_4'){ ## this flow cell was sequenced by the core facility, the path and files names are special
        print_fh($closest, '+');
      }
      else {
        print_fh($closest, $read_name);
      }
      ## print the quality of the read
      $nextline = <IN>;
      chomp($nextline);
      print_fh($closest, $nextline);
    }
  }
  close IN;

  ## close all out files
  foreach my $sample (keys %{$samples{$lane}}){
    close_fh($samples{$lane}->{$sample}->{'index'});
  }

  ## some stats
  print "\tAnalysis complete for samples in $lane. $count_rescued reads rescued out of $count_reads parsed.\n";
}
exit;

#subroutines
sub seqdist	{
  #reports the number of differences in sequence between a template and query sequence. 
  #USAGE: seqdist($query, $template)

	my($query, $template) = @_;

	my @q = split(//, $query);
	my @t = split(//, $template);

	my $diffs = 0;
	for(my $i = 0; $i < scalar(@t); $i++){
		if($q[$i] ne $t[$i]){ $diffs++ }
  }
	return($diffs);
}

sub checkseqs	{
  #uses seqdist to return the closest adapter ID, or "NA" if match is ambiguous.
  #USAGE: checkseqs($query, \%hash_with_samples_ref_indexes, $maxdist)

	my($query, $hash, $maxdist) = @_;
  my %ref_indexes = %{$hash};

	my $minimum_diff = ($maxdist + 1);
	my $out = 'NA';
  my $multiple = 0; # if multiple ref indexes have the same distance to the analyzed sequence

	foreach my $sample (keys %ref_indexes){
    my $ref_index = $ref_indexes{$sample}->{'index'};
		my $output = seqdist($query, $ref_index);
    ## if we previosuly seen a good match and this is another one: multiple matches!
    if (($output == $minimum_diff) and ($minimum_diff <= $maxdist)){
			$multiple = 1;
    }
		if ($output < $minimum_diff){
			$minimum_diff = $output;
			$out = $ref_index;
      $multiple = 0;
    }
  }
  if ($multiple eq 1){
    $out = 'NA';
  }
  return($out);
}

sub print_fh{	# print fastq elements to a specified filehandle
	local *FH = shift;
	my $line = shift;
	print FH "$line\n";
}

sub open_fh{	# open a specified file to append, with the supplied variable as the filehandle. USAGE: open_fh($file, $filehandle);
  my $outfile = shift;
	local *FH = shift;

  # open the output (compressed) file
	open(*FH, "| gzip - > $outfile") or die ("Cannot write rescued sequences to file\n");
}

sub close_fh{	# open a specified file to print with the supplied variable as the filehandle. USAGE: close_fh($filehandle);
	local *FH = shift;
	close(*FH);
}



















