#!perl -w
use strict;
# Julien Roux, Apr 28 2015
# generate BEDgraph files for each sample
# - see http://genome.ucsc.edu/goldenPath/help/bedgraph.html for definition of format
# - zero-based half-open format
# - 19 49303800 49303801 10 means coverage of 10 on the base number 49303801 (1-based system)
# - regions with 0 coverage are not represented
# Also generate BigWig files from BEDgraph files for each sample
# https://genome.ucsc.edu/goldenPath/help/wiggle.html

## get the path to all samples to convert
my @samples = glob("./*/accepted_hits.bam");
# push(@samples, glob("./*/accepted_hits_subsample_p*.bam"));
# TO DO? Add those?

my %samples;
foreach my $path (@samples){
  $path =~ m/^\.\/(.+)\/accepted_hits.+$/;
  my $sample = $1;

  next if ((defined $ARGV[0]) and ($sample ne $ARGV[0])); ## if we want to relaunch only 1 sample

  $samples{$sample}->{$path}++;
  # print "$sample $path \n";
}

# launch BEDtools and script bedGraphToBigWig
SAMPLE:
foreach my $sample (sort keys %samples){
  print "Processing $sample\n";
  foreach my $input (sort keys %{$samples{$sample}}){


    my $output_bed = $input;
    $output_bed =~ s/accepted_hits(.*)\.bam/genome_coverage$1\.bed/;

    my $output_wig = $input;
    $output_wig =~ s/accepted_hits(.*)\.bam/genome_coverage$1\.bw/;

    my $outdir = $input;
    $outdir =~ s/accepted_hits.*\.bam//;

    ## check that the output file doesn't already exist and is not empty
    if ((-s $output_bed) or (-s $output_bed.".gz")){
      warn "The bedGraph files seem to have already been generated!\n";
      next SAMPLE;
    }

    my $genome;
    if ($sample =~ m/^H/){
      $genome = "/mnt/lustre/home/jroux/data/bedtools/hg19_chrs_length.genome";
    }
    elsif ($sample =~ m/^C/){
      $genome = "/mnt/lustre/home/jroux/data/bedtools/panTro3_chrs_length.genome";
    }
    elsif ($sample =~ m/^R/){
      $genome = "/mnt/lustre/home/jroux/data/bedtools/rheMac2_chrs_length.genome";
    }

    ## we could use the option -trackline of bedtools, but we want to control what is the name of the track, so we put it in the file first hand
    ## we use option -bga to output also positions with 0 coverage. Without this, many downstream programs tought there was no data at these locations ans did not consider them :(
    my $command = "echo track type=bedGraph name=\'\"".$sample." coverage\"\' description=\'\"".$sample." coverage\"\' priority=20 > $output_bed; echo \"/data/tools/bedtools-2.17/genomeCoverageBed -bga -split -ibam $input >> $output_bed; /mnt/lustre/home/jroux/bin/bedGraphToBigWig $output_bed $genome $output_wig; gzip -9 $output_bed\" | qsub -l h_vmem=8g -v PATH -cwd -N bed_$sample -o ".$outdir."genome_coverage.out -e ".$outdir."genome_coverage.err\n";
    #print "$command\n";

    ## Launch the command
    system($command)==0
        or warn "Failed to submit the job to the cluster\n";
  }
}

## Note: the export of the track header doesn't seem to be reliable
## Also check that the first line start at chromosome 1 position 0
##      for i in `find . -name genome_coverage.bed.gz`; do echo $i; zcat $i | head -n2; done

## check the end of the file
## for i in `find . -name genome_coverage.bed.gz`; do echo $i; zcat $i | tail -n1; done

## check the size of all files and the good completion of the jobs
# for i in `find . -name genome_coverage.err`; do ll $i; done
# for i in `find . -name genome_coverage.bed.gz`; do ll $i; done
# for i in `find . -name genome_coverage.bw`; do ll $i; done
