#!perl -w
use strict;
# Julien Roux, March 25 2014
# generate BAM files that are subsamples of full BAM
# downsampling to the minimum number of reads mapped in sample R1H: 28544039 reads

# Problem:
#  In samtools the sampling is not based on a number of reads, but on a probability, so it's hard to have an exact number of reads. Also the results do not seem to be reliable...
#  Bamtools random would be nice because it can choose a given number of reads, but it is too slow!
#  -> we use picard (probability too :( )

## get the path to all samples that were mapped
my @samples = glob("./*/accepted_hits.bam");

my %samples;
foreach my $path (@samples){
  $path =~ m/^\.\/(.*)\/accepted_hits.bam$/;
  my $sample = $1;
  $samples{$sample} = $path;
}

foreach my $sample (sort keys %samples){
  print "Processing $sample\n";
  my $outdir = $samples{$sample};
  $outdir =~ s/accepted_hits\.bam//;

  ## if smallest sample: just copy the file
  if ($sample eq "R1H"){
    system("cp $samples{$sample} $outdir/accepted_hits_subsample_p=1.bam");
  }
  else {
    ## read report file and extract number of mapped reads
    my $prob = '28544039'; ## minimum number of mapped reads in R1H
    open(IN, $outdir.$sample.".report") or die "Cannot read report file\n";
    while (defined (my $line = <IN>)) {
      chomp $line;
      if ($line =~ m/Total number of reads mapped: (\d+)/){
        $prob /= $1;
      }
    }
    print "\tPortion of kept reads: $prob\n";

    my $output = $samples{$sample};
    $output =~ s/accepted_hits\.bam/accepted_hits_subsample_p=$prob\.bam/;

    ## Does not seem to work on cluster: launched on a supdling node
    my $command = "java -Xmx512m -Xms512m -mx512m -jar /mnt/lustre/home/jroux/bin/picard-tools-1.88/DownsampleSam.jar INPUT=$samples{$sample} OUTPUT=$output PROBABILITY=$prob";
    # print "$command\n";
    ## Launch the command
    system($command)==0
        or warn "Failed to submit the job\n";
  }
}

