#!/usr/bin/perl
use strict;
use warnings;

## Julien Roux
## July 6, 2012
## Recursively execute launch_tophat.pl for all samples listed on the file "samples_to_map.txt"
## (mapping them with tophat on the genome)

open(IN, '<', 'samples_to_map.txt') or die ("Cannot open info file\n");
while (defined (my $line = <IN>)) {
  chomp $line;
  next if $line eq '';
  next if $line =~ /^#/;

  # launch the job
  if ($line =~ m/^H/){ ## human
    system("perl launch_tophat.pl $line ../raw_data/SamplesDirectories.txt /mnt/lustre/home/jroux/data/indexes/Bowtie2/hg19_norandom_noM_noY ~/data/indexes/Bowtie2/refSeq_hg19.gtf 2")==0
        or warn "Failed to launch sample $line.\n";
  }
  elsif ($line =~ m/^C/){ ## chimp
    system("perl launch_tophat.pl $line ../raw_data/SamplesDirectories.txt /mnt/lustre/home/jroux/data/indexes/Bowtie2/panTro3_nonrandom ~/data/indexes/Bowtie2/refSeq_panTro3.gtf 2")==0
        or warn "Failed to launch sample $line.\n";
  }
  elsif ($line =~ m/^R/){ ## rhesus macaque
    system("perl launch_tophat.pl $line ../raw_data/SamplesDirectories.txt /mnt/lustre/home/jroux/data/indexes/Bowtie2/rheMac2_norandom_noUr ~/data/indexes/Bowtie2/refSeq_rheMac2.gtf 2")==0
        or warn "Failed to launch sample $line.\n";
  }
}
close IN;
exit;

