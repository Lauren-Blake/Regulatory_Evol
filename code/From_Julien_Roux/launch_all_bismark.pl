#!/usr/bin/perl
use strict;
use warnings;

## Julien Roux
## May 31, 2012
## adapted fron launch_all_tophat.pl
## Recursively execute launch_bismark.pl for all samples listed on the file SamplesDirectories.txt that have not been mapped before

# look which samples were already mapped
my %mapped;
opendir(DIR, '.') or die "can't opendir: $!";
while (defined(my $sample_dir = readdir(DIR))) {
  if ($sample_dir =~ m/(Methylation\_\d+)\_(.+)/){
    my $flowcell = $1;
    my $sample = $2;
    $mapped{$sample}->{$flowcell}++;
    print "Sample $sample in flow-cell $flowcell seem to be already mapped. It will not be remapped here.\n";
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
  next if (exists $mapped{$tmp[3]}->{$tmp[0]});
  ## sample -> flow cell = fastq file path
  $samples{$tmp[3]}->{$tmp[0]} = "../raw_data/".$tmp[0].'/s_'.$tmp[2].'_sequence.txt.gz';
}
close IN;


foreach my $sample (sort keys %samples){
  foreach my $flowcell (sort keys %{$samples{$sample}}){
    # launch the job
    if ($sample =~ m/^H/) {     ## human
      # print "perl launch_bismark.pl $sample $flowcell /home/jroux/Methylation/raw_data/SamplesDirectories.txt /home/jroux/data/indexes/Bowtie/hg19_norandom_noM_noY+Lambda_prepared_bismark/\n";
      system("perl launch_bismark.pl $sample $flowcell /home/jroux/Methylation/raw_data/SamplesDirectories.txt /home/jroux/data/indexes/Bowtie/hg19_norandom_noM_noY+Lambda_prepared_bismark/")==0
          or warn "Failed to launch sample $sample.\n";
    }
    elsif ($sample =~ m/^C/) {  ## chimp
      # print "perl launch_bismark.pl $sample $flowcell /home/jroux/Methylation/raw_data/SamplesDirectories.txt /home/jroux/data/indexes/Bowtie/panTro3_nonrandom+Lambda_prepared_bismark/\n";
      system("perl launch_bismark.pl $sample $flowcell /home/jroux/Methylation/raw_data/SamplesDirectories.txt /home/jroux/data/indexes/Bowtie/panTro3_nonrandom+Lambda_prepared_bismark/")==0
          or warn "Failed to launch sample $sample.\n";
    }
    elsif ($sample =~ m/^R/) {  ## rhesus macaque
      # print "perl launch_bismark.pl $sample $flowcell /home/jroux/Methylation/raw_data/SamplesDirectories.txt /home/jroux/data/indexes/Bowtie/rheMac2_norandom_noUr+Lambda_prepared_bismark/\n";
      system("perl launch_bismark.pl $sample $flowcell /home/jroux/Methylation/raw_data/SamplesDirectories.txt /home/jroux/data/indexes/Bowtie/rheMac2_norandom_noUr+Lambda_prepared_bismark/")==0
          or warn "Failed to launch sample $sample.\n";
    }
  }
}
exit;


