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
    print "Launching sample $sample in flow-cell $flowcell.\n";

    my $command = "echo \"perl check_mapping.pl $sample_dir\" | qsub -l h_vmem=1g -v PATH -q batch64.q -cwd -N check_$sample -o $sample_dir/$sample.out -e $sample_dir/$sample.err";
    #print "$command\n";

    system($command)==0
        or warn "Failed to launch sample $sample.\n";
  }
}
closedir(DIR);
exit;


