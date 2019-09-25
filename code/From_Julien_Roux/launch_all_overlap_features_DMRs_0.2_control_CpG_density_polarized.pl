#!/usr/bin/perl
use strict;
use warnings;

## Julien Roux
## Oct 30 2014
# - Run on spudhead
# - Launches overlap_features_DMRs_0.2.pl for all DMRs, and all the random DMRs files of same length / CpG density
# - takes into account the direction of the methylation change: only possible for tissue-specific and species-specific comparisons.

scalar @ARGV >= 1 or die<<USAGE;
Usage:
perl launch_all_overlap_features_DMRs_0.2_control_CpG_density.pl <DMR type>
USAGE
my $DMRtype = $ARGV[0];
print "Launching overlap of genomic features for $DMRtype DMRs\n";

# first create BED files of DMR positions
my @files = glob("../bsseq/DMRs/$DMRtype/*Specific*DMRs.txt");
### my @files = glob("../bsseq/DMRs/$DMRtype/Chimp_heartSpecific*DMRs.txt");

foreach my $file (@files){
  foreach my $direction (('hypo', 'hyper')){
    print "$file $direction\n";
    $file =~ m/\.\.\/bsseq\/DMRs\/$DMRtype\/(.+)\_DMRs\.txt/;
    my $baseName = $1;
    my $DMRfile = $baseName.'_'.$direction.'_DMRs';

    # Export BED file
    # BED end coordinate are 0-based exclusive (== 1-based)

    # 'sed -n -e 2p -e 4p $file' Selects 2nd and 4th lines
    # system("sed 1d $file | sed -n". $positions{$direction}." | awk 'BEGIN{OFS=\"\t\";}{ print \$1,\$2,\$3+1 }' > ./DMRs/$DMRtype/$DMRfile.bed");
    ## Doesn't work when number of arguments is too big! Use perl instead:

    open(OUT, ">./DMRs/$DMRtype/$DMRfile.bed") or die ("Cannot open BED file\n");
    open(IN, $file) or die ("Cannot open DMRs file\n");
    my %positions;
    my $header = <IN>;
    my $i = 1;
    while (defined (my $line = <IN>)) {
      chomp $line;
      my @tmp = split("\t", $line);
      if ($tmp[15] eq $direction){
        $positions{$direction}->{$i}++;
        print OUT "$tmp[0]\t$tmp[1]\t", $tmp[2]+1, "\n";
      }
      $i++;
    }
    close IN;
    close OUT;

    # launch script overlap_features_DMRs_0.2.pl
    system("echo \"perl overlap_features_DMRs_0.2.pl ./DMRs/$DMRtype/$DMRfile.bed\" | qsub -l h_vmem=4g -v PATH -cwd -N $DMRfile -o ./DMRs/$DMRtype/$DMRfile.out -e ./DMRs/$DMRtype/$DMRfile.err");

    # launch script overlap_features_DMRs_0.2.pl for 100 sets of random DMRs of same length/CpG density
    my $randomDMRfile = $baseName.'_'.$direction.'_randomDMRs_';

    for (my $i=1; $i <= 100; $i++) {
      # Export BED file
      ## system("sed -n". $positions{$direction}." ./DMRs/$DMRtype/randomized_control_CpG_density/$baseName\_randomDMRs\_$i.bed > ./DMRs/$DMRtype/randomized_control_CpG_density/$randomDMRfile$i.bed");
      ## Doesn't work when number of arguments is too big! This should work better:

      open(OUT, ">./DMRs/$DMRtype/randomized_control_CpG_density/$randomDMRfile$i.bed") or die ("Cannot open BED file\n");
      open(IN, "./DMRs/$DMRtype/randomized_control_CpG_density/$baseName\_randomDMRs\_$i.bed") or die ("Cannot open DMRs file\n");
      my $j = 1;
      while (defined (my $line = <IN>)) {
        chomp $line;
        if (exists $positions{$direction}->{$j}){
          my @tmp = split("\t", $line);
          print OUT "$tmp[0]\t$tmp[1]\t$tmp[2]\n";
        }
        $j++;
      }
      close IN;
      close OUT;

      system("echo \"perl overlap_features_DMRs_0.2.pl ./DMRs/$DMRtype/randomized_control_CpG_density/$randomDMRfile$i.bed\" | qsub -l h_vmem=4g -v PATH -cwd -N $randomDMRfile$i -o ./DMRs/$DMRtype/randomized_control_CpG_density/$randomDMRfile$i.out -e ./DMRs/$DMRtype/randomized_control_CpG_density/$randomDMRfile$i.err");
    }
  }
  # last; # testing
}

