#!/usr/bin/perl
use strict;
use warnings;
$| = 1;

## Julien Roux
## Oct 7, 2014
# - Run with lot of memory
# - only one DMR type at a time
# - Collect all DMRs files, export them in BED format and create BED files of random DMRs on same chromosome, of same length and CpG density

scalar @ARGV >= 1 or die<<USAGE;
Usage:
perl generate_random_DMRs_control_CpG_density.pl <DMR type>
USAGE
my $DMRtype = $ARGV[0];
print "Generating random DMRs for $DMRtype DMRs\n";

# first record the position of all CpGs in human genome (per chromosome). 35M CpGs in total
print "Recording CpG map...\n";
my %CpG_map;
open(IN, "zcat CpG_map/CpG_map.txt.gz |") or die "Can't read input file\n";
my $count = 0;
while (defined (my $line = <IN>)) {
  $count++;
  ## print the script progress
  if ($count % 1000000 eq 0) {
    print "  ".$count." CpGs parsed\n";
  }
  chomp $line;
  my @tmp = split("\t", $line);
  # chr -> array of sorted positions
  push(@{$CpG_map{$tmp[0]}}, $tmp[1]);
}
close IN;

# no need to generate randomDMRs that are more or less dense than any DMR in real dataset!
# So we record the range of density for all real DMRs (hash: all numbers of CpGs -> all lengths)
my %densityRange;
my $maxLength = 0;

# Create BED files of DMRs and record their length and number of CpG
my @files = glob("../bsseq/DMRs/$DMRtype/*_DMRs.txt");
my %DMRs;
foreach my $file (@files){
  print "Treating file: ", $file, "\n";
  # Export of BED file was already done by script launch_all_overlap_features_DMRs.pl

  # record DMR length and density
  open(IN, $file) or die ("Cannot open DMRs file\n");
  my $header = <IN>;
  while (defined (my $line = <IN>)) {
    chomp $line;
    my @tmp = split("\t", $line);
    my $length = $tmp[7];

    my $chr = $tmp[0];
    # In contrast Rhesus_1_3_DMRs.txt, one DMR was mapped to chrY in human: this is problematic since we don't have this chromosome in the reference genome. So we take random DMRs from chromosome X instead
    if ($chr eq 'chrY'){
      $chr = 'chrX'
    }
    # some chimp or rhesus DMRs are mapped to chr17_gl000204_random or chr11_gl000204_random
    if ($chr =~ m/(chr\d+)\_gl\d+\_random/){
      $chr = $1;
    }

    # record number of CpGs
    # the "n" column ($tmp[6]) in DMR file gives the number of CpG after filtering: not a good indicator for CpG density.
    # we need to use the CpG map to get the number of CpGs:
    my $numberCpG = 1; # at least 1 CpG on starting position
  CPG:
    for (my $i = 0; $i < scalar(@{$CpG_map{$chr}}); $i++) {
      # looking for starting position: stop as soon as positions >= starting position
      # (some starting positions could be on - strand, so we subtract 1)
      # (some starting positions mapped from chimp or rhesus could not be CpGs in human, so we start counting on the next CpG)
      if (${$CpG_map{$chr}}[$i] >= ($tmp[1]-1)){
        # print "Starting position found: index $i and position: ", ${$CpG_map{$chr}}[$i], "\n";

        # go to next positions until last position (just check that there is indeed a next position
        while ((($i + $numberCpG) <= $#{$CpG_map{$chr}}) and (${$CpG_map{$chr}}[$i + $numberCpG] <= $tmp[2])){
            #print "Passing through position :", ${$CpG_map{$chr}}[$i + $numberCpG], "\n";
          $numberCpG++;
        }
        last CPG; # no need to iterate through the end of the array
      }
    }
    # If the chimp or rhesus is mapped to a sequence that was strongly deaminated in human, the density will be underestimated. Thus is this case we use the number of CpGs that are inside DMR
    if ($tmp[6] > $numberCpG){
      $numberCpG = $tmp[6];
    }
    #print "DMR of length: ", $length, " and number of CpGs: ", $numberCpG, "\n";

    # record density for generation of randomDMRs
    # add to hash: number of CpG -> length 
    $densityRange{$numberCpG}->{$length} = ();
    if ($length > $maxLength){
      $maxLength = $length;
    }

    push(@{$DMRs{$file}}, [$chr, $length, $numberCpG]); ## build an hash of array
  }
  close IN;

  # last; #testing
}
print "Max length: $maxLength\n";

# Generate list of all possible randomDMRs
# For each CpG in genome, record coordinates of all DNA strings extending to CpGs on the right
# length < $maxLength
# only obsreved values of length and number of CpGs are stored
print "Creating list of all potential randomDMRs...\n";
my %randomDMRs;
foreach my $chr (keys %CpG_map){
  print "  $chr\n";

  # for each CpG position
  for (my $i = 0; $i < scalar(@{$CpG_map{$chr}}); $i++) {
    ## print the script progress
    if ($i % 100000 eq 0) {
      print "    ".$i." CpGs positions done\n";
    }
    #print ".";

    my $j = $i;
  EXTEND:
    # no need to extend longer than max length of any DMR
    # j is incremented in the loop
    while (($j < scalar(@{$CpG_map{$chr}})) and ((${$CpG_map{$chr}}[$j] - ${$CpG_map{$chr}}[$i] + 1) <= $maxLength)){
      #print ${$CpG_map{$chr}}[$j] - ${$CpG_map{$chr}}[$i] + 1, " ";
      # no need to generate random strings with a number of CpGs that is not seen in real DMRs
      next EXTEND unless (exists $densityRange{$j - $i + 1});

      # no need to generate random strings with a length that is not seen in real DMRs
      next EXTEND unless (exists $densityRange{$j - $i + 1}->{${$CpG_map{$chr}}[$j] - ${$CpG_map{$chr}}[$i] + 1});

      #chromosome -> number of CpGs -> length -> start = end
      $randomDMRs{$chr}->{$j - $i + 1}->{${$CpG_map{$chr}}[$j] - ${$CpG_map{$chr}}[$i] + 1}->{${$CpG_map{$chr}}[$i]} = ${$CpG_map{$chr}}[$j];
    }
    # Executed after each loop iteration and before the conditional statement is evaluated. A good place to increment counters.
    continue {
      $j++;
    };
  }
  #last; #testing
}
## TO DO? read through all hash and trim number of strings for some values or numberCpGs and length which have throusands of items (are there some?)

foreach my $file (sort keys %DMRs){
  print "Creating files of random DMRs for: ", $file, "\n";
  # Create 100 BED files of random DMRs of same CpG density

#   $file =~ m/\.\.\/bsseq\/DMRs\/(.+)\/(.+\_DMRs)\.txt/;
#   # DMR type -> file
#   my $DMRtype = $1;
#   my $DMRfile = $2;

  $file =~ m/\.\.\/bsseq\/DMRs\/$DMRtype\/(.+\_DMRs)\.txt/;
  my $DMRfile = $1;

  # base name for randomDMRs files
  my $randomDMRfile = $DMRfile;
  $randomDMRfile=~ s/DMRs/randomDMRs_/;

  # 100 randomizations
  for (my $i=1; $i <= 100; $i++) {
    print "  Randomization $i\n";
    open(OUT, '>', "./DMRs/$DMRtype/randomized_control_CpG_density/$randomDMRfile".$i.".bed") or die ("Cannot open OUT random DMRs file\n");

    # For each DMR, pick a randomDMR of same length and same number of CpG on same chromosome
    for (my $j=0; $j < scalar(@{$DMRs{$file}}); $j++) {
      my $chr = ${$DMRs{$file}}[$j][0];
      my $DMRlength = ${$DMRs{$file}}[$j][1];
      my $DMRnumberCpG = ${$DMRs{$file}}[$j][2];
      #print "  Generating random DMR on chromosome $chr, of length $DMRlength and with $DMRnumberCpG CpGs...\n";

      # check that there are pre-generated random string of same number of CpGs and length
      unless (exists $randomDMRs{$chr}->{$DMRnumberCpG}){
        # sort by minimal absolute difference between existing keys and actual number of DMRs
        # see http://www.perlmonks.org/?node_id=884064
        my $newDMRnumberCpG = (sort { abs( $a - $DMRnumberCpG ) <=> abs( $b - $DMRnumberCpG ) } keys %{$randomDMRs{$chr}})[0];
        print "    DMR on $chr of length $DMRlength and with $DMRnumberCpG CpGs: no random string with same number of CpGs: changing to closest value: ", $newDMRnumberCpG, "\n";
        $DMRnumberCpG =  $newDMRnumberCpG;
      }
      unless (exists $randomDMRs{$chr}->{$DMRnumberCpG}->{$DMRlength}){
        # sort by minimal absolute difference between existing keys and actual length
        my $newDMRlength = (sort { abs( $a - $DMRlength ) <=> abs( $b - $DMRlength ) } keys %{$randomDMRs{$chr}->{$DMRnumberCpG}})[0];
        print "    DMR on $chr of length $DMRlength and with $DMRnumberCpG CpGs: no random string with same length: changing to closest value: ", $newDMRlength, "\n";
        $DMRlength =  $newDMRlength;
      }

      # pick random string in hash
      # $hash{(keys %hash)[rand keys %hash]};
      my $randomStart = (keys %{$randomDMRs{$chr}->{$DMRnumberCpG}->{$DMRlength}})[rand keys %{$randomDMRs{$chr}->{$DMRnumberCpG}->{$DMRlength}}];
      my $randomEnd = $randomDMRs{$chr}->{$DMRnumberCpG}->{$DMRlength}->{$randomStart};

      print OUT "$chr\t$randomStart\t$randomEnd\n";
    }
    close OUT;
  }
  # last; # testing
}

# TO DO? 
# - it is strange when no fragment with same number of DMRs and similar size in randomDMRs. There should be at least one!
# - When it is the case we take the closest values, but when we do so for number of DMRs, we should decrease the length (sometimes we increase the length if the closest value is bigger)

