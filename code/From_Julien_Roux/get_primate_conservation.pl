#!/usr/bin/env perl -w
use strict;
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::AlignIO;
use Getopt::Long;
use Bio::EnsEMBL::Registry;
use List::Util qw(sum);

if (scalar(@ARGV) ne 1){
  die "USAGE: perl getprimate_conservation.pl DMR_file\n"
}
## read the DMR file
my $DMR_file = $ARGV[0];

# Connect to the Ensembl API
Bio::EnsEMBL::Registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous',
    -port => 5306);

# which Ensembl version is used?
use Bio::EnsEMBL::ApiVersion;
printf( "The API version used is %s. Please check it is appropriate.\n", software_version() );

# Getting the GenomicAlignBlock adaptor
my $genomic_align_block_adaptor = Bio::EnsEMBL::Registry->get_adaptor(
    'Multi', 'compara', 'GenomicAlignBlock');
# Getting the GenomeDB adaptor
my $genome_db_adaptor = Bio::EnsEMBL::Registry->get_adaptor(
    'Multi', 'compara', 'GenomeDB');
# Getting the MethodLinkSpeciesSet adaptor
my $method_link_species_set_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi', 'compara', 'MethodLinkSpeciesSet');

## 8-way primate EPO alignments
my $method_link_species_set = $method_link_species_set_adaptor->fetch_by_method_link_type_species_set_name("EPO", "primates");

# # Function calculating conservation at a region
# # take % conservation of most present nucleotide. Not related to human conservation
# # Inspired by function at https://www.biostars.org/p/9173/
# sub conservation {
#   my ($self) = @_;
#   my @cons;
#   my $num_sequences = $self->no_sequences; #num_sequences in more recent Bioperl versions
#   my $num_species = 0;
#   foreach my $seq ( $self->each_seq() ) {
#     # print $seq->display_id(), "\n";
#     $num_species++ if ($seq->display_id() =~ m/homo_sapiens/);
#     $num_species++ if ($seq->display_id() =~ m/pan_troglodytes/);
#     $num_species++ if ($seq->display_id() =~ m/macaca_mulatta/);
#   }
#   # print "Num species: $num_species\n";
#   if ($num_species < 3){ ## we want at least human, chimp and macaque
#     return 'NA';
#   }
#   else {
#     foreach my $point (0..$self->length-1) { 
#       my %hash;
#       foreach my $seq ( $self->each_seq() ) {
#         my $letter = substr($seq->seq, $point, 1);
#         ($letter eq $self->gap_char || $letter =~ /\./) && next; ## if it is a gap: next
#         $hash{$letter}++;       ## otherwise: record letter
#       }
#       # max frequency of a non-gap letter
#       my $max = (sort {$b<=>$a} values %hash )[0];
#       push @cons, 100 * $max / $num_sequences;
#     }
#     sub mean {
#       return sum(@_)/@_;
#     }
#     return mean(@cons);
#   }
# }

# Function calculating conservation of human sequence at a region
# take % conservation of most present nucleotide. Not related to human conservation
# Inspired by function at https://www.biostars.org/p/9173/
sub humanConservation {
  my ($self) = @_;
  my @cons;
  my $num_sequences = $self->no_sequences; #num_sequences in more recent Bioperl versions
  my $num_species = 0;
  foreach my $seq ( $self->each_seq() ) {
    # print $seq->display_id(), "\n";
    $num_species++ if ($seq->display_id() =~ m/homo_sapiens/);
    $num_species++ if ($seq->display_id() =~ m/pan_troglodytes/);
    $num_species++ if ($seq->display_id() =~ m/macaca_mulatta/);
  }
  # print "Num species: $num_species\n";
  if ($num_species < 3){ # we want at least human, chimp and macaque
    return 'NA';
  }
  else {
    foreach my $point (0..$self->length-1) {
      my %hash;
      my $ref_letter; 
      foreach my $seq ( $self->each_seq() ) {
        my $letter = substr($seq->seq, $point, 1);
        if ($seq->display_id() =~ m/homo_sapiens/){
          $ref_letter = $letter;
        }
        if ($letter eq $self->gap_char || $letter =~ /\./){
          next; # gaps
        }
        $hash{$letter}++; # record all letters (and their counts) at this position
      }
      # frequency of human letter
      if (defined $hash{$ref_letter}){
        push @cons, 100 * $hash{$ref_letter} / $num_sequences;
      } else {
        push @cons, 0; ## if human has a gap 
      } 

    }
    sub mean {
      return sum(@_)/@_;
    }
    return mean(@cons);
  }
}

my $query_species = 'human';
my $outfile = $ARGV[0];
# $outfile =~ s/\.txt/\_meanSeqConservation\.txt/;
$outfile =~ s/\.txt/\_humanSeqConservation\.txt/;
open(OUT, '>', './'.$outfile) or die ("Cannot open OUT file\n");
select( (select(OUT), $| = 1)[0] ); # disables buffering to output file
print OUT "chr\tstart\tend\tmeanGERPscore\thumanPrimateConservation\n"; # header

## read each DMR and calculate conservation
open(IN, $DMR_file) or die "Can't read file";
my $header = <IN>;
while (defined (my $line = <IN>)) {
  chomp $line;
  my @tmp = split(/\t/, $line);
  next unless ($tmp[0] =~ m/chr/);

  my $seq_region = $tmp[0];
  my $seq_region_start = $tmp[1];
  my $seq_region_end = $tmp[2];
  print "Analyzing $seq_region: $seq_region_start-$seq_region_end\n";
  print OUT $seq_region, "\t", $seq_region_start, "\t", $seq_region_end;

  # Getting the Slice adaptor:
  my $slice_adaptor = Bio::EnsEMBL::Registry->get_adaptor($query_species, 'core', 'Slice');
  # Fetching a Slice object:
  my $seq_region_ensembl = $seq_region;
  $seq_region_ensembl =~ s/chr//;
  my $query_slice = $slice_adaptor->fetch_by_region('toplevel', $seq_region_ensembl, $seq_region_start, $seq_region_end);

  ################################################################
  ## Get the GERP score for this slice (mammalian conservation) ##
  ################################################################
  # See https://github.com/Ensembl/ensembl-compara/blob/release/80/scripts/examples/dna_getConservationScores.pl

  #get method_link_species_set object for gerp conservation scores for mammals
  my $mlss = $method_link_species_set_adaptor->fetch_by_method_link_type_species_set_name("GERP_CONSERVATION_SCORE", "mammals");

  my $cs_adaptor = Bio::EnsEMBL::Registry->get_adaptor("Multi", 'compara', 'ConservationScore');
  #To get one score per base in the slice, must set display_size to the size of the slice.
  my $display_size = $query_slice->end - $query_slice->start + 1;
  my $scores = $cs_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $query_slice, $display_size);

  #Meaningful: difference between observed and expected scores.
  my $num_scores = 0;
  my $sum_scores;
  foreach my $score (@$scores) {
    if (defined $score->diff_score) {
      $num_scores++;
      $sum_scores += $score->diff_score;
      # print $score->diff_score, " \n";
    }
  }
  my $mean_score;
  if ($num_scores ne 0){
    $mean_score = $sum_scores / $num_scores;
  }
  if (defined $mean_score) {
    print OUT "\t", $mean_score;
  }
  else {
    print OUT "\tNA";
  }

  ################################################################################
  ## Get the sequence conservation in primates (compared to ref human sequence) ##
  ################################################################################
  # Fetching all the GenomicAlignBlock corresponding to this Slice from the primate alignment
  my $genomic_align_blocks = $genomic_align_block_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($method_link_species_set, $query_slice);

  # We will then (usually) need to restrict the blocks to the required positions in the reference sequence 
  # ($seq_region_start and $seq_region_end), and print all the genomic alignments using a Bio::AlignIO object
  my $alignIO = Bio::AlignIO->newFh(
    -interleaved => 0,
    -fh => \*STDOUT,
    -format => "clustalw",      ## [clustalw|fasta|...]
    -idlength => 10 );

  my $mean_conservation;
  foreach my $genomic_align_block ( @{ $genomic_align_blocks }) {
    my $restricted_gab = $genomic_align_block->restrict_between_reference_positions($seq_region_start, $seq_region_end);

    my $simple_align = $restricted_gab->get_SimpleAlign;
    # print $alignIO $simple_align;

    ## get mean conservation score for this DMR
    ## retunrs NA if on species in human/chimp/macaque missing
    # $mean_conservation = conservation($simple_align);
    $mean_conservation = humanConservation($simple_align);
  }
  if (defined $mean_conservation) {
    print OUT "\t", $mean_conservation, "\n";
  }
  else {
    print OUT "\tNA\n";
  }

  #########################################
  ## Get the mean age of the human bases ##
  #########################################
  ## TO DO!

}
close IN;
close OUT;
exit;
