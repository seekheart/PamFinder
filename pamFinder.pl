#!/usr/bin/perl
=head
Author: Mike Tung, Meryl STav
Pam Finder 1.00
http://www.github.com/seekheart/PamFinder/
=cut

# Load packages/libraries
use warnings;
use strict;
use feature qw(say);
use Bio::SeqIO;
use Getopt::Long;
use Pod::Usage;

# Make command line options
my $fastaFile;
my $strain;
my $guideFile;
my $outfile;
my $usage = "\nPam Finder\n
Options     Description
-fasta      Fasta File to search
-strain     Specify Strain(s) cas9 to Use
-guide      Guide RNA File
-outfile    Output name
-help       Show this message
";

GetOptions(
  "fasta=s"     =>    \$fastaFile,
  "strain=s"    =>    \$strain,
  "guide=s"     =>    \$guideFile,
  "outfile=s"   =>    \$outfile,
  "help"        =>    sub{pod2usage($usage);},
) or die "$usage";

# Check for required CL options
unless($fastaFile and $strain and $guideFile and $outfile){
  die "Error!!\n$usage"
}

################################################################################
# Subroutines Here:
################################################################################

# Subroutine that takes a FASTA file as an argument and returns the sequence as
# a scalar.
sub processFasta{
  my ($fasta) = @_;
  my @sequence = ();

  # Make bioseqIO object to house fasta then get_sequence()
  my $seqIn = Bio::SeqIO->new(
          -file   =>   $fasta,
          -format =>   'Fasta',
          );

  # Process the seq object
  while (my $seq = $seqIn->next_seq()){
    @sequence = split("", $seq->seq()) ;
  }
  my $seq = join('', @sequence);
  return $seq;
}

# Reverse Complement Sub that takes as input a scalar and outputs an array
sub reverseComplement{
  my ($seq) = @_;
  my @sequence = split("", $seq);
  my @newSeq = ();
  foreach my $nucleotide (@sequence){
    chomp $nucleotide;
    if ($nucleotide eq "A"){
      push @newSeq, "T";
    }
    elsif($nucleotide eq "T"){
      push @newSeq, "A";
    }
    elsif($nucleotide eq "C"){
      push @newSeq, "G";
    }
    else{
      push @newSeq, "C";
    }
  }
  my $newSeq = reverse(@newSeq);
  return $newSeq;
}

sub outputFile{
  my @line = @_;
  open(my $fileHandle, ">>", $outfile) or die "Can't open $outfile! $!";
  say $fileHandle join("\t", @line);
  close $fileHandle;
  say "Writing to $outfile!";
}

# Main Subroutine executes entire script.
sub main{
  # Setup hashes
  # Hash of cas9 Variants and their respective PAM Sites
  # In Future Versions more Variants will be added.
  my %cas9Variants = (
  	# SP => NGG, NAG, NGCG, NGAN, NGGG, NGTG 
  	"SP"  =>   [
          "AGG", "GGG","CGG", "TGG",
          "AAG", "CAG", "TAG", "GAG",
          "AGCG", "CGCG", "TGCG", "GGCG",
          "AGAG", "CGAG", "TGAG", "GGAG",
          "AGAA", "CGAA", "TGAA", "GGAA",
          "AGAC", "CGAC", "TGAC", "GGAC",
          "AGAT", "CGAT", "TGAT", "GGAT",
          "AGGG", "CGGG", "TGGG", "GGGG",
          "AGTG", "CGTG", "TGTG", "GGTG"
          ],
  	# SpCas9 => NGG, NAG, NGCG, NGAN, NGGG, NGTG
    "SpCas9"  =>   [
          "AGG", "GGG","CGG", "TGG",
          "AAG", "CAG", "TAG", "GAG",
          "AGCG", "CGCG", "TGCG", "GGCG",
          "AGAG", "CGAG", "TGAG", "GGAG",
          "AGAA", "CGAA", "TGAA", "GGAA",
          "AGAC", "CGAC", "TGAC", "GGAC",
          "AGAT", "CGAT", "TGAT", "GGAT",
          "AGGG", "CGGG", "TGGG", "GGGG",
          "AGTG", "CGTG", "TGTG", "GGTG"
          ],
  	# SpCas9_D1135E => NGG (reduced binding), NAG, NGCG, NGAN, NGGG, NGTG
    "SpCas9_D1135E"  =>   [
          "AGG", "GGG","CGG", "TGG",
          "AAG", "CAG", "TAG", "GAG",
          "AGCG", "CGCG", "TGCG", "GGCG",
          "AGAG", "CGAG", "TGAG", "GGAG",
          "AGAA", "CGAA", "TGAA", "GGAA",
          "AGAC", "CGAC", "TGAC", "GGAC",
          "AGAT", "CGAT", "TGAT", "GGAT",
          "AGGG", "CGGG", "TGGG", "GGGG",
          "AGTG", "CGTG", "TGTG", "GGTG"
          ],
  	# SpCas9_VRER => NGCG, NAG, NGAN, NGGG, NGTG
    "SpCas9_VRER"  =>   [
          "AGCG", "CGCG", "TGCG", "GGCG",
          "AAG", "CAG", "TAG", "GAG",
          "AGAG", "CGAG", "TGAG", "GGAG",
          "AGAA", "CGAA", "TGAA", "GGAA",
          "AGAC", "CGAC", "TGAC", "GGAC",
          "AGAT", "CGAT", "TGAT", "GGAT",
          "AGGG", "CGGG", "TGGG", "GGGG",
          "AGTG", "CGTG", "TGTG", "GGTG"
          ],
  	# SpCas9_EQR => NGAG, NAG, NGCG, NGAN, NGGG, NGTG
    "SpCas9_EQR"  =>   [
          "AGAG", "CGAG", "TGAG", "GGAG",
          "AAG", "CAG", "TAG", "GAG",
          "AGCG", "CGCG", "TGCG", "GGCG",
          "AGAA", "CGAA", "TGAA", "GGAA",
          "AGAC", "CGAC", "TGAC", "GGAC",
          "AGAT", "CGAT", "TGAT", "GGAT",
          "AGGG", "CGGG", "TGGG", "GGGG",
          "AGTG", "CGTG", "TGTG", "GGTG"
          ],
  	# SpCas9_VQR => NGAN, NGNG, NAG, (repeated NGCG, NGGG, NGTG)
    "SpCas9_VQR"  =>   [
          "AGAA", "CGAA", "TGAA", "GGAA",
          "AGAC", "CGAC", "TGAC", "GGAC",
          "AGAT", "CGAT", "TGAT", "GGAT",
          "AGAG", "CGAG", "TGAG", "GGAG",
          "AGCG", "CGCG", "TGCG", "GGCG",
          "AGTG", "CGTG", "TGTG", "GGTG",
          "AGGG", "CGGG", "TGGG", "GGGG",
          "AAG", "CAG", "TAG", "GAG",
          ],
       );
  # Process the sGRNA file in to array of guides (already rev complemented)
  my @guides;
  my $guide;
  open(my $fileHandle, "<", $guideFile) or die "Couldn't Open File!";
  while (<$fileHandle>) {
    chomp;
    $guide = join('', $_);
    $guide = reverseComplement($guide);
    push @guides, $guide;
  }
  close $fileHandle;

  # Load the FASTA File's Sequence
  my $referenceSeq = processFasta($fastaFile);
  my $reverseRefSeq = reverseComplement($referenceSeq);

  # For each of the guides attach PAM site and check ref seq for hit can find location
  # using @- @+ special vars in perl
  my @line;
  my $count = 0;
  my $targetSite;
  foreach $guide (@guides){
    foreach my $pam (@{$cas9Variants{$strain}}){
      $targetSite = $guide . $pam;
      if ($count == 0){
      #if the filename already exists delete it.
        if (-e $outfile){
          system("rm $outfile");
        }
        say "Making Header";
        @line = ("GUIDE RNA", "START_POS", "END_POS", "REF_FILE", "STRAND");
        outputFile(@line);
      }
      if ($referenceSeq =~ m/$targetSite/i){
        say "Guide: $guide Matched at position [@-, @+] in Sense strand!";
        @line = ($guide, @-, @+, $fastaFile, "+");
        outputFile(@line);
      }
      if ($reverseRefSeq =~ m/$targetSite/i){
        say "Guide: $guide Matched at position [@-, @+] in Anti-Sense Strand!";
        @line = ($guide, @-, @+, $fastaFile, "-");
        outputFile(@line);
      }
      $count++;
    }
  }
  say "Done!";


}

#Run program here through main sub
say "";
main();
