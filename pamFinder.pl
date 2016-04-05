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
  	# Streptococcus pyogenes (SP) => NGG, NAG (D1135E var), NGCG (VRER var), NGAN (VQR var
  	# and EQR var (NGAG)), NGNG (VQR var)
  	"SP"  =>   [
          "AGG", "CGG", "TGG", "GGG",
          
          "AAG", "CAG", "TAG", "GAG",
          
          "AGCG", "CGCG", "TGCG", "GGCG",
          
          "AGAA", "CGAA", "TGAA", "GGAA",
          "AGAC", "CGAC", "TGAC", "GGAC",
          "AGAT", "CGAT", "TGAT", "GGAT",
          "AGAG", "CGAG", "TGAG", "GGAG",
          
          "AGTG", "CGTG", "TGTG", "GGTG",
          "AGGG", "CGGG", "TGGG", "GGGG"
          ],
  	# Staphylococcus aureus (SA) => NNGRR, NNGRRN (includes NNGRRT)
  	"SA"  =>   [
  	      "AAGRR", "CAGRR", "TAGRR", "GAGRR",
  	      "ACGRR", "CCGRR", "TCGRR", "GCGRR",
  	      "ATGRR", "CTGRR", "TTGRR", "GTGRR",
  	      "AGGRR", "CGGRR", "TGGRR", "GGGRR",
  	      
  	      "AAGRRA", "CAGRRA", "TAGRRA", "GAGRRA",
  	      "ACGRRA", "CCGRRA", "TCGRRA", "GCGRRA",
  	      "ATGRRA", "CTGRRA", "TTGRRA", "GTGRRA",
  	      "AGGRRA", "CGGRRA", "TGGRRA", "GGGRRA",
  	      
  	      "AAGRRC", "CAGRRC", "TAGRRC", "GAGRRC",
  	      "ACGRRC", "CCGRRC", "TCGRRC", "GCGRRC",
  	      "ATGRRC", "CTGRRC", "TTGRRC", "GTGRRC",
  	      "AGGRRC", "CGGRRC", "TGGRRC", "GGGRRC",
  	      
  	      "AAGRRT", "CAGRRT", "TAGRRT", "GAGRRT",
  	      "ACGRRT", "CCGRRT", "TCGRRT", "GCGRRT",
  	      "ATGRRT", "CTGRRT", "TTGRRT", "GTGRRT",
  	      "AGGRRT", "CGGRRT", "TGGRRT", "GGGRRT",
  	      
  	      "AAGRRG", "CAGRRG", "TAGRRG", "GAGRRG",
  	      "ACGRRG", "CCGRRG", "TCGRRG", "GCGRRG",
  	      "ATGRRG", "CTGRRG", "TTGRRG", "GTGRRG",
  	      "AGGRRG", "CGGRRG", "TGGRRG", "GGGRRG"
          ],
  	# Neisseria meningitidis (NM) => NNNNGATT
    "NM"  =>   [
          "AAAAGATT", "CAAAGATT", "TAAAGATT", "GAAAGATT",
          "ACAAGATT", "CCAAGATT", "TCAAGATT", "GCAAGATT",
          "ATAAGATT", "CTAAGATT", "TTAAGATT", "GTAAGATT",
          "AGAAGATT", "CGAAGATT", "TGAAGATT", "GGAAGATT",
          
          "AACAGATT", "CACAGATT", "TACAGATT", "GACAGATT",
          "ACCAGATT", "CCCAGATT", "TCCAGATT", "GCCAGATT",
          "ATCAGATT", "CTCAGATT", "TTCAGATT", "GTCAGATT",
          "AGCAGATT", "CGCAGATT", "TGCAGATT", "GGCAGATT",
          
          "AATAGATT", "CATAGATT", "TATAGATT", "GATAGATT",
          "ACTAGATT", "CCTAGATT", "TCTAGATT", "GCTAGATT",
          "ATTAGATT", "CTTAGATT", "TTTAGATT", "GTTAGATT",
          "AGTAGATT", "CGTAGATT", "TGTAGATT", "GGTAGATT",
          
          "AAGAGATT", "CAGAGATT", "TAGAGATT", "GAGAGATT",
          "ACGAGATT", "CCGAGATT", "TCGAGATT", "GCGAGATT",
          "ATGAGATT", "CTGAGATT", "TTGAGATT", "GTGAGATT",
          "AGGAGATT", "CGGAGATT", "TGGAGATT", "GGGAGATT",
          
          
          "AAACGATT", "CAACGATT", "TAACGATT", "GAACGATT",
          "ACACGATT", "CCACGATT", "TCACGATT", "GCACGATT",
          "ATACGATT", "CTACGATT", "TTACGATT", "GTACGATT",
          "AGACGATT", "CGACGATT", "TGACGATT", "GGACGATT",
          
          "AACCGATT", "CACCGATT", "TACCGATT", "GACCGATT",
          "ACCCGATT", "CCCCGATT", "TCCCGATT", "GCCCGATT",
          "ATCCGATT", "CTCCGATT", "TTCCGATT", "GTCCGATT",
          "AGCCGATT", "CGCCGATT", "TGCCGATT", "GGCCGATT",
          
          "AATCGATT", "CATCGATT", "TATCGATT", "GATCGATT",
          "ACTCGATT", "CCTCGATT", "TCTCGATT", "GCTCGATT",
          "ATTCGATT", "CTTCGATT", "TTTCGATT", "GTTCGATT",
          "AGTCGATT", "CGTCGATT", "TGTCGATT", "GGTCGATT",
          
          "AAGCGATT", "CAGCGATT", "TAGCGATT", "GAGCGATT",
          "ACGCGATT", "CCGCGATT", "TCGCGATT", "GCGCGATT",
          "ATGCGATT", "CTGCGATT", "TTGCGATT", "GTGCGATT",
          "AGGCGATT", "CGGCGATT", "TGGCGATT", "GGGCGATT",
          
          
          "AAATGATT", "CAATGATT", "TAATGATT", "GAATGATT",
          "ACATGATT", "CCATGATT", "TCATGATT", "GCATGATT",
          "ATATGATT", "CTATGATT", "TTATGATT", "GTATGATT",
          "AGATGATT", "CGATGATT", "TGATGATT", "GGATGATT",
          
          "AACTGATT", "CACTGATT", "TACTGATT", "GACTGATT",
          "ACCTGATT", "CCCTGATT", "TCCTGATT", "GCCTGATT",
          "ATCTGATT", "CTCTGATT", "TTCTGATT", "GTCTGATT",
          "AGCTGATT", "CGCTGATT", "TGCTGATT", "GGCTGATT",
          
          "AATTGATT", "CATTGATT", "TATTGATT", "GATTGATT",
          "ACTTGATT", "CCTTGATT", "TCTTGATT", "GCTTGATT",
          "ATTTGATT", "CTTTGATT", "TTTTGATT", "GTTTGATT",
          "AGTTGATT", "CGTTGATT", "TGTTGATT", "GGTTGATT",
          
          "AAGTGATT", "CAGTGATT", "TAGTGATT", "GAGTGATT",
          "ACGTGATT", "CCGTGATT", "TCGTGATT", "GCGTGATT",
          "ATGTGATT", "CTGTGATT", "TTGTGATT", "GTGTGATT",
          "AGGTGATT", "CGGTGATT", "TGGTGATT", "GGGTGATT",
          
          
          "AAAGGATT", "CAAGGATT", "TAAGGATT", "GAAGGATT",
          "ACAGGATT", "CCAGGATT", "TCAGGATT", "GCAGGATT",
          "ATAGGATT", "CTAGGATT", "TTAGGATT", "GTAGGATT",
          "AGAGGATT", "CGAGGATT", "TGAGGATT", "GGAGGATT",
          
          "AACGGATT", "CACGGATT", "TACGGATT", "GACGGATT",
          "ACCGGATT", "CCCGGATT", "TCCGGATT", "GCCGGATT",
          "ATCGGATT", "CTCGGATT", "TTCGGATT", "GTCGGATT",
          "AGCGGATT", "CGCGGATT", "TGCGGATT", "GGCGGATT",
          
          "AATGGATT", "CATGGATT", "TATGGATT", "GATGGATT",
          "ACTGGATT", "CCTGGATT", "TCTGGATT", "GCTGGATT",
          "ATTGGATT", "CTTGGATT", "TTTGGATT", "GTTGGATT",
          "AGTGGATT", "CGTGGATT", "TGTGGATT", "GGTGGATT",
          
          "AAGGGATT", "CAGGGATT", "TAGGGATT", "GAGGGATT",
          "ACGGGATT", "CCGGGATT", "TCGGGATT", "GCGGGATT",
          "ATGGGATT", "CTGGGATT", "TTGGGATT", "GTGGGATT",
          "AGGGGATT", "CGGGGATT", "TGGGGATT", "GGGGGATT"
          ],
  	# Streptococcus thermophilus (ST)	=> NNAGAAW
    "ST"  =>   [
          "AAAGAAA", "CAAGAAA", "TAAGAAA", "GAAGAAA",
          "ACAGAAA", "CCAGAAA", "TCAGAAA", "GCAGAAA",
          "ATAGAAA", "CTAGAAA", "TTAGAAA", "GTAGAAA",
          "AGAGAAA", "CGAGAAA", "TGAGAAA", "GGAGAAA",
          
          "AAAGAAT", "CAAGAAT", "TAAGAAT", "GAAGAAT",
          "ACAGAAT", "CCAGAAT", "TCAGAAT", "GCAGAAT",
          "ATAGAAT", "CTAGAAT", "TTAGAAT", "GTAGAAT",
          "AGAGAAT", "CGAGAAT", "TGAGAAT", "GGAGAAT"
          ],
  	# Treponema denticola (TD)	=> NAAAAC
    "TD"  =>   [
          "AAAAAC", "CAAAAC", "TAAAAC", "GAAAAC"
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
