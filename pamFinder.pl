#!/usr/bin/perl
=head
Author: Mike Tung
Pam Finder 1.00
http://www.github.com/seekheart/PamFinder/
=cut

#load packages/libraries
use warnings;
use strict;
use feature qw(say);
use Bio::SeqIO;
use Getopt::Long;
use Pod::Usage;

#Make command line options
my $fastaFile;
my $strain;
my $guideFile;
my $outfile;
my $usage = "\nPam Finder\n
Options     Description
-fasta      Fasta File to search
-strain     Specify Strain(s) Cas9 to Use
-guide      Guide RNA File
-outfile    Output name
-help       Show this message
";

GetOptions(
	"fasta=s"       =>      \$fastaFile,
	"strain=s"        =>      \$strain,
	"guide=s"      =>      \$guideFile,
	"outfile=s"	=>	\$outfile,
	"help"              =>      sub{pod2usage($usage);},
	  ) or die "$usage";

#check for required CL options
unless($fastaFile and $strain and $guideFile and $outfile){
	die "Error!!\n$usage"
}

#Subroutines Here

#Subroutine that takes a FASTA file as an argument and returns
#the sequence as a scalar.
sub processFasta{
	my ($fasta) = @_;
	my @sequence = ();

#make bioseqIO object to house fasta then get_sequence()
	my $seqIn = Bio::SeqIO->new(
					-file   => 	$fasta,
					-format => 	'Fasta',
					);

#process the seq object
	while (my $seq = $seqIn->next_seq()){
		@sequence = split("", $seq->seq()) ;
	}
	my $seq = join('', @sequence);
	return $seq;
}

#Reverse Complement Sub that takes as input a scalar and outputs an array
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
	open(my $fh, ">>", $outfile) or die "Can't open $outfile! $!";
	say $fh join("\t", @line);
	close $fh;
	say "Writing to $outfile!";
}

#Main Subroutine executes entire script.
sub main{
	#Setup hashes
	#Hash of Cas9 Variants and their respective PAM Sites
	#In Future Versions more Variants will be added.
	my %cas9 = ("SP"	=> 	[
					"AGG", "GGG","CGG", "TGG",
					"AAG", "CAG", "TAG", "GAG",
					"AGCG", "CGCG", "TGCG", "GGCG",
					"AGAG", "CGAG", "TGAG", "GGAG",
					"AGAA", "AGAG", "AGAC", "AGAT",
					"GGAA", "GGAG", "GGAC", "GGAT",
					"CGAA", "CGAG", "CGAC", "CGAT",
					"TGAA", "TGAG", "TGAC", "TGAT",
					"AGAG", "AGGG", "AGCG", "AGTG",
					"GGAG", "GGGG", "GGCG", "GGTG",
					"CGAG", "CGGG", "CGCG", "CGTG",
					"TGAG", "TGGG", "TGCG", "TGTG"
					],
		   );
	#Process the sGRNA file in to array of guides (already rev complemented)
	my @guides;
	my $guide;
	open(my $fh, "<", $guideFile) or die "Couldn't Open File!";
	while (<$fh>) {
		chomp;
		$guide = join('', $_);
		$guide = reverseComplement($guide);
		push @guides, $guide;
	}
	close $fh;


	#Load the FASTA File's Sequence
	my $referenceSeq = processFasta($fastaFile);
	my $reverseRefSeq = reverseComplement($referenceSeq);

	#For each of the guides attach PAM site and check ref seq for hit
	#can find location using @- @+ special vars in perl
	my @line;
	my $count = 0;
	my $targetSite;
	foreach $guide (@guides){
		foreach my $pam (@{$cas9{$strain}}){
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
