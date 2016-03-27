#!/usr/bin/perl
#@AUTHORS - Mike Tung, Meryl Stav
use warnings;
use strict;
use feature qw(say);
use Bio::SeqIO;
use Getopt::Long;
use Pod::Usage;


#Make Options
my $fastaFile;
my $cas9;
my $guideFile;
my $usage = "\n$0\n
Options     Description
-fasta      Fasta File to search
-cas9       Specify Strain(s) Cas9 to Use
-guide      Guide RNA File
-help       Show this message
";

GetOptions(
                    "fasta=s"       =>      \$fastaFile,
                    "cas9=s"         =>      \$cas9,
                    "guide=s"      =>      \$guideFile,
                    "help"              =>      sub{pod2usage($usage);},
    ) or die "$usage";

unless($fastaFile and $cas9 and $guideFile){
    die "Error!!\n$usage"
}

#Subroutines
sub processFasta{
    my ($fasta) = @_;
    my @sequence = ();

    #make bioseqIO object to house fasta then get_sequence()
    my $seqIn = Bio::SeqIO->new(   -file => $fasta,
                                                    -format => 'Fasta');
    #process the seq object
    while (my $seq = $seqIn->next_seq()){
        @sequence = split("", $seq->seq()) ;
    }
    return @sequence;
}

sub main{
    my @sequence = processFasta($fastaFile);
    my %cas9 = (    "SP"                => ("NGG"),
                             "SP D1135E"    => ("NGG", "NAG"),
                             "SP VRER"        => ("NGCG"),
                             "SP EQR"          => ("NGAG"),
                             "SP VQR"          => ("NGAN", "NGNG"),
                             "SA"                 => ("NNGRRT", "NNGRR", "NNGRRN"),
                             "NM"                 => ("NNNNGATT"),
                             "ST"                   => ("NNAGAAW"),
                             "TD"                   => ("NAAAAC"),
                        );
    my %code = (
                            "A"         =>          "A",
                            "C"         =>          "C",
                            "G"         =>          "G",
                            "T"         =>          "T",
                            "R"         =>          "[AG]",
                            "Y"         =>          "[CT]",
                            "S"         =>          "[CG]",
                            "W"         =>          "[AT]",
                            "K"         =>          "[GT]",
                            "M"         =>          "[AC]",
                            "B"         =>          "[CGT]",
                            "D"         =>          "[AGT]",
                            "H"         =>          "[ACT]",
                            "V"         =>          "[ACG]",
                            "N"         =>          "[AGCT]",
                        );
    say $cas9{$cas9};

}


###TODO###
# [X] Start
# [X] Hashtable of {Cas9 Strain : PAM}
# [X] IUPAC CODE {S: REGEX PATTERN ([GC])}

#run
main();
