#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::DB::EUtilities;
use Bio::SeqIO;
use feature qw(say);

# Make the options
my $infile;
my $outfile;
my $email;
my $location;
my $usage = "\n$0\n

Usage:
    -infile         Input File with Genbank/RefSeq Ids
    -email          Your Email
    -outfile        Output Base File Name
    -location       Folder name for all outputs
    -verbose        Display output, off by default
    -help           Display This

";

GetOptions(
    "infile=s"      =>      \$infile,
    "email=s"       =>      \$email,
    "outfile=s"     =>      \$outfile,
    "location=s"    =>      \$location,
    help            =>      sub{pod2usage($usage);},
);

# check required options are met
unless($infile && $email && $outfile){
    die "$usage";
}

# Make output directory
my $currentDir = `pwd`;
chomp $currentDir;
my $dirName = $currentDir . "/" . $location;
unless (-e $dirName){
    `mkdir -p $dirName`;
}

# Process and extract out the RefSEq IDs
open(my $fh, '<', $infile) or die "Invalid File!";
my @ids = ();
while (<$fh>){
    chomp;
    if ($_ =~ m/(NC.{9})/) {
        push @ids, $1;
    }
}
close($fh);

# Make a factory to fetch fastas
foreach my $id (@ids){
    my $factory = Bio::DB::EUtilities->new(
        -eutil      =>      'efetch',
        -rettype    =>      'fasta',
        -db         =>      'nucleotide',
        -id         =>      $id,
        -email      =>      $email,
        -tool       =>      'get_refseq_id'
    );
    # Output file to directory
    my $file = $dirName . '/' . $outfile . '.' . $id . '.fasta';
    if (-s $file){
        say STDERR "File Exists Skipping...";
    }
    else{
        $factory->get_Response( -file => $file);
        say STDERR "Fetching data from NCBI for $id";
        sleep(3);
    }
}