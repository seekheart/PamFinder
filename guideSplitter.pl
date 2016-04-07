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
my $usage = "\n$0\n

Usage:
-infile           Input file of guide strand that is longer than 20 characters
-outfile          Output file of guide converted to multiple 14-20 character wide lines
-verbose          Display output, off by default
-help             Display This

";

GetOptions(
    "infile=s"    =>    \$infile,
    "outfile=s"   =>    \$outfile,
    help          =>    sub{pod2usage($usage);},
);

# check required options are met
unless($infile && $outfile){
    die "$usage";
}

# Process out the full guide sequence
open(my $fileHandle, '<', $infile) or die "Invalid File!";
my $guide = "";
while (<$fileHandle>){
    chomp;
    if ($_ =~ m/^([^>].*)/){
        $guide = $guide . $1;
    }
}
close ($fileHandle);

# Figure out what is the optimal size of each line (must be minimum 14 and maximum 20)
my $guideSize = length($guide);
my $guideStep = 20;

my $option2 = 0;
while ($guideSize % $guideStep != 0 && $guideSize % $guideStep < 14){
    $guideStep = $guideStep - 1;
#    print "New step size = $guideStep\n";
#    print $guideSize % $guideStep, "\n";
    if ($guideStep == 13)
    {
        $guideStep = 20;
        $option2 = 1;
        last;
    }
}

for (my $count = 14; $count < 100; $count++){
    $guideStep = 20;
    $option2 = 0;
    while($count % $guideStep != 0 && $count % $guideStep < 14){
        $guideStep = $guideStep - 1;
        if ($guideStep == 13)
        {
            $option2 = 1;
            last;
        }
    }
    if ($option2 == 0)
    {
        print "$count\t$guideStep\n";
    }
    else
    {
        print "******\t$count is not easily convertible\t******\n"
    }
}
