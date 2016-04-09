#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

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
# Option 1: All lines (but last one) are size Step (and last is between 14 and Step)
# Option 2: All line but last are size Step and last is between Step and 20
# Option 3: No easy conversion possible, so cut off minimum possible
my $guideLength = length($guide);
my $stepSize = 20;
my $option1 = 0;
my $option2 = 0;
my $option3 = 1;
my $remainder = -1;

# In the case that the guide is already within the correct boundaries, there is no need to
# readjust; only to copy the string to a new file.
if ($guideLength <= 20 && $guideLength >= 14)
{
    $stepSize = $guideLength;
    $remainder = 0;
}
elsif ($guideLength < 14)
{
    print "********ERROR: GUIDE STRING IS TOO SHORT************************\n";
    print "********ERROR: PLEASE RE-RUN WITH DIFFERENT GUIDE STRING********\n";
    $stepSize = $guideLength;
    $remainder = $guideLength;
}
else
{
    while($stepSize > 13)
    {
        $remainder = $guideLength % $stepSize;
        if ($remainder == 0 || $remainder > 13)
        {
            $option1 = 1;
            $option3 = 0;
            last;
        }
        elsif ($remainder + $stepSize >=14 && $remainder + $stepSize <= 20)
        {
            $option2 = 1;
            $option3 = 0;
            last;
        }       
        else
        {
            $stepSize--;
        }
    }

    # If string needs to be cut, find the permutation with least to remove
    if ($option3 == 1)
    {
        $stepSize = 20;
        my $minStep = $stepSize;
        my $minRemainder = $guideLength % $stepSize;
        while ($stepSize > 13)
        {
            $stepSize = $stepSize - 1;
            if ($guideLength % $stepSize < $minRemainder)
            {
                $minStep = $stepSize;
                $minRemainder = $guideLength % $stepSize;
            }
        }
        # Readjust the length of the guide rena
        $guide = substr($guide, 0, -1*$minRemainder);
        $guideLength = length($guide);
        $remainder = $minRemainder;
        $option1 = 1;
        $option3 = 0;
    }
}

# Write rearranged guide to new file
open($fileHandle, '>', $outfile) or die "Invalid File!";
my $windowEnd = $guideLength - $remainder;
my $windowStart;
my $line;
for ($windowStart = 0; $windowStart < $windowEnd; $windowStart = $windowStart + $stepSize)
{
    $line = substr($guide, $windowStart, $stepSize);
    print $fileHandle $line, "\n";
}
if ($windowStart != $windowEnd)
{
    $line = substr($guide, $windowStart, $remainder);
    print $fileHandle $line, "\n";
}

close ($fileHandle);

## Unit Test for parsing guide RNA
# for (my $count = 14; $count <= 100; $count++){
#     $stepSize = 20;
#     $remainder = $count % $stepSize;
#     $option1 = 0;
#     $option2 = 0;
#     $option3 = 1;
#     while($stepSize > 13)
#     {
#         $remainder = $count % $stepSize;
#         if ($remainder == 0 || $remainder > 13)
#         {
#             $option1 = 1;
#             $option3 = 0;
#             last;
#         }
#         elsif ($remainder + $stepSize >=14 && $remainder + $stepSize <= 20)
#         {
#             $option2 = 1;
#             $option3 = 0;
#             last;
#         }       
#         else
#         {
#             $stepSize--;
#         }
#     }
#     if ($option1 == 1)
#     {
#         print "1\t$count\t$stepSize\t$remainder\n";
#     }
#     if ($option2 == 1)
#     {
#         print "2\t$count\t$stepSize\t", $remainder+$stepSize,"\n";
#     }
#     if ($option3 == 1)
#     {
#         print "******\t$count is not easily convertible\t******\n";
#         $stepSize = 20;
#         my $minStep = $stepSize;
#         my $minRemainder = $guideLength % $stepSize;
#         while ($stepSize > 13)
#         {
#             $stepSize = $stepSize - 1;
#             if ($guideLength % $stepSize < $minRemainder)
#             {
#                 $minStep = $stepSize;
#                 $minRemainder = $guideLength % $stepSize;
#             }
#         }
#         Readjust the length of the guide rna
#         $guide = repeat($count, "a");
#         $guide = substr($guide, 0, -1*$minRemainder);
#         $guideLength = length($guide);
#         $remainder = $minRemainder;
#         $option1 = 1;
#         $option3 = 0;
#         print "******\t$count is now $guideLength\t$stepSize\t$remainder\t******\n";
#     }
# }