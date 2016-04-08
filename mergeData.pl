#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

# Make the options
my $dir;
my $prefix;
my $outfile;
my $usage = "\n$0\n

Usage:
-dir            Input directory containing all files to be merged
-prefix         Prefix of all files in directory to be merged
-outfile        Output merged file
-verbose        Display output, off by default
-help           Display This

";

GetOptions(
    "dir=s"         =>      \$dir,
    "prefix=s"      =>      \$prefix,
    "outfile=s"     =>      \$outfile,
    help            =>      sub{pod2usage($usage);},
);

# check required options are met
unless($dir && $prefix && $outfile){
    die "$usage";
}

# Extract out the files that need to be merged
opendir (my $directoryHandle, $dir) or die $!;
my @filesToMerge = ();
while (my $file = readdir($directoryHandle)) {
    if ($file =~ m/(^\Q$prefix\E.*)/) {
        push @filesToMerge, $dir . "/" . $file;
    }
}
closedir($directoryHandle);

# Extract information from files and put into merged file
open(my $outFileHandle, '>', $outfile) or die "Invalid File!";

print $outFileHandle "GUIDE_RNA\tSTART_POS\tEND_POS\tREF_FILE\tSTRAND\tOUT_FILE\n";

foreach my $infile (@filesToMerge){
    open(my $inFileHandle, '<', $infile) or die "Invalid File!";
    while (<$inFileHandle>){
        chomp;
        if ($_ !~ m/^GUIDE.*/) {
            print $outFileHandle $_, "\t", $infile, "\n";
        }
    }
    close($inFileHandle);
}
close($outFileHandle);