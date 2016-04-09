#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use POSIX qw(strftime);

# Get timestamp
my $timestamp = strftime "%F_%H-%M-%S", localtime;

##########################################################################################
# User, please fill in information below. Default is currently entered.
##########################################################################################

my $infile          =   "yeastDataID";          # file of information to parse from NCBI
my $email           =   'merylstav@gmail.com';  # your email to use for searching NCHBI
my $prefix          =   "TEST";                 # prefix for all intermediate files
my $intermediates   =   "output";               # location of where to put intermediates

my $strain          =   "SP";                   # which cas9 do you want to use
my $guideFile       =   "S288C_ARS103_ARS103_genomic.fsa";  # name of your guide RNA file

my $finalOutput     =   $timestamp . "_yeastResults.txt";

##########################################################################################
# User, please fill in information above. Default is currently entered.
##########################################################################################
print "Step 1 Start: Downloaded all necessary data from NCBI database\n";
# Download files from NCBI database
my @ARGS_download = (
    '-infile    yeastDataID',
    "-email     $email",
    "-outfile   $prefix",
    "-location  $intermediates"
);
system("perl downloadData.pl @ARGS_download");

print "Step 1 Complete: All data downloaded from NCBI database\n";
print "Step 2 Start: Readjust shape of guide RNA file if necessary\n";

my $dir             =   $intermediates;
my $guideAdjusted   =   $dir . "/" . "guide_adjusted.txt";

# Make sure guideRNA file is set so that each line is between 14-20 characters only
my @ARGS_guideSplitter = (
    "-infile    $guideFile",
    "-outfile   $guideAdjusted"
);
system ("perl guideSplitter.pl @ARGS_guideSplitter");

print "Step 2 Complete: Guide RNA file has been readjusted to lines of length 14-20\n";
print "Step 3 Start: Parse each fasta file for possible PAM sites\n";

# Create list of all fasta files downloaded
opendir (my $directoryHandle, $dir) or die $!;
my @fastaFiles = ();
while (my $file = readdir($directoryHandle)) {
    if ($file =~ m/(^\Q$prefix\E.*)/) {
        push @fastaFiles, $file;
    }
}
closedir($directoryHandle);

# Run pamFinder on each fasta file
foreach my $fastaFile (@fastaFiles){
    my $PAMresults = $dir . "/" . $timestamp . "_output_" . $fastaFile;
    my @ARGS_pamFinder = (
        "-fasta     $dir/$fastaFile",
        "-strain    $strain",
        "-guide     $guideAdjusted",
        "-outfile   $PAMresults"
    );
    system("perl pamFinder.pl @ARGS_pamFinder");
}
print "Step 3 Complete: All fasta files have been parsed for possible PAM sites\n";
print "Step 4 Start: Merge all intermediate output files into one file\n";

# Merge all results into 1 handy to read file
my @ARGS_merge = (
    "-dir       $dir",
    "-prefix    $timestamp",
    "-outfile   $finalOutput"
);
system("perl mergeData.pl @ARGS_merge");
print "Step 4 Complete: All files merged into $finalOutput in current directory\n";

print "Automator complete\n";