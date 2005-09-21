#!/usr/local/bin/perl
# $Id: asm2TampaResults.pl,v 1.2 2005-09-21 20:12:30 catmandew Exp $
#
# Wrapper to run and post-process results from TAMPA
# (Tool for Analyzing Mate Pairs in Assemblies)
#
# Ian M. Dew, Brian Walenz, Granger Sutton. "A Tool for Analyzing Mate
#   Pairs in Assemblies (TAMPA)", J Comp Biol 2005 12(5):497-513.
#
# Written by Ian Dew
#

use Carp;
# use strict;
use FileHandle;
use Getopt::Long;

my $MY_VERSION = " Version 1.01 (Build " . (qw/$Revision: 1.2 $/ )[1]. ")";

my $HELPTEXT = qq~
Produce TAMPA results from an assembly

    asm2TampaResults  [options]  <-a assembly>

    <-a assembly>  The prefix of the input filenames such that
                   assembly.gkpStore = gatekeeper store
                   assembly.frgStore = fragment store
                   assembly.asm = assembly file
  
    options:
      -h               Print help.

$MY_VERSION

~;


my $assemblyPrefix = "";
my $helpRequested;
GetOptions("a=s" => \$assemblyPrefix,
           "h|help" => \$helpRequested);

if($helpRequested)
{
  print STDERR "Help requested:\n\n";
  print STDERR $HELPTEXT;
  exit 0;
}

if(!$assemblyPrefix)
{
  print STDERR "Please specify an assembly\n\n";
  print STDERR $HELPTEXT;
  exit 1;
}

# build an assembly store
my $gkpStorename = $assemblyPrefix . ".gkpStore";
my $frgStorename = $assemblyPrefix . ".frgStore";
my $asmFilename = $assemblyPrefix . ".asm";
my $asmStorename = $assemblyPrefix . ".asmStore";
my $command = "asm2asmStore -g $gkpStorename -f $frgStorename -a $asmFilename -s $asmStorename";
print "Running $command\n";
system($command) == 0
  or die "\nFailed to run command\n$command\n\n";

# create a tampa subdirectory
my $subdir = "tampa";
print "Creating & cd'ing to $subdir subdirectory\n";
mkdir($subdir) or die "Failed to create $subdir dir";
chdir($subdir) or die "Failed to chdir to $subdir";

# populate tampa subdirectory with files
$command = "dumpForTampa -a $assemblyPrefix -s ../$asmStorename";
print "Running $command\n";
system($command) == 0
  or die "\nFailed to run command\n$command\n\n";

# run TAMPA
my $libFilename = $assemblyPrefix . "Libs.txt";
$command = "runTampa.pl -a $assemblyPrefix -l $libFilename";
print "Running $command\n";
system($command) == 0
  or die "\nFailed to run command\n$command\n\n";

# copy spreadsheet file
my $spreadsheetFilename = $assemblyPrefix . ".intra.spreadsheet.txt";
$command = "cp $spreadsheetFilename ..";
print "Running $command\n";
system($command) == 0
  or die "\nFailed to run command\n$command\n\n";

# get back to original directory
chdir("../") or die "Failed to chdir to $subdir";
