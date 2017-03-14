#!/usr/bin/env perl

###############################################################################
 #
 #  This file is part of canu, a software program that assembles whole-genome
 #  sequencing reads into contigs.
 #
 #  This software is based on:
 #    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 #    the 'kmer package' (http://kmer.sourceforge.net)
 #  both originally distributed by Applera Corporation under the GNU General
 #  Public License, version 2.
 #
 #  Canu branched from Celera Assembler at its revision 4587.
 #  Canu branched from the kmer project at its revision 1994.
 #
 #  Modifications by:
 #
 #    Brian P. Walenz beginning on 2015-MAR-14
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

use strict;
use File::Basename;

#  A simple wrapper to emulate an object store.  Implements:
#    describe --name <name>          - prints <name> if the store has the named object, nothing otherwise
#    upload --path <name> <file>     - uploads local <file> to the store as object <name>
#    download --output <file> <name> - downloads object <name> into <file>
#
#  As seen immediately below, requires a magic hardcoded path.  Files are copied to this directory
#  to become the objects.

my $STASH;

$STASH = "/assembly/STASH"          if (-d "/assembly/STASH");
$STASH = "/Users/walenzbp/STASH"    if (-d "/Users/walenzbp/STASH");

die "No STASH found\n"  if (!defined($STASH));

my $task = shift @ARGV;



if ($task eq "describe") {
    my $file;
    my $path = "";
    my $wait = 0;

    while (scalar(@ARGV) > 0) {
        my $arg = shift @ARGV;

        if      ($arg eq "--name") {
            $path = shift @ARGV;
        }

        else {
            die "Unknown option $arg\n";
        }
    }

    if (-e "$STASH/$path") {
        print "$path\n";
    }
}



if ($task eq "upload") {
    my $file;
    my $path;
    my $wait = 0;

    while (scalar(@ARGV) > 0) {
        my $arg = shift @ARGV;

        if      ($arg eq "--path") {
            $path = shift @ARGV;
        }

        elsif ($arg eq "--wait") {
            $wait = 1;
        }

        elsif (!defined($file)) {
            $file = $arg;
        }

        else {
            die "Unknown option $arg\n";
        }
    }

    #  Copy local file $file, assumed to be in this directory, to the stash as $path.

    die "dx upload - no stash path supplied.\n"      if (!defined($path));
    die "dx upload - no input file supplied.\n"      if (!defined($file));
    die "dx upload - input file $file not found.\n"  if (($file ne "-") && (! -e $file));

    system("mkdir -p $STASH/" . dirname($path));

    if ($file eq "-") {
        system("dd status=none of=$STASH/$path");
    } else {
        system("cp -fp $file $STASH/$path");
    }
}



if ($task eq "download") {
    my $file;
    my $path;
    my $wait = 0;

    while (scalar(@ARGV) > 0) {
        my $arg = shift @ARGV;

        if      ($arg eq "--output") {
            $file = shift @ARGV;
        }

        elsif (!defined($path)) {
            $path = $arg;
        }

        else {
            die "Unknown option $arg\n";
        }
    }

    #  Copy local file $file, assumed to be in this directory, to the stash as $path.

    die "dx download - no stash path supplied.\n"       if (!defined($path));
    die "dx download - no output file supplied.\n"      if (!defined($file));
    #die "dx download - stash path $path not found.\n"   if (! -e "$STASH/$path");

    exit(0)  if (! -e "$STASH/$path");

    if ($file eq "-") {
        system("dd status=none if=$STASH/$path");
    } else {
        system("cp -fp $STASH/$path $file");
    }
}


exit(0);
