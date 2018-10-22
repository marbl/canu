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
 #    Brian P. Walenz beginning on 2018-MAY-04
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

use strict;
use File::Basename;

my $STASH;

$STASH = "/assembly/objectstore"    if (-d "/assembly/objectstore");

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
            die "Unknown 'describe' option $arg\n";
        }
    }

    if (-e "$STASH/$path") {
        print "$path\n";
    }
}



if ($task eq "mv") {
    my $file;
    my $path = "";
    my $wait = 0;

    my $oldname = shift @ARGV;
    my $newname = shift @ARGV;

    print STDERR "DX:  STASH '$STASH'\n";
    print STDERR "DX:  oldname '$oldname'\n";
    print STDERR "DX:  newname '$newname'\n";
    print STDERR "DX: '$STASH/$oldname' -> '$STASH/$newname'\n";

    rename("$STASH/$oldname", "$STASH/$newname");
}



if ($task eq "rm") {
    my $file;
    my $path = "";
    my $recursive = 0;

    while (scalar(@ARGV) > 0) {
        my $arg = shift @ARGV;

        if      ($arg eq "--recursive") {    #  NOT SUPPORTED!
            $recursive = 1;
            next;
        }

        unlink("$STASH/$arg");
    }
}



if ($task eq "upload") {
    my $file;
    my $path;
    my $wait    = 0;
    my $parents = 0;

    while (scalar(@ARGV) > 0) {
        my $arg = shift @ARGV;

        if      ($arg eq "--path") {
            $path = shift @ARGV;
        }

        elsif ($arg eq "--wait") {
            $wait = 1;
        }

        elsif ($arg eq "--parents") {
            $parents = 1;
        }

        elsif ($arg eq "--no-progress") {
        }

        elsif (!defined($file)) {
            $file = $arg;
        }

        else {
            die "Unknown 'upload' option $arg\n";
        }
    }

    #  Copy local file $file, assumed to be in this directory, to the stash as $path.

    die "dx upload - no stash path supplied.\n"      if (!defined($path));
    die "dx upload - no input file supplied.\n"      if (!defined($file));
    die "dx upload - input file $file not found.\n"  if (($file ne "-") && (! -e $file));

    system("mkdir -p $STASH/" . dirname($path) . " 2> /dev/null")   if ($parents);

    if ($file eq "-") {
        system("dd status=none \"of=$STASH/$path\" 2> /dev/null");
    } else {
        system("cp -fp \"$file\" \"$STASH/$path\" 2> /dev/null");
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

        elsif ($arg eq "--no-progress") {
        }

        elsif ($arg eq "--overwrite") {
        }

        elsif (!defined($path)) {
            $path = $arg;
        }

        else {
            die "Unknown 'download' option $arg\n";
        }
    }

    #  Copy local file $file, assumed to be in this directory, to the stash as $path.

    die "dx download - no stash path supplied.\n"       if (!defined($path));
    die "dx download - no output file supplied.\n"      if (!defined($file));
    #die "dx download - stash path $path not found.\n"   if (! -e "$STASH/$path");

    print STDERR "FETCH 'cp -fp \"$STASH/$path\" \"$file\" 2> /dev/null'\n";

    exit(0)  if (! -e "$STASH/$path");

    if ($file eq "-") {
        system("dd status=none \"if=$STASH/$path\" 2> /dev/null");
    } else {
        system("cp -fp \"$STASH/$path\" \"$file\" 2> /dev/null");
    }
}


