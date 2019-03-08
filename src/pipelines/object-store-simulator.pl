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


my $mode = $0;

if    ($mode =~ m/ua$/) {
    if (rand() < 0.25) {          #  Simulate errors.
        print STDERR "Fail!\n";
        exit(1);
    }

    upload(@ARGV);
}

elsif ($mode =~ m/dx$/) {
    my $task = shift @ARGV;

    if      ($task eq "describe") {  describe(@ARGV);
    } elsif ($task eq "mv")       {  mv(@ARGV);
    } elsif ($task eq "rm")       {  rm(@ARGV);
    } elsif ($task eq "upload")   {  upload(@ARGV);
    } elsif ($task eq "download") {  download(@ARGV);
    }
}

else {
    die "Unknown mode '$mode'.\n";
}




sub describe (@) {
    my @args = @_;
    my $file;
    my $path = "";

    while (scalar(@args) > 0) {
        my $arg = shift @args;

        if      ($arg eq "--name") {
            $path = shift @args;
        }

        else {
            die "Unknown 'describe' option $arg\n";
        }
    }

    if (-e "$STASH/$path") {
        print "$path\n";
    }
}



sub mv (@) {
    my @args = @_;
    my $file;
    my $path = "";

    my $oldname = shift @args;
    my $newname = shift @args;

    print STDERR "DX:  STASH '$STASH'\n";
    print STDERR "DX:  oldname '$oldname'\n";
    print STDERR "DX:  newname '$newname'\n";
    print STDERR "DX: '$STASH/$oldname' -> '$STASH/$newname'\n";

    rename("$STASH/$oldname", "$STASH/$newname");
}



sub rm (@) {
    my @args = @_;
    my $file;
    my $path = "";
    my $recursive = 0;

    while (scalar(@args) > 0) {
        my $arg = shift @args;

        if      ($arg eq "--recursive") {    #  NOT SUPPORTED!
            $recursive = 1;
            next;
        }

        unlink("$STASH/$arg");
    }
}



#  Handles two variants.
#
#  dx upload  --path PR:NS/NA
#
#  dx-ua      --project PR --folder FL --name NA <path-to-file>
#                PR
#                FL must begin with a /
#                NA is just text, any 'directories' implied in it are part of the name
#
#  Canu uploads to
#    --project getGlobal("objectStoreProject");
#    --folder  getGlobal("objectStoreNameSpace");
#
sub upload (@) {
    my @args = @_;

    my $path;

    my $project;
    my $folder;
    my $name;

    my $file;

    while (scalar(@args) > 0) {
        my $arg = shift @args;

        if    ($arg eq "--path") {
            $path = shift @args;
        }

        elsif ($arg eq "--project") {
            $project = shift @args;
        }

        elsif ($arg eq "--folder") {
            $folder = shift @args;
        }

        elsif ($arg eq "--name") {
            $name = shift @args;
        }

        elsif (($arg eq "--wait") ||
               ($arg eq "--wait-on-close") ||
               ($arg eq "--do-not-compress") ||
               ($arg eq "--no-progress")) {
        }

        elsif (!defined($file)) {
            $file = $arg;
        }

        else {
            die "Unknown 'upload' option $arg\n";
        }
    }

    #  Check that the input file exists.

    die "dx upload - no input file supplied.\n"      if (!defined($file));
    die "dx upload - input file $file not found.\n"  if (! -e $file);

    #  If path exists, we're pretending to be 'dx upload'.

    if (defined($path)) {
        exit(1);
        system("mkdir -p $STASH/" . dirname($path) . " 2> /dev/null");
        system("cp -fp \"$file\" \"$STASH/$path\" 2> /dev/null");
    }

    #  Otherwise, we're pretending to be 'ua'.

    else {
        #print STDERR "dx-ua '$file' -> '$STASH' / '$project' : '$folder' / '$name'\n";

        system("mkdir -p $STASH/$project:$folder 2> /dev/null");
        system("cp -fp \"$file\" \"$STASH/$project:$folder/$name\" 2> /dev/null");
    }
}



sub download (@) {
    my @args = @_;
    my $file;
    my $path;

    while (scalar(@args) > 0) {
        my $arg = shift @args;

        if      ($arg eq "--output") {
            $file = shift @args;
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
    die "dx download - don't want to use stdout.\n"     if ($file eq "-");

    exit(0)  if (! -e "$STASH/$path");

    print STDERR "dx download '$STASH' / '$path' -> '$file'\n";

    system("cp -fp \"$STASH/$path\" \"$file\" 2> /dev/null");
}


