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
 #  This file is derived from:
 #
 #    src/AS_RUN/replaceUIDwithName-simple.pl
 #
 #  Modifications by:
 #
 #    Brian P. Walenz on 2013-SEP-22
 #      are Copyright 2013 J. Craig Venter Institute, and
 #      are subject to the GNU General Public License version 2
 #
 #    Brian P. Walenz on 2014-OCT-01
 #      are Copyright 2014 Battelle National Biodefense Institute, and
 #      are subject to the BSD 3-Clause License
 #
 #    Brian P. Walenz beginning on 2015-OCT-12
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

use strict;

my %UIDtoNAME;

#  Usage
#  fastqUIDmap posmnap posmap posmap

my $fastqUIDmap = shift @ARGV;

print STDERR "Loading UID map from '$fastqUIDmap'.\n";

open (F, "< $fastqUIDmap") or die "Failed to open '$fastqUIDmap'\n";
while (<F>) {
    chomp;

    my @v = split '\s+', $_;

    if      (scalar(@v) == 3) {
        $UIDtoNAME{$v[0]} = $v[2];

    } elsif (scalar(@v) == 6) {
        $UIDtoNAME{$v[0]} = $v[2];
        $UIDtoNAME{$v[3]} = $v[5];

    } else {
        die "unknown format '$_'\n";
    }
}
close(F);

my $inFile;
my $otFile;

while (scalar(@ARGV)) {
    $inFile = shift @ARGV;

    print STDERR "Renaming '$inFile' to '$inFile.UID'.\n";

    rename "$inFile", "$inFile.UID";

    open(F, "< $inFile.UID") or die "Failed to open '$inFile.UID' for reading\n";
    open(O, "> $inFile")     or die "Failed to open '$inFile' for writing\n";

    while (!eof(F)) {
        my $a = <F>;
        my @a = split '\s+', $a;

        foreach my $a (@a) {
            if (exists($UIDtoNAME{$a})) {
                $a = $UIDtoNAME{$a};
            }
        }

        print O join "\t", @a;
        print O "\n";
    }

    close(F);
    close(O);
}
