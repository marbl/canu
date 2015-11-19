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
 #    Brian P. Walenz on 2013-SEP-22
 #      are Copyright 2013 J. Craig Venter Institute, and
 #      are subject to the GNU General Public License version 2
 #
 #    Brian P. Walenz beginning on 2015-OCT-12
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

use strict;

my %IIDtoNAME;

my $fastqUIDmap = shift @ARGV;

print STDERR "Loading UID map from '$fastqUIDmap'.\n";

open (F, "< $fastqUIDmap") or die "Failed to open '$fastqUIDmap'\n";
while (<F>) {
    my @v = split '\s+', $_;

    if      (scalar(@v) == 3) {
        $IIDtoNAME{$v[1]} = $v[2];

    } elsif (scalar(@v) == 6) {
        $IIDtoNAME{$v[1]} = $v[2];
        $IIDtoNAME{$v[4]} = $v[5];

    } else {
        die "unknown format '$_'\n";
    }
}
close(F);


while (<STDIN>) {
    $_ =~ s/^\s+//;
    $_ =~ s/\s+$//;

    my @v = split '\s+', $_;

    die "Didn't find IID '$v[0]' in overlap '$_'.\n"  if (!exists($IIDtoNAME{$v[0]}));
    die "Didn't find IID '$v[0]' in overlap '$_'.\n"  if (!exists($IIDtoNAME{$v[1]}));

    $v[0] = $IIDtoNAME{$v[0]};
    $v[1] = $IIDtoNAME{$v[1]};

    print join("\t", @v), "\n";
}
