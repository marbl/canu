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
 #    Brian P. Walenz beginning on 2015-OCT-12
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

use strict;

#  Convert BLASR default output to OBT overlaps:  aIID bIID [f|r] aBgn aEnd bBgn bEnd error

my %IDmap;

open(F, "< $ARGV[0]") or die "Failed to open '$ARGV[0]' for reading sequence names.\n";
while (<F>) {
    my ($uid, $iid, $name) = split '\s+', $_;

    $IDmap{$name} = $iid;
}
close(F);


while (<STDIN>) {
    my @v = split '\s+', $_;

    #  The first read seems to have a sub-read range appended.

    if ($v[0] =~ m/^(.*)\/\d+_\d+$/) {
        $v[0] = $1;
    }

    my $aiid  = $IDmap{$v[0]};
    my $biid  = $IDmap{$v[1]};
    my $fr    = ($v[8] == 0) ? "f" : "r";
    my $error = 100.0 - $v[3];

    ($v[9], $v[10]) = ($v[10], $v[9])  if ($v[8] == 1);

    die "No A iid found for '$v[0]'\n" if (!defined($aiid));
    die "No B iid found for '$v[1]'\n" if (!defined($biid));

    next if ($aiid == $biid);

    #  HACK!
    #$error = 0.01;

    print "$aiid\t$biid\t$fr\t$v[5]\t$v[6]\t$v[9]\t$v[10]\t$error\n";
}
