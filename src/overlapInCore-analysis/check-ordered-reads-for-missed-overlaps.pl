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

my $dataset;

$dataset = "BL";
$dataset = "CAordered";

my $maxB = 0;
my $maxA = 0;

sub processIIDs ($@) {
    my $thisIID = shift @_;
    my @iids    = sort { $a <=> $b } @_;

    return if (scalar(@iids) == 0);

    my $minIID = $iids[0];
    my $curIID = 0;
    my $maxIID = $iids[scalar(@iids) - 1];

    my $expectedB = 0;
    my $expectedA = 0;

    my $foundB = 0;
    my $foundA = 0;

    my $missingB = 0;
    my $missingA = 0;

    foreach my $iid (@iids) {
        #print "$iid -- $curIID\n";

        next if ($iid < $thisIID - 30);  #  Ignore false overlaps
        next if ($iid > $thisIID + 30);

        $expectedB = $thisIID - $iid   if ($expectedB < $thisIID - $iid);
        $expectedA = $iid - $thisIID   if ($expectedA < $iid - $thisIID);

        $curIID = $iid - 1  if ($curIID == 0);  #  If curIID==0, first time here, set to make no missing.

        if ($iid < $thisIID) {
            $foundB++;
            if ($curIID + 1 < $iid) {
                #print STDERR "  missingB $curIID - $iid\n";
                $missingB += $iid - 1 - $curIID;
            }

        } else {
            $foundA++;
            if ($curIID + 1 < $iid) {
                print STDERR "  missingA $curIID - $iid\n";
                $missingA += $iid - 1 - $curIID;
            }
        }

        $curIID   = $iid;

        $curIID++  if ($curIID + 1 == $thisIID);
    }

    my $found    = $foundB    + $foundA;
    my $expected = $expectedB + $expectedA;
    my $missing  = $missingB  + $missingA;

    if ($missingA > 0) {
        my $fracB  = int(10000 * $missingB / $expectedB) / 100;
        my $fracA  = int(10000 * $missingA / $expectedA) / 100;

        print "$thisIID range $expectedB/$expectedA found $foundB/$foundA missing $missingB/$missingA $fracB%/$fracA%\n";
    }
}






my $lastIID = 0;
my @iids;

open(F, "overlapStore -d JUNKTEST3/$dataset/test.ovlStore |");
while (<F>) {
    s/^\s+//;
    s/\s+$//;

    #  v[0] = iid
    #  v[1] = iid
    #  v[2] = orient
    #  v[3] = hang
    #  v[4] = hang
    #  v[5] = ident
    #  v[6] = ident

    my @v = split '\s+', $_;

    if ($v[0] != $lastIID) {
        #print STDERR "PROCESS $lastIID ", scalar(@iids), "\n";

        processIIDs($lastIID, @iids);

        undef @iids;
        $lastIID = $v[0];
    }

    #print STDERR "$v[0] $v[1]\n";

    push @iids, $v[1];
}
close(F);
