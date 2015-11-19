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

se strict;

#
###########################################################################
#
# This file is part of Celera Assembler, a software program that
# assembles whole-genome shotgun reads into contigs and scaffolds.
# Copyright (C) 2009, J. Craig Venter Instititue.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received (LICENSE.txt) a copy of the GNU General Public
# License along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
###########################################################################
#
#  Figures out what pieces ended up after trimming.  Use with generate-random-chimeric-fragments.pl
#
#  User unfriendly.

my $uid;
my $iid;

open(F, "/work/wgs/FreeBSD-amd64/bin/gatekeeper -dumpfragments /work/test/test.gkpStore |") or die;
while (<F>) {
    if (m/^fragmentIdent\s+=\s+(.*),(\d+)$/) {
        $uid = $1;
        $iid = $2;
    }
    if (m/^fragmentIsDeleted\s+=\s+1$/) {
        undef $uid;
        undef $iid;
    }

    if (m/^fragmentClear\s+=\s+OBTCHIMERA,(\d+),(\d+)$/) {
        my $bgn = $1;
        my $end = $2;

        my $b = 0;
        my $e = 0;

        if ($uid =~ m/(\d+)-/) {
            my @v = split '-', $uid;

            shift @v;  #  iid of read

            my $numMatches = 0;
            my $strMatches = "";
            my $badMatches = 0;

            while (scalar(@v) > 0) {
                my $typ = shift @v;
                my $len = shift @v;

                $e = $b + $len;

                if (($bgn < $e - 5) && ($b + 5 < $end)) {
                    $numMatches++;
                    $strMatches .= "$b-$e-$typ ";

                    if (($typ ne "seq") && ($typ ne "rev")) {
                        $badMatches++;
                    }
                }

                $b = $e;
            }

            if (($numMatches > 1) || ($badMatches > 0)) {
                my $ovl = `/work/wgs/FreeBSD-amd64/bin/overlapStore -g /work/test/test.gkpStore -dp $iid /work/test/0-overlaptrim/test.obtStore`;

                if (length($ovl) > 0) {
                    print "$uid\n";
                    print "$bgn-$end\n";
                    print "$strMatches\n";
                    print $ovl;
                    print "========================================================================================================================\n";
                }
            }
        }
    }
}
close(F);
