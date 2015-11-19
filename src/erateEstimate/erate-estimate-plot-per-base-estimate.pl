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
 #    src/AS_BAT/erate-estimate-plot-per-base-estimate.pl
 #
 #  Modifications by:
 #
 #    Brian P. Walenz from 2014-NOV-15 to 2015-AUG-07
 #      are Copyright 2014-2015 Battelle National Biodefense Institute, and
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

# 279 280
# 281 284 285 interesting

my @iid;  #= ( 281, 284, 285, 3360, 4206, 4944, 5118 );

push @iid, 2276;  #  normal
push @iid, 281;
push @iid, 284;
push @iid, 285;
push @iid, 3360;
push @iid, 4206;
push @iid, 4944;
push @iid, 5118;

foreach my $iid (@iid) {
    my $prefix = "iid$iid";
    my @erates;

    open(F, "overlapStore -p $iid test.ovlStore test.gkpStore CLR |");

    open(D, "> $prefix.stats.dat");
    open(E, "> $prefix.erate.dat");

    #  Ignore the header.
    $_ = <F>;

    while (<F>) {
        s/^\s+//;
        s/\s+$//;

        #my @v = split '\s+', $_;

        #my $biid = $v[0];
        #my $abgn = $v[2];
        #my $aend = $v[3];

        my ($biid, $abgn, $aend, $erate);

        if (m/(\d+)\s+A:\s+(\d+)\s+(\d+)\s+\(\s*\d+\)\s+B:\s+\d+\s+\d+\s+\(\s*\d+\)\s+(\d+.\d+)%/) {
            $biid  = $1;
            $abgn  = $2;
            $aend  = $3;
            $erate = $4;
        } else {
            print "NOPE $_\n";
        }

        for (my $ii=$abgn; $ii<=$aend; $ii++) {
            if (!defined($erates[$ii])) {
                $erates[$ii]  = "$erate";
            } else {
                $erates[$ii] .= ":$erate";
            }
        }

        #print E "$abgn\t$aend\t$erate\n";
    }

    close(F);


    my $maxPoint = scalar(@erates);

    for (my $pp=0; $pp<$maxPoint; $pp++) {
        my @vals = split ':', $erates[$pp];
        my $depth = scalar (@vals);

        @vals = sort {$a <=> $b } @vals;

        my $min = $vals[0];
        my $max = $vals[$depth-1];

        my $sum = 0;
        foreach my $v (@vals) {
            $sum += $v;
        }

        my $ave = $sum / $depth;

        print D "$pp\t$min\t$ave\t$max\t$depth\n";

        foreach my $v (@vals) {
            print E "$pp\t$v\n";
        }
    }

    close(D);
    close(E);

    system("overlapStore -p $iid test.ovlStore test.gkpStore CLR > $prefix.ovlPicture");

    open(F, "| gnuplot");
    print F "set terminal 'png' size 1280,800\n";
    print F "set output   '$prefix.png'\n";
    print F "plot [][0:30]";
    print F "  '$prefix.stats.dat' using 1:2 with lines title 'min % identity',";
    print F "  '$prefix.stats.dat' using 1:3 with lines title 'ave % identity',";
    print F "  '$prefix.stats.dat' using 1:4 with lines title 'max % identity',";
    print F "  '$prefix.stats.dat' using 1:5 with lines title 'depth',";
    #print F "  '$prefix.erate.dat' using 1:2 with points ps 0.2 title 'erate'\n";
    print F "  '$prefix.erate.dat' using 1:3 with points ps 0.3 title 'bgnerates',";
    print F "  '$prefix.erate.dat' using 2:3 with points ps 0.3 title 'enderates'\n";
    close(F);
}
