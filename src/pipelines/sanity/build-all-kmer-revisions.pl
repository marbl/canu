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

#  Given a local rsync'd copy of the repository, check out all versions and compile.

my $kmersvn = "/work/NIGHTLY/kmer-svn"

open(F, "< $kmersvn/db/current");
my $latest = <F>;
chomp $latest;
close(F);


for (my $i=1917; $i<=$latest; $i++) {
    if (! -d "kmer$i") {
        print "Check out r$i\n";
        system("mkdir kmer$i");
        system("cd kmer$i && svn co -r $i file://$kmersvn/trunk . > kmer-checkout.err 2>&1");
    }

    if (! -d "kmer$i/FreeBSD-amd64") {
        print "Compile r$i\n";
        system("cd kmer$i && gmake install > kmer-build.err 2>& 1 &");
    }
}
