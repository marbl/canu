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
 #    src/AS_MER/gkrpt.pl
 #
 #  Modifications by:
 #
 #    Brian P. Walenz from 2008-JAN-03 to 2013-AUG-01
 #      are Copyright 2008,2013 J. Craig Venter Institute, and
 #      are subject to the GNU General Public License version 2
 #
 #    Brian P. Walenz on 2015-APR-10
 #      are Copyright 2015 Battelle National Biodefense Institute, and
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

my %kmers;

while (!eof(STDIN)){
    my $count = <STDIN>;  chomp $count;
    my $kmer  = <STDIN>;  chomp $kmer;
    $count =~ s/>//g;
    $count = int($count);
    $kmer =~ tr/acgt/ACGT/;
    my $numbad = $kmer =~ tr/ACGT//c;
    if ($numbad > 0){
        die "$numbad nonacgtACGT characters in kmer $kmer\n";
    }
    if (defined($kmers{$kmer})){
        die "kmer repeated in input $kmer\n";
    }
    $kmers{$kmer} = $count;
}

foreach my $kmer (sort { $kmers{$b} <=> $kmers{$a} } (keys %kmers)){
    my $startkmer = my $curkmer = $kmer;
    my $startcount = my $curcount = $kmers{$kmer};
    if ($curcount < 0) {
        next;
    } else {
        $kmers{$kmer} = - $curcount;
    }
    my $currpt = $curkmer;
    for ( ; ; ) {
        my $nextcount = -1;
        my $realcount = -1;
        my $maxkmer;
        my $realkmer;
        my $tmpkmer;
        my $maxcount = -1;
        $curkmer = (substr $curkmer, 1) . "A";
        my $rckmer = reverse $curkmer;
        $rckmer =~ tr/ACGT/TGCA/;
        if (defined($kmers{$curkmer})){
            $nextcount = $kmers{$curkmer};
            $tmpkmer = $curkmer;
        } else {
            if (defined($kmers{$rckmer})){
                $nextcount = $kmers{$rckmer};
                $tmpkmer = $rckmer;
            } else {
                $nextcount = -1;
            }
        }
        if ((abs $nextcount) > $maxcount){
            $maxcount = abs $nextcount;
            $realcount = $nextcount;
            $realkmer = $tmpkmer;
            $maxkmer = $curkmer;
        }
        substr($curkmer, -1) = "C";
        if (defined($kmers{$curkmer})){
            $nextcount = $kmers{$curkmer};
            $tmpkmer = $curkmer;
        } else {
            substr($rckmer, 0, 1) = "G";
            if (defined($kmers{$rckmer})){
                $nextcount = $kmers{$rckmer};
                $tmpkmer = $rckmer;
            } else {
                $nextcount = -1;
            }
        }
        if ((abs $nextcount) > $maxcount){
            $maxcount = abs $nextcount;
            $realcount = $nextcount;
            $realkmer = $tmpkmer;
            $maxkmer = $curkmer;
        }
        substr($curkmer, -1) = "G";
        if (defined($kmers{$curkmer})){
            $nextcount = $kmers{$curkmer};
            $tmpkmer = $curkmer;
        } else {
            substr($rckmer, 0, 1) = "C";
            if (defined($kmers{$rckmer})){
                $nextcount = $kmers{$rckmer};
                $tmpkmer = $rckmer;
            } else {
                $nextcount = -1;
            }
        }
        if ((abs $nextcount) > $maxcount){
            $maxcount = abs $nextcount;
            $realcount = $nextcount;
            $realkmer = $tmpkmer;
            $maxkmer = $curkmer;
        }
        substr($curkmer, -1) = "T";
        if (defined($kmers{$curkmer})){
            $nextcount = $kmers{$curkmer};
            $tmpkmer = $curkmer;
        } else {
            substr($rckmer, 0, 1) = "A";
            if (defined($kmers{$rckmer})){
                $nextcount = $kmers{$rckmer};
                $tmpkmer = $rckmer;
            } else {
                $nextcount = -1;
            }
        }
        if ((abs $nextcount) > $maxcount){
            $maxcount = abs $nextcount;
            $realcount = $nextcount;
            $realkmer = $tmpkmer;
            $maxkmer = $curkmer;
        }
        if (($realcount < 0) || ($realcount < ($curcount / 2))) {
            last;
        } else {
            $curkmer = $maxkmer;
            $curcount = $realcount;
            $kmers{$realkmer} = - $realcount;
            $currpt .= (substr $curkmer, -1);
        }
    }
    $curcount = $startcount;
    $curkmer = $startkmer;

    for ( ; ; ) {
        my $nextcount = -1;
        my $realcount = -1;
        my $maxkmer;
        my $realkmer;
        my $tmpkmer;
        my $maxcount = -1;
        $curkmer = "A" . (substr $curkmer, 0, -1);
        my $rckmer = reverse $curkmer;
        $rckmer =~ tr/ACGT/TGCA/;
        if (defined($kmers{$curkmer})){
            $nextcount = $kmers{$curkmer};
            $tmpkmer = $curkmer;
        } else {
            if (defined($kmers{$rckmer})){
                $nextcount = $kmers{$rckmer};
                $tmpkmer = $rckmer;
            } else {
                $nextcount = -1;
            }
        }
        if ((abs $nextcount) > $maxcount){
            $maxcount = abs $nextcount;
            $realcount = $nextcount;
            $realkmer = $tmpkmer;
            $maxkmer = $curkmer;
        }
        substr($curkmer, 0, 1) = "C";
        if (defined($kmers{$curkmer})){
            $nextcount = $kmers{$curkmer};
            $tmpkmer = $curkmer;
        } else {
            substr($rckmer, -1) = "G";
            if (defined($kmers{$rckmer})){
                $nextcount = $kmers{$rckmer};
                $tmpkmer = $rckmer;
            } else {
                $nextcount = -1;
            }
        }
        if ((abs $nextcount) > $maxcount){
            $maxcount = abs $nextcount;
            $realcount = $nextcount;
            $realkmer = $tmpkmer;
            $maxkmer = $curkmer;
        }
        substr($curkmer, 0, 1) = "G";
        if (defined($kmers{$curkmer})){
            $nextcount = $kmers{$curkmer};
            $tmpkmer = $curkmer;
        } else {
            substr($rckmer, -1) = "C";
            if (defined($kmers{$rckmer})){
                $nextcount = $kmers{$rckmer};
                $tmpkmer = $rckmer;
            } else {
                $nextcount = -1;
            }
        }
        if ((abs $nextcount) > $maxcount){
            $maxcount = abs $nextcount;
            $realcount = $nextcount;
            $realkmer = $tmpkmer;
            $maxkmer = $curkmer;
        }
        substr($curkmer, 0, 1) = "T";
        if (defined($kmers{$curkmer})){
            $nextcount = $kmers{$curkmer};
            $tmpkmer = $curkmer;
        } else {
            substr($rckmer, -1) = "A";
            if (defined($kmers{$rckmer})){
                $nextcount = $kmers{$rckmer};
                $tmpkmer = $rckmer;
            } else {
                $nextcount = -1;
            }
        }
        if ((abs $nextcount) > $maxcount){
            $maxcount = abs $nextcount;
            $realcount = $nextcount;
            $realkmer = $tmpkmer;
            $maxkmer = $curkmer;
        }
        if (($realcount < 0) || ($realcount < ($curcount / 2))) {
            last;
        } else {
            $curkmer = $maxkmer;
            $curcount = $realcount;
            $kmers{$realkmer} = - $realcount;
            $currpt = (substr $curkmer, 0, 1) . $currpt;
        }
    }
    if ((my $lenrpt = length $currpt) > $ARGV[0]) {
        print ">$startkmer $startcount $lenrpt\n$currpt\n";
    }
}
