#!/usr/bin/perl

#
###########################################################################
#
# This file is part of Celera Assembler, a software program that
# assembles whole-genome shotgun reads into contigs and scaffolds.
# Copyright (C) 2007, J. Craig Venter Instititue.
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
#  Read in fragments from Solexa fastq-format sequence and quality files,
#  writes version 2 fragments.
#

use strict;

my $seqL;
my $seqR;

my $lib    = undef;
my $mean   = 0.0;
my $stddev = 0.0;

my $srcstr;

{
    local $, = " ";
    $srcstr = "$0 @ARGV";
}

my $err = 0;
while (scalar(@ARGV) > 0) {
    my $arg = shift @ARGV;
    if      ($arg eq "-sl") {
        $seqL = shift @ARGV;
    } elsif ($arg eq "-sr") {
        $seqR = shift @ARGV;

    } elsif ($arg eq "-l") {
        $lib = shift @ARGV;
    } elsif ($arg eq "-mean") {
        $mean = shift @ARGV;
    } elsif ($arg eq "-stddev") {
        $stddev = shift @ARGV;

    } else {
        $err++;
    }
}
if (($err) ||
    (!defined($seqL)) ||
    (defined($seqL) && defined($seqR) && !defined($lib))) {
    print STDERR "usage: $0 [options] -l libraryname -s seq.fasta -q qlt.fasta > new.frg\n";
    print STDERR "  -sl seq                  Fasta file of sequences, left side.\n";
    print STDERR "  -sr seq                  Fasta file of sequences, right side.\n";
    print STDERR "  -l libraryname           Name of the library; freeformat text.\n";
    print STDERR "  -mean m                  Insert has mean size of m.\n";
    print STDERR "  -stddev s                Insert has std dev of s.\n";
    exit(1);
}

print "{VER\n";
print "ver:2\n";
print "}\n";

print "{LIB\n";
print "act:A\n";
print "acc:$lib\n";
if ($mean > 0) {
    print "ori:I\n";
    print "mea:$mean\n";
    print "std:$stddev\n";
} else {
    print "ori:U\n";
    print "mea:0.0\n";
    print "std:0.0\n";
}
print "src:\n";
print "$srcstr\n";
print ".\n";
print "nft:9\n";
print "fea:\n";
print "doNotTrustHomopolymerRuns=0\n";
print "discardReadsWithNs=0\n";
print "doNotQVTrim=1\n";
print "deletePerfectPrefixes=0\n";
print "doNotOverlapTrim=1\n";
print "hpsIsFlowGram=0\n";
print "hpsIsPeakSpacing=0\n";
print "doNotOverlapTrim=1\n";
print "isNotRandom=0\n";
print ".\n";
print "}\n";

open(SEQL, "< $seqL") or die "Failed to open '$seqL'\n";
open(SEQR, "< $seqR") or die "Failed to open '$seqR'\n";

my ($seqLid, $seqL, $qltL) = readSeq(*SEQL);
my ($seqRid, $seqR, $qltR) = readSeq(*SEQR);

while (defined($seqL)) {
    print "{FRG\n";
    print "act:A\n";
    print "acc:$seqLid\n";
    print "rnd:1\n";
    print "sta:G\n";
    print "lib:$lib\n";
    print "pla:0\n";
    print "loc:0\n";
    print "src:\n";

    print ".\n";
    print "seq:\n";
    print "$seqL\n";
    print ".\n";
    print "qlt:\n";
    print "$qltL\n";
    print ".\n";
    print "hps:\n";
    print ".\n";
    print "clr:0,", length($seqL), "\n";
    print "}\n";

    print "{FRG\n";
    print "act:A\n";
    print "acc:$seqRid\n";
    print "rnd:1\n";
    print "sta:G\n";
    print "lib:$lib\n";
    print "pla:0\n";
    print "loc:0\n";
    print "src:\n";
    print ".\n";
    print "seq:\n";
    print "$seqR\n";
    print ".\n";
    print "qlt:\n";
    print "$qltR\n";
    print ".\n";
    print "hps:\n";
    print ".\n";
    print "clr:0,", length($seqR), "\n";
    print "}\n";

    print "{LKG\n";
    print "act:A\n";
    print "frg:$seqLid\n";
    print "frg:$seqRid\n";
    print "}\n";

    ($seqLid, $seqL, $qltL) = readSeq(*SEQL);
    ($seqRid, $seqR, $qltR) = readSeq(*SEQR);
}

close(SEQL);
close(SEQR);

print "{VER\n";
print "ver:1\n";
print "}\n";


sub readSeq(*) {
    local *F = shift;

    my ($ids, $seq, $idq, $qlt);

    $ids = <F>;  chomp $ids;
    $seq = <F>;  chomp $seq;
    $idq = <F>;  chomp $idq;
    $qlt = <F>;  chomp $qlt;

    if (!defined($seq) || !defined($qlt)) {
        return(undef, undef, undef);
    }

    if ($ids =~ m/.(.*):(\d+):(\d+):(\d+):(\d+)\/([0-9])/) {
        my $E = ($6 == 1) ? "L" : "R";
        $ids = "$1-$2-$3-$4-$5-$E";
    }

    if ($idq =~ m/.(.*):(\d+):(\d+):(\d+):(\d+)\/([0-9])/) {
        my $E = ($6 == 1) ? "L" : "R";
        $idq = "$1-$2-$3-$4-$5-$E";
    }

    if ($ids ne $idq) {
        print STDERR "WARNING: $ids != $idq\n";
    }

    return($ids, $seq, $qlt);
}
