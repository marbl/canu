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
#  Read in fragments from fasta-format sequence and quality files,
#  writes version 2 fragments.
#
#  Assumes that ALL READS in the input are in the same library.  If
#  not, well, there is then no way to tell which library the unmated
#  reads are really from.  That library is assumed to be in a DST
#  record that appears BEFORE all reads.
#
#  If you give it "-v vecFile" it'll also populate the vector clear
#  range.

use strict;

my $vec;
my %clv;
my %clq;  #  Currently, we never have this info

my $clvFound = 0;
my $clvNotFound = 0;
my $clqFound = 0;
my $clqNotFound = 0;

my $noOBT = 0;
my $isNCBIqv = 0;

my $seqfile;
my $qltfile;
my $matefile;

my $is454  = 0;

my $lib    = undef;
my $mean   = 0.0;
my $stddev = 0.0;

my $idregex = "^>(\\S+)";

my $srcstr;

{
    local $, = " ";
    $srcstr = "$0 @ARGV";
}

my $err = 0;
while (scalar(@ARGV) > 0) {
    my $arg = shift @ARGV;
    if      ($arg eq "-v") {
        $vec = shift @ARGV;

    } elsif ($arg eq "-noobt") {
        $noOBT = 1;

    } elsif ($arg eq "-s") {
        $seqfile = shift @ARGV;
    } elsif ($arg eq "-q") {
        $qltfile = shift @ARGV;

    } elsif ($arg eq "-m") {
        $matefile = shift @ARGV;

    } elsif ($arg eq "-454") {
        $is454 = 1;

    } elsif ($arg eq "-l") {
        $lib = shift @ARGV;
    } elsif ($arg eq "-mean") {
        $mean = shift @ARGV;
    } elsif ($arg eq "-stddev") {
        $stddev = shift @ARGV;

    } elsif ($arg eq "-idregex") {
        $idregex = shift @ARGV;

    } else {
        $err++;
    }
}
if (defined($matefile) && ($mean == 0 || $stddev == 0)){
   print STDERR "A mean and standard deviation are required when a matepairing is provided.\n";
   $err++;
}
if (($err) || (!defined($seqfile)) || (!defined($qltfile)) || (!defined($lib))) {
    print STDERR "usage: $0 [options] -l libraryname -s seq.fasta -q qlt.fasta > new.frg\n";
    print STDERR "  -v vector-clear-file     A file of 'readUID vecBeg vecEnd', one per line, that is the vector clear range.\n";
    print STDERR "  -noobt                   Set the 'doNotOverlapTrim' library feature.\n";
    print STDERR "  -454                     Set library features appropriate for 454 reads (see also sffToCA).\n";
    print STDERR "  -idregex pattern         Use this perl regex to extract the read name from the seq defline.\n";
    print STDERR "  -l libraryname           Name of the library; freeformat text.\n";
    print STDERR "  -mean m                  Insert has mean size of m.\n";
    print STDERR "  -stddev s                Insert has std dev of s.\n";
    print STDERR "  -s seq                   Fasta file of sequences.\n";
    print STDERR "  -q qual                  Fasta file of quality values.\n";
    print STDERR "  -m matepairing           A file of pairs of read UIDs for mated reads, one pair per line, whitespace separated.\n";
    exit(1);
}

if (defined($vec)) {
    open(F, "< $vec") or die "Failed to open '$vec'\n";
    while (<F>) {
        s/^\s+//;
        s/\s$//;
        my @v = split '\s+', $_;
        $clv{$v[0]} = "$v[1],$v[2]";
    }
    close(F);
    #print STDERR "Read vector info for ", scalar(keys %clv), " reads.\n";
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
if ($is454) {
    print "doNotTrustHomopolymerRuns=1\n";
    print "discardReadsWithNs=1\n";
    print "doNotQVTrim=1\n";
    print "deletePerfectPrefixes=1\n";
} else {
    print "doNotTrustHomopolymerRuns=0\n";
    print "discardReadsWithNs=0\n";
    print "doNotQVTrim=0\n";
    print "deletePerfectPrefixes=0\n";
}
print "doNotOverlapTrim=$noOBT\n";
print "hpsIsFlowGram=0\n";
print "hpsIsPeakSpacing=0\n";
print "doNotOverlapTrim=$noOBT\n";
print "isNotRandom=0\n";
print ".\n";
print "}\n";

open(SEQ, "< $seqfile") or die "Failed to open '$seqfile'\n";
open(QLT, "< $qltfile") or die "Failed to open '$qltfile'\n";

my ($seqid, $seq) = readFasta();
my ($qltid, $qlt) = readQual();

while (defined($seq) && defined($qlt)) {

    if ($seqid ne $qltid) {
        die "Misordered seq/qual; got seq='$seqid' and qual='$qltid'\n";
    }

    print "{FRG\n";
    print "act:A\n";
    print "acc:$seqid\n";
    print "rnd:1\n";
    print "sta:G\n";
    print "lib:$lib\n";
    print "pla:0\n";
    print "loc:0\n";
    print "src:\n";
    print ".\n";
    print "seq:\n";
    print "$seq\n";
    print ".\n";
    print "qlt:\n";
    print "$qlt\n";
    print ".\n";
    print "hps:\n";
    print ".\n";
    if (defined($clv{$seqid})) {
        $clvFound++;
        print "clv:$clv{$seqid}\n";
    } else {
        $clvNotFound++;
    }
    if (defined($clq{$seqid})) {
        $clqFound++;
        print "clq:$clq{$seqid}\n";
    } else {
        $clqNotFound++;
    }
    print "clr:0,", length($seq), "\n";
    print "}\n";

    if (length($seq) != length($qlt)) {
        die "Length of sequence (" . length($seq) . ") and quality (" . length($qlt) . ") do not agree for $seqid.\n";
    }

    ($seqid, $seq) = readFasta();
    ($qltid, $qlt) = readQual();
}

close(SEQ);
close(QLT);

my $mateline = 1;
my $matefail = 0;

if (defined($matefile)) {
    open(F, "< $matefile") or die "Failed to open '$matefile'\n";
    while (<F>) {
        $_ =~ s/^\s+//;
        $_ =~ s/\s+$//;

        if ($_ ne "") {
            my ($a, $b, $c) = split '\s+', $_;

            if ((!defined($a)) || (!defined($b)) || (defined($c))) {
                print STDERR "Invalid mate pair '$_' on line $mateline -- NO MATE CREATED.\n";
                $matefail++;
            }

            print "{LKG\n";
            print "act:A\n";
            print "frg:$a\n";
            print "frg:$b\n";
            print "}\n";
        }

        $mateline++;
    }
    close(F);
}

print "{VER\n";
print "ver:1\n";
print "}\n";


if (defined(%clv) && ($clvNotFound > 0)) {
    print STDERR "Updated $clvFound vector clear ranges ($clvNotFound NOT updated).\n";
}

#print STDERR "Updated $clqFound quality clear ranges ($clqNotFound NOT updated).\n";


#
#  Both routines stolen from AS_MSG/tracedb-to-frg.pl
#

my $fhdr;
my $qhdr;

sub readFasta {
    my $fstr;

    while (!eof(SEQ)) {
        $_ = <SEQ>;
        chomp;

        if (m/^>/) {
            my $ret = $fhdr;

            if      (m/$idregex/) {
                $fhdr = $1;
            } else {
                die "Failed to parse an ID out of the sequence defline '$_'\n";
            }

            if (defined($ret)) {
                return($ret, $fstr);
            }
        } else {
            my $q = $_;
            $q =~ s/^\s+//;
            $q =~ s/\s+$//;
            $fstr .= $q;
        }
    }

    #   Got eof, return whatever
    my $ret = $fhdr;
    undef $fhdr;
    return($ret, $fstr);
}


sub readQual {
    my $qstr;

    while (!eof(QLT)) {
        $_ = <QLT>;
        chomp;

        if (m/^>/) {
            my $ret = $qhdr;

            if      (m/$idregex/) {
                $qhdr = $1;
            } else {
                die "Failed to parse an ID out of the quality defline '$_'\n";
            }

            if (defined($ret)) {
                return($ret, $qstr);
            }
        } else {
            my $q = $_;
            $q =~ s/^\s+//;
            $q =~ s/\s+$//;

            if ($q =~ m/\s+/) {
                $isNCBIqv = 1;

                foreach my $qv (split '\s+', $q) {
                    if ($qv > 60) {$qv = 60;}
                    $qstr .= chr(ord('0') + $qv);
                }
            } elsif ($isNCBIqv) {
                $qstr .= chr(ord('0') + $q);
            } else {
                $qstr .= $q;
            }
        }
    }

    #   Got eof, return whatever
    my $ret = $qhdr;
    undef $qhdr;
    return($ret, $qstr);
}

if ($matefail > 0) {
    print STDERR "\nThere were $matefail errors in the matefile.\n";
} else {
    print STDERR "Finished successfully.\n";
}

exit($matefail != 0);
