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
#  Read in version 1 format fragment file, write version 2.
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
my $lib;
my $clvFound = 0;
my $clvNotFound = 0;
my $clqFound = 0;
my $clqNotFound = 0;

my $noOBT = 0;
my $mateStatus = "I";

my $err = 0;
while (scalar(@ARGV) > 0) {
    my $arg = shift @ARGV;
    if      ($arg eq "-v") {
        $vec = shift @ARGV;
    } elsif ($arg eq "-noobt") {
        $noOBT = 1;
    } elsif ($arg eq "-unmated") {
        $mateStatus = "U";
    } else {
        $err++;
    }
}
if ($err) {
    print STDERR "usage: $0 [-v vector-clear-file] [-noobt] [-unmated] < old.frg > new.frg\n";
    print STDERR "-v vector-clear-file     A file of 'readUID vecBeg vecEnd', one per line, that is the vector clear range.\n";
    print STDERR "-noobt                   Set the 'doNotOverlapTrim' library feature (used for 454 reads; see sffToCA instead).\n";
    print STDERR "-unmated                 Make all reads unmated.\n";
    exit(1);
}

if (defined($vec)) {
    open(F, "< $vec") or die "Failed to open '$vec'\n";
    while (<F>) {
        s/^\s+//;
        s/\s$//;

        my @v = split '\s+', $_;
        my @t = split '\t', $_;

        if      (scalar(@v) == 3) {
            $clv{$v[0]} = "$v[1],$v[2]";

        } elsif (scalar(@t) == 13) {
            $clv{$t[0]} = "$t[11],$t[12]";
        } elsif (scalar(@t) == 14) {
            $clv{$t[0]} = "$t[11],$t[12]";

        } else {
            die "Unknown metadata format in line '$_'\n";
        }
    }
    close(F);

    #print STDERR "$vec: Read vector info for ", scalar(keys %clv), " reads.\n";
}


sub readMultiLineDot {
    #$_ = <STDIN>;               #  Read and discard the tag
    my $save = $/;  $/ = ".\n";  #  Prepare to read the whole thing
    my $src = <STDIN>;           #  Read it
    $/ = $save;                  #  Reset our end of line marker
    $src =~ s/\s+//g;            #  Replace spaces and newlines
    chop $src;                   #  Get rid of the trailing .
    return($src);
}

print "{VER\n";
print "ver:2\n";
print "}\n";

while (!eof(STDIN)) {
    my $line = <STDIN>; chomp $line;

    if      ($line =~ m/^{BAT$/) {
        #  Discard BAT
        while ($line ne "}") {
            $line = <STDIN>; chomp $line;
        }
    } elsif ($line =~ m/^{ADT$/) {
        #  Discard ADT/ADL
        $line = <STDIN>; chomp $line;
        while ($line =~ m/^{ADL$/) {
            while ($line ne "}") {
                $line = <STDIN>; chomp $line;
            }
        }
        $line = <STDIN>; chomp $line;
        $line = <STDIN>; chomp $line;
    } elsif ($line =~ m/^{DST$/) {
        #  Convert a DST into a LIB.
        my $acc;
        my $mea;
        my $std;
        while ($line ne "}") {
            if ($line =~ m/^acc:(\S+)$/) {
                $acc = $1;
            }
            if ($line =~ m/^mea:(\d+\.*\d*)$/) {
                $mea = $1;
            }
            if ($line =~ m/^std:(\d+\.*\d*)$/) {
                $std = $1;
            }
            $line = <STDIN>; chomp $line;
        }
        print "{LIB\n";
        print "act:A\n";
        print "acc:$acc\n";
        print "ori:$mateStatus\n";
        print "mea:$mea\n";
        print "std:$std\n";
        print "src:\n";
        print "convert-v1-to-v2\n";
        print ".\n";
        print "nft:3\n";
        print "fea:\n";
        print "doNotTrustHomopolymerRuns=0\n";
        print "doNotOverlapTrim=$noOBT\n";
        print "isNotRandom=0\n";
        print ".\n";
        print "}\n";
        $lib = $acc;
    } elsif ($line =~ m/^{FRG$/) {
        my $acc;
        my $src;
        my $seq;
        my $qlt;
        my $clr;
        while ($line ne "}") {

            if ($line =~ m/^acc:(\S+)$/) {
                $acc = $1;
            }
            if ($line =~ m/^clr:(\d+,\d+)$/) {
                $clr = $1;
            }

            if ($line =~ m/^src:$/) {
                $src = readMultiLineDot();
            }
            if ($line =~ m/^seq:$/) {
                $seq = readMultiLineDot();
            }
            if ($line =~ m/^qlt:$/) {
                $qlt = readMultiLineDot();
            }

            $line = <STDIN>; chomp $line;
        }

        #  A bit of a hack, but needed for our existing shreds; if no
        #  library has been printed, print one.
        if (!defined($lib)) {
            print "{LIB\n";
            print "act:A\n";
            print "acc:454\n";
            print "ori:U\n";
            print "mea:0.0\n";
            print "std:0.0\n";
            print "src:\n";
            print "convert-v1-to-v2; forced 454\n";
            print ".\n";
            print "nft:5\n";
            print "fea:\n";
            print "doNotTrustHomopolymerRuns=0\n";
            print "hpsIsFlowGram=0\n";
            print "hpsIsPeakSpacing=0\n";
            print "doNotOverlapTrim=1\n";
            print "isNotRandom=0\n";
            print ".\n";
            print "}\n";
            $lib = "454";
        }
        print "{FRG\n";
        print "act:A\n";
        print "acc:$acc\n";
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
        if (defined($clv{$acc})) {
            $clvFound++;
            print "clv:$clv{$acc}\n";
        } else {
            $clvNotFound++;
        }
        if (defined($clq{$acc})) {
            $clqFound++;
            print "clq:$clq{$acc}\n";
        } else {
            $clqNotFound++;
        }
        print "clr:$clr\n";
        print "}\n";
    } elsif ($line =~ m/^{LKG$/) {
        my $fg1;
        my $fg2;
        while ($line ne "}") {
            if ($line =~ m/^fg1:(\S+)$/) {
                $fg1 = $1;
            }
            if ($line =~ m/^fg2:(\S+)$/) {
                $fg2 = $1;
            }
            $line = <STDIN>; chomp $line;
        }
        print "{LKG\n";
        print "act:A\n";
        print "frg:$fg1\n";
        print "frg:$fg2\n";
        print "}\n";
    } else {
        die "Unsupported line '$_'\n";
    }
}

print "{VER\n";
print "ver:1\n";
print "}\n";

if (defined(%clv) && ($clvNotFound > 0)) {
    print STDERR "$vec: Updated $clvFound vector clear ranges ($clvNotFound NOT updated).\n";
}

#print STDERR "Updated $clqFound quality clear ranges ($clqNotFound NOT updated).\n";
