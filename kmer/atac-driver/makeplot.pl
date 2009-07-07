#!/usr/bin/env perl
#
# This file is part of A2Amapper.
# Copyright (c) 2008-2009 J. Craig Venter Institute
# Author: Brian Walenz
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

use strict;
use FindBin;

my $mm = shift @ARGV;
my $in = shift @ARGV;
my $ot = shift @ARGV;

if (($mm ne "u") && ($mm ne "r")) {
    die "First arg must be 'u' (ungapped matches) or 'r' (runs).\n";
}

if (!defined($ot) && ($in =~ m/(.*).atac/)) {
    $ot = $1;
}
if ($ot =~ m/^(.*).png/) {
    $ot = $1;
}
if ($ot =~ m/^(.*).ps/) {
    $ot = $1;
}

die if (!defined($in));
die if (!defined($ot));

my $version = `gnuplot -V`;

if ($version =~ m/gnuplot\s+(\d+\.\d+)\s+/) {
    $version = $1;
} else {
    chomp $version;
    print STDERR "WARNING:  Unknown gnuplot version '$version'\n";
    $version = 0;
}

if ($version < 4.2) {
    print STDERR "gnuplot version 4.2 is needed for plots.\n";
    exit(0);
}

open(FD, "> $ot.fdat");
open(RD, "> $ot.rdat");
open(GP, "> $ot.gp");

print GP "set size 1,1\n";
print GP "set grid\n";
print GP "unset key\n";
print GP "set border 10\n";
print GP "set tics scale 0\n";
print GP "set xlabel \"REF\"\n";
print GP "set ylabel \"ASM\"\n";
print GP "set title \"$ot\"\n";
print GP "set format \"%.0f\"\n";
#print GP "set mouse format \"%.0f\"\n";
#print GP "set mouse mouseformat \"[%.0f, %.0f]\"\n";
#print GP "set mouse clipboardformat \"[%.0f, %.0f]\"\n";
print GP "set style line 1  lt 1 lw 1 pt 6 ps 1\n";
print GP "set style line 2  lt 3 lw 1 pt 6 ps 1\n";
print GP "set style line 3  lt 2 lw 1 pt 6 ps 1\n";

#  We need to know the length of the reference so we can cycle
#  the coordinates
#
my $refLength = 0;

#  And we need to know the lengths of the scaffolds in the assembly.
#
my $asmFile1;
my $asmId1;
my $asmFile2;
my $asmId2;

open(IN, "< $in") or die;
while (<IN>) {
    if (m/assemblyId1=(.*)/) {
        $asmId1 = $1;
    }
    if (m/assemblyFile1=(.*)/) {
        $asmFile1 = $1;
    }
    if (m/assemblyFile2=(.*)/) {
        $asmFile2 = $1;
    }
    if (m/assemblyId2=(.*)/) {
        $asmId2 = $1;
    }
}
close(IN);


#  Figure out which scaffolds are reversed.

my %reversed;
my %goofed;

{
    open(IN, "< $in") or die;
    while (<IN>) {
        my @v = split '\s+', $_;
        if (($v[0] eq "M") && ($v[1] eq $mm)) {
            $reversed{$v[8]} += $v[10] * $v[11];  #  length * orientation

            #  Remember if there was a goofy inversion in this
            #  scaffold.  First time through we remember the
            #  orientation of the first match, then later matches
            #  check it is the same, reseting to a token value (2) if
            #  they differ.

            if (!defined($goofed{$v[8]})) {
                $goofed{$v[8]} = $v[11];
            } else {
                if ($goofed{$v[8]} != $v[11]) {
                    $goofed{$v[8]} = 2;
                }
            }
        }
    }
    close(IN);
    foreach my $k (keys %reversed) {
        if ($reversed{$k} < 0) {
            $reversed{$k} = 1;
        } else {
            $reversed{$k} = 0;
        }
        if ($goofed{$k} == 2) {
            $goofed{$k} = 1;
        } else {
            $goofed{$k} = 0;
        }
    }
}



#  Find the reference length

my $refLength;

open(F, "< $asmFile1") or die "Failed to open genome reference '$asmFile1'\n";
$_ = <F>;  #  defline
while (<F>) {
    s/^\s+//;
    s/\s+$//;
    $refLength += length ($_);
}
close(F);


#  Find the assembly lengths -- ignore anything without a match.

my %asmLength;
my $asmLength;

open(F, "$FindBin::Bin/leaff -F $asmFile2 -i $asmId2 |");
while (<F>) {
    my @v = split '\s+', $_;
    if ($v[0] eq "G") {
        $asmLength{$v[2]} = $v[3] if (exists($reversed{$v[2]}));
        $asmLength       += $v[3];
    }
}
close(F);

print GP "set xrange [0:", $refLength+1, "]\n";
print GP "set yrange [0:", $asmLength+1, "]\n";


#  Figure out where to split the reference - we shift both X and Y, so
#  we arbitrarily pick any assembly sequence and anchor it at the
#  origin.  Well, not arbitrary.  We pick the longest only so that
#  we're pretty sure we didn't pick some crappy tiny contig.

my $refSplit  = 0;

{
    my $asmAnchorSequence;
    my $minAsm = 999999999;

    #  Pick the longest sequence as the anchor -- skipping this block
    #  will instead pick the first thing in the file.
    #
    #  Don't pick anything with reversed crap on the ends; this
    #  greatly screws up.
    #
    foreach my $k (keys %asmLength) {
        if ($goofed{$k}) {
            next;
        }
        if (!defined($asmAnchorSequence) || ($asmLength{$asmAnchorSequence} < $asmLength{$k})) {
            $asmAnchorSequence = $k;
        }
    }

    #  Except when all scaffolds are goofed.  *sigh*

    if (!defined($asmAnchorSequence)) {
        foreach my $k (keys %asmLength) {
            if (!defined($asmAnchorSequence) || ($asmLength{$asmAnchorSequence} < $asmLength{$k})) {
                $asmAnchorSequence = $k;
            }
        }
    }

    open(IN, "< $in") or die;
    while (<IN>) {
        my @v = split '\s+', $_;

        if (($v[0] eq "M") && ($v[1] eq $mm)) {
            if ($reversed{$v[8]}) {
                $v[9]   = $asmLength{$v[8]} - ($v[9] + $v[10]);
                $v[11] *= -1;
            }

            if (!defined($asmAnchorSequence)) {
                $asmAnchorSequence = $v[8];
            }

            if ($v[8] eq $asmAnchorSequence) {
                if ($v[9] < $minAsm) {
                    $minAsm   = $v[9];
                    $refSplit = $v[5];
                }
            }
        }
    }
    close(IN);
}

#  Figure out how to place the assembly sequences
#
#  $refSplit  controls where we shift the reference origin.
#  %offsetRef controls where we place the assembly in the Y axis.
#
#  We want to rotate the reference so the largest scaffold is placed
#  at the origin.

my %offsetRef;

{
    my %minY;
    my %minYbackup;

    open(IN, "< $in") or die;
    while (<IN>) {
        my @v = split '\s+', $_;
        if (($v[0] eq "M") && ($v[1] eq $mm)) {
            if ($reversed{$v[8]}) {
                $v[9]   = $asmLength{$v[8]} - ($v[9] + $v[10]);
                $v[11] *= -1;
            }

            #  Confusing.  Rotate the reference coordinate, then
            #  remember the smallest for each scaffold.

            $v[5] -= $refSplit;
            $v[5] += $refLength if ($v[5] < 0);

            my $d = $v[5];

            #  Ignore if this is a tiny crappy little thing.  This
            #  allows us to place most of the real matches on the
            #  diagonal, showing obvious small chimers.

            if ($v[10] >= 2000) {
                if (!exists($minY{$v[8]}) || ($d < $minY{$v[8]})) {
                    $minY{$v[8]} = $d;
                }
            } else {
                if (!exists($minY{$v[8]}) || ($d < $minYbackup{$v[8]})) {
                    $minYbackup{$v[8]} = $d;
                }
            }
        }
    }
    close(IN);

    #  If we never found a large block, use the biggest we did find.

    foreach my $k (keys %minYbackup) {
        if (!exists($minY{$k})) {
            $minY{$k} = $minYbackup{$k};
        }
    }

    my @sortme;
    my $lengthsum = 0;
    foreach my $k (keys %minY) {
        push @sortme, "$minY{$k}\0$k";
    }
    @sortme = sort { $a <=> $b } @sortme;
    foreach my $v (@sortme) {
        my ($p, $k) = split '\0', $v;
        $offsetRef{$k} = $lengthsum;
        $lengthsum += $asmLength{$k};
    }
}
    
print GP "set ytics ( \\\n";

{
    #  Ugh, gross.
    my @keys = keys %offsetRef;
    while (scalar(@keys)) {
        my $k = shift @keys;

        if (scalar(@keys) > 0) {
            print GP "  \"$k\" $offsetRef{$k},\\\n";
        } else {
            print GP "  \"$k\" $offsetRef{$k}\\\n";
        }
    }
}
print GP ")\n";

print GP "set xtics ( \\\n";
#for (my $p=500000; $p<$refLength; $p += 500000) {
#    my $i = $p - $refSplit;
#    if ($i < 0) {
#        $i += $refLength;
#    }
#    print GP "  \"$p\" $i,\\\n";
#}
print GP "  \"origin\" ", $refLength - $refSplit, "\\\n";
print GP ")\n";




my $hasFdat = 0;
my $hasRdat = 0;

open(IN, "< $in") or die;
while (<IN>) {
    my @v = split '\s+', $_;

    if (($v[0] eq "M") && ($v[1] eq $mm)) {
        if ($reversed{$v[8]}) {
            $v[9]   = $asmLength{$v[8]} - ($v[9] + $v[10]);
            $v[11] *= -1;
        }

        my $abeg = $v[5];
        my $aend = $v[5] + $v[6];
        my $bbeg = $v[9];
        my $bend = $v[9] + $v[10];

        $abeg -= $refSplit;
        $aend -= $refSplit;

        if (($abeg < 0) || ($aend < 0)) {
            $abeg += $refLength;
            $aend += $refLength;
        }

        $bbeg += $offsetRef{$v[8]};
        $bend += $offsetRef{$v[8]};

        if ($v[11] == 1) {
            $hasFdat++;
            print FD "$abeg $bbeg\n";
            print FD "$aend $bend\n";
            print FD "\n\n";
        } else {
            $hasRdat++;
            print RD "$abeg $bend\n";
            print RD "$aend $bbeg\n";
            print RD "\n\n";
        }
    }
}
close(IN);

close(FD);
close(RD);

print GP "set terminal png tiny size 800,800\n";
print GP "set output \"$ot.png\"\n";

if ($hasFdat && $hasRdat) {
    print GP "plot \\\n";
    print GP "  \"$ot.fdat\" w lp ls 1, \\\n";
    print GP "  \"$ot.rdat\" w lp ls 2\n";
    #print GP "pause -1\n";
} elsif ($hasFdat) {
    print GP "plot \\\n";
    print GP "  \"$ot.fdat\" w lp ls 1\n";
    #print GP "pause -1\n";
} elsif ($hasRdat) {
    print GP "plot \\\n";
    print GP "  \"$ot.rdat\" w lp ls 2\n";
    #print GP "pause -1\n";
} else {
    #  No matches??
    #die;
}

print GP "set terminal postscript color\n";
print GP "set output \"$ot.ps\"\n";
print GP "replot\n";

close(GP);

system("gnuplot $ot.gp");

#unlink "$ot.fdat";
#unlink "$ot.rdat";
#unlink "$ot.gp";
