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

#  Simple hack to merge multiple QC reports into one.  Uses the first
#  QC on the command line to figure out what fields exist, then writes
#  a multiple column QC using those fields and data from all QC files
#  on the command line.

my @labels;
my %values;

my $s;
my $l;
my $v;

my $firstFile = 1;

while (scalar(@ARGV) > 0) {
    if (open(F, "< $ARGV[0]")) {
        if ($ARGV[0] =~ m!.*/NIGHTLY/(\d\d\d\d-\d\d-\d\d-\d\d\d\d)/(.*)/9-terminator/.*.qc$!) {
            print "$1/$2\n";
        } else {
            print "$ARGV[0]\n";
        }

        while (<F>) {
            $_ =~ s/^\s+//;
            $_ =~ s/\s+$//;

            if ($_ =~ m/^\s+$/) {
                next;
            }

            if ($_ =~ m/^\[(.*)\]$/) {
                $s = $1;
                next;
            }

            if ($_ =~ m/^(.*)=(.*)$/) {
                $l = "$s\0$1";
                $v = $2;
            } else {
                next;
            }

            if ($firstFile) {
                push @labels, $l;
            }

            $values{$l} .= substr(" $v                ", 0, 16);
            $values{$l} .= "BRIWASHERE";
        }
        close(F);

        my @k = keys %values;
        foreach my $l (@k) {
            if ($values{$l} =~ m/^(.*)BRIWASHERE$/) {
                $values{$l}  = $1;
            } else {
                $values{$l} .= substr(" N/A                ", 0, 16);
            }
        }

        $firstFile = 0;
    }

    shift @ARGV;
}


my $lasts;
foreach my $xx (@labels) {
    ($s, $l) = split '\0', $xx;

    if ($s ne $lasts) {
        print "\n[$s]\n";
        $lasts = $s;
    }

    $l = substr("$l                    ", 0, 20);

    my $d = "    ";

    my @v = split '\s+', $values{$xx};
    foreach my $a (@v) {
        foreach my $b (@v) {
            #  Because $xx is whitespace justified, the first
            #  thing from the split can be empty.
            if (($a ne "") && ($b ne "") && ($a ne $b)) {
                $d = " ** ";
            }
        }
    }

    print "$l$d$values{$xx}\n";
}
