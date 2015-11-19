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

my $FILE = "test.004.buildUnitigs";
$FILE = shift @ARGV if (scalar(@ARGV) > 0);

if ((! -e "$FILE.tigStore") ||
    (! -e "$FILE.fasta")) {
    die "Missing tigStore or fasta.  Run build-fasta.pl\n";
}

########################################

print STDERR "Aligning unitigs to unitigs\n";
if (! -e "$FILE.allVall.sim4db" ) {
    my $cmd;

    $cmd .= "/work/kmer/FreeBSD-amd64/bin/snapper2 ";
    $cmd .= "  -queries $FILE.fasta ";
    $cmd .= "  -genomic $FILE.fasta ";
    $cmd .= "  -mersize 16 ";
    $cmd .= "  -minmatchidentity 94 ";
    $cmd .= "  -minmatchcoverage 94 ";
    $cmd .= "  -aligns ";
    $cmd .= "  -verbose > $FILE.allVall.sim4db";

    system($cmd);
}

print STDERR "Filtering alignments\n";
if (! -e "$FILE.bubbles.sim4db") {
    my $cmd;

    $cmd .= "/work/kmer/FreeBSD-amd64/bin/filterPolishes ";
    $cmd .= "  -selfhits ";
    $cmd .= "  -D ";
    $cmd .= " < $FILE.allVall.sim4db ";
    $cmd .= " > $FILE.bubbles.sim4db";

    system($cmd);
}

########################################

my %detected1;
my %detected2;

my $utg;
my $ref;
my $ide;
my $cov;

print STDERR "Reading bubble mapping\n";
open(F, "< $FILE.bubbles.sim4db") or die;
while (<F>) {
    chomp;

    if (m/^edef=utg(\d+)/) {
        $utg = $1;
    }
    if (m/^ddef=utg(\d+)/) {
        $ref = $1;
    }
    if (m/^\d+\[(\d+)-\d+-\d+\]\s\d+\[\d+-\d+\]\s+<(\d+)-\d+-(\d+)-\w+-\w+>$/) {
        $ide = $3;
        $cov = $2;
    }
    if (m/^(\d+)-(\d+)\s\((\d+)-(\d+)\)\s/) {
        $cov = 100 * $cov / ($2 - $1);
    }
    if (m/^sim4end$/) {
        die if (!defined($utg));
        die if (!defined($ref));

        if (($ide >= 96) && ($cov >= 98)) {
            $detected1{"$utg"}++;
            $detected2{"$utg-$ref"}++;
        }

        #if ($cov < 100) {
        #    print STDERR "LOW COV $utg - $ref -- $cov\n";
        #}

        undef $utg;
        undef $ref;
    }
}
close(F);

print STDERR "Found ", scalar(keys %detected2), " actual bubble instances, from ", scalar(keys %detected1), " unitigs.\n";

my $numUnique = 0;

foreach my $k (keys %detected1) {
    if ($detected1{$k} == 1) {
        #print STDERR "UNIQUE $k\n";
        $numUnique++;
    }
}

print STDERR "Found $numUnique uniquely placeable bubbles.\n";

#
#  Scan the log.  Discover if we successfully merged a bubble.  Report bubbles that we merged that
#  were not verified by mapping.  Remove bubbles we are successful on from the list of alignments
#  (so we can next report what we failed on).
#

my $numMerged = 0;
my $numMergedUnique = 0;
my $numMergedMultiple = 0;
my $numMergedExtra = 0;
my $numMergedExtraMultiple = 0;

open(F, "< test.005.bubblePopping.log") or die;
while (<F>) {
    chomp;

    if (m/merged\sbubble\sunitig\s(\d+)\sinto\sunitig\s(\d+)$/) {
        $numMerged++;

        if (!exists($detected1{"$1"})) {
            #print STDERR "EXTRA $_ (no mapping at all)\n";
            $numMergedExtra++;

        } elsif (!exists($detected2{"$1-$2"})) {
            #print STDERR "EXTRA $_ (no mapping for this pair out of $detected1{\"$1\"} alignments - possibly incorrectly placed!)\n";
            $numMergedExtraMultiple++;

        } elsif ($detected1{"$1"} == 1) {
            $numMergedUnique++;

        } elsif ($detected1{"$1"} > 1) {
            #print STDERR "MULTI $_ ($detected1{\"$1\"} alignments - possibly incorrectly placed!)\n";
            $numMergedMultiple++;

        } else {
            die;
        }

        delete $detected2{"$1-$2"};
        delete $detected1{"$1"};
    }

}
close(F);

print STDERR "Merged $numMerged bubbles:\n";
print STDERR "       $numMergedUnique were uniquely alignable.\n";
print STDERR "       $numMergedMultiple were multiply alignable.\n";
print STDERR "       $numMergedExtra were not alignable at all.\n";
print STDERR "       $numMergedExtraMultiple were not alignable at the spot it was merged.\n";

#
#  Scan the layouts.  For each un-popped bubble (BOG didn't pop, but it aligned to a unitig) count
#  the number of fragments in the unitig, and length of the tig.  Report that and the number of
#  places the bubble mapped.
#

my $numMergedMissed = 0;
my $numMergedMissedUnique = 0;

my $id = 0;  #  Unitig ID we're looking at
my $nf = 0;  #  Number of fragments
my $ln = 0;  #  Length of the unitig

open(F, "< $FILE.layout");
open(O, "> $FILE.bubblesMISSED");
while (<F>) {
    chomp;

    if (m/^unitig\s+(\d+)$/) {
        if ($nf > 0) {
            $numMergedMissed++       if ($detected1{$id} > 0);
            $numMergedMissedUnique++ if ($detected1{$id} == 1);
            print O "Missed merge unitig $id of length $ln with $nf frags -- mapped to $detected1{$id} places.\n";
        }

        $id = $1;
        $nf = 0;
        $ln = 0;
    }

    if (($detected1{"$id"} > 0) &&
        (m/^FRG.*position\s+(\d+)\s+(\d+)$/)) {
        $nf++;

        $ln = $1 if ($ln < $1);
        $ln = $2 if ($ln < $2);
    }
}

close(F);

if ($nf > 0) {
    $numMergedMissed++       if ($detected1{$id} > 0);
    $numMergedMissedUnique++ if ($detected1{$id} == 1);
    print O "Missed merge unitig $id of length $ln with $nf frags -- mapped to $detected1{$id} places.\n";
}
close(O);

print STDERR "Failed to merge $numMergedMissed unitigs.\n";
print STDERR "                $numMergedMissedUnique were uniquely placeable.\n";
