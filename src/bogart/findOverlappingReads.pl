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
 #    Brian P. Walenz beginning on 2016-AUG-05
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

use strict;

if (scalar(@ARGV) == 0) {
    die "usage: $0 assembly-prefix readID [tigStore-bogart-stage]\n";
}

my $prefix = shift @ARGV;
my $readID = shift @ARGV;
my $store  = shift @ARGV;

my $gkpStore = "$prefix.gkpStore";
my $ovlStore = "$prefix.ovlStore";
my $tigStore = "$prefix.tigStore";
my $tigVers  = 1;

$tigStore = "$prefix.003.buildUnitigs.tigStore"    if ($store eq "003");
$tigStore = "$prefix.004.placeContains.tigStore"   if ($store eq "004");
$tigStore = "$prefix.005.mergeOrphans.tigStore"    if ($store eq "005");
$tigStore = "$prefix.007.breakRepeats.tigStore"    if ($store eq "007");
$tigStore = "$prefix.009.generateOutputs.tigStore" if ($store eq "009");

$gkpStore = "../$gkpStore"  if (! -d $gkpStore);
$gkpStore = "../$gkpStore"  if (! -d $gkpStore);
$gkpStore = "../$gkpStore"  if (! -d $gkpStore);

$ovlStore = "../$ovlStore"  if (! -d $ovlStore);
$ovlStore = "../$ovlStore"  if (! -d $ovlStore);
$ovlStore = "../$ovlStore"  if (! -d $ovlStore);

$tigStore = "../$tigStore"  if (! -d $tigStore);
$tigStore = "../$tigStore"  if (! -d $tigStore);
$tigStore = "../$tigStore"  if (! -d $tigStore);

die "failed to find gkpStore $prefix.gkpStore"  if (! -d $gkpStore);
die "failed to find ovlStore $prefix.ovlStore"  if (! -d $ovlStore);
die "failed to find tigStore $prefix.tigStore"  if (! -d $tigStore);


my %readOvl;

my $nOvl = 0;
my $nTig = 0;

open(F, "ovStoreDump -G $gkpStore -O $ovlStore -p $readID |");
while (<F>) {
    chomp;

    #  For -d dumps
    if (m/^\s*\d+\s+(\d+)\s+/) {
        $nOvl++;
        $readOvl{$1} = $_;
    }

    #  For -p dumps
    if (m/^\s*(\d+)\s+A:\s+\d+\s+/) {
        $nOvl++;
        $readOvl{$1} = $_;
    }
}
close(F);


system("ovStoreDump -G $gkpStore -O $ovlStore -p $readID");
print "\n";


my $tig;
my $len;
my $num;

open(F, "tgStoreDump -G $gkpStore -T $tigStore $tigVers -layout |");
while (<F>) {
    chomp;

    $tig = $1  if (m/^tig\s+(\d+)$/);
    $len = $1  if (m/^len\s+(\d+)$/);
    $num = $1  if (m/^numChildren\s+(\d+)$/);

    if (m/^read\s+(\d+)\s+/) {
        my $r = $1;

        if      ($r == $readID) {
            print "tig $tig len $len -- $_\n";
        } elsif (exists($readOvl{$1})) {
            $nTig++;
            printf "tig %6d len %8d -- %s -- %s\n", $tig, $len, $_, $readOvl{$1};
        }
    }
}
close(F);

print "\n";
print STDERR "Found $nOvl overlaps in $ovlStore (with $gkpStore).\n";
print STDERR "Found $nTig placements in $tigStore.\n";
