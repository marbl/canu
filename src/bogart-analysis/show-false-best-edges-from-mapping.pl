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

my $asm = shift @ARGV;

#  Reads:
#    $asm//best.edges
#    filtered-overlaps.true.ova
#    filtered-overlaps.false.ova
#
#  Reports which best.edges are false

my %true;
my %false;

open(F, "< filtered-overlaps.true.ova") or die;
while (<F>) {
    s/^\s+//;
    s/\s+$//;

    my @v = split '\s+', $_;

    $true{"$v[0]-$v[1]"}++;
    $true{"$v[1]-$v[0]"}++;
}
close(F);


open(F, "< filtered-overlaps.false.ova") or die;
while (<F>) {
    s/^\s+//;
    s/\s+$//;

    my @v = split '\s+', $_;

    $false{"$v[0]-$v[1]"}++;
    $false{"$v[1]-$v[0]"}++;
}
close(F);


    my $true5  = 0;
    my $true3  = 0;
    my $false5 = 0;
    my $false3 = 0;
    my $novel5 = 0;
    my $novel3 = 0;

open(F, "< $asm/best.edges") or die;
while (<F>) {
    s/^\s+//;
    s/\s+$//;

    my @v = split '\s+', $_;

    my $p5 = "$v[0]-$v[2]";
    my $p3 = "$v[0]-$v[4]";

    if      (exists($true{$p5})) {
        $true5++;
    } elsif (exists($false{$p5})) {
        print STDERR "FALSE 5 $_\n";
        $false5++;
    } else {
        print STDERR "NOVEL 5 $_\n";
        $novel5++;
    }

    if      (exists($true{$p3})) {
        $true3++;
    } elsif (exists($false{$p3})) {
        print STDERR "FALSE 3 $_\n";
        $false3++;
    } else {
        print STDERR "NOVEL 3 $_\n";
        $novel3++;
    }
}
close(F);

print "true5   $true5\n";
print "true3   $true3\n";
print "false5  $false5\n";
print "false3  $false3\n";
print "novel5  $novel5\n";
print "novel3  $novel3\n";
