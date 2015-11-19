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

my %readLen;

open(F, "gatekeeper -dumpfragments -tabular test.gkpStore |");
while (<F>) {
    my @v = split '\s+', $_;
    $readLen{$v[1]} = $v[9];
}
close(F);


my %readHash;
my %readConf;
my $lastRead = 0;
my $nOvl     = 0;


my $potSub = 0;
my $proSub = 0;

open(F, "overlapStore -d test.ovlStore |");
while (<F>) {
    s/^\s+//;
    s/\s+$//;

    my @v = split '\s+', $_;

    if ($lastRead != $v[0]) {
        my $mCount = scalar(keys %readConf);

        if ($mCount > 0) {
            print "$lastRead\t$mCount\n";
            $potSub++;
        }

        if ($mCount > 5) {
            $proSub++;
        }

        $lastRead = $v[0];
        $nOvl     = 0;
        undef %readHash;
        undef %readConf;
    }

    next  if (($v[3] < 0) && ($v[4] > 0));  #  read is contained in other

    $nOvl++;

    if (! exists($readHash{$v[1]})) {
        $readHash{$v[1]} = "$v[3] $v[4]";
        next;
    }

    my ($na, $nb) = ($v[3], $v[4]);
    my ($oa, $ob) = split '\s+', $readHash{$v[1]};

    my $len = $readLen{$v[0]};

    my $nbgn = ($na > 0) ?         $na  : 0;
    my $nend = ($nb < 0) ? ($len + $nb) : $len;

    my $obgn = ($oa > 0) ?         $oa  : 0;
    my $oend = ($ob < 0) ? ($len + $ob) : $len;

    #print "$_ -- $nbgn $nend - $obgn $oend\n";

    goto saveConf   if (($nbgn <= $obgn) && ($nend <= $obgn));  # n before o
    goto saveConf   if (($nbgn >= $oend) && ($nend >= $oend));  # n after o

    my $min = ($nbgn < $obgn) ? $obgn : $nbgn;  #  Largest bgn coord
    my $max = ($nend < $oend) ? $nend : $oend;  #  Smallest end coord

    my $nlen = $nend - $nbgn;
    my $olen = $oend - $obgn;

    my $olap = $max - $min;

    #print "$_ -- $olap $nlen $olen\n";

    next   if (($olap / $nlen > 0.5) || ($olap / $olen > 0.5));  #  reads overlap by lots

    #  different regions, probably sub read detected
  saveConf:
    $readConf{$v[1]}++;

    #print "$_\n";
}
close(F);

print STDERR "potential unsplit subreads:  $potSub\n";
print STDERR "probable  unsplit subreads:  $proSub\n";
