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

#  cd /work/blasr-overlapper/JUNKTEST3/CAordered
#  time sh /work/scripts/blasr-overlaps.sh BLASROVERLAP.fasta BLASROVERLAP.fasta BLASROVERLAP
#7043.059u 16.347s 22:09.56 530.9%       6573+816k 11+192523io 11881pf+0w

#
#  NMAME to IID MAP
#

my %NAMEtoIID;

open(F, "< CAordered/test.gkpStore.fastqUIDmap") or die "Failed to open 'CAordered/test.gkpStore.fastqUIDmap'\n";
while (<F>) {
    my @v = split '\s+', $_;

    #$IIDtoNAME{$v[1]} = $v[2];
    $NAMEtoIID{$v[2]} = $v[1];
}
close(F);

#
#  READ LENGTHS
#

my @readLength;
my %contained;

open(F, "gatekeeper -dumpfragments -tabular BL/test.gkpStore |");
while (<F>) {
    my @v = split '\s+', $_;

    $readLength[$v[1]] = $v[9];
}
close(F);

open(F, "gatekeeper -dumpfragments -tabular CAordered/test.gkpStore |");
while (<F>) {
    my @v = split '\s+', $_;

    die "IID $v[1] BL $readLength[$v[1]] != CA $v[9]\n"   if ($readLength[$v[1]] != $v[9]);
}
close(F);

#
#
#

my %reported;


open(F, "< BO/PAIRS.blasr.sam.coords") or die "Failed to open 'BO/PAIRS.blasr.sam.coords'\n";
open(O, "> BO/PAIRS.blasr.sam.ova");

$_ = <F>;
$_ = <F>;
$_ = <F>;
$_ = <F>;
while (<F>) {
    my @v = split '\s+', $_;

    my $aIID = $NAMEtoIID{$v[9]};
    my $bIID = $NAMEtoIID{$v[10]};

    next   if (!defined($aIID));  #  Extra overlaps are not an error.
    next   if (!defined($bIID));

    die "undef for a $v[9]\n"   if (!defined($aIID));  #  Extra overlaps are an error.
    die "undef for b $v[10]\n"  if (!defined($bIID));

    next   if ($aIID == $bIID);

    next   if (exists($reported{"$aIID-$bIID"}));

    $reported{"$aIID-$bIID"}++;
    $reported{"$bIID-$aIID"}++;


    my $l1 = $readLength[$aIID];
    my $a1 = $v[0] - 1;       #  Amount unaligned on left of first
    my $b1 = $l1 - $v[1];     #  Amount unaligned on right of first

    my $ori;

    my $l2 = $readLength[$bIID];
    my $a2;
    my $b2;

    if ($v[2] < $v[3]) {
        $ori = "N";
        $a2  = $v[2] - 1;
        $b2  = $l2 - $v[3];
    } else {
        $ori = "I";
        $a2  = $v[3] - 1;
        $b2  = $l2 - $v[2];
    }

    #  Extend near global to be global.

    my $maxTol = 10;

    $a1 = 0    if ($a1 < $maxTol);
    $b1 = $l1  if ($l1 < $b1 + $maxTol);

    $a2 = 0    if ($a2 < $maxTol);
    $b2 = $l2  if ($l2 < $b2 + $maxTol);

    my $fl = ($ori eq "I");

    my $ahang = 0;
    my $bhang = 0;
    my $label = "";

    #  Handle A contained in B
    if    (($a1 == 0) && ($b1 == 0)) {
        if ($fl == 0) {
            $ahang = -$a2;
            $bhang =  $b2;
            $label = "AinBf";
        } else {
            $ahang = -$b2;
            $bhang =  $a2;
            $label = "AinBr";
        }
    }

    #  Handle B contained in A
    elsif (($a2 == 0) && ($b2 == 0)) {
        if ($fl == 0) {
            $ahang =  $a1;
            $bhang = -$b1;
            $label = "BinAf";
        } else {
            $ahang =  $b1;
            $bhang = -$a1;
            $label = "BinAr";
        }
    }

    #  Handle a dovetail off the left end of A
    elsif (($a1 == 0) && ($b2 == 0) && ($fl == 0)) {
        $ahang = -$a2;
        $bhang = -$b1;
        $label = "BdoveAf";
    }
    elsif (($a1 == 0) && ($a2 == 0) && ($fl == 1)) {
        $ahang = -$b2;
        $bhang = -$b1;
        $label = "BdoveAr";
    }

    #  Handle dovetail off the right end of A
    elsif (($b1 == 0) && ($a2 == 0) && ($fl == 0)) {
        $ahang =  $a1;
        $bhang =  $b2;
        $label = "AdoveBf";
    }
    elsif (($b1 == 0) && ($b2 == 0) && ($fl == 1)) {
        $ahang =  $a1;
        $bhang =  $a2;
        $label = "AdoveBr";
    }

    #  All the rest aren't valid overlaps.
    else {
        next;
        $label = "INVALID";
    }

    my $error = 100 - $v[8];
    my $ecorr = 100 - $v[8];

    #print "$aIID\t$bIID\t$ni\t$ahang\t$bhang\t$error\t$ecorr\t$label\t$_";
    print O "$aIID\t$bIID\t$ori\t$ahang\t$bhang\t$error\t$ecorr\n";
}
close(F);
close(O);

system("convertOverlap -ovl < BO/PAIRS.blasr.sam.ova  > BO/PAIRS.blasr.sam.ovb");
system("overlapStoreBuild -g CAordered/test.gkpStore -o blasr-pairs.ovlStore -F 1 BO/PAIRS.blasr.sam.ovb");

