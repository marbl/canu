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

#  Convert BLASR default output to OVL overlaps:  aIID bIID [I|N] aHang bHang error error_corrected

my %IDmap;

#  BLASR is not reporting symmetric overlaps.  We need to keep track of which overlaps we have reported.

my %reported;


open(F, "< $ARGV[0]") or die "Failed to open '$ARGV[0]' for reading sequence names.\n";
while (<F>) {
    my ($uid, $iid, $name) = split '\s+', $_;

    $IDmap{$name} = $iid;
}
close(F);

while (<STDIN>) {
    my @v = split '\s+', $_;

    $_ = join '\t', $_;

    #  The first read seems to have a sub-read range appended.
    if ($v[0] =~ m/^(.*)\/\d+_\d+$/) {
        $v[0] = $1;
    }

    my $aiid;
    my $biid;

    #  Or it might have UID,IID.  If this is ALWAYS the case, we don't need the uid map, yay!
    if ($v[0] =~ m/^\d+,(\d+)$/) {
        $aiid = $1;
    } else {
        $aiid  = $IDmap{$v[0]};
    }

    if ($v[1] =~ m/^\d+,(\d+)$/) {
        $biid = $1;
    } else {
        $biid  = $IDmap{$v[1]};
    }

    my $ni    = ($v[8] == 0) ? "N" : "I";
    my $error = 100.0 - $v[3];
    my $ecorr = 100.0 - $v[3];

    if ((exists($reported{"$aiid-$biid"}) || (exists($reported{"$biid-$aiid"})))) {
        next;
    }

    $reported{"$aiid-$biid"} = 0;

    #  Argh!  Do NOT flip coords if reversed.
    #($v[9], $v[10]) = ($v[10], $v[9])  if ($v[8] == 1);

    die "No A iid found for '$v[0]'\n" if (!defined($aiid));
    die "No B iid found for '$v[1]'\n" if (!defined($biid));

    next if ($aiid == $biid);
    #next if ($aiid <  $biid);

    die "First read is flipped\n$_\n" if ($v[4] != 0);

    my $a1 = $v[5];  #  Amount unaligned on left of first
    my $b1 = $v[7] - $v[6];

    my $a2 = $v[9];
    my $b2 = $v[11] - $v[10];

    #print "'$v[7]' '$v[8]' '$v[9]' '$v[10]'\n";

    #print "a1 $a1 from [5] $v[5]\n";
    #print "b1 $b1 from [7]-[6] $v[7] - $v[6]\n";
    #print "a2 $a2 from [9] $v[9]\n";
    #print "b2 $b2 from [11]-[10] $v[11] - $v[10]\n";

    my $fl = $v[8];

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

    #  HACK!
    $error = 0.01;
    $ecorr = 0.01;

    #print "$aiid\t$biid\t$ni\t$ahang\t$bhang\t$error\t$ecorr\t$label\t$_";
    print "$aiid\t$biid\t$ni\t$ahang\t$bhang\t$error\t$ecorr\n";
}
