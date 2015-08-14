#!/usr/bin/perl

use strict;

#  Given:
#    an assembly of true inferred overlaps
#    an assembly of computed overlaps
#  report the number of true overlaps that are missed in the computed overlaps

my %truthOverlaps;
my @readLengths;


open(F, "gatekeeper -dumpfragments -tabular JUNKTEST3/BL/test.gkpStore |");
while (<F>) {
    my @v = split '\s+', $_;

    $readLengths[$v[1]] = $v[9];
}
close(F);

open(F, "gatekeeper -dumpfragments -tabular JUNKTEST3/CAordered/test.gkpStore |");
while (<F>) {
    my @v = split '\s+', $_;

    die "IID $v[1] BL $readLengths[$v[1]] != CA $v[9]\n"   if ($readLengths[$v[1]] != $v[9]);
}
close(F);




open(F, "overlapStore -d JUNKTEST3/BL/test.ovlStore |");
while (<F>) {
    s/^\s+//;
    s/\s+$//;

    my ($aIID, $bIID, $orient, $aHang, $bHang, $origID, $corrID) = split '\s+', $_;

    my $aLen = $readLengths[$aIID];
    my $bLen = $readLengths[$bIID];

    my $aOvl = 0;
    my $bOvl = 0;
    my $oLen = 0;

    #  Swiped from bogart

    if ($aHang < 0) {
      #  bHang < 0      ?     ----------  :     ----
      #                 ?  ----------     :  ----------
      #
      $aOvl = ($bHang < 0) ? ($aLen + $bHang) : ($aLen);
      $bOvl = ($bHang < 0) ? ($bLen + $aHang) : ($bLen + $aHang - $bHang);
    } else {
      #  bHang < 0      ?  ----------              :  ----------
      #                 ?     ----                 :     ----------
      #
      $aOvl = ($bHang < 0) ? ($aLen - $aHang + $bHang) : ($aLen - $aHang);
      $bOvl = ($bHang < 0) ? ($bLen)                   : ($bLen - $bHang);
    }

    $oLen = ($aOvl + $bOvl) / 2;

    $truthOverlaps{"$aIID-$bIID"}++;
}
close(F);

print STDERR "truthOverlaps:        ", scalar(keys %truthOverlaps), "\n";

open(F, "overlapStore -d JUNKTEST3/CAordered/test.ovlStore |");
while (<F>) {
    s/^\s+//;
    s/\s+$//;

    #  v[0] = iid
    #  v[1] = iid
    #  v[2] = orient
    #  v[3] = hang
    #  v[4] = hang
    #  v[5] = ident
    #  v[6] = ident

    my ($aIID, $bIID, $orient, $aHang, $bHang, $origID, $corrID) = split '\s+', $_;


    delete $truthOverlaps{"$aIID-$bIID"};
}
close(F);

print STDERR "truthOverlaps missed: ", scalar(keys %truthOverlaps), "\n";
