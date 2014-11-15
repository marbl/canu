#!/usr/bin/perl

use strict;

my $log1 = shift @ARGV;
my $log2 = shift @ARGV;

if (!defined($log1) || !defined($log2)) {
    die "usage: $0 run1.log run2.log\n";
}

open(L1, "< $log1") or die;
open(L2, "< $log1") or die;

my ($a1, $a2, $a3, $a4, $a5, $a6, $a7, $a8, $a9);
my ($b1, $b2, $b3, $b4, $b5, $b6, $b7, $b8, $b9);

while (!eof(L1) && !eof(L2)) {

  anotherA:
    do {
        $a1 = <L1>;  #  FINAL, or "Correct" lines
    } while ($a1 !~ m/^FINAL/);
    $a2 = <L1>;  #  ORI seq
    $a3 = <L1>;  #  COR seq
    $a4 = <L1>;  #  COR qlt
    $a5 = <L1>;  #  COVERAGE
    $a6 = <L1>;  #  CORRECTIONS
    $a7 = <L1>;  #  DISCONNECTION
    $a8 = <L1>;  #  ADAPTER
    $a9 = <L1>;  #  RESULT

    #if ($a1 =~ m/^ADAPTERSEARCH/) {
    #    goto anotherA;
    #}

  anotherB:
    do {
        $b1 = <L2>;  #  FINAL, or "Correct" lines
    } while ($b1 !~ m/^FINAL/);
    $b2 = <L2>;  #  ORI seq
    $b3 = <L2>;  #  COR seq
    $b4 = <L2>;  #  COR qlt
    $b5 = <L2>;  #  COVERAGE
    $b6 = <L2>;  #  CORRECTIONS
    $b7 = <L2>;  #  DISCONNECTION
    $b8 = <L2>;  #  ADAPTER
    $b9 = <L2>;  #  RESULT

    #if ($b1 =~ m/^ADAPTERSEARCH/) {
    #    goto anotherB;
    #}

    my ($aID, $aLen, $aBgn, $aEnd);
    my ($bID, $bLen, $bBgn, $bEnd);

    #  FINAL or ADAPTERSEARCH
    if ($a1 =~ m/^\w+\sread\s(\d+)\slen\s(\d+)\s\(trim\s(\d+)-(\d+)\)$/) {
        $aID  = $1;
        $aLen = $2;
        $aBgn = $3;
        $aEnd = $4;
    } else {
        die "Nope a1 $a1";
    }

    if ($b1 =~ m/^\w+\sread\s(\d+)\slen\s(\d+)\s\(trim\s(\d+)-(\d+)\)$/) {
        $bID  = $1;
        $bLen = $2;
        $bBgn = $3;
        $bEnd = $4;
    } else {
        die "Nope b1 $b1";
    }

    die "ID mismatch $aID $bID\n" if ($aID != $bID);

    if (($aBgn != $bBgn) || ($aEnd != $bEnd)) {
        print "$aID/$bID $aLen/$bLen $aBgn-$aEnd $bBgn-$bEnd\n";
    }

    if (($aID % 10000) == 0) {
        print STDERR "$aID\n";
    }
}

