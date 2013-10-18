#!/usr/bin/perl

use strict;

#  Computes various stats on the final trim log.  Should fold into the program itself.

die "usage: $0 prefix.finalTrim.log\n"  if (! -e "$ARGV[0]");

open(F, "< $ARGV[0]") or die "Failed to open finalTrim.log '$ARGV[0]'.\n";

my $header = <F>;

my $numReadsIn = 0;
my $numBasesIn = 0;

my $numReadsOut = 0;
my $numBasesOut = 0;

my $trimL = 0;
my $trimR = 0;

my $noTrimL = 0;
my $noTrimR = 0;

my $bigL = 0;
my $bigR = 0;
my $bigE = 0;
my $bigB = 0;

my $noOvl      = 0;
my $noOvlBases = 0;

my $del      = 0;
my $delBases = 0;

while (<F>) {
    chomp;

    my @v = split '\s+', $_;

    $numReadsIn++;
    $numBasesIn += $v[2] - $v[1];

    $numReadsOut++                  if ($v[5] ne "DEL");
    $numBasesOut += $v[4] - $v[3]   if ($v[5] ne "DEL");

    die "trimL got larger '$_'\n" if ($v[3] < $v[1]);
    die "trimR got larger '$_'\n" if ($v[2] < $v[4]);

    my $origLen = $v[2] - $v[1];
    my $trimLen = $v[4] - $v[3];

    my $tl = $v[3] - $v[1];
    my $tr = $v[2] - $v[4];

    $trimL    += $tl;
    $trimR    += $tr;

    $bigL++  if ($tl >= $trimLen / 2);
    $bigR++  if ($tr >= $trimLen / 2);

    if (($tl >= $trimLen / 2) ||
        ($tr >= $trimLen / 2)) {
        $bigE++;
    }

    if (($tl >= $trimLen / 2) &&
        ($tr >= $trimLen / 2)) {
        $bigB++;
    }

    $noTrimL++  if ($v[3] == $v[1]);
    $noTrimR++  if ($v[4] == $v[2]);

    if ($_ =~ m/no\soverlaps/) {
        $noOvl++;
        $noOvlBases += $origLen;
    }
    if ($v[5] eq "DEL") {
        $del++;
        $delBases += $origLen;
    }
}
close(F);

print "numReads    $numReadsIn reads in input\n";
print "numBases    $numBasesIn bases in input\n";
print "\n";
print "trimL       $trimL bases trimmed from left; ", $numReadsIn - $noTrimL, " reads trimmed\n";
print "trimR       $trimR bases trimmed from right; ", $numReadsIn - $noTrimR, " reads trimmed\n";
print "noTrimL     $noTrimL reads with no trimming on left\n";
print "noTrimR     $noTrimR reads with no trimming on right\n";
print "\n";
print "bigL        $bigL reads trimmed more than trimmedLength/2 from left\n";
print "bigR        $bigR reads trimmed more than trimmedLength/2 from right\n";
print "bigE        $bigE reads trimmed more than trimmedLength/2 from either side\n";
print "bigB        $bigB reads trimmed more than trimmedLength/2 from both sides\n";
print "\n";
print "numReads    $numReadsOut reads in output\n";
print "numBases    $numBasesOut bases in output\n";
print "noOvl       $noOvl reads $noOvlBases bp with no overlaps and no trimming done\n";
print "del         $del reads $delBases bp deleted\n";
