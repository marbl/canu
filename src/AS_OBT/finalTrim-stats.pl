#!/usr/bin/perl

use strict;

#  Computes various stats on the final trim log.  Should fold into the program itself.

open(F, "< $ARGV[0]") or die "Failed to open final.trim.log '$ARGV[0]'.\n";

my $header = <F>;

my $numReadsIn;
my $numBasesIn;

my $numReadsOut;
my $numBasesOut;

my $trimL;
my $trimR;

my $noTrimL++;
my $noTrimR++;

my $bigL;
my $bigR;
my $bigE;
my $bigB;

my $noOvl++;
my $del++;

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

    $noOvl++    if ($_ =~ m/no\soverlaps/);
    $del++      if ($v[5] eq "DEL");
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
print "noOvl       $noOvl reads with no overlaps and no trimming done\n";
print "del         $del reads deleted\n";
