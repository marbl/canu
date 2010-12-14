#!/usr/bin/perl

use strict;

my $prefix = "test.004.buildUnitigs";
my $IDEAL  = "../../ideal/porphyromonas_gingivalis_w83.flx.fragment.E9T0MN001.ideal.intervals";

$prefix = $ARGV[0]  if (scalar(@ARGV) > 0);
$IDEAL  = $ARGV[1]  if (scalar(@ARGV) > 1);

if (! -e "$prefix.bogus.out") {
    my $cmd;

    $cmd  = "bogusness \\\n";
    $cmd .= " -nucmer $prefix.coords \\\n";
    $cmd .= " -ideal  $IDEAL \\\n";
    $cmd .= " > $prefix.bogusness.out";

    system($cmd);
}

my $lastUtg;

open(F, "sort -k10n -k18n < $prefix.bogusness.out |") or die "Failed to read sort output.\n";
open(O, "> $prefix.bogusness.wiki") or die "Failed to open '$prefix.bogusness.wiki' for writing.\n";

print O "{| class=\"wikitable\" border=1\n";
print O "! unitig !! align num !! utg coords !! gen coords !! status !! ideal type !! ideal index !! ideal coords !! length !! utg cov !! ideal cov !! annotation\n";
print O "|-\n";

while (<F>) {
    chomp;

    if (m/^\|\s(utg\d+)\s/) {
        if (!defined($lastUtg)) {
            $lastUtg = $1;
        }
        if ($lastUtg ne $1) {
            print O "| colspan=12 bgcolor=#666666 |\n";
            print O "|-\n";
            $lastUtg = $1;
        }

        s!BEGINSin  \|\| UNIQ!bgcolor="FireBrick" \| BEGINSin  \|\| bgcolor="FireBrick" \| UNIQ!;
        s!ENDSin    \|\| UNIQ!bgcolor="FireBrick" \| ENDSin    \|\| bgcolor="FireBrick" \| UNIQ!;

        s!BEGINSin  \|\| REPT!bgcolor="ForestGreen" \| BEGINSin  \|\| bgcolor="ForestGreen" \| REPT!;
        s!ENDSin    \|\| REPT!bgcolor="ForestGreen" \| ENDSin    \|\| bgcolor="ForestGreen" \| REPT!;

        s!CONTAINED \|\| UNIQ!bgcolor="Indigo" \| CONTAINED \|\| bgcolor="Indigo" \| UNIQ!;

        s!CONTAINED \|\| REPT!bgcolor="SteelBlue" \| CONTAINED \|\| bgcolor="SteelBlue" \| REPT!;

        print O "$_\n";
        print O "|-\n";
    }
}

print O "|}\n";

close(O);
close(F);
