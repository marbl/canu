#!/usr/local/bin/perl

#                    Confidential -- Do Not Distribute
#   Copyright (c) 2002 PE Corporation (NY) through the Celera Genomics Group
#                           All Rights Reserved.

$| = 1;

use strict;

my $path = $ARGV[0];

$path = "." if ($path eq "");

open(F, "< $path/0-input/scaffolds-list");
my @scafList = <F>;
close(F);
chomp @scafList;


my $sysTimeS = 0;
my $usrTimeS = 0;
my $clkTimeS = 0;
my $sysTimeH = 0;
my $usrTimeH = 0;
my $clkTimeH = 0;
my $sysTimeM = 0;
my $usrTimeM = 0;
my $clkTimeM = 0;
foreach my $s (@scafList) {
    if (-e "$path/1-search/$s.stats") {
        open(F, "< $path/1-search/$s.stats");
        while (<F>) {
            if (m/^systemTime:\s+(\d+\.\d+)$/) {
                $sysTimeS += $1;
            }
            if (m/^userTime:\s+(\d+\.\d+)$/) {
                $usrTimeS += $1;
            }
            if (m/^total:\s+(\d+.\d+)$/) {
                $clkTimeS += $1;
            }
        }
        close(F);
    }
}
$clkTimeH = $clkTimeS / 3600;
$sysTimeH = $sysTimeS / 3600;
$usrTimeH = $usrTimeS / 3600;
$clkTimeM = int(($clkTimeH - int($clkTimeH)) * 60);
$sysTimeM = int(($sysTimeH - int($sysTimeH)) * 60);
$usrTimeM = int(($usrTimeH - int($usrTimeH)) * 60);
$clkTimeH = int($clkTimeH);
$sysTimeH = int($sysTimeH);
$usrTimeH = int($usrTimeH);
$clkTimeS = int($clkTimeS);
$sysTimeS = int($sysTimeS);
$usrTimeS = int($usrTimeS);
printf STDOUT "ESTmapper: search required %7d seconds wall   (%3d:%02d).\n", $clkTimeS, $clkTimeH, $clkTimeM;
printf STDOUT "ESTmapper: search required %7d seconds system (%3d:%02d).\n", $sysTimeS, $sysTimeH, $sysTimeM;
printf STDOUT "ESTmapper: search required %7d seconds user   (%3d:%02d).\n", $usrTimeS, $usrTimeH, $usrTimeM;

$sysTimeS = 0;
$usrTimeS = 0;
$clkTimeS = 0;
open(F, "< $path/3-polish/run-script");
while (!eof(F)) {
    my $idx = <F>;  chomp $idx;
    my $cmd = <F>;

    #  Fix for reading old formats
    #
    if ($idx =~ m/(\d+).touch$/) {
        $idx = $1;
    }

    if (-e "$path/3-polish/$idx.stats") {
        open(X, "< $path/3-polish/$idx.stats");
        while (<X>) {
            if (m/^clockTime:\s+(\d+\.\d+)$/) {
                $clkTimeS += $1;
            }
            if (m/^systemTime:\s+(\d+\.\d+)$/) {
                $sysTimeS += $1;
            }
            if (m/^userTime:\s+(\d+\.\d+)$/) {
                $usrTimeS += $1;
            }
        }
        close(X);
    }
}
close(F);
$clkTimeH = $clkTimeS / 3600;
$sysTimeH = $sysTimeS / 3600;
$usrTimeH = $usrTimeS / 3600;
$clkTimeM = int(($clkTimeH - int($clkTimeH)) * 60);
$sysTimeM = int(($sysTimeH - int($sysTimeH)) * 60);
$usrTimeM = int(($usrTimeH - int($usrTimeH)) * 60);
$clkTimeH = int($clkTimeH);
$sysTimeH = int($sysTimeH);
$usrTimeH = int($usrTimeH);
$clkTimeS = int($clkTimeS);
$sysTimeS = int($sysTimeS);
$usrTimeS = int($usrTimeS);
printf STDOUT "ESTmapper: polish required %7d seconds wall   (%3d:%02d).\n", $clkTimeS, $clkTimeH, $clkTimeM;
printf STDOUT "ESTmapper: polish required %7d seconds system (%3d:%02d).\n", $sysTimeS, $sysTimeH, $sysTimeM;
printf STDOUT "ESTmapper: polish required %7d seconds user   (%3d:%02d).\n", $usrTimeS, $usrTimeH, $usrTimeM;


