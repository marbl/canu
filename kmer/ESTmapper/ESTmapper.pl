#!/usr/bin/perl

#  Copyright (c) 2002 PE Corporation (NY) through the Celera Genomics Group
#  Copyright (c) 2003, 2004 Applied Biosystems
#  Copyright (c) 2004, 2005, 2006 Brian Walenz

$| = 1;

# Perl version 5.005_03 is too old, it requires two args to mkdir.

use strict;

use vars "%prog", "%args";

use FindBin;
use lib "$FindBin::Bin/util";
use scheduler;
use fasta;

require "run.pl";
require "1-configure.pl";
require "2-search.pl";
require "3-filter.pl";
require "4-polish.pl";
require "5-assemble.pl";

setExecutables();
parseArgs(@ARGV);

if      ($args{'runstyle'} eq "est") {
    configure();
    search();
    filter();
    polish();
    assembleOutput();
} elsif ($args{'runstyle'} eq "mrna") {
    $args{'relink'} = 1000;
    $args{'abort'}  = 1;

    configure();
    search();
    filter();
    polish();
    assembleOutput();
} elsif ($args{'runstyle'} eq "snp") {
    $args{'minidentity'} = 95;
    $args{'mincoverage'} = 80;

    configure();
    search();
    filter();
    polish();
    assembleOutput();
    parseSNP();
} elsif ($args{'runstyle'} eq "-time") {
    summarizeTime();
    exit(0);
} else {
    print STDERR "Basic help N/A.\n";
}

print STDERR "ESTmapper: script finished everything in ", time() - $args{'startTime'}, " wall-clock seconds.\n" if (time() != $args{'startTime'});


if (-e $args{'runInforFile'}) {
    my $time = time();

    open(F, ">> $args{'runInforFile'}");
    print F "endTime:  $time (", scalar(localtime($time)), ")\n";
    close(F);
}

exit(0);

















sub parseSNP {
    #  Parse the SNPs out
    #
    if (! -e "$args{'path'}/snps-parsed") {
        print STDERR "ESTmapper--  Parsing the SNPs\n";

        #  Sort, if needed.
        #
        if (! -e "$args{'path'}/polishes-good.sorted") {
            print STDERR "ESTmapper--  Sorting polishes by sequence ID; using 2GB memory maximum.\n";
            if (runCommand("$prog{'sortPolishes'} -m 2000 -c < $args{'path'}/polishes-good > $args{'path'}/polishes-good.sorted")) {
                unlink "$args{'path'}/polishes-good.sorted";
                die "Failed to sort the polishes.\n";
            }
        }

        #  Parse the options, looking for SNP specific ones
        #
        my @ARGS = @ARGV;
        my $snpdelimiter = "";
        my $snpsizetag   = "";
        my $snppostag    = "";
        my $snpoffset    = "";
        my $snpoutformat = "";

        while (scalar @ARGS > 0) {
            my $arg = shift @ARGS;

            if ($arg eq "-snpdelimiter") {
                $arg = shift @ARGS;
                $snpdelimiter = "-d \"$arg\"";
            } elsif ($arg eq "-snpsizetag") {
                $arg = shift @ARGS;
                $snpsizetag = "-s \"$arg\"";
            } elsif ($arg eq "-snppostag") {
                $arg = shift @ARGS;
                $snppostag = "-p \"$arg\"";
            } elsif ($arg eq "-snpoffset") {
                $arg = shift @ARGS;
                $snpoffset = "-o $arg";
            } elsif ($arg eq "-snpoutformat") {
                $arg = shift @ARGS;
                $snpoutformat = "-format $arg";
            }
        }

        #  PARSE!
        #
        if (runCommand("$prog{'parseSNPs'} $snpdelimiter $snpsizetag $snppostag $snpoffset $snpoutformat -F $args{'path'}/snps-failed -O $args{'path'}/snps-parsed < $args{'path'}/polishes-good.sorted > $args{'path'}/summary-snps")) {
            unlink "$args{'path'}/snps-failed";
            unlink "$args{'path'}/snps-parsed";
            unlink "$args{'path'}/summary-snps";
            die "Failed to parse SNP locations from polishes.\n";
        }
    }
}




sub sTOhms ($) {
    my ($s, $m, $h) = @_;
    $h = $s / 3600;
    $m = int(($h - int($h)) * 60);
    $h = int($h);
    $s = int($s);
    return($h,$m,$s);
}

sub summarizeTime {
    my ($sysTimeS, $usrTimeS, $clkTimeS);
    my ($sysTimeH, $usrTimeH, $clkTimeH);
    my ($sysTimeM, $usrTimeM, $clkTimeM);

    $sysTimeS = $usrTimeS = $clkTimeS = 0;
    open(F, "find $args{'path'}/1-search -name *.stats -print |");
    while (<F>) {
        chomp;
        open(G, "< $_");
        while (<G>) {
            $sysTimeS += $1 if (m/^systemTime:\s+(\d+\.\d+)$/);
            $usrTimeS += $1 if (m/^userTime:\s+(\d+\.\d+)$/);
            $clkTimeS += $1 if (m/^total:\s+(\d+.\d+)$/);
        }
        close(G);
    }
    close(F);

    ($clkTimeH,$clkTimeM,$clkTimeS) = sTOhms($clkTimeS);
    ($sysTimeH,$sysTimeM,$sysTimeS) = sTOhms($sysTimeS);
    ($usrTimeH,$usrTimeM,$usrTimeS) = sTOhms($usrTimeS);
    printf STDOUT "ESTmapper: search required %7d seconds wall   (%3d:%02d).\n", $clkTimeS, $clkTimeH, $clkTimeM;
    printf STDOUT "ESTmapper: search required %7d seconds system (%3d:%02d).\n", $sysTimeS, $sysTimeH, $sysTimeM;
    printf STDOUT "ESTmapper: search required %7d seconds user   (%3d:%02d).\n", $usrTimeS, $usrTimeH, $usrTimeM;

    $sysTimeS = $usrTimeS = $clkTimeS = 0;
    open(F, "find $args{'path'}/3-polish -name *.stats -print |");
    while (<F>) {
        chomp;
        open(G, "< $_");
        while (<G>) {
            $clkTimeS += $1 if (m/^clockTime:\s+(\d+\.\d+)$/);
            $sysTimeS += $1 if (m/^systemTime:\s+(\d+\.\d+)$/);
            $usrTimeS += $1 if (m/^userTime:\s+(\d+\.\d+)$/);
        }
        close(G);
    }
    close(F);

    ($clkTimeH,$clkTimeM,$clkTimeS) = sTOhms($clkTimeS);
    ($sysTimeH,$sysTimeM,$sysTimeS) = sTOhms($sysTimeS);
    ($usrTimeH,$usrTimeM,$usrTimeS) = sTOhms($usrTimeS);
    printf STDOUT "ESTmapper: polish required %7d seconds wall   (%3d:%02d).\n", $clkTimeS, $clkTimeH, $clkTimeM;
    printf STDOUT "ESTmapper: polish required %7d seconds system (%3d:%02d).\n", $sysTimeS, $sysTimeH, $sysTimeM;
    printf STDOUT "ESTmapper: polish required %7d seconds user   (%3d:%02d).\n", $usrTimeS, $usrTimeH, $usrTimeM;
}
