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

if (scalar(@ARGV) != 3) {
    die "wrong args.\n";
}

my $wrkdir   = shift @ARGV;
my $prefix   = shift @ARGV;
my $thisdate = shift @ARGV;

#  Special hack to add group rw

system("chmod -R ug+rw $wrkdir/$thisdate/$prefix");
system("chgrp -R atg   $wrkdir/$thisdate/$prefix");

my $lastdate = undef;
my $gooddate = undef;

my $resultlast = "none";
my $resultgood = "none";

my $diffs = "$wrkdir/$thisdate/$prefix/9-terminator/$prefix.qc";

#  Special case; don't set last/ref to the current date.  Useful in
#  testing (and restarting in general).

if (-e "$wrkdir/POINTERS/$prefix.last") {
    open(F, "< $wrkdir/POINTERS/$prefix.last");
    while (<F>) {
        chomp;
        $lastdate = $_ if ($_ lt $thisdate);
    }
    close(F);
}
if (-e "$wrkdir/POINTERS/$prefix.reference") {
    open(F, "< $wrkdir/POINTERS/$prefix.reference");
    while (<F>) {
        chomp;
        $gooddate = $_ if ($_ lt $thisdate);
    }
    close(F);
}

open(RESULT, "> $wrkdir/$thisdate/$prefix/sanity-result.out") or die;
open(ERROR,  "> $wrkdir/$thisdate/$prefix/sanity-error.out")  or die;

if (! -e "$wrkdir/$thisdate/$prefix/9-terminator/$prefix.asm") {
    print RESULT "Assembly result for $thisdate/$prefix: FAILURE (no asm file)\n";

    open(F, "ls -ltr $wrkdir/$thisdate/$prefix |");
    while (<F>) {
        print ERROR $_;
    }
    close(F);

    print ERROR "\n";
    print ERROR "\n";
    print ERROR "\n";

    my $logfile;

    open(F, "ls $wrkdir/$thisdate/$prefix/runCA.sge.out.[0-9][0-9] |");
    while (<F>) {
        chomp;
        $logfile = $_;
    }
    close(F);

    open(F, "< $logfile");
    while (<F>) {
        print ERROR $_;
    }
    close(F);

    close(RESULT);
    close(ERROR);
    exit(0);
} else {
    if (defined($lastdate) && (-e "$wrkdir/$lastdate/$prefix/9-terminator/$prefix.asm")) {
        $diffs = "$wrkdir/$lastdate/$prefix/9-terminator/$prefix.qc $diffs";
        if (system("diff -q $wrkdir/$thisdate/$prefix/9-terminator/$prefix.asm $wrkdir/$lastdate/$prefix/9-terminator/$prefix.asm > /dev/null 2>&1") == 0) {
            $resultlast = "$lastdate same";
        } else {
            $resultlast = "$lastdate differs";
        }
    }

    if (defined($gooddate) && (-e "$wrkdir/$gooddate/$prefix/9-terminator/$prefix.asm")) {
        $diffs = "$wrkdir/$gooddate/$prefix/9-terminator/$prefix.qc $diffs";
        if (system("diff -q $wrkdir/$thisdate/$prefix/9-terminator/$prefix.asm $wrkdir/$gooddate/$prefix/9-terminator/$prefix.asm > /dev/null 2>&1") == 0) {
            $resultgood = "$gooddate same";
        } else {
            $resultgood = "$gooddate differs";
        }
    }

    print RESULT "Assembly result for $thisdate/$prefix: SUCCESS (last: $resultlast) (reference: $resultgood)\n";

    system("perl $wrkdir/sanity-merge-qc.pl $diffs > $wrkdir/$thisdate/$prefix/sanity-qc.out");

    open(F, ">> $wrkdir/POINTERS/$prefix.last");
    print F "$thisdate\n";
    close(F);
}


close(RESULT);
close(ERROR);
