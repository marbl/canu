#!/usr/bin/perl

if (scalar(@ARGV) != 6) {
    die "wrong args.\n";
}

my $bindir   = shift @ARGV;
my $wrkdir   = shift @ARGV;
my $prefix   = shift @ARGV;
my $thisdate = shift @ARGV;
my $lastdate = shift @ARGV;
my $gooddate = shift @ARGV;

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

    exit(0);
}

my $resultlast = "none";
my $resultgood = "none";


if (-e "$wrkdir/$lastdate/$prefix/9-terminator/$prefix.asm") {
    if (system("diff -q $wrkdir/$thisdate/$prefix/9-terminator/$prefix.asm $wrkdir/$lastdate/$prefix/9-terminator/$prefix.asm > /dev/null 2>&1") == 0) {
        $resultlast = "same";
    } else {
        $resultlast = "differs";
    }
}

if (-e "$wrkdir/$gooddate/$prefix/9-terminator/$prefix.asm") {
    if (system("diff -q $wrkdir/$thisdate/$prefix/9-terminator/$prefix.asm $wrkdir/$gooddate/$prefix/9-terminator/$prefix.asm > /dev/null 2>&1") == 0) {
        $resultgood = "same";
    } else {
        $resultgood = "differs";
    }
}

print RESULT "Assembly result for $thisdate/$prefix: SUCCESS (last: $resultlast) (reference: $resultgood)\n";

system("perl $wrkdir/mergeqc.pl $wrkdir/$gooddate/$prefix/$prefix.qc $wrkdir/$lastdate/$prefix/$prefix.qc $wrkdir/$thisdate/$prefix/$prefix.qc > $wrkdir/$thisdate/$prefix/sanity-qc.out");

close(RESULT);
close(ERROR);
