#!/usr/bin/perl

use strict;
use Config;  #  for @signame

#  The only two globals:
#
my $wrkdir = "/home/work/nightly";
my $cvsdir = "/home/work/nightly/wgs-assembler-cvs";


#  Command line options are
#
#    'oper'   -- thing to do.
#
#    'ddir' -- properly formatted directory name to do it in.  op ==
#    checkout creates this directory.  If still supplied, it will
#    force a checkout from that date.
#
#  rsync
#     update the local repository in 'wgs-assembler-cvs'.
#
#  checkout date
#     Date must be of the form yyyy-mm-dd-hhmm.  This will checkout
#     the assembler current at that time.  If not present, the current
#     time is assumed.
#
#  build date
#     Builds the assembler checked out into 'date'.
#
#  assemble date specFile ....
#     Launches runCA on each specFile, running it in directory 'date'.
#     Each specFile must have a unique name; the assembly is named
#     after the specFile.
#
#  summarize date
#     Build a summary message for the directory 'date'.
#

{
    my $oper   = shift @ARGV;
    my $ddir   = shift @ARGV;

    my ($thisdate, $lastdate) = parseDate($ddir);

    if ($oper eq "rsync") {
        system("cd $cvsdir && rsync -av rsync://wgs-assembler.cvs.sourceforge.net/cvsroot/wgs-assembler/\* . > rsync.out 2>&1");
    }

    if ($oper eq "checkout") {
        checkoutAndLog($thisdate, $lastdate);
    }

    if ($oper eq "build") {
        build($thisdate);
    }

    if ($oper eq "assemble") {
        assemble($thisdate, @ARGV);
    }

    if (!defined($oper)) {
    }

    exit;
}


sub parseDate ($) {
    my $pathdate = shift @_;
    my $thisdate;
    my $lastdate;

    #  Until our Makefile is fixed, the ":" cannot be used, and so we need
    #  to adopt a slightly goofy timestamp name.

    if (defined($pathdate)) {
        if ($pathdate =~ m/^\d\d\d\d-\d\d-\d\d-\d\d\d\d$/) {
            $thisdate = $pathdate;
        } else {
            die "Malformed ddir '$pathdate'\n";
        }
    } else {
        my @v = localtime;
        $v[5] += 1900;
        $v[4]++;

        $v[5] = substr("0000$v[5]", -4);
        $v[4] = substr("0000$v[4]", -2);
        $v[3] = substr("0000$v[3]", -2);
        $v[2] = substr("0000$v[2]", -2);
        $v[1] = substr("0000$v[1]", -2);

        $thisdate = "$v[5]-$v[4]-$v[3]-$v[2]$v[1]";
    }

    open(F, "ls $wrkdir |");
    while (<F>) {
        chomp;
        if (m/^\d\d\d\d-\d\d-\d\d-\d\d\d\d$/) {
            $lastdate = $_;
        }
    }
    close(F);

    print STDERR "Working on '$thisdate' -- last found is '$lastdate'.\n";

    return($thisdate, $lastdate);
}


sub checkoutAndLog ($$) {
    my $thisdate = shift @_;
    my $lastdate = shift @_;

    print STDERR "Checking out $wrkdir/$thisdate\n";

    if (-d "$wrkdir/$thisdate") {
        print STDERR "$wrkdir/$thisdate already exists.  Please remove to rerun.\n";
        return;
    }

    system("mkdir $wrkdir/$thisdate");
    system("mkdir $wrkdir/$thisdate/wgs");


    #  Convert time-names to dates that cvs can use.

    my $thisdatecvs;
    my $lastdatecvs;

    my $tz = "EDT";

    if ($thisdate =~ m/(\d\d\d\d)-(\d\d)-(\d\d)-(\d\d)(\d\d)/) {
        $thisdatecvs = "$1-$2-$3 $4$5 $tz";
    } else {
    }

    if (defined($lastdate)) {
        if ($lastdate =~ m/(\d\d\d\d)-(\d\d)-(\d\d)-(\d\d)(\d\d)/) {
            $lastdatecvs = "$1-$2-$3 $4$5 $tz";
        } else {
        }
    }

    system("cd $wrkdir/$thisdate/wgs && cvs -r -R -d $cvsdir -z3 co  -N    -D '$thisdatecvs' src > checkout.err 2>&1");
    system("cd $wrkdir/$thisdate/wgs && cvs -r -R -d $cvsdir -z3 log -N -S -d '$lastdatecvs<$thisdatecvs' src > updates-raw");

    my $log;
    my $revs;
    my $date;
    my $auth;
    my $chng;
    my $file;
    my $inlog;
    my %logs;
    my %logdate;

    open(F, "< $wrkdir/$thisdate/wgs/updates-raw");
    while (<F>) {
        if (m/^==========.*========$/) {
            $inlog = 0;
            if (!defined($logs{$log})) {
                $logs{$log}  = "$date\t$revs\t$chng\t$auth\t$file\n";
            } else {
                $logs{$log} .= "$date\t$revs\t$chng\t$auth\t$file\n";
            }
            $logdate{"$date\0$log"}++;
            undef $log;
        }
        if (m/^----------------------------$/) {
            $inlog = 0;
            if (!defined($logs{$log})) {
                $logs{$log}  = "$date\t$revs\t$chng\t$auth\t$file\n";
            } else {
                $logs{$log} .= "$date\t$revs\t$chng\t$auth\t$file\n";
            }
            $logdate{"$date\0$log"}++;
            undef $log;
        }
        if ($inlog) {
            $log .= $_;
        }
        if (m/^RCS\s+file:\s+(.*),v/) {
            $file = $1;
            $file =~ s!$cvsdir/src/!!g;
        }
        if (m/^revision\s+(\d+.\d+)/) {
            $revs = $1;
        }
        if (m/^date:\s+(.*);\s+author:\s+(.*);\s+state.*lines:\s+(.*)\s*$/) {
            $date = $1;
            $auth = $2;
            $chng = $3;
            $inlog = 1;
        }
    }
    close(F);


    open(F, "> $wrkdir/$thisdate/wgs/updates");
    my @keys = sort keys %logdate;
    foreach my $l (@keys) {
        my ($d, $l) = split '\0', $l;
    
        if ((defined($logs{$l})) && (length($logs{$l} > 0))) {
            print F "----------------------------------------\n";
            print F "$logs{$l}\n";
            print F "$l\n";
            undef $logs{$l};
        }
    }
    close(F);

    print STDERR "$thisdate checked out!\n";
}


sub build ($) {
    my $thisdate = shift @_;

    print STDERR "Building $wrkdir/$thisdate\n";

    if (-e "$wrkdir/$thisdate/wgs/src/make.err") {
        print STDERR "$wrkdir/$thisdate was already built once (successuflly or not).  Please cleanup first.\n";
        return;
    }

    system("cd $wrkdir/$thisdate/wgs/src && gmake > make.out.raw 2> make.err.raw");

    my %lines;

    open(F, "< $wrkdir/$thisdate/wgs/src/make.err.raw");
    open(G, "> $wrkdir/$thisdate/wgs/src/make.err");
    while (<F>) {
        chomp;
        next if (m/^ar:\s+creating/);
        if (!exists($lines{$_})) {
            $lines{$_}++;
            print G "$_\n";
        }
    }
    close(F);
    close(G);
}



sub assemble ($@) {
    my $thisdate  = shift @_;
    my $cwd       = $ENV{'PWD'};
    my $holds;

    #  Script that will run (in the assembly directory) when an assembly is finished.
    #
    open(F, "> $wrkdir/$thisdate/assembly-done.sh");
    print F "#!/bin/sh\n";
    print F "ls -l\n";
    print F "ls -ld */*\n";
    close(F);

    #  Script that will run (in the main directory) when all assemblies are finished.
    #
    open(F, "> $wrkdir/$thisdate/summarize.sh");
    print F "#!/bin/sh\n";
    print F "ls -l\n";
    print F "ls -ld */*\n";
    close(F);

    foreach my $s (@_) {
        my ($n, undef) = split '\.', $s;

        print STDERR "----------------------------------------\n";
        print STDERR "Submitting assembly '$n'.\n";

        my $arch;

        #  Ripped from runCA/util.pl
        my $syst = `uname -s`;    chomp $syst;  #  OS implementation
        my $arch = `uname -m`;    chomp $arch;  #  Hardware platform
        my $name = `uname -n`;    chomp $name;  #  Name of the system

        $arch = "amd64"  if ($arch eq "x86_64");
        $arch = "ppc"    if ($arch eq "Power Macintosh");
        #  end of rip.

        system("mkdir $wrkdir/$thisdate/$n");
        system("cd $wrkdir/$thisdate    && perl $wrkdir/$thisdate/wgs/$syst-$arch/bin/runCA -p $n -d $n -s $cwd/$n.specFile sgePropagateHold=ca$n");
        system("cd $wrkdir/$thisdate/$n && qsub -b n -cwd -j y -o assembly-done.out -hold_jid runCA_$n -N ca$n ../assembly-done.sh");

        if (defined($holds)) {
            $holds .= ",ca$n";
        } else {
            $holds  = "ca$n";
        }
    }

    system("cd $wrkdir/$thisdate && qsub -b n -cwd -j y -o summarize.out -hold_jid $holds -N ca$thisdate summarize.sh");
}


exit(0);
