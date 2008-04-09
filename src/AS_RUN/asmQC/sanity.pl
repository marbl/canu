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

{
    my $oper   = shift @ARGV;
    my $ddir   = shift @ARGV;

    my ($thisdate, $lastdate) = parseDate($ddir);

    if ($oper eq "rsync") {
        print STDERR "NOT RSYNCing UNTIL DEVELOPMENT FINISHED.\n";
        #runCommand($cvsdir, "rsync -av rsync://wgs-assembler.cvs.sourceforge.net/cvsroot/wgs-assembler/\* .");
    }

    if ($oper eq "checkout") {
        checkoutAndLog($thisdate, $lastdate);
    }

    if ($oper eq "build") {
        build($thisdate);
    }

    if ($oper eq "daily") {
        daily($thisdate);
    }

    if ($oper eq "weekly") {
        weekly($thisdate);
    }

    if ($oper eq "monthly") {
        monthly($thisdate);
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

    runCommand($wrkdir, "mkdir $wrkdir/$thisdate");
    runCommand($wrkdir, "mkdir $wrkdir/$thisdate/wgs");


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

    runCommand("$wrkdir/$thisdate/wgs", "cvs -r -R -d $cvsdir -z3 co  -N    -D '$thisdatecvs' src > checkout.err 2>&1");
    runCommand("$wrkdir/$thisdate/wgs", "cvs -r -R -d $cvsdir -z3 log -N -S -d '$lastdatecvs<$thisdatecvs' src > updates-raw");

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

    runCommand("$wrkdir/$thisdate/wgs/src", "gmake > make.out.raw 2> make.err.raw");

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



sub daily ($) {
    my $thisdate = shift @_;

    print STDERR "Submitting daily for $wrkdir/$thisdate\n";
}



sub weekly ($) {
    my $thisdate = shift @_;

    print STDERR "Submitting weekly for $wrkdir/$thisdate\n";
}



sub monthly ($) {
    my $thisdate = shift @_;

    print STDERR "Submitting monthly for $wrkdir/$thisdate\n";
}


exit(0);



#  Utility to run a command and check the exit status, report time used.
#
sub runCommand ($$) {
    my $dir = shift @_;
    my $cmd = shift @_;

    my $t = localtime();
    my $d = time();
    print STDERR "----------------------------------------START $t\n$cmd\n";

    my $rc = 0xffff & system("cd $dir && $cmd");

    $t = localtime();
    print STDERR "----------------------------------------END $t (", time() - $d, " seconds)\n";

    #  Pretty much copied from Programming Perl page 230

    return(0) if ($rc == 0);

    #  Bunch of busy work to get the names of signals.  Is it really worth it?!
    #
    my @signame;
    if (defined($Config{sig_name})) {
        my $i = 0;
        foreach my $n (split('\s+', $Config{sig_name})) {
            $signame[$i] = $n;
            $i++;
        }
    }

    my $error = "ERROR: $cmd\nERROR: Command failed with ";

    if ($rc == 0xff00) {
        $error .= "$!\n";
    } else {
        if ($rc & 0x80) {
            $error .= "coredump from ";
        }
    
        if ($rc > 0x80) {
            $rc >>= 8;
        }
        $rc &= 127;

        if (defined($signame[$rc])) {
            $error .= "signal $signame[$rc] ($rc)\n";
        } else {
            $error .= "signal $rc\n";
        }
    }

    print STDERR $error;

    return(1);
}
