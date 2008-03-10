#!/usr/bin/perl

use strict;
use Config;  #  for @signame

my %logs;
my $file;
my $revs;
my $date;
my $auth;
my $inlog = 0;
my $log;
my %logdate;

my $daten   = `date -u`;
my $datep   = `tail -1 update-history`;

chomp $daten;
chomp $datep;

my $message = "log-$daten";

my $srcdir  = "src";


#  Check out the source if we don't have it.
#
if (! -d $srcdir) {
    #  Anonymous cvs doesn't work for this -- there is apparently a lag behing reality.
    #runCommand(".", "cvs -z3 -d:pserver:anonymous\@wgs-assembler.cvs.sourceforge.net:/cvsroot/wgs-assembler co -P src");
    #runCommand(".", "cvs co -P test");
}


#  Build logs
#
runCommand($srcdir, "cvs -z3 log -N -S -d '$datep<now' > \"../$message\"");

open(F, "< $message");
while (<F>) {
    if (m/^==========.*========$/) {
        $inlog = 0;
        if (!defined($logs{$log})) {
            $logs{$log}  = "$date\t$revs\t$auth\t$file\n";
        } else {
            $logs{$log} .= "$date\t$revs\t$auth\t$file\n";
        }
        $logdate{"$date\0$log"}++;
        undef $log;
    }
    if (m/^----------------------------$/) {
        $inlog = 0;
        if (!defined($logs{$log})) {
            $logs{$log}  = "$date\t$revs\t$auth\t$file\n";
        } else {
            $logs{$log} .= "$date\t$revs\t$auth\t$file\n";
        }
        $logdate{"$date\0$log"}++;
        undef $log;
    }
    if ($inlog) {
        $log .= $_;
    }
    if (m/^RCS\s+file:\s+(.*),v/) {
        $file = $1;
    }
    if (m/^revision\s+(\d+.\d+)/) {
        $revs = $1;
    }
    if (m/^date:\s+(.*);\s+author:\s+(.*);\s+state/) {
        $date = $1;
        $auth = $2;
        $inlog = 1;
    }
}
close(F);


#  Get outta here if there are no updates
#
if (scalar(keys %logdate) == 0) {
    unlink $message;
    exit(0);
}


#  Update our chatter.
#
open(F, ">> update-history");
print F "$daten\n";
close(F);


#  Report what changed.
#
open(F, "> $message");
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


#  Update the source
#
runCommand($srcdir, "cvs update -A -D '$daten' -d");


#  Clean and build it.
#
runCommand(".", "rm -rf FreeBSD* Darwin* Linux*");
runCommand($srcdir, "gmake > /dev/null 2>> \"../$message\"");


#  Run our sanity checks.
#
system("rm -rf f1 f1.out f2 f2.out");

if (runCommand(".", "perl */bin/runCA-OBT.pl -p f1 -d f1 fakeUIDs=1 vectorIntersect=HUREF-D-01-1P5-3KB.vector HUREF-D-01-1P5-3KB.version1.frg > f1.out 2>&1")) {
    system("cat $message f1.out");
}

if (runCommand(".", "perl */bin/runCA-OBT.pl -p f2 -d f2 fakeUIDs=1 HUREF-D-02-1P5-3KB.version2.frg > f2.out 2>&1")) {
    system("cat $message f2.out");
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
