#!/usr/bin/perl
#
# A file with all the (recent) commits to the source tree.  Annotated
# with the time when any tags were placed, search for TAG: or BRANCH:
# to see when.
#
# One can use the following script to generate an update though you'll
# need to edit the $date var to be appropriate for your time frame
# then perl COMMITS >> COMMITS should do the update Might want to
# prune some of the early stuff at some point
#
#--------------------------------------------------------------------------------

use strict;

my %logs;
my $file;
my $revs;
my $auth;
my $inlog = 0;
my $log;
my %logdate;

my $date;
$date = "-d '1 day ago<now'";
$date = "-d '1 week ago<now'";
$date = "-d '5 month ago<now'";
$date = "-d '1 year ago<now'";

open(F, "/usr/bin/cvs -z3 log -N -S $date |");
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

my @keys = sort keys %logdate;

foreach my $l (@keys) {
    my ($d, $l) = split '\0', $l;

    if ((defined($logs{$l})) && (length($logs{$l} > 0))) {
        print "----------------------------------------\n";
        $logs{$l} =~ s!/cvsroot/wgs-assembler/src/!!g;
        print "$logs{$l}\n";
        print "$l\n";
        undef $logs{$l};
    }
}
