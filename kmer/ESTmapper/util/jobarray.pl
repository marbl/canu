#!/usr/local/bin/perl

#  ESTmapper utility script that runs heavy computes on an LSF farm.
#  This script is designed to be run as a job array.  It takes two
#  arguments:
#    path -- the path to the ESTmapper base directory.
#    type -- "search" or "polish", denoting which stage should be run.

use strict;

my $path = shift @ARGV;
my $type = shift @ARGV;
my $jobi = $ENV{'LSB_JOBINDEX'} - 1;

if (!defined($ENV{'LSB_JOBINDEX'})) {
    print STDERR "ERROR: LSB_JOBINDEX not defined.\n";
    print STDERR "       path = '$path'\n";
    print STDERR "       type = '$type'\n";
    die;
}


if (! -d "$path") {
    print STDERR "ERROR: path does not exist.\n";
    print STDERR "       path = '$path'\n";
    print STDERR "       type = '$type'\n";
    die;
}

#  The big if block below will figure out which jobs need to run, and
#  put the command to run in @jobsToRun.
#
my @jobsToRun;

if ($type eq "search") {

    #  Read the list of segments, and see if the search for that segment is done.
    #
    open(F, "< $path/0-input/scaffolds-list");
    while (<F>) {
        chomp;
        push @jobsToRun, "$path/1-search/$_.cmd" if (! -e "$path/1-search/$_.count");
    }

} elsif ($type eq "polish") {

    #  ESTmapper.pl has already figured out what jobs need to be run.
    #
    open(F, "$path/3-polish/run.sh");
    while (<F>) {
        chomp;
        push @jobsToRun, $_;
    }
    close(F);

} else {
    print STDERR "ERROR: type is not 'search' or 'polish'.\n";
    print STDERR "       path = '$path'\n";
    print STDERR "       type = '$type'\n";
    die;
}


#  If the $jobi'th job isn't finished, run it.  No error checking is
#  really needed (except for diagnostics) beacuse ESTmapper.pl
#  verifies the job completed successfully.
#
if (defined($jobsToRun[$jobi])) {
    system("$jobsToRun[$jobi]");
}
