#!/usr/bin/perl

use strict;
use Config;  #  for @signame
use FindBin;

use vars qw($bin $gin $wrk $asm);

use vars qw($numFrags);

#  Set some not reasonable defaults.
$wrk = undef;
$asm = "asm";

my %global;
my $specFile = undef;
my @specOpts;
my @fragFiles;

while (scalar(@ARGV)) {
    my $arg = shift @ARGV;

    if      ($arg =~ m/^-d/) {
        $wrk = shift @ARGV;
        if ($wrk !~ m!^/!) {
            $wrk = "$ENV{'PWD'}/$wrk";
            print STDERR "WARNING:  Constructing full path to work directory: $wrk\n";
        }
    } elsif ($arg =~ m/^-p/) {
        $asm = shift @ARGV;
    } elsif ($arg =~ m/^-s/) {
        $specFile = shift @ARGV;
    } elsif (-e $arg) {
        push @fragFiles, $arg;
    } elsif ($arg =~ m/=/) {
        push @specOpts, $arg;
    } else {
        print STDERR "$0: Unknown argument '$arg'\n";
    }
}

setDefaults();
setParameters($specFile, @specOpts);
checkDirectories();

#  Begin

preoverlap(@fragFiles);
overlapTrim();
meryl();
createOverlapJobs("normal");                    #  parallel, run manually
checkOverlap("normal");
createOverlapStore();
createFragmentCorrectionJobs();                 #  parallel, run manually
mergeFragmentCorrection();
createOverlapCorrectionJobs();                  #  parallel, run manually
applyOverlapCorrection();
if ($global{'useBogUnitig'}) {
    bogUnitigger();
} else {
    unitigger();
}
createPostUnitiggerConsensusJobs();             #  parallel, run manually
checkPostUnitiggerConsensus();
scaffolder();
createConsensusJobs();                          #  parallel, run manually
terminate();

exit(0);

