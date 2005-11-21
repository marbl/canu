#!/usr/bin/perl

use strict;
use Config;  #  for @signame
use FindBin;

use vars qw($bin $gin $wrk $asm);

use vars qw($numFrags);
use vars qw($useGrid);

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

die "ERROR: I need a directory to run the assembly in (-d option).\n" if (!defined($wrk));
system("mkdir -p $wrk") if (! -d $wrk);
die "ERROR: Directory '$wrk' doesn't exist (-d option).  Fail.\n" if (! -d $wrk);

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
unitigger();
createPostUnitiggerConsensusJobs();             #  parallel, run manually
checkPostUnitiggerConsensus();
scaffolder();
createConsensusJobs();                          #  parallel, run manually
terminate();

exit(0);

