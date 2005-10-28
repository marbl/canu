#!/usr/bin/perl

use strict;
use Config;  #  for @signame

use vars qw($bin $gin $scratch $wrk $asm @fragFiles);
use vars qw($vectorIntersect $doOverlapTrimming);
use vars qw($ovlThreads $ovlHashBlockSize $ovlRefBlockSize $ovlMemory $ovlStoreMemory);
use vars qw($frgCorrBatchSize $frgCorrThreads $frgCorrOnGrid);
use vars qw($ovlCorrBatchSize $ovlCorrOnGrid);
use vars qw($numFrags);
use vars qw($useGrid);
use vars qw($processStats);
use vars qw($stoneLevel);
use vars qw($doExtendClearRanges);
use vars qw($doBackups);

use FindBin;
#use lib "$FindBin::Bin/runCAutil";

#  In alphabetical order, instead of execution order.
#require "applyOverlapCorrection.pl";
#require "checkOverlap.pl";
#require "checkPostUnitiggerConsensus.pl";
#require "createConsensusJobs.pl";
#require "createFragmentCorrectionJobs.pl";
#require "createOverlapCorrectionJobs.pl";
#require "createOverlapJobs.pl";
#require "createOverlapStore.pl";
#require "createPostUnitiggerConsensus.pl";
#require "mergeFragmentCorrection.pl";
#require "meryl.pl";
#require "overlapTrim.pl";
#require "preoverlap.pl";
#require "scaffolder.pl";
#require "terminate.pl";
#require "unitigger.pl";
#require "util.pl";


$doBackups = 1;


my $specFile = undef;

while (scalar(@ARGV)) {
    my $arg = shift @ARGV;

    if      ($arg eq "-d") {
        $wrk = shift @ARGV;
    } elsif ($arg eq "-p") {
        $asm = shift @ARGV;
    } elsif ($arg eq "-s") {
        $specFile = shift @ARGV;
    } elsif (-e $arg) {
        push @fragFiles, $arg;
    } else {
        print STDERR "$0: Unknown argument '$arg'\n";
    }
}

setParameters($specFile);

#  We should see assembler binaries in both $bin and $gin, check that
#  it is so.
#
die "Can't find local bin/gatekeeper in $bin\n" if (! -e "$bin/gatekeeper");
die "Can't find grid bin/gatekeeper in $gin\n"  if (! -e "$gin/gatekeeper");


#  For Justin - if we are given no -d, assume the cwd, search for frag
#  files in this directory.
#
if (!defined($wrk)) {
    $wrk = `pwd`;
    chomp $wrk;

    if (scalar(@fragFiles) == 0) {
        open(F, "ls *.frg |");
        @fragFiles = <F>;
        chomp @fragFiles;
        close(F);
    }
}


system("mkdir $wrk") if (! -d $wrk);
die "Directory '$wrk' doesn't exist (-d option).  Fail.\n" if (!defined($wrk) || (! -d $wrk));
die "Run prefix not supplied (-p option).\n" if (!defined($asm));


#  Begin

checkProcessLimits();

preoverlap();
overlapTrim() if ($doOverlapTrimming);
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
extendClearRanges() if ($doExtendClearRanges);
scaffolder("7-CGW");                            #  ignored, if $doExtendClearRanges
createConsensusJobs();                          #  parallel, run manually
terminate();

exit(0);

########################################

