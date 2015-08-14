#!/usr/bin/perl

use strict;

use FindBin;
use Cwd;

use lib "$FindBin::RealBin";
use lib "$FindBin::RealBin/ca3g/lib/perl5";
use lib "$FindBin::RealBin/ca3g/lib64/perl5";

use File::Path qw(make_path remove_tree);

use Carp;
#use FileHandle;

use POSIX "ceil";
use POSIX "floor";
use POSIX "sys_wait_h";  #  waitpid()

use ca3g::Defaults;
use ca3g::Execution;

use ca3g::Gatekeeper;
use ca3g::Meryl;
use ca3g::OverlapInCore;
use ca3g::OverlapMhap;
use ca3g::OverlapStore;

use ca3g::CorrectReads;

use ca3g::OverlapBasedTrimming;

use ca3g::OverlapErrorAdjustment;
use ca3g::Unitig;
use ca3g::Consensus;
use ca3g::Output;


my $bin = undef;  #  Path to binaries, set once in main.
my $cmd = undef;  #  Temporary string passed to system().
my $wrk = undef;  #  Path to our assembly directory.
my $asm = undef;  #  Name of our assembly.

my %global;       #  Global parameters
my %synops;       #  Global parameters - description
my %synnam;       #  Global parameters - case sensitive name

my @specFiles;    #  Files of specs
my @specOpts;     #  Command line specs
my @inputFiles;   #  Command line inputs, later inputs in spec files are added

print STDERR "-- Detected ", getNumberOfCPUs(), " CPUs and ", getPhysicalMemorySize(), " gigabytes of memory.\n";

#  Initialize our defaults.

setDefaults();

#  Check for the presence of a -options switch BEFORE we do any work.
#  This lets us print the default values of options.

foreach my $arg (@ARGV) {
    if (($arg eq "-options") ||
        ($arg eq "-defaults")) {
        setGlobal("options", 1);
        printHelp($bin);
    }
}

#  The first arg must be the mode of operation.  Currently three modes
#  are supported: correct, trim, assemble.  The modes only affect the
#  steps taken in the 'pipeline' at the bottom of this file.

my $mode = shift @ARGV;
my $step = $mode;

if (($mode ne "run") &&
    ($mode ne "correct") &&
    ($mode ne "trim") &&
    ($mode ne "assemble")) {
    setGlobal("help", "ERROR: first parameter must be the mode of operation (run, correct, trim or assemble).\n");
    printHelp($bin);
    exit(1);
}

addCommandLineOption($mode);

#  At some pain, we stash the original options for later use.  We need
#  to use these when we resubmit ourself to the grid.  We can't simply dump
#  all of @ARGV into here, because we need to fix up relative paths.


#  Enabling/disabling algorithm features is done through library features
#  set in the input gkp files.  This is inconvenient, as you cannot easily
#  change the algorithm without rebuilding gkpStore.  This is flexible, letting
#  you disable an algorithm, or use different parameters for different reads.


while (scalar(@ARGV)) {
    my $arg = shift @ARGV;

    if      ($arg =~ m/^-d/) {
        $wrk = shift @ARGV;
        $wrk = "$ENV{'PWD'}/$wrk" if ($wrk !~ m!^/!);
        addCommandLineOption("-d \"$wrk\"");

    } elsif ($arg eq "-p") {
        $asm = shift @ARGV;
        addCommandLineOption("-p \"$asm\"");

    } elsif ($arg eq "-s") {
        my $spec = shift @ARGV;
        $spec = "$ENV{'PWD'}/$spec" if ($spec !~ m!^/!);

        push @specFiles, $spec;

        addCommandLineOption("-s \"$spec\"");


    } elsif ($arg eq "-version") {
        setGlobal("version", 1);

    } elsif (($arg eq "-options") ||
             ($arg eq "-defaults")) {
        #  Do nothing.  Handled above, but we still need to process it here.
        #setGlobal("options", 1);

    } elsif (($arg eq "-pacbio-raw")       ||
             ($arg eq "-pacbio-corrected") ||
             ($arg eq "-nanopore-raw")     ||
             ($arg eq "-nanopore-corrected")) {
        
        setGlobal("help",
                  getGlobal("help") . "ERROR:  $arg file '$ARGV[0]' not found\n")
            if (! -e $ARGV[0]);

        while (-e $ARGV[0]) {
            my $file = shift @ARGV;

            $file = "$ENV{'PWD'}/$file"  if ($file !~ m!^/!);

            push @inputFiles, "$arg:$file";
            addCommandLineOption("$arg \"$file\"");
        }

    } elsif (-e $arg) {
        setGlobal("help",
                  getGlobal("help") . "ERROR:  File supplied on command line; use -s, -pacbio-raw, -pacbio-corrected, -nanopore-raw, or -nanopore-corrected.\n");

    } elsif ($arg =~ m/=/) {
        push @specOpts, $arg;
        addCommandLineOption("\"$arg\"");

    } else {
        setGlobal("help",
                  getGlobal("help") . "ERROR:  Invalid command line option '$arg'\n");
    }
}


setGlobal("help", getGlobal("help") . "ERROR:  Assembly name prefix not supplied with -p.\n") if (!defined($asm));
setGlobal("help", getGlobal("help") . "ERROR:  Directory not supplied with -d.\n")            if (!defined($wrk));


$bin = getBinDirectory();

#@inputFiles = setParametersFromFile("$bin/spec/runCA.default.specFile", @inputFiles)   if (-e "$bin/spec/runCA.default.specFile");
#@inputFiles = setParametersFromFile("$ENV{'HOME'}/.runCA",              @inputFiles)   if (-e "$ENV{'HOME'}/.runCA");


#  For each of the spec files, parse it, setting parameters and remembering any input files discovered.

foreach my $specFile (@specFiles) {
    @inputFiles = setParametersFromFile($specFile, @inputFiles);
}

#  Set parameters from the command line.

setParametersFromCommandLine(@specOpts);

#  Finish setting parameters.

checkParameters($bin);

#  If anything complained, global{'help'} will be defined, and we'll print help (and the error) and
#  stop.

printHelp($bin);

#  Fail immediately if we run the script on the grid, and the gkpStore directory doesn't exist and
#  we have no input files.  Without this check we'd fail only after being scheduled on the grid.

my $cor = (-e "$wrk/correction/$asm.gkpStore") || (-e "$wrk/$asm.correctedReads.fastq") || (-e "$wrk/$asm.correctedReads.gkp");
my $obt = (-e "$wrk/trimming/$asm.gkpStore")   || (-e "$wrk/$asm.trimmedReads.fastq")   || (-e "$wrk/$asm.trimmedReads.gkp");
my $utg = (-e "$wrk/$asm.gkpStore");

if (($cor + $obt + $utg == 0) &&
    (scalar(@inputFiles) == 0)) {
    caExit("no input files specified, and store not already created", undef);
}

#  Check that we were supplied a work directory, and that it exists, or we can create it.

caExit("no run directory (-d option) specified", undef)  if (!defined($wrk));

make_path("$wrk")             if (! -d "$wrk");
make_path("$wrk/runCA-logs")  if (! -d "$wrk/runCA-logs");

caExit("run directory (-d option) '$wrk' doesn't exist and couldn't be created", undef)  if (! -d $wrk);

#  This environment variable tells the binaries to log their execution in runCA-logs/

$ENV{'AS_RUNCA_DIRECTORY'} = $wrk;

#  Report the parameters used.

writeLog($wrk);

#  Submit ourself for grid execution?  If not grid enabled, or already running on the grid, this
#  call just returns.  The arg MUST be undef.

submitScript($wrk, $asm, undef);





#
#  When doing 'run', this sets options for each stage.
#    - overlapper 'mhap' for correction, 'ovl' for trimming and assembly.
#    - consensus 'falconpipe' for correction, 'utgcns' for assembly.  No consensus in trimming.
#    - errorRates 15% for correction and 2% for trimming and assembly.  Internally, this is
#      multiplied by three for obt, ovl, cns, etc.
#

sub setOptions ($$) {
    my $mode = shift @_;
    my $step = shift @_;

    make_path("$wrk/correction")  if ((! -d "$wrk/correction") && ($step eq "correct"));
    make_path("$wrk/trimming")    if ((! -d "$wrk/trimming")   && ($step eq "trim"));

    return($mode)  if ($mode ne "run");

    if ($step eq "correct") {
        setErrorRate(0.15);

        return($step);
    }

    if ($step eq "trim") {
        setErrorRate(0.01);

        return($step);
    }

    if ($step eq "assemble") {
        setErrorRate(0.01);

        return($step);
    }
}

#
#  Pipeline piece
#

sub overlap ($$$) {
    my $wrk  = shift @_;
    my $asm  = shift @_;
    my $tag  = shift @_;

    my $ovlType = ($tag eq "utg") ? "normal" : "partial";

    if (getGlobal("${tag}overlapper") eq "mhap") {
        mhapConfigure($wrk, $asm, $tag, $ovlType);

        mhapPrecomputeCheck($wrk, $asm, $tag, $ovlType, 1);
        mhapPrecomputeCheck($wrk, $asm, $tag, $ovlType, 2);
        mhapPrecomputeCheck($wrk, $asm, $tag, $ovlType, 3);

        #  this also does mhapReAlign

        mhapCheck($wrk, $asm, $tag, $ovlType, 1);
        mhapCheck($wrk, $asm, $tag, $ovlType, 2);
        mhapCheck($wrk, $asm, $tag, $ovlType, 3);

    } else {
        overlapConfigure($wrk, $asm, $tag, $ovlType);

        overlapCheck($wrk, $asm, $tag, $ovlType, 1);
        overlapCheck($wrk, $asm, $tag, $ovlType, 2);
        overlapCheck($wrk, $asm, $tag, $ovlType, 3);
    }

    createOverlapStore($wrk, $asm, $tag, getGlobal("ovlStoreMethod"));
}

#
#  Begin pipeline
#

print STDERR "--\n";
print STDERR "-- This is ca3g parallel iteration #" . getGlobal("ca3gIteration") . ", out of a maximum of " . getGlobal("ca3gIterationMax") . " attempts.\n";
print STDERR "--\n";
print STDERR "-- Final error rates before starting pipeline:\n";

showErrorRates("-- ");

if (setOptions($mode, "correct") eq "correct") {
    print STDERR "----------------------------------------BEGIN CORRECTION\n";

    gatekeeper($wrk, $asm, "cor", @inputFiles);
    meryl($wrk, $asm, "cor");
    overlap($wrk, $asm, "cor");

    buildCorrectionLayouts($wrk, $asm);

    generateCorrectedReads($wrk, $asm, 1);
    generateCorrectedReads($wrk, $asm, 2);
    generateCorrectedReads($wrk, $asm, 3);

    dumpCorrectedReads($wrk, $asm);

    undef @inputFiles;
    push  @inputFiles, "-pacbio-corrected:$wrk/correction/$asm.correctedReads.fastq";
}


if (setOptions($mode, "trim") eq "trim") {
    print STDERR "----------------------------------------BEGIN TRIMMING\n";

    gatekeeper($wrk, $asm, "obt", @inputFiles);
    meryl($wrk, $asm, "obt");
    overlap($wrk, $asm, "obt");

    trimReads  ($wrk, $asm);
    splitReads ($wrk, $asm);
    dumpReads  ($wrk, $asm);
    #summarizeReads($wrk, $asm);

    undef @inputFiles;
    push  @inputFiles, "-pacbio-corrected:$wrk/trimming/$asm.trimmedReads.fastq";
}


if (setOptions($mode, "assemble") eq "assemble") {
    print STDERR "----------------------------------------BEGIN ASSEMBLY\n";

    gatekeeper($wrk, $asm, "utg", @inputFiles);
    meryl($wrk, $asm, "utg");
    overlap($wrk, $asm, "utg");

    #readErrorDetection($wrk, $asm);

    readErrorDetectionConfigure($wrk, $asm);

    readErrorDetectionCheck($wrk, $asm, 1);
    readErrorDetectionCheck($wrk, $asm, 2);
    readErrorDetectionCheck($wrk, $asm, 3);

    overlapErrorAdjustmentConfigure($wrk, $asm);

    overlapErrorAdjustmentCheck($wrk, $asm, 1);
    overlapErrorAdjustmentCheck($wrk, $asm, 2);
    overlapErrorAdjustmentCheck($wrk, $asm, 3);

    updateOverlapStore($wrk, $asm);

    unitig($wrk, $asm);

    consensusConfigure($wrk, $asm);
    consensusCheck($wrk, $asm, 1);
    consensusCheck($wrk, $asm, 2);
    consensusCheck($wrk, $asm, 3);

    consensusLoad($wrk, $asm);
    consensusFilter($wrk, $asm);

    outputGraph($wrk, $asm);
    outputLayout($wrk, $asm);
    outputSequence($wrk, $asm);
}

exit(0);
