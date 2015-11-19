#!/usr/bin/env perl

###############################################################################
 #
 #  This file is part of canu, a software program that assembles whole-genome
 #  sequencing reads into contigs.
 #
 #  This software is based on:
 #    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 #    the 'kmer package' (http://kmer.sourceforge.net)
 #  both originally distributed by Applera Corporation under the GNU General
 #  Public License, version 2.
 #
 #  Canu branched from Celera Assembler at its revision 4587.
 #  Canu branched from the kmer project at its revision 1994.
 #
 #  This file is derived from:
 #
 #    src/pipelines/ca3g.pl
 #
 #  Modifications by:
 #
 #    Brian P. Walenz beginning on 2015-FEB-27
 #      are Copyright 2015 Battelle National Biodefense Institute, and
 #      are subject to the BSD 3-Clause License
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

use strict;

use FindBin;
use Cwd;

use lib "$FindBin::RealBin/lib";
use lib "$FindBin::RealBin/lib/canu/lib/perl5";
use lib "$FindBin::RealBin/lib/canu/lib64/perl5";

use File::Path qw(make_path remove_tree);

use Carp;

use canu::Defaults;
use canu::Execution;

use canu::Gatekeeper;
use canu::Meryl;
use canu::OverlapInCore;
use canu::OverlapMhap;
use canu::OverlapStore;

use canu::CorrectReads;

use canu::OverlapBasedTrimming;

use canu::OverlapErrorAdjustment;
use canu::Unitig;
use canu::Consensus;
use canu::Output;

use canu::HTML;


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

#  Initialize our defaults.  Must be done before defaults are reported in printHelp() below.

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

#  Print usage if no arguments

if (scalar(@ARGV) == 0) {
    setGlobal("help", 1);
    printHelp($bin);
}

#  By default, all three steps are run.  Options -correct, -trim and -assemble
#  can limit the pipeline to just that stage.

#  At some pain, we stash the original options for later use.  We need
#  to use these when we resubmit ourself to the grid.  We can't simply dump
#  all of @ARGV into here, because we need to fix up relative paths first.

my $mode = "run";
my $step = "run";

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


    } elsif ($arg eq "-correct") {
        $mode = $step = "correct";
        addCommandLineOption("-correct");

    } elsif ($arg eq "-trim") {
        $mode = $step = "trim";
        addCommandLineOption("-trim");

    } elsif ($arg eq "-assemble") {
        $mode = $step = "assemble";
        addCommandLineOption("-assemble");

    } elsif ($arg eq "-trim-assemble") {
        $mode = $step = "trim-assemble";
        addCommandLineOption("-trim-assemble");


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

        $mode = "trim-assemble"  if ($arg =~ m/corrected/);

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
my $utg = (-e "$wrk/unitigging/$asm.gkpStore");

if (($cor + $obt + $utg == 0) &&
    (scalar(@inputFiles) == 0)) {
    caExit("no input files specified, and store not already created", undef);
}

#  Check that we were supplied a work directory, and that it exists, or we can create it.

caExit("no run directory (-d option) specified", undef)  if (!defined($wrk));

make_path("$wrk")               if (! -d "$wrk");
make_path("$wrk/canu-logs")     if (! -d "$wrk/canu-logs");
make_path("$wrk/canu-scripts")  if (! -d "$wrk/canu-scripts");

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
    my $mode = shift @_;  #  E.g,. "run" or "trim-assemble" or just plain ol' "trim"
    my $step = shift @_;  #  Step we're setting options for.

    #  Decide if we care about running this step in this mode.  I almost applied
    #  De Morgan's Laws to this.  I don't think it would have been any clearer.

    if (($mode eq $step) ||
        ($mode eq "run") ||
        (($mode eq "trim-assemble") && ($step eq "trim")) ||
        (($mode eq "trim-assemble") && ($step eq "assemble"))) {
        #  Do run this.
    } else {
        return("don't run this");
    }

    #  Create directories for the step, if needed.

    make_path("$wrk/correction")  if ((! -d "$wrk/correction") && ($step eq "correct"));
    make_path("$wrk/trimming")    if ((! -d "$wrk/trimming")   && ($step eq "trim"));
    make_path("$wrk/unitigging")  if ((! -d "$wrk/unitigging") && ($step eq "assemble"));

    #  Now, set options for this step.

    if ($step eq "correct") {
        setErrorRate(0.15);
    }
    if ($step eq "trim") {
        setErrorRate(0.01);
    }
    if ($step eq "assemble") {
        setErrorRate(0.01);
    }

    #  And return that we want to run this step.

    return($step);
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

if (getGlobal("canuIteration") > 0) {
    print STDERR "--\n";
    print STDERR "-- This is canu parallel iteration #" . getGlobal("canuIteration") . ", out of a maximum of " . getGlobal("canuIterationMax") . " attempts.\n";
}

print STDERR "--\n";
print STDERR "-- Final error rates before starting pipeline:\n";

showErrorRates("--   ");

if (setOptions($mode, "correct") eq "correct") {
    print STDERR "--\n";
    print STDERR "--\n";
    print STDERR "-- BEGIN CORRECTION\n";
    print STDERR "--\n";

    gatekeeper($wrk, $asm, "cor", @inputFiles);
    meryl($wrk, $asm, "cor");
    overlap($wrk, $asm, "cor");

    buildCorrectionLayouts($wrk, $asm);

    generateCorrectedReads($wrk, $asm, 1);
    generateCorrectedReads($wrk, $asm, 2);
    generateCorrectedReads($wrk, $asm, 3);

    dumpCorrectedReads($wrk, $asm);

    buildHTML($wrk, $asm, "cor");

    undef @inputFiles;
    push  @inputFiles, "-pacbio-corrected:$wrk/$asm.correctedReads.fastq";
}


if (setOptions($mode, "trim") eq "trim") {
    print STDERR "--\n";
    print STDERR "--\n";
    print STDERR "-- BEGIN TRIMMING\n";
    print STDERR "--\n";

    gatekeeper($wrk, $asm, "obt", @inputFiles);
    meryl($wrk, $asm, "obt");
    overlap($wrk, $asm, "obt");

    trimReads ($wrk, $asm);
    splitReads($wrk, $asm);
    dumpReads ($wrk, $asm);
    #summarizeReads($wrk, $asm);

    buildHTML($wrk, $asm, "obt");

    undef @inputFiles;
    push  @inputFiles, "-pacbio-corrected:$wrk/$asm.trimmedReads.fastq";
}


if (setOptions($mode, "assemble") eq "assemble") {
    print STDERR "--\n";
    print STDERR "--\n";
    print STDERR "-- BEGIN ASSEMBLY\n";
    print STDERR "--\n";

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
    consensusAnalyze($wrk, $asm);
    consensusFilter($wrk, $asm);

    outputGraph($wrk, $asm);
    outputLayout($wrk, $asm);
    outputSequence($wrk, $asm);
    outputSummary($wrk, $asm);
}

exit(0);
