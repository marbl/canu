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
 #    Brian P. Walenz from 2015-FEB-27 to 2015-AUG-26
 #      are Copyright 2015 Battelle National Biodefense Institute, and
 #      are subject to the BSD 3-Clause License
 #
 #    Brian P. Walenz beginning on 2015-NOV-03
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #    Sergey Koren beginning on 2015-NOV-19
 #      are a 'United States Government Work', and
 #      are released in the public domain
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

use canu::Configure;

use canu::Grid;
use canu::Grid_SGE;
use canu::Grid_Slurm;
use canu::Grid_PBSTorque;
use canu::Grid_LSF;

use canu::Gatekeeper;
use canu::Meryl;
use canu::OverlapInCore;
use canu::OverlapMhap;
use canu::OverlapMMap;
use canu::OverlapStore;

use canu::CorrectReads;
use canu::ErrorEstimate;

use canu::OverlapBasedTrimming;

use canu::OverlapErrorAdjustment;
use canu::Unitig;
use canu::Consensus;
use canu::Output;

use canu::HTML;

my @specFiles;    #  Files of specs
my @specOpts;     #  Command line specs
my @inputFiles;   #  Command line inputs, later inputs in spec files are added

#  Initialize our defaults.  Must be done before defaults are reported in printOptions() below.

setDefaults();

#  The bin directory is needed for -version, can only be set after setDefaults(), but really should be
#  set after checkParameters() so it can know pathMap.

my $bin = getBinDirectory();  #  Path to binaries, reset later.
my $cmd = undef;              #  Temporary string passed to system().
my $wrk = undef;              #  Path to our assembly directory.
my $asm = undef;              #  Name of our assembly.

#  Check for the presence of a -options switch BEFORE we do any work.
#  This lets us print the default values of options.

foreach my $arg (@ARGV) {
    if (($arg eq "-options") ||
        ($arg eq "-defaults")) {
        printOptions();
        exit(0);
    }

    if (($arg eq "-version") ||
        ($arg eq "--version")) {
        system("$bin/gatekeeperCreate --version");
        exit(0);
    }
}

#  By default, all three steps are run.  Options -correct, -trim and -assemble
#  can limit the pipeline to just that stage.

#  At some pain, we stash the original options for later use.  We need
#  to use these when we resubmit ourself to the grid.  We can't simply dump
#  all of @ARGV into here, because we need to fix up relative paths first.

my $mode = undef;
my $step = "run";
my $haveRaw = 0;
my $haveCorrected = 0;

while (scalar(@ARGV)) {
    my $arg = shift @ARGV;

    if      ($arg =~ m/^-d/) {
        $wrk = shift @ARGV;
        $wrk = "$ENV{'PWD'}/$wrk" if ($wrk !~ m!^/!);
        addCommandLineOption("-d '$wrk'");
        setGlobal("onExitDir", $wrk);

    } elsif ($arg eq "-p") {
        $asm = shift @ARGV;
        addCommandLineOption("-p '$asm'");
        setGlobal("onExitNam", $asm);

    } elsif ($arg eq "-s") {
        my $spec = shift @ARGV;
        $spec = "$ENV{'PWD'}/$spec" if ($spec !~ m!^/!);

        push @specFiles, $spec;

        addCommandLineOption("-s '$spec'");

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

    } elsif (($arg eq "-pacbio-raw")       ||    #  File handling is also present in
             ($arg eq "-pacbio-corrected") ||    #  Defaults.pm around line 438
             ($arg eq "-nanopore-raw")     ||
             ($arg eq "-nanopore-corrected")) {
        addCommandLineError("ERROR:  File '$ARGV[0]' not found.\n")   if (! -e $ARGV[0]);

        while (-e $ARGV[0]) {
            my $file = shift @ARGV;

            $file = "$ENV{'PWD'}/$file"  if ($file !~ m!^/!);

            push @inputFiles, "$arg\0$file";
            addCommandLineOption("$arg '$file'");
        }

    } elsif (-e $arg) {
        addCommandLineError("ERROR:  File supplied on command line; use -s, -pacbio-raw, -pacbio-corrected, -nanopore-raw, or -nanopore-corrected.\n");

    } elsif ($arg =~ m/=/) {
        push @specOpts, $arg;
        addCommandLineOption("'$arg'");

    } else {
        addCommandLineError("ERROR:  Invalid command line option '$arg'.  Did you forget quotes around options with spaces?\n");
    }
}

#  Fail if some obvious things aren't set.

addCommandLineError("ERROR:  Assembly name prefix not supplied with -p.\n")   if (!defined($asm));
addCommandLineError("ERROR:  Directory not supplied with -d.\n")              if (!defined($wrk));

#  If the mode isn't set - which is allowed only if a gkpStore exists somewhere - be a little smart
#  and figure out which store exists.

$mode = "run"            if (!defined($mode) && (-d "$wrk/correction/$asm.gkpStore"));
$mode = "trim-assemble"  if (!defined($mode) && (-d "$wrk/trimming/$asm.gkpStore"));
$mode = "assemble"       if (!defined($mode) && (-d "$wrk/unitigging/$asm.gkpStore"));

#  Load paramters from the defaults files

@inputFiles = setParametersFromFile("$bin/canu.defaults",   @inputFiles)   if (-e "$bin/canu.defaults");
@inputFiles = setParametersFromFile("$ENV{'HOME'}/.canu",   @inputFiles)   if (-e "$ENV{'HOME'}/.canu");

#  For each of the spec files, parse it, setting parameters and remembering any input files discovered.

foreach my $specFile (@specFiles) {
    @inputFiles = setParametersFromFile($specFile, @inputFiles);
}

#  Set parameters from the command line.

setParametersFromCommandLine(@specOpts);

#  Set parameters based on file types supplied.

foreach my $typefile (@inputFiles) {
    my ($type, $file) = split '\0', $typefile;

    $mode = "trim-assemble"             if (!defined($mode) && ($type =~ m/corrected/));
    $mode = "run"                       if (!defined($mode) && ($type =~ m/raw/));

    $haveCorrected = 1                  if ($type =~ m/corrected/);
    $haveRaw = 1                        if ($type =~ m/raw/);

    setErrorRate(0.015, 0)              if ($type =~ m/pacbio/);
    setGlobal("corErrorRate", "0.30")   if ($type =~ m/pacbio/);

    setErrorRate(0.048, 0)              if ($type =~ m/nanopore/);
    setGlobal("corErrorRate", "0.50")   if ($type =~ m/nanopore/);
}

#  Fail if both raw and corrected are supplied.

addCommandLineError("ERROR:  Canu does not currently support mixing raw and corrected sequences.\n")   if ($haveRaw && $haveCorrected);

#  When resuming a run without input files, set the error rates based on library type in the
#  gkpStore.  If the user set the error rate already, do nothing.
#
#  Also, check if we have gkpStores but no input files and reset error rates based on gkpStore.

if (scalar(@inputFiles) == 0 && ! defined(getGlobal("errorRate"))) {
    my $gkpStore = undef;
    $gkpStore = "$wrk/correction/$asm.gkpStore" if -e "$wrk/correction/$asm.gkpStore/libraries.txt";
    $gkpStore = "$wrk/trimming/$asm.gkpStore"   if -e "$wrk/trimming/$asm.gkpStore/libraries.txt";
    $gkpStore = "$wrk/unitigging/$asm.gkpStore" if -e "$wrk/unitigging/$asm.gkpStore/libraries.txt";

    # set to the default if we can't find anything
    if (!defined($gkpStore)) {
        setErrorRate(0.01, 0);
    } else {
        my $numPacBioRaw         = 0;
        my $numPacBioCorrected   = 0;
        my $numNanoporeRaw       = 0;
        my $numNanoporeCorrected = 0;

        open(L, "< $gkpStore/libraries.txt") or caExit("can't open '$gkpStore/libraries.txt' for reading: $!", undef);
        while (<L>) {
            $numPacBioRaw++           if (m/pacbio-raw/);
            $numPacBioCorrected++     if (m/pacbio-corrected/);
            $numNanoporeRaw++         if (m/nanopore-raw/);
            $numNanoporeCorrected++   if (m/nanopore-corrected/);
        }
        if ($numPacBioRaw > 0 || $numPacBioCorrected > 0) {
            setErrorRate(0.015, 0);
            setGlobal("corErrorRate", "0.30");
            setGlobal("cnsMaxCoverage", 40);
        }
        if ($numNanoporeRaw > 0 || $numNanoporeCorrected > 0) {
            setErrorRate(0.048, 0);
            setGlobal("corErrorRate", "0.50");
            setGlobal("cnsMaxCoverage", 40);
        }
    }
}

#  Finish setting parameters, then reset the bin directory using pathMap.

checkParameters();

$bin = getBinDirectory();

#  If anything complained (invalid option, missing file, etc) printHelp() will trigger and exit.

printHelp();

#  Now that we know the bin directory, print the version so those pesky users
#  will (hopefully) include it when they paste in logs.

printVersion($bin);

#  Check java and gnuplot.

checkJava();
checkGnuplot();

#  And one last chance to fail - because java and gnuplot both can set an error.

printHelp();

#  Detect grid support.  If 'gridEngine' isn't set, the execution methods submitScript() and
#  submitOrRunParallelJob() will return without submitting, or run locally (respectively).  This
#  means that we can leave the default of 'useGrid' to 'true', and execution will do the right thing
#  when there isn't a grid.

print STDERR "-- Detected ", getNumberOfCPUs(), " CPUs and ", getPhysicalMemorySize(), " gigabytes of memory.\n";

detectSGE();
detectSlurm();
detectPBSTorque();
detectLSF();

#  Report if no grid engine found, or if the user has disabled grid support.

if (!defined(getGlobal("gridEngine"))) {
    print STDERR "-- No grid engine detected, grid disabled.\n";
}

if ((getGlobal("useGrid") eq "0") && (defined(getGlobal("gridEngine")))) {
    print STDERR "-- Grid engine disabled per useGrid=false option.\n";
    setGlobal("gridEngine", undef);
}

#  Finish setting up the grid.  This is done AFTER parameters are set from the command line, to
#  let the user override any of our defaults.

configureSGE();
configureSlurm();
configurePBSTorque();
configureLSF();
configureRemote();

#  Based on genomeSize, configure the execution of every component.  This needs to be done AFTER the grid is setup!

configureAssembler();

#  Fail immediately if we run the script on the grid, and the gkpStore directory doesn't exist and
#  we have no input files.  Without this check we'd fail only after being scheduled on the grid.

my $cor = (-e "$wrk/correction/$asm.gkpStore") || sequenceFileExists("$wrk/$asm.correctedReads") || (-e "$wrk/$asm.correctedReads.gkp");
my $obt = (-e "$wrk/trimming/$asm.gkpStore")   || sequenceFileExists("$wrk/$asm.trimmedReads")   || (-e "$wrk/$asm.trimmedReads.gkp");
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

#  This environment variable tells the binaries to log their execution in canu-logs/

$ENV{'CANU_DIRECTORY'} = $wrk;

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

    #  Return that we want to run this step.

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

        mhapPrecomputeCheck($wrk, $asm, $tag, $ovlType)  foreach (1..getGlobal("canuIterationMax") + 1);

        #  this also does mhapReAlign

        mhapCheck($wrk, $asm, $tag, $ovlType)  foreach (1..getGlobal("canuIterationMax") + 1);

   } elsif (getGlobal("${tag}overlapper") eq "minimap") {
        mmapConfigure($wrk, $asm, $tag, $ovlType);

        mmapPrecomputeCheck($wrk, $asm, $tag, $ovlType)  foreach (1..getGlobal("canuIterationMax") + 1);

        mmapCheck($wrk, $asm, $tag, $ovlType)   foreach (1..getGlobal("canuIterationMax") + 1);

    } else {
        overlapConfigure($wrk, $asm, $tag, $ovlType);

        overlapCheck($wrk, $asm, $tag, $ovlType)  foreach (1..getGlobal("canuIterationMax") + 1);
    }

    createOverlapStore($wrk, $asm, $tag, getGlobal("ovsMethod"));
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

    merylConfigure($wrk, $asm, "cor");
    merylCheck($wrk, $asm, "cor")  foreach (1..getGlobal("canuIterationMax") + 1);
    merylProcess($wrk, $asm, "cor");

    overlap($wrk, $asm, "cor");

    buildCorrectionLayouts($wrk, $asm);
    generateCorrectedReads($wrk, $asm)  foreach (1..getGlobal("canuIterationMax") + 1);
    dumpCorrectedReads($wrk, $asm);

    estimateCorrectedError($wrk, $asm, "cor");

    buildHTML($wrk, $asm, "cor");

    my $correctedReads = sequenceFileExists("$wrk/$asm.correctedReads");

    caExit("can't find corrected reads in '$wrk/$asm.correctedReads*'", undef)  if (!defined($correctedReads));

    undef @inputFiles;
    push  @inputFiles, "-pacbio-corrected\0$correctedReads";
}


if (setOptions($mode, "trim") eq "trim") {
    print STDERR "--\n";
    print STDERR "--\n";
    print STDERR "-- BEGIN TRIMMING\n";
    print STDERR "--\n";

    gatekeeper($wrk, $asm, "obt", @inputFiles);

    merylConfigure($wrk, $asm, "obt");
    merylCheck($wrk, $asm, "obt")  foreach (1..getGlobal("canuIterationMax") + 1);
    merylProcess($wrk, $asm, "obt");

    overlap($wrk, $asm, "obt");

    trimReads ($wrk, $asm);
    splitReads($wrk, $asm);
    dumpReads ($wrk, $asm);
    #summarizeReads($wrk, $asm);

    buildHTML($wrk, $asm, "obt");

    my $trimmedReads = sequenceFileExists("$wrk/$asm.trimmedReads");

    caExit("can't find trimmed reads in '$wrk/$asm.trimmedReads*'", undef)  if (!defined($trimmedReads));

    undef @inputFiles;
    push  @inputFiles, "-pacbio-corrected\0$trimmedReads";
}


if (setOptions($mode, "assemble") eq "assemble") {
    print STDERR "--\n";
    print STDERR "--\n";
    print STDERR "-- BEGIN ASSEMBLY\n";
    print STDERR "--\n";

    gatekeeper($wrk, $asm, "utg", @inputFiles);

    merylConfigure($wrk, $asm, "utg");
    merylCheck($wrk, $asm, "utg")  foreach (1..getGlobal("canuIterationMax") + 1);
    merylProcess($wrk, $asm, "utg");

    overlap($wrk, $asm, "utg");

    #readErrorDetection($wrk, $asm);

    readErrorDetectionConfigure($wrk, $asm);
    readErrorDetectionCheck($wrk, $asm)  foreach (1..getGlobal("canuIterationMax") + 1);

    overlapErrorAdjustmentConfigure($wrk, $asm);
    overlapErrorAdjustmentCheck($wrk, $asm)  foreach (1..getGlobal("canuIterationMax") + 1);

    updateOverlapStore($wrk, $asm);

    unitig($wrk, $asm);
    unitigCheck($wrk, $asm)  foreach (1..getGlobal("canuIterationMax") + 1);

    consensusConfigure($wrk, $asm);
    consensusCheck($wrk, $asm)  foreach (1..getGlobal("canuIterationMax") + 1);

    consensusLoad($wrk, $asm);
    consensusAnalyze($wrk, $asm);

    generateOutputs($wrk, $asm);
}

exit(0);
