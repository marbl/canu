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
use Cwd qw(getcwd abs_path);

use lib "$FindBin::RealBin/lib";
use lib "$FindBin::RealBin/lib/canu/lib/perl5";
use lib "$FindBin::RealBin/lib/canu/lib64/perl5";

use File::Path 2.08 qw(make_path remove_tree);

use Carp;

use canu::Defaults;
use canu::Execution;

use canu::Configure;

use canu::Grid;
use canu::Grid_Cloud;
use canu::Grid_SGE;
use canu::Grid_Slurm;
use canu::Grid_PBSTorque;
use canu::Grid_LSF;
use canu::Grid_DNANexus;

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

my $bin     = getBinDirectory();  #  Path to binaries, reset later.
my $cmd     = undef;              #  Temporary string passed to system().
my $asm     = undef;              #  Name of our assembly.
my $asmAuto = undef;              #  If set, the name was auto-discovered.


#  What a mess.  We can't set the version string until after we have a bin directory, and
#  Defaults.pm can't call stuff in Execution.pm.  So, we need to special case setting the version
#  string.

setVersion($bin);

#  Check for the presence of a -options switch BEFORE we do any work.
#  This lets us print the default values of options.

if (scalar(@ARGV) == 0) {
    printHelp(1);
}

foreach my $arg (@ARGV) {
    if (($arg eq "-options") ||
        ($arg eq "-defaults")) {
        printOptions();
        exit(0);
    }

    if (($arg eq "-version") ||
        ($arg eq "--version")) {
        print getGlobal("version") . "\n";
        exit(0);
    }
}

#  By default, all three steps are run.  Options -correct, -trim and -assemble
#  can limit the pipeline to just that stage.

#  At some pain, we stash the original options for later use.  We need
#  to use these when we resubmit ourself to the grid.  We can't simply dump
#  all of @ARGV into here, because we need to fix up relative paths first.

my $rootdir       = undef;
my $readdir       = undef;
my $mode          = undef;   #  "correct", "trim", "trim-assemble" or "assemble"
my $type          = undef;   #  "pacbio" or "nanopore"
my $step          = "run";
my $haveRaw       = 0;
my $haveCorrected = 0;

while (scalar(@ARGV)) {
    my $arg = shift @ARGV;

    if     (($arg eq "-h") || ($arg eq "-help") || ($arg eq "--help")) {
        printHelp(1);

    } elsif (($arg eq "-citation") || ($arg eq "--citation")) {
        print STDERR "\n";
        printCitation(undef);
        exit(0);

    } elsif ($arg eq "-d") {
        $rootdir = shift @ARGV;

    } elsif ($arg eq "-p") {
        $asm = shift @ARGV;
        addCommandLineOption("-p '$asm'");

    } elsif ($arg eq "-s") {
        my $spec = shift @ARGV;
        $spec = abs_path($spec);

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

    } elsif ($arg eq "-readdir") {
        $readdir = shift @ARGV;
        addCommandLineOption("-readdir '$readdir'");

    } elsif (($arg eq "-pacbio") ||
             ($arg eq "-nanopore")) {
        $type = "pacbio"    if ($arg eq "-pacbio");
        $type = "nanopore"  if ($arg eq "-nanopore");

    } elsif (($arg eq "-pacbio-raw")       ||    #  File handling is also present in
             ($arg eq "-pacbio-corrected") ||    #  Defaults.pm around line 438
             ($arg eq "-nanopore-raw")     ||
             ($arg eq "-nanopore-corrected")) {

        my $file = $ARGV[0];
        my $fopt = addSequenceFile($readdir, $file, 1);

        while (defined($fopt)) {
            push @inputFiles, "$arg\0$fopt";
            addCommandLineOption("$arg '$fopt'");

            shift @ARGV;

            $file = $ARGV[0];
            $fopt = addSequenceFile($readdir, $file);
        }

    } elsif (-e $arg) {
        addCommandLineError("ERROR:  File '$arg' supplied on command line; use -s, -pacbio-raw, -pacbio-corrected, -nanopore-raw, or -nanopore-corrected.\n");

    } elsif ($arg =~ m/=/) {
        push @specOpts, $arg;
        addCommandLineOption("'$arg'");

    } else {
        addCommandLineError("ERROR:  Invalid command line option '$arg'.  Did you forget quotes around options with spaces?\n");
    }
}

#  If no $asm or $dir, see if there is an assembly here.  If so, set $asm to what was found.

if (!defined($asm)) {
    $asmAuto = 1;   #  If we don't actually find a prefix, we'll fail right after this, so OK to set blindly.

    open(F, "ls -d */*gkpStore |");
    while (<F>) {
        $asm = $1   if (m/^correction\/(.*).gkpStore$/);
        $asm = $1   if (m/^trimming\/(.*).gkpStore$/);
        $asm = $1   if (m/^unitigging\/(.*).gkpStore$/);
    }
    close(F);
}

#  Fail if some obvious things aren't set.

addCommandLineError("ERROR:  Assembly name prefix not supplied with -p.\n")   if (!defined($asm));

#  Load paramters from the defaults files

@inputFiles = setParametersFromFile("$bin/canu.defaults", $readdir, @inputFiles)   if (-e "$bin/canu.defaults");
@inputFiles = setParametersFromFile("$ENV{'HOME'}/.canu", $readdir, @inputFiles)   if (-e "$ENV{'HOME'}/.canu");

#  For each of the spec files, parse it, setting parameters and remembering any input files discovered.

foreach my $specFile (@specFiles) {
    @inputFiles = setParametersFromFile($specFile, $readdir, @inputFiles);
}

#  Set parameters from the command line.

setParametersFromCommandLine(@specOpts);

#  Reset $bin, now that all options, specifically the pathMap, are set.

$bin = getBinDirectory();

#  If anything complained (invalid option, missing file, etc) printHelp() will trigger and exit.

printHelp();

#  Now that we know the bin directory, print the version so those pesky users
#  will (hopefully) include it when they paste in logs.

print STDERR "-- " . getGlobal("version") . "\n";
print STDERR "--\n";
print STDERR "-- CITATIONS\n";
print STDERR "--\n";
printCitation("-- ");
print STDERR "-- CONFIGURE CANU\n";
print STDERR "--\n";

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
print STDERR "-- Limited to ", getGlobal("maxMemory"), " gigabytes from maxMemory option.\n"  if (defined(getGlobal("maxMemory")));
print STDERR "-- Limited to ", getGlobal("maxThreads"), " CPUs from maxThreads option.\n"     if (defined(getGlobal("maxThreads")));

detectSGE();
detectSlurm();
detectPBSTorque();
detectLSF();
detectDNANexus();

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
configureDNANexus();

#  Based on genomeSize, configure the execution of every component.
#  This needs to be done AFTER the grid is setup!

configureAssembler();

#  And, finally, move to the assembly directory, finish setting things up, and report the critical
#  parameters.

setWorkDirectory();

if (defined($rootdir)) {
    make_path($rootdir)  if (! -d $rootdir);
    chdir($rootdir);
}

setGlobal("onExitDir", getcwd());
setGlobal("onExitNam", $asm);

setGlobalIfUndef("objectStoreNameSpace", $asm);   #  No good place to put this.

#  Figure out read inputs.  From an existing store?  From files?  Corrected?  Etc, etc.

my $haveCorrected    = 0;
my $haveRaw          = 0;

my $setUpForPacBio   = 0;
my $setUpForNanopore = 0;

#  If we're a cloud run, fetch the store we expect to be working with.

fetchStore("unitigging/$asm.gkpStore")    if ((! -e "unitigging/$asm.gkpStore") && (fileExists("unitigging/$asm.gkpStore.tar")));
fetchStore("trimming/$asm.gkpStore")      if ((! -e "trimming/$asm.gkpStore")   && (fileExists("trimming/$asm.gkpStore.tar"))   && (! -e "unitigging/$asm.gkpStore"));
fetchStore("correction/$asm.gkpStore")    if ((! -e "correction/$asm.gkpStore") && (fileExists("correction/$asm.gkpStore.tar")) && (! -e "trimming/$asm.gkpStore"));

#  Scan for an existing gkpStore.  If the output from that stage exists, ignore the store there.

my $gkp;

$gkp = "correction/$asm.gkpStore"   if ((-e "correction/$asm.gkpStore/libraries.txt")  && (sequenceFileExists("$asm.correctedReads") eq undef));
$gkp = "trimming/$asm.gkpStore"     if ((-e "trimming/$asm.gkpStore/libraries.txt")    && (sequenceFileExists("$asm.trimmedReads")   eq undef));
$gkp = "unitigging/$asm.gkpStore"   if ((-e "unitigging/$asm.gkpStore/libraries.txt"));

#  Scan for existing stage outputs.  These only get used if there isn't a gkpStore found above.

my $reads;

$reads = sequenceFileExists("$asm.correctedReads")  if (!defined($reads));
$reads = sequenceFileExists("$asm.trimmedReads")    if (!defined($reads));

#  A handy function for reporting what reads we found.

sub reportReadsFound ($$$$) {
    my ($setUpForPacBio, $setUpForNanopore, $haveRaw, $haveCorrected) = @_;

    my $rt;
    my $ct;

    $rt = "both PacBio and Nanopore"    if (($setUpForPacBio  > 0) && ($setUpForNanopore  > 0));
    $rt = "PacBio"                      if (($setUpForPacBio  > 0) && ($setUpForNanopore == 0));
    $rt = "Nanopore"                    if (($setUpForPacBio == 0) && ($setUpForNanopore  > 0));
    $rt = "unknown"                     if (($setUpForPacBio == 0) && ($setUpForNanopore == 0));

    $ct = "uncorrected"                 if (($haveRaw         > 0) && ($haveCorrected    == 0));
    $ct = "corrected"                   if (($haveRaw        == 0) && ($haveCorrected     > 0));
    $ct = "uncorrected AND corrected"   if (($haveRaw         > 0) && ($haveCorrected     > 0));

    return("$rt $ct");
}

#  If a gkpStore was found, scan the reads in it to decide what we're working with.

if (defined($gkp)) {
    my $numPacBioRaw         = 0;
    my $numPacBioCorrected   = 0;
    my $numNanoporeRaw       = 0;
    my $numNanoporeCorrected = 0;

    open(L, "< $gkp/libraries.txt") or caExit("can't open '$gkp/libraries.txt' for reading: $!", undef);
    while (<L>) {
        $numPacBioRaw++           if (m/pacbio-raw/);
        $numPacBioCorrected++     if (m/pacbio-corrected/);
        $numNanoporeRaw++         if (m/nanopore-raw/);
        $numNanoporeCorrected++   if (m/nanopore-corrected/);
    }
    close(L);

    $setUpForPacBio++      if ($numPacBioRaw       + $numPacBioCorrected   > 0);
    $setUpForNanopore++    if ($numNanoporeRaw     + $numNanoporeCorrected > 0);

    $haveRaw++             if ($numPacBioRaw       + $numNanoporeRaw       > 0);
    $haveCorrected++       if ($numPacBioCorrected + $numNanoporeCorrected > 0);

    my $rtct = reportReadsFound($setUpForPacBio, $setUpForNanopore, $haveRaw, $haveCorrected);

    print STDERR "--\n";
    print STDERR "-- Found $rtct reads in '$gkp'.\n";
}

#  Like above, scan the gkpStore to decide what we're working with.  The catch here is that
#  we scan the previous store, and all reads are corrected.

elsif (defined($reads)) {

    $gkp = "correction/$asm.gkpStore"   if ((-e "correction/$asm.gkpStore/libraries.txt")  && (sequenceFileExists("$asm.correctedReads")));
    $gkp = "trimming/$asm.gkpStore"     if ((-e "trimming/$asm.gkpStore/libraries.txt")    && (sequenceFileExists("$asm.trimmedReads")));

    my $numPacBio   = 0;
    my $numNanopore = 0;

    if (defined($gkp)) {
        open(L, "< $gkp/libraries.txt") or caExit("can't open '$gkp/libraries.txt' for reading: $!", undef);
        while (<L>) {
            $numPacBio++           if (m/pacbio/);
            $numNanopore++         if (m/nanopore/);
        }
        close(L);

        $setUpForPacBio++      if ($numPacBio    > 0);
        $setUpForNanopore++    if ($numNanopore  > 0);

        $haveCorrected++;
    } else {
        #$setUpForPacBio++;   #  Leaving both setUp's as zero reports 'unknown' and
        $haveCorrected++;     #  defaults to Pacbio below (search for setUpForNanopore).
    }

    #  Regardless of what the user gave us, we always want to restart with these reads.

    undef @inputFiles;
    push  @inputFiles, (($setUpForNanopore == 0) ? "-pacbio" : "-nanopore") . "-corrected\0$reads";

    my $rtct = reportReadsFound($setUpForPacBio, $setUpForNanopore, $haveRaw, $haveCorrected);

    print STDERR "--\n";
    print STDERR "-- Found $rtct reads in '$reads'.\n";
}

#  Scan input files, counting the different types of libraries we have.

elsif (scalar(@inputFiles) > 0) {
    foreach my $typefile (@inputFiles) {
        my ($type, $file) = split '\0', $typefile;

        $haveCorrected++         if ($type =~ m/corrected/);
        $haveRaw++               if ($type =~ m/raw/);

        $setUpForPacBio++        if ($type =~ m/pacbio/);
        $setUpForNanopore++      if ($type =~ m/nanopore/);
    }

    my $rtct = reportReadsFound($setUpForPacBio, $setUpForNanopore, $haveRaw, $haveCorrected);

    print STDERR "--\n";
    print STDERR "-- Found $rtct reads in the input files.\n";
}

#  Set an initial run mode, based on the libraries we have found, or the stores that exist (unless
#  it was set on the command line).

if (!defined($mode)) {
    $mode = "run"            if ($haveRaw       > 0);
    $mode = "trim-assemble"  if ($haveCorrected > 0);

    $mode = "run"            if (-e "correction/$asm.gkpStore/libraries.txt");
    $mode = "trim-assemble"  if (-e "trimming/$asm.gkpStore/libraries.txt");
    $mode = "assemble"       if (-e "unitigging/$asm.gkpStore/libraries.txt");
}

#  Set the type of the reads.  A command line option could force the type, e.g., "-pacbio" or
#  "-nanopore", to let you do cRaZy stuff like "-nanopore -pacbio-raw *fastq".

if (!defined($type)) {
    $type = "pacbio"        if ($setUpForPacBio   > 0);
    $type = "nanopore"      if ($setUpForNanopore > 0);
}

#  Now set error rates (if not set already) based on the dominant read type.

if ($type eq"nanopore") {
    setGlobalIfUndef("corOvlErrorRate",  0.320);
    setGlobalIfUndef("obtOvlErrorRate",  0.144);
    setGlobalIfUndef("utgOvlErrorRate",  0.144);
    setGlobalIfUndef("corErrorRate",     0.500);
    setGlobalIfUndef("obtErrorRate",     0.144);
    setGlobalIfUndef("utgErrorRate",     0.144);
    setGlobalIfUndef("cnsErrorRate",     0.192);
}

if ($type eq"pacbio") {
    setGlobalIfUndef("corOvlErrorRate",  0.240);
    setGlobalIfUndef("obtOvlErrorRate",  0.045);
    setGlobalIfUndef("utgOvlErrorRate",  0.045);
    setGlobalIfUndef("corErrorRate",     0.300);
    setGlobalIfUndef("obtErrorRate",     0.045);
    setGlobalIfUndef("utgErrorRate",     0.045);
    setGlobalIfUndef("cnsErrorRate",     0.075);
}

#  Check for a few errors:
#    no mode                -> don't have any reads or any store to run from.
#    both raw and corrected -> don't know how to process these

caExit("ERROR: No reads supplied, and can't find any reads in any gkpStore", undef)   if (!defined($mode));
caExit("ERROR: Failed to determine the sequencing technology of the reads", undef)    if (!defined($type));
caExit("ERROR: Can't mix uncorrected and corrected reads", undef)                     if ($haveRaw && $haveCorrected);

#  Do a final check on parameters, cleaning up paths and case, and failing on invalid stuff.

checkParameters();

#  And one final last chance to fail - because java and gnuplot both can set an error.

printHelp();

#  Go!

printf STDERR "--\n";
printf STDERR "-- Generating assembly '$asm' in '" . getcwd() . "'\n";
printf STDERR "--\n";
printf STDERR "-- Parameters:\n";
printf STDERR "--\n";
printf STDERR "--  genomeSize        %s\n", getGlobal("genomeSize");
printf STDERR "--\n";
printf STDERR "--  Overlap Generation Limits:\n";
printf STDERR "--    corOvlErrorRate %6.4f (%6.2f%%)\n", getGlobal("corOvlErrorRate"), getGlobal("corOvlErrorRate") * 100.0;
printf STDERR "--    obtOvlErrorRate %6.4f (%6.2f%%)\n", getGlobal("obtOvlErrorRate"), getGlobal("obtOvlErrorRate") * 100.0;
printf STDERR "--    utgOvlErrorRate %6.4f (%6.2f%%)\n", getGlobal("utgOvlErrorRate"), getGlobal("utgOvlErrorRate") * 100.0;
printf STDERR "--\n";
printf STDERR "--  Overlap Processing Limits:\n";
printf STDERR "--    corErrorRate    %6.4f (%6.2f%%)\n", getGlobal("corErrorRate"), getGlobal("corErrorRate") * 100.0;
printf STDERR "--    obtErrorRate    %6.4f (%6.2f%%)\n", getGlobal("obtErrorRate"), getGlobal("obtErrorRate") * 100.0;
printf STDERR "--    utgErrorRate    %6.4f (%6.2f%%)\n", getGlobal("utgErrorRate"), getGlobal("utgErrorRate") * 100.0;
printf STDERR "--    cnsErrorRate    %6.4f (%6.2f%%)\n", getGlobal("cnsErrorRate"), getGlobal("cnsErrorRate") * 100.0;

#  Check that we were supplied a work directory, and that it exists, or we can create it.

make_path("canu-logs")     if (! -d "canu-logs");
make_path("canu-scripts")  if (! -d "canu-scripts");

#  This environment variable tells the binaries to log their execution in canu-logs/

$ENV{'CANU_DIRECTORY'} = getcwd();

#  Report the parameters used.

writeLog();

#  Submit ourself for grid execution?  If not grid enabled, or already running on the grid, this
#  call just returns.  The arg MUST be undef.

submitScript($asm, undef);

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

    make_path("correction")  if ((! -d "correction") && ($step eq "correct"));
    make_path("trimming")    if ((! -d "trimming")   && ($step eq "trim"));
    make_path("unitigging")  if ((! -d "unitigging") && ($step eq "assemble"));

    #  Return that we want to run this step.

    return($step);
}

#
#  Pipeline piece
#

sub overlap ($$) {
    my $asm  = shift @_;
    my $tag  = shift @_;

    my $ovlType = ($tag eq "utg") ? "normal" : "partial";

    if (getGlobal("${tag}overlapper") eq "mhap") {
        mhapConfigure($asm, $tag, $ovlType);

        mhapPrecomputeCheck($asm, $tag, $ovlType)  foreach (1..getGlobal("canuIterationMax") + 1);

        #  this also does mhapReAlign

        mhapCheck($asm, $tag, $ovlType)  foreach (1..getGlobal("canuIterationMax") + 1);

   } elsif (getGlobal("${tag}overlapper") eq "minimap") {
        mmapConfigure($asm, $tag, $ovlType);

        mmapPrecomputeCheck($asm, $tag, $ovlType)  foreach (1..getGlobal("canuIterationMax") + 1);

        mmapCheck($asm, $tag, $ovlType)   foreach (1..getGlobal("canuIterationMax") + 1);

    } else {
        overlapConfigure($asm, $tag, $ovlType);

        overlapCheck($asm, $tag, $ovlType)  foreach (1..getGlobal("canuIterationMax") + 1);
    }

    createOverlapStore($asm, $tag, getGlobal("ovsMethod"));
}

#
#  Begin pipeline
#
#  The checks for sequenceFileExists() at the start aren't needed except for
#  object storage mode.  Gatekeeper has no way of knowing, inside
#  gatekeeper(), that this stage is completed and it shouldn't fetch the
#  store.  In 'normal' operation, the store exists already, and we just
#  return.
#

if (setOptions($mode, "correct") eq "correct") {
    if (sequenceFileExists("$asm.correctedReads") eq undef) {
        print STDERR "--\n";
        print STDERR "--\n";
        print STDERR "-- BEGIN CORRECTION\n";
        print STDERR "--\n";

        gatekeeper($asm, "cor", @inputFiles);

        merylConfigure($asm, "cor");
        merylCheck($asm, "cor")  foreach (1..getGlobal("canuIterationMax") + 1);
        merylProcess($asm, "cor");

        overlap($asm, "cor");

        setupCorrectionParameters($asm);

        buildCorrectionLayoutsConfigure($asm);
        buildCorrectionLayoutsCheck($asm)      foreach (1..getGlobal("canuIterationMax") + 1);

        filterCorrectionLayouts($asm);

        generateCorrectedReadsConfigure($asm);
        generateCorrectedReadsCheck($asm)      foreach (1..getGlobal("canuIterationMax") + 1);

        dumpCorrectedReads($asm);

        buildHTML($asm, "cor");
    }

    my $correctedReads = sequenceFileExists("$asm.correctedReads");

    caExit("can't find corrected reads '$asm.correctedReads*' in directory '" . getcwd() . "'", undef)  if (!defined($correctedReads));

    undef @inputFiles;
    push  @inputFiles, "-$type-corrected\0$correctedReads";
}


if (setOptions($mode, "trim") eq "trim") {
    if (sequenceFileExists("$asm.trimmedReads") eq undef) {
        print STDERR "--\n";
        print STDERR "--\n";
        print STDERR "-- BEGIN TRIMMING\n";
        print STDERR "--\n";

        gatekeeper($asm, "obt", @inputFiles);

        merylConfigure($asm, "obt");
        merylCheck($asm, "obt")  foreach (1..getGlobal("canuIterationMax") + 1);
        merylProcess($asm, "obt");

        overlap($asm, "obt");

        trimReads ($asm);
        splitReads($asm);
        dumpReads ($asm);
        #summarizeReads($asm);

        buildHTML($asm, "obt");
    }

    my $trimmedReads = sequenceFileExists("$asm.trimmedReads");

    caExit("can't find trimmed reads '$asm.trimmedReads*' in directory '" . getcwd() . "'", undef)  if (!defined($trimmedReads));

    undef @inputFiles;
    push  @inputFiles, "-$type-corrected\0$trimmedReads";
}


if (setOptions($mode, "assemble") eq "assemble") {
    if (sequenceFileExists("$asm.contigs") eq undef) {
        print STDERR "--\n";
        print STDERR "--\n";
        print STDERR "-- BEGIN ASSEMBLY\n";
        print STDERR "--\n";

        gatekeeper($asm, "utg", @inputFiles);

        merylConfigure($asm, "utg");
        merylCheck($asm, "utg")  foreach (1..getGlobal("canuIterationMax") + 1);
        merylProcess($asm, "utg");

        overlap($asm, "utg");

        #readErrorDetection($asm);

        readErrorDetectionConfigure($asm);
        readErrorDetectionCheck($asm)  foreach (1..getGlobal("canuIterationMax") + 1);

        overlapErrorAdjustmentConfigure($asm);
        overlapErrorAdjustmentCheck($asm)  foreach (1..getGlobal("canuIterationMax") + 1);

        updateOverlapStore($asm);

        unitig($asm);
        unitigCheck($asm)  foreach (1..getGlobal("canuIterationMax") + 1);

        foreach (1..getGlobal("canuIterationMax") + 1) {   #  Consensus wants to change the script between the first and
            consensusConfigure($asm);                      #  second iterations.  The script is rewritten in
            consensusCheck($asm);                          #  consensusConfigure(), so we need to add that to the loop.
        }

        consensusLoad($asm);
        consensusAnalyze($asm);

        alignGFA($asm)  foreach (1..getGlobal("canuIterationMax") + 1);

        generateOutputs($asm);
    }
}

print STDERR "--\n";
print STDERR "-- Bye.\n";

exit(0);
