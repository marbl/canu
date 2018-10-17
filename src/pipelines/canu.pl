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
use warnings "all";
no  warnings "uninitialized";

use FindBin;
use Cwd qw(getcwd abs_path);

use lib "$FindBin::RealBin/../lib/site_perl";
#use lib "$FindBin::RealBin/lib/canu/lib/perl5";
#use lib "$FindBin::RealBin/lib/canu/lib64/perl5";

use File::Path 2.08 qw(make_path remove_tree);

use Carp;

use canu::Defaults;
use canu::Execution;

use canu::Configure;

use canu::HaplotypeReads;

use canu::SequenceStore;
use canu::Meryl;
use canu::OverlapInCore;
use canu::OverlapMhap;
use canu::OverlapMMap;
use canu::OverlapStore;

use canu::CorrectReads;

use canu::OverlapBasedTrimming;

use canu::OverlapErrorAdjustment;
use canu::Unitig;
use canu::Consensus;
use canu::Output;

use canu::Grid;
use canu::Grid_Cloud;
use canu::Grid_SGE;
use canu::Grid_Slurm;
use canu::Grid_PBSTorque;
use canu::Grid_LSF;
use canu::Grid_DNANexus;

my @specFiles;       #  Files of specs
my @specOpts;        #  Command line specs
my @inputFiles;      #  Command line inputs, later inputs in spec files are added

my %haplotypeReads;  #  Inpout reads for haplotypes; each element is a NUL-delimited list of files

#  Initialize our defaults.  Must be done before defaults are reported in printOptions() below.

setDefaults();

my $bin     = getBinDirectory();  #  Path to binaries, must be set after setDefaults().
my $cmd     = undef;              #  Temporary string passed to system().
my $asm     = undef;              #  Name of our assembly.
my $asmAuto = undef;              #  If set, the name was auto-discovered.


#  What a mess.  We can't set the version string until after we have a bin directory, and
#  Defaults.pm can't call stuff in Execution.pm.  So, we need to special case setting the version
#  string.

setVersion($bin);

#  Check for the presence of -option switches BEFORE we do any work.
#  This lets us print the default values of options (which we don't do anymore, because
#  too many are set later), and change the defaults (-fast and -slow) before
#  other options are applied.

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

    if ($arg eq "-fast") {
        #  All defaults, unless noted.
        setGlobal("corOverlapper",  "mhap");
        setGlobal("obtOverlapper",  "mhap");    # Changed
        setGlobal("utgOverlapper",  "mhap");    # Changed
        setGlobal("utgReAlign",        "true");    # Changed

        if (-e "$bin/wtdbg-1.2.8") {
           setGlobal("unitigger",      "wtdbg");   # Changed
        }
    }

    if ($arg eq "-accurate") {
        #  All defaults, unless noted.
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

while (scalar(@ARGV)) {
    my $arg = shift @ARGV;

    if     (($arg eq "-h") || ($arg eq "-help") || ($arg eq "--help")) {
        printHelp(1);

    } elsif (($arg eq "-fast") ||
             ($arg eq "-accurate")) {
        addCommandLineOption($arg);

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

    } elsif (($arg eq "-pacbio-raw")         ||  #  File handling is also present in Defaults.pm,
             ($arg eq "-pacbio-corrected")   ||  #  look for addSequenceFile().
             ($arg eq "-nanopore-raw")       ||
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

    } elsif ($arg =~ m/^-haplotype(\w+)$/) {
        my $hapn = $1;
        my $file = $ARGV[0];
        my $fopt = addSequenceFile($readdir, $file, 1);

        while (defined($fopt)) {
            $haplotypeReads{$hapn} .= "$fopt\0";

            addCommandLineOption("$arg '$fopt'");

            shift @ARGV;

            $file = $ARGV[0];
            $fopt = addSequenceFile($readdir, $file);
        }

    } elsif (-e $arg) {
        addCommandLineError("ERROR:  File '$arg' supplied on command line, don't know what to do with it.\n");

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

    open(F, "ls -d . *seqStore |");
    while (<F>) {
        $asm = $1   if (m/^(.*).seqStore$/);
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
checkMinimap($bin);
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
configureCloud($asm);
configureDNANexus();

#  Set jobs sizes based on genomeSize and available hosts;
#  Check that parameters (except error rates) are valid and consistent;
#  Fail if any thing flagged an error condition;

configureAssembler();  #  Set job sizes and etc bases on genomeSize and hosts available.
checkParameters();     #  Check all parameters (except error rates) are valid and consistent.
printHelp();           #  And one final last chance to fail.

#  Make space for us to work in, and move there.

setWorkDirectory($asm, $rootdir);

#  Figure out read inputs.  From an existing store?  From files?  Corrected?  Etc, etc.

my $haveRaw          = 0;
my $haveCorrected    = 0;

my $setUpForPacBio   = 0;
my $setUpForNanopore = 0;

#  If we're a cloud run, fetch the store.

fetchSeqStore($asm);

#  Scan for an existing seqStore and decide what reads and types are in it.  This is
#  pretty horrible and needs to be consolidated into a single report.

my $nCor = getNumberOfReadsInStore($asm, "cor");   #  Number of raw reads ready for correction.
my $nOBT = getNumberOfReadsInStore($asm, "obt");   #  Number of corrected reads ready for OBT.
my $nAsm = getNumberOfReadsInStore($asm, "utg");   #  Number of trimmed reads ready for assembly.

#  If a seqStore was found, scan the reads in it to decide what we're working with.
#  The sqStoreDumpMetaData here isn't used normally; it's only used when canu runs
#  on a seqStore created by hand.

if ($nCor + $nOBT + $nAsm > 0) {
    my $numPacBioRaw         = 0;
    my $numPacBioCorrected   = 0;
    my $numNanoporeRaw       = 0;
    my $numNanoporeCorrected = 0;

    if (! -e "./$asm.seqStore/libraries.txt") {
        if (runCommandSilently(".", "$bin/sqStoreDumpMetaData -S ./$asm.seqStore -libs > ./$asm.seqStore/libraries.txt 2> /dev/null", 1)) {
            caExit("failed to generate list of libraries in store", undef);
        }
    }

    open(L, "< ./$asm.seqStore/libraries.txt") or caExit("can't open './$asm.seqStore/libraries.txt' for reading: $!", undef);
    while (<L>) {
        my @v = split '\s+', $_;

        $numPacBioRaw++           if ($v[2] eq "pacbio-raw");
        $numPacBioCorrected++     if ($v[2] eq "pacbio-corrected");
        $numNanoporeRaw++         if ($v[2] eq "nanopore-raw");
        $numNanoporeCorrected++   if ($v[2] eq "nanopore-corrected");
    }
    close(L);

    $setUpForPacBio++      if ($numPacBioRaw       + $numPacBioCorrected   > 0);
    $setUpForNanopore++    if ($numNanoporeRaw     + $numNanoporeCorrected > 0);

    $haveRaw++             if ($numPacBioRaw       + $numNanoporeRaw       > 0);
    $haveCorrected++       if ($numPacBioCorrected + $numNanoporeCorrected > 0);

    my $rt;

    $rt = "both PacBio and Nanopore"    if (($setUpForPacBio  > 0) && ($setUpForNanopore  > 0));
    $rt = "PacBio"                      if (($setUpForPacBio  > 0) && ($setUpForNanopore == 0));
    $rt = "Nanopore"                    if (($setUpForPacBio == 0) && ($setUpForNanopore  > 0));

    #my $rtct = reportReadsFound($setUpForPacBio, $setUpForNanopore, $haveRaw, $haveCorrected);

    print STDERR "--\n";
    print STDERR "-- In '$asm.seqStore', found $rt reads:\n";
    print STDERR "--   Raw:        $nCor\n";
    print STDERR "--   Corrected:  $nOBT\n";
    print STDERR "--   Trimmed:    $nAsm\n";
}

#  Otherwise, scan input files, counting the different types of libraries we have.

elsif (scalar(@inputFiles) > 0) {
    foreach my $typefile (@inputFiles) {
        my ($type, $file) = split '\0', $typefile;

        $haveCorrected++         if ($type =~ m/corrected/);
        $haveRaw++               if ($type =~ m/raw/);

        $setUpForPacBio++        if ($type =~ m/pacbio/);
        $setUpForNanopore++      if ($type =~ m/nanopore/);
    }

    my $rt;
    my $ct;

    $rt = "both PacBio and Nanopore"    if (($setUpForPacBio  > 0) && ($setUpForNanopore  > 0));
    $rt = "PacBio"                      if (($setUpForPacBio  > 0) && ($setUpForNanopore == 0));
    $rt = "Nanopore"                    if (($setUpForPacBio == 0) && ($setUpForNanopore  > 0));
    $rt = "unknown"                     if (($setUpForPacBio == 0) && ($setUpForNanopore == 0));

    $ct = "uncorrected"                 if (($haveRaw         > 0) && ($haveCorrected    == 0));
    $ct = "corrected"                   if (($haveRaw        == 0) && ($haveCorrected     > 0));
    $ct = "uncorrected AND corrected"   if (($haveRaw         > 0) && ($haveCorrected     > 0));

    #my $rtct = reportReadsFound($setUpForPacBio, $setUpForNanopore, $haveRaw, $haveCorrected);

    print STDERR "--\n";
    print STDERR "-- Found $rt $ct reads in the input files.\n";
}

#  Otherwise, no reads found in a store, and no input files.

else {
    caExit("ERROR: No reads supplied, and can't find any reads in any seqStore", undef);
}

#  Set an initial run mode, based on the libraries we have found, or the stores that exist (unless
#  it was set on the command line).

if (!defined($mode)) {
    $mode = "run"            if ($haveRaw       > 0);   #  If no seqStore, these are set based
    $mode = "trim-assemble"  if ($haveCorrected > 0);   #  on flags describing the input files.

    $mode = "run"            if ($nCor > 0);            #  If a seqStore, these are set based
    $mode = "trim-assemble"  if ($nOBT > 0);            #  on the reads present in the stores.
    $mode = "assemble"       if ($nAsm > 0);
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
    setGlobalIfUndef("obtOvlErrorRate",  0.120);
    setGlobalIfUndef("utgOvlErrorRate",  0.120);
    setGlobalIfUndef("corErrorRate",     0.500);
    setGlobalIfUndef("obtErrorRate",     0.120);
    setGlobalIfUndef("utgErrorRate",     0.120);
    setGlobalIfUndef("cnsErrorRate",     0.200);
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

caExit("ERROR: No reads supplied, and can't find any reads in any seqStore", undef)   if (!defined($mode));
caExit("ERROR: Failed to determine the sequencing technology of the reads", undef)    if (!defined($type));
caExit("ERROR: Can't mix uncorrected and corrected reads", undef)                     if ($haveRaw && $haveCorrected);

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

    make_path("haplotype")    if ((! -d "haplotype")   && ($step eq "haplotype"));
    make_path("correction")   if ((! -d "correction")  && ($step eq "correct"));
    make_path("trimming")     if ((! -d "trimming")    && ($step eq "trim"));
    make_path("unitigging")   if ((! -d "unitigging")  && ($step eq "assemble"));

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
        mhapCheck($asm, $tag, $ovlType)            foreach (1..getGlobal("canuIterationMax") + 1);          #  this also does mhapReAlign

   } elsif (getGlobal("${tag}overlapper") eq "minimap") {
        mmapConfigure($asm, $tag, $ovlType);
        mmapPrecomputeCheck($asm, $tag, $ovlType)  foreach (1..getGlobal("canuIterationMax") + 1);
        mmapCheck($asm, $tag, $ovlType)            foreach (1..getGlobal("canuIterationMax") + 1);

    } else {
        overlapConfigure($asm, $tag, $ovlType);
        overlapCheck($asm, $tag, $ovlType)         foreach (1..getGlobal("canuIterationMax") + 1);
    }

    createOverlapStore($asm, $tag);
}

#
#  Begin pipeline
#

my @haplotypes = sort keys %haplotypeReads;

if (setOptions($mode, "haplotype") eq "haplotype") {
    if ((! -e "./haplotype/haplotyping.success") &&
        (haplotypeReadsExist($asm, @haplotypes) eq "no")) {

        submitScript($asm, undef);   #  See comments there as to why this is safe.

        print STDERR "--\n";
        print STDERR "--\n";
        print STDERR "-- BEGIN HAPLOTYPING\n";
        print STDERR "--\n";

        haplotypeCountConfigure($asm, %haplotypeReads);

        haplotypeCountCheck($asm)                   foreach (1..getGlobal("canuIterationMax") + 1);
        haplotypeMergeCheck($asm, @haplotypes)      foreach (1..getGlobal("canuIterationMax") + 1);
        haplotypeSubtractCheck($asm, @haplotypes)   foreach (1..getGlobal("canuIterationMax") + 1);

        haplotypeReadsConfigure($asm, \@haplotypes, \@inputFiles);
        haplotypeReadsCheck($asm)                   foreach (1..getGlobal("canuIterationMax") + 1);
    }
}

#  If haplotype reads exist, bootstrap the assemblies.
#
#  I tried to use submitScript() to launch these, but that didn't work
#  so nicely - if not on grid, it wouldn't do anything.

if (haplotypeReadsExist($asm, @haplotypes) eq "yes") {
    my $techtype = removeHaplotypeOptions();
    my @options  = getCommandLineOptions();

    #  Find the maximum length of haplotype names, to make the output pretty.

    my $displLen = 0;

    foreach my $haplotype (@haplotypes) {
        my $hapLen = length($haplotype);
        $displLen = ($displLen < $hapLen) ? $hapLen : $displLen;
    }

    #  Decide if we should use or ignore the unassigned reads, and if we should
    #  even bother assembling.

    fetchFile("");

    my %hapReads;
    my %hapBases;

    my $totReads = 0;
    my $totBases = 0;

    open(F, "< haplotype/haplotype.log") or caExit("can't open 'haplotype/haplotype.log' for reading: $!", undef);
    while (<F>) {
        if (m/(\d+)\s+reads\s+(\d+)\s+bases\s+written\s+to\s+haplotype\s+file\s+.*haplotype-(\w+).fasta.gz/) {
            $hapReads{$3} = $1;     $totReads += $1;
            $hapBases{$3} = $2;     $totBases += $2;
        }
    }
    close(F);

    print STDERR "--\n";
    foreach my $haplotype (@haplotypes) {
        printf STDERR "-- Found %8d reads and %12d bases for haplotype $haplotype.\n", $hapReads{$haplotype}, $hapBases{$haplotype};
    }
    printf STDERR "-- Found %8d reads and %12d bases assigned to no haplotype.\n", $hapReads{"unknown"}, $hapBases{"unknown"};

    #  Ignore the unknown reads if there aren't that many.

    my $withUnknown = ($hapBases{"unknown"} / $totBases < 0.02) ? 0 : 1;

    if ($withUnknown == 0) {
        print STDERR "--\n";
        print STDERR "-- Too few bases in unassigned reads to care; don't use them in assemblies.\n";
    }


    #  For each haplotype, emit a script to run the assembly.

    print STDERR "--\n";

    foreach my $haplotype (@haplotypes) {
        my $hs = substr("$haplotype" . " " x $displLen, 0, $displLen);

        print STDERR "-- Assemble haplotype $hs with command './$asm-haplotype$haplotype.sh'.\n";

        open(F, "> ./$asm-haplotype$haplotype.sh");
        print F "#!/bin/sh\n";
        print F "\n";

        if (defined(getGlobal("objectStore"))) {
            print F "\n";
            print F "#  Fetch the haplotyped reads.  This is just a bit weird.\n";
            print F "#  The fetch (boilerplate) only works from within a subdirectory,\n";
            print F "#  so we must cd into it first, fetch, the go back to the root.\n";
            print F "\n";
            print F "mkdir -p haplotype\n";
            print F "cd       haplotype\n";
            print F fetchFileShellCode("haplotype", "haplotype-$haplotype.fasta.gz", "");
            print F fetchFileShellCode("haplotype", "haplotype-unnown.fasta.gz", "")      if ($withUnknown);
            print F "cd ..\n";
        }

        print F "\n";
        print F "$bin/canu \\\n";
        print F "  -p $asm-haplotype$haplotype \\\n";
        print F "  -d $asm-haplotype$haplotype \\\n";
        print F "  $_ \\\n"   foreach (@options);
        print F "  $techtype ./haplotype/haplotype-$haplotype.fasta.gz \\\n";
        print F "  $techtype ./haplotype/haplotype-unknown.fasta.gz \\\n"     if ($withUnknown);
        print F "> ./$asm-haplotype$haplotype.out 2>&1\n";
        print F "\n";
        print F "exit 0\n";
        print F "\n";
        close(F);

        makeExecutable("./$asm-haplotype$haplotype.sh");
    }

    #  Fail if too many unassigned reads.

    print STDERR "--\n";

    if ($hapBases{"unknown"} / $totBases > 0.50) {
        caExit("too many unassigned reads", undef);
    }

    #  Then run the scripts.

    foreach my $haplotype (@haplotypes) {
        if (getGlobal("useGrid") ne "1") {
            print STDERR "-- Starting haplotype assembly $haplotype with script '$rootdir/$asm-haplotype$haplotype.sh'.\n";
        } else {
            print STDERR "-- Submitting '$rootdir/$asm-haplotype$haplotype.sh' to generate haplotype assembly $haplotype.\n";
        }

        runCommand(".", "./$asm-haplotype$haplotype.sh");
    }

    if (getGlobal("useGrid") ne "1") {
        print STDERR "--\n";
        print STDERR "-- Haplotype assemblies completed (or failed).\n";
    } else {
        print STDERR "--\n";
        print STDERR "-- Haplotype assemblies submitted.\n";
    }

    exit(0);
}

if (setOptions($mode, "correct") eq "correct") {
    if ((getNumberOfBasesInStore($asm, "obt") == 0) &&
        (! fileExists("$asm.correctedReads.fasta.gz")) &&
        (! fileExists("$asm.correctedReads.fastq.gz"))) {

        submitScript($asm, undef);   #  See comments there as to why this is safe.

        print STDERR "--\n";
        print STDERR "--\n";
        print STDERR "-- BEGIN CORRECTION\n";
        print STDERR "--\n";

        if (checkSequenceStore($asm, "cor", @inputFiles)) {
            merylConfigure($asm, "cor");
            merylCountCheck($asm, "cor")          foreach (1..getGlobal("canuIterationMax") + 1);
            merylProcessCheck($asm, "cor")        foreach (1..getGlobal("canuIterationMax") + 1);

            overlap($asm, "cor");

            setupCorrectionParameters($asm);

            buildCorrectionLayoutsConfigure($asm);
            buildCorrectionLayoutsCheck($asm)     foreach (1..getGlobal("canuIterationMax") + 1);

            filterCorrectionLayouts($asm);

            generateCorrectedReadsConfigure($asm);
            generateCorrectedReadsCheck($asm)     foreach (1..getGlobal("canuIterationMax") + 1);

            loadCorrectedReads($asm);
        }
    }
}

dumpCorrectedReads($asm);

if ((setOptions($mode, "trim") eq "trim") &&
    (getGlobal("unitigger") ne "wtdbg")) {
    if ((getNumberOfBasesInStore($asm, "utg") == 0) &&
        (! fileExists("$asm.trimmedReads.fasta.gz")) &&
        (! fileExists("$asm.trimmedReads.fastq.gz"))) {

        submitScript($asm, undef);   #  See comments there as to why this is safe.

        print STDERR "--\n";
        print STDERR "--\n";
        print STDERR "-- BEGIN TRIMMING\n";
        print STDERR "--\n";

        if (checkSequenceStore($asm, "obt", @inputFiles)) {
            merylConfigure($asm, "obt");
            merylCountCheck($asm, "obt")     foreach (1..getGlobal("canuIterationMax") + 1);
            merylProcessCheck($asm, "obt")   foreach (1..getGlobal("canuIterationMax") + 1);

            overlap($asm, "obt");

            trimReads($asm);
            splitReads($asm);

            loadTrimmedReads($asm);
        }
    }
}

dumpTrimmedReads ($asm);

if (setOptions($mode, "assemble") eq "assemble") {
    if ((! fileExists("$asm.contigs.fasta")) &&
        (! fileExists("$asm.contigs.fastq"))) {

        submitScript($asm, undef);   #  See comments there as to why this is safe.

        print STDERR "--\n";
        print STDERR "--\n";
        print STDERR "-- BEGIN ASSEMBLY\n";
        print STDERR "--\n";

        if (checkSequenceStore($asm, "utg", @inputFiles)) {
            if (getGlobal("unitigger") ne "wtdbg") {
                merylConfigure($asm, "utg");
                merylCountCheck($asm, "utg")       foreach (1..getGlobal("canuIterationMax") + 1);
                merylProcessCheck($asm, "utg")     foreach (1..getGlobal("canuIterationMax") + 1);

                overlap($asm, "utg");

                #readErrorDetection($asm);

                readErrorDetectionConfigure($asm);
                readErrorDetectionCheck($asm)      foreach (1..getGlobal("canuIterationMax") + 1);

                overlapErrorAdjustmentConfigure($asm);
                overlapErrorAdjustmentCheck($asm)  foreach (1..getGlobal("canuIterationMax") + 1);

                updateOverlapStore($asm);
            }

            unitig($asm);
            unitigCheck($asm)  foreach (1..getGlobal("canuIterationMax") + 1);

            foreach (1..getGlobal("canuIterationMax") + 1) {   #  Consensus wants to change the script between the first and
                consensusConfigure($asm);                      #  second iterations.  The script is rewritten in
                consensusCheck($asm);                          #  consensusConfigure(), so we need to add that to the loop.
            }

            consensusLoad($asm);
            consensusAnalyze($asm);

            if (getGlobal("unitigger") ne "wtdbg") {
                alignGFA($asm)  foreach (1..getGlobal("canuIterationMax") + 1);
            }
            generateOutputs($asm);
        }
    }
}

print STDERR "--\n";
print STDERR "-- Bye.\n";

exit(0);
