#!/usr/bin/env perl

###############################################################################
 #
 #  This file is part of canu, a software program that assembles whole-genome
 #  sequencing reads into contigs.
 #
 #  This software is based on:
 #    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 #    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 #
 #  Except as indicated otherwise, this is a 'United States Government Work',
 #  and is released in the public domain.
 #
 #  File 'README.licenses' in the root directory of this distribution
 #  contains full conditions and disclaimers.
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

my $citationMode = "canu";        #  For printing the citation(s).

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

    if (($arg eq "-citation")  || ($arg eq "--citation") ||
        ($arg eq "-citations") || ($arg eq "--citations")) {
        print STDERR "\n";
        printCitation(undef, "all");
        exit(0);
    }

    if ($arg =~ m/-haplotype\w+/) {
        $citationMode = "trio";
    }

    if ($arg eq "-pacbio-hifi") {
        $citationMode = "hicanu";
    }

    if ($arg eq "-fast") {
        #  All defaults, unless noted.
        setGlobal("corOverlapper",  "mhap");
        setGlobal("obtOverlapper",  "mhap");
        setGlobal("utgOverlapper",  "mhap");
        setGlobal("utgReAlign",     "true");
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

my $rootdir          = undef;
my $readdir          = undef;

my $mode             = undef;   #  "haplotype", "correct", "trim", "trim-assemble" or "assemble"
my $step             = "run";   #  Step to start at (?)

#  If no seqStore exists, we depend on the user setting command line
#  options to tell us the status of the reads.
my $readsAreRaw      = 1;       #  They're either raw or corrected.
my $readsAreTrimmed  = 0;       #  They're either untrimmed or trimmed.

while (scalar(@ARGV)) {
    my $arg = shift @ARGV;

    if     (($arg eq "-h") || ($arg eq "-help") || ($arg eq "--help")) {
        printHelp(1);
    }

    elsif (($arg eq "-fast") ||
           ($arg eq "-accurate")) {
        addCommandLineOption($arg);
    }

    elsif ($arg eq "-d") {
        $rootdir = shift @ARGV;
    }

    elsif ($arg eq "-p") {
        $asm = shift @ARGV;
        addCommandLineOption("-p '$asm'");
    }

    elsif ($arg eq "-s") {
        my $spec = shift @ARGV;
        $spec = abs_path($spec);

        push @specFiles, $spec;

        addCommandLineOption("-s '$spec'");
    }

    elsif ($arg eq "-readdir") {
        $readdir = shift @ARGV;
        addCommandLineOption("-readdir '$readdir'");
    }

    #  Set the type of read we're getting.

    elsif ($arg eq "-corrected") {
        $readsAreRaw     = 0;
        addCommandLineOption($arg);
    }

    elsif ($arg eq "-trimmed") {
        $readsAreTrimmed = 1;
        addCommandLineOption($arg);
    }

    #  Remember the read files we're given.

    elsif (($arg eq "-pacbio") ||
           ($arg eq "-nanopore") ||
           ($arg eq "-pacbio-raw") ||           # Compatibility with v1.x.
           ($arg eq "-nanopore-raw") ||         # Compatibility with v1.x.
           ($arg eq "-pacbio-corrected") ||     # Compatibility with v1.x.
           ($arg eq "-nanopore-corrected") ||   # Compatibility with v1.x.
           ($arg eq "-pacbio-hifi")) {

      anotherFile:
        my $file = $ARGV[0];
        my $fopt = addSequenceFile($readdir, $file);

        if (defined($fopt)) {
            push @inputFiles, "$arg\0$fopt";
            addCommandLineOption("$arg '$fopt'");
            shift @ARGV;
            goto anotherFile;
        }
    }

    #  Trio reads

    elsif ($arg =~ m/^-haplotype(\w+)$/) {
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
    }

    #  Running only a specific stage.

    elsif ($arg eq "-haplotype") {
        $mode = $step = "haplotype";
        addCommandLineOption("-haplotype");
    }

    elsif ($arg eq "-correct") {
        $mode = $step = "correct";
        addCommandLineOption("-correct");
    }

    elsif ($arg eq "-trim") {
        $mode = $step = "trim";
        addCommandLineOption("-trim");
    }

    elsif ($arg eq "-trim-assemble") {
        $mode = $step = "trim-assemble";
        addCommandLineOption("-trim-assemble");
    }

    elsif ($arg eq "-assemble") {
        $readsAreTrimmed = 1;
        $mode = $step = "assemble";
        addCommandLineOption("-assemble");
    }

    #  General options and errors.

    elsif (-e $arg) {
        addCommandLineError("ERROR:  File '$arg' supplied on command line, don't know what to do with it.\n");
    }

    elsif ($arg =~ m/=/) {
        push @specOpts, $arg;
        addCommandLineOption("'$arg'");
    }

    else {
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

addCommandLineError("ERROR:  Assembly name prefix (-p) not supplied.\n")   if (!defined($asm));
addCommandLineError("ERROR:  Assembly name prefix (-p) cannot contain the path delimiter '/'.\n")   if ($asm =~ m!/!);
addCommandLineError("ERROR:  Assembly name prefix (-p) cannot contain spaces.\n")   if ($asm =~ m!\s!);

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
printCitation("-- ", $citationMode);
print STDERR "-- CONFIGURE CANU\n";
print STDERR "--\n";

#  Check java and gnuplot.

checkJava();
checkMinimap($bin);
checkGnuplot();

#  And one last chance to fail - because java and gnuplot both can set an error.

printHelp();


################################################################################



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
#    If we're off grid (disabled by user, ignore stage directory
#    for remote, we leave it alone since we still want the user-submitted jobs to use staging
if (!defined(getGlobal("gridEngine"))) {
    print STDERR "-- No grid engine detected, grid and staging disabled.\n";
    setGlobal("stageDirectory", undef);
}

if ((getGlobal("useGrid") eq "0") && (defined(getGlobal("gridEngine")))) {
    print STDERR "-- Grid engine and staging disabled per useGrid=false option.\n";
    setGlobal("gridEngine", undef);
    setGlobal("stageDirectory", undef);
}

#  Finish setting up the grid.  This is done AFTER parameters are set from the command line, to
#  let the user override any of our defaults.

configureSGE();
configureSlurm();
configurePBSTorque();
configureLSF();
configureRemote();
configureCloud($asm, $rootdir);
configureDNANexus();

#  Set jobs sizes based on genomeSize and available hosts;
#  Check that parameters (except error rates) are valid and consistent;
#  Fail if any thing flagged an error condition;

configureAssembler();  #  Set job sizes and etc bases on genomeSize and hosts available.

#  Make space for us to work in, and move there.

setWorkDirectory($asm, $rootdir);


################################################################################


#  If we're a cloud run, fetch the store.

fetchSeqStore($asm);

#  Scan for an existing seqStore and decide what reads and types are in it.  This is
#  pretty horrible and needs to be consolidated into a single report.

#my $nCor = getNumberOfReadsInStore($asm, "cor");   #  Number of raw reads ready for correction.
#my $nOBT = getNumberOfReadsInStore($asm, "obt");   #  Number of corrected reads ready for OBT.
#my $nAsm = getNumberOfReadsInStore($asm, "utg");   #  Number of trimmed reads ready for assembly.

my $numRaw    = 0;
my $numRawTri = 0;
my $numCor    = 0;
my $numCorTri = 0;

my $numPacBio     = 0;
my $numNanopore   = 0;
my $numHiFi       = 0;

#  If a seqStore was found, scan the reads in it to decide what we're working with.

#if ($nCor + $nOBT + $nAsm > 0) {
if (-e "./$asm.seqStore/info") {
    ($numRaw,
     $numRawTri,
     $numCor,
     $numCorTri,
     $numPacBio,
     $numNanopore,
     $numHiFi) = getSequenceStoreStats($asm);

    my $rt;

    if      (($numPacBio == 0) && ($numNanopore == 0) && ($numHiFi == 0)) {
        $rt = "NO READS AT ALL";
    } elsif (($numPacBio == 0) && ($numNanopore == 0) && ($numHiFi != 0)) {
        $rt = "PacBio HiFi";
    } elsif (($numPacBio == 0) && ($numNanopore != 0) && ($numHiFi == 0)) {
        $rt = "Nanopore";
    } elsif (($numPacBio == 0) && ($numNanopore != 0) && ($numHiFi != 0)) {
        $rt = "Nanopore and PacBio HiFi";
    } elsif (($numPacBio != 0) && ($numNanopore == 0) && ($numHiFi == 0)) {
        $rt = "PacBio CLR";
    } elsif (($numPacBio != 0) && ($numNanopore == 0) && ($numHiFi != 0)) {
        $rt = "PacBio CLR and PacBio HiFi";
    } elsif (($numPacBio != 0) && ($numNanopore != 0) && ($numHiFi == 0)) {
        $rt = "PacBio CLR and Nanopore";
    } elsif (($numPacBio != 0) && ($numNanopore != 0) && ($numHiFi != 0)) {
        $rt = "PacBio CLR, Nanopore and PacBio HiFi";
    }

    print STDERR "--\n";
    print STDERR "-- In '$asm.seqStore', found $rt reads:\n";
    print STDERR "--   PacBio CLR:               $numPacBio\n"     if ($numPacBio   > 0);
    print STDERR "--   Nanopore:                 $numNanopore\n"   if ($numNanopore > 0);
    print STDERR "--   PacBio HiFi:              $numHiFi\n"       if ($numHiFi     > 0);
    print STDERR "--\n";
    print STDERR "--   Raw:                      $numRaw\n"        if ($numRaw      > 0);
    print STDERR "--   Raw and Trimmed:          $numRawTri\n"     if ($numRawTri   > 0);
    print STDERR "--   Corrected:                $numCor\n"        if ($numCor      > 0);
    print STDERR "--   Corrected and Trimmed:    $numCorTri\n"     if ($numCorTri   > 0);
}

#  Otherwise, scan input files, counting the different types of libraries we have.
#  At the same time, rewrite the sqStoreCreate option to accound for read status.

elsif (scalar(@inputFiles) > 0) {

    foreach my $tf (@inputFiles) {
        my ($t, $f) = split '\0', $tf;

        #  See setParametersFromFile() in Defaults.pm for more
        #  COMPATIBILITY MODE stuff.

        if ($t =~ m/^(-.*)-raw/) {
            print STDERR "-- WARNING:\n";
            print STDERR "-- WARNING:  Option '$t <files>' is deprecated.\n";
            print STDERR "-- WARNING:  Use option '$1 <files>' in the future.\n";
            print STDERR "-- WARNING:\n";
            $t = $1;
            $readsAreRaw = 1;
        }

        if ($t =~ m/^(-.*)-corrected/) {
            print STDERR "-- WARNING:\n";
            print STDERR "-- WARNING:  Option '$t' is deprecated.\n";
            print STDERR "-- WARNING:  Use options '-corrected $1 <files>' in the future.\n";
            print STDERR "-- WARNING:\n";
            $t = $1;
            $readsAreRaw = 0;
        }

        if      ($t eq "-pacbio") {
            $numPacBio++;
        } elsif ($t eq "-nanopore") {
            $numNanopore++;
        } elsif ($t eq "-pacbio-hifi") {
            $numHiFi++;
        } else {
            caExit("ERROR: Unknown type '$t' for file '$f'\n", undef);
        }

        $tf = "$t\0$f";
    }

    my $ct;
    my $rt;
    my $st;

    #  If no mode set, default -pacbio-hifi to trimmed status.

    if (($mode eq "") && ($numHiFi > 0)) {
        $readsAreTrimmed = 1;
    }

    #  Figure out a human description of the reads, and set
    #  flags to pass to sqStoreCreate.

    if      (($readsAreRaw == 0) && ($readsAreTrimmed == 0)) {
        $numCor = 1;
        $ct = "untrimmed corrected";
        $st = "-corrected";
    } elsif (($readsAreRaw == 0) && ($readsAreTrimmed != 0)) {
        $numCorTri = 1;
        $ct = "trimmed corrected";
        $st = "-corrected -trimmed";
    } elsif (($readsAreRaw != 0) && ($readsAreTrimmed == 0)) {
        $numRaw = 1;
        $ct = "untrimmed raw";
        $st = "-raw";
    } elsif (($readsAreRaw != 0) && ($readsAreTrimmed != 0)) {
        $numRawTri = 1;
        $ct = "trimmed raw";
        $st = "-raw -trimmed";
    }

    if      (($numPacBio == 0) && ($numNanopore == 0) && ($numHiFi == 0)) {
        $rt = "NO READS AT ALL";
    } elsif (($numPacBio == 0) && ($numNanopore == 0) && ($numHiFi != 0)) {
        $rt = "PacBio HiFi";
    } elsif (($numPacBio == 0) && ($numNanopore != 0) && ($numHiFi == 0)) {
        $rt = "Nanopore";
    } elsif (($numPacBio == 0) && ($numNanopore != 0) && ($numHiFi != 0)) {
        $rt = "Nanopore and PacBio HiFi";
    } elsif (($numPacBio != 0) && ($numNanopore == 0) && ($numHiFi == 0)) {
        $rt = "PacBio CLR";
    } elsif (($numPacBio != 0) && ($numNanopore == 0) && ($numHiFi != 0)) {
        $rt = "PacBio CLR and PacBio HiFi";
    } elsif (($numPacBio != 0) && ($numNanopore != 0) && ($numHiFi == 0)) {
        $rt = "PacBio CLR and Nanopore";
    } elsif (($numPacBio != 0) && ($numNanopore != 0) && ($numHiFi != 0)) {
        $rt = "PacBio CLR, Nanopore and PacBio HiFi";
    }

    foreach my $tf (@inputFiles) {
        $tf = "$st $tf";
    }

    print STDERR "--\n";
    print STDERR "-- Found $ct $rt reads in the input files.\n";

}

#  Otherwise, no reads found in a store, and no input files.

else {
    caExit("ERROR: No reads supplied, and can't find a seqStore", undef);
}

#  Set the type of the reads, and basic parameters tuned for that type.

if (($numPacBio   == 0) &&
    ($numNanopore == 0) &&
    ($numHiFi     == 0)) {
    caExit("ERROR: Failed to determine the sequencing technology of the reads", undef);
}

if (($numPacBio   >= 0) &&
    ($numNanopore  > 0) &&
    ($numHiFi     == 0)) {
    setGlobalIfUndef("corOvlErrorRate",  0.320);
    setGlobalIfUndef("obtOvlErrorRate",  0.120);
    setGlobalIfUndef("utgOvlErrorRate",  0.120);
    setGlobalIfUndef("corErrorRate",     0.500);
    setGlobalIfUndef("obtErrorRate",     0.120);
    setGlobalIfUndef("utgErrorRate",     0.120);
    setGlobalIfUndef("cnsErrorRate",     0.200);
}

if (($numPacBio    > 0) &&
    ($numNanopore >= 0) &&
    ($numHiFi     == 0)) {
    setGlobalIfUndef("corOvlErrorRate",  0.240);
    setGlobalIfUndef("obtOvlErrorRate",  0.045);
    setGlobalIfUndef("utgOvlErrorRate",  0.045);
    setGlobalIfUndef("corErrorRate",     0.300);
    setGlobalIfUndef("obtErrorRate",     0.045);
    setGlobalIfUndef("utgErrorRate",     0.045);
    setGlobalIfUndef("cnsErrorRate",     0.075);
}

if (($numPacBio   == 0) &&
    ($numNanopore == 0) &&
    ($numHiFi      > 0)) {
    setGlobalIfUndef("corOvlErrorRate",  0.000);
    setGlobalIfUndef("obtOvlErrorRate",  0.025);
    setGlobalIfUndef("utgOvlErrorRate",  0.010);
    setGlobalIfUndef("corErrorRate",     0.000);
    setGlobalIfUndef("obtErrorRate",     0.025);
    setGlobalIfUndef("utgErrorRate",     0.010);
    setGlobalIfUndef("cnsErrorRate",     0.050);
    setGlobalIfUndef("homoPolyCompress", 1);
    setGlobalIfUndef("maxInputCoverage", 50);
    setGlobalIfUndef("batOptions",       "-eg 0.0003 -sb 0.01 -dg 0 -db 3 -dr 0 -ca 1000 -cp 10");
}

if (($numPacBio > 0 || $numNanopore >0) && $numHiFi > 0) {
   caExit("ERROR: HiFi data cannot currently be combined with another read type", undef);
}

#  Set an initial run mode based on what we discovered above.

if (!defined($mode)) {
    $mode = "run"            if ($numRaw    > 0);
    $mode = "trim-assemble"  if ($numCor    > 0);
    $mode = "assemble"       if ($numHiFi   > 0);
    $mode = "assemble"       if ($numCorTri > 0);
}

checkParameters();     #  Check all parameters (except error rates) are valid and consistent.
printHelp();           #  And one final last chance to fail.

################################################################################

#  Go!

printf STDERR "--\n";
printf STDERR "-- Generating assembly '$asm' in '" . getcwd() . "':\n";

if ($mode eq "run") {
    print STDERR "--    - separate reads into haplotypes.\n"   if (scalar(keys %haplotypeReads) > 0);
    print STDERR "--    - correct raw reads.\n";
    print STDERR "--    - trim corrected reads.\n";
    print STDERR "--    - assemble corrected and trimmed reads.\n";
}

if ($mode eq "haplotype") {
    print STDERR "--    - only separate reads into haplotypes.\n";
}

if ($mode eq "correct") {
    print STDERR "--    - only correct raw reads.\n";
}

if ($mode eq "trim") {
    print STDERR "--    - only trim corrected reads.\n";
}

if ($mode eq "trim-assemble") {
   print STDERR "--    - trim corrected reads.\n";
   print STDERR "--    - assemble corrected and trimmed reads.\n";
}

if ($mode eq "assemble") {
    if ($numHiFi  > 0) {
       print STDERR "--    - assemble HiFi reads.\n";
    } else {
       print STDERR "--    - assemble corrected and trimmed reads.\n"    if ($numHiFi == 0);
   }
}

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

if ((scalar(@haplotypes) > 0) &&
    (setOptions($mode, "haplotype") eq "haplotype")) {
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
            $hapReads{$3} = $1;
            $hapBases{$3} = $2;

            $totReads += $1   if ($3 ne "unknown");
            $totBases += $2   if ($3 ne "unknown");
        }
        if (m/(\d+)\s+reads\s+(\d+)\s+bases\s+filtered\s+for\s+being\s+too\s+short/) {
            $hapReads{"short"} = $1;
            $hapBases{"short"} = $2;
        }
    }
    close(F);

    print STDERR "--\n";
    foreach my $haplotype (@haplotypes) {
        printf STDERR "-- Found   %8d reads and %12d bases for haplotype $haplotype.\n", $hapReads{$haplotype}, $hapBases{$haplotype};
    }
    printf STDERR "-- Found   %8d reads and %12d bases assigned to no haplotype.\n", $hapReads{"unknown"}, $hapBases{"unknown"};
    printf STDERR "-- Ignored %8d reads and %12d bases because they were short.\n",  $hapReads{"short"},   $hapBases{"short"};

    #  Ignore the unknown reads if there aren't that many.

    my $unknownFraction =  getGlobal("hapUnknownFraction");
    my $withUnknown = (($totBases > 0) && ($hapBases{"unknown"} / $totBases < $unknownFraction)) ? 0 : 1;

    if ($withUnknown == 0) {
        print STDERR "--\n";
        print STDERR "-- Fewer than " . $unknownFraction*100 . " % of bases in unassigned reads; don't use them in assemblies.\n";
    } else {
        print STDERR "--\n";
        print STDERR "-- More than " .  $unknownFraction*100 . " % of bases in unassigned reads; including them in assemblies.\n";
    }

    #  For each haplotype, emit a script to run the assembly.

    print STDERR "--\n";
    print STDERR "-- Haplotype assembly commands:\n";

    foreach my $haplotype (@haplotypes) {
        my $hs = substr("$haplotype" . " " x $displLen, 0, $displLen);

        print STDERR "--   $rootdir/$asm-haplotype$haplotype.sh\n";

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

    if ($totBases == 0) {
        print STDERR "--\n";
        print STDERR "-- ERROR:\n";
        print STDERR "-- ERROR:  No reads assigned to haplotypes.  Assemblies not started.\n";
        print STDERR "-- ERROR:\n";
    }

    elsif ($hapBases{"unknown"} / $totBases > 0.50) {
        print STDERR "--\n";
        print STDERR "-- ERROR:\n";
        print STDERR "-- ERROR:  Too many bases in unassigned reads.  Assemblies not started.\n";
        print STDERR "-- ERROR:\n";
        print STDERR "-- ERROR:  If you run them manually, note that the unassigned reads\n";
        print STDERR "-- ERROR:  are included in ALL assemblies.\n";
        print STDERR "-- ERROR:\n";
    }

    #  Or stop if we're not running assemblies.

    elsif (setOptions($mode, "run") ne "run") {
        print STDERR "--\n";
        print STDERR "-- Assemblies not started per '-haplotype' option.\n";
    }

    #  Or run the assemblies.

    else {
        print STDERR "--\n";

        foreach my $haplotype (@haplotypes) {
            if (getGlobal("useGrid") ne "1") {
                print STDERR "-- Starting haplotype assembly $haplotype.\n";
            } else {
                print STDERR "-- Submitting haplotype assembly $haplotype.\n";
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
    }

    #  And now we're done.

    print STDERR "--\n";
    print STDERR "-- Bye.\n";

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

if (setOptions($mode, "trim") eq "trim") {
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
}

#  User-supplied termination command.

if (defined(getGlobal("onSuccess"))) {
    print STDERR "--\n";
    print STDERR "-- Running user-supplied termination command.\n";
    runCommand(getGlobal("onExitDir"), getGlobal("onSuccess") . " $asm");
}


print STDERR "--\n";
print STDERR "-- Bye.\n";

exit(0);
