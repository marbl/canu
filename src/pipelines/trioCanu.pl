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
 #    src/pipelines/canu.pl
 #
 #  Modifications by:
 #
 #    Sergey Koren beginning on 2018-FEB-08
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #    Brian P. Walenz beginning on 2018-FEB-08
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

use strict;

use FindBin;
use Cwd qw(getcwd abs_path);
use POSIX qw(ceil);

use lib "$FindBin::RealBin/../lib/site_perl";

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
use canu::HaplotypeReads;

sub getHaplotypeReadInfo($$) {
   my $asm = shift @_;
   my $h   = shift @_;
   my $mean = 0;
   my $countH = 0;
   my $totalH = 0;

   open(F, "< $asm.${h}Reads.length") or caFailure("failed to read length information '$asm.${h}Reads.length'", undef);
   while (<F>) {
      my @l = split '\s+', $_;
      $countH++;
      $totalH += $l[-1];
   }
   close(F);
   $mean = sprintf "%.2f", ($totalH/$countH) if $countH > 0;
   return ($countH, $totalH, $mean);
}

my @specFiles;    #  Files of specs
my @specOpts;     #  Command line specs
my @inputFiles;   #  Command line inputs, later inputs in spec files are added
my %haplotypes;   #  Command line haplotype inputs

#  Initialize our defaults.  Must be done before defaults are reported in printOptions() below.

setDefaults();

#  The bin directory is needed for -version, can only be set after setDefaults(), but really should be
#  set after checkParameters() so it can know pathMap.

my $bin     = getBinDirectory();  #  Path to binaries, reset later.
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
}

#  At some pain, we stash the original options for later use.  We need
#  to use these when we resubmit ourself to the grid.  We can't simply dump
#  all of @ARGV into here, because we need to fix up relative paths first.

my $rootdir       = undef;
my $readdir       = undef;
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
        addCommandLineOption("-correct");

    } elsif ($arg eq "-trim") {
        addCommandLineOption("-trim");

    } elsif ($arg eq "-assemble") {
        addCommandLineOption("-assemble");

    } elsif ($arg eq "-trim-assemble") {
        addCommandLineOption("-trim-assemble");

    } elsif ($arg eq "-readdir") {
        addCommandLineOption("-readdir '$readdir'");

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

    } elsif ($arg =~ m/haplotype(\S*)$/) {
        my $hapID  = $1;

        my $file = $ARGV[0];
        my $fopt = addSequenceFile($readdir, $file, 1);

        while (defined($fopt)) {
            push @{ $haplotypes{$hapID} }, "-pacbio-raw\0$fopt";
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

#  Check for a few errors:
#    no mode                -> don't have any reads or any store to run from.
#    both raw and corrected -> don't know how to process these

caExit("ERROR: Can't mix uncorrected and corrected reads", undef)                     if ($haveRaw && $haveCorrected);

#  Go!
# we duplicate some parameters from assembly so we want to record them before our run and set them back at the end
#
my $asmReadLength=getGlobal("minReadLength");
my $genomesize = getGlobal("genomesize");

# compute kmer size given genome size and error rate
my $erate = 0.001;
my $kmer = log($genomesize * (1-$erate)/$erate) / log(4);
$kmer = int(ceil($kmer));

#initialize meryl params
setGlobal("hapOverlapper", "ovl");
setGlobal("hapovlmerthreshold", "auto");
setGlobal("hapovlmerdistinct", undef);
setGlobal("hapovlmertotal", undef);
setGlobal("hapovlfrequentmers", undef);
setGlobal("hapovlmertotal", undef);
setGlobal("hapovlMerSize", $kmer);

printf STDERR "--\n";
printf STDERR "-- Generating assembly '$asm' in '" . getcwd() . "'\n";
printf STDERR "--\n";
printf STDERR "-- Parameters:\n";
printf STDERR "--\n";
printf STDERR "--  genomeSize        %s\n", getGlobal("genomeSize");
printf STDERR "--  merSize           %s\n", getGlobal("hapovlMerSize");
printf STDERR "--\n";

#  Check that we were supplied a work directory, and that it exists, or we can create it.

make_path("canu-logs")     if (! -d "canu-logs");
make_path("canu-scripts")  if (! -d "canu-scripts");

#  This environment variable tells the binaries to log their execution in canu-logs/

$ENV{'CANU_DIRECTORY'} = getcwd();

#  Report the parameters used.

writeLog();

# check some params
caExit("ERROR: No reads supplied, and can't find any reads in any gkpStore", undef)   if (scalar(@inputFiles) == 0);
caExit("ERROR: Need at least two haplotypes", undef) if (scalar(keys %haplotypes) < 2);
foreach my $h (keys(%haplotypes)) {
   caExit("ERROR: No haplotype reads supplied for haplotype", undef) if (scalar(@{ $haplotypes{$h} }) == 0);
}

#  Submit ourself for grid execution?  If not grid enabled, or already running on the grid, this
#  call just returns.  The arg MUST be undef.
#
# no resuming for now so no grid
submitScript($asm, undef);

#
#  Begin pipeline
#
#  The checks for sequenceFileExists() at the start aren't needed except for
#  object storage mode.  Gatekeeper has no way of knowing, inside
#  gatekeeper(), that this stage is completed and it shouldn't fetch the
#  store.  In 'normal' operation, the store exists already, and we just
#  return.
#
print STDERR "--\n";
print STDERR "--\n";
print STDERR "-- BEGIN HAPLOTYPING\n";
print STDERR "--\n";

make_path("haplotype")  if (! -d "haplotype");

if (checkHaplotypeReads($asm, "haplotype") != 1) {
   # reset read length
   setGlobal("minReadLength", 50);

   # gkpStore for each haplotype and setup meryl
   gatekeeper($asm, "hap", @inputFiles);
   foreach my $h (keys(%haplotypes)) {
      gatekeeper("haplotype$h", "hap", @{ $haplotypes{$h} });
      merylConfigure("haplotype$h", "hap");
      merylCheck("haplotype$h", "hap")  foreach (1..getGlobal("canuIterationMax") + 1);
   }
   foreach my $h (keys(%haplotypes)) {
      merylSubtract("haplotype$h", "hap");
   }
   foreach my $h (keys(%haplotypes)) {
      merylFinishSubtraction("haplotype$h", "hap");
   }

   # now that we have haplotype information, classify the reads and dump them
   setGlobal("minReadLength", $asmReadLength);

   haplotypeConfigure($asm);
   haplotypeCheck($asm)      foreach (1..getGlobal("canuIterationMax") + 1);
   dumpHaplotypeReads($asm);
}

# now launch canu
#
# remove any read/haplotype options but keep anything else
my @commandOptions = split /\s+/, getCommandLineOptions();
my @fixedOptions;
my $setUpForPacBio = 0;
my $setUpForNanopore = 0;
my $readType = undef;
my $i = 0;
my $haveCorrected = 0;
my $haveRaw = 0;

while ($i < scalar(@commandOptions)) {
   my $option = $commandOptions[$i];

   # skip options we don't want
   if ($option =~ m/^-pacbio/) {
      if ($option =~ m/^-pacbio-raw/) {
         $haveRaw++;
      } else {
         $haveCorrected++;
      }

      if (!defined($readType)) {  # only set if not set first, if it gets set to nanopore, don't edit it
         $readType = $option;
      }
      $setUpForPacBio++;
      $i ++;
   } elsif ($option =~ m/^-nanopore/) {
      if ($option =~ m/^-nanopore-raw/) {
         $haveRaw++;
      } else {
         $haveCorrected++;
      }

      $readType = $option; # this one always wins
      $setUpForNanopore++;
      $i ++;
   } elsif ($option =~ m/^-haplotype/) {
      $i ++;
   } elsif ($option =~ m/^-d/) {
      $i ++;
   }
   else {
      push @fixedOptions, $option;
   }
   $i++;
}
caExit("ERROR: Can't mix uncorrected and corrected reads", undef)                     if ($haveRaw && $haveCorrected);

# report some stats
print STDERR "--\n";

my $totalBases = 0;
foreach my $h (keys(%haplotypes)) {
   my ($countH, $totalH, $mean) = getHaplotypeReadInfo($asm, "haplotype$h");
   $totalBases += $totalH;
   print STDERR "-- Haplotype $h has $countH sequences with $totalH bp (mean = $mean).\n";
}
my ($unknownC, $unknownB, $mean) = getHaplotypeReadInfo($asm, "unknown");
$totalBases += $unknownB;
my $unknownPercent = sprintf "%3.2f", ($unknownB / $totalBases * 100);
my $unknownReads = "";

# if we have almost no unclassified, don't bother with them
if ($unknownPercent <= 2) {
   $unknownReads = "";
} elsif ($unknownPercent > 50) {
   print STDERR "-- Haplotype unknown has $unknownPercent\% of total, not auto-assembling\n";
   caExit("Failed to haplotype reads, majority is unclassified", undef);
} else {
   $unknownReads = " $readType $asm.unknownReads.fasta.gz"
}
print STDERR "-- Haplotype unknown has $unknownC sequences with $unknownB bp (% of total = $unknownPercent\%, mean = $mean).\n";

foreach my $h (keys(%haplotypes)) {
   my $qcmd = "$bin/canu " . (join " ", @fixedOptions) .  " -d haplotype$h $readType $asm.haplotype${h}Reads.fasta.gz $unknownReads canuIteration=0 stopOnReadQuality=false\n";
   runCommand(getcwd(), $qcmd) and caFailure("Failed to run canu  on haplotype $h", undef);
}

print STDERR "--\n";
print STDERR "-- Bye.\n";

exit(0);
