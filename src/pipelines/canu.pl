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
use canu::Grid_Local;
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
        setGlobal("corOverlapper",  "mhap");
        setGlobal("obtOverlapper",  "mhap");
        setGlobal("utgOverlapper",  "mhap");
        setGlobal("utgReAlign",     "true");   #  Enabled by default, too.
    }

    if ($arg eq "-accurate") {
        #  All defaults, unless noted.
    }
}


#  By default, all three steps are run.  Options -correct, -trim and -assemble
#  can limit the pipeline to just that stage.

my $rootdir            = undef;
my $readdir            = undef;

my $mode               = undef;   #  "haplotype", "correct", "trim", "trim-assemble" or "assemble"
my $step               = "run";   #  Step to start at (?)

#  The filesAre variable decides what we should do with files.  It remembers
#  the last '-pacbio' et al. option given, and applies that to any file later
#  supplied on the command line.  Note that:
#    -pacbio -haplotypeA A.fasta pacbio.fasta
#  will treat the two files as belonging to haplotype "A" and result with no
#  files for assembly.  A more conventional parsing would notice that there
#  are no read files given for the -pacbio option and fail.

my $filesAre           = "unknown"; #  Treat files as 'pacbio', 'nanopore', 'pacbio-hifi' or 'haplotype' reads.
my $filesAreHap        = undef;     #  Treat files as haplotype short reads.

#  If no seqStore exists, we depend on the user setting command line options
#  to tell us the status of the reads.  The use of four readsAre variables is
#  so we can detect invalid cases (-raw -corrected -pacbio) instead of just
#  using the last one supplied.
#
#  If no options are supplied (canu -pacbio file.fasta) we'll later default
#  to 'raw' and 'untrimmed'.

my $readsAreRaw        = 0;         #  They're either raw or corrected.
my $readsAreCorrected  = 0;         #    If neither is set, we'll set to raw later.
my $readsAreUntrimmed  = 0;         #  They're either trimmed or not.
my $readsAreTrimmed    = 0;         #    If neither is set, we'll set to untrimmed later.

while (scalar(@ARGV)) {
    my $arg  = shift @ARGV;
    my $file = $arg;

    #  Decide if this argument is a file of reads.
    #   - Append any -readdir path.
    #   - Convert to an absolute path.
    #   - Unset if the path doesn't exist.
    #   - But then re-set if it looks like a DNA Nexus link (dnanexus:file-Ac36P534JbZvV2Gd1979x5Qv=reads.fasta.gz)
    #   - Complain if it's a file we can't support.

    $file = "$readdir/$arg"   if (defined($readdir));

    $file = abs_path($file)   if (  -e $file);
    $file = undef             if (! -e $file);

    $file = "$arg"            if ($arg =~ m/^dnanexus:.*=.*$/);

    if ((-e $file) && ($file =~ m/bam$/)) {
        addCommandLineError("ERROR: BAM input not supported: file '$arg'.\n");
    }

    #  Now just run through all the valid options.
    #
    #  addCommandLineOption() is just saving a (slightly modified) copy of
    #  the command line so we can resubmit ourself to the grid.

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
    #   - 'trimmed' reads have historically implied they are also corrected.

    elsif ($arg eq "-raw") {
        $readsAreRaw = 1;
        addCommandLineOption($arg);
    }

    elsif ($arg eq "-corrected") {
        $readsAreCorrected = 1;
        addCommandLineOption($arg);
    }

    elsif ($arg eq "-untrimmed") {
        $readsAreUntrimmed = 1;
        addCommandLineOption($arg);
    }

    elsif ($arg eq "-trimmed") {
        $readsAreCorrected = 1;
        $readsAreTrimmed = 1;
        addCommandLineOption($arg);
    }

    #  Set the technology of read we're getting.  This one block handles both
    #  the new style ('-pacbio') and old style ('-pacbio-raw') options.  An
    #  unfortunate side effect is that '-pacbio-hifi-raw' is allowed.
    #
    #  Search for the next (or around the next) instance of readsAreRaw to
    #  see where the rest of the settings are made.  There we set trimmed
    #  status and make sure that both raw and corrected aren't set.

    elsif ($arg =~ m/^-(pacbio|nanopore|pacbio-hifi)(-raw|-corrected){0,1}$/) {
        $filesAre          = $1;
        $readsAreRaw       = 1  if ($2 eq "-raw");
        $readsAreCorrected = 1  if ($2 eq "-corrected");
        $readsAreRaw       = 1  if ($1 eq "pacbio-hifi");

        addCommandLineOption($2)  if (defined($2));
    }

    elsif ($arg =~ m/^-haplotype(\w+)$/) {
        $filesAre    = "haplotype";
        $filesAreHap = $1;
    }

    #  File of reads.

    elsif (defined($file)) {
        if    ($filesAre eq "haplotype") {
            $haplotypeReads{$filesAreHap} .= "$file\0";
            addCommandLineOption("-haplotype$filesAreHap '$file'");
        }

        elsif ($filesAre ne "unknown") {
            push @inputFiles, "-$filesAre\0$file";
            addCommandLineOption("-$filesAre '$file'");
        }

        else {
            addCommandLineError("ERROR:  File '$arg' supplied on command line, don't know what to do with it.\n");
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

    open(F, "ls -d . *seqStore 2> /dev/null|");
    while (<F>) {
        $asm = $1   if (m/^(.*).seqStore$/);
    }
    close(F);
}

#  Fail if some obvious things aren't set.

addCommandLineError("ERROR:  Assembly name prefix (-p) not supplied.\n")   if (!defined($asm));
addCommandLineError("ERROR:  Assembly name prefix (-p) cannot contain the path delimiter '/'.\n")   if ($asm =~ m!/!);
addCommandLineError("ERROR:  Assembly name prefix (-p) cannot contain spaces.\n")   if ($asm =~ m!\s!);

#  Load paramters from the defaults files, then from any spec files on the
#  command line, then from the command line itself.

setParametersFromFile("$bin/canu.defaults")   if (-e "$bin/canu.defaults");
setParametersFromFile("$ENV{'HOME'}/.canu")   if (-e "$ENV{'HOME'}/.canu");

foreach my $specFile (@specFiles) {
    setParametersFromFile($specFile);
}

setParametersFromCommandLine(@specOpts);

checkParameters();
printHelp();           #  Fail now if any complaints so far.

#  OK so far, print citations and start run-time configuration by locating
#  extra binaries.

print STDERR "-- " . getGlobal("version") . "\n";
print STDERR "--\n";
print STDERR "-- CITATIONS\n";
print STDERR "--\n";
printCitation("-- ", $citationMode);
print STDERR "-- CONFIGURE CANU\n";
print STDERR "--\n";

checkJava();
checkMinimap($bin);
checkGnuplot();

checkParameters();
printHelp();           #  Fail now if any complaints so far.


################################################################################
#
#  Detect grid support.  If 'gridEngine' isn't set, the execution methods submitScript() and
#  submitOrRunParallelJob() will return without submitting, or run locally (respectively).  This
#  means that we can leave the default of 'useGrid' to 'true', and execution will do the right thing
#  when there isn't a grid.

detectLocal();
detectSGE();
detectSlurm();
detectPBSTorque();
detectLSF();
detectDNANexus();

#  Finish setting up the grid.  This is done AFTER parameters are set from the command line, to
#  let the user override any of our defaults.

configureLocal();      #  Also handles useGrid=false; needs to be after the detect() above.
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

checkParameters();
configureAssembler();  #  Set job sizes and etc bases on genomeSize and hosts available.

#  Make space for us to work in, and move there.
#
#  Unfortunately, if we fail later on because of inconsistent parameters
#  (e.g., -assemble -raw X.fasta) we end up with an empty output directory.
#  There is no easy way to get around this; fetchSeqStore() just below needs
#  the directory to exist before it can fetch, and we don't know enough yet
#  if the parameters are inconsistent.

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
    print STDERR "-- Found $rt reads in '$asm.seqStore':\n";
    print STDERR "--   Libraries:\n";
    print STDERR "--     PacBio CLR:            $numPacBio\n"     if ($numPacBio   > 0);
    print STDERR "--     Nanopore:              $numNanopore\n"   if ($numNanopore > 0);
    print STDERR "--     PacBio HiFi:           $numHiFi\n"       if ($numHiFi     > 0);
    print STDERR "--   Reads:\n";
    print STDERR "--     Raw:                   $numRaw\n"        if ($numRaw      > 0);
    print STDERR "--     Raw and Trimmed:       $numRawTri\n"     if ($numRawTri   > 0);
    print STDERR "--     Corrected:             $numCor\n"        if ($numCor      > 0);
    print STDERR "--     Corrected and Trimmed: $numCorTri\n"     if ($numCorTri   > 0);
    print STDERR "--\n";
}

#  Otherwise, scan input files, counting the different types of libraries we have.
#  At the same time, rewrite the sqStoreCreate option to accound for read status.

elsif (scalar(@inputFiles) > 0) {
    foreach my $tf (@inputFiles) {
        my ($t, $f) = split '\0', $tf;

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

    #  If no read type is set, default to 'raw' and 'untrimmed', unless
    #  they're HiFi reads, then they're 'raw' and 'trimmed'.

    if (($readsAreRaw == 0) && ($readsAreCorrected == 0)) {
        $readsAreRaw = 1;
    }

    if (($readsAreUntrimmed == 0) && ($readsAreTrimmed == 0)) {
        if ($numHiFi == 0) {
            $readsAreUntrimmed = 1;
        } else {
            $readsAreTrimmed = 1;
        }
    }

    #  If the user told us our HiFi reads are corrected (using something like
    #  "-corrected -pacbio-hifi"), silently undo that.  We need them to be
    #  called "raw" so they're stored in the store correctly.

    if (($readsAreRaw == 1) && ($readsAreCorrected == 1) && ($numHiFi > 0)) {
        $readsAreCorrected = 0;
    }

    #  Figure out a human description of the reads, and set
    #  flags to pass to sqStoreCreate.

    my $ct;   #  correction/trimmed status
    my $rt;   #  read tech
    my $st;   #  sqStore flags for the reads

    if      (($readsAreRaw == 0) && ($readsAreTrimmed == 0)) {
        $numCor = 1;
        $ct = "untrimmed corrected";
        $st = "-corrected";
    } elsif (($readsAreRaw == 0) && ($readsAreTrimmed == 1)) {
        $numCorTri = 1;
        $ct = "trimmed corrected";
        $st = "-corrected -trimmed";
    } elsif (($readsAreRaw == 1) && ($readsAreTrimmed == 0)) {
        $numRaw = 1;
        $ct = "untrimmed raw";
        $st = "-raw";
    } elsif (($readsAreRaw == 1) && ($readsAreTrimmed == 1)) {
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
    setGlobalIfUndef("corErrorRate",     0.300);
    setGlobalIfUndef("obtErrorRate",     0.120);
    setGlobalIfUndef("oeaErrorRate",     getGlobal("utgOvlErrorRate"));
    setGlobalIfUndef("oeaHaploConfirm",  5);
    setGlobalIfUndef("utgErrorRate",     0.120);
    setGlobalIfUndef("cnsErrorRate",     0.200);
}

if (($numPacBio    > 0) &&
    ($numNanopore >= 0) &&
    ($numHiFi     == 0)) {
    setGlobalIfUndef("corOvlErrorRate",  0.240);
    setGlobalIfUndef("obtOvlErrorRate",  0.045);
    setGlobalIfUndef("utgOvlErrorRate",  0.045);
    setGlobalIfUndef("corErrorRate",     0.250);
    setGlobalIfUndef("obtErrorRate",     0.045);
    setGlobalIfUndef("oeaErrorRate",     getGlobal("utgOvlErrorRate"));
    setGlobalIfUndef("oeaHaploConfirm",  5);
    setGlobalIfUndef("utgErrorRate",     0.045);
    setGlobalIfUndef("cnsErrorRate",     0.075);
}

if (($numPacBio   == 0) &&
    ($numNanopore == 0) &&
    ($numHiFi      > 0)) {
    setGlobalIfUndef("corOvlErrorRate",    0.000);
    setGlobalIfUndef("obtOvlErrorRate",    0.025);
    setGlobalIfUndef("utgOvlErrorRate",    0.010);
    setGlobalIfUndef("corErrorRate",       0.000);
    setGlobalIfUndef("obtErrorRate",       0.025);
    setGlobalIfUndef("utgErrorRate",       0.0003);
    setGlobalIfUndef("cnsErrorRate",       0.050);
    setGlobalIfUndef("homoPolyCompress",   1);
    setGlobalIfUndef("maxInputCoverage",   50);
    setGlobalIfUndef("utgGraphDeviation",  4);
    setGlobalIfUndef("utgRepeatDeviation", 0);
    setGlobalIfUndef("oeaErrorRate",       0.003);
    setGlobalIfUndef("oeaMaskTrivial",     1);
    setGlobalIfUndef("oeaHaploConfirm",    5);
}

if (($numPacBio > 0 || $numNanopore >0) && $numHiFi > 0) {
   caExit("ERROR: HiFi data cannot currently be combined with another read type", undef);
}

checkParameters();     #  Check all parameters (except error rates) are valid and consistent.
printHelp();           #  And one final last chance to fail.

################################################################################

#  Go!

#  Set an initial run mode based on what we discovered above.

#  SHOULD THIS BE SETTING THE MODE, OR LET IT DEFAULT TO RUN?
#    needs to skip correction for hifi
#
if (!defined($mode)) {
    $mode = "run"            if  ($numRaw    > 0);
    $mode = "trim-assemble"  if  ($numCor    > 0);
    $mode = "assemble"       if (($numHiFi   > 0) && ($readsAreTrimmed   == 1));
    $mode = "trim-assemble"  if (($numHiFi   > 0) && ($readsAreUntrimmed == 1));
    $mode = "assemble"       if  ($numCorTri > 0);
}

#  Print a message about what we're going to do.

printf STDERR "--\n";
printf STDERR "-- Generating assembly '$asm' in '" . getcwd() . "':\n";
printf STDERR "--   genomeSize:\n";
printf STDERR "--     %s\n", getGlobal("genomeSize");
printf STDERR "--\n";
printf STDERR "--   Overlap Generation Limits:\n";
printf STDERR "--     corOvlErrorRate %6.4f (%6.2f%%)\n", getGlobal("corOvlErrorRate"), getGlobal("corOvlErrorRate") * 100.0;
printf STDERR "--     obtOvlErrorRate %6.4f (%6.2f%%)\n", getGlobal("obtOvlErrorRate"), getGlobal("obtOvlErrorRate") * 100.0;
printf STDERR "--     utgOvlErrorRate %6.4f (%6.2f%%)\n", getGlobal("utgOvlErrorRate"), getGlobal("utgOvlErrorRate") * 100.0;
printf STDERR "--\n";
printf STDERR "--   Overlap Processing Limits:\n";
printf STDERR "--     corErrorRate    %6.4f (%6.2f%%)\n", getGlobal("corErrorRate"), getGlobal("corErrorRate") * 100.0;
printf STDERR "--     obtErrorRate    %6.4f (%6.2f%%)\n", getGlobal("obtErrorRate"), getGlobal("obtErrorRate") * 100.0;
printf STDERR "--     utgErrorRate    %6.4f (%6.2f%%)\n", getGlobal("utgErrorRate"), getGlobal("utgErrorRate") * 100.0;
printf STDERR "--     cnsErrorRate    %6.4f (%6.2f%%)\n", getGlobal("cnsErrorRate"), getGlobal("cnsErrorRate") * 100.0;
print  STDERR "--\n";
print  STDERR "--   Stages to run:\n";
print  STDERR "--     separate reads into haplotypes.\n"         if (($mode eq "run") && (scalar(keys %haplotypeReads) > 0));
print  STDERR "--     correct raw reads.\n"                      if (($mode eq "run"));
print  STDERR "--     trim corrected reads.\n"                   if (($mode eq "run"));
print  STDERR "--     assemble corrected and trimmed reads.\n"   if (($mode eq "run"));
print  STDERR "--     only separate reads into haplotypes.\n"    if (($mode eq "haplotype"));
print  STDERR "--     only correct raw reads.\n"                 if (($mode eq "correct"));
print  STDERR "--     only trim corrected reads.\n"              if (($mode eq "trim"));
print  STDERR "--     trim corrected reads.\n"                   if (($mode eq "trim-assemble"));
print  STDERR "--     assemble corrected and trimmed reads.\n"   if (($mode eq "trim-assemble"));
print  STDERR "--     assemble HiFi reads.\n"                    if (($mode eq "assemble") && ($numHiFi  > 0));
print  STDERR "--     assemble corrected and trimmed reads.\n"   if (($mode eq "assemble") && ($numHiFi == 0));
print  STDERR "--\n";

#  Fail with an error if the options don't make sense or if the types are inconsistent.

if (! -e "./$asm.seqStore/info") {
    my $inconsistent = 0;

    if (($readsAreRaw == 1) && ($numHiFi == 0) &&
        (($mode eq "trim") || ($mode eq "trim-assemble") || ($mode eq "assemble"))) {
        print STDERR "--\n";
        print STDERR "-- ERROR:\n";
        print STDERR "-- ERROR:  Cannot trim or assemmble raw (uncorrected) reads.\n";
        print STDERR "-- ERROR:   - If the reads are corrected, add '-corrected'.\n";
        print STDERR "-- ERROR:   - If the reads are raw, remove '-trim' or '-trim-assemble'.\n";
        print STDERR "-- ERROR:\n";
        print STDERR "-- ERROR:  readsAreRaw=true\n"         if ($readsAreRaw);
        print STDERR "-- ERROR:  readsAreCorrected=true\n"   if ($readsAreCorrected);
        print STDERR "-- ERROR:  mode=$mode\n";
        print STDERR "-- ERROR:\n";
        $inconsistent = 1;
    }

    if ((($readsAreCorrected == 1) || ($readsAreTrimmed == 1)) &&
        (($mode eq "correct"))) {
        print STDERR "--\n";
        print STDERR "-- ERROR:\n";
        print STDERR "-- ERROR:  Cannot correct already corrected or already trimmed reads.\n";
        print STDERR "-- ERROR:   - If the reads are raw (uncorrected), or to run\n";
        print STDERR "-- ERROR:     correction again, remove '-corrected'.\n";
        print STDERR "-- ERROR:   - Otherwise, remove the '-correct' option.\n";
        print STDERR "-- ERROR:\n";
        $inconsistent = 1;
    }

    if (($readsAreTrimmed == 1) &&
        (($mode eq "trim"))) {
        print STDERR "--\n";
        print STDERR "-- ERROR:\n";
        print STDERR "-- ERROR:  Cannot trim already trimmed reads.\n";
        print STDERR "-- ERROR:   - If the reads are untrimmed, or to run trimming again (*)n";
        print STDERR "-- ERROR:     replace '-trimmed' with '-corrected'.\n";
        print STDERR "-- ERROR:       * - Double trimming probably doesn't do much.\n";
        print STDERR "-- ERROR:   - Otherwise, remove the '-trim' option.\n";
        print STDERR "-- ERROR:\n";
        $inconsistent = 1;
    }

    if ($inconsistent) {
        caExit("inconsistent action for read types supplied", undef);
    }

    if (($readsAreRaw       == 1) &&
        ($readsAreCorrected == 1)) {
        print STDERR "--\n";
        print STDERR "-- ERROR:\n";
        print STDERR "-- ERROR:  Reads specified as 'raw' and 'corrected'.\n";
        print STDERR "-- ERROR:  Reads must be either all raw or all corrected.\n";
        print STDERR "-- ERROR:\n";
        $inconsistent = 1;
    }

    if (($readsAreRaw     == 1) &&
        ($readsAreTrimmed == 1) &&
        ($numHiFi         == 0)) {
        print STDERR "--\n";
        print STDERR "-- ERROR:\n";
        print STDERR "-- ERROR:  Reads specified as 'raw' and 'trimmed'.\n";
        print STDERR "-- ERROR:  Canu doesn't support this.  Remove the -trimmed option.\n";
        print STDERR "-- ERROR:  To not trim after correction, run correction separately,\n";
        print STDERR "-- ERROR:  then assemble.\n";
        print STDERR "-- ERROR:\n";
        $inconsistent = 1;
    }

    if (($readsAreUntrimmed == 1) &&
        ($readsAreTrimmed   == 1)) {
        print STDERR "--\n";
        print STDERR "-- ERROR:\n";
        print STDERR "-- ERROR:  Reads specified as both 'untrimmed' and 'trimmed'.\n";
        print STDERR "-- ERROR:  Please submit an issue for this; it shouldn't happen.\n";
        print STDERR "-- ERROR:  Try removing the -trimmed option, if present.\n";
        print STDERR "-- ERROR:\n";
        $inconsistent = 1;
    }

    if ($inconsistent) {
        caExit("inconsistent read types", undef);
    }
}

#  Make space for logs, and tell binaries where to write their execution
#  logging, then dump the parameters given to canu.

make_path("canu-logs")     if (! -d "canu-logs");
make_path("canu-scripts")  if (! -d "canu-scripts");

$ENV{'CANU_DIRECTORY'} = getcwd();

writeLog();


#  Stop if coverage is too low, unless no seqStore exists.
#
sub stopOnLowCoverage ($$) {
    my $asm = shift @_;
    my $tag = shift @_;

    return   if (! -e "./$asm.seqStore/info.txt");

    my $mincov = getGlobal("stopOnLowCoverage");
    my $curcov = getExpectedCoverage($asm, $tag);

    if ($curcov < $mincov) {
        print STDERR "--\n";
        print STDERR "-- ERROR:  Read coverage ($curcov) lower than allowed.\n";
        print STDERR "-- ERROR:    stopOnLowCoverage = $mincov\n";
        print STDERR "-- ERROR:\n";
        print STDERR "-- ERROR:  This could be caused by an incorrect genomeSize or poor\n";
        print STDERR "-- ERROR:  quality reads that cound not be sufficiently corrected.\n";
        print STDERR "-- ERROR:\n";
        print STDERR "-- ERROR:  You can force Canu to continue by decreasing parameter\n";
        print STDERR "-- ERROR:  stopOnLowCoverage (and possibly minInputCoverage too).\n";
        print STDERR "-- ERROR:  Be warned that the quality of corrected reads and/or\n";
        print STDERR "-- ERROR:  contiguity of contigs will be poor.\n";
        print STDERR "--\n";

        caExit("", undef);
    }
}



#  Decide if we need to enter the haplotyping pipeline
#   - no if there is no haplotype data supplied
#   - no if the mode isn't haplotype or run (and if it is run, the first 'no'
#     will obviously catch when there is no data to process
#   - no if the haplotyping is finsihed
#
sub doHaplotyping ($$@) {
    my $asm        = shift @_;
    my $mode       = shift @_;
    my @haplotypes =       @_;
    my $reason     = undef;

    #  Silently return if no haplotypes exist and/or it's just not enabled.

    return(0)   if (haplotypeReadsExist($asm, @haplotypes) eq "no-haplotypes");
    return(0)   if (($mode ne "haplotype") && ($mode ne "run"));

    $reason = "Haplotyping finished"   if (haplotypeReadsExist($asm, @haplotypes) eq "yes");

    if (defined($reason)) {
        print STDERR "--\n";
        print STDERR "-- $reason.\n";
        return(0);
    }

    submitScript($asm, undef);   #  See comments there as to why this is safe.

    print STDERR "--\n";
    print STDERR "-- BEGIN HAPLOTYPING\n";

    return(1);
}



#  Decide if we need to enter the correction pipeline.
#   - no if the mode tells us not to
#   - no if the output of correction is present
#   - no if there are corrected or trimmed reads available in the store
#
sub doCorrection ($$) {
    my $asm      = shift @_;
    my $mode     = shift @_;
    my $reason   = undef;

    my $cr = "corrected reads exist";

    my $exist1 = (fileExists("$asm.correctedReads.fasta.gz") ||
                  fileExists("$asm.correctedReads.fastq.gz"));
    my $exist2 = ((getNumberOfBasesInStore($asm, "obt") > 0) ||
                  (getNumberOfBasesInStore($asm, "utg") > 0));

    $reason = "no raw reads exist in $asm.seqStore"                 if ((fileExists("$asm.seqStore/info")) &&
                                                                        (getNumberOfBasesInStore($asm, "cor") == 0));
    $reason = "$cr in $asm.correctedReads.*.gz"                     if ($exist1);
    $reason = "$cr in $asm.seqStore"                                if ($exist2);
    $reason = "$cr in $asm.correctedReads.*.gz and $asm.seqStore"   if ($exist1 && $exist2);
    $reason = "not enabled"                                         if (($mode ne "correct") &&
                                                                        ($mode ne "run"));

    if (defined($reason)) {
        print STDERR "--\n";
        print STDERR "-- Correction skipped; $reason.\n";
        return(0);
    }

    stopOnLowCoverage($asm, "cor");
    submitScript($asm, undef);   #  See comments there as to why this is safe.

    print STDERR "--\n";
    print STDERR "-- BEGIN CORRECTION\n";

    return(1);
}



#  Decide if we need to enter the trimming pipeline.
#   - no if the mode tells us not to
#   - no if the output of trimming exists.
#   - no if there are trimmed reads available in the store
#
sub doTrimming ($$) {
    my $asm    = shift @_;
    my $mode   = shift @_;
    my $reason = undef;

    my $tr = "trimmed reads exist";

    my $exist1 = (fileExists("$asm.trimmedReads.fasta.gz") ||
                  fileExists("$asm.trimmedReads.fastq.gz"));
    my $exist2 = (getNumberOfBasesInStore($asm, "utg") > 0);

    $reason = "no corrected reads exist in $asm.seqStore"         if ((fileExists("$asm.seqStore/info")) &&
                                                                      (getNumberOfBasesInStore($asm, "obt") == 0));
    $reason = "$tr in $asm.trimmedReads.*.gz"                     if ($exist1);
    $reason = "$tr in $asm.seqStore"                              if ($exist2);
    $reason = "$tr in $asm.trimmedReads.*.gz and $asm.seqStore"   if ($exist1 && $exist2);
    $reason = "not enabled"                                       if (($mode ne "trim") &&
                                                                      ($mode ne "trim-assemble") &&
                                                                      ($mode ne "run"));

    if (defined($reason)) {
        print STDERR "--\n";
        print STDERR "-- Trimming skipped; $reason.\n";
        return(0);
    }

    stopOnLowCoverage($asm, "obt");
    submitScript($asm, undef);   #  See comments there as to why this is safe.

    print STDERR "--\n";
    print STDERR "-- BEGIN TRIMMING\n";
                    
    return(1);
}



#  Decide if we need to enter the unitigging pipeline.
#   - no if the mode tells us not to
#   - no if the output of trimming exists.
#   - no if there are trimmed reads available in the store
#
sub doUnitigging ($$) {
    my $asm    = shift @_;
    my $mode   = shift @_;
    my $reason = undef;

    $reason = "no corrected reads to assemble"        if (fileExists("$asm.seqStore/info") &&
                                                          (getNumberOfBasesInStore($asm, "obt") == 0) &&
                                                          (getNumberOfBasesInStore($asm, "utg") == 0));
    $reason = "contigs exist in $asm.contigs.fasta"   if (fileExists("$asm.contigs.fasta"));
    $reason = "not enabled"                           if (($mode ne "assemble") &&
                                                          ($mode ne "trim-assemble") &&
                                                          ($mode ne "run"));

    if (defined($reason)) {
        print STDERR "--\n";
        print STDERR "-- Unitigging skipped; $reason.\n";
        return(0);
    }

    stopOnLowCoverage($asm, "utg");
    submitScript($asm, undef);   #  See comments there as to why this is safe.

    print STDERR "--\n";
    print STDERR "-- BEGIN ASSEMBLY\n";

    return(1);
}



#  Run overlap jobs.
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
#  The start of haplotyping.  Separate reads into haplotypes, then maybe assemble.
#

my @haplotypes = sort keys %haplotypeReads;

if (doHaplotyping($asm, $mode, @haplotypes)) {
    haplotypeCountConfigure($asm, %haplotypeReads);

    haplotypeCountCheck($asm)                   foreach (1..getGlobal("canuIterationMax") + 1);
    haplotypeMergeCheck($asm, @haplotypes)      foreach (1..getGlobal("canuIterationMax") + 1);
    haplotypeSubtractCheck($asm, @haplotypes)   foreach (1..getGlobal("canuIterationMax") + 1);

    haplotypeReadsConfigure($asm, \@haplotypes, \@inputFiles);
    haplotypeReadsCheck($asm)                   foreach (1..getGlobal("canuIterationMax") + 1);
}

if (haplotypeReadsExist($asm, @haplotypes) eq "yes") {
    bootstrapHaplotypeAssemblies($asm, @haplotypes);

    #  If the mode lets us, start the assemblies.

    print STDERR "--\n";

    if ($mode eq "run") {
        foreach my $haplotype (@haplotypes) {
            print STDERR "-- Starting haplotype assembly $haplotype.\n"     if (getGlobal("useGrid") ne "1");
            print STDERR "-- Submitting haplotype assembly $haplotype.\n"   if (getGlobal("useGrid") eq "1");

            runCommand(".", "./$asm-haplotype$haplotype.sh");
        }

        print STDERR "--\n";
        print STDERR "-- Haplotype assemblies completed (or failed).\n"     if (getGlobal("useGrid") ne "1");
        print STDERR "-- Haplotype assemblies submitted.\n"                 if (getGlobal("useGrid") eq "1");
    } else {
        print STDERR "-- Not running assemblies.  Haplotyping finished.\n";
    }

    #  And now we're done.

    print STDERR "--\n";
    print STDERR "-- Bye.\n";

    exit(0);
}

#
#  The start of the pipeline.
#

if (doCorrection($asm, $mode)) {
    fetchSeqStore($asm);
    createSequenceStore($asm, @inputFiles);

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
    dumpCorrectedReads($asm);
}

if (doTrimming($asm, $mode)) {
    fetchSeqStore($asm);
    createSequenceStore($asm, @inputFiles);

    merylConfigure($asm, "obt");
    merylCountCheck($asm, "obt")     foreach (1..getGlobal("canuIterationMax") + 1);
    merylProcessCheck($asm, "obt")   foreach (1..getGlobal("canuIterationMax") + 1);

    overlap($asm, "obt");

    trimReads($asm);
    splitReads($asm);

    loadTrimmedReads($asm);
    dumpTrimmedReads($asm);
}

if (doUnitigging($asm, $mode)) {
    fetchSeqStore($asm);
    createSequenceStore($asm, @inputFiles);

    merylConfigure($asm, "utg");
    merylCountCheck($asm, "utg")       foreach (1..getGlobal("canuIterationMax") + 1);
    merylProcessCheck($asm, "utg")     foreach (1..getGlobal("canuIterationMax") + 1);

    overlap($asm, "utg");

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

    generateOutputs($asm);
}

stopAfter("theEnd");   #  Never returns.
exit(0);
