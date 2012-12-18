#!/usr/bin/env perl

use strict;

use Config;  #  for @signame
use FindBin;
use Cwd;
use Carp;
use FileHandle;

use Sys::Hostname;

use POSIX qw(ceil floor sys_wait_h);

my $bin = undef;  #  Path to binaries, set once in main.
my $cmd = undef;  #  Temporary string passed to system().
my $wrk = undef;  #  Path to our assembly directory.
my $asm = undef;  #  Name of our assembly.

my $numFrags;
my %global;
my %synops;

my $commandLineOptions = "";
my $specLog            = "";

my @specFiles;
my @specOpts;
my @fragFiles;


sub submitBatchJobs($$) {
   my $SGE = shift @_;
   my $TAG = shift @_;

   if (runningOnGrid()) {
       runCommand($wrk, $SGE) and caFailure("Failed to submit batch jobs.");
       submitScript($TAG);
   } else {
       pleaseExecute($SGE);
   }
}


#  Decide what bin directory to use.
#
#  When we are running on SGE, the path of this perl script is NOT
#  always the correct architecture.  If the submission host is
#  FreeBSD, but the grid is Linux, the BSD box will submit
#  FreeBSD/bin/runCA.pl to the grid -- unless it knows in advance,
#  there is no way to pick the correct one.  The grid host then has to
#  have enough smarts to choose the correct binaries, and that is what
#  we're doing here.
#
#  To make it more trouble, shell scripts need to do all this by
#  themselves.
#
sub getBinDirectory () {
    my $installDir;

    ###
    ### CODE DUPLICATION WITH getBinDirectoryShellCode
    ###

    #  Assume the current binary path is the path to the global CA
    #  install directory.

    #  CODE DUPLICATION!!!
    my @t = split '/', "$FindBin::RealBin";
    pop @t;                      #  bin
    pop @t;                      #  arch, e.g., FreeBSD-amd64
    my $installDir = join '/', @t;  #  path to the assembler
    #  CODE DUPLICATION!!!

    #  Guess what platform we are currently running on.

    my $syst = `uname -s`;    chomp $syst;  #  OS implementation
    my $arch = `uname -m`;    chomp $arch;  #  Hardware platform
    my $name = `uname -n`;    chomp $name;  #  Name of the system

    $arch = "amd64"  if ($arch eq "x86_64");
    $arch = "ppc"    if ($arch eq "Power Macintosh");

    my $path = "$installDir/$syst-$arch/bin";

    my $pathMap = getGlobal("pathMap");
    if (defined($pathMap)) {
        open(F, "< $pathMap") or caFailure("failed to open pathMap '$pathMap'", undef);
        while (<F>) {
            my ($n, $b) = split '\s+', $_;
            $path = $b if ($name eq $n);
        }
        close(F);
    }

    return($path);
}

sub getBinDirectoryShellCode () {
    my $string;

    #  CODE DUPLICATION!!!
    my @t = split '/', "$FindBin::RealBin";
    pop @t;                      #  bin
    pop @t;                      #  arch, e.g., FreeBSD-amd64
    my $installDir = join '/', @t;  #  path to the assembler
    #  CODE DUPLICATION!!!

    $string  = "\n";
    $string .= "syst=`uname -s`\n";
    $string .= "arch=`uname -m`\n";
    $string .= "name=`uname -n`\n";
    $string .= "\n";
    $string .= "if [ \"\$arch\" = \"x86_64\" ] ; then\n";
    $string .= "  arch=\"amd64\"\n";
    $string .= "fi\n";
    $string .= "if [ \"\$arch\" = \"Power Macintosh\" ] ; then\n";
    $string .= "  arch=\"ppc\"\n";
    $string .= "fi\n";
    $string .= "\n";
    $string .= "bin=\"$installDir/\$syst-\$arch/bin\"\n";
    $string .= "\n";

    my $pathMap = getGlobal("pathMap");
    if (defined($pathMap)) {
        open(PM, "< $pathMap") or caFailure("failed to open pathMap '$pathMap'", undef);
        while (<PM>) {
            my ($n, $b) = split '\s+', $_;
            $string .= "if [ \"\$name\" = \"$n\" ] ; then\n";
            $string .= "  bin=\"$b\"\n";
            $string .= "fi\n";
        }
        close(PM);
        $string .= "\n";
    }

    return($string);
}





#  Return the second argument, unless the first argument is found in
#  %global, in which case return that.
#
sub getGlobal ($) {
    my $var = shift @_;
    caFailure("script error -- $var has no defined value", undef) if (!exists($global{$var}));
    return($global{$var});
}

sub setGlobal ($$) {
    my $var = shift @_;
    my $val = shift @_;

    #  If no value, set the field to undefined, the default for many of the options.

    $val = undef  if ($val eq "");

    #  Handle special cases.

    if ($var eq "merSize") {
        setGlobal("obtMerSize", $val);
        setGlobal("ovlMerSize", $val);
        return;
    }

    if ($var eq "merThreshold") {
        setGlobal("obtMerThreshold", $val);
        setGlobal("ovlMerThreshold", $val);
        return;
    }

    if ($var eq "overlapper") {
        setGlobal("obtOverlapper", $val);
        setGlobal("ovlOverlapper", $val);
        return;
    }

    #  Update obsolete usage.

    if ($var eq "doOverlapTrimming") {
        print STDERR "WARNING:  option doOverlapTrimming deprecated.  Use doOverlapBasedTrimming in the future.\n";
        $var = "doOverlapBasedTrimming";
    }

    if ($var eq "ovlMemory") {
        print STDERR "ERROR:  the runCA option ovlMemory is not recognized.\n";
        print STDERR "This runCA functionality changed in CABOG version 7.\n";
        print STDERR "Here are suggestions for new runCA option values:\n";
        print STDERR " ovlHashBits=23\n ovlHashBlockLength= 30000000\n (replaces ovlMemory=2GB)\n";
        print STDERR " ovlHashBits=24\n ovlHashBlockLength=110000000\n (replaces ovlMemory=4GB)\n";
        print STDERR " ovlHashBits=25\n ovlHashBlockLength=180000000\n (replaces ovlMemory=8GB)\n";
	exit(1);
    }

    if ($var eq "ovlHashBlockSize") {
        print STDERR "ERROR:  the runCA option ovlHashBlockSize is not recognized.\n";
        print STDERR "This runCA functionality changed in CABOG version 7.\n";
        print STDERR "Try using ovlHashBlockLength instead.\n";
        exit(1);
    }

    #  Update aliases.

    $var = "doOverlapBasedTrimming"    if ($var eq "doOBT");
    $var = "doExtendClearRanges"       if ($var eq "doECR");

    #  If "help" exists, we're parsing command line options, and will catch this failure in
    #  printHelp().  Otherwise, this is an internal error, and we should bomb now.
    #
    if (!exists($global{$var})) {
        if (exists($global{"help"})) {
            setGlobal("help", getGlobal("help") . "'$var' is not a valid option; see 'runCA -options' for a list of valid options.\n");
        } else {
            caFailure("'$var' is not a valid option Global variable.", undef);
        }
    }

    $global{$var} = $val;
}

sub setDefaults () {

    #  The rules:
    #
    #  1) Before changing these defaults, read the (printed) documentation.
    #  2) After changing, update the documentation.
    #  3) Add new defaults in the correct section.
    #  4) Keep defaults in the same order as the documentation.
    #  5) UPDATE THE DOCUMENTATION.
    #

    #####  Special Case for runCA-dedupe.pl

    $global{"outputPrefix"}               = undef;
    $synops{"outputPrefix"}               = "FASTQ output file prefix";

    #####  General Configuration Options (aka miscellany)

    $global{"pathMap"}                     = undef;
    $synops{"pathMap"}                     = "File with a hostname to binary directory map";

    $global{"shell"}                       = "/bin/sh";
    $synops{"shell"}                       = "Command interpreter to use; sh-compatible (e.g., bash), NOT C-shell (csh or tcsh)";

    #####  Error Rates

    $global{"ovlErrorRate"}                = 0.06;
    $synops{"ovlErrorRate"}                = "Overlaps above this error rate are not computed";

    $global{"utgErrorRate"}                = 0.015;
    $synops{"utgErrorRate"}                = "Overlaps at or below this error rate are used to construct unitigs (BOG and UTG)";

    $global{"utgErrorLimit"}               = 2.5;
    $synops{"utgErrorLimit"}               = "Overlaps at or below this number of errors are used to construct unitigs (BOG and UTG)";

    $global{"utgGraphErrorRate"}           = 0.030;
    $synops{"utgGraphErrorRate"}           = "Overlaps at or below this error rate are used to construct unitigs (BOGART)";

    $global{"utgGraphErrorLimit"}          = 3.25;
    $synops{"utgGraphErrorLimit"}          = "Overlaps at or below this number of errors are used to construct unitigs (BOGART)";

    $global{"utgMergeErrorRate"}           = 0.045;
    $synops{"utgMergeErrorRate"}           = "Overlaps at or below this error rate are used to construct unitigs (BOGART)";

    $global{"utgMergeErrorLimit"}          = 5.25;
    $synops{"utgMergeErrorLimit"}          = "Overlaps at or below this number of errors are used to construct unitigs (BOGART)";

    $global{"cnsErrorRate"}                = 0.06;
    $synops{"cnsErrorRate"}                = "Consensus expects alignments at about this error rate";

    $global{"cgwErrorRate"}                = 0.10;
    $synops{"cgwErrorRate"}                = "Unitigs/Contigs are not merged if they align above this error rate";

    #####  Minimums

    $global{"frgMinLen"}                   = 64;
    $synops{"frgMinLen"}                   = "Fragments shorter than this length are not loaded into the assembler";

    $global{"ovlMinLen"}                   = 40;
    $synops{"ovlMinLen"}                   = "Overlaps shorter than this length are not computed";

    #####  Stopping conditions

    $global{"stopBefore"}                  = undef;
    $synops{"stopBefore"}                  = "Tell runCA when to halt execution";

    $global{"stopAfter"}                   = undef;
    $synops{"stopAfter"}                   = "Tell runCA when to halt execution";

    #####  Sun Grid Engine

    $global{"useGrid"}                     = 0;
    $synops{"useGrid"}                     = "Enable SGE globally";

    $global{"scriptOnGrid"}                = 0;
    $synops{"scriptOnGrid"}                = "Enable SGE for runCA (and unitigger, scaffolder, other sequential phases)";

    $global{"mbtOnGrid"}                   = 1;
    $synops{"mbtOnGrid"}                   = "Enable SGE for mer-based trimming computations";

    $global{"ovlOnGrid"}                   = 1;
    $synops{"ovlOnGrid"}                   = "Enable SGE for overlap computations";

    $global{"frgCorrOnGrid"}               = 0;
    $synops{"frgCorrOnGrid"}               = "Enable SGE for the fragment error correction";

    $global{"ovlCorrOnGrid"}               = 0;
    $synops{"ovlCorrOnGrid"}               = "Enable SGE for the overlap error correction";

    $global{"cnsOnGrid"}                   = 1;
    $synops{"cnsOnGrid"}                   = "Enable SGE for consensus";

    $global{"sge"}                         = undef;
    $synops{"sge"}                         = "SGE options applied to all SGE jobs";

    $global{"sgeName"}                     = undef;
    $synops{"sgeName"}                     = "SGE jobs name suffix";

    $global{"sgeScript"}                   = undef;
    $synops{"sgeScript"}                   = "SGE options applied to runCA jobs (and unitigger, scaffolder, other sequential phases)";

    $global{"sgeMerTrim"}                  = undef;
    $synops{"sgeMerTrim"}                  = "SGE options applied to mer-based trimming jobs";

    $global{"sgeOverlap"}                  = undef;
    $synops{"sgeOverlap"}                  = "SGE options applied to overlap computation jobs";

    $global{"sgeMerOverlapSeed"}           = undef;
    $synops{"sgeMerOverlapSeed"}           = "SGE options applied to mer overlap seed (overmerry) jobs";

    $global{"sgeMerOverlapExtend"}         = undef;
    $synops{"sgeMerOverlapExtend"}         = "SGE options applied to mer overlap extend (olap-from-seeds) jobs";

    $global{"sgeConsensus"}                = undef;
    $synops{"sgeConsensus"}                = "SGE options applied to consensus jobs";

    $global{"sgeFragmentCorrection"}       = undef;
    $synops{"sgeFragmentCorrection"}       = "SGE options applied to fragment error correction jobs";

    $global{"sgeOverlapCorrection"}        = undef;
    $synops{"sgeOverlapCorrection"}        = "SGE options applied to overlap error correction jobs";

    $global{"sgePropagateHold"}            = undef;
    $synops{"sgePropagateHold"}            = undef;  #  Internal option

    #####  Preoverlap

    $global{"gkpFixInsertSizes"}           = 1;
    $synops{"gkpFixInsertSizes"}           = "Update stddev to 0.10 * mean if it is too large";

    $global{"gkpAllowInefficientStorage"}  = 0;
    $synops{"gkpAllowInefficientStorage"}  = "Allow mis-ordered reads in gkpStore; storage is inefficient and memory consuming";

    #####  Vector Trimming

    $global{"vectorIntersect"}             = undef;
    $synops{"vectorIntersect"}             = "File of vector clear ranges";

    $global{"vectorTrimmer"}               = "ca";
    $synops{"vectorTrimmer"}               = "Use the CA default vector trimmer, or figaro";

    $global{"figaroFlags"}                 = "-T 30 -M 100 -E 500 -V f";
    $synops{"figaroFlags"}                 = "Options to the figaro vector trimmer";

    #####  Overlap Based Trimming

    $global{"doOverlapBasedTrimming"}      = 1;
    $synops{"doOverlapBasedTrimming"}      = "Enable the Overlap Based Trimming module (doOBT and doOverlapTrimming are aliases)";

    $global{"doDeDuplication"}             = 1;
    $synops{"doDeDuplication"}             = "Enable the OBT duplication detection and cleaning module for 454 reads, enabled automatically";

    $global{"doChimeraDetection"}          = "normal";
    $synops{"doChimeraDetection"}          = "Enable the OBT chimera detection and cleaning module; 'off', 'normal' or 'aggressive'";

    #####  Mer Based Trimming

    $global{"mbtBatchSize"}                = 1000000;
    $synops{"mbtBatchSize"}                = "Process this many fragments per merTrim batch";

    $global{"mbtThreads"}                  = 4;
    $synops{"mbtThreads"}                  = "Number of threads to use when computing mer-based trimming";

    $global{"mbtConcurrency"}              = 1;
    $synops{"mbtConcurrency"}              = "If not SGE, number of mer-based trimming processes to run at the same time";

    #####  Overlapper

    $global{"obtOverlapper"}               = "ovl";
    $synops{"obtOverlapper"}               = "Which overlap algorithm to use for OBT overlaps";

    $global{"ovlOverlapper"}               = "ovl";
    $synops{"ovlOverlapper"}               = "Which overlap algorithm to use for OVL (unitigger) overlaps";

    $global{"ovlStoreMemory"}              = 1024;
    $synops{"ovlStoreMemory"}              = "How much memory (MB) to use when constructing overlap stores";

    $global{"ovlThreads"}                  = 2;
    $synops{"ovlThreads"}                  = "Number of threads to use when computing overlaps";

    $global{"ovlConcurrency"}              = 1;
    $synops{"ovlConcurrency"}              = "If not SGE, number of overlapper processes to run at the same time";

    $global{"ovlHashBlockLength"}          = 100000000;
    $synops{"ovlHashBlockLength"}          = "Amount of sequence (bp) to load into the overlap hash table";

    $global{"ovlRefBlockSize"}             = 2000000;
    $synops{"ovlRefBlockSize"}             = "Number of fragments to search against the hash table per batch";

    $global{"ovlRefBlockLength"}           = 0;
    $synops{"ovlRefBlockLength"}           = "Amount of sequence (bp) to search against the hash table per batch";

    $global{"ovlHashBits"}                 = "22";
    $synops{"ovlHashBits"}                 = "Width of the kmer hash.  Width 22=1gb, 23=2gb, 24=4gb, 25=8gb.  Plus 10b per ovlHashBlockLength";

    $global{"ovlHashLoad"}                 = "0.75";
    $synops{"ovlHashLoad"}                 = "Maximum hash table load.  If set too high, table lookups are inefficent; if too low, search overhead dominates run time";

    $global{"ovlMerSize"}                  = 22;
    $synops{"ovlMerSize"}                  = "K-mer size for seeds in overlaps";

    $global{"ovlMerThreshold"}             = "auto";
    $synops{"ovlMerThreshold"}             = "K-mer frequency threshold; mers more frequent than this are ignored";

    $global{"ovlFrequentMers"}              = undef;
    $synops{"ovlFrequentMers"}              = "Do not seed overlaps with these kmers (fasta format)";

    $global{"obtMerSize"}                  = 22;
    $synops{"obtMerSize"}                  = "K-mer size";

    $global{"obtMerThreshold"}             = "auto";
    $synops{"obtMerThreshold"}             = "K-mer frequency threshold; mers more frequent than this are ignored";

    $global{"obtFrequentMers"}              = undef;
    $synops{"obtFrequentMers"}              = "Do not seed overlaps with these kmers (fasta format)";

    $global{"ovlHashLibrary"}              = "0";
    $synops{"ovlHashLibrary"}              = "Only load hash fragments from specified lib, 0 means all";

    $global{"ovlRefLibrary"}                = 0;
    $synops{"ovlRefLibrary"}                = "Only load ref fragments from specified lib, 0 means all";

    $global{"obtHashLibrary"}              = "0";
    $synops{"obtHashLibrary"}              = "Only load hash fragments from specified lib, 0 means all";

    $global{"obtRefLibrary"}                = 0;
    $synops{"obtRefLibrary"}                = "Only load ref fragments from specified lib, 0 means all";


    $global{"merCompression"}              = 1;
    $synops{"merCompression"}              = "K-mer size";

    $global{"merOverlapperThreads"}        = 2;
    $synops{"merOverlapperThreads"}        = "Number of threads to use for both mer overlapper seed finding and extension jobs";

    $global{"merOverlapperSeedBatchSize"}  = 100000;
    $synops{"merOverlapperSeedBatchSize"}  = "Number of fragments in a mer overlapper seed finding batch; directly affects memory usage";

    $global{"merOverlapperExtendBatchSize"}= 75000;
    $synops{"merOverlapperExtendBatchSize"}= "Number of fragments in a mer overlapper seed extension batch; directly affects memory usage";

    $global{"merOverlapperCorrelatedDiffs"}= 0;
    $synops{"merOverlapperCorrelatedDiffs"}= "EXPERIMENTAL!";

    $global{"merOverlapperSeedConcurrency"}= 1;
    $synops{"merOverlapperSeedConcurrency"}= "If not SGE, number of mer overlapper seed finding processes to run at the same time";

    $global{"merOverlapperExtendConcurrency"}= 1;
    $synops{"merOverlapperExtendConcurrency"}= "If not SGE, number of mer overlapper seed extension processes to run at the same time";

    $global{"umdOverlapperFlags"}          = "-use-uncleaned-reads -trim-error-rate 0.03 -max-minimizer-cutoff 150";
    $synops{"umdOverlapperFlags"}          = "Options for the UMD overlapper";

    $global{"saveOverlaps"}                = 0;
    $synops{"saveOverlaps"}                = "Save intermediate overlap files";

    #####  Mers

    $global{"merylMemory"}                 = 800;
    $synops{"merylMemory"}                 = "Amount of memory, in MB, to use for mer counting";

    $global{"merylThreads"}                = 1;
    $synops{"merylThreads"}                = "Number of threads to use for mer counting";

    #####  Fragment/Overlap Error Correction

    $global{"frgCorrBatchSize"}            = 200000;
    $synops{"frgCorrBatchSize"}            = "Number of fragments per fragment error detection batch, directly affects memory usage";

    $global{"doFragmentCorrection"}        = 1;
    $synops{"doFragmentCorrection"}        = "Do overlap error correction";

    $global{"frgCorrThreads"}              = 2;
    $synops{"frgCorrThreads"}              = "Number of threads to use while computing fragment errors";

    $global{"frgCorrConcurrency"}          = 1;
    $synops{"frgCorrConcurrency"}          = "If not SGE, number of fragment error detection processes to run at the same time";

    $global{"ovlCorrBatchSize"}            = 200000;
    $synops{"ovlCorrBatchSize"}            = "Number of fragments per overlap error correction batch";

    $global{"ovlCorrConcurrency"}          = 4;
    $synops{"ovlCorrConcurrency"}          = "If not SGE, number of overlap error correction processes to run at the same time";

    #####  Illumina MP Classification

    $global{"dncMPlibraries"}              = undef;
    $synops{"dncMPlibraries"}              = "List of library names to run the 'de novo classifier' on";

    $global{"dncBBlibraries"}              = undef;
    $synops{"dncBBlibraries"}              = "List of library names to use as the bacbone in the 'de novo classifier'";

    #####  Unitigger & BOG & bogart Options

    $global{"unitigger"}                   = undef;
    $synops{"unitigger"}                   = "Which unitig algorithm to use; utg (if no SFF files) or bog (Best Overlap Graph, if SFF files)";

    $global{"utgGenomeSize"}               = 0;
    $synops{"utgGenomeSize"}               = "An estimate of the size of the genome; decides if unitigs are unique or repeats";

    $global{"utgBubblePopping"}            = 1;
    $synops{"utgBubblePopping"}            = "Smooth polymorphic regions";

    $global{"utgRecalibrateGAR"}           = 1;
    $synops{"utgRecalibrateGAR"}           = "Use an experimental algorithm to decide unique/repeat";

    $global{"bogBreakAtIntersections"}     = 1;
    $synops{"bogBreakAtIntersections"}     = "EXPERT!";

    $global{"bogBadMateDepth"}             = 7;
    $synops{"bogBadMateDepth"}             = "EXPERT!";

    $global{"batRebuildRepeats"}           = 0;
    $synops{"batRebuildRepeats"}           = "Shatter repeats, rebuild more stringently";

    $global{"batMateExtension"}            = 0;
    $synops{"batMateExtension"}            = "Shatter repeats, extend unique unitigs with mates (leaves repeats shattered)";

    $global{"batMemory"}                   = undef;
    $synops{"batMemory"}                   = "Approximate maximum memory usage for loading overlaps, in gigabytes, default is unlimited";

    #####  Scaffolder Options

    $global{"cgwPurgeCheckpoints"}         = 1;
    $synops{"cgwPurgeCheckpoints"}         = "Remove cgw checkpoint files when a scaffolding step finishes successfully";

    $global{"cgwDemoteRBP"}                = 1;
    $synops{"cgwDemoteRBP"}                = "EXPERT!";

    $global{"cgwUseUnitigOverlaps"}        = 0;
    $synops{"cgwUseUnitigOverlaps"}        = "Use unused best overlaps (from BOG) in scaffolder (EXPERIMENTAL)";

    $global{"astatLowBound"}               = 1;
    $synops{"astatLowBound"}               = "EXPERT!";

    $global{"astatHighBound"}              = 5;
    $synops{"astatHighBound"}              = "EXPERT!";

    $global{"stoneLevel"}                  = 2;
    $synops{"stoneLevel"}                  = "EXPERT!";

    $global{"computeInsertSize"}           = undef;
    $synops{"computeInsertSize"}           = "Compute a scratch scaffolding to estimate insert sizes; default: do only if less than 1 million reads";

    $global{"cgwDistanceSampleSize"}       = 100;
    $synops{"cgwDistanceSampleSize"}       = "Require N mates to reestimate insert sizes";

    $global{"doResolveSurrogates"}         = 1;
    $synops{"doResolveSurrogates"}         = "Place fragments in surrogates in the final assembly";

    $global{"doExtendClearRanges"}         = 2;
    $synops{"doExtendClearRanges"}         = "Enable the clear range extension heuristic";

    $global{"extendClearRangesStepSize"}   = undef;
    $synops{"extendClearRangesStepSize"}   = "Batch N scaffolds per ECR run";

    $global{"kickOutNonOvlContigs"}        = 0;
    $synops{"kickOutNonOvlContigs"}        = "Allow kicking out a contig placed in a scaffold by mate pairs that has no overlaps to both its left and right neighbor contigs. EXPERT!\n";

    $global{"doUnjiggleWhenMerging"}       = 0;
    $synops{"doUnjiggleWhenMerging"}       = "After inserting rocks/stones try shifting contig positions back to their original location when computing overlaps to see if they overlap with the rock/stone and allow them to merge if they do. EXPERT!\n";
    
    $global{"cgwContigShatterWeight"}      = 0;
    $synops{"cgwContigShatterWeight"}      = "When starting from a checkpoint, for any contig connected to its scaffold by a link with less than cgwContigShatterWeight, remove it and place it into a singleton scaffold. EXPERT!\n";

    $global{"cgwMergeMissingThreshold"}    = 0;
    $synops{"cgwMergeMissingThreshold"}    = "When merging scaffolds, missing mates are those mates that should fall within the merged scaffold but do not. In metagenomics, this may not be the case for a conserved region within strains as the mates are missing because they are in a different strain. This is a value between 0 and 1 to specify the percentage of missing mates (relative to good mates) to ignore. A value of -1 means ignore all missing mates when merging. EXPERT!\n";
    
    #####  Consensus Options

    $global{"cnsPartitions"}               = 128;
    $synops{"cnsPartitions"}               = "Partition consensus into N jobs";

    $global{"cnsMinFrags"}                 = 75000;
    $synops{"cnsMinFrags"}                 = "Don't make a consensus partition with fewer than N fragments";

    $global{"cnsConcurrency"}              = 2;
    $synops{"cnsConcurrency"}              = "If not SGE, number of consensus jobs to run at the same time";

    $global{"cnsPhasing"}                  = 0;
    $synops{"cnsPhasing"}                  = "Options for consensus phasing of SNPs\n\t0 - Do not phase SNPs to be consistent.\n\t1 - If two SNPs are joined by reads, phase them to be consistent.";

    $global{"consensus"}                   = "cns";
    $synops{"consensus"}                   = "Which consensus algorithm to use; currently only 'cns' is supported";

    #####  Terminator Options

    $global{"fakeUIDs"}                    = 0;
    $synops{"fakeUIDs"}                    = "Don't query a UID server, use UIDs specific to this assembly";

    $global{"uidServer"}                   = undef;
    $synops{"uidServer"}                   = "EXPERT!";

    $global{"createAGP"}                   = 0;
    $synops{"createAGP"}                   = "Create an AGP file for the assembly";

    $global{"createACE"}                   = 0;
    $synops{"createACE"}                   = "Create an ACE file for the assembly";

    $global{"createPosMap"}                = 1;
    $synops{"createPosMap"}                = "Create the POSMAP files for the assembly";

    $global{"merQC"}                       = 0;
    $synops{"merQC"}                       = "Compute a mer-based QC for the assembly";

    $global{"merQCmemory"}                 = 1024;
    $synops{"merQCmemory"}                 = "Memory to use for the mer-based QC";

    $global{"merQCmerSize"}                = 22;
    $synops{"merQCmerSize"}                = "Mer size to use for the mer-based QC";

    $global{"cleanup"}                     = "none";
    $synops{"cleanup"}                     = "At the end of a successful assembly, remove none/some/many/all of the intermediate files";

    #####  Options for toggling assembly. 

    $global{"doToggle"}                     = 0;
    $synops{"doToggle"}                     = "At the end of a successful assembly, search for placed surrogates and toggle them to be unique unitigs. Re-run the assembly starting from scaffolder";

    $global{"toggleUnitigLength"}           = 2000;
    $synops{"toggleUnitigLength"}           = "Minimum length for a surrogate to be toggled.";

    $global{"toggleNumInstances"}            = 1;
    $synops{"toggleNumInstances"}            = "Number of instances for a surrogate to be toggled. If 0 is specified, all non-singleton unitigs are toggled to unique status.";

    $global{"toggleMaxDistance"}             = 1000;
    $synops{"toggleMaxDistance"}             = "Toggling will look for surrogates that appear exactly twice, both at the end of a scaffold. This parameter specifies how close to the scaffold end the surrogate must be.";
    
    $global{"toggleDoNotDemote"}             = 0;
    $synops{"toggleDoNotDemote"}            = "Do not allow CGW to demote toggled unitigs based on branching patterns.";

    #### Closure Options

    $global{"closureOverlaps"}            = 0;
    $synops{"closureOverlaps"}             = "Option for handling overlaps involving closure reads.\n\t0 - Treat them just like regular reads, \n\t1 - Do not allow any overlaps (i.e. closure reads will stay as singletons until scaffolding), \n\t2 - allow overlaps betweeen closure reads and non-closure reads only";

    $global{"closurePlacement"}           = 2;
    $synops{"closurePlacement"}           = "Option for placing closure reads using the constraints.\n\t0 - Place at the first location found\n\t1 - Place at the best location (indicated by most constraints)\n\t2 - Place at multiple locations as long as the closure read/unitig in question is not unique";

    #####  Ugly, command line options passed to printHelp()

    $global{"help"}                        = "";
    $synops{"help"}                        = undef;

    $global{"version"}                     = 0;
    $synops{"version"}                     = undef;

    $global{"options"}                     = 0;
    $synops{"options"}                     = undef;



    if (exists($ENV{'AS_OVL_ERROR_RATE'})) {
        setGlobal("ovlErrorRate", $ENV{'AS_OVL_ERROR_RATE'});
        print STDERR "ENV: ovlErrorRate $ENV{'AS_OVL_ERROR_RATE'} (from \$AS_OVL_ERROR_RATE)\n";
    }

    if (exists($ENV{'AS_CGW_ERROR_RATE'})) {
        setGlobal("cgwErrorRate", $ENV{'AS_CGW_ERROR_RATE'});
        print STDERR "ENV: cgwErrorRate $ENV{'AS_CGW_ERROR_RATE'} (from \$AS_CGW_ERROR_RATE)\n";
    }

    if (exists($ENV{'AS_CNS_ERROR_RATE'})) {
        setGlobal("cnsErrorRate", $ENV{'AS_CNS_ERROR_RATE'});
        print STDERR "ENV: cnsErrorRate $ENV{'AS_CNS_ERROR_RATE'} (from \$AS_CNS_ERROR_RATE)\n";
    }

    if (exists($ENV{'AS_READ_MIN_LEN'})) {
        setGlobal("frgMinLen", $ENV{'AS_READ_MIN_LEN'});
        print STDERR "ENV: frgMinLen $ENV{'AS_READ_MIN_LEN'} (from \$AS_READ_MIN_LEN)\n";
    }

    if (exists($ENV{'AS_OVERLAP_MIN_LEN'})) {
        setGlobal("ovlMinLen", $ENV{'AS_OVERLAP_MIN_LEN'});
        print STDERR "ENV: ovlMinLen $ENV{'AS_OVERLAP_MIN_LEN'} (from \$AS_OVERLAP_MIN_LEN)\n";
    }
}

sub makeAbsolute ($) {
    my $var = shift @_;
    my $val = getGlobal($var);
    if (defined($val) && ($val !~ m!^/!)) {
        $val = "$ENV{'PWD'}/$val";
        setGlobal($var, $val);
        $commandLineOptions .= " \"$var=$val\" ";
    }
}

sub fixCase ($) {
    my $var = shift @_;
    my $val = getGlobal($var);
    $val =~ tr/A-Z/a-z/;
    setGlobal($var, $val);
}

sub setParametersFromFile ($@) {
    my $specFile  = shift @_;
    my @fragFiles = @_;

    #  Client should be ensuring that the file exists before calling this function.
    die "specFile '$specFile' not found.\n"  if (! -e "$specFile");

    $specLog .= "\n";
    $specLog .= "###\n";
    $specLog .= "###  Reading options from '$specFile'\n";
    $specLog .= "###\n";
    $specLog .= "\n";

    open(F, "< $specFile") or caFailure("Couldn't open '$specFile'", undef);

    while (<F>) {
        $specLog .= $_;

        s/^\s+//;
        s/\s+$//;

        next if (m/^\s*\#/);
        next if (m/^\s*$/);

        if (-e $_) {
            my $xx = $_;
            $xx = "$ENV{'PWD'}/$xx" if ($xx !~ m!^/!);
            if (-e $xx) {
                push @fragFiles, $xx;
            } else {
                setGlobal("help", getGlobal("help") . "File not found '$_' after appending absolute path.\n");
            }
        } elsif (m/\s*(\w*)\s*=([^#]*)#*.*$/) {
            my ($var, $val) = ($1, $2);
            $var =~ s/^\s+//; $var =~ s/\s+$//;
            $val =~ s/^\s+//; $val =~ s/\s+$//;
            undef $val if ($val eq "undef");
            setGlobal($var, $val);
        } else {
            setGlobal("help", getGlobal("help") . "File not found or unknown specFile option line '$_'.\n");
        }
    }
    close(F);

    return(@fragFiles);
}


sub setParametersFromCommandLine(@) {
    my @specOpts = @_;

    if (scalar(@specOpts) > 0) {
        $specLog .= "\n";
        $specLog .= "###\n";
        $specLog .= "###  Reading options from the command line.\n";
        $specLog .= "###\n";
        $specLog .= "\n";
    }

    foreach my $s (@specOpts) {
        $specLog .= "$s\n";

        if ($s =~ m/\s*(\w*)\s*=(.*)/) {
            my ($var, $val) = ($1, $2);
            $var =~ s/^\s+//; $var =~ s/\s+$//;
            $val =~ s/^\s+//; $val =~ s/\s+$//;
            setGlobal($var, $val);
        } else {
            setGlobal("help", getGlobal("help") . "Misformed command line option '$s'.\n");
        }
    }
}


sub setParameters () {

    #  Update obsolete usages.
    #
    if (getGlobal("doChimeraDetection") eq "1") {
        print STDERR "WARNING: 'doChimeraDetection=1' is obsolete; use 'doChimeraDetection=normal' in the future.\n";
        setGlobal("doChimeraDetection", "normal");
    }
    if (getGlobal("doChimeraDetection") eq "0") {
        print STDERR "WARNING: 'doChimeraDetection=0' is obsolete; use 'doChimeraDetection=off' in the future.\n";
        setGlobal("doChimeraDetection", "off");
    }

    #  Fiddle with filenames to make them absolute paths.
    #
    makeAbsolute("vectorIntersect");
    makeAbsolute("pathMap");

    #  Adjust case on some of them
    #
    fixCase("doChimeraDetection");
    fixCase("obtOverlapper");
    fixCase("ovlOverlapper");
    fixCase("unitigger");
    fixCase("vectorTrimmer");
    fixCase("stopBefore");
    fixCase("stopAfter");
    fixCase("consensus");
    fixCase("cleanup");

    if ((getGlobal("doChimeraDetection") ne "off") && (getGlobal("doChimeraDetection") ne "normal") && (getGlobal("doChimeraDetection") ne "aggressive")) {
        caFailure("invalid doChimeraDetection specified (" . getGlobal("doChimeraDetection") . "); must be 'off', 'normal', or 'aggressive'", undef);
    }
    if ((getGlobal("obtOverlapper") ne "mer") && (getGlobal("obtOverlapper") ne "ovl") && (getGlobal("obtOverlapper") ne "ovm")) {
        caFailure("invalid obtOverlapper specified (" . getGlobal("obtOverlapper") . "); must be 'mer' or 'ovl' (or DEVEL_ONLY 'ovm')", undef);
    }
    if ((getGlobal("ovlOverlapper") ne "mer") && (getGlobal("ovlOverlapper") ne "ovl") && (getGlobal("ovlOverlapper") ne "ovm")) {
        caFailure("invalid ovlOverlapper specified (" . getGlobal("ovlOverlapper") . "); must be 'mer' or 'ovl' (or DEVEL_ONLY 'ovm')", undef);
    }
    if (defined(getGlobal("unitigger")) && (getGlobal("unitigger") ne "utg") && (getGlobal("unitigger") ne "bog") && (getGlobal("unitigger") ne "bogart")) {
        caFailure("invalid unitigger specified (" . getGlobal("unitigger") . "); must be 'utg' or 'bog' or 'bogart'", undef);
    }
    if ((getGlobal("vectorTrimmer") ne "ca") && (getGlobal("vectorTrimmer") ne "figaro")) {
        caFailure("invalid vectorTrimmer specified (" . getGlobal("vectorTrimmer") . "); must be 'ca' or 'figaro'", undef);
    }
    if ((getGlobal("consensus") ne "cns") && (getGlobal("consensus") ne "seqan")) {
        caFailure("invalid consensus specified (" . getGlobal("consensus") . "); must be 'cns' or 'seqan'", undef);
    }
    if ((getGlobal("cnsPhasing") ne "0") && (getGlobal("cnsPhasing") ne "1")) {
       caFailure("invalid cnsPhasing specified (" . getGlobal("cnsPhasing") . "); must be '0' or '1'", undef);
    }
    if ((getGlobal("cleanup") ne "none") &&
        (getGlobal("cleanup") ne "light") &&
        (getGlobal("cleanup") ne "heavy") &&
        (getGlobal("cleanup") ne "aggressive")) {
        caFailure("invalid cleaup specified (" . getGlobal("cleanup") . "); must be 'none', 'light', 'heavy' or 'aggressive'", undef);
    }

    if (defined(getGlobal("stopBefore"))) {
        my $ok = 0;
        my $st = getGlobal("stopBefore");
        $st =~ tr/A-Z/a-z/;

        my $failureString = "Invalid stopBefore specified (" . getGlobal("stopBefore") . "); must be one of:\n";

        my @stopBefore = ("meryl",
                          "initialTrim",
                          "deDuplication",
                          "finalTrimming",
                          "chimeraDetection",
                          "classifyMates",
                          "unitigger",
                          "scaffolder",
                          "CGW",
                          "eCR",
                          "extendClearRanges",
                          "eCRPartition",
                          "extendClearRangesPartition",
                          "terminator");

        foreach my $sb (@stopBefore) {
            $failureString .= "    '$sb'\n";
            $sb =~ tr/A-Z/a-z/;
            if ($st eq $sb) {
                $ok++;
                setGlobal('stopBefore', $st);
            }
        }

        caFailure($failureString, undef) if ($ok == 0);
    }


    if (defined(getGlobal("stopAfter"))) {
        my $ok = 0;
        my $st = getGlobal("stopAfter");
        $st =~ tr/A-Z/a-z/;

        my $failureString = "Invalid stopAfter specified (" . getGlobal("stopAfter") . "); must be one of:\n";

        my @stopAfter = ("initialStoreBuilding",
                         "overlapper",
                         "OBT",
                         "overlapBasedTrimming",
                         "unitigger",
                         "classifyMates",
                         "utgcns",
                         "consensusAfterUnitigger",
                         "scaffolder",
                         "ctgcns",
                         "consensusAfterScaffolder");

        foreach my $sa (@stopAfter) {
            $failureString .= "    '$sa'\n";
            $sa =~ tr/A-Z/a-z/;
            if ($st eq $sa) {
                $ok++;
                setGlobal('stopAfter', $st);
            }
        }

        caFailure($failureString, undef) if ($ok == 0);
    }


    #  PIck a nice looking set of binaries, and check them.
    #
    {
        caFailure("can't find 'gatekeeper' program in $bin.  Possibly incomplete installation", undef) if (! -x "$bin/gatekeeper");
        caFailure("can't find 'meryl' program in $bin.  Possibly incomplete installation", undef)      if (! -x "$bin/meryl");
        caFailure("can't find 'overlap' program in $bin.  Possibly incomplete installation", undef)    if (! -x "$bin/overlapInCore");
        caFailure("can't find 'unitigger' program in $bin.  Possibly incomplete installation", undef)  if (! -x "$bin/unitigger");
        caFailure("can't find 'cgw' program in $bin.  Possibly incomplete installation", undef)        if (! -x "$bin/cgw");
        caFailure("can't find 'utgcns' program in $bin.  Possibly incomplete installation", undef)     if (! -x "$bin/utgcns");
        caFailure("can't find 'ctgcns' program in $bin.  Possibly incomplete installation", undef)     if (! -x "$bin/ctgcns");
        caFailure("can't find 'terminator' program in $bin.  Possibly incomplete installation", undef) if (! -x "$bin/terminator");

        if ((getGlobal("obtOverlapper") eq "mer") || (getGlobal("ovlOverlapper") eq "mer")) {
            caFailure("can't find 'overmerry' program in $bin.  Possibly incomplete installation", undef) if (! -x "$bin/overmerry");
        }
    }

    #  Set the globally accessible error rates.  Adjust them if they
    #  look strange.
    #
    #  We must have:     ovl <= cns <= cgw
    #  We usually have:  ovl == cns <= cgw
    #
    my $ovlER = getGlobal("ovlErrorRate");
    my $utgER = getGlobal("utgErrorRate");
    my $cgwER = getGlobal("cgwErrorRate");
    my $cnsER = getGlobal("cnsErrorRate");

    if (($ovlER < 0.0) || (0.25 < $ovlER)) {
        caFailure("ovlErrorRate is $ovlER, this MUST be between 0.00 and 0.25", undef);
    }
    if (($utgER < 0.0) || (0.25 < $utgER)) {
        caFailure("utgErrorRate is $utgER, this MUST be between 0.00 and 0.25", undef);
    }
    if (($cgwER < 0.0) || (0.25 < $cgwER)) {
        caFailure("cgwErrorRate is $cgwER, this MUST be between 0.00 and 0.25", undef);
    }
    if (($cnsER < 0.0) || (0.25 < $cnsER)) {
        caFailure("cnsErrorRate is $cnsER, this MUST be between 0.00 and 0.25", undef);
    }
    if ($utgER > $ovlER) {
        caFailure("utgErrorRate is $utgER, this MUST be <= ovlErrorRate ($ovlER)", undef);
    }
    if ($ovlER > $cnsER) {
        caFailure("ovlErrorRate is $ovlER, this MUST be <= cnsErrorRate ($cnsER)", undef);
    }
    if ($ovlER > $cgwER) {
        caFailure("ovlErrorRate is $ovlER, this MUST be <= cgwErrorRate ($cgwER)", undef);
    }
    if ($cnsER > $cgwER) {
        caFailure("cnsErrorRate is $cnsER, this MUST be <= cgwErrorRate ($cgwER)", undef);
    }

    $ENV{'AS_OVL_ERROR_RATE'} = $ovlER;
    $ENV{'AS_CGW_ERROR_RATE'} = $cgwER;
    $ENV{'AS_CNS_ERROR_RATE'} = $cnsER;

    $ENV{'AS_READ_MIN_LEN'}    = getGlobal("frgMinLen");
    $ENV{'AS_OVERLAP_MIN_LEN'} = getGlobal("ovlMinLen");
}

sub logVersion() {
        system("$bin/gatekeeper   --version");
        system("$bin/overlap      --version");
        system("$bin/unitigger    --version");
        system("$bin/buildUnitigs --version");
        system("$bin/cgw          --version");
        system("$bin/consensus    --version");
        system("$bin/terminator   --version");
}

sub printHelp () {

    if (getGlobal("version")) {
        logVersion();
        exit(0);
    }

    if (getGlobal("options")) {
        foreach my $k (sort keys %global) {
            my $o = substr("$k                                    ", 0, 35);
            my $d = substr(getGlobal($k) . "                      ", 0, 20);
            my $u = $synops{$k};

            if (!defined(getGlobal($k))) {
                $d = substr("<unset>                    ", 0, 20);
            }

            print "$o$d($u)\n";
        }
        exit(0);
    }

    if (getGlobal("help") ne "") {
        print "usage: runCA-dedupe -d <dir> -p <prefix> [options] <frg> ...\n";
        print "  -d <dir>          Use <dir> as the working directory.  Required\n";
        print "  -p <prefix>       Use <prefix> as the output prefix.  Required\n";
        print "\n";
        print "  -s <specFile>     Read options from the specifications file <specfile>.\n";
        print "                      <specfile> can also be one of the following key words:\n";
        print "                      [no]OBT - run with[out] OBT\n";
        print "                      noVec   - run with OBT but without Vector\n";
        print "\n";
        print "  -version          Version information\n";
        print "  -help             This information\n";
        print "  -options          Describe specFile options, and show default values\n";
        print "\n";
        print "  <frg>             CA formatted fragment file\n";
        print "\n";
        print "Complete documentation at http://wgs-assembler.sourceforge.net/\n";
        print "\n";
        print $global{"help"};
        exit(0);
    }

    undef $global{"version"};
    undef $global{"options"};
    undef $global{"help"};
}



sub checkDirectories () {

    #  Check that we were supplied a work directory, and that it
    #  exists, or we can create it.
    #
    die "ERROR: I need a directory to run the assembly in (-d option).\n" if (!defined($wrk));

    system("mkdir -p $wrk") if (! -d $wrk);
    chmod 0755, "$wrk";

    system("mkdir -p $wrk/runCA-logs") if (! -d "$wrk/runCA-logs");
    chmod 0755, "$wrk/runCA-logs";

    $ENV{'AS_RUNCA_DIRECTORY'} = $wrk;

    caFailure("directory '$wrk' doesn't exist (-d option) and couldn't be created", undef) if (! -d $wrk);
}


sub outputSpecLog () {
    my $time = time();
    my $host = hostname();
    my $pid  = $$;

    open(F, "> $wrk/runCA-logs/${time}_${host}_${pid}_runCA");
    print F $specLog;
    close(F);
}



sub findFirstCheckpoint ($) {
    my $dir      = shift @_;
    my $firstckp = 0;

    $dir = "$wrk/$dir" if (! -d $dir);

    open(F, "ls -1 $dir/$asm.ckp.[0-9]* |");
    while (<F>) {
        chomp;

        if (m/ckp.(\d+)$/) {
            $firstckp = $1 if ($1 < $firstckp);
        }
    }
    close(F);

    return($firstckp);
}

sub findLastCheckpoint ($) {
    my $dir     = shift @_;
    my $lastckp = 0;

    $dir = "$wrk/$dir" if (-d "$wrk/$dir");

    open(F, "ls -1 $dir/$asm.ckp.[0-9]* |");
    while (<F>) {
        chomp;

        if (m/ckp.(\d+)$/) {
            $lastckp = $1 if ($1 > $lastckp);
        }
    }
    close(F);

    return($lastckp);
}

sub findNumScaffoldsInCheckpoint ($$) {
    my $dir     = shift @_;
    my $lastckp = shift @_;

    open(F, "cd $wrk/$dir && $bin/getNumScaffolds ../$asm.gkpStore $asm $lastckp 2> /dev/null |");
    my $numscaf = <F>;  chomp $numscaf;
    close(F);
    $numscaf = int($numscaf);

    return($numscaf);
}


sub getNumberOfFragsInStore ($$) {
    my $wrk = shift @_;
    my $asm = shift @_;

    $numFrags = 0;

    if (-e "$wrk/$asm.gkpStore/inf") {
        open(F, "$bin/gatekeeper -lastfragiid $wrk/$asm.gkpStore 2> /dev/null |") or caFailure("failed to run gatekeeper to get the number of frags in the store", undef);
        $_ = <F>;    chomp $_;
        close(F);

        $numFrags = $1 if (m/^Last frag in store is iid = (\d+)$/);

        #  Check that the store is efficient

        my $typeMapSize = (-s "$wrk/$asm.gkpStore/inf") / 1024 / 1024;

        if (($typeMapSize > 128) && (getGlobal("gkpAllowInefficientStorage") == 0)) {
            $typeMapSize = int($typeMapSize * 3);

            my $err;
            $err .= "reads ordered inefficiently.  reorder (illumina first) or use\n";
            $err .= "gkpAllowInefficientStorage=1.  note that overlap jobs will use\n";
            $err .= "an additional $typeMapSize mb memory per job.";

            caFailure($err, "$wrk/$asm.gkpStore.err");
        }
    }

    return($numFrags);
}


#  Decide if we have the CA meryl or the Mighty one.
#
sub merylVersion () {
    my $ver = "unknown";

    open(F, "$bin/meryl -V |");
    while (<F>) {
        $ver = "CA"     if (m/CA/);
        $ver = "Mighty" if (m/Mighty/);
    }
    close(F);
    return($ver);
}

sub stopBefore ($$) {
    my $stopBefore = shift @_;  $stopBefore =~ tr/A-Z/a-z/;
    my $cmd        = shift @_;
    if (defined($stopBefore) &&
        defined(getGlobal('stopBefore')) &&
        (getGlobal('stopBefore') eq $stopBefore)) {
        print STDERR "Stop requested before '$stopBefore'.\n";
        print STDERR "Command:\n$cmd\n" if (defined($cmd));
        exit(0);
    }
}

sub stopAfter ($) {
    my $stopAfter = shift @_;  $stopAfter =~ tr/A-Z/a-z/;
    if (defined($stopAfter) &&
        defined(getGlobal('stopAfter')) &&
        (getGlobal('stopAfter') eq $stopAfter)) {
        print STDERR "Stop requested after '$stopAfter'.\n";
        exit(0);
    }
}

sub runningOnGrid () {
    return(defined($ENV{'SGE_TASK_ID'}));
}

sub findNextScriptOutputFile () {
    my $idx = "00";
    while (-e "$wrk/runCA.sge.out.$idx") {
        $idx++;
    }
    return("$wrk/runCA.sge.out.$idx");
}

sub submitScript ($) {
    my $waitTag = shift @_;

    return if (getGlobal("scriptOnGrid") == 0);

    my $output = findNextScriptOutputFile();
    my $script = "$output.sh";

    open(F, "> $script") or caFailure("failed to open '$script' for writing", undef);
    print F "#!" . getGlobal("shell") . "\n";
    print F "#\n";
    print F "#  Attempt to (re)configure SGE.  For reasons Bri doesn't know,\n";
    print F "#  jobs submitted to SGE, and running under SGE, fail to read his\n";
    print F "#  .tcshrc (or .bashrc, limited testing), and so they don't setup\n";
    print F "#  SGE (or ANY other paths, etc) properly.  For the record,\n";
    print F "#  interactive SGE logins (qlogin, etc) DO set the environment.\n";
    print F "\n";
    print F ". \$SGE_ROOT/\$SGE_CELL/common/settings.sh\n";
    print F "\n";
    print F "#  On the off chance that there is a pathMap, and the host we\n";
    print F "#  eventually get scheduled on doesn't see other hosts, we decide\n";
    print F "#  at run time where the binary is.\n";

    print F getBinDirectoryShellCode();

    print F "/usr/bin/env perl \$bin/runCA-dedupe $commandLineOptions\n";
    close(F);

    system("chmod +x $script");

    my $sge         = getGlobal("sge");
    my $sgeName     = getGlobal("sgeName");
    my $sgeScript   = getGlobal("sgeScript");
    my $sgePropHold = getGlobal("sgePropagateHold");

    $sgeName = "_$sgeName"              if (defined($sgeName));
    $waitTag = "-hold_jid \"$waitTag\"" if (defined($waitTag));

    my $qcmd = "qsub $sge $sgeScript -cwd -N \"rCA_$asm$sgeName\" -j y -o $output $waitTag $script";

    runCommand($wrk, $qcmd) and caFailure("Failed to submit script.\n");

    if (defined($sgePropHold)) {
        my $acmd = "qalter -hold_jid \"rCA_$asm$sgeName\" \"$sgePropHold\"";
        system($acmd) and print STDERR "WARNING: Failed to reset hold_jid trigger on '$sgePropHold'.\n";
    }

    exit(0);
}


sub caFailure ($$) {
    my  $msg = shift @_;
    my  $log = shift @_;

    print STDERR "================================================================================\n";
    print STDERR "\n";
    print STDERR "runCA-dedupe failed.\n";
    print STDERR "\n";

    print STDERR "----------------------------------------\n";
    print STDERR "Stack trace:\n";
    print STDERR "\n";
    carp;

    if (-e $log) {
        print STDERR "\n";
        print STDERR "----------------------------------------\n";
        print STDERR "Last few lines of the relevant log file ($log):\n";
        print STDERR "\n";
        system("tail -n 20 $log");
    }

    print STDERR "\n";
    print STDERR "----------------------------------------\n";
    print STDERR "Failure message:\n";
    print STDERR "\n";
    print STDERR "$msg\n";
    print STDERR "\n";

    exit(1);
}



#  Bit of a wierd one here; assume path are supplied relative to $wrk.
#  Potentially gives us a bit of safety.
#
sub rmrf (@) {
    foreach my $f (@_) {
        unlink("$wrk/$f")         if (-f "$wrk/$f");
        system("rm -rf $wrk/$f")  if (-d "$wrk/$f");
    }
}


#  Create an empty file.  Much faster than system("touch ...").
#
sub touch ($) {
    open(F, "> $_[0]") or caFailure("failed to touch file '$_[0]'", undef);
    close(F);
}


sub pleaseExecute ($) {
    my $file = shift @_;

    print STDERR "Please execute:\n";
    print STDERR "  $file\n";
    print STDERR "to submit jobs to the grid, then restart this script when all\n";
    print STDERR "jobs finish.  I'll make sure all jobs finished properly.\n";
}


#  Utility to run a command and check the exit status, report time used.
#
sub runCommand ($$) {
    my $dir = shift @_;
    my $cmd = shift @_;

    if (! -d $dir) {
        caFailure("Directory '$dir' doesn't exist, can't run command.\n", "");
    }

    my $t = localtime();
    my $d = time();
    print STDERR "----------------------------------------START $t\n$cmd\n";

    my $rc = 0xffff & system("cd $dir && $cmd");

    $t = localtime();
    print STDERR "----------------------------------------END $t (", time() - $d, " seconds)\n";

    #  Pretty much copied from Programming Perl page 230

    return(0) if ($rc == 0);

    #  Bunch of busy work to get the names of signals.  Is it really worth it?!
    #
    my @signame;
    if (defined($Config{sig_name})) {
        my $i = 0;
        foreach my $n (split('\s+', $Config{sig_name})) {
            $signame[$i] = $n;
            $i++;
        }
    }

    my $error = "ERROR: Failed with ";

    if ($rc == 0xff00) {
        $error .= "$!\n";
    } else {
        if ($rc & 0x80) {
            $error .= "coredump from ";
        }

        if ($rc > 0x80) {
            $rc >>= 8;
        }
        $rc &= 127;

        if (defined($signame[$rc])) {
            $error .= "signal $signame[$rc] ($rc)\n";
        } else {
            $error .= "signal $rc\n";
        }
    }

    print STDERR $error;

    return(1);
}

################################################################################
################################################################################
################################################################################

#  Functions for running multiple processes at the same time.

my $numberOfProcesses       = 0;     #  Number of jobs concurrently running
my $numberOfProcessesToWait = 0;     #  Number of jobs we can leave running at exit
my @processQueue            = ();
my @processesRunning        = ();
my $printProcessCommand     = 1;     #  Show commands as they run

sub schedulerSetNumberOfProcesses {
    $numberOfProcesses = shift @_;
}

sub schedulerSubmit {
    chomp @_;
    push @processQueue, @_;
}

sub schedulerForkProcess {
    my $process = shift @_;
    my $pid;

    #  From Programming Perl, page 167
  FORK: {
      if ($pid = fork) {
          # Parent
          #
          return($pid);
     } elsif (defined $pid) {
         # Child
         #
         exec($process);
      } elsif ($! =~ /No more processes/) {
          # EAGIN, supposedly a recoverable fork error
          sleep 1;
          redo FORK;
      } else {
          die "Can't fork: $!\n";
      }
  }
}

sub schedulerReapProcess {
    my $pid = shift @_;

    if (waitpid($pid, &WNOHANG) > 0) {
        return(1);
    } else {
        return(0);
    }
}

sub schedulerRun {
    my @newProcesses;

    #  Reap any processes that have finished
    #
    undef @newProcesses;
    foreach my $i (@processesRunning) {
        if (schedulerReapProcess($i) == 0) {
            push @newProcesses, $i;
        }
    }
    undef @processesRunning;
    @processesRunning = @newProcesses;

    #  Run processes in any available slots
    #
    while ((scalar(@processesRunning) < $numberOfProcesses) &&
           (scalar(@processQueue) > 0)) {
        my $process = shift @processQueue;
        print STDERR "$process\n";
        push @processesRunning, schedulerForkProcess($process);
    }
}

sub schedulerFinish {
    my $child;
    my @newProcesses;
    my $remain;

    my $t = localtime();
    my $d = time();
    print STDERR "----------------------------------------START CONCURRENT $t\n";

    $remain = scalar(@processQueue);

    #  Run all submitted jobs
    #
    while ($remain > 0) {
        schedulerRun();

        $remain = scalar(@processQueue);

        if ($remain > 0) {
            $child = waitpid -1, 0;

            undef @newProcesses;
            foreach my $i (@processesRunning) {
                push @newProcesses, $i if ($child != $i);
            }
            undef @processesRunning;
            @processesRunning = @newProcesses;
        }
    }

    #  Wait for them to finish, if requested
    #
    while (scalar(@processesRunning) > $numberOfProcessesToWait) {
        waitpid(shift @processesRunning, 0);
    }

    $t = localtime();
    print STDERR "----------------------------------------END CONCURRENT $t (", time() - $d, " seconds)\n";
}

################################################################################
################################################################################
################################################################################


sub preoverlap {
    my @fragFiles = @_;

    $numFrags = getNumberOfFragsInStore($wrk, $asm);

    #  Return if there are fragments in the store, and die if there
    #  are no fragments and no source files.
    #
    if ($numFrags > 0) {
        goto stopafter;
    }

    caFailure("no fragment files specified, and stores not already created", undef)
    	if (scalar(@fragFiles) == 0);

    if ((! -d "$wrk/$asm.gkpStore") ||
        (! -e "$wrk/$asm.gkpStore/inf")) {

        #  Make sure all the inputs are here.  We also shred any
        #  supplied ace files, and convert the sff's to frg's.
        #
        my $failedFiles = undef;
        my $gkpInput = "";
        foreach my $frg (@fragFiles) {
            if (! -e $frg) {
                if (defined($failedFiles)) {
                    $failedFiles .= "; '$frg' not found";
                } else {
                    $failedFiles = "'$frg' not found";
                }
            }

            if ($frg =~ m/^(.*)\.ace$/) {
                my @fff = split '/', $1;
                my $ace = $frg;
                my $nam = pop @fff;

                $frg = "$wrk/$nam.shred.frg";

                if (! -e "$frg") {
                    print STDERR "Shredding '$ace' -> '$frg'\n";
                    shredACE($ace, $frg);
                }
            }

            if (($frg =~ m/^(.*)\.sff$/) ||
                ($frg =~ m/^(.*)\.sff.gz$/) ||
                ($frg =~ m/^(.*)\.sff.bz2$/)) {
                my @fff = split '/', $1;
                my $sff = $frg;
                my $nam = pop @fff;
                my $log = "$wrk/$nam.sff.log";

                $frg = "$wrk/$nam.sff.frg";

                if (! -e "$frg") {
                    print STDERR "Converting '$sff' -> '$frg'\n";

                    if (runCommand($wrk, "$bin/sffToCA -libraryname $nam -linker flx -linker titanium -insertsize 3000 300 -output $frg $sff > $frg.err 2>&1")) {
                        unlink "$wrk/$frg";
                        caFailure("sffToCA failed", "$frg.err");
                    }
                }
            }

            $gkpInput .= " $frg";
        }
        caFailure($failedFiles, undef) if defined($failedFiles);

        $cmd  = "$bin/gatekeeper ";
        $cmd .= " -o $wrk/$asm.gkpStore.BUILDING ";
        $cmd .= " -T " if (getGlobal("doOverlapBasedTrimming"));
        $cmd .= " -F " if (getGlobal("gkpFixInsertSizes"));
        $cmd .= "$gkpInput ";
        $cmd .= "> $wrk/$asm.gkpStore.err 2>&1";

        if (runCommand($wrk, $cmd)) {
            caFailure("gatekeeper failed", "$wrk/$asm.gkpStore.err");
        }

        rename "$wrk/$asm.gkpStore.BUILDING",             "$wrk/$asm.gkpStore";
        rename "$wrk/$asm.gkpStore.BUILDING.errorLog",    "$wrk/$asm.gkpStore.errorLog";
        rename "$wrk/$asm.gkpStore.BUILDING.fastqUIDmap", "$wrk/$asm.gkpStore.fastqUIDmap";
    }

    $numFrags = getNumberOfFragsInStore($wrk, $asm);

    print STDERR "numFrags = $numFrags\n";

    caFailure("gatekeeper failed to add fragments", "$wrk/$asm.gkpStore.err")  if ($numFrags == 0);

  stopafter:
    stopAfter("initialStoreBuilding");
}

################################################################################
################################################################################
################################################################################


sub runMeryl ($$$$$$) {
    my $merSize      = shift @_;
    my $merComp      = shift @_;
    my $merCanonical = shift @_;
    my $merThresh    = shift @_;
    my $merScale     = 1.0;
    my $merType      = shift @_;
    my $merDump      = shift @_;

    system("mkdir $wrk/0-mercounts") if (! -d "$wrk/0-mercounts");

    #  The fasta file we should be creating.
    my $ffile = "$wrk/0-mercounts/$asm.nmers.$merType.fasta";

    if ($merThresh =~ m/auto\s*\*\s*(\S+)/) {
        $merThresh = "auto";
        $merScale  = $1;
    }

    if ($merThresh =~ m/auto\s*\/\s*(\S+)/) {
        $merThresh = "auto";
        $merScale  = 1.0 / $1;
    }

    if (($merThresh ne "auto") && ($merThresh == 0)) {
        touch $ffile;
        return;
    }

    #if (-e $ffile) {
    #    print STDERR "runMeryl() would have returned.\n";
    #}
    #print STDERR "$ffile\n";

    if (merylVersion() eq "Mighty") {

        #  Use the better meryl!  This is straightforward.  We count,
        #  then we dump.

        #  Intermediate file
        my $ofile = "$wrk/0-mercounts/$asm$merCanonical-ms$merSize-cm$merComp";

        if (! -e "$ofile.mcdat") {
            my $merylMemory  = getGlobal("merylMemory");
            my $merylThreads = getGlobal("merylThreads");

            if ($merylMemory !~ m/^-/) {
                $merylMemory = "-memory $merylMemory";
            }

            #  A small optimization we could do if (a) not mer
            #  overlapper, (b) not auto threshold: only save mer
            #  counts above the smaller (of obt & ovl thresholds).
            #  It's complicated, and potentially screws up restarts
            #  (if the threshold is changed after meryl is finished,
            #  for example).  It's only useful on large assemblies,
            #  which we usually assume you know what's going on
            #  anyway.
            #
            #  N.B. the mer overlapper NEEDS all mer counts 2 and
            #  higher.

            $cmd  = "$bin/meryl ";
            $cmd .= " -B $merCanonical -v -m $merSize $merylMemory -threads $merylThreads -c $merComp ";
            $cmd .= " -L 2 ";
            $cmd .= " -s $wrk/$asm.gkpStore:chain ";  #  Or 'latest' to get the original version
            $cmd .= " -o $ofile ";
            $cmd .= "> $wrk/0-mercounts/meryl.err 2>&1";

            stopBefore("meryl", $cmd);

            if (runCommand("$wrk/0-mercounts", $cmd)) {
                caFailure("meryl failed", "$wrk/0-mercounts/meryl.err");
            }
            unlink "$wrk/0-mercounts/meryl.err";
        }

        if ($merThresh eq "auto") {
            if (! -e "$ofile.estMerThresh.out") {
                $cmd  = "$bin/estimate-mer-threshold ";
                $cmd .= " -m $ofile ";
                $cmd .= " > $ofile.estMerThresh.out ";
                $cmd .= "2> $ofile.estMerThresh.err";

                stopBefore("meryl", $cmd);

                if (runCommand("$wrk/0-mercounts", $cmd)) {
                    rename "$ofile.estMerThresh.out", "$ofile.estMerThresh.out.FAILED";
                    caFailure("estimate-mer-threshold failed", "$ofile.estMerThresh.err");
                }
            }

            open(F, "< $ofile.estMerThresh.out") or caFailure("failed to read estimated mer threshold from '$ofile.estMerThresh.out'", undef);
            $merThresh = <F>;
            $merThresh = int($merThresh * $merScale);
            close(F);

            if ($merThresh == 0) {
                caFailure("failed to estimate a mer threshold", "$ofile.estMerThresh.err");
            }
        }

        #  We only need the ascii dump if we're doing overlapper, mer
        #  overlapper reads meryl directly.
        #
        if ($merDump) {
            if (! -e $ffile) {
                $cmd  = "$bin/meryl ";
                $cmd .= "-Dt -n $merThresh ";
                $cmd .= "-s $ofile ";
                $cmd .= "> $ffile ";
                $cmd .= "2> $ffile.err ";

                if (runCommand("$wrk/0-mercounts", $cmd)) {
                    unlink $ffile;
                    caFailure("meryl failed to dump frequent mers", "$ffile.err");
                }
                unlink "$ffile.err";
            }
        }
    } elsif (merylVersion() eq "CA") {

        #  Sigh.  The old meryl.  Not as easy.  If we assume the
        #  process, in particular that the Ovl threshold is less than
        #  the Obt threshold, and that we have already computed the
        #  Ovl mers, we could filter the Ovl mers to get the Obt mers.
        #  But that's tough, especially if we allow mer compression.

        my $merSkip = 10;

        #  Intermediate file
        my $ofile = "$wrk/0-mercounts/$asm-ms$merSize-mt$merThresh-mk$merSkip.$merType.fasta";

        if ($merComp > 0) {
            print STDERR "ERROR!  merCompression not supported without installing kmer\n";
            print STDERR "        (http://sourceforge.net/projects/kmer/).\n";
            print STDERR "If you have installed kmer, then your build is broken, as I\n";
            print STDERR "did not find the correct 'meryl' (meryl -V should have said Mighty).\n";
            die;
        }

        if ($merCanonical ne "-C") {
            print STDERR "ERROR!  mer overlapper not supported without installing kmer\n";
            print STDERR "        (http://sourceforge.net/projects/kmer/).\n";
            print STDERR "If you have installed kmer, then your build is broken, as I\n";
            print STDERR "did not find the correct 'meryl' (meryl -V should have said Mighty).\n";
            die;
        }

        if ($merThresh eq "auto") {
            print STDERR "WARNING!  auto picking a mer threshold not supported without installing kmer\n";
            print STDERR "          (http://sourceforge.net/projects/kmer/).\n";
            print STDERR "Using historical defaults.\n";

            if ($merType eq "obt") {
                $merThresh = 1000;
            } else {
                $merThresh = 500;
            }
        }

        if (! -e $ofile) {
            my $mt = $merThresh / $merSkip;

            $cmd  = "$bin/meryl ";
            $cmd .= "-s $wrk/$asm.gkpStore -m $merSize -n $mt -K $merSkip ";
            $cmd .= " -o $ofile";
            $cmd .= "> $wrk/0-mercounts/meryl.err 2>&1";

            stopBefore("meryl", $cmd);

            if (runCommand("$wrk/0-mercounts", $cmd)) {
                unlink $ofile;
                caFailure("meryl failed to dump frequent mers", "$wrk/0-mercounts/meryl.err");
            }
            unlink "$wrk/0-mercounts/meryl.err";
        }

        symlink($ofile, $ffile) if (! -e $ffile);
    } else {
        caFailure("unknown meryl version '" . merylVersion() . "'", "");
    }

    return($merThresh);
}

sub meryl {

    if (getGlobal("ovlOverlapper") eq "umd") {
        caFailure("meryl attempted to compute mer counts for the umd overlapper", undef);
    }

    my $ovlc = 0;  #  No compression, unless we're the mer overlapper
    my $obtc = 0;

    my $ovlC = "-C";  #  Canonical, unless we're the mer overlapper
    my $obtC = "-C";  #  (except the mer overlapper now wants canonical)

    my $ovlD = 1;  #  Dump, unless we're the mer overlapper
    my $obtD = 1;

    my $obtT = 0;  #  New threshold
    my $ovlT = 0;

    if (getGlobal("ovlOverlapper") eq "mer") {
        $ovlc = getGlobal("merCompression");
        $ovlC = "-C";
        $ovlD = 0;
    }
    if (getGlobal("obtOverlapper") eq "mer") {
        $obtc = getGlobal("merCompression");
        $obtC = "-C";
        $obtD = 0;
    }


    system("mkdir $wrk/0-mercounts") if (! -d "$wrk/0-mercounts");

    if (defined(getGlobal("obtFrequentMers"))) {
        my $ffile = "$wrk/0-mercounts/$asm.nmers.obt.fasta";
        my $sfile = getGlobal("obtFrequentMers");

        if (! -e $ffile) {
            caFailure("obtFrequentMers '$sfile' not found", undef)  if (! -e $sfile);
            print STDERR "Using obt frequent mers in '$sfile'\n";
            symlink $sfile, $ffile;
        }
    }


    if (defined(getGlobal("ovlFrequentMers"))) {
        my $ffile = "$wrk/0-mercounts/$asm.nmers.ovl.fasta";
        my $sfile = getGlobal("ovlFrequentMers");

        if (! -e $ffile) {
            caFailure("ovlFrequentMers '$sfile' not found", undef)  if (! -e $sfile);
            print STDERR "Using ovl frequent mers in '$sfile'\n";
            symlink $sfile, $ffile;
        }
    }

    #  We're only here to compute mercounts for the overlappers (all of mer, obt and ovl).  If the
    #  frequent mer file(s) exist, don't bother running meryl (which runs to build mcidx/mcdat even
    #  if the fasta is there).

    if ((getGlobal("doOverlapBasedTrimming")) &&
        ((getGlobal('obtOverlapper') eq 'mer') || (! -e "$wrk/0-mercounts/$asm.nmers.obt.fasta"))) {
        $obtT = runMeryl(getGlobal('obtMerSize'), $obtc, $obtC, getGlobal("obtMerThreshold"), "obt", $obtD) if (getGlobal("doOverlapBasedTrimming"));
    } else {
        print STDERR "No need to run meryl for OBT (OBT is disabled).\n"             if (getGlobal("doOverlapBasedTrimming") == 0);
        print STDERR "No need to run meryl for OBT ($asm.nmers.obt.fasta exists).\n" if (-e "$wrk/0-mercounts/$asm.nmers.obt.fasta");
    }

    if ((getGlobal('ovlOverlapper') eq 'mer') || (! -e "$wrk/0-mercounts/$asm.nmers.ovl.fasta")) {
        $ovlT = runMeryl(getGlobal('ovlMerSize'), $ovlc, $ovlC, getGlobal("ovlMerThreshold"), "ovl", $ovlD);
    } else {
        print STDERR "No need to run meryl for OVL ($asm.nmers.ovl.fasta exists).\n";
    }

    if (($obtT > 0) && (getGlobal("obtMerThreshold") ne $obtT) && (getGlobal("doOverlapBasedTrimming"))) {
        print STDERR "Reset OBT mer threshold from ", getGlobal("obtMerThreshold"), " to $obtT.\n";
        setGlobal("obtMerThreshold", $obtT);
    }
    
    if (($ovlT > 0) && (getGlobal("ovlMerThreshold") ne $ovlT)) {
        print STDERR "Reset OVL mer threshold from ", getGlobal("ovlMerThreshold"), " to $ovlT.\n";
        setGlobal("ovlMerThreshold", $ovlT);
    }
}

################################################################################
################################################################################
################################################################################



sub createOverlapJobs($) {
    my $isTrim = shift @_;

    return if (-d "$wrk/$asm.ovlStore");

    caFailure("overlapper detected no fragments", undef) if ($numFrags == 0);
    caFailure("overlapper needs to know if trimming or assembling", undef) if (!defined($isTrim));

    my $ovlThreads        = getGlobal("ovlThreads");
    my $ovlHashBits       = getGlobal("ovlHashBits");
    my $ovlHashLoad       = getGlobal("ovlHashLoad");

    my $outDir  = "1-overlapper";
    my $ovlOpt  = "";
    my $merSize = getGlobal("ovlMerSize");
    my $merComp = getGlobal("merCompression");
    my $overlap = "overlapInCore";

    if ($isTrim eq "trim") {
        $outDir  = "0-overlaptrim-overlap";
        $ovlOpt  = "-G";
        $merSize = getGlobal("obtMerSize");
        $overlap = "overlapInCore";
    }

    system("mkdir $wrk/$outDir") if (! -d "$wrk/$outDir");

    return if (-e "$wrk/$outDir/overlap.sh");

    #  umd overlapper here
    #
    if (getGlobal("ovlOverlapper") eq "umd") {
        #  For Sergey:
        #
        #  UMDoverlapper() needs to dump the gkpstore, run UMD, build
        #  the ovlStore and update gkpStore with new clear ranges.
        #  The explicit call to UMDoverlapper in main() can then go away.
        #  OBT is smart enough to disable itself if umd is enabled.
        #
        UMDoverlapper();
        return;
    }

    #  mer overlapper here
    #
    if ((($isTrim eq "trim") && (getGlobal("obtOverlapper") eq "mer")) ||
        (($isTrim ne "trim") && (getGlobal("ovlOverlapper") eq "mer"))) {
        merOverlapper($isTrim);
        return;
    }

    my $hashLibrary = 0;
    my $refLibrary = 0;
    if ($isTrim eq "trim") {
       $hashLibrary = getGlobal("obtHashLibrary");
       $refLibrary = getGlobal("obtRefLibrary");
    } else {
       $hashLibrary = getGlobal("ovlHashLibrary");
       $refLibrary = getGlobal("ovlRefLibrary");
    }

    #  To prevent infinite loops -- stop now if the overlap script
    #  exists.  This will unfortunately make restarting from transient
    #  failures non-trivial.
    #
    #  FAILUREHELPME
    #
    caFailure("overlapper failed\nmanual restart needed to prevent infinite loops\nremove file '$wrk/$outDir/overlap.sh'", undef) if (-e "$wrk/$outDir/overlap.sh");

    meryl();

    #  We make a giant job array for this -- we need to know hashBeg,
    #  hashEnd, refBeg and refEnd -- from that we compute batchName
    #  and jobName.
    #
    #  ovlopts.pl returns the batch name ($batchName), the job name
    #  ($jobName) and options to pass to overlap (-h $hashBeg-$hashEnd
    #  -r $refBeg-$refEnd).  From those, we can construct the command
    #  to run.
    #
    open(F, "> $wrk/$outDir/overlap.sh") or caFailure("can't open '$wrk/$outDir/overlap.sh'", undef);
    print F "#!" . getGlobal("shell") . "\n";
    print F "\n";
    print F "perl='/usr/bin/env perl'\n";
    print F "\n";
    print F "jobid=\$SGE_TASK_ID\n";
    print F "if [ x\$jobid = x -o x\$jobid = xundefined ]; then\n";
    print F "  jobid=\$1\n";
    print F "fi\n";
    print F "if [ x\$jobid = x ]; then\n";
    print F "  echo Error: I need SGE_TASK_ID set, or a job index on the command line.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "\n";
    print F "bat=`head -n \$jobid $wrk/$outDir/ovlbat | tail -n 1`\n";
    print F "job=`head -n \$jobid $wrk/$outDir/ovljob | tail -n 1`\n";
    print F "opt=`head -n \$jobid $wrk/$outDir/ovlopt | tail -n 1`\n";
    print F "jid=\$\$\n";
    print F "\n";
    print F "if [ ! -d $wrk/$outDir/\$bat ]; then\n";
    print F "  mkdir $wrk/$outDir/\$bat\n";
    print F "fi\n";
    print F "\n";
    print F "if [ -e $wrk/$outDir/\$bat/\$job.ovb.gz ]; then\n";
    print F "  echo Job previously completed successfully.\n";
    print F "  exit\n";
    print F "fi\n";
    print F "\n";
    print F "if [ x\$bat = x ]; then\n";
    print F "  echo Error: Job index out of range.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "\n";
    print F "AS_OVL_ERROR_RATE=" , getGlobal("ovlErrorRate"), "\n";
    print F "AS_CNS_ERROR_RATE=" , getGlobal("cnsErrorRate"), "\n";
    print F "AS_CGW_ERROR_RATE=" , getGlobal("cgwErrorRate"), "\n";
    print F "AS_OVERLAP_MIN_LEN=", getGlobal("ovlMinLen"),    "\n";
    print F "AS_READ_MIN_LEN="   , getGlobal("frgMinLen"),    "\n";
    print F "export AS_OVL_ERROR_RATE AS_CNS_ERROR_RATE AS_CGW_ERROR_RATE AS_OVERLAP_MIN_LEN AS_READ_MIN_LEN\n";

    print F getBinDirectoryShellCode();

    print F "\$bin/$overlap $ovlOpt --hashbits $ovlHashBits --hashload $ovlHashLoad -t $ovlThreads \\\n";
    print F "  \$opt \\\n";
    print F "  -k $merSize \\\n";
    print F "  -k $wrk/0-mercounts/$asm.nmers.obt.fasta \\\n" if ($isTrim eq "trim");
    print F "  -k $wrk/0-mercounts/$asm.nmers.ovl.fasta \\\n" if ($isTrim ne "trim");
    print F "  -o $wrk/$outDir/\$bat/\$job.ovb.WORKING.gz \\\n";
    print F " -H $hashLibrary -R $refLibrary \\\n";
    print F "  $wrk/$asm.gkpStore \\\n";
    print F "&& \\\n";
    print F "mv $wrk/$outDir/\$bat/\$job.ovb.WORKING.gz $wrk/$outDir/\$bat/\$job.ovb.gz\n";
    print F "\n";
    print F "exit 0\n";
    close(F);

    system("chmod +x $wrk/$outDir/overlap.sh");


    my $jobs      = 0;
    my $batchName = "";
    my $jobName   = "";

    {
        my $cmd;

        my $ovlHashBlockLength = getGlobal("ovlHashBlockLength");
        my $ovlHashBlockSize   = 0;
        my $ovlRefBlockSize    = getGlobal("ovlRefBlockSize");
        my $ovlRefBlockLength  = getGlobal("ovlRefBlockLength");

        if (($ovlRefBlockSize > 0) && ($ovlRefBlockLength > 0)) {
            caFailure("can't set both ovlRefBlockSize and ovlRefBlockLength", undef);
        }

        $cmd  = "$bin/overlap_partition \\\n";
        $cmd .= " -g  $wrk/$asm.gkpStore \\\n";
        $cmd .= " -bl $ovlHashBlockLength \\\n";
        $cmd .= " -bs $ovlHashBlockSize \\\n";
        $cmd .= " -rs $ovlRefBlockSize \\\n";
        $cmd .= " -rl $ovlRefBlockLength \\\n";
        $cmd .= " -o  $wrk/$outDir";

        if (runCommand($wrk, $cmd)) {
            caFailure("failed partition for overlapper", undef);
        }

        open(F, "< $wrk/$outDir/ovlbat") or caFailure("failed partition for overlapper: no ovlbat file found", undef);
        my @bat = <F>;
        close(F);

        open(F, "< $wrk/$outDir/ovljob") or caFailure("failed partition for overlapper: no ovljob file found", undef);
        my @job = <F>;
        close(F);

        $jobs      = scalar(@job);
        $batchName = $bat[$jobs-1];  chomp $batchName;
        $jobName   = $job[$jobs-1];  chomp $jobName;
    }

    print STDERR "Created $jobs overlap jobs.  Last batch '$batchName', last job '$jobName'.\n";

    #  Submit to the grid (or tell the user to do it), or just run
    #  things here
    #
    if (getGlobal("useGrid") && getGlobal("ovlOnGrid")) {
        my $sge        = getGlobal("sge");
        my $sgeName    = getGlobal("sgeName");
        my $sgeOverlap = getGlobal("sgeOverlap");

        $sgeName = "_$sgeName" if (defined($sgeName));

        my $SGE;
        $SGE  = "qsub $sge $sgeOverlap -cwd -N ovl_$asm$sgeName \\\n";
        $SGE .= "  -t 1-$jobs \\\n";
        $SGE .= "  -j y -o $wrk/$outDir/\\\$TASK_ID.out \\\n";
        $SGE .= "  $wrk/$outDir/overlap.sh\n";

	submitBatchJobs($SGE, "ovl_$asm$sgeName");
        exit(0);
    } else {
        for (my $i=1; $i<=$jobs; $i++) {
            my $out = substr("000000" . $i, -6);
            schedulerSubmit("$wrk/$outDir/overlap.sh $i > $wrk/$outDir/$out.out 2>&1");
        }

        schedulerSetNumberOfProcesses(getGlobal("ovlConcurrency"));
        schedulerFinish();
    }
}

################################################################################
################################################################################
################################################################################

#  Check that the overlapper jobs properly executed.  If not,
#  complain, but don't help the user fix things.


sub checkOverlapper ($) {
    my $isTrim = shift @_;

    my $outDir = "1-overlapper";
    my $ovlOpt = "";

    if ($isTrim eq "trim") {
        $outDir = "0-overlaptrim-overlap";
        $ovlOpt = "-G";
    }

    my $failedJobs = 0;
    my $failureMessage = "";

    open(B, "< $wrk/$outDir/ovlbat") or caFailure("failed to open '$wrk/$outDir/ovlbat'", undef);
    open(J, "< $wrk/$outDir/ovljob") or caFailure("failed to open '$wrk/$outDir/ovljob'", undef);

    while (!eof(B) && !eof(J)) {
        my $b = <B>;  chomp $b;
        my $j = <J>;  chomp $j;

        if ((! -e "$wrk/$outDir/$b/$j.ovb.gz") &&
            (! -e "$wrk/$outDir/$b/$j.ovb")) {
            $failureMessage .= "ERROR:  Overlap job $wrk/$outDir/$b/$j FAILED.\n";
            $failedJobs++;
        }
    }

    if (!eof(B) || !eof(J)) {
        print STDERR "Partitioning error; '$wrk/$outDir/ovlbat' and '$wrk/$outDir/ovljob' have extra lines.\n";
    }

    #  FAILUREHELPME
    #
    $failureMessage .= "\n$failedJobs overlapper jobs failed";
    caFailure($failureMessage, undef) if ($failedJobs);
}


sub checkMerOverlapper ($) {
    my $isTrim = shift @_;

    my $outDir = "1-overlapper";

    if ($isTrim eq "trim") {
        $outDir = "0-overlaptrim-overlap";
    }

    my $batchSize  = getGlobal("merOverlapperExtendBatchSize");
    my $jobs       = int($numFrags / $batchSize) + (($numFrags % $batchSize == 0) ? 0 : 1);
    my $failedJobs = 0;

    for (my $i=1; $i<=$jobs; $i++) {
        my $job = substr("0000" . $i, -4);

        if ((! -e "$wrk/$outDir/olaps/$job.ovb.gz") &&
            (! -e "$wrk/$outDir/olaps/$job.ovb")) {
            print STDERR "$wrk/$outDir/olaps/$job failed.\n";
            $failedJobs++;
        }
    }
    
    caFailure("$failedJobs overlapper jobs failed", undef) if ($failedJobs);
}


sub checkOverlap {
    my $isTrim = shift @_;

    caFailure("overlap checker needs to know if trimming or assembling", undef) if (!defined($isTrim));

    if ($isTrim eq "trim") {
        return if (-d "$wrk/$asm.obtStore");
        return if (-d "$wrk/$asm.dupStore");
        if      (getGlobal("obtOverlapper") eq "ovl") {
            checkOverlapper($isTrim);
        } elsif (getGlobal("obtOverlapper") eq "ovm") {
            checkOverlapper($isTrim);
        } elsif (getGlobal("obtOverlapper") eq "mer") {
            checkMerOverlapper($isTrim);
        } elsif (getGlobal("obtOverlapper") eq "umd") {
            caFailure("checkOverlap() wanted to check umd overlapper for obt?\n", undef);
        } else {
            caFailure("checkOverlap() unknown obt overlapper?\n", undef);
        }
    } else {
        return if (-d "$wrk/$asm.ovlStore");
        if      (getGlobal("ovlOverlapper") eq "ovl") {
            checkOverlapper($isTrim);
        } elsif (getGlobal("ovlOverlapper") eq "ovm") {
            checkOverlapper($isTrim);
        } elsif (getGlobal("ovlOverlapper") eq "mer") {
            checkMerOverlapper($isTrim);
        } elsif (getGlobal("ovlOverlapper") eq "umd") {
            #  Nop.
        } else {
            caFailure("checkOverlap() unknown ovl overlapper?\n", undef);
        }
    }
}

################################################################################
################################################################################
################################################################################

sub createOverlapStore {

    goto alldone if (-d "$wrk/$asm.ovlStore");

    if (runCommand($wrk, "find -L $wrk/1-overlapper \\( -name \\*ovb.gz -or -name \\*ovb \\) -print > $wrk/$asm.ovlStore.list")) {
        caFailure("failed to generate a list of all the overlap files", undef);
    }

    $cmd  = "$bin/overlapStoreBuild ";
    $cmd .= " -o $wrk/$asm.ovlStore.BUILDING ";
    $cmd .= " -g $wrk/$asm.gkpStore ";
    $cmd .= " -M " . getGlobal("ovlStoreMemory");
    $cmd .= " -L $wrk/$asm.ovlStore.list ";
    $cmd .= " > $wrk/$asm.ovlStore.err 2>&1";

    if (runCommand($wrk, $cmd)) {
        caFailure("failed to create the overlap store", "$wrk/$asm.ovlStore.err");
    }

    rename "$wrk/$asm.ovlStore.BUILDING", "$wrk/$asm.ovlStore";

    if (getGlobal("saveOverlaps") == 0) {
        open(F, "< $wrk/$asm.ovlStore.list");
        while (<F>) {
            chomp;
            unlink $_;
        }
        close(F);
    }

    rmrf("$wrk/$asm.ovlStore.list");
    rmrf("$wrk/$asm.ovlStore.err");

  alldone:
    stopAfter("overlapper");
}

################################################################################
################################################################################
################################################################################

sub overlapTrim {

    return if (getGlobal("doOverlapBasedTrimming") == 0);
    return if (getGlobal("ovlOverlapper") eq "umd");

    #  Skip overlap based trimming if it is done, or if the ovlStore already exists.
    #
    goto alldone if (-e "$wrk/0-overlaptrim/overlaptrim.success");
    goto alldone if (-d "$wrk/$asm.ovlStore");

    system("mkdir $wrk/0-overlaptrim")         if (! -d "$wrk/0-overlaptrim");
    system("mkdir $wrk/0-overlaptrim-overlap") if (! -d "$wrk/0-overlaptrim-overlap");


    #
    #  Compute overlaps, if we don't have them already -- the obtStore isn't needed for dedupe
    #

    if (! -e "$wrk/0-overlaptrim/$asm.dupStore") {
        createOverlapJobs("trim");
        checkOverlap("trim");
    }

    #
    #  Deduplicate?
    #

    if (! -e "$wrk/0-overlaptrim/$asm.deduplicate.summary") {

        if (! -e "$wrk/0-overlaptrim/$asm.dupStore") {
            if (runCommand("$wrk/0-overlaptrim",
                           "find -L $wrk/0-overlaptrim-overlap -follow \\( -name \\*ovb.gz -or -name \\*ovb \\) -print > $wrk/0-overlaptrim/$asm.dupStore.list")) {
                caFailure("failed to generate a list of all the overlap files", undef);
            }

            $cmd  = "$bin/overlapStoreBuild \\\n";
            $cmd .= " -dup \\\n";
            $cmd .= " -o $wrk/0-overlaptrim/$asm.dupStore.BUILDING \\\n";
            $cmd .= " -g $wrk/$asm.gkpStore \\\n";
            $cmd .= " -M \\\n" . getGlobal('ovlStoreMemory');
            $cmd .= " -L $wrk/0-overlaptrim/$asm.dupStore.list \\\n";
            $cmd .= " > $wrk/0-overlaptrim/$asm.dupStore.err 2>&1";

            if (runCommand("$wrk/0-overlaptrim", $cmd)) {
                caFailure("failed to build the dup store", "$wrk/0-overlaptrim/$asm.dupStore.err");
            }

            rename "$wrk/0-overlaptrim/$asm.dupStore.BUILDING", "$wrk/0-overlaptrim/$asm.dupStore";

            #  Delete overlaps unless we're told to save them
            if (getGlobal("saveOverlaps") == 0) {
                open(F, "< $wrk/0-overlaptrim/$asm.dupStore.list");
                while (<F>) {
                    chomp;
                    unlink $_;
                }
                close(F);
            }

            rmrf("$asm.dupStore.list");
            rmrf("$asm.dupStore.err");
        }

        $cmd  = "$bin/deduplicate \\\n";
        $cmd .= "-gkp     $wrk/$asm.gkpStore \\\n";
        #$cmd .= "-ovs     $wrk/0-overlaptrim/$asm.obtStore \\\n";  #  Not needed
        $cmd .= "-ovs     $wrk/0-overlaptrim/$asm.dupStore \\\n";
        $cmd .= "-report  $wrk/0-overlaptrim/$asm.deduplicate.log \\\n";
        $cmd .= "-summary $wrk/0-overlaptrim/$asm.deduplicate.summary \\\n";
        $cmd .= "> $wrk/0-overlaptrim/$asm.deduplicate.err 2>&1";

        stopBefore("deDuplication", $cmd);

        if (runCommand("$wrk/0-overlaptrim", $cmd)) {
            unlink "$wrk/0-overlaptrim/$asm.deduplicate.summary";
            caFailure("failed to deduplicate the reads", "$wrk/0-overlaptrim/$asm.deduplicate.err");
        }
    }

    touch("$wrk/0-overlaptrim/overlaptrim.success");

  alldone:
    stopAfter("overlapBasedTrimming");
    stopAfter("OBT");
}


################################################################################
################################################################################
################################################################################


setDefaults();

#  Check for the presence of a -options switch BEFORE we do any work.
#  This lets us print the default values of options.

foreach my $arg (@ARGV) {
    if ($arg eq "-options") {
        setGlobal("options", 1);
        printHelp();
    }
}

#  At some pain, we stash the original options for later use.  We need
#  to use these when we resubmit ourself to SGE.  We can't simply dump
#  all of @ARGV into here, because we need to fix up relative paths.

while (scalar(@ARGV)) {
    my $arg = shift @ARGV;

    if      ($arg =~ m/^-d/) {
        $wrk = shift @ARGV;
        $wrk = "$ENV{'PWD'}/$wrk" if ($wrk !~ m!^/!);
        $commandLineOptions .= " -d \"$wrk\"";

    } elsif ($arg eq "-p") {
        $asm = shift @ARGV;
        $commandLineOptions .= " -p \"$asm\"";

    } elsif ($arg eq "-s") {
        push @specFiles, shift @ARGV;

    } elsif ($arg eq "-version") {
        setGlobal("version", 1);

    } elsif ($arg eq "-options") {
        #  Do nothing.  Handled above, but we still need to process it here.
        #setGlobal("options", 1);

    } elsif (($arg =~ /\.frg$|frg\.gz$|frg\.bz2$/i) && (-e $arg)) {
        $arg = "$ENV{'PWD'}/$arg" if ($arg !~ m!^/!);
        push @fragFiles, $arg;
        $commandLineOptions .= " \"$arg\"";

    } elsif (($arg =~ /\.sff$|sff\.gz$|sff\.bz2$/i) && (-e $arg)) {
        $arg = "$ENV{'PWD'}/$arg" if ($arg !~ m!^/!);
        push @fragFiles, $arg;
        $commandLineOptions .= " \"$arg\"";

    } elsif (($arg =~ /\.ace$/i) && (-e $arg)) {
        $arg = "$ENV{'PWD'}/$arg" if ($arg !~ m!^/!);
        push @fragFiles, $arg;
        $commandLineOptions .= " \"$arg\"";

    } elsif ($arg =~ m/=/) {
        push @specOpts, $arg;
        $commandLineOptions .= " \"$arg\"";

    } else {
        setGlobal("help",
                  getGlobal("help") . "File not found or invalid command line option '$arg'\n");
    }
}


setGlobal("help", getGlobal("help") . "Assembly name prefix not supplied with -p.\n") if (!defined($asm));
setGlobal("help", getGlobal("help") . "Directory not supplied with -d.\n")            if (!defined($wrk));


$bin = getBinDirectory();

@fragFiles = setParametersFromFile("$bin/spec/runCA.default.specFile", @fragFiles)   if (-e "$bin/spec/runCA.default.specFile");
@fragFiles = setParametersFromFile("$ENV{'HOME'}/.runCA",              @fragFiles)   if (-e "$ENV{'HOME'}/.runCA");


#  For each of the specfiles on the command line, find the actual file and make it an absolute path.
#  These can be in the current directory (e.g., 'my.spec'), or in the installed directory ('$bin/spec').
#
foreach my $specFile (@specFiles) {

    if ((-e "$specFile") && (! -d "$specFile")) {
        $specFile = "$ENV{'PWD'}/$specFile" if ($specFile !~ m!^/!);

    } elsif ((-e "$bin/spec/$specFile") && (! -d "$bin/spec/$specFile")) {
        $specFile = "$bin/spec/$specFile";

    } elsif ((-e "$bin/spec/$specFile.specFile") && (! -d "$bin/spec/$specFile.specFile")) {
        $specFile = "$bin/spec/$specFile.specFile";

    } else {
        die "specFile '$specFile' not found.\n";
    }

    $commandLineOptions .= " -s \"$specFile\"";

    @fragFiles = setParametersFromFile($specFile, @fragFiles);
}

setParametersFromCommandLine(@specOpts);

setParameters();

printHelp();

#  Fail immediately if we run the script on the grid, and the gkpStore
#  directory doesn't exist and we have no input files.  Without this
#  check we'd fail only after being scheduled on the grid.
#
if ((getGlobal("scriptOnGrid") == 1) &&
    (! -d "$wrk/$asm.gkpStore") &&
    (scalar(@fragFiles) == 0)) {
    caFailure("no fragment files specified, and stores not already created", undef);
}

checkDirectories();
outputSpecLog();

#  If not already on the grid, see if we should be on the grid.
#  N.B. the arg MUST BE undef.
#
submitScript(undef) if (!runningOnGrid());

#  Begin

if (defined(getGlobal("outputPrefix"))) {
    die "'outputPrefix' obsolete; output now written to work directory (-d) using assembly name (-p).\n";
}

preoverlap(@fragFiles);
overlapTrim();

#  End

if (! -e "$wrk/$asm.gkpStore.dump.err") {
    $cmd = "$bin/gatekeeper -dumpfastq $wrk/$asm $wrk/$asm.gkpStore > $wrk/$asm.gkpStore.dump.err 2>&1";

    if (runCommand($wrk, $cmd)) {
        caFailure("failed to dump deduplicated reads", "$wrk/$asm.gkpStore.dump.err");
    }
}


if (! -e "$wrk/$asm.gkpStore.dumpRename.err") {
    $cmd = "$bin/replaceUIDwithName $wrk/$asm.gkpStore.fastqUIDmap $wrk/$asm.1.fastq $wrk/$asm.2.fastq $wrk/$asm.paired.fastq $wrk/$asm.unmated.fastq > $wrk/$asm.gkpStore.dumpRename.err 2>&1";

    if (runCommand($wrk, $cmd)) {
        caFailure("failed to restore names in deduplicated reads", "$wrk/$asm.$asm.gkpStore.dumpRename.err");
    }

    rename "$wrk/0-overlaptrim/$asm.deduplicate.summary", "$wrk/$asm.deduplicate.summary";
    rename "$wrk/0-overlaptrim/$asm.deduplicate.log",     "$wrk/$asm.deduplicate.log";
}

rmrf("$wrk/0-mercounts");
rmrf("$wrk/0-overlaptrim");
rmrf("$wrk/0-overlaptrim-overlap");

rmrf("$wrk/$asm.gkpStore");
rmrf("$wrk/$asm.gkpStore.err");
rmrf("$wrk/runCA-logs");

exit(0);
