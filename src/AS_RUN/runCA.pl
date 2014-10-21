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


sub findBlasr($) {
   my $cns = shift @_;

   if ($cns ne "pbdagcon" && $cns ne "pbutgcns") {
      return undef;
   }

   my $CA = getBinDirectory();
   my $BLASR = "$CA/../../../smrtanalysis/current/analysis/bin/";

   if (! -e "$BLASR/blasr" && $cns eq "pbdagcon") {
      if (-e "$CA/blasr") {
         $BLASR = $CA;
      } else {
         # try to use path
         my $amosPath = `which blasr`;
         chomp $amosPath;
         my @t = split '/', "$amosPath";
         pop @t;                      #  blasr
         $BLASR = join '/', @t;  #  path to the assembler
      }
      # if we really can't find it just give up
      if (! -e "$BLASR/blasr") {
         return undef;
      }

      # check for consensus too
      # make sure we have the pb consensus module available if it was requested
      if (! -e "$BLASR/blasr" || ! -e "$BLASR/pbdagcon") {
         return undef;
      }
   } elsif (! -e "$BLASR/pbutgcns" && $cns eq "pbutgcns") {
      if (-e "$CA/pbutgcns") {
         $BLASR = $CA;
      } else {
         # try to use path
         my $amosPath = `which pbutgcns`;
         chomp $amosPath;
         my @t = split '/', "$amosPath";
         pop @t;                      #  pbutgcns
         $BLASR = join '/', @t;  #  path to the assembler
      }
      # if we really can't find it just give up
      if (! -e "$BLASR/pbutgcns") {
         return undef;
      }
   }
   return $BLASR;
}

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
sub getInstallDirectory () {
    my @t = split '/', "$FindBin::RealBin";
    pop @t;                         #  bin
    pop @t;                         #  arch, e.g., FreeBSD-amd64
    my $installDir = join '/', @t;  #  path to the assembler

    return($installDir);
}

sub getBinDirectory () {
    my $installDir = getInstallDirectory();

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
    my $installDir = getInstallDirectory();
    my $string;

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

    if ($var eq "merDistinct") {
        setGlobal("obtMerDistinct", $val);
        setGlobal("ovlMerDistinct", $val);
        return;
    }

    if ($var eq "merTotal") {
        setGlobal("obtMerTotal", $val);
        setGlobal("ovlMerTotal", $val);
        return;
    }

    if ($var eq "overlapper") {
        setGlobal("obtOverlapper", $val);
        setGlobal("ovlOverlapper", $val);
        return;
    }

    if (($var eq "gridEngine") && ($val eq "SGE")) {
        setGlobal("gridSubmitCommand",      "qsub");
        setGlobal("gridHoldOption",         "-hold_jid \"WAIT_TAG\"");
        setGlobal("gridHoldOptionNoArray",  undef);
        setGlobal("gridSyncOption",         "-sync y");
        setGlobal("gridNameOption",         "-cwd -N");
        setGlobal("gridArrayOption",        "-t ARRAY_JOBS");
        setGlobal("gridArrayName",          "ARRAY_NAME");
        setGlobal("gridOutputOption",       "-j y -o");
        setGlobal("gridPropagateCommand",   "qalter -hold_jid \"WAIT_TAG\"");
        setGlobal("gridNameToJobIDCommand", undef);
        setGlobal("gridTaskID",             "SGE_TASK_ID");
        setGlobal("gridArraySubmitID",      "\\\$TASK_ID");
        setGlobal("gridJobID",              "JOB_ID");
    }

    if (($var eq "gridEngine") && ($val eq "LSF")) {
        setGlobal("gridSubmitCommand",      "bsub");
        setGlobal("gridHoldOption",         "-w \"numended\(\"WAIT_TAG\", \*\)\"");
        setGlobal("gridHoldOptionNoArray",  "-w \"done\(\"WAIT_TAG\"\)\"");
        setGlobal("gridSyncOption",         "-K");
        setGlobal("gridNameOption",         "-J");
        setGlobal("gridArrayOption",        "");
        setGlobal("gridArrayName",          "ARRAY_NAME\[ARRAY_JOBS\]");
        setGlobal("gridOutputOption",       "-o");
        setGlobal("gridPropagateCommand",   "bmodify -w \"done\(\"WAIT_TAG\"\)\"");
        setGlobal("gridNameToJobIDCommand", "bjobs -A -J \"WAIT_TAG\" | grep -v JOBID");
        setGlobal("gridNameToJobIDCommandNoArray", "bjobs -J \"WAIT_TAG\" | grep -v JOBID");
        setGlobal("gridTaskID",             "LSB_JOBINDEX");
        setGlobal("gridArraySubmitID",      "%I");
        setGlobal("gridJobID",              "LSB_JOBID");
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

    $global{"overlapStoreOnGrid"}          = 0;
    $synops{"overlapStoreOnGrid"}          = "Enable OverlapStore Build on Grid";


    #####  General Configuration Options (aka miscellany)

    $global{"showNext"}                    = undef;
    $synops{"showNext"}                    = "Don't run any commands, just report what would run";

    $global{"pathMap"}                     = undef;
    $synops{"pathMap"}                     = "File with a hostname to binary directory map";

    $global{"shell"}                       = "/bin/sh";
    $synops{"shell"}                       = "Command interpreter to use; sh-compatible (e.g., bash), NOT C-shell (csh or tcsh)";

    #####  Error Rates

    $global{"ovlErrorRate"}                = 0.06;
    $synops{"ovlErrorRate"}                = "Overlaps above this error rate are not computed";

    $global{"obtErrorRate"}                = undef;
    $synops{"obtErrorRate"}                = "Overlaps at or below this error rate are used for Overlap Based Trimming (OBT)";

    $global{"obtErrorLimit"}               = undef;
    $synops{"obtErrorLimit"}               = "Overlaps at or below this number of errors are used for Overlap Based Trimming (OBT)";

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

    #####  Grid Engine configuration, how to submit jobs, etc

    $global{"gridEngine"}		   = "SGE";
    $global{"gridSubmitCommand"}           = "qsub";
    $global{"gridHoldOption"}              = "-hold_jid \"WAIT_TAG\"";         # for lsf: -w "done("WAIT_TAG")"
    $global{"gridHoldOptionNoArray"}       = undef;
    $global{"gridSyncOption"}              = "-sync y";                        # for lsf: -K
    $global{"gridNameOption"}              = "-cwd -N";                        # for lsf: -J
    $global{"gridArrayOption"}             = "-t ARRAY_JOBS";                  # for lsf: empty ("")
    $global{"gridArrayName"}               = "ARRAY_NAME";                     # for lsf: ARRAY_NAME[ARRAY_JOBS]
    $global{"gridOutputOption"}            = "-j y -o";                        # for lsf: -o
    $global{"gridPropagateCommand"}        = "qalter -hold_jid \"WAIT_TAG\"";  # for lsf: bmodify -w "done(WAIT_TAG)"
    $global{"gridNameToJobIDCommand"}      = undef;                            # for lsf: bjobs -J "WAIT_TAG" | grep -v JOBID
    $global{"gridNameToJobIDCommandNoArray"} = undef;
    $global{"gridTaskID"}                  = "SGE_TASK_ID";                    # for lsf: LSB_JOBINDEX
    $global{"gridArraySubmitID"}           = "\\\$TASK_ID";                    # for lsf: %I
    $global{"gridJobID"}                   = "JOB_ID";                         # for lsf: LSB_JOBID

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

    $global{"mbtIlluminaAdapter"}          = 0;
    $synops{"mbtIlluminaAdapter"}          = "Remove Illumina adapter sequence during merTrim";

    $global{"mbt454Adapter"}               = 0;
    $synops{"mbt454Adapter"}               = "Remove 454 adapter sequence during merTrim";

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

    $global{"ovlMerThreshold"}             = undef;
    $synops{"ovlMerThreshold"}             = "K-mer frequency threshold; mers more frequent than this count are ignored";

    $global{"ovlMerDistinct"}              = undef;
    $synops{"ovlMerDistinct"}              = "K-mer frequency threshold; the least frequent fraction of distinct mers can seed overlaps";

    $global{"ovlMerTotal"}                 = undef;
    $synops{"ovlMerTotal"}                 = "K-mer frequency threshold; the least frequent fraction of all mers can seed overlaps";

    $global{"ovlFrequentMers"}             = undef;
    $synops{"ovlFrequentMers"}             = "Do not seed overlaps with these kmers (fasta format)";

    $global{"obtMerSize"}                  = 22;
    $synops{"obtMerSize"}                  = "K-mer size";

    $global{"obtMerThreshold"}             = undef;
    $synops{"obtMerThreshold"}             = "K-mer frequency threshold; mers more frequent than this are ignored";

    $global{"obtMerDistinct"}              = undef;
    $synops{"obtMerDistinct"}              = "K-mer frequency threshold; the least frequent fraction of distinct mers can seed overlaps";

    $global{"obtMerTotal"}                 = undef;
    $synops{"obtMerTotal"}                 = "K-mer frequency threshold; the least frequent fraction of all mers can seed overlaps";

    $global{"obtFrequentMers"}              = undef;
    $synops{"obtFrequentMers"}              = "Do not seed overlaps with these kmers (fasta format)";

    $global{"ovlHashLibrary"}               = "0";
    $synops{"ovlHashLibrary"}               = "For ovl overlaps, only load hash fragments from specified lib, 0 means all";

    $global{"ovlRefLibrary"}                = "0";
    $synops{"ovlRefLibrary"}                = "For ovl overlaps, only load ref fragments from specified lib, 0 means all";

    $global{"obtHashLibrary"}               = "0";
    $synops{"obtHashLibrary"}               = "For obt overlaps, only load hash fragments from specified lib, 0 means all";

    $global{"obtRefLibrary"}                = "0";
    $synops{"obtRefLibrary"}                = "For obt overlaps, only load ref fragments from specified lib, 0 means all";

    $global{"obtCheckLibrary"}              = 1;
    $synops{"obtCheckLibrary"}              = "Check that all libraries are used during obt overlaps";

    $global{"ovlCheckLibrary"}              = 1;
    $synops{"ovlCheckLibrary"}              = "Check that all libraries are used during ovl overlaps";


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

    $global{"batOptions"}                  = undef;
    $synops{"batOptions"}                  = "Advanced options to bogart";

    $global{"batMemory"}                   = undef;
    $synops{"batMemory"}                   = "Approximate maximum memory usage for loading overlaps, in gigabytes, default is unlimited";

    $global{"batThreads"}                  = undef;
    $synops{"batThreads"}                  = "Number of threads to use in the Merge/Split/Join phase; default is whatever OpenMP wants";

    $global{"doUnitigSplitting"}           = 1;
    $synops{"doUnitigSplitting"}           = "Split unitigs based on low coverage and high bad mate evidence";

    #####  Unitig Repeat/Unique Options (formerly in scaffolder)

    $global{"astatLowBound"}               = 1;
    $synops{"astatLowBound"}               = "EXPERT!";

    $global{"astatHighBound"}              = 5;
    $synops{"astatHighBound"}              = "EXPERT!";

    $global{"maxSingleReadSpan"}           = undef;
    $synops{"maxSingleReadSpan"}           = "Unitigs with a single read spanning more than this fraction of the unitig are never labeled unique";

    $global{"lowCoverageDepth"}            = undef;
    $synops{"lowCoverageDepth"}            = "See lowCoverageAllowed";

    $global{"lowCoverageAllowed"}          = undef;
    $synops{"lowCoverageAllowed"}          = "Unitigs with more than this fraction lowCoverageDepth bases are never labeled unique";

    $global{"minReadsUnique"}              = undef;
    $synops{"minReadsUnique"}              = "Unitigs with fewer reads that this are never labeled unique";

    $global{"maxRepeatLength"}             = undef;
    $synops{"maxRepeatLength"}             = "Unitigs longer than this are always labeled unique";

    #####  Scaffolder Options

    $global{"cgwPurgeCheckpoints"}         = 1;
    $synops{"cgwPurgeCheckpoints"}         = "Remove cgw checkpoint files when a scaffolding step finishes successfully";

    $global{"cgwCompressTigStore"}         = 0;
    $synops{"cgwCompressTigStore"}         = "Remove tigStore versions when a scaffolding step finishes successfully";

    $global{"cgwDemoteRBP"}                = 1;
    $synops{"cgwDemoteRBP"}                = "EXPERT!";

    $global{"cgwUseUnitigOverlaps"}        = 0;
    $synops{"cgwUseUnitigOverlaps"}        = "Use unused best overlaps (from BOG) in scaffolder (EXPERIMENTAL)";

    $global{"cgwReloadMates"}              = 0;
    $synops{"cgwReloadMates"}              = "Load new mate pairs from gkpStore after ckp is loaded (EXPERIMENTAL)";

    $global{"stoneLevel"}                  = 2;
    $synops{"stoneLevel"}                  = "EXPERT!";

    $global{"cgwMergeFilterLevel"}         = 1;
    $synops{"cgwMergeFilterLevel"}         = "Apply no (0), classic (1, default), suggested (2), or stringent (5) criteria before considering a scaffold merge";

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

    $global{"cgwMinMergeWeight"}           = 2;
    $synops{"cgwMinMergeWeight"}           = "When merging scaffolds, do not use edges with weight below this.\n";

    $global{"cgwPreserveConsensus"}        = 0;
    $synops{"cgwPreserveConsensus"}        = "Do not remove the contig consensus sequence at the end of scaffolder; ctgcns will be skipped (faster) but quality will be lower.\n";

    #####  Consensus Options

    $global{"cnsPartitions"}               = 128;
    $synops{"cnsPartitions"}               = "Partition consensus into N jobs";

    $global{"cnsMinFrags"}                 = 75000;
    $synops{"cnsMinFrags"}                 = "Don't make a consensus partition with fewer than N fragments";

    $global{"cnsConcurrency"}              = 2;
    $synops{"cnsConcurrency"}              = "If not SGE, number of consensus jobs to run at the same time";

    $global{"cnsPhasing"}                  = 0;
    $synops{"cnsPhasing"}                  = "Options for consensus phasing of SNPs\n\t0 - Do not phase SNPs to be consistent.\n\t1 - If two SNPs are joined by reads, phase them to be consistent.";

    $global{"cnsMaxCoverage"}              = 0;
    $synops{"cnsMaxCoverage"}              = "Limit unitig consensus to to at most this coverage";

    $global{"cnsReuseUnitigs"}             = 0;
    $synops{"cnsReuseUnitigs"}             = "Do not compute single-unitig contigs again, just reuse the unitig.";

    #$global{"cnsRecycleUnitigs"}          = 0;
    #$synops{"cnsRecycleUnitigs"}          = "At some point, we'll come up with something for this.";

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

    $global{"dumpFASTQ"}                   = undef;
    $synops{"dumpFASTQ"}                   = "At the end of a successful assembly, output final reads in FASTQ format";

    #####  Options for toggling assembly.

    $global{"doToggle"}                     = 0;
    $synops{"doToggle"}                     = "At the end of a successful assembly, search for placed surrogates and toggle them to be unique unitigs. Re-run the assembly starting from scaffolder";

    $global{"toggleUnitigLength"}           = 2000;
    $synops{"toggleUnitigLength"}           = "Minimum length for a surrogate to be toggled.";

    $global{"toggleNumInstances"}           = 1;
    $synops{"toggleNumInstances"}           = "Number of instances for a surrogate to be toggled. If 0 is specified, all non-singleton unitigs are toggled to unique status.";

    $global{"toggleMaxDistance"}            = 1000;
    $synops{"toggleMaxDistance"}            = "Toggling will look for surrogates that appear exactly twice, both at the end of a scaffold. This parameter specifies how close to the scaffold end the surrogate must be.";

    $global{"toggleDoNotDemote"}            = 0;
    $synops{"toggleDoNotDemote"}            = "Do not allow CGW to demote toggled unitigs based on branching patterns.";

    #### Closure Options

    $global{"closureOverlaps"}              = undef;
    $synops{"closureOverlaps"}              = "Option for handling overlaps involving closure reads.\n\t0 - Treat them just like regular reads, \n\t1 - Do not allow any overlaps (i.e. closure reads will stay as singletons until scaffolding), \n\t2 - allow overlaps betweeen closure reads and non-closure reads only";

    $global{"closurePlacement"}             = 2;
    $synops{"closurePlacement"}             = "Option for placing closure reads using the constraints.\n\t0 - Place at the first location found\n\t1 - Place at the best location (indicated by most constraints)\n\t2 - Place at multiple locations as long as the closure read/unitig in question is not unique";

    #####  Ugly, command line options passed to printHelp()

    $global{"help"}                        = "";
    $synops{"help"}                        = undef;

    $global{"version"}                     = 0;
    $synops{"version"}                     = undef;

    $global{"options"}                     = 0;
    $synops{"options"}                     = undef;


    #  If this is set, it breaks the consensus.sh and overlap.sh scripts.  Good grief!  Why
    #  are you running runCA in a task array!?
    #
    if (exists($ENV{$global{"gridTaskID"}})) {
        undef $ENV{$global{"gridTaskID"}};
        print STDERR "ENV: $global{'gridTaskID'} needs to be unset, done.\n";
    }


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
        $val =~ s/\\\"/\"/g;
        $val =~ s/\"/\\\"/g;
        $val =~ s/\\\$/\$/g;
        $val =~ s/\$/\\\$/g;
        $commandLineOptions .= " \"$var=$val\" ";
    }
}

sub fixCase ($) {
    my $var = shift @_;
    my $val = getGlobal($var);

    if (defined($val)) {
        $val =~ tr/A-Z/a-z/;
        setGlobal($var, $val);
    }
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

    makeAbsolute("obtFrequentMers");
    makeAbsolute("ovlFrequentMers");

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
    if ((getGlobal("consensus") ne "cns") && (getGlobal("consensus") ne "seqan") && (getGlobal("consensus") ne "pbdagcon") && (getGlobal("consensus") ne "pbutgcns")) {
        caFailure("invalid consensus specified (" . getGlobal("consensus") . "); must be 'cns' or 'seqan' or 'pbdagcon' or 'pbutgcns'", undef);
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
                         "meryl",
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

    if (($ovlER < 0.0) || (0.40 < $ovlER)) {
        caFailure("ovlErrorRate is $ovlER, this MUST be between 0.00 and 0.40", undef);
    }
    if (($utgER < 0.0) || (0.40 < $utgER)) {
        caFailure("utgErrorRate is $utgER, this MUST be between 0.00 and 0.40", undef);
    }
    if (($cgwER < 0.0) || (0.40 < $cgwER)) {
        caFailure("cgwErrorRate is $cgwER, this MUST be between 0.00 and 0.40", undef);
    }
    if (($cnsER < 0.0) || (0.40 < $cnsER)) {
        caFailure("cnsErrorRate is $cnsER, this MUST be between 0.00 and 0.40", undef);
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
        print "usage: runCA -d <dir> -p <prefix> [options] <frg> ...\n";
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

    open(F, "find -L $dir -type f -name $asm.ckp.\[0-9\]\* |");
    while (<F>) {
        if (m/ckp.(\d+)\s*$/) {
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

    open(F, "find -L $dir -type f -name $asm.ckp.\[0-9\]\* |");
    while (<F>) {
        if (m/ckp.(\d+)\s*$/) {
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




#  Default to unitigger, unless the gkpStore says otherwise.
#
sub getUnitigger () {
    my $unitigger = getGlobal("unitigger");

    if (!defined($unitigger)) {
        setGlobal("unitigger", "utg");

        if (system("$bin/gatekeeper -isfeatureset 0 forceBOGunitigger $wrk/$asm.gkpStore") == 0) {
            setGlobal("unitigger", "bog");
        }
    }

    $unitigger = getGlobal("unitigger");

    die "Unitigger not defined.\n" if (!defined($unitigger));

    return($unitigger);
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
    my $jobID = getGlobal("gridJobID");
    return(exists($ENV{$jobID}));
}

sub findNextScriptOutputFile () {
    my $idx = "00";
    while (-e "$wrk/runCA.sge.out.$idx.sh") {
        $idx++;
    }
    return("$wrk/runCA.sge.out.$idx");
}

sub buildGridArray($$$) {
    my $name = shift @_;
    my $maxLimit = shift @_;
    my $globalValue = shift @_;

    my $arrayJobName = getGlobal($globalValue);
    $arrayJobName =~ s/ARRAY_NAME/$name/g;
    $arrayJobName =~ s/ARRAY_JOBS/1-$maxLimit/g;

    return $arrayJobName;
}

sub getGridArrayName($$) {
    my $name = shift @_;
    my $maxLimit = shift @_;
    return buildGridArray($name, $maxLimit, "gridArrayName");
}

sub getGridArrayOption($$) {
    my $name = shift @_;
    my $maxLimit = shift @_;
    return buildGridArray($name, $maxLimit, "gridArrayOption");
}

sub submitScript ($) {
    my $waitTag = shift @_;

    return if (getGlobal("scriptOnGrid") == 0);

    my $output = findNextScriptOutputFile();
    my $script = "$output.sh";

    open(F, "> $script") or caFailure("failed to open '$script' for writing", undef);
    print F "#!" . getGlobal("shell") . "\n";
    print F "#\n";
    print F "if [ \"x\$SGE_ROOT\" != \"x\" ]; then \n";
    print F "   #  Attempt to (re)configure SGE.  For reasons Bri doesn't know,\n";
    print F "   #  jobs submitted to SGE, and running under SGE, fail to read his\n";
    print F "   #  .tcshrc (or .bashrc, limited testing), and so they don't setup\n";
    print F "   #  SGE (or ANY other paths, etc) properly.  For the record,\n";
    print F "   #  interactive SGE logins (qlogin, etc) DO set the environment.\n";
    print F "   \n";
    print F "   . \$SGE_ROOT/\$SGE_CELL/common/settings.sh\n";
    print F "fi\n";
    print F "\n";
    print F "#  On the off chance that there is a pathMap, and the host we\n";
    print F "#  eventually get scheduled on doesn't see other hosts, we decide\n";
    print F "#  at run time where the binary is.\n";

    print F getBinDirectoryShellCode();

    print F "/usr/bin/env perl \$bin/runCA $commandLineOptions\n";
    close(F);

    system("chmod +x $script");

    my $sge         = getGlobal("sge");
    my $sgeName     = getGlobal("sgeName");
    my $sgeScript   = getGlobal("sgeScript");
    my $sgePropHold = getGlobal("sgePropagateHold");

    my $submitCommand  = getGlobal("gridSubmitCommand");
    my $holdOption   = getGlobal("gridHoldOption");
    my $nameOption   = getGlobal("gridNameOption");
    my $outputOption  = getGlobal("gridOutputOption");
    my $holdPropagateCommand  = getGlobal("gridPropagateCommand");

    $sgeName = "_$sgeName"              if (defined($sgeName));
    my $jobName = "rCA_$asm$sgeName";

    if (defined($waitTag)) {
        my $hold = $holdOption;
        if (getGlobal("gridEngine") eq "LSF"){
           my $tcmd = getGlobal("gridNameToJobIDCommand");
           $tcmd =~ s/WAIT_TAG/$waitTag/g;
           my $propJobCount = `$tcmd |wc -l`;
           chomp $propJobCount;
           if ($propJobCount == 0) {
              $tcmd = getGlobal("gridNameToJobIDCommandNoArray");
              $tcmd =~ s/WAIT_TAG/$waitTag/g;
              $hold = getGlobal("gridHoldOptionNoArray");
              $propJobCount = `$tcmd |wc -l`;
           }
           if ($propJobCount != 1) {
              print STDERR "Warning: multiple IDs for job $sgePropHold got $propJobCount and should have been 1.\n";
           }
           my $jobID = `$tcmd |tail -n 1 |awk '{print \$1}'`;
           chomp $jobID;
           $hold =~ s/WAIT_TAG/$jobID/g;
        } else{
           $hold =~ s/WAIT_TAG/$waitTag/g;
        }
        $waitTag = $hold;
    }
    my $qcmd = "$submitCommand $sge $sgeScript $nameOption \"$jobName\" $waitTag $outputOption $output  $script";
    runCommand($wrk, $qcmd) and caFailure("Failed to submit script.\n");

    if (defined($sgePropHold)) {
        if (defined($holdPropagateCommand)) {
            my $translateCmd = getGlobal("gridNameToJobIDCommandNoArray");

            # translate hold option to job id if necessary
            if (defined($translateCmd) && $translateCmd ne "") {
                my $tcmd = $translateCmd;
                $tcmd =~ s/WAIT_TAG/$sgePropHold/g;
                my $propJobCount = `$tcmd |wc -l`;
                chomp $propJobCount;
                if ($propJobCount != 1) {
                    print STDERR "Warning: multiple IDs for job $sgePropHold got $propJobCount and should have been 1.\n";
                }
                #my $jobID = `$tcmd |head -n 1 |awk '{print \$1}'`;
                #chomp $jobID;
                #print STDERR "Translated job ID $sgePropHold to be job $jobID\n";
                #$sgePropHold = $jobID;
                open(PROPS, "$tcmd |awk '{print \$1}' | ") or die("Couldn't get list of jobs that need to hold", undef);

                # now we can get the job we are holding for
                $tcmd = $translateCmd;
                $tcmd =~ s/WAIT_TAG/$jobName/g;
                my $holdJobCount = `$tcmd |wc -l`;
                chomp $propJobCount;
                if ($propJobCount != 1) {
                    print STDERR "Warning: multiple IDs for job $jobName got $propJobCount and should have been 1.\n";
                }
                #$jobID = `$tcmd |head -n 1 |awk '{print \$1}'`;
                #chomp $jobID;
                #print STDERR "Translated job ID $sgePropHold to be job $jobID\n";
                #$jobName = $jobID;
                open(HOLDS, "$tcmd |awk '{print \$1}' | ") or die("Couldn't get list of jobs that should be held for", undef);

                # loop over all jobs and all sge hold commands to modify the jobs. We have no way to know which is the right one unfortunately
                while (my $prop = <PROPS>) {
                    while (my $hold = <HOLDS>) {
                        chomp $hold;
                        chomp $prop;
                        my $hcmd = $holdPropagateCommand;
                        $hcmd =~ s/WAIT_TAG/$hold/g;
                        my $acmd = "$hcmd $prop";
                        print STDERR "Propagating hold to $prop to wait for job $hold\n";
                        system($acmd) and print STDERR "WARNING: Failed to reset hold_jid trigger on '$prop'.\n";
                    }
                }
                close(HOLDS);
                close(PROPS);
            } else {
                $sgePropHold = "\"$sgePropHold\"";
                $holdPropagateCommand =~ s/WAIT_TAG/$jobName/g;
                my $acmd = "$holdPropagateCommand $sgePropHold";
                system($acmd) and print STDERR "WARNING: Failed to reset hold_jid trigger on '$sgePropHold'.\n";
            }
        } else {
            print STDERR "WARNING: Failed to reset hold '$sgePropHold', not supported on current grid environment.\n";
        }
    }

    exit(0);
}


sub caFailure ($$) {
    my  $msg = shift @_;
    my  $log = shift @_;

    print STDERR "================================================================================\n";
    print STDERR "\n";
    print STDERR "runCA failed.\n";
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
        system("tail -n 50 $log");
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

    if (getGlobal('showNext')) {
        print STDERR "----------------------------------------NEXT-COMMAND\n$cmd\n";
        exit(0);
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

sub dumpInfo {
    if (! -e "$wrk/$asm.gkpStore.info") {
        $cmd  = "$bin/gatekeeper ";
        $cmd .= " -dumpinfo -tabular $wrk/$asm.gkpStore ";
        $cmd .= ">  $wrk/$asm.gkpStore.info ";
        $cmd .= "2> $wrk/$asm.gkpStore.info.err";

        if (runCommand($wrk, $cmd)) {
            caWarn("gatekeeper -dumpinfo failed", "$wrk/$asm.gkpStore.info.err");
        }
    }
}

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
                ($frg =~ m/^(.*)\.sff.bz2$/) ||
                ($frg =~ m/^(.*)\.sff.xz$/)) {
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
        rename "$wrk/$asm.gkpStore.BUILDING.info",        "$wrk/$asm.gkpStore.info";
        rename "$wrk/$asm.gkpStore.BUILDING.errorLog",    "$wrk/$asm.gkpStore.errorLog";
        rename "$wrk/$asm.gkpStore.BUILDING.fastqUIDmap", "$wrk/$asm.gkpStore.fastqUIDmap";
    }

    #  This should never run; the info is now created by gatekeeper.
    dumpInfo();

    generateVectorTrim();

    my $vi = getGlobal("vectorIntersect");

    if ((defined($vi)) && (! -e "$wrk/$asm.gkpStore/$asm.vectorClearLoaded.log")) {
        $cmd  = "$bin/gatekeeper -a -v $vi -o $wrk/$asm.gkpStore ";
        $cmd .= "  > $wrk/$asm.gkpStore/$asm.vectorClearLoaded.log";
        $cmd .= " 2> $wrk/$asm.gkpStore/$asm.vectorClearLoaded.err";

        if (runCommand($wrk, $cmd)) {
            rename "$wrk/$asm.gkpStore/$asm.vectorClearLoaded.log", "$wrk/$asm.gkpStore/$asm.vectorClearLoaded.log.FAILED";
            caFailure("gatekeeper failed to update clear ranges", "$wrk/$asm.gkpStore/$asm.vectorClearLoaded.err");
        }

        unlink "$wrk/$asm.gkpStore/$asm.vectorClearLoaded.err";
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

#
#  Parameters
#

my $MIN_COVERAGE      = 1;  #  Should be 2 if there are "fake" reads in ace file

my $MIN_READS         = 4;
my $MIN_CONTIG_SIZE   = 600;

my $SHRED_READ_LENGTH = 600;

my $LOW_QUAL_DIVISOR  = 4;
my $DEFAULT_QUAL      = 3;

#
#  Methods for reading an ACE file.
#

sub read_AS{
    my $fh=shift;

    while(<$fh>){
        chomp;
        my ($id, $num_contigs, $num_reads)=split /\s+/;
        if($id eq "AS"){
            return ($num_contigs, $num_reads);
        }
    }
    die "Could not find AS to read.\n";
}


sub read_CO{
    my $fh=shift;

    while(<$fh>){
        chomp;
        my ($id, $contig_id, $num_bases, $num_reads, $num_segments, $complementation, $sequence)=split /\s+/;

        if($id eq "CO"){
            while(<$fh>){
                chomp;
                if($_ eq ""){
                    last;
                }else{
                    $sequence.=$_;
                }
            }
            return($contig_id, $num_bases, $num_reads, $num_segments, $complementation, $sequence);
        }
    }
    die "Could not find CO to read.\n";
}


sub read_BQ{
    my $fh=shift;

    my ($id, $sequence);

    while(<$fh>){
        chomp;
        ($id)=split /\s+/;

        if($id eq "BQ"){
            while(<$fh>){
                chomp;
                if($_ eq ""){
                    last;
                }else{
                    $sequence.=$_;
                }
            }
            return($sequence);
        }
    }
    die "Could not find BQ to read.\n";
}


sub read_AF{
    my $fh=shift;

    while(<$fh>){
        chomp;
        my ($id, $read_id, $complementation, $start)=split /\s+/;
        if($id eq "AF"){
            return($read_id, $complementation, $start);
        }
    }
    die "Could not find AF to read.\n";
}


sub read_BS{
    my $fh=shift;

    while(<$fh>){
        chomp;
        my ($id, $start, $end, $read_id)=split /\s+/;
        if($id eq "BS"){
            return($start, $end, $read_id);
        }
    }
    die "Could not find BS to read.\n";
}


sub read_RD{
    my $fh=shift;

    while(<$fh>){
        chomp;
        my ($id, $read_id, $num_bases, $num_read_info_items, $num_read_tags)=split /\s+/;
        my $sequence;
        if($id eq "RD"){
            while(<$fh>){
                chomp;
                if($_ eq ""){
                    last;
                }else{
                    $sequence.=$_;
                }
            }
            return($read_id, $num_bases, $num_read_info_items, $num_read_tags, $sequence);
        }
    }
    die "Could not find RD to read.\n";
}


sub read_QA{
    my $fh=shift;

    while(<$fh>){
        chomp;
        my ($id, $qual_start, $qual_end, $clip_start, $clip_end)=split /\s+/;
        if($id eq "QA"){
            return($qual_start, $qual_end, $clip_start, $clip_end);
        }
    }
    die "Could not find QA to read.\n";
}


sub read_DS{
    my $fh=shift;
    my $id;
    while(<$fh>){
        chomp;
        my ($id)=split /\s+/;
        if($id eq "DS"){
            return("not implemented");
        }
    }
    die "Could not find DS to read.\n";
}

#
#
#

sub emitFragment ($$$$) {
    my $uid = shift;
    my $lid = shift;
    my $seq = shift;
    my $oh  = shift;

    my $len = length($seq);

    my $qvs = $seq;

    my $q = chr($DEFAULT_QUAL                   + ord("0"));
    my $l = chr($DEFAULT_QUAL/$LOW_QUAL_DIVISOR + ord("0"));

    $qvs =~ s/[^ACGT]/$l/og;
    $qvs =~ s/[ACGT]/$q/og;

    print $oh "{FRG\n";
    print $oh "act:A\n";
    print $oh "acc:$uid\n";
    print $oh "rnd:1\n";
    print $oh "sta:G\n";
    print $oh "lib:$lid\n";
    print $oh "pla:0\n";
    print $oh "loc:0\n";
    print $oh "src:\n.\n";
    print $oh "seq:\n$seq\n.\n";
    print $oh "qlt:\n$qvs\n.\n";
    print $oh "hps:\n.\n";
    print $oh "clr:0,$len\n";
    print $oh "}\n";
}

#
#
#

sub shredContig ($$$$$) {
    my $ctgId       = shift;
    my $avgCoverage = shift;
    my $sequence    = shift;
    my $libId       = shift;
    my $oh          = shift;

    my $seq_len=length($sequence);

    my @begin_shred;
    my @end_shred;

    {
        #
        #                  |*******|
        #                  |###############|
        # |-------------------------------------------------|
        #  ----------------1----------------
        #          ----------------2----------------
        #                  ----------------3----------------
        #
        # #### represents the distance between center of read 1 and read 3
        #            [$center_range_width]
        #       **** represents the distance between centers of consective reads
        #            [$center_increments]
        #

        my $shred_len = $SHRED_READ_LENGTH;
        $shred_len = $seq_len - 50 if $seq_len < $SHRED_READ_LENGTH;

        my $num_reads=int($seq_len * $avgCoverage / $shred_len);
        my $center_range_width = $seq_len - $shred_len;

        if($num_reads==1){
            push @begin_shred, 0;
            push @end_shred, $shred_len;
        }else{
            my $center_increments = $center_range_width / ($num_reads-1);

            # Cap the number of reads we will make so that we don't get
            # redundant reads

            my $i;
            my ($prev_begin, $prev_end)=(-1,-1);
            for($i=0; $i<$num_reads; $i++){
                my $begin=$center_increments*$i;
                my $end=$begin+$shred_len;

                $begin=int($begin);
                $end=int($end);

                if($begin!=$prev_begin || $end!=$prev_end){
                    push @begin_shred, $begin;
                    push @end_shred, $end;
                    $prev_begin=$begin;
                    $prev_end=$end;
                }
            }
        }

    }

    my $num_shreds = scalar(@begin_shred);

    my $accomplished_coverage = $num_shreds * $SHRED_READ_LENGTH / $seq_len;

    # Output sequence after it has been formatted to the specified width
    my $shred_idx;
    for($shred_idx=0; $shred_idx<$num_shreds; $shred_idx++){
        my $shredded_sequence=substr($sequence,
                                     $begin_shred[$shred_idx],
                                     $end_shred[$shred_idx]-$begin_shred[$shred_idx]);

        #"/contig=$contigID\.$shred_idx " ,
        #"/target_coverage=$avgCoverage " ,
        #"/accomplished_coverage=$accomplished_coverage " ,
        #"/input_length=$seq_len " ,
        #"/range=${$begin_shred_ref}[$shred_idx]-" ,
        #       "${$end_shred_ref}[$shred_idx]\n";

        emitFragment("$libId.$ctgId.frag$shred_idx.$begin_shred[$shred_idx]-$end_shred[$shred_idx]", $libId, $shredded_sequence, $oh);
    }
}

#
#  Main
#

sub shredACE ($$) {
    my $aceFile = shift;
    my $outFile = shift;
    my $libId   = $aceFile;

    if ($aceFile =~ m/^.*\/(.*).ace/) {
        $libId = $1;
    }

    my $fh = new FileHandle "< $aceFile";
    my $oh = new FileHandle "> $outFile";

    print $oh "{VER\n";
    print $oh "ver:2\n";
    print $oh "}\n";
    print $oh "{LIB\n";
    print $oh "act:A\n";
    print $oh "acc:$libId\n";
    print $oh "ori:U\n";
    print $oh "mea:0.0\n";
    print $oh "std:0.0\n";
    print $oh "src:\n";
    print $oh ".\n";
    print $oh "nft:1\n";
    print $oh "fea:\n";
    print $oh "doNotOverlapTrim=1\n";
    print $oh ".\n";
    print $oh "}\n";

    my ($num_contigs, $num_reads)=read_AS($fh);

    my $contig_idx;
    for($contig_idx=0; $contig_idx<$num_contigs; $contig_idx++){

        my %read_position_hash;

        my ($contig_id, $num_consensus_bases, $num_reads, $num_segments, $complementation, $consensus_sequence) = read_CO($fh);

        my @coverage_array;
        my $i;

        # Initialize Coverage Array
        for($i=0; $i<$num_consensus_bases; $i++){
            $coverage_array[$i]=0;
        }

        my $quality=read_BQ($fh);

        my $read_idx;
        for($read_idx=0; $read_idx<$num_reads; $read_idx++){
            my ($read_id, $complementation, $consensus_start_pos)=read_AF($fh);
            $read_position_hash{$read_id}=$consensus_start_pos;
        }

        my ($base_line_start, $base_line_end, $base_line_read_id)=read_BS($fh);

        for($read_idx=0; $read_idx<$num_reads; $read_idx++){
            my ($read_id, $num_padded_bases, $num_read_info_items, $num_read_tags, $read_sequence)= read_RD($fh);
            my ($qual_start, $qual_end, $align_start, $align_end)=read_QA($fh);
            my $startPos = $read_position_hash{$read_id};

            my $begin = $align_start + $startPos - 1;
            my $end   = $align_end   + $startPos - 1;

            for($i=$begin; $i<$end; $i++){
                $coverage_array[$i]++;
            }
            my ($null)=read_DS($fh);
        }


        my $in_deep_enough=0;
        my @sub_contig_begin_arr;
        my @sub_contig_end_arr;

        # Keep track of where we go into deep coverage region from low coverage regions
        for($i=0; $i<$num_consensus_bases; $i++){
            if($coverage_array[$i]>$MIN_COVERAGE && !$in_deep_enough){
                push @sub_contig_begin_arr, $i;
                $in_deep_enough=1;
            }
            if($coverage_array[$i]<=$MIN_COVERAGE && $in_deep_enough){
                push @sub_contig_end_arr, ($i);
                $in_deep_enough=0;
            }
        }

        if($in_deep_enough){
            push @sub_contig_end_arr, ($i);
        }

        for($i=0; $i<=$#sub_contig_begin_arr; $i++){
            # Sum up coverage for each sub contig
            my $cov_idx;
            my $cov_sum=0;
            for($cov_idx=$sub_contig_begin_arr[$i];
                $cov_idx<$sub_contig_end_arr[$i];
                $cov_idx++){
                $cov_sum+=$coverage_array[$cov_idx];
            }

            # Compute average coverage depth

            my $sub_seq_len=$sub_contig_end_arr[$i]-$sub_contig_begin_arr[$i];
            my $avg_cov = $cov_sum / $sub_seq_len;

            if($num_reads > $MIN_READS && $sub_seq_len>=$MIN_CONTIG_SIZE){
                my $sub_contig_seq  = substr($consensus_sequence,
                                             $sub_contig_begin_arr[$i],
                                             $sub_seq_len);

                # Remove padding
                $sub_contig_seq=~s/\*//g;

                shredContig($contig_id, $avg_cov, $sub_contig_seq, $libId, $oh);
            }
        }
    }

    print $oh "{VER\n";
    print $oh "ver:1\n";
    print $oh "}\n";
}

#
#  For standalone use
#

#die "usage: $0 file.ace > file.frg\n" if (scalar(@ARGV) == 0);
#shredACE($ARGV[0], "a.frg");
#exit();

################################################################################
################################################################################
################################################################################


sub runMeryl ($$$$$$$$) {
    my $merSize      = shift @_;
    my $merComp      = shift @_;
    my $merCanonical = shift @_;
    my $merThresh    = shift @_;
    my $merScale     = 1.0;
    my $merDistinct  = shift @_;
    my $merTotal     = shift @_;
    my $merType      = shift @_;
    my $merDump      = shift @_;

    system("mkdir $wrk/0-mercounts") if (! -d "$wrk/0-mercounts");

    #  This is the historical default.
    if (!defined($merThresh) && !defined($merDistinct) && !defined($merTotal)) {
        $merThresh = "auto";
    }

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

    #  And a special case; if the threshold is zero, we can skip the rest.
    if (($merThresh ne "auto") && ($merThresh == 0) && (!defined($merDistinct)) && (!defined($merTotal))) {
        touch $ffile;
        return;
    }

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
        }

        if (defined($merDistinct) || defined($merTotal)) {
            if (! -e "$ofile.histogram") {
                $cmd  = "$bin/meryl ";
                $cmd .= " -Dh -s $ofile ";
                $cmd .= " > $ofile.histogram ";
                $cmd .= "2> $ofile.histogram.err";

                stopBefore("meryl", $cmd);

                if (runCommand("$wrk/0-mercounts", $cmd)) {
                    rename "$ofile.histogram", "$ofile.histogram.FAILED";
                    caFailure("meryl histogram failed", "$ofile.histogram.err");
                }
            }

            open(F, "< $ofile.histogram") or caFailure("failed to read mer histogram from '$ofile.histogram'", undef);
            while (<F>) {
                my ($threshold, $num, $distinct, $total) = split '\s+', $_;

                if (($merThresh > 0) && ($merThresh < $threshold)) {
                    print STDERR "Supplied merThreshold $merThresh is the smallest.\n";
                    last;
                }

                if ((defined($merDistinct)) && ($merDistinct <= $distinct)) {
                    $merThresh = (($merThresh > 0) && ($merThresh < $threshold)) ? $merThresh : $threshold;
                    print STDERR "Supplied merDistinct $merDistinct with threshold $threshold is the smallest.\n";
                    last;
                }

                if ((defined($merTotal)) && ($merTotal <= $total)) {
                    $merThresh = (($merThresh > 0) && ($merThresh < $threshold)) ? $merThresh : $threshold;
                    print STDERR "Supplied merTotal $merTotal with threshold $threshold is the smallest.\n";
                    last;
                }
            }
            close(F);
        }

        #if ($merThresh == 0) {
        #    caFailure("failed to estimate a mer threshold", "$ofile.estMerThresh.err");
        #}

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
        $obtT = runMeryl(getGlobal('obtMerSize'), $obtc, $obtC, getGlobal("obtMerThreshold"), getGlobal("obtMerDistinct"), getGlobal("obtMerTotal"), "obt", $obtD) if (getGlobal("doOverlapBasedTrimming"));
    } else {
        print STDERR "No need to run meryl for OBT (OBT is disabled).\n"             if (getGlobal("doOverlapBasedTrimming") == 0);
        print STDERR "No need to run meryl for OBT ($asm.nmers.obt.fasta exists).\n" if (-e "$wrk/0-mercounts/$asm.nmers.obt.fasta");
    }

    if ((getGlobal('ovlOverlapper') eq 'mer') || (! -e "$wrk/0-mercounts/$asm.nmers.ovl.fasta")) {
        $ovlT = runMeryl(getGlobal('ovlMerSize'), $ovlc, $ovlC, getGlobal("ovlMerThreshold"), getGlobal("ovlMerDistinct"), getGlobal("ovlMerTotal"), "ovl", $ovlD);
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

    stopAfter("meryl");
}

################################################################################
################################################################################
################################################################################

sub getUMDOverlapperClearRange ($) {
    my $dir     = shift @_;
    my $fileName = "$asm.obtClrRange";

    open(F, "ls -1 -d $wrk/$dir/*overlapperRunDir* |");
    open(G, ">$wrk/$dir/$fileName") or caFailure("failed to write '$wrk/$dir/$fileName'", undef);
    while (<F>) {
        chomp;

        open(T, "< $_/revisedOrigTrimsForReads.txt") or caFailure("failed to open '$_/revisedOrigTrimsForReads.txt'", undef);
        while (<T>) {
            my @trimData = split(/\s+/,$_);
            my $uid = $trimData[0];
            my $bgn = $trimData[1];
            my $end = $trimData[2];

            if ($bgn < $end) {
                print G "frg uid $uid obt all $bgn $end\n";
            } else {
                print G "frg uid $uid obt all $end $bgn\n";
            }
        }
        close(T);
    }
    close(F);
    close(G);

    return $fileName;
}

sub UMDoverlapper () {
    goto alldone if (-d "$wrk/$asm.ovlStore");
    goto alldone if (getGlobal("ovlOverlapper") ne "umd");

    my $outDir  = "1-overlapper";
    system("mkdir $wrk/$outDir") if (! -d "$wrk/$outDir");

    my $jobID = "0000001";
    system("mkdir $wrk/$outDir/$jobID") if (! -d "$wrk/$outDir/$jobID");

    my $vi = getGlobal("vectorIntersect");

    #dump the frag file from gkp if it does not exist already
    # should check if vector clear then dump vec range else dump this range
    if (defined($vi)) {
        if (runCommand($wrk, "$bin/gatekeeper -clear VEC -dumpfrg $wrk/$asm.gkpStore 2> $wrk/gatekeeper.err | grep -v 'No source' > $wrk/$asm.vec.frg")) {
            caFailure("failed to dump gatekeeper store for UMD overlapper", "$wrk/gatekeeper.err");
        }
    }
    elsif ( ! -s "$wrk/$asm.frg" ) {
        if (runCommand($wrk, "$bin/gatekeeper -dumpfrg $wrk/$asm.gkpStore 2> $wrk/gatekeeper.err | grep -v 'No source' > $wrk/$asm.frg")) {
            caFailure("failed to dump gatekeeper store for UMD overlapper", "$wrk/gatekeeper.err");
        }
    }

    # create a job list (we have only one job for right now)
    open(SUB, "> $wrk/$outDir/ovljobs.dat") or caFailure("failed to open '$wrk/$outDir/ovljobs.dat'", undef);
    print SUB "$jobID ";   print SUB "\n";
    print SUB "$jobID ";   print SUB "\n";
    close(SUB);

    # run frg file command
    #
    $cmd  = "$bin/runUMDOverlapper ";
    $cmd .= getGlobal("umdOverlapperFlags") . " ";

    # when we have vector clear, pass it to the overlapper, otherwise tell the overlapper to figure it out
    if (defined($vi)) {
        $cmd .= "-vector-trim-file $wrk/$asm.vec.frg $wrk/$asm.vec.frg ";
    } else {
        $cmd .= "-calculate-trims $wrk/$asm.frg ";
    }

    $cmd .= "$wrk/$outDir/$jobID/$asm.umd.frg ";
    $cmd .= " > $wrk/$outDir/$jobID/overlapper.out 2>$wrk/$outDir/$jobID/overlapper.err";

    if (runCommand("$wrk/$outDir", $cmd)) {
        caFailure("failed to run UMD overlapper", "$wrk/$outDir/$jobID/overlapper.err");
    }

    my $trimFile = getUMDOverlapperClearRange($outDir);
    $cmd = "";
    $cmd .= "$bin/gatekeeper --edit ";
    $cmd .= "$wrk/$outDir/$trimFile $wrk/$asm.gkpStore";
    if (runCommand("$wrk/$outDir", $cmd)) {
        caFailure("failed to update OBT trims", "undef");
    }

    # now create the binary overlaps
    $cmd = "";
    $cmd .= "cat $wrk/$outDir/$jobID/$asm.umd.reliable.overlaps | ";
    $cmd .= "awk '{print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6\"\\t\"\$7}' | ";
    $cmd .= "$bin/convertOverlap ";
    $cmd .= "-b -ovldump ";
    $cmd .= " -out $wrk/$outDir/$jobID/$jobID.ovb";
    if (runCommand("$wrk/$outDir", $cmd)) {
        caFailure("failed to create overlaps", undef);
    }

    #cleanup
    rmrf("$asm.vec.frg");

    touch("$wrk/$outDir/$jobID/$jobID.success");
    stopAfter("overlapper");

  alldone:
}

################################################################################
################################################################################
################################################################################

sub getFigaroClearRange ($) {
    my $outDir     = shift @_;
    my $fileName = "$asm.clv";

    # the figaro output is UID,IID CLR_BGN
    # first reformat is as UID CLR_BGN
    runCommand("$wrk/$outDir", "awk '{print substr(\$1, 1, index(\$1, \",\")-1)\" \"\$2}' $wrk/$outDir/$asm.vectorcuts > $wrk/$outDir/$asm.clrBgn");

    # sort by UID and join it together with the read end to form the full vector clear range
    runCommand("$wrk/$outDir", "sort -nk 1 -T $wrk/$outDir $wrk/$outDir/$asm.clrBgn > $wrk/$outDir/$asm.clrBgn.sorted");
    runCommand("$wrk/$outDir", "join $wrk/$outDir/$asm.clrBgn.sorted $wrk/$asm.untrimmed -o 1.1,1.2,2.3 > $wrk/$outDir/$fileName");

    # clean up
    rmrf("$outDir/$asm.clrBgn");
    rmrf("$outDir/$asm.clrBgn.sorted");

    return $fileName;
}

sub generateFigaroTrim($) {
    my $outDir = shift @_;

    return if (-e "$wrk/$outDir/trim.success");

    # run command
    #
    $cmd  = "$bin/figaro ";
    $cmd .= getGlobal("figaroFlags") . " ";
    $cmd .= "-F $wrk/$asm.fasta -P $asm ";
    $cmd .= " > $wrk/$outDir/figaro.out 2>$wrk/$outDir/figaro.err";

    if (runCommand("$wrk/$outDir", $cmd)) {
        caFailure("figaro died", "$wrk/$outDir/figaro.err");
    }

    # update the gkpStore with newly computed clear ranges
    return getFigaroClearRange($outDir);
}

sub getUMDTrimClearRange($) {
    my $outDir = shift @_;
    my $fileName = "$asm.clv";

    # the umd output is CLR_BGN (in the same order as the input)
    # to join it with the UID we first number both the list of UIDs in the fasta file and the CLR_BGN
    runCommand("$wrk/$outDir", "cat $wrk/$asm.fasta | grep \">\" | awk '{print NR\" \"substr(\$1, 2, index(\$1, \",\")-2)}' > $wrk/$outDir/$asm.numberedUids");
    runCommand("$wrk/$outDir", "awk '{print NR\" \"\$0}' $asm.vectorcuts > $asm.numberedCuts");

    # now we join them together
    runCommand("$wrk/$outDir", "join $wrk/$outDir/$asm.numberedUids $wrk/$outDir/$asm.numberedCuts -o 1.2,2.2 > $wrk/$outDir/$asm.clrBgn");

    # now we can join together the UID CLR_BGN with the read-end information for the full clear range
    runCommand("$wrk/$outDir", "sort -nk 1 -T $wrk/$outDir $wrk/$outDir/$asm.clrBgn > $wrk/$outDir/$asm.clrBgn.sorted");
    runCommand("$wrk/$outDir", "join $wrk/$outDir/$asm.clrBgn.sorted $wrk/$asm.untrimmed -o 1.1,1.2,2.3 > $wrk/$outDir/$fileName");

    # clean up
    rmrf("$outDir/$asm.numberedUids");
    rmrf("$outDir/$asm.numberedCuts");
    rmrf("$outDir/$asm.clrBgn");
    rmrf("$outDir/$asm.clrBgn.sorted");
    rmrf("$outDir/vectorTrimIntermediateFile001.*");

    return $fileName;
}

sub generateUMDTrim($) {
    my $outDir = shift @_;

    return if (-e "$wrk/$outDir/trim.success");

    # run command
    #
    $cmd  = "$bin/dataWorkReduced/findVectorTrimPoints.perl ";
    $cmd .= "$wrk/$asm.fasta $wrk/$outDir/$asm.vectorcuts ";
    $cmd .= " > $wrk/$outDir/umd.out 2>$wrk/$outDir/umd.err";

    if (runCommand("$wrk/$outDir", $cmd)) {
        caFailure("UMD overlapper dataWorkReduced/findVectorTrimPoints.perl died",
                  "$wrk/$outDir/umd.err");
    }

    return getUMDTrimClearRange($outDir);
}

sub generateVectorTrim ($) {
    my $vi = getGlobal("vectorIntersect");
    my $trimmer = getGlobal("vectorTrimmer");
    my $outDir  = "0-preoverlap";
    my $trimFile = undef;

    # when vector insersect is specified or no external trimming is requested, do nothing
    return if (defined($vi));
    return if ($trimmer eq "ca");
    return if (-e "$wrk/$outDir/trim.success");

    #dump the fasta file from gkp
    if ( ! -e "$wrk/$asm.fasta" ) {
        if (runCommand($wrk, "$bin/gatekeeper -dumpfastaseq -clear UNTRIM $wrk/$asm.gkpStore 2> $wrk/$outDir/gatekeeper.err > $wrk/$asm.fasta")) {
            caFailure("failed to dump gatekeeper store for figaro trimmer",
                      "$wrk/$outDir/gatekeeper.err");
        }
    }
    #dump the clr range
    if ( ! -e "$wrk/$asm.untrimmed" ) {
        if (runCommand($wrk, "$bin/gatekeeper -dumpfragments -tabular -clear UNTRIM $wrk/$asm.gkpStore 2> $wrk/$outDir/gatekeeper.err | grep -v 'UID' |awk '{print \$1\" \"\$12\" \"\$13}' | sort -nk 1 -T $wrk/ > $wrk/$asm.untrimmed")) {
            caFailure("failed to dump gatekeeper quality trim points for figaro trimmer",
                      "$wrk/$outDir/gatekeeper.err");
        }
    }

    if ($trimmer eq "figaro") {
        $trimFile = generateFigaroTrim($outDir);
    } elsif($trimmer eq "umd") {
        $trimFile = generateUMDTrim($outDir);
    } else {
        caFailure("unknown vector trimmer $trimmer", undef);
    }

    # set the global vector trim file so that the subsequent code will update the gkp for us
    setGlobal("vectorIntersect", "$wrk/$outDir/$trimFile");

    #cleanup
    rmrf("$asm.fasta");
    rmrf("$asm.untrimmed");

    touch("$wrk/$outDir/trim.success");

    return;
}

################################################################################
################################################################################
################################################################################

sub findMBTFailures ($) {
    my $mbtJobs  = shift @_;
    my $failures = 0;

    for (my $i=1; $i<=$mbtJobs; $i++) {
        my $jobid = substr("0000" . $i, -4);
        if (-e "$wrk/0-mertrim/$asm.$jobid.merTrim.WORKING") {
            #print STDERR "FAILURE $wrk/0-mertrim/$asm.$jobid.merTrim.WORKING\n";
            $failures++;
        }
    }

    return $failures;
}

sub findMBTSuccess ($) {
    my $mbtJobs  = shift @_;
    my $successes= 0;

    for (my $i=1; $i<=$mbtJobs; $i++) {
        my $jobid = substr("0000" . $i, -4);
        if (-e "$wrk/0-mertrim/$asm.$jobid.merTrim") {
            #print STDERR "SUCCESS $wrk/0-mertrim/$asm.$jobid.merTrim\n";
            $successes++;
        }
    }

    return($successes == $mbtJobs);
}


sub merTrim {

    return if (getGlobal("doOverlapBasedTrimming") == 0);
    return if (getGlobal("ovlOverlapper") eq "umd");

    #  Skip mer based trimming if it is done, or if the ovlStore already exists.
    #
    goto alldone if (-e "$wrk/0-mertrim/mertrim.success");
    goto alldone if (-d "$wrk/$asm.ovlStore");

    system("mkdir $wrk/0-mertrim")         if (! -d "$wrk/0-mertrim");

    #  Decide if any libraries request mer based trimming.  There is a simpler method to get
    #  this (see unitigger.pl), but for unity with overlapTrim.pl, we count the number
    #  of libraries that want MBT.
    #
    my $mbtNeeded = 0;

    open(F, "$bin/gatekeeper -nouid -dumplibraries $wrk/$asm.gkpStore |");
    while (<F>) {
        $mbtNeeded++ if (m/doTrim_initialMerBased.*=.*1/);
    }
    close(F);

    if ($mbtNeeded == 0) {
        touch("$wrk/0-mertrim/mertrim.success");
        goto alldone;
    }

    #
    #  Run mer trim on the grid.
    #

    my $mbtBatchSize = getGlobal("mbtBatchSize");
    my $mbtJobs      = int($numFrags / $mbtBatchSize) + (($numFrags % $mbtBatchSize == 0) ? 0 : 1);

    my $merSize      = getGlobal("obtMerSize");
    my $merComp      = 0;  # getGlobal("merCompression");  --  DOES NOT WORK WITH merTrim

    my $mbtThreads   = getGlobal("mbtThreads");

    my $taskID       = getGlobal("gridTaskID");
    my $submitTaskID = getGlobal("gridArraySubmitID");

    runMeryl($merSize, $merComp, "-C", "auto", undef, undef, "mbt", 0);

    if (! -e "$wrk/0-mertrim/mertrim.sh") {
        open(F, "> $wrk/0-mertrim/mertrim.sh") or caFailure("can't open '$wrk/0-mertrim/mertrim.sh'", undef);
        print F "#!" . getGlobal("shell") . "\n";
        print F "\n";
        print F "perl='/usr/bin/env perl'\n";
        print F "\n";
        print F "jobid=\$$taskID\n";
        print F "if [ x\$jobid = x -o x\$jobid = xundefined -o x\$jobid = x0 ]; then\n";
        print F "  jobid=\$1\n";
        print F "fi\n";
        print F "if [ x\$jobid = x ]; then\n";
        print F "  echo Error: I need $taskID set, or a job index on the command line.\n";
        print F "  exit 1\n";
        print F "fi\n";
        print F "\n";
        print F "jobid=`printf %04d \$jobid`\n";
        print F "minid=`expr \$jobid \\* $mbtBatchSize - $mbtBatchSize + 1`\n";
        print F "maxid=`expr \$jobid \\* $mbtBatchSize`\n";
        print F "\n";
        print F "if [ \$maxid -gt $numFrags ] ; then\n";
        print F "  maxid=$numFrags\n";
        print F "fi\n";
        print F "\n";
        print F "if [ \$minid -gt \$maxid ] ; then\n";
        print F "  echo Job partitioning error -- minid=\$minid maxid=\$maxid.\n";
        print F "  exit\n";
        print F "fi\n";
        print F "\n";
        print F "if [ -e $wrk/0-mertrim/$asm.\$jobid.merTrim ]; then\n";
        print F "  echo Job previously completed successfully.\n";
        print F "  exit\n";
        print F "fi\n";
        print F "\n";

        print F "AS_OVL_ERROR_RATE=" , getGlobal("ovlErrorRate"), "\n";
        print F "AS_CNS_ERROR_RATE=" , getGlobal("cnsErrorRate"), "\n";
        print F "AS_CGW_ERROR_RATE=" , getGlobal("cgwErrorRate"), "\n";
        print F "AS_OVERLAP_MIN_LEN=", getGlobal("ovlMinLen"),    "\n";
        print F "AS_READ_MIN_LEN="   , getGlobal("frgMinLen"),    "\n";
        print F "export AS_OVL_ERROR_RATE AS_CNS_ERROR_RATE AS_CGW_ERROR_RATE AS_OVERLAP_MIN_LEN AS_READ_MIN_LEN\n";

        print F getBinDirectoryShellCode();

        print F "\$bin/merTrim \\\n";
        print F " -b  \$minid \\\n";
        print F " -e  \$maxid \\\n";
        print F " -g  $wrk/$asm.gkpStore \\\n";
        print F " -mCillumina \\\n"  if (getGlobal("mbtIlluminaAdapter") ne "0");
        print F " -mC454 \\\n"       if (getGlobal("mbt454Adapter")      ne "0");
        print F " -t  $mbtThreads \\\n";
        print F " -m  $merSize \\\n";
        print F " -mc $wrk/0-mercounts/$asm-C-ms$merSize-cm$merComp \\\n";
        print F " -o  $wrk/0-mertrim/$asm.\$jobid.merTrim.WORKING \\\n";
        print F " >   $wrk/0-mertrim/$asm.\$jobid.err 2>&1 \\\n";
        print F "&& \\\n";
        print F "mv $wrk/0-mertrim/$asm.\$jobid.merTrim.WORKING $wrk/0-mertrim/$asm.\$jobid.merTrim\n";
        print F "\n";
        print F "exit 0\n";
        close(F);

        system("chmod +x $wrk/0-mertrim/mertrim.sh");
    }

    stopBefore("initialTrim", undef);

    #  Don't try to rerun failures.
    #
    #  FAILUREHELPME
    #
    if (findMBTFailures($mbtJobs) > 0) {
        caFailure("merTrim failed.  See *.err in $wrk/0-mertrim", undef);
    }

    #  Submit to the grid (or tell the user to do it), or just run
    #  things here
    #
    if (findMBTSuccess($mbtJobs) == 0) {
        if (getGlobal("useGrid") && getGlobal("mbtOnGrid")) {
            my $submitCommand  = getGlobal("gridSubmitCommand");
            my $nameOption   = getGlobal("gridNameOption");
            my $outputOption  = getGlobal("gridOutputOption");

            my $sge        = getGlobal("sge");
            my $sgeName    = getGlobal("sgeName");
            my $sgeMerTrim = getGlobal("sgeMerTrim");

            $sgeName = "_$sgeName" if (defined($sgeName));
            my $jobName = getGridArrayName("mbt_$asm$sgeName", $mbtJobs);
            my $arrayOpt = getGridArrayOption("mbt_$asm$sgeName", $mbtJobs);

            my $SGE;
            $SGE  = "$submitCommand $sge $sgeMerTrim $nameOption \"$jobName\" $arrayOpt \\\n";
            $SGE .= "  $outputOption $wrk/0-mertrim/$asm.merTrim.$submitTaskID.sge.err \\\n";
            $SGE .= "  $wrk/0-mertrim/mertrim.sh\n";

            submitBatchJobs($SGE, $jobName);
            exit(0);
        } else {
            for (my $i=1; $i<=$mbtJobs; $i++) {
                my $out = substr("0000" . $i, -4);
                schedulerSubmit("$wrk/0-mertrim/mertrim.sh $i > $wrk/0-mertrim/$asm.$out.err 2>&1");
            }

            schedulerSetNumberOfProcesses(getGlobal("mbtConcurrency"));
            schedulerFinish();
        }
    }

    #  Make sure everything finished ok.
    #
    #  FAILUREHELPME
    #
    if (findMBTFailures($mbtJobs) > 0) {
        caFailure("merTrim failed.  See *.err in $wrk/0-mertrim", undef);
    }


    if (runCommand($wrk, "find -L $wrk/0-mertrim -name \\*.merTrim -print | sort > $wrk/0-mertrim/$asm.merTrim.list")) {
        caFailure("failed to generate a list of all the merTrim results", undef);
    }

    {
        $cmd  = "$bin/merTrimApply \\\n";
        $cmd .= " -g $wrk/$asm.gkpStore \\\n";
        $cmd .= " -L $wrk/0-mertrim/$asm.merTrim.list \\\n";
        $cmd .= " -l $wrk/0-mertrim/$asm.merTrim.log \\\n";
        $cmd .= " > $wrk/0-mertrim/$asm.merTrimApply.err 2>&1";

        if (runCommand($wrk, $cmd)) {
            rename "$wrk/0-mertrim/$asm.merTrim.log", "$wrk/0-mertrim/$asm.merTrim.log.FAILED";
            caFailure("merTrimApply failed", "$wrk/0-mertrim/$asm.merTrimApply.err");
        }
    }

    touch("$wrk/0-mertrim/mertrim.success");

  alldone:
    #stopAfter("overlapBasedTrimming");
    #stopAfter("OBT");
}

################################################################################
################################################################################
################################################################################


sub findOvermerryFailures ($$) {
    my $outDir   = shift @_;
    my $ovmJobs  = shift @_;
    my $failures = 0;

    for (my $i=1; $i<=$ovmJobs; $i++) {
        my $out = substr("0000" . $i, -4);
        if (-e "$wrk/$outDir/seeds/$out.ovm.WORKING.gz") {
            $failures++;
        }
    }

    return $failures;
}

sub findOvermerrySuccess ($$) {
    my $outDir   = shift @_;
    my $ovmJobs  = shift @_;
    my $successes= 0;

    for (my $i=1; $i<=$ovmJobs; $i++) {
        my $out = substr("0000" . $i, -4);
        if (-e "$wrk/$outDir/seeds/$out.ovm.gz") {
            $successes++;
        }
    }

    return($successes == $ovmJobs);
}

sub findOlapFromSeedsFailures ($$) {
    my $outDir   = shift @_;
    my $olpJobs  = shift @_;
    my $failures = 0;

    for (my $i=1; $i<=$olpJobs; $i++) {
        my $out = substr("0000" . $i, -4);
        if (-e "$wrk/$outDir/olaps/$out.ovb.WORKING.gz") {
            $failures++;
        }
    }

    return $failures;
}

sub findOlapFromSeedsSuccess ($$) {
    my $outDir   = shift @_;
    my $olpJobs  = shift @_;
    my $successes= 0;

    for (my $i=1; $i<=$olpJobs; $i++) {
        my $out = substr("0000" . $i, -4);
        if (-e "$wrk/$outDir/olaps/$out.ovb.gz") {
            $successes++;
        }
    }

    return $successes == $olpJobs;
}


sub merOverlapper($) {
    my $isTrim = shift @_;

    return if (-d "$wrk/$asm.ovlStore");
    return if (-d "$wrk/$asm.obtStore") && ($isTrim eq "trim");
    return if (-d "$wrk/$asm.tigStore");

    caFailure("mer overlapper detected no fragments", undef) if ($numFrags == 0);
    caFailure("mer overlapper doesn't know if trimming or assembling", undef) if (!defined($isTrim));

    my ($outDir, $ovlOpt, $merSize, $merComp, $merType, $merylNeeded);

    #  Set directories and parameters for either 'trimming' or 'real'
    #  overlaps.

    if ($isTrim eq "trim") {
        $outDir      = "0-overlaptrim-overlap";
        $ovlOpt      = "-G";
        $merSize     = getGlobal("obtMerSize");
        $merComp     = getGlobal("merCompression");
        $merType     = "obt";
        $merylNeeded = (getGlobal("obtMerThreshold") =~ m/auto/) ? 1 : 0;
    } else {
        $outDir      = "1-overlapper";
        $ovlOpt      = "";
        $merSize     = getGlobal("ovlMerSize");
        $merComp     = getGlobal("merCompression");
        $merType     = "ovl";
        $merylNeeded = (getGlobal("ovlMerThreshold") =~ m/auto/) ? 1 : 0;
    }

    system("mkdir $wrk/$outDir")       if (! -d "$wrk/$outDir");
    system("mkdir $wrk/$outDir/seeds") if (! -d "$wrk/$outDir/seeds");
    system("mkdir $wrk/$outDir/olaps") if (! -d "$wrk/$outDir/olaps");

    #  Make the directory (to hold the corrections output) and claim
    #  that fragment correction is all done.  after this, the rest of
    #  the fragment/overlap correction pipeline Just Works.
    #
    system("mkdir $wrk/3-overlapcorrection") if ((! -d "$wrk/3-overlapcorrection") && ($isTrim ne "trim"));

    my $ovmBatchSize = getGlobal("merOverlapperSeedBatchSize");
    my $ovmJobs      = int($numFrags / $ovmBatchSize) + (($numFrags % $ovmBatchSize == 0) ? 0 : 1);

    my $olpBatchSize = getGlobal("merOverlapperExtendBatchSize");
    my $olpJobs      = int($numFrags / $olpBatchSize) + (($numFrags % $olpBatchSize == 0) ? 0 : 1);

    #  Need mer counts, unless there is only one partition.
    meryl() if (($ovmJobs > 1) || ($merylNeeded));

    my $taskID       = getGlobal("gridTaskID");
    my $submitTaskID = getGlobal("gridArraySubmitID");

    #  Create overmerry and olap-from-seeds jobs
    #
    open(F, "> $wrk/$outDir/overmerry.sh") or caFailure("can't open '$wrk/$outDir/overmerry.sh'", undef);
    print F "#!" . getGlobal("shell") . "\n";
    print F "\n";
    print F "jobid=\$$taskID\n";
    print F "if [ x\$jobid = x -o x\$jobid = xundefined -o x\$jobid = x0 ]; then\n";
    print F "  jobid=\$1\n";
    print F "fi\n";
    print F "if [ x\$jobid = x ]; then\n";
    print F "  echo Error: I need $taskID set, or a job index on the command line.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "\n";
    print F "jobid=`printf %04d \$jobid`\n";
    print F "minid=`expr \$jobid \\* $ovmBatchSize - $ovmBatchSize + 1`\n";
    print F "maxid=`expr \$jobid \\* $ovmBatchSize`\n";
    print F "runid=\$\$\n";
    print F "\n";
    print F "if [ \$maxid -gt $numFrags ] ; then\n";
    print F "  maxid=$numFrags\n";
    print F "fi\n";
    print F "if [ \$minid -gt \$maxid ] ; then\n";
    print F "  echo Job partitioning error -- minid=\$minid maxid=\$maxid.\n";
    print F "  exit\n";
    print F "fi\n";
    print F "\n";
    print F "AS_OVL_ERROR_RATE=", getGlobal("ovlErrorRate"), "\n";
    print F "AS_CNS_ERROR_RATE=", getGlobal("cnsErrorRate"), "\n";
    print F "AS_CGW_ERROR_RATE=", getGlobal("cgwErrorRate"), "\n";
    print F "AS_OVERLAP_MIN_LEN=", getGlobal("ovlMinLen"),    "\n";
    print F "AS_READ_MIN_LEN="   , getGlobal("frgMinLen"),    "\n";
    print F "export AS_OVL_ERROR_RATE AS_CNS_ERROR_RATE AS_CGW_ERROR_RATE AS_OVERLAP_MIN_LEN AS_READ_MIN_LEN\n";
    print F "\n";
    print F "if [ ! -d $wrk/$outDir/seeds ]; then\n";
    print F "  mkdir $wrk/$outDir/seeds\n";
    print F "fi\n";
    print F "\n";
    print F "if [ -e $wrk/$outDir/seeds/\$jobid.ovm.gz ]; then\n";
    print F "  echo Job previously completed successfully.\n";
    print F "  exit\n";
    print F "fi\n";
    print F getBinDirectoryShellCode();
    print F "\$bin/overmerry \\\n";
    print F " -g  $wrk/$asm.gkpStore \\\n";
    if ($ovmJobs > 1) {
        print F " -mc $wrk/0-mercounts/$asm-C-ms$merSize-cm$merComp \\\n";
        print F " -tb \$minid -te \$maxid \\\n";
        print F " -qb \$minid \\\n";
    }
    print F " -m $merSize \\\n";
    print F " -c $merComp \\\n";
    print F " -T ", getGlobal("obtMerThreshold"), " \\\n" if ($isTrim eq "trim");
    print F " -T ", getGlobal("ovlMerThreshold"), " \\\n" if ($isTrim ne "trim");
    print F " -t " . getGlobal("merOverlapperThreads") . "\\\n";
    print F " -o $wrk/$outDir/seeds/\$jobid.ovm.WORKING.gz \\\n";
    print F "&& \\\n";
    print F "mv $wrk/$outDir/seeds/\$jobid.ovm.WORKING.gz $wrk/$outDir/seeds/\$jobid.ovm.gz\n";
    close(F);

    system("chmod +x $wrk/$outDir/overmerry.sh");




    open(F, "> $wrk/$outDir/olap-from-seeds.sh") or caFailure("can't open '$wrk/$outDir/olap-from-seeds.sh'", undef);
    print F "#!" . getGlobal("shell") . "\n";
    print F "\n";
    print F "jobid=\$$taskID\n";
    print F "if [ x\$jobid = x -o x\$jobid = xundefined -o x\$jobid = x0 ]; then\n";
    print F "  jobid=\$1\n";
    print F "fi\n";
    print F "if [ x\$jobid = x ]; then\n";
    print F "  echo Error: I need $taskID set, or a job index on the command line.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "\n";
    print F "jobid=`printf %04d \$jobid`\n";
    print F "minid=`expr \$jobid \\* $olpBatchSize - $olpBatchSize + 1`\n";
    print F "maxid=`expr \$jobid \\* $olpBatchSize`\n";
    print F "runid=\$\$\n";
    print F "\n";
    print F "if [ \$maxid -gt $numFrags ] ; then\n";
    print F "  maxid=$numFrags\n";
    print F "fi\n";
    print F "if [ \$minid -gt \$maxid ] ; then\n";
    print F "  echo Job partitioning error -- minid=\$minid maxid=\$maxid.\n";
    print F "  exit\n";
    print F "fi\n";
    print F "\n";
    print F "AS_OVL_ERROR_RATE=", getGlobal("ovlErrorRate"), "\n";
    print F "AS_CNS_ERROR_RATE=", getGlobal("cnsErrorRate"), "\n";
    print F "AS_CGW_ERROR_RATE=", getGlobal("cgwErrorRate"), "\n";
    print F "AS_OVERLAP_MIN_LEN=", getGlobal("ovlMinLen"),    "\n";
    print F "AS_READ_MIN_LEN="   , getGlobal("frgMinLen"),    "\n";
    print F "export AS_OVL_ERROR_RATE AS_CNS_ERROR_RATE AS_CGW_ERROR_RATE AS_OVERLAP_MIN_LEN AS_READ_MIN_LEN\n";
    print F "\n";
    print F "if [ ! -d $wrk/$outDir/olaps ]; then\n";
    print F "  mkdir $wrk/$outDir/olaps\n";
    print F "fi\n";
    print F "\n";
    print F "if [ -e $wrk/$outDir/olaps/\$jobid.ovb.gz ]; then\n";
    print F "  echo Job previously completed successfully.\n";
    print F "  exit\n";
    print F "fi\n";
    print F getBinDirectoryShellCode();
    print F "\$bin/olap-from-seeds \\\n";
    print F " -a -b \\\n";
    print F " -t " . getGlobal("merOverlapperThreads") . "\\\n";
    print F " -S $wrk/$outDir/$asm.merStore \\\n";
    if ($isTrim eq "trim") {
        print F " -G \\\n";  #  Trim only
        print F " -o $wrk/$outDir/olaps/\$jobid.ovb.WORKING.gz \\\n";
        print F " $wrk/$asm.gkpStore \\\n";
        print F " \$minid \$maxid \\\n";
        print F "&& \\\n";
        print F "mv $wrk/$outDir/olaps/\$jobid.ovb.WORKING.gz $wrk/$outDir/olaps/\$jobid.ovb.gz\n";
    } else {
        print F " -w \\\n" if (getGlobal("merOverlapperCorrelatedDiffs"));
        print F " -c $wrk/3-overlapcorrection/\$jobid.frgcorr.WORKING \\\n";
        print F " -o $wrk/$outDir/olaps/\$jobid.ovb.WORKING.gz \\\n";
        print F " $wrk/$asm.gkpStore \\\n";
        print F " \$minid \$maxid \\\n";
        print F "&& \\\n";
        print F "mv $wrk/$outDir/olaps/\$jobid.ovb.WORKING.gz $wrk/$outDir/olaps/\$jobid.ovb.gz \\\n";
        print F "&& \\\n";
        print F "mv $wrk/3-overlapcorrection/\$jobid.frgcorr.WORKING $wrk/3-overlapcorrection/\$jobid.frgcorr\n";
    }
    close(F);

    system("chmod +x $wrk/$outDir/olap-from-seeds.sh");




    if (! -e "$wrk/$outDir/$asm.merStore") {

        #  To prevent infinite loops -- stop now if the overmerry script
        #  exists.  This will unfortunately make restarting from transient
        #  failures non-trivial.
        #
        #  FAILUREHELPME
        #
        if (findOvermerryFailures($outDir, $ovmJobs) > 0) {
            caFailure("overmerry failed.  See *.err in $wrk/$outDir", undef);
        }

        #  Submit to the grid (or tell the user to do it), or just run
        #  things here
        #
        if (findOvermerrySuccess($outDir, $ovmJobs) == 0) {
            if (getGlobal("useGrid") && getGlobal("ovlOnGrid")) {
                my $submitCommand  = getGlobal("gridSubmitCommand");
                my $nameOption   = getGlobal("gridNameOption");
                my $outputOption  = getGlobal("gridOutputOption");

                my $sge        = getGlobal("sge");
                my $sgeName    = getGlobal("sgeName");
                my $sgeOverlap = getGlobal("sgeMerOverlapSeed");

                $sgeName = "_$sgeName" if (defined($sgeName));
                my $jobName = getGridArrayName("mer_$asm$sgeName", $ovmJobs);
                my $arrayOpt = getGridArrayOption("mer_$asm$sgeName", $ovmJobs);

                my $SGE;
                $SGE  = "$submitCommand $sge $sgeOverlap $nameOption \"$jobName\" $arrayOpt \\\n";
                $SGE .= "  $outputOption $wrk/$outDir/seeds/$submitTaskID.err \\\n";
                $SGE .= "  $wrk/$outDir/overmerry.sh\n";

                submitBatchJobs($SGE, $jobName);
                exit(0);
            } else {
                for (my $i=1; $i<=$ovmJobs; $i++) {
                    my $out = substr("0000" . $i, -4);
                    schedulerSubmit("$wrk/$outDir/overmerry.sh $i > $wrk/$outDir/seeds/$out.err 2>&1");
                }

                schedulerSetNumberOfProcesses(getGlobal("merOverlapperSeedConcurrency"));
                schedulerFinish();
            }
        }

        #  Make sure everything finished ok.
        #
        #  FAILUREHELPME
        #
        if (findOvermerryFailures($outDir, $ovmJobs) > 0) {
            caFailure("overmerry failed.  See *.err in $wrk/$outDir", undef);
        }


        if (runCommand($wrk, "find -L $wrk/$outDir/seeds \\( -name \\*ovm.gz -or -name \\*ovm \\) -print > $wrk/$outDir/$asm.merStore.list")) {
            caFailure("failed to generate a list of all the overlap files", undef);
        }

        #$cmd  = "$bin/overlapStore";
        #$cmd .= " -c $wrk/$outDir/$asm.merStore.WORKING";
        #$cmd .= " -g $wrk/$asm.gkpStore";
        #$cmd .= " -M " . getGlobal("ovlStoreMemory");
        #$cmd .= " -L $wrk/$outDir/$asm.merStore.list";
        #$cmd .= " > $wrk/$outDir/$asm.merStore.err 2>&1";

        $cmd  = "$bin/overlapStoreBuild";
        $cmd .= " -o $wrk/$outDir/$asm.merStore.WORKING";
        $cmd .= " -g $wrk/$asm.gkpStore";
        $cmd .= " -M " . getGlobal("ovlStoreMemory");
        $cmd .= " -L $wrk/$outDir/$asm.merStore.list";
        $cmd .= " > $wrk/$outDir/$asm.merStore.err 2>&1";

        if (runCommand($wrk, $cmd)) {
            caFailure("overlap store building failed", "$wrk/$outDir/$asm.merStore.err");
        }

        rename "$wrk/$outDir/$asm.merStore.WORKING", "$wrk/$outDir/$asm.merStore";

        if (getGlobal("saveOverlaps") == 0) {
            open(F, "< $wrk/$outDir/$asm.merStore.list");
            while (<F>) {
                chomp;
                unlink $_;
            }
            close(F);
        }

        rmrf("$outDir/$asm.merStore.list");
        rmrf("$outDir/$asm.merStore.err");
    }


    #  To prevent infinite loops -- stop now if the overmerry script
    #  exists.  This will unfortunately make restarting from transient
    #  failures non-trivial.
    #
    #  FAILUREHELPME
    #
    if (findOlapFromSeedsFailures($outDir, $olpJobs) > 0) {
        caFailure("olap-from-seeds failed.  See *.err in $wrk/$outDir.", undef);
    }

    #  Submit to the grid (or tell the user to do it), or just run
    #  things here
    #
    if (findOlapFromSeedsSuccess($outDir, $olpJobs) == 0) {
        if (getGlobal("useGrid") && getGlobal("ovlOnGrid")) {
            my $submitCommand  = getGlobal("gridSubmitCommand");
            my $nameOption   = getGlobal("gridNameOption");
            my $outputOption  = getGlobal("gridOutputOption");

            my $sge        = getGlobal("sge");
            my $sgeName    = getGlobal("sgeName");
            my $sgeOverlap = getGlobal("sgeMerOverlapExtend");

            $sgeName = "_$sgeName" if (defined($sgeName));
            my $jobName = getGridArrayName("olp_$asm$sgeName", $olpJobs);
            my $arrayOpt = getGridArrayOption("mer_$asm$sgeName", $ovmJobs);

            my $SGE;
            $SGE  = "$submitCommand $sge $sgeOverlap $nameOption \"$jobName\" $arrayOpt \\\n";
            $SGE .= "  $outputOption $wrk/$outDir/olaps/$submitTaskID.err \\\n";
            $SGE .= "  $wrk/$outDir/olap-from-seeds.sh\n";

            submitBatchJobs($SGE, $jobName);
            exit(0);
        } else {
            for (my $i=1; $i<=$olpJobs; $i++) {
                my $out = substr("0000" . $i, -4);
                schedulerSubmit("$wrk/$outDir/olap-from-seeds.sh $i > $wrk/$outDir/olaps/$out.err 2>&1");
            }

            schedulerSetNumberOfProcesses(getGlobal("merOverlapperExtendConcurrency"));
            schedulerFinish();
        }
    }

    #  Make sure everything finished ok.
    #
    #  FAILUREHELPME
    #
    if (findOlapFromSeedsFailures($outDir, $olpJobs) > 0) {
        caFailure("olap-from-seeds failed.  See *.err in $wrk/$outDir.", undef);
    }
}

################################################################################
################################################################################
################################################################################



sub createOverlapJobs($) {
    my $isTrim = shift @_;

    stopAfter("meryl")  if (-d "$wrk/$asm.ovlStore");
    return              if (-d "$wrk/$asm.ovlStore");
    return              if (-d "$wrk/$asm.tigStore");

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

    if (-e "$wrk/$outDir/overlap.sh") {
        stopAfter("meryl");
    }
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

    my $hashLibrary = ($isTrim eq "trim") ? getGlobal("obtHashLibrary") : getGlobal("ovlHashLibrary");
    my $refLibrary  = ($isTrim eq "trim") ? getGlobal("obtRefLibrary")  : getGlobal("ovlRefLibrary");
    my $checkLibrary = ($isTrim eq "trim") ? getGlobal("obtCheckLibrary") : getGlobal("ovlCheckLibrary");

    #  To prevent infinite loops -- stop now if the overlap script
    #  exists.  This will unfortunately make restarting from transient
    #  failures non-trivial.
    #
    #  FAILUREHELPME
    #
    caFailure("overlapper failed\nmanual restart needed to prevent infinite loops\nremove file '$wrk/$outDir/overlap.sh'", undef) if (-e "$wrk/$outDir/overlap.sh");

    meryl();

    my $taskID       = getGlobal("gridTaskID");
    my $submitTaskID = getGlobal("gridArraySubmitID");

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
    print F "jobid=\$$taskID\n";
    print F "if [ x\$jobid = x -o x\$jobid = xundefined -o x\$jobid = x0 ]; then\n";
    print F "  jobid=\$1\n";
    print F "fi\n";
    print F "if [ x\$jobid = x ]; then\n";
    print F "  echo Error: I need $taskID set, or a job index on the command line.\n";
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
    print F "  -H $hashLibrary \\\n" if ($hashLibrary ne "0");
    print F "  -R $refLibrary \\\n"  if ($refLibrary ne "0");
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
        $cmd .= " -H $hashLibrary \\\n" if ($hashLibrary ne "0");
        $cmd .= " -R $refLibrary \\\n"  if ($refLibrary ne "0");
        $cmd .= " -C \\\n" if (!$checkLibrary);
        $cmd .= " -o  $wrk/$outDir \\\n";
        $cmd .= "> $wrk/$outDir/overlap_partition.err 2>&1";

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
        my $submitCommand  = getGlobal("gridSubmitCommand");
        my $nameOption   = getGlobal("gridNameOption");
        my $outputOption  = getGlobal("gridOutputOption");

        my $sge        = getGlobal("sge");
        my $sgeName    = getGlobal("sgeName");
        my $sgeOverlap = getGlobal("sgeOverlap");

        $sgeName = "_$sgeName" if (defined($sgeName));
        my $jobName = getGridArrayName("ovl_$asm$sgeName", $jobs);
        my $arrayOpt = getGridArrayOption("ovl_$asm$sgeName", $jobs);

        my $SGE;
        $SGE  = "$submitCommand $sge $sgeOverlap $nameOption \"$jobName\" $arrayOpt \\\n";
        $SGE .= "  $outputOption $wrk/$outDir/$submitTaskID.out \\\n";
        $SGE .= "  $wrk/$outDir/overlap.sh\n";

        submitBatchJobs($SGE, $jobName);
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

    return  if (-d "$wrk/$asm.tigStore");
    return  if (($isTrim eq "trim") && (-d "$wrk/0-overlaptrim/$asm.obtStore"));
    return  if (($isTrim ne "trim") && (-d "$wrk/$asm.ovlStore"));

    my $outDir = ($isTrim ne "trim") ? "1-overlapper" : "0-overlaptrim-overlap";

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
    goto alldone if (-d "$wrk/$asm.tigStore");

    if (runCommand($wrk, "find -L $wrk/1-overlapper \\( -name \\*ovb.gz -or -name \\*ovb \\) -print > $wrk/$asm.ovlStore.list")) {
        caFailure("failed to generate a list of all the overlap files", undef);
    }

    if (getGlobal("overlapStoreOnGrid") == 1) {
        caFailure("on grid overlapStore construction not supported yet; use runCA-overlapStoreBuild for now", undef);
    }

    $cmd  = "$bin/overlapStoreBuild ";
    $cmd .= " -o $wrk/$asm.ovlStore.BUILDING ";
    $cmd .= " -g $wrk/$asm.gkpStore ";

    if (defined(getGlobal("closureOverlaps"))){
        $cmd .= " -i " . getGlobal("closureOverlaps");
    }

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
    goto alldone if (-d "$wrk/$asm.tigStore");

    system("mkdir $wrk/0-overlaptrim")         if (! -d "$wrk/0-overlaptrim");
    system("mkdir $wrk/0-overlaptrim-overlap") if (! -d "$wrk/0-overlaptrim-overlap");

    #  Disable dedup, unless reads request it.  This avoids an expensive ovlStore build.
    #
    if (getGlobal("doDeDuplication") != 0) {
        setGlobal("doDeDuplication", 0);

        if (system("$bin/gatekeeper -isfeatureset 0 doRemoveDuplicateReads $wrk/$asm.gkpStore") == 0) {
            setGlobal("doDeDuplication", 1);
        }
    }

    #
    #  Do an initial overly-permissive quality trimming, intersected with any known vector trimming.
    #  This step also applies any clear range computed by MBT.
    #

    if ((! -e "$wrk/0-overlaptrim/$asm.initialTrim.log") &&
        (! -e "$wrk/0-overlaptrim/$asm.initialTrim.log.bz2")) {
        $cmd  = "$bin/initialTrim \\\n";
        $cmd .= " -log $wrk/0-overlaptrim/$asm.initialTrim.log \\\n";
        $cmd .= " -frg $wrk/$asm.gkpStore \\\n";
        $cmd .= " >  $wrk/0-overlaptrim/$asm.initialTrim.summary \\\n";
        $cmd .= " 2> $wrk/0-overlaptrim/$asm.initialTrim.err ";

        stopBefore("initialTrim", $cmd);

        if (runCommand("$wrk/0-overlaptrim", $cmd)) {
            rename "$wrk/0-overlaptrim/$asm.initialTrim.log", "$wrk/0-overlaptrim/$asm.initialTrim.log.FAILED";
            caFailure("initial trimming failed", "$wrk/0-overlaptrim/$asm.initialTrim.err");
        }

        unlink "0-overlaptrim/$asm.initialTrim.err";
    }

    #
    #  Decide if any libraries request overlap based trimming -- if all libraries are
    #  asking for mer based trimming, we can skip OBT.
    #

    my $obtNeeded = 0;

    open(F, "$bin/gatekeeper -nouid -dumplibraries $wrk/$asm.gkpStore |");
    while (<F>) {
        $obtNeeded++ if (m/doRemoveDuplicateReads.*=.*1/);
        $obtNeeded++ if (m/doTrim_finalLargestCovered.*=.*1/);
        $obtNeeded++ if (m/doTrim_finalEvidenceBased.*=.*1/);
        $obtNeeded++ if (m/doRemoveSpurReads.*=.*1/);
        $obtNeeded++ if (m/doRemoveChimericReads.*=.*1/);
    }
    close(F);

    if ($obtNeeded == 0) {
        touch("$wrk/0-overlaptrim/overlaptrim.success");
        goto alldone;
    }

    #
    #  Compute overlaps, if we don't have them already
    #

    if (! -e "$wrk/0-overlaptrim/$asm.obtStore") {
        createOverlapJobs("trim");
        checkOverlap("trim");

        #  Sort the overlaps -- this also duplicates each overlap so that
        #  all overlaps for a fragment A are localized.

        if (runCommand("$wrk/0-overlaptrim",
                       "find -L $wrk/0-overlaptrim-overlap \\( -name \\*ovb.gz -or -name \\*ovb \\) -print > $wrk/0-overlaptrim/$asm.obtStore.list")) {
            caFailure("failed to generate a list of all the overlap files", undef);
        }

        $cmd  = "$bin/overlapStoreBuild ";
        $cmd .= " -obt ";
        $cmd .= " -o $wrk/0-overlaptrim/$asm.obtStore.BUILDING ";
        $cmd .= " -g $wrk/$asm.gkpStore ";
        $cmd .= " -M " . getGlobal('ovlStoreMemory');
        $cmd .= " -L $wrk/0-overlaptrim/$asm.obtStore.list";
        $cmd .= " > $wrk/0-overlaptrim/$asm.obtStore.err 2>&1";

        if (runCommand("$wrk/0-overlaptrim", $cmd)) {
            caFailure("failed to build the obt store", "$wrk/0-overlaptrim/$asm.obtStore.err");
        }

        rename "$wrk/0-overlaptrim/$asm.obtStore.BUILDING", "$wrk/0-overlaptrim/$asm.obtStore";

        #  Delete overlaps unless we're told to save them, or we need to dedup.
        if ((getGlobal("saveOverlaps") == 0) && (getGlobal("doDeDuplication") == 0)) {
            open(F, "< $wrk/0-overlaptrim/$asm.obtStore.list");
            while (<F>) {
                chomp;
                unlink $_;
            }
            close(F);
        }

        rmrf("$wrk/0-overlaptrim/$asm.obtStore.list");
        rmrf("$wrk/0-overlaptrim/$asm.obtStore.err");
    }

    #
    #  Deduplicate?
    #

    if ((getGlobal("doDeDuplication") != 0) &&
        (! -e "$wrk/0-overlaptrim/$asm.deduplicate.summary")) {

        if (! -e "$wrk/0-overlaptrim/$asm.dupStore") {
            if (runCommand("$wrk/0-overlaptrim",
                           "find -L $wrk/0-overlaptrim-overlap \\( -name \\*ovb.gz -or -name \\*ovb \\) -print > $wrk/0-overlaptrim/$asm.dupStore.list")) {
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
        $cmd .= "-ovs     $wrk/0-overlaptrim/$asm.obtStore \\\n";
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

    if ((! -e "$wrk/0-overlaptrim/$asm.finalTrim.log") &&
        (! -e "$wrk/0-overlaptrim/$asm.finalTrim.log.bz2")) {
        my $erate  = 0.03;  #  Used for non-Sanger 'largest covered' style
        my $elimit = 4.5;

        my $utg = getUnitigger();

        if      (defined(getGlobal("obtErrorRate"))) {
            $erate  = getGlobal("obtErrorRate");
            $elimit = getGlobal("obtErrorLimit");

        } elsif      ($utg eq "utg") {
            $erate  = getGlobal("utgErrorRate");
            $elimit = getGlobal("utgErrorLimit");

        } elsif ($utg eq "bog") {
            $erate  = getGlobal("utgErrorRate");
            $elimit = getGlobal("utgErrorLimit");

        } elsif ($utg eq "bogart") {
            $erate  = getGlobal("utgGraphErrorRate");
            $elimit = getGlobal("utgGraphErrorLimit");

        } else {
            caFailure("unknown unitigger '$utg' during finalTrim", undef);
        }

        $cmd  = "$bin/finalTrim \\\n";
        $cmd .= "  -G $wrk/$asm.gkpStore \\\n";
        $cmd .= "  -O $wrk/0-overlaptrim/$asm.obtStore \\\n";
        $cmd .= "  -e $erate \\\n";
        $cmd .= "  -E $elimit \\\n"  if (defined($elimit));
        $cmd .= "  -o $wrk/0-overlaptrim/$asm.finalTrim \\\n";
        $cmd .= "> $wrk/0-overlaptrim/$asm.finalTrim.err 2>&1";

        stopBefore("finalTrimming", $cmd);

        if (runCommand("$wrk/0-overlaptrim", $cmd)) {
            unlink "$wrk/0-overlaptrim/$asm.finalTrim.log";
            unlink "$wrk/0-overlaptrim/$asm.finalTrim.stats";
            caFailure("failed to compute final trimming", "$wrk/0-overlaptrim/$asm.finalTrim.err");
        }
    }


    if (getGlobal("doChimeraDetection") ne 'off') {
        if ((! -e "$wrk/0-overlaptrim/$asm.chimera.log") &&
            (! -e "$wrk/0-overlaptrim/$asm.chimera.log.bz2")) {
            my $erate  = getGlobal("ovlErrorRate");
            my $elimit = getGlobal("ovlErrorRate") * 10 + 0.5;

            $cmd  = "$bin/chimera \\\n";
            $cmd .= " -G $wrk/$asm.gkpStore \\\n";
            $cmd .= " -O $wrk/0-overlaptrim/$asm.obtStore \\\n";
            $cmd .= " -e $erate \\\n";
            $cmd .= " -E $elimit \\\n";
            $cmd .= " -o $wrk/0-overlaptrim/$asm.chimera \\\n";
            $cmd .= " -mininniepair 0 -minoverhanging 0 \\\n" if (getGlobal("doChimeraDetection") eq "aggressive");
            $cmd .= " > $wrk/0-overlaptrim/$asm.chimera.err 2>&1";

            stopBefore("chimeraDetection", $cmd);

            if (runCommand("$wrk/0-overlaptrim", $cmd)) {
                rename "$wrk/0-overlaptrim/$asm.chimera.log", "$wrk/0-overlaptrim/$asm.chimera.log.FAILED";
                caFailure("chimera cleaning failed", "$wrk/0-overlaptrim/$asm.chimera.err");
            }
        }
    }

    #rmrf("$asm.obtStore");

    touch("$wrk/0-overlaptrim/overlaptrim.success");

  alldone:
    stopAfter("overlapBasedTrimming");
    stopAfter("OBT");
}

################################################################################
################################################################################
################################################################################



sub overlapCorrection {
    my $cleanup = 1;

    return if (getGlobal("doFragmentCorrection") == 0);

    return if (-e "$wrk/3-overlapcorrection/$asm.erates.updated");
    return if (-e "$wrk/$asm.ovlStore/corrected");
    return if (-d "$wrk/$asm.tigStore");

    system("mkdir $wrk/3-overlapcorrection") if (! -e "$wrk/3-overlapcorrection");

    if (((getGlobal("ovlOverlapper") eq "ovl") &&
         (getGlobal("ovlOverlapper") eq "ovm")) ||
        (! -e "$wrk/3-overlapcorrection/frgcorr.sh")) {
        my $batchSize   = getGlobal("frgCorrBatchSize");
        my $numThreads  = getGlobal("frgCorrThreads");
        my $jobs        = int($numFrags / $batchSize) + (($numFrags % $batchSize == 0) ? 0 : 1);

        my $taskID       = getGlobal("gridTaskID");
        my $submitTaskID = getGlobal("gridArraySubmitID");

        open(F, "> $wrk/3-overlapcorrection/frgcorr.sh") or caFailure("failed to write to '$wrk/3-overlapcorrection/frgcorr.sh'", undef);
        print F "#!" . getGlobal("shell") . "\n\n";
        print F "jobid=\$$taskID\n";
        print F "if [ x\$jobid = x -o x\$jobid = xundefined -o x\$jobid = x0 ]; then\n";
        print F "  jobid=\$1\n";
        print F "fi\n";
        print F "if [ x\$jobid = x ]; then\n";
        print F "  echo Error: I need $taskID set, or a job index on the command line.\n";
        print F "  exit 1\n";
        print F "fi\n";
        print F "\n";
        print F "jobid=`printf %04d \$jobid`\n";
        print F "minid=`expr \$jobid \\* $batchSize - $batchSize + 1`\n";
        print F "maxid=`expr \$jobid \\* $batchSize`\n";
        print F "runid=\$\$\n";
        print F "\n";
        print F "if [ \$maxid -gt $numFrags ] ; then\n";
        print F "  maxid=$numFrags\n";
        print F "fi\n";
        print F "if [ \$minid -gt \$maxid ] ; then\n";
        print F "  echo Job partitioning error -- minid=\$minid maxid=\$maxid.\n";
        print F "  exit\n";
        print F "fi\n";
        print F "\n";
        print F "AS_OVL_ERROR_RATE=", getGlobal("ovlErrorRate"), "\n";
        print F "AS_CNS_ERROR_RATE=", getGlobal("cnsErrorRate"), "\n";
        print F "AS_CGW_ERROR_RATE=", getGlobal("cgwErrorRate"), "\n";
        print F "AS_OVERLAP_MIN_LEN=", getGlobal("ovlMinLen"),    "\n";
        print F "AS_READ_MIN_LEN="   , getGlobal("frgMinLen"),    "\n";
        print F "export AS_OVL_ERROR_RATE AS_CNS_ERROR_RATE AS_CGW_ERROR_RATE AS_OVERLAP_MIN_LEN AS_READ_MIN_LEN\n";
        print F "\n";
        print F "if [ -e $wrk/3-overlapcorrection/\$jobid.frgcorr ] ; then\n";
        print F "  echo Job previously completed successfully.\n";
        print F "  exit\n";
        print F "fi\n";

        print F getBinDirectoryShellCode();

        print F "\$bin/correct-frags \\\n";
        print F "  -t $numThreads \\\n";
        print F "  -S $wrk/$asm.ovlStore \\\n";
        print F "  -o $wrk/3-overlapcorrection/\$jobid.frgcorr.WORKING \\\n";
        print F "  $wrk/$asm.gkpStore \\\n";
        print F "  \$minid \$maxid \\\n";
        print F "&& \\\n";
        print F "mv $wrk/3-overlapcorrection/\$jobid.frgcorr.WORKING $wrk/3-overlapcorrection/\$jobid.frgcorr\n";

        close(F);

        chmod 0755, "$wrk/3-overlapcorrection/frgcorr.sh";

        if (getGlobal("frgCorrOnGrid") && getGlobal("useGrid")) {
            #  Run the correction job on the grid.
            my $submitCommand  = getGlobal("gridSubmitCommand");
            my $nameOption   = getGlobal("gridNameOption");
            my $outputOption  = getGlobal("gridOutputOption");

            my $sge                   = getGlobal("sge");
            my $sgeName               = getGlobal("sgeName");
            my $sgeFragmentCorrection = getGlobal("sgeFragmentCorrection");

            $sgeName = "_$sgeName" if (defined($sgeName));
            my $jobName = getGridArrayName("frg_$asm$sgeName", $jobs);
            my $arrayOpt = getGridArrayOption("ovl_$asm$sgeName", $jobs);

            my $SGE;
            $SGE  = "$submitCommand $sge $sgeFragmentCorrection $nameOption \"$jobName\" $arrayOpt ";
            $SGE .= " $outputOption $wrk/3-overlapcorrection/$submitTaskID.err ";
            $SGE .= "$wrk/3-overlapcorrection/frgcorr.sh\n";

            submitBatchJobs($SGE, $jobName);
            exit(0);
        } else {
            #  Run the correction job right here, right now.

            for (my $i=1; $i<=$jobs; $i++) {
                my $out = substr("0000" . $i, -4);
                schedulerSubmit("$wrk/3-overlapcorrection/frgcorr.sh $i > $wrk/3-overlapcorrection/$out.err 2>&1");
            }

            schedulerSetNumberOfProcesses($global{"frgCorrConcurrency"});
            schedulerFinish();
        }
    }

    #
    #  MERGE CORRECTION
    #

    if (! -e "$wrk/3-overlapcorrection/$asm.frgcorr") {
        my $batchSize  = (getGlobal("ovlOverlapper") eq "mer") ? getGlobal("merOverlapperExtendBatchSize") : getGlobal("frgCorrBatchSize");
        my $jobs       = int($numFrags / $batchSize) + (($numFrags % $batchSize == 0) ? 0 : 1);
        my $failedJobs = 0;

        open(F, "> $wrk/3-overlapcorrection/cat-corrects.frgcorrlist");
        for (my $i=1; $i<=$jobs; $i++) {
            my $jobid = substr("0000" . $i, -4);

            if (! -e "$wrk/3-overlapcorrection/$jobid.frgcorr") {
                print STDERR "Fragment correction job $jobid failed.\n";
                $failedJobs++;
            }

            print F "$wrk/3-overlapcorrection/$jobid.frgcorr\n";
        }
        close(F);

        #  FAILUREHELPME

        if ($failedJobs) {
            if ((getGlobal("ovlOverlapper") eq "ovl") || (getGlobal("ovlOverlapper") eq "ovm")) {
                caFailure("$failedJobs overlap jobs failed; remove $wrk/3-overlapcorrection/frgcorr.sh to try again", undef);
            } else {
                caFailure("$failedJobs overlap jobs failed due to mer overlap seed extension", undef);
            }
        }

        $cmd  = "$bin/cat-corrects ";
        $cmd .= "-L $wrk/3-overlapcorrection/cat-corrects.frgcorrlist ";
        $cmd .= "-o $wrk/3-overlapcorrection/$asm.frgcorr ";
        $cmd .= "> $wrk/3-overlapcorrection/cat-corrects.err 2>&1";

        if (runCommand("$wrk/3-overlapcorrection", $cmd)) {
            rename "$wrk/3-overlapcorrection/$asm.frgcorr", "$wrk/3-overlapcorrection/$asm.frgcorr.FAILED";
            caFailure("failed to concatenate the fragment corrections", "$wrk/3-overlapcorrection/cat-corrects.err");
        }

        if ($cleanup) {
            open(F, "< $wrk/3-overlapcorrection/cat-corrects.frgcorrlist");
            while (<F>) {
                if (m/^(.*)\/([0-9]*).frgcorr/) {
                    #unlink "$1/$2.frgcorr";
                    #unlink "$1/$2.err";
                    my $sge = int($2);
                    #unlink "$1/$sge.err";
                }
            }
            close(F);
            #unlink "$wrk/3-overlapcorrection/cat-corrects.frgcorrlist";
            #unlink "$wrk/3-overlapcorrection/cat-corrects.err";
        }
    }

    #
    #  CREATE OVERLAP CORRECTION
    #

    if (! -e "$wrk/3-overlapcorrection/ovlcorr.sh") {
        my $batchSize  = getGlobal("ovlCorrBatchSize");
        my $jobs       = int($numFrags / $batchSize) + (($numFrags % $batchSize == 0) ? 0 : 1);
        my $taskID       = getGlobal("gridTaskID");
        my $submitTaskID = getGlobal("gridArraySubmitID");

        open(F, "> $wrk/3-overlapcorrection/ovlcorr.sh") or caFailure("failed to write '$wrk/3-overlapcorrection/ovlcorr.sh'", undef);
        print F "jobid=\$$taskID\n";
        print F "if [ x\$jobid = x -o x\$jobid = xundefined -o x\$jobid = x0 ]; then\n";
        print F "  jobid=\$1\n";
        print F "fi\n";
        print F "if [ x\$jobid = x ]; then\n";
        print F "  echo Error: I need $taskID set, or a job index on the command line.\n";
        print F "  exit 1\n";
        print F "fi\n";
        print F "\n";
        print F "if [ \$jobid -gt $jobs ] ; then\n";
        print F "  exit\n";
        print F "fi\n";
        print F "\n";
        print F "jobid=`printf %04d \$jobid`\n";
        print F "frgBeg=`expr \$jobid \\* $batchSize - $batchSize + 1`\n";
        print F "frgEnd=`expr \$jobid \\* $batchSize`\n";
        print F "if [ \$frgEnd -ge $numFrags ] ; then\n";
        print F "  frgEnd=$numFrags\n";
        print F "fi\n";
        print F "frgBeg=`printf %08d \$frgBeg`\n";
        print F "frgEnd=`printf %08d \$frgEnd`\n";

        print F getBinDirectoryShellCode();

        print F "if [ ! -e $wrk/3-overlapcorrection/\$jobid.erate ] ; then\n";

        print F "  \$bin/correct-olaps \\\n";
        print F "    -S $wrk/$asm.ovlStore \\\n";
        print F "    -e $wrk/3-overlapcorrection/\$jobid.erate.WORKING \\\n";
        print F "    $wrk/$asm.gkpStore \\\n";
        print F "    $wrk/3-overlapcorrection/$asm.frgcorr \\\n";
        print F "    \$frgBeg \$frgEnd \\\n";
        print F "  &&  \\\n";
        print F "  mv $wrk/3-overlapcorrection/\$jobid.erate.WORKING $wrk/3-overlapcorrection/\$jobid.erate\n";
        print F "fi\n";
        close(F);

        chmod 0755, "$wrk/3-overlapcorrection/ovlcorr.sh";

        if (getGlobal("ovlCorrOnGrid") && getGlobal("useGrid")) {
            #  Run the correction job on the grid.
            my $submitCommand  = getGlobal("gridSubmitCommand");
            my $nameOption   = getGlobal("gridNameOption");
            my $outputOption  = getGlobal("gridOutputOption");

            my $sge                   = getGlobal("sge");
            my $sgeName               = getGlobal("sgeName");
            my $sgeOverlapCorrection  = getGlobal("sgeOverlapCorrection");

            $sgeName = "_$sgeName" if (defined($sgeName));
            my $jobName = getGridArrayName("ovc_$asm$sgeName", $jobs);
            my $arrayOpt = getGridArrayOption("ovc_$asm$sgeName", $jobs);

            my $SGE;
            $SGE  = "$submitCommand $sge $sgeOverlapCorrection $nameOption \"$jobName\" $arrayOpt ";
            $SGE .= " $outputOption $wrk/3-overlapcorrection/$submitTaskID.err ";
            $SGE .= "$wrk/3-overlapcorrection/ovlcorr.sh\n";

            submitBatchJobs($SGE, $jobName);
            exit(0);
        } else {
            #  Run the correction job right here, right now.

            for (my $i=1; $i<=$jobs; $i++) {
                my $out = substr("0000" . $i, -4);
                schedulerSubmit("$wrk/3-overlapcorrection/ovlcorr.sh $i > $wrk/3-overlapcorrection/$out.err 2>&1");
            }

            schedulerSetNumberOfProcesses($global{"ovlCorrConcurrency"});
            schedulerFinish();
        }
    }

    #
    #  APPLY OVERLAP CORRECTION
    #

    if (! -e "$wrk/3-overlapcorrection/$asm.erates.updated") {
        my $batchSize   = getGlobal("ovlCorrBatchSize");
        my $failedJobs  = 0;
        my $jobs        = int($numFrags / $batchSize) + (($numFrags % $batchSize == 0) ? 0 : 1);

        open(F, "> $wrk/3-overlapcorrection/cat-erates.eratelist");
        for (my $i=1; $i<=$jobs; $i++) {
            my $jobid = substr("0000" . $i, -4);

            if (! -e "$wrk/3-overlapcorrection/$jobid.erate") {
                print STDERR "Overlap correction job $i ($wrk/3-overlapcorrection/$jobid) failed.\n";
                $failedJobs++;
            }

            print F "$wrk/3-overlapcorrection/$jobid.erate\n";
        }
        close(F);

        #  FAILUREHELPME

        if ($failedJobs) {
            caFailure("$failedJobs overlap correction jobs failed; remove $wrk/3-overlapcorrection/ovlcorr.sh (or run by hand) to try again", undef);
        }

        #unlink "$wrk/3-overlapcorrection/$asm.frgcorr" if ($cleanup);

        $cmd  = "$bin/cat-erates ";
        $cmd .= "-L $wrk/3-overlapcorrection/cat-erates.eratelist ";
        $cmd .= "-o $wrk/3-overlapcorrection/$asm.erates ";
        $cmd .= "> $wrk/3-overlapcorrection/cat-erates.err 2>&1";
        if (runCommand("$wrk/3-overlapcorrection", $cmd)) {
            rename "$wrk/3-overlapcorrection/$asm.erates", "$wrk/3-overlapcorrection/$asm.erates.FAILED";
            caFailure("failed to concatenate the overlap erate corrections", "$wrk/3-overlapcorrection/cat-erates.err");
        }

        $cmd  = "$bin/overlapStore ";
        $cmd .= " -u $wrk/$asm.ovlStore ";
        $cmd .= " $wrk/3-overlapcorrection/$asm.erates ";
        $cmd .= "> $wrk/3-overlapcorrection/overlapStore-update-erates.err 2>&1";
        if (runCommand("$wrk/3-overlapcorrection", $cmd)) {
            caFailure("failed to apply the overlap corrections", "$wrk/3-overlapcorrection/overlapStore-update-erates.err");
        }

        touch("$wrk/3-overlapcorrection/$asm.erates.updated");
        touch("$wrk/$asm.ovlStore/corrected");

        if ($cleanup) {
            open(F, "< $wrk/3-overlapcorrection/cat-erates.eratelist");
            while (<F>) {
                if (m/^(.*)\/([0-9]*).erate/) {
                    #unlink "$1/$2.erate";
                    #unlink "$1/$2.err";
                    my $sge = int($2);
                    #unlink "$1/$sge.err";
                }
            }
            close(F);

            #unlink "$wrk/3-overlapcorrection/overlapStore-update-erates.err";
            #unlink "$wrk/3-overlapcorrection/$asm.erates";

            #unlink "$wrk/3-overlapcorrection/cat-erates.err";
            #unlink "$wrk/3-overlapcorrection/cat-erates.eratelist";

            #unlink "$wrk/3-overlapcorrection/frgcorr.sh";
            #unlink "$wrk/3-overlapcorrection/ovlcorr.sh";
        }
    }
}

################################################################################
################################################################################
################################################################################

sub classifyMates () {

    #  This is configuration hell.  For now, we just ask for the library to classify.
    #  We could later extend with parallelization, run time, memory, etc, etc.
    #
    #  classifyIlluminaBB=<libraryname>,<libraryname>,...
    #  classifyIlluminaMP=<libraryname>,<libraryname>,...

    my @libsToClassify = split ',', getGlobal("dncMPlibraries");
    my @backboneToUse  = split ',', getGlobal("dncBBlibraries");

    if ((scalar(@libsToClassify) == 0) ||
        (scalar(@backboneToUse)  == 0) ||
        (-e "$wrk/2-classifyMates/classify.success")) {
        return;
    }

    system("mkdir $wrk/2-classifyMates") if (! -e "$wrk/2-classifyMates");

    #  Load a map of library name to IID.

    my %libToIID;

    open(F, "$bin/gatekeeper -dumplibraries -tabular $wrk/$asm.gkpStore |");
    while (<F>) {
        my @v = split '\s+', $_;
        $libToIID{$v[0]} = $v[1];
    }
    close(F);

    #  Map the backbone library names to a IIDs.

    my $bbIID;

    foreach my $lib (@backboneToUse) {
        if (!exists($libToIID{$lib})) {
            caFailure("Backbone library '$lib' doesn't exist in the assembly, classifyMates failed\n", undef);
        }
        if (defined($bbIID)) {
            $bbIID .= " -bl $libToIID{$lib}";
        } else {
            $bbIID  = "$libToIID{$lib}";
        }
    }

    #  Build jobs to run.  These are NOT optimal, but building optimal configurations is probably
    #  impossible.  You can save memory by classifying each MP library seperately, but doing all at
    #  once is easier.
    #
    #  The jobs are built, but they are only templates.  It is up to the user to run them, and
    #  user could modify the way they are run.

    my @classifyJobs;
    my $classifiedOutputs = "";

    foreach my $lib (@libsToClassify) {
        my $mpIID = $libToIID{$lib};

        if (!exists($libToIID{$lib})) {
            caFailure("Mate pair library '$lib' doesn't exist in the assembly, classifyMates failed\n", undef);
        }

        open(F, "> $wrk/2-classifyMates/classify-$lib.sh");
        print F "#!/bin/sh\n";
        print F "\n";
        print F "\n";
        print F "if [ ! -e \"$wrk/2-classifyMates/classifyMates.sl$mpIID.O.0100.1500.bfs.100000\" ] ; then\n";
        print F "  $bin/classifyMates \\\n";
        print F "    -G $wrk/$asm.gkpStore \\\n";
        print F "    -O $wrk/$asm.ovlStore \\\n";
        print F "    -t 8 \\\n";
        print F "    -m 128 \\\n";
        print F "    -sl $mpIID \\\n";
        print F "    -bl $bbIID \\\n";
        print F "    -outtie -min  100 -max 1500 -bfs 100000 \\\n";
        print F "    -o $wrk/2-classifyMates/classifyMates.sl$mpIID.O.0100.1500.bfs.100000\n";
        print F "fi\n";
        print F "\n";
        print F "#  classifyMates internally appends .WORKING, then renames just before it exits.\n";
        print F "if [ ! -e \"$wrk/2-classifyMates/classifyMates.sl$mpIID.O.0100.1500.bfs.100000\" ] ; then\n";
        print F "  exit 1\n";
        print F "fi\n";
        print F "\n";
        print F "exit 0\n";
        close(F);

        chmod 0755, "$wrk/2-classifyMates/classify-$lib.sh";

        push @classifyJobs, "$wrk/2-classifyMates/classify-$lib.sh";
        $classifiedOutputs .= "  -r $wrk/2-classifyMates/classifyMates.sl$mpIID.O.0100.1500.bfs.100000 \\\n";
    }

    open(F, "> $wrk/2-classifyMates/classify-apply.sh");
    print F "#!/bin/sh\n";
    print F "\n";
    print F "$bin/classifyMatesApply \\\n";
    print F "  -G $wrk/$asm.gkpStore \\\n";
    print F "  -p \\\n";
    print F "$classifiedOutputs";  #  NO INDENTATION!  NO LINE CONTINUATION!  Both are in the string directly.
    print F "  -o $wrk/2-classifyMates/$asm.classified.gkpStore.edit \\\n";
    print F ">  $wrk/2-classifyMates/$asm.classified.log \\\n";
    print F "2> $wrk/2-classifyMates/$asm.classified.summary \\\n";
    print F "&& \\\n";
    print F "$bin/gatekeeper --edit \\\n";
    print F "  $wrk/2-classifyMates/$asm.classified.gkpStore.edit \\\n";
    print F "  $wrk/$asm.gkpStore \\\n";
    print F "> $wrk/2-classifyMates/$asm.classified.gkpStore.edit.out \\\n";
    print F "&& \\\n";
    print F "touch $wrk/2-classifyMates/classify.success\n";
    print F "\n";
    print F "exit 0\n";
    close(F);

    chmod 0755, "$wrk/2-classifyMates/classify-apply.sh";

    #  Force the stop.  User must run scripts by hand.

    stopBefore("classifyMates", undef);

    foreach my $j (@classifyJobs) {
        if (runCommand("$wrk/2-classifyMates", $j)) {
            caFailure("failed to run classify mates command $j", undef);
        }
    }

    if (runCommand("$wrk/2-classifyMates", "$wrk/2-classifyMates/classify-apply.sh > $wrk/2-classifyMates/classify-apply.err 2>&1")) {
        caFailure("failed to apply classify mates results", undef);
    }
    if (! -e "$wrk/2-classifyMates/classify.success") {
        caFailure("classifyMatesApply failed.", "");
    }

    stopAfter("classifyMates");
}




################################################################################
################################################################################
################################################################################

sub unitigger () {

    if (0) {
        $cmd = "$bin/removeMateOverlap -gkp $wrk/$asm.gkpStore -ovl $wrk/$asm.ovlStore";
        if (runCommand("$wrk", $cmd)) {
            caFailure("failed to remove mate overlaps", undef);
        }
    }

    if (-e "$wrk/4-unitigger/unitigger.success") {
        goto alldone;
    }

    if (-e "$wrk/$asm.tigStore") {
        print STDERR "Skipping unitigger because tigStore exists at $wrk/$asm.tigStore\n";
        goto alldone;
    }

    system("mkdir $wrk/4-unitigger") if (! -e "$wrk/4-unitigger");

    my $e    = getGlobal("utgErrorRate");        #  Unitigger and BOG
    my $E    = getGlobal("utgErrorLimit");
    my $eg   = getGlobal("utgGraphErrorRate");   #  BOGART
    my $Eg   = getGlobal("utgGraphErrorLimit");
    my $em   = getGlobal("utgMergeErrorRate");
    my $Em   = getGlobal("utgMergeErrorLimit");

    my $B = int($numFrags / getGlobal("cnsPartitions"));
    $B = getGlobal("cnsMinFrags") if ($B < getGlobal("cnsMinFrags"));

    my $u = getGlobal("utgBubblePopping");

    my $unitigger = getUnitigger();

    if ($unitigger eq "bogart") {
        my $th   = getGlobal("batThreads");
        my $mem  = getGlobal("batMemory");
        my $opts = getGlobal("batOptions");

        $cmd  = "$bin/bogart ";
        $cmd .= " -O $wrk/$asm.ovlStore ";
        $cmd .= " -G $wrk/$asm.gkpStore ";
        $cmd .= " -T $wrk/$asm.tigStore ";
        $cmd .= " -B $B ";
        $cmd .= " -eg $eg ";
        $cmd .= " -Eg $Eg ";
        $cmd .= " -em $em ";
        $cmd .= " -Em $Em ";
        $cmd .= " -threads $th " if (defined($th));
        $cmd .= " -M $mem "      if (defined($mem));
        $cmd .= " $opts "        if (defined($opts));
        $cmd .= " -o $wrk/4-unitigger/$asm ";
        $cmd .= " > $wrk/4-unitigger/unitigger.err 2>&1";
    } elsif ($unitigger eq "bog") {
        my $bmd = getGlobal("bogBadMateDepth");

        $cmd  = "$bin/buildUnitigs ";
        $cmd .= " -O $wrk/$asm.ovlStore ";
        $cmd .= " -G $wrk/$asm.gkpStore ";
        $cmd .= " -T $wrk/$asm.tigStore ";
        $cmd .= " -B $B ";
        $cmd .= " -e $e ";
        $cmd .= " -E $E ";
        $cmd .= " -b "      if (getGlobal("bogBreakAtIntersections") == 1);
        $cmd .= " -m $bmd " if (defined($bmd));
        $cmd .= " -U "      if ($u == 1);
        $cmd .= " -o $wrk/4-unitigger/$asm ";
        $cmd .= " > $wrk/4-unitigger/unitigger.err 2>&1";
    } elsif ($unitigger eq "utg") {
        $cmd  = "$bin/unitigger ";
        $cmd .= " -I $wrk/$asm.ovlStore ";
        $cmd .= " -F $wrk/$asm.gkpStore ";
        $cmd .= " -T $wrk/$asm.tigStore ";
        $cmd .= " -B $B ";
        $cmd .= " -e $e ";
        $cmd .= " -k " if (getGlobal("utgRecalibrateGAR") == 1);
        $cmd .= " -d 1 -x 1 -z 10 -j 5 -U $u ";
        $cmd .= " -o $wrk/4-unitigger/$asm ";
        $cmd .= " > $wrk/4-unitigger/unitigger.err 2>&1";
    } else {
        caFailure("unknown unitigger '$unitigger'; must be 'bog' or 'utg'", undef);
    }

    stopBefore("unitigger", $cmd);

    if (runCommand("$wrk/4-unitigger", $cmd)) {
        caFailure("failed to unitig", "$wrk/4-unitigger/unitigger.err");
    }

    touch("$wrk/4-unitigger/unitigger.success");

  alldone:
    stopAfter("unitigger");
}

################################################################################
################################################################################
################################################################################

sub createPostUnitiggerConsensusJobs (@) {
    my $blasr = findBlasr(getGlobal("consensus"));
    if (!defined($blasr)) { setGlobal("consensus", "cns"); }
    my $consensusType = getGlobal("consensus");

    return if (-e "$wrk/5-consensus/consensus.sh");

    if (! -e "$wrk/5-consensus/$asm.partitioned") {
        $cmd  = "$bin/gatekeeper ";
        $cmd .= " -P $wrk/4-unitigger/$asm.partitioning ";
        $cmd .= " $wrk/$asm.gkpStore ";
        $cmd .= "> $wrk/5-consensus/$asm.partitioned.err 2>&1";
        if (runCommand("$wrk/5-consensus", $cmd)) {
            caFailure("failed to partition the fragStore", "$wrk/5-consensus/$asm.partitioned.err");
        }

        touch "$wrk/5-consensus/$asm.partitioned";
    }

    my $jobs = 0;

    open(F, "< $wrk/4-unitigger/$asm.partitioningInfo") or caFailure("can't open '$wrk/4-unitigger/$asm.partitioningInfo'", undef);
    while (<F>) {
        if (m/Partition\s+(\d+)\s+has\s+(\d+)\s+unitigs\sand\s+(\d+)\s+fragments./) {
            $jobs = $1;
        }
    }
    close(F);

    my $taskID       = getGlobal("gridTaskID");
    open(F, "> $wrk/5-consensus/consensus.sh") or caFailure("can't open '$wrk/5-consensus/consensus.sh'", undef);
    print F "#!" . getGlobal("shell") . "\n";
    print F "\n";
    print F "jobid=\$$taskID\n";
    print F "if [ x\$jobid = x -o x\$jobid = xundefined -o x\$jobid = x0 ]; then\n";
    print F "  jobid=\$1\n";
    print F "fi\n";
    print F "if [ x\$jobid = x ]; then\n";
    print F "  echo Error: I need $taskID set, or a job index on the command line.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "\n";
    print F "if [ \$jobid -gt $jobs ]; then\n";
    print F "  echo Error: Only $jobs partitions, you asked for \$jobid.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "\n";
    print F "jobid=`printf %03d \$jobid`\n";
    print F "\n";
    print F "if [ -e $wrk/5-consensus/${asm}_\$jobid.success ] ; then\n";
    print F "  exit 0\n";
    print F "fi\n";
    print F "\n";
    print F "AS_OVL_ERROR_RATE=", getGlobal("ovlErrorRate"), "\n";
    print F "AS_CNS_ERROR_RATE=", getGlobal("cnsErrorRate"), "\n";
    print F "AS_CGW_ERROR_RATE=", getGlobal("cgwErrorRate"), "\n";
    print F "AS_OVERLAP_MIN_LEN=", getGlobal("ovlMinLen"),    "\n";
    print F "AS_READ_MIN_LEN="   , getGlobal("frgMinLen"),    "\n";
    print F "export AS_OVL_ERROR_RATE AS_CNS_ERROR_RATE AS_CGW_ERROR_RATE AS_OVERLAP_MIN_LEN AS_READ_MIN_LEN\n";

    print F getBinDirectoryShellCode();

    if ($consensusType eq "cns") {
        my $maxCov = getGlobal('cnsMaxCoverage');

        print F "\$bin/utgcns \\\n";
        print F "  -g $wrk/$asm.gkpStore \\\n";
        print F "  -t $wrk/$asm.tigStore 1 \$jobid \\\n";
        print F "  -maxcoverage $maxCov \\\n";
        print F "> $wrk/5-consensus/${asm}_\$jobid.cns.err 2>&1 \\\n";
        print F "&& \\\n";
        print F "\$bin/utgcnsfix \\\n";
        print F "  -g $wrk/$asm.gkpStore \\\n";
        print F "  -t $wrk/$asm.tigStore 2 \$jobid \\\n";
        print F "  -o $wrk/5-consensus/${asm}_\$jobid.fixes \\\n";
        print F "> $wrk/5-consensus/${asm}_\$jobid.fix.err 2>&1";
        print F "&& \\\n";
        print F "touch $wrk/5-consensus/${asm}_\$jobid.success\n";
    } elsif ($consensusType eq "seqan") {
        print F "\$bin/SeqAn_CNS \\\n";
        print F "  -G $wrk/$asm.gkpStore \\\n";
        print F "  -c \$cgbfile \\\n";
        print F "  -s \$bin/graph_consensus \\\n";
        print F "  -w $wrk/5-consensus/ \\\n";
        print F "  -o $wrk/5-consensus/${asm}_\$jobid.cgi \\\n";
        print F " > $wrk/5-consensus/${asm}_\$jobid.cns.err 2>&1 \\\n";
        print F "&& \\\n";
        print F "touch $wrk/5-consensus/${asm}_\$jobid.success\n";
    } elsif ($consensusType eq "pbdagcon" || $consensusType eq "pbutgcns") {
        print F "cat $wrk/4-unitigger/$asm.partitioning |awk -v JOB=\$jobid '{if (\$1 == JOB) print \$NF}' > $wrk/5-consensus/$asm.\$jobid.iid\n";
        print F "$bin/gatekeeper -dumpfasta $wrk/5-consensus/$asm.\$jobid -iid $wrk/5-consensus/$asm.\$jobid.iid $wrk/$asm.gkpStore\n";
        print F "rm -f $wrk/5-consensus/$asm.\$jobid.q*\n";
        print F "\$bin/tigStore -d layout -U -t $wrk/$asm.tigStore 1 -up \$jobid -g $wrk/$asm.gkpStore > $wrk/5-consensus/$asm.\$jobid.lay\n";
        print F "\$bin/convertToPBCNS -path $blasr -consensus $consensusType -coverage 1 -threads " . getGlobal("cnsConcurrency") . " -prefix $wrk/5-consensus/$asm.\$jobid.tmp -length 500 -sequence $wrk/5-consensus/$asm.\$jobid.fasta -input $wrk/5-consensus/$asm.\$jobid.lay -output $wrk/5-consensus/$asm.\$jobid.fa\n";
        print F "\$bin/addCNSToStore -path \$bin -input $wrk/5-consensus/$asm.\$jobid.fa -lay $wrk/5-consensus/$asm.\$jobid.lay -output $wrk/5-consensus/$asm.\$jobid.cns -prefix $wrk/$asm -sequence $wrk/5-consensus/$asm.\$jobid.fasta -partition \$jobid && \$bin/utgcnsfix -g $wrk/$asm.gkpStore  -t $wrk/$asm.tigStore 2 \$jobid -o $wrk/5-consensus/${asm}_\$jobid.fixes > $wrk/5-consensus/${asm}_\$jobid.fix.err 2>&1 && touch $wrk/5-consensus/${asm}_\$jobid.success\n";
        print F "if [ -e $wrk/5-consensus/${asm}_\$jobid.success ]; then\n";
        print F "   rm -f $wrk/5-consensus/${asm}.\$jobid.fasta*\n";
        print F "   rm -f $wrk/5-consensus/${asm}.\$jobid.lay\n";
        print F "fi\n";
        setGlobal("cnsConcurrency", 1);

    } else {
        caFailure("unknown consensus type $consensusType; should be 'cns' or 'seqan'", undef);
    }
    close(F);

    chmod 0755, "$wrk/5-consensus/consensus.sh";

    if (getGlobal("useGrid") && getGlobal("cnsOnGrid")) {
        my $submitCommand  = getGlobal("gridSubmitCommand");
        my $nameOption   = getGlobal("gridNameOption");
        my $outputOption  = getGlobal("gridOutputOption");

        my $sge          = getGlobal("sge");
        my $sgeName      = getGlobal("sgeName");
        my $sgeConsensus = getGlobal("sgeConsensus");

        $sgeName = "_$sgeName" if (defined($sgeName));
        my $jobName = getGridArrayName("utg_$asm$sgeName", $jobs);
        my $arrayOpt = getGridArrayOption("ovc_$asm$sgeName", $jobs);

        my $SGE;
        $SGE  = "$submitCommand $sge $sgeConsensus $nameOption \"$jobName\" $arrayOpt ";
        $SGE .= "$outputOption /dev/null ";
        $SGE .= "$wrk/5-consensus/consensus.sh\n";

        submitBatchJobs($SGE, $jobName);
        exit(0);
    } else {
        for (my $i=1; $i<=$jobs; $i++) {
            schedulerSubmit("$wrk/5-consensus/consensus.sh $i > /dev/null 2>&1");
        }

        schedulerSetNumberOfProcesses(getGlobal("cnsConcurrency"));
        schedulerFinish();
    }
}



sub postUnitiggerConsensus () {
    my $cmd;

    return if (-e "$wrk/5-consensus/consensus.success");

    if (-e "$wrk/$asm.tigStore/seqDB.v005.dat") {
        print STDERR "Skipping consensus because tigStore version 5 exists at $wrk/$asm.tigStore\n";
        return;
    }

    system("mkdir $wrk/5-consensus") if (! -d "$wrk/5-consensus");

    #  NOT IDEAL.  We assume that jobs are finished if the script exists.  We need
    #  better logic in createPostUnitiggerConsensusJobs() to tell if the job crashed.
    #
    if (! -e "$wrk/5-consensus/consensus.sh") {
        createPostUnitiggerConsensusJobs();
    }

    #
    #  Check that consensus finished properly
    #

    my $failedJobs = 0;

    open(F, "< $wrk/4-unitigger/$asm.partitioningInfo") or caFailure("can't open '$wrk/4-unitigger/$asm.partitioningInfo'", undef);
    while (<F>) {
        if (m/Partition\s+(\d+)\s+has\s+(\d+)\s+unitigs\sand\s+(\d+)\s+fragments./) {
            my $id = substr("000" . $1, -3);

            if (! -e "$wrk/5-consensus/${asm}_$id.success") {
                print STDERR "$wrk/5-consensus/${asm}_$id failed -- no .success.\n";
                $failedJobs++;
            }
        }
    }
    close(F);

    #  FAILUREHELPME
    #
    caFailure("$failedJobs unitig consensus jobs failed; remove $wrk/5-consensus/consensus.sh to try again", undef) if ($failedJobs);

    #
    #  Summarize the sampling done by reading logs
    #

    if (! -e "$wrk/5-consensus/$asm.sampling") {
        open(O, "> $wrk/5-consensus/$asm.sampling");
        open(D, "> $wrk/5-consensus/$asm.sampling.dat");

        print D "utgIID\t#_removed\tcov_removed\t#_saved\tcov_saved\t#_kept\tcov_kept\n";

        open(F, "< $wrk/4-unitigger/$asm.partitioningInfo") or caFailure("can't open '$wrk/4-unitigger/$asm.partitioningInfo'", undef);
        while (<F>) {
            if (m/Partition\s+(\d+)\s+has\s+(\d+)\s+unitigs\sand\s+(\d+)\s+fragments./) {
                my $id = substr("000" . $1, -3);

                open(G, "< $wrk/5-consensus/${asm}_$id.cns.err") or warn "Failed to open '$wrk/5-consensus/${asm}_$id.cns.err'\n";
                while (<G>) {
                    s/^\s+//;
                    s/\s+$//;

                    if (m/unitig\s+(\d+)\s+removing\s+(\d+)\s+\((\d+.\d+)x\)\s+contained\s+reads;\s+processing\s+only\s+(\d+)\s+contained\s+\((\d+.\d+)x\)\s+and\s+(\d+)\s+dovetail\s+\((\d+.\d+)x\)\s+reads/) {
                        print O "$_\n";
                        print D "$1\t$2\t$3\t$4\t$5\t$6\t$7\n";
                    }
                }
                close(G);
            }
        }
        close(F);
        close(D);
        close(O);
    }

    #
    #  Apply the utgcnsfix changes
    #

    if (! -e "$wrk/5-consensus/$asm.fixes") {
        open(O, "> $wrk/5-consensus/$asm.fixes");
        open(F, "< $wrk/4-unitigger/$asm.partitioningInfo") or caFailure("can't open '$wrk/4-unitigger/$asm.partitioningInfo'", undef);
        while (<F>) {
            if (m/Partition\s+(\d+)\s+has\s+(\d+)\s+unitigs\sand\s+(\d+)\s+fragments./) {
                my $id = substr("000" . $1, -3);

                open(G, "< $wrk/5-consensus/${asm}_$id.fixes") or die "Failed to open '$wrk/5-consensus/${asm}_$id.fixes'\n";
                while (<G>) {
                    print O $_;
                }
                close(G);
            }
        }
        close(F);
        close(O);

        $cmd  = "$bin/tigStore";
        $cmd .= " -g $wrk/$asm.gkpStore";
        $cmd .= " -t $wrk/$asm.tigStore 2";
        $cmd .= " -N";
        $cmd .= " -R $wrk/5-consensus/$asm.fixes";
        $cmd .= " > $asm.fixes.err 2>&1";

        if (runCommand("$wrk/5-consensus", $cmd)) {
            caFailure("unitig utgcnsfix failed", "$wrk/5-consensus/$asm.fixes.err");
            rename "$wrk/5-consensus/$asm.fixes", "$wrk/5-consensus/$asm.fixes.FAILED";
        }
    }

    #
    #  Estimate insert size estimates
    #
    #  tigStore writes output files using the tigStore as a prefix.  We thus need to temporarily symlink
    #  The tigStore to the proper directory.
    #

    if (! -e "$wrk/5-consensus-insert-sizes/estimates.out") {
        system("mkdir $wrk/5-consensus-insert-sizes") if (! -d "$wrk/5-consensus-insert-sizes");

        system("ln -s ../$asm.tigStore $wrk/5-consensus-insert-sizes/$asm.tigStore");

        $cmd  = "$bin/tigStore \\\n";
        $cmd .= " -g $wrk/$asm.gkpStore \\\n";
        $cmd .= " -t $wrk/5-consensus-insert-sizes/$asm.tigStore 3 \\\n";
        $cmd .= " -d matepair -U \\\n";
        $cmd .= "> $wrk/5-consensus-insert-sizes/estimates.out 2>&1";

        if (runCommand("$wrk/5-consensus-insert-sizes", $cmd)) {
            caFailure("Insert size estimation failed", "$wrk/5-consensus-insert-sizes/estimates.out");
        }

        unlink("$wrk/5-consensus-insert-sizes/$asm.tigStore");
    }

    #
    #  Update estimates in gatekeeper
    #

    if (! -e "$wrk/5-consensus-insert-sizes/updates.err") {
        if (! -e "$wrk/5-consensus-insert-sizes/$asm.tigStore.distupdate") {
            rename "$wrk/5-consensus-insert-sizes/estimates.out", "$wrk/5-consensus-insert-sizes/estimates.out.FAILED";
            caFailure("Failed to find insert size estimates", "$wrk/5-consensus-insert-sizes/estimates.out.FAILED");
        }

        $cmd  = "$bin/gatekeeper \\\n";
        $cmd .= " --edit $wrk/5-consensus-insert-sizes/$asm.tigStore.distupdate \\\n";
        $cmd .= "        $wrk/$asm.gkpStore \\\n";
        $cmd .= "> $wrk/5-consensus-insert-sizes/updates.err 2>&1";

        if (runCommand("$wrk/5-consensus-insert-sizes", $cmd)) {
            rename "$wrk/5-consensus-insert-sizes/updates.err", "$wrk/5-consensus-insert-sizes/updates.err.FAILED";
            caFailure("Insert size updates failed", "$wrk/5-consensus-insert-sizes/updates.err.FAILED");
        }
    }

    #
    #  Run the chimeric unitig splitter, and then the fixer again.
    #

    if (getGlobal("doUnitigSplitting")) {
        if (! -e "$wrk/5-consensus-split/splitUnitigs.out") {
            system("mkdir $wrk/5-consensus-split") if (! -d "$wrk/5-consensus-split");

            $cmd  = "$bin/splitUnitigs \\\n";
            $cmd .= " -g $wrk/$asm.gkpStore \\\n";
            $cmd .= " -t $wrk/$asm.tigStore 3 \\\n";
            $cmd .= "> $wrk/5-consensus-split/splitUnitigs.out 2>&1";

            if (runCommand("$wrk/5-consensus-split", $cmd)) {
                rename "$wrk/5-consensus-split/splitUnitigs.out", "$wrk/5-consensus-split/splitUnitigs.out.FAILED";
                caFailure("Unitig splitter failed", "$wrk/5-consensus-split/splitUnitigs.out.FAILED");
            }
        }

        if (! -e "$wrk/5-consensus-split/consensus-fix.out") {
            $cmd  = "$bin/utgcnsfix \\\n";
            $cmd .= " -g $wrk/$asm.gkpStore \\\n";
            $cmd .= " -t $wrk/$asm.tigStore 4 . \\\n";
            $cmd .= "> $wrk/5-consensus-split/consensus-fix.out 2>&1";

            if (runCommand("$wrk/5-consensus-split", $cmd)) {
                rename "$wrk/5-consensus-split/consensus-fix.out", "$wrk/5-consensus-split/consensus-fix.out.FAILED";
                caFailure("post split utgcnsfix failed", "$wrk/5-consensus-split/consensus-fix.out.FAILED");
            }
        }
    }

    #
    #  And finally, compute the coverage stat for all unitigs
    #
    my $l = getGlobal("utgGenomeSize");

    if (! -e "$wrk/5-consensus-coverage-stat/computeCoverageStat.err") {
        system("mkdir $wrk/5-consensus-coverage-stat") if (! -d "$wrk/5-consensus-coverage-stat");

        $cmd  = "$bin/computeCoverageStat \\\n";
        $cmd .= " -g $wrk/$asm.gkpStore \\\n";
        $cmd .= " -t $wrk/$asm.tigStore 5 \\\n";
        $cmd .= " -s $l \\\n" if defined($l);
        $cmd .= " -o $wrk/5-consensus-coverage-stat/$asm \\\n";
        $cmd .= "> $wrk/5-consensus-coverage-stat/computeCoverageStat.err 2>&1";

        if (runCommand("$wrk/5-consensus-coverage-stat", $cmd)) {
            rename "$wrk/5-consensus-coverage-stat/computeCoverageStat.err", "$wrk/5-consensus-coverage-stat/computeCoverageStat.err.FAILED";
            caFailure("Unitig coverage stat computation failed", "$wrk/5-consensus-coverage-stat/computeCoverageStat.err.FAILED");
        }

        if (-d "$wrk/7-0-CGW") {
            caFailure("Unitig coverage stat updated, but there is a scaffolding already started.  Remove old scaffold directories to proceed.", "");
        }
    }

    if (! -e "$wrk/5-consensus-coverage-stat/markRepeatUnique.err") {
        my $astatLow       = getGlobal("astatLowBound");
        my $astatHigh      = getGlobal("astatHighBound");

        my $maxSingleSpan  = getGlobal("maxSingleReadSpan");
        my $lowCovDepth    = getGlobal("lowCoverageDepth");
        my $lowCovAllowed  = getGlobal("lowCoverageAllowed");
        my $minReadsUnique = getGlobal("minReadsUnique");
        my $maxRepeatLen   = getGlobal("maxRepeatLength");

        system("mkdir $wrk/5-consensus-coverage-stat") if (! -d "$wrk/5-consensus-coverage-stat");

        $cmd  = "$bin/markRepeatUnique \\\n";
        $cmd .= " -g $wrk/$asm.gkpStore \\\n";
        $cmd .= " -t $wrk/$asm.tigStore 5 \\\n";
        $cmd .= " -j $astatLow \\\n";
        $cmd .= " -k $astatHigh \\\n";
        $cmd .= " -span   $maxSingleSpan \\\n"              if (defined($maxSingleSpan));
        $cmd .= " -lowcov $lowCovDepth $lowCovAllowed \\\n" if (defined($lowCovAllowed));
        $cmd .= " -reads  $minReadsUnique \\\n"             if (defined($minReadsUnique));
        $cmd .= " -length $maxRepeatLen \\\n"               if (defined($maxRepeatLen));
        $cmd .= " -o $wrk/5-consensus-coverage-stat/$asm.markRepeatUnique \\\n";
        $cmd .= "> $wrk/5-consensus-coverage-stat/markRepeatUnique.err 2>&1";

        if (runCommand("$wrk/5-consensus-coverage-stat", $cmd)) {
            rename "$wrk/5-consensus-coverage-stat/markRepeatUnique.err", "$wrk/5-consensus-coverage-stat/markRepeatUnique.err.FAILED";
            caFailure("Unitig repeat/unique marking failed", "$wrk/5-consensus-coverage-stat/markRepeatUnique.err.FAILED");
        }

        if (-d "$wrk/7-0-CGW") {
            caFailure("Unitig repeat/unique markings updated, but there is a scaffolding already started.  Remove old scaffold directories to proceed.", "");
        }
    }

    #  All jobs finished.  Remove the partitioning from the gatekeeper store.
    #
    #system("rm -f $wrk/$asm.gkpStore/???.[0-9][0-9][0-9]");

    touch("$wrk/5-consensus/consensus.success");

  alldone:
    stopAfter("utgcns");
    stopAfter("consensusAfterUnitigger");
}

################################################################################
################################################################################
################################################################################


#  Don't do interleaved merging unless we are throwing stones.

sub CGW ($$$$$$) {
    my $thisDir     = shift @_;
    my $lastDir     = shift @_;
    my $tigStore    = shift @_;
    my $stoneLevel  = shift @_;
    my $logickp     = shift @_;
    my $finalRun    = shift @_;
    my $lastckp     = undef;
    my $ckp         = undef;

    return($thisDir) if (-e "$wrk/$thisDir/cgw.success");

    if (defined($lastDir)) {
        $lastckp = findLastCheckpoint($lastDir);
    }
    if (defined($lastckp) && defined($logickp)) {
        $ckp     = "-R $lastckp -N $logickp"
    }

    #  If there is a timing file here, assume we are restarting.  Not
    #  all restarts are possible, but we try hard to make it so.
    #
    if (-e "$wrk/$thisDir/$asm.timing") {
        my $restartckp = undef;

        open(F, "< $wrk/$thisDir/$asm.timing");
        while (<F>) {
            if (m/Writing.*ckp.(\d+)\s\(logical\s(.+)\)/) {
                $restartckp = "-R $1 -N $2";
            }
        }
        close(F);

        if (!defined($restartckp)) {
            print STDERR "Found an empty timing file, starting from the beginning: $ckp\n";
        } else {
            $ckp = $restartckp;
            print STDERR "Found a timing file, restarting: $ckp\n";
        }
    }

    system("mkdir $wrk/$thisDir")               if (! -d "$wrk/$thisDir");

    system("ln -s ../$lastDir/$asm.ckp.$lastckp $wrk/$thisDir/$asm.ckp.$lastckp") if (defined($lastDir));

    if (-e "$wrk/$thisDir/cgw.out") {
        my $ckp = findLastCheckpoint($thisDir);
        my $ver = "00";
        while (-e "$wrk/$thisDir/cgw.out.$ver.ckp.$ckp") {
            $ver++;
        }
        rename "$wrk/$thisDir/cgw.out", "$wrk/$thisDir/cgw.out.$ver.ckp.$ckp"
    }

    my $sampleSize = getGlobal("cgwDistanceSampleSize");

    my $astatLow = getGlobal("astatLowBound");
    my $astatHigh = getGlobal("astatHighBound");

    my $B = int($numFrags / getGlobal("cnsPartitions"));
    $B = getGlobal("cnsMinFrags") if ($B < getGlobal("cnsMinFrags"));

    my $P = getGlobal("closurePlacement");

    my $shatterLevel  = getGlobal("cgwContigShatterWeight");
    my $missingMate   = getGlobal("cgwMergeMissingThreshold");
    my $minWeight     = getGlobal("cgwMinMergeWeight");
    my $filterLevel   = getGlobal("cgwMergeFilterLevel");

    $cmd  = "$bin/cgw $ckp \\\n";
    $cmd .= "  -j $astatLow -k $astatHigh \\\n";
    $cmd .= "  -r 5 \\\n";
    $cmd .= "  -s $stoneLevel \\\n";
    $cmd .= "  -filter $filterLevel \\\n"                if ($filterLevel != 0);
    $cmd .= "  -minmergeweight $minWeight \\\n";
    $cmd .= "  -S 0 \\\n"                                if (($finalRun == 0)   || (getGlobal("doResolveSurrogates") == 0));
    $cmd .= "  -G \\\n"                                  if ($finalRun == 0);
    $cmd .= "  -GG \\\n"                                 if (getGlobal("cgwPreserveConsensus") == 1);
    $cmd .= "  -z \\\n"                                  if (getGlobal("cgwDemoteRBP") == 1);
    $cmd .= "  -P $P \\\n"                               if (defined($P));
    $cmd .= "  -K \\\n"                                  if (getGlobal("kickOutNonOvlContigs") != 0);
    $cmd .= "  -U \\\n"                                  if (getGlobal("doUnjiggleWhenMerging") != 0);
    $cmd .= "  -F \\\n"                                  if (getGlobal("toggleDoNotDemote") != 0);
    $cmd .= "  -B $B \\\n";
    $cmd .= "  -u $wrk/4-unitigger/$asm.unused.ovl \\\n" if (getGlobal("cgwUseUnitigOverlaps") != 0);
    $cmd .= "  -reloadmates \\\n"                        if (getGlobal("cgwReloadMates") != 0);
    $cmd .= "  -shatter $shatterLevel \\\n";
    $cmd .= "  -missingMate $missingMate \\\n";
    $cmd .= "  -m $sampleSize \\\n";
    $cmd .= "  -g $wrk/$asm.gkpStore \\\n";
    $cmd .= "  -t $tigStore \\\n";
    $cmd .= "  -o $wrk/$thisDir/$asm \\\n";
    $cmd .= " > $wrk/$thisDir/cgw.out 2>&1";

    stopBefore("CGW", $cmd);

    if (runCommand("$wrk/$thisDir", $cmd)) {
        caFailure("scaffolder failed", "$wrk/$thisDir/cgw.out");
    }

    if (getGlobal("cgwCompressTigStore") != 0) {
        my $f = findFirstCheckpoint($thisDir);
        my $l = findLastCheckpoint($thisDir);

        if ($f < $l) {
            $cmd  = "$bin/tigStore \\\n";
            $cmd .= " -g $wrk/$asm.gkpStore \\\n";
            $cmd .= " -t $tigStore $l \\\n";
            $cmd .= " -compress \\\n";
            $cmd .= "> $wrk/$thisDir/tigStore-compress.out 2>&1";

            if (runCommand("$wrk/$thisDir", $cmd)) {
                caFailure("tigStore compression failed", "$wrk/$thisDir/tigStore-compress.out");
            }
        }
    }

    if (getGlobal("cgwPurgeCheckpoints") != 0) {
        my $f = findFirstCheckpoint($thisDir);
        my $l = findLastCheckpoint($thisDir);

        while ($f < $l) {
            #print STDERR "Purging $wrk/$thisDir/$asm.ckp.$f\n";
            unlink "$wrk/$thisDir/$asm.ckp.$f";
            $f++;
        }
    }

    touch("$wrk/$thisDir/cgw.success");

    return $thisDir;
}


sub eCR ($$$) {
    my $thisDir = shift @_;
    my $lastDir = shift @_;
    my $iter    = shift @_;

    return $thisDir if (-e "$wrk/$thisDir/extendClearRanges.success");

    my $lastckp = findLastCheckpoint($lastDir);

    system("mkdir $wrk/$thisDir") if (! -d "$wrk/$thisDir");

    system("ln -s ../$lastDir/$asm.ckp.$lastckp $wrk/$thisDir/$asm.ckp.$lastckp")  if (! -e "$wrk/$thisDir/$asm.ckp.$lastckp");

    #  Partition eCR.

    if (! -e "$wrk/$thisDir/extendClearRanges.partitionInfo") {
        $cmd  = "$bin/extendClearRangesPartition ";
        $cmd .= " -g $wrk/$asm.gkpStore ";
        $cmd .= " -t $wrk/$asm.tigStore ";
        $cmd .= " -n $lastckp ";
        $cmd .= " -c $asm ";
        $cmd .= " -N 4 ";
        $cmd .= " -p $wrk/$thisDir/extendClearRanges.partitionInfo";
        $cmd .= "  > $wrk/$thisDir/extendClearRanges.partitionInfo.err 2>&1";

        stopBefore("eCRPartition", $cmd);
        stopBefore("extendClearRangesPartition", $cmd);

        if (runCommand("$wrk/$thisDir", $cmd)) {
            caFailure("extendClearRanges partitioning failed", "$wrk/$thisDir/extendClearRanges.partitionInfo.err");
        }

        #  Remove any existing eCR scripts -- possibly left behind by the user deleting
        #  the partitioinInfo and restarting.

        open(F, "find -L $wrk/$thisDir -name 'extendClearRanges-scaffold.*' -print |");
        while (<F>) {
            chomp;
            print STDERR "Repartitioned; remove extraneous file $_\n";
            unlink $_;
        }
        close(F);
    }

    #  Read the partitioning info, create jobs.  No partitions?  No ECR jobs.

    my @jobs;

    open(P, "< $wrk/$thisDir/extendClearRanges.partitionInfo") or caFailure("failed to find extendClearRanges partitioning file $wrk/$thisDir/extendClearRanges.partitionInfo", undef);
    while (<P>) {
        #  Fields are: partitionNum BgnScf partitionNum EndScf NumFrags

        my @v = split '\s+', $_;

        my $curScaffold = substr("000000000$v[1]", -7);
        my $endScaffold = substr("000000000$v[3]", -7);

        my $j = "$wrk/$thisDir/extendClearRanges-scaffold.$curScaffold";

        if (! -e "$j.success") {
            if (! -e "$j.sh") {
                open(F, "> $j.sh");
                print F "#!" . getGlobal("shell") . "\n\n";
                print F "\n";
                print F "AS_OVL_ERROR_RATE=", getGlobal("ovlErrorRate"), "\n";
                print F "AS_CNS_ERROR_RATE=", getGlobal("cnsErrorRate"), "\n";
                print F "AS_CGW_ERROR_RATE=", getGlobal("cgwErrorRate"), "\n";
                print F "AS_OVERLAP_MIN_LEN=", getGlobal("ovlMinLen"),    "\n";
                print F "AS_READ_MIN_LEN="   , getGlobal("frgMinLen"),    "\n";
                print F "export AS_OVL_ERROR_RATE AS_CNS_ERROR_RATE AS_CGW_ERROR_RATE AS_OVERLAP_MIN_LEN AS_READ_MIN_LEN\n";
                print F "\n";
                print F "$bin/extendClearRanges \\\n";
                print F " -g $wrk/$asm.gkpStore \\\n";
                print F " -t $wrk/$asm.tigStore \\\n";
                print F " -n $lastckp \\\n";
                print F " -c $asm \\\n";
                print F " -b $curScaffold -e $endScaffold \\\n";
                print F " -i $iter \\\n";
                print F " > $j.err 2>&1\n";
                close(F);

                system("chmod +x $j.sh");
            }

            push @jobs, "$j";

            $lastckp++;
        }
    }
    close(P);

    #  Run jobs.

    stopBefore("eCR", undef);
    stopBefore("extendClearRanges", undef);

    foreach my $j (@jobs) {
        if (runCommand("$wrk/$thisDir", "$j.sh")) {
            caFailure("extendClearRanges failed", "$j.err");
        }
        touch("$j.success");
    }

    touch("$wrk/$thisDir/extendClearRanges.success");

    return $thisDir;
}


sub updateDistanceRecords ($) {
    my $thisDir = shift @_;

    return $thisDir   if (-e "$wrk/$thisDir/$asm.distupdate.success");

    open(F, "> $wrk/$thisDir/$asm.distupdate");

    print F "## Scaffold Based Estimates\n";
    if (-e "$wrk/$thisDir/stat/scaffold_final.distupdate") {
        open(A, "< $wrk/$thisDir/stat/scaffold_final.distupdate");
        while (<A>) {
            print F $_;
        }
        close(A);
    }

    print F "## Contig Based Estimates\n";
    if (-e "$wrk/$thisDir/stat/contig_final.distupdate") {
        open(A, "< $wrk/$thisDir/stat/contig_final.distupdate");
        while (<A>) {
            print F $_;
        }
        close(A);
    }

    close(F);

    $cmd  = "$bin/gatekeeper \\\n";
    $cmd .= " --edit $wrk/$thisDir/$asm.distupdate \\\n";
    $cmd .= "        $wrk/$asm.gkpStore \\\n";
    $cmd .= "> $wrk/$thisDir/$asm.distupdate.err 2>&1";

    if (runCommand("$wrk/$thisDir", $cmd)) {
        caFailure("Insert size updates failed", "$wrk/$thisDir/$asm.distupdate.err");
    }

    touch("$wrk/$thisDir/$asm.distupdate.success");

    return $thisDir;
}


sub scaffolder () {
    my $lastDir    = undef;
    my $thisDir    = 0;
    my $stoneLevel = getGlobal("stoneLevel");

    stopBefore("scaffolder", undef);

    goto alldone if (-e "$wrk/7-CGW/cgw.success");

    #  Do an initial CGW to update distances, then update the
    #  gatekeeper.  This initial run shouldn't be used for later
    #  CGW'ing.
    #
    my $cis = getGlobal("computeInsertSize");
    if (!defined($cis) && ($numFrags < 1000000)) {
        $cis = 1;
    }
    if ($cis == 1) {
        if (! -e "$wrk/6-clonesize/$asm.tigStore") {
            system("mkdir -p $wrk/6-clonesize/$asm.tigStore");
            system("cd $wrk/6-clonesize/$asm.tigStore && ln -s ../../$asm.tigStore/* .");
        }
        updateDistanceRecords(CGW("6-clonesize", undef, "$wrk/6-clonesize/$asm.tigStore", $stoneLevel, undef, 0));
    }


    #  If we're not doing eCR, we just do a single scaffolder run, and
    #  get the heck outta here!
    #
    if (getGlobal("doExtendClearRanges") == 0) {
        $lastDir = CGW("7-$thisDir-CGW", $lastDir, "$wrk/$asm.tigStore", $stoneLevel, undef, 1);
        $thisDir++;
    } else {

        #  Do the initial CGW, making sure to not throw stones.
        #
        $lastDir = updateDistanceRecords(CGW("7-$thisDir-CGW", $lastDir, "$wrk/$asm.tigStore", 0, undef, 0));
        $thisDir++;

        #  Followed by at least one eCR
        #
        $lastDir = eCR("7-$thisDir-ECR", $lastDir, 1);
        $thisDir++;

        #  Iterate eCR: do another scaffolder still without stones,
        #  then another eCR.  Again, and again, until we get dizzy and
        #  fall over.
        #
        my $iterationMax = getGlobal("doExtendClearRanges") + 1;
        for (my $iteration = 2; $iteration < $iterationMax; $iteration++) {
            $lastDir = updateDistanceRecords(CGW("7-$thisDir-CGW", $lastDir, "$wrk/$asm.tigStore", 0, "ckp03-SCF-partial", 0));
            $thisDir++;

            $lastDir = eCR("7-$thisDir-ECR", $lastDir, $iteration);
            $thisDir++;
        }

        #  Then another scaffolder, chucking stones into the big holes,
        #  filling in surrogates, and writing output.
        #
        $lastDir = updateDistanceRecords(CGW("7-$thisDir-CGW", $lastDir, "$wrk/$asm.tigStore", $stoneLevel, "ckp03-SCF-partial", 1));
        $thisDir++;
    }


    #  And, finally, hold on, we're All Done!  Point to the correct output directory.
    #
    system("ln -s $lastDir $wrk/7-CGW") if (! -d "$wrk/7-CGW");

  alldone:
    stopAfter("scaffolder");
}


################################################################################
################################################################################
################################################################################


#  Prepare for consensus on the grid
#    Partition the contigs
#    Repartition the frag store

sub createPostScaffolderConsensusJobs () {
    my $blasr = findBlasr(getGlobal("consensus"));
    if (!defined($blasr)) { setGlobal("consensus", "cns"); }
    my $consensusType = getGlobal("consensus");

    return if (-e "$wrk/8-consensus/consensus.sh");

    caFailure("contig consensus didn't find '$wrk/$asm.tigStore'", undef)    if (! -d "$wrk/$asm.tigStore");

    my $tigVersion = findLastCheckpoint("$wrk/7-CGW");
    caFailure("contig consensus didn't find any checkpoints in '$wrk/7-CGW'", undef) if (!defined($tigVersion));

    ########################################
    #
    #  Partition the gkpStore for consensus.
    #
    if (! -e "$wrk/8-consensus/$asm.partitioned") {
        $cmd  = "$bin/gatekeeper -P $wrk/7-CGW/$asm.partitioning $wrk/$asm.gkpStore ";
        $cmd .= "> $wrk/8-consensus/$asm.partitioned.err 2>&1";

        caFailure("gatekeeper partitioning failed", "$wrk/8-consensus/$asm.partitioned.err") if (runCommand("$wrk/8-consensus", $cmd));
        touch("$wrk/8-consensus/$asm.partitioned");
    }

    ########################################
    #
    #  Build consensus jobs for the grid -- this is very similar to that in createPostUnitiggerConsensus.pl
    #
    my $jobs = 0;

    open(F, "< $wrk/7-CGW/$asm.partitionInfo") or caFailure("can't open '$wrk/7-CGW/$asm.partitionInfo'", undef);
    while (<F>) {
        if (m/Partition\s+(\d+)\s+has\s+(\d+)\s+contigs\sand\s+(\d+)\s+fragments./) {
            $jobs = $1;
        }
    }
    close(F);

    my $taskID       = getGlobal("gridTaskID");
    open(F, "> $wrk/8-consensus/consensus.sh") or caFailure("can't open '$wrk/8-consensus/consensus.sh'", undef);
    print F "#!" . getGlobal("shell") . "\n";
    print F "\n";
    print F "jobid=\$$taskID\n";
    print F "if [ x\$jobid = x -o x\$jobid = xundefined -o x\$jobid = x0 ]; then\n";
    print F "  jobid=\$1\n";
    print F "fi\n";
    print F "if [ x\$jobid = x ]; then\n";
    print F "  echo Error: I need $taskID set, or a job index on the command line.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "if [ \$jobid -gt $jobs ]; then\n";
    print F "  echo Error: Only $jobs partitions, you asked for \$jobid.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "\n";
    print F "jobid=`printf %03d \$jobid`\n";
    print F "\n";
    print F "if [ -e $wrk/8-consensus/${asm}_\$jobid.success ] ; then\n";
    print F "  exit 0\n";
    print F "fi\n";
    print F "\n";
    print F "AS_OVL_ERROR_RATE=", getGlobal("ovlErrorRate"), "\n";
    print F "AS_CNS_ERROR_RATE=", getGlobal("cnsErrorRate"), "\n";
    print F "AS_CGW_ERROR_RATE=", getGlobal("cgwErrorRate"), "\n";
    print F "AS_OVERLAP_MIN_LEN=", getGlobal("ovlMinLen"),    "\n";
    print F "AS_READ_MIN_LEN="   , getGlobal("frgMinLen"),    "\n";
    print F "export AS_OVL_ERROR_RATE AS_CNS_ERROR_RATE AS_CGW_ERROR_RATE AS_OVERLAP_MIN_LEN AS_READ_MIN_LEN\n";

    print F getBinDirectoryShellCode();

    if ($consensusType eq "cns") {
        print F "\$bin/ctgcns \\\n";
        print F "  -g $wrk/$asm.gkpStore \\\n";
        print F "  -t $wrk/$asm.tigStore $tigVersion \$jobid \\\n";
        print F "  -P ", getGlobal("cnsPhasing"), " \\\n";
        print F "  -U \\\n"  if (getGlobal("cnsReuseUnitigs") != 0);
        print F " > $wrk/8-consensus/${asm}_\$jobid.err 2>&1 \\\n";
        print F "&& \\\n";
        print F "touch $wrk/8-consensus/${asm}_\$jobid.success\n";
    } elsif ($consensusType eq "seqan") {
        print F "\$bin/SeqAn_CNS \\\n";
        print F "  -G $wrk/$asm.gkpStore \\\n";
        print F "  -u $wrk/$asm.SeqStore \\\n";
        print F "  -V $tigVersion \\\n";
        print F "  -p \$jobid \\\n";
        print F "  -S \$jobid \\\n";
        #print F "  -c $cgwDir/$asm.cgw_contigs.\$jobid \\\n";
        print F "  -s \$bin/graph_consensus \\\n";
        print F "  -w $wrk/8-consensus/ \\\n";
        print F "  -o $wrk/8-consensus/$asm.cns_contigs.\$jobid \\\n";
        print F " > $wrk/8-consensus/$asm.cns_contigs.\$jobid.err 2>&1 \\\n";
        print F "&& \\\n";
        print F "touch $wrk/8-consensus/$asm.cns_contigs.\$jobid.success\n";
    } elsif ($consensusType eq "pbdagcon" || $consensusType eq "pbutgcns") {
        print F "\$bin/tigStore -d layout -C -t $wrk/$asm.tigStore $tigVersion -cp \$jobid -g $wrk/$asm.gkpStore > $wrk/8-consensus/$asm.\$jobid.lay\n";
        print F "\$bin/tigStore -d consensus -U -t $wrk/$asm.tigStore 2 -cp \$jobid -g $wrk/$asm.gkpStore > $wrk/8-consensus/$asm.\$jobid.fasta\n";
        print F "\$bin/convertToPBCNS -path $blasr -consensus $consensusType -coverage 1 -threads " . getGlobal("cnsConcurrency") . " -prefix $wrk/8-consensus/$asm.\$jobid.tmp -length 500 -sequence $wrk/8-consensus/$asm.\$jobid.fasta -input $wrk/8-consensus/$asm.\$jobid.lay -output $wrk/8-consensus/$asm.\$jobid.fa\n";
        print F "\$bin/addCNSToStore -path \$bin -version $tigVersion -input $wrk/8-consensus/$asm.\$jobid.fa -lay $wrk/8-consensus/$asm.\$jobid.lay -output $wrk/8-consensus/$asm.\$jobid.cns -prefix $wrk/$asm -sequence $wrk/8-consensus/$asm.\$jobid.fasta -partition \$jobid && touch $wrk/8-consensus/${asm}_\$jobid.success\n";
        print F "if [ -e $wrk/8-consensus/${asm}_\$jobid.success ]; then\n";
        print F "   rm -f $wrk/8-consensus/${asm}.\$jobid.fasta*\n";
        print F "   rm -f $wrk/8-consensus/${asm}.\$jobid.lay\n";
        print F "fi\n";
        setGlobal("cnsConcurrency", 1);
    } else {
        caFailure("unknown consensus type $consensusType; must be 'cns' or 'seqan'", undef);
    }
    print F "exit 0\n";
    close(F);

    chmod 0755, "$wrk/8-consensus/consensus.sh";

    if (getGlobal("cnsOnGrid") && getGlobal("useGrid")) {
        my $submitCommand  = getGlobal("gridSubmitCommand");
        my $nameOption   = getGlobal("gridNameOption");
        my $outputOption  = getGlobal("gridOutputOption");

        my $sge          = getGlobal("sge");
        my $sgeName      = getGlobal("sgeName");
        my $sgeConsensus = getGlobal("sgeConsensus");

        $sgeName = "_$sgeName" if (defined($sgeName));
        my $jobName = getGridArrayName("ctg_$asm$sgeName", $jobs);
        my $arrayOpt = getGridArrayOption("ctg_$asm$sgeName", $jobs);

        my $SGE;
        $SGE  = "$submitCommand $sge $sgeConsensus $nameOption \"$jobName\" $arrayOpt ";
        $SGE .= "$outputOption /dev/null ";
        $SGE .= "$wrk/8-consensus/consensus.sh\n";

        submitBatchJobs($SGE, $jobName);
        exit(0);
    } else {
        for (my $i=1; $i<=$jobs; $i++) {
            schedulerSubmit("$wrk/8-consensus/consensus.sh $i > /dev/null 2>&1");
        }

        schedulerSetNumberOfProcesses(getGlobal("cnsConcurrency"));
        schedulerFinish();
    }
}


sub postScaffolderConsensus () {

    system("mkdir $wrk/8-consensus") if (! -d "$wrk/8-consensus");

    goto alldone if (-e "$wrk/8-consensus/consensus.success");

    createPostScaffolderConsensusJobs();

    #
    #  Check that consensus finished properly
    #
    my $failedJobs = 0;

    open(F, "< $wrk/7-CGW/$asm.partitionInfo") or caFailure("can't open '$wrk/7-CGW/$asm.partitionInfo'", undef);
    while (<F>) {
        if (m/Partition\s+(\d+)\s+has\s+(\d+)\s+contigs\sand\s+(\d+)\s+fragments./) {
            my $id = substr("000" . $1, -3);

            if (! -e "$wrk/8-consensus/${asm}_$id.success") {
                print STDERR "$wrk/8-consensus/${asm}_$id failed -- no .success.\n";
                $failedJobs++;
            }
        }
    }
    close(F);

    #  FAILUREHELPME
    #
    caFailure("$failedJobs consensusAfterScaffolder jobs failed; remove $wrk/8-consensus/consensus.sh to try again", undef) if ($failedJobs);

    #  All jobs finished.  Remove the partitioning from the gatekeeper store.
    #
    #system("rm -f $wrk/$asm.gkpStore/???.[0-9][0-9][0-9]");

    touch("$wrk/8-consensus/consensus.success");

  alldone:
    stopAfter("ctgcns");
    stopAfter("consensusAfterScaffolder");
}

################################################################################
################################################################################
################################################################################


sub summarizeConsensusStatistics ($) {
    my $dir = shift @_;

    if (! -e "$dir/consensus.stats.summary") {
        my $NumColumnsInUnitigs           = 0;
        my $NumGapsInUnitigs              = 0;
        my $NumRunsOfGapsInUnitigReads    = 0;
        my $NumColumnsInContigs           = 0;
        my $NumGapsInContigs              = 0;
        my $NumRunsOfGapsInContigReads    = 0;
        my $NumAAMismatches               = 0;
        my $NumFAMismatches               = 0;
        my $NumVARRecords                 = 0;
        my $NumVARStringsWithFlankingGaps = 0;
        my $NumUnitigRetrySuccess         = 0;

        open(F, "ls $dir/$asm*.err |");
        my @files = <F>;
        chomp @files;
        close(F);

        foreach my $f (@files) {
            open(F, "< $f");
            while (<F>) {
                $NumColumnsInUnitigs += $1           if (m/NumColumnsInUnitigs\s+=\s+(\d+)/);
                $NumGapsInUnitigs += $1              if (m/NumGapsInUnitigs\s+=\s+(\d+)/);
                $NumRunsOfGapsInUnitigReads += $1    if (m/NumRunsOfGapsInUnitigReads\s+=\s+(\d+)/);
                $NumColumnsInContigs += $1           if (m/NumColumnsInContigs\s+=\s+(\d+)/);
                $NumGapsInContigs += $1              if (m/NumGapsInContigs\s+=\s+(\d+)/);
                $NumRunsOfGapsInContigReads += $1    if (m/NumRunsOfGapsInContigReads\s+=\s+(\d+)/);
                $NumAAMismatches += $1               if (m/NumAAMismatches\s+=\s+(\d+)/);
                $NumFAMismatches += $1               if (m/NumFAMismatches\s+=\s+(\d+)/);
                $NumVARRecords += $1                 if (m/NumVARRecords\s+=\s+(\d+)/);
                $NumVARStringsWithFlankingGaps += $1 if (m/NumVARStringsWithFlankingGaps\s+=\s+(\d+)/);
                $NumUnitigRetrySuccess += $1         if (m/NumUnitigRetrySuccess\s+=\s+(\d+)/);
            }
            close(F);
        }

        open(F, "> $dir/consensus.stats.summary");
        print F "NumColumnsInUnitigs=$NumColumnsInUnitigs\n"                     if ($NumColumnsInUnitigs > 0);
        print F "NumGapsInUnitigs=$NumGapsInUnitigs\n"                           if ($NumGapsInUnitigs > 0);
        print F "NumRunsOfGapsInUnitigReads=$NumRunsOfGapsInUnitigReads\n"       if ($NumRunsOfGapsInUnitigReads > 0);
        print F "NumColumnsInContigs=$NumColumnsInContigs\n"                     if ($NumColumnsInContigs > 0);
        print F "NumGapsInContigs=$NumGapsInContigs\n"                           if ($NumGapsInContigs > 0);
        print F "NumRunsOfGapsInContigReads=$NumRunsOfGapsInContigReads\n"       if ($NumRunsOfGapsInContigReads > 0);
        print F "NumAAMismatches=$NumAAMismatches\n"                             if ($NumAAMismatches > 0);
        print F "NumFAMismatches=$NumFAMismatches\n"                             if ($NumFAMismatches > 0);
        print F "NumVARRecords=$NumVARRecords\n"                                 if ($NumVARRecords > 0);
        print F "NumVARStringsWithFlankingGaps=$NumVARStringsWithFlankingGaps\n" if ($NumVARStringsWithFlankingGaps > 0);
        print F "NumUnitigRetrySuccess=$NumUnitigRetrySuccess\n"                 if ($NumUnitigRetrySuccess > 0);
        close(F);
    }
}



sub terminate () {
    my $perl = "/usr/bin/env perl";

    my $termDir = "$wrk/9-terminator";
    system("mkdir $termDir") if (! -e "$termDir");

    stopBefore("terminator", undef);

    if (! -e "$termDir/$asm.asm") {
        my $uidServer = getGlobal("uidServer");
        my $fakeUIDs  = getGlobal("fakeUIDs");

        my $ckpVersion = findLastCheckpoint("$wrk/7-CGW");
        my $tigVersion = $ckpVersion + 1;

        caFailure("contig consensus didn't find any checkpoints in '$wrk/7-CGW'", undef) if (!defined($tigVersion));

        $cmd  = "$bin/terminator";
        $cmd .= " -E $uidServer" if (defined($uidServer));
        $cmd .= " -s $fakeUIDs" if ($fakeUIDs > 0);
        $cmd .= " -g $wrk/$asm.gkpStore";
        $cmd .= " -t $wrk/$asm.tigStore $tigVersion";
        $cmd .= " -c $wrk/7-CGW/$asm $ckpVersion";
        $cmd .= " -o $termDir/$asm";
        $cmd .= " > $termDir/$asm.asm.err 2>&1";

        if (runCommand("$termDir", $cmd)) {
            rename "$termDir/$asm.asm", "$termDir/$asm.asm.FAILED";
            rename "$termDir/$asm.map", "$termDir/$asm.map.FAILED";
            caFailure("terminator failed", "$termDir/$asm.asm.err");
        }
    }


    my $asmOutputFasta = "$bin/asmOutputFasta";
    if (! -e "$termDir/$asm.scf.fasta") {
        $cmd  = "$asmOutputFasta -p $termDir/$asm $termDir/$asm.asm > $termDir/asmOutputFasta.err 2>&1";
        if (runCommand("$termDir", $cmd)) {
            rename "$termDir/$asm.scfcns.fasta", "$termDir/$asm.scfcns.fasta.FAILED";
            caFailure("fasta output failed", "$termDir/asmOutputFasta.err");
        }
        unlink "$termDir/asmOutputFasta.err";
    }


    if (! -e "$termDir/$asm.singleton.fasta") {
        my $ckpVersion = findLastCheckpoint("$wrk/7-CGW");
        my $tigVersion = $ckpVersion + 1;

        $cmd  = "$bin/dumpSingletons ";
        $cmd .= " -g $wrk/$asm.gkpStore ";
        $cmd .= " -t $wrk/$asm.tigStore ";
        $cmd .= " -c $wrk/7-CGW/$asm -n $ckpVersion -S ";
        $cmd .= "> $termDir/$asm.singleton.fasta ";
        $cmd .= "2> $termDir/dumpSingletons.err ";
        if (runCommand("$termDir", $cmd)) {
            print STDERR "Failed.\n";
            rename "$termDir/$asm.singleton.fasta", "$termDir/$asm.singleton.fasta.FAILED";
        }
        unlink "$termDir/dumpSingletons.err";
    }


    ########################################
    #
    #  Generate fragment/unitig/contig/scaffold mappings
    #
    ########################################


    if (getGlobal("createPosMap") > 0) {
        if (! -e "$termDir/$asm.posmap.frgscf") {
            if (runCommand("$termDir", "$bin/buildPosMap -o $asm -g $wrk/$asm.gkpStore < $termDir/$asm.asm > $termDir/buildPosMap.err 2>&1")) {
                rename "$termDir/$asm.posmap.frgscf", "$termDir/$asm.posmap.frgscf.FAILED";
                caFailure("buildPosMap failed", "$termDir/buildPosMap.err");
            }
            unlink "$termDir/buildPosMap.err";
        }
    }

    ########################################
    #
    #  Generate a read depth histogram
    #
    ########################################
    if ((  -e "$termDir/$asm.posmap.frgscf") &&
        (! -e "$termDir/$asm.qc.readdepth") &&
        (! -e "$termDir/$asm.qc")) {

        #  Youch.  Run five commands, do something if all are successful.

        $cmd  = "sort -k2n -k3n -T $termDir $termDir/$asm.posmap.frgscf > $termDir/$asm.posmap.frgscf.sorted &&";
        $cmd .= "$bin/fragmentDepth -min       0 -max    3000 < $termDir/$asm.posmap.frgscf.sorted > $termDir/$asm.posmap.frgscf.histogram1 && ";
        $cmd .= "$bin/fragmentDepth -min    3001 -max   10000 < $termDir/$asm.posmap.frgscf.sorted > $termDir/$asm.posmap.frgscf.histogram2 && ";
        $cmd .= "$bin/fragmentDepth -min   10001 -max 1000000 < $termDir/$asm.posmap.frgscf.sorted > $termDir/$asm.posmap.frgscf.histogram3 && ";
        $cmd .= "$bin/fragmentDepth -min 1000001              < $termDir/$asm.posmap.frgscf.sorted > $termDir/$asm.posmap.frgscf.histogram4 ";

        if (runCommand("$termDir", $cmd) == 0) {
            my @H1;
            my @H2;
            my @H3;
            my @H4;
            my $histMax = 0;

            open(G, "<  $termDir/$asm.posmap.frgscf.histogram1") or caFailure("failed to open '$termDir/$asm.posmap.frgscf.histogram1'", undef);
            while (<G>) {
                my ($v, $s) = split '\s+', $_;
                $H1[$v] = $s;
                $histMax = $v if ($histMax < $v);
            }
            close(G);

            open(G, "<  $termDir/$asm.posmap.frgscf.histogram2") or caFailure("failed to open '$termDir/$asm.posmap.frgscf.histogram2'", undef);
            while (<G>) {
                my ($v, $s) = split '\s+', $_;
                $H2[$v] = $s;
                $histMax = $v if ($histMax < $v);
            }
            close(G);

            open(G, "<  $termDir/$asm.posmap.frgscf.histogram3") or caFailure("failed to open '$termDir/$asm.posmap.frgscf.histogram3'", undef);
            while (<G>) {
                my ($v, $s) = split '\s+', $_;
                $H3[$v] = $s;
                $histMax = $v if ($histMax < $v);
            }
            close(G);

            open(G, "<  $termDir/$asm.posmap.frgscf.histogram4") or caFailure("failed to open '$termDir/$asm.posmap.frgscf.histogram4'", undef);
            while (<G>) {
                my ($v, $s) = split '\s+', $_;
                $H4[$v] = $s;
                $histMax = $v if ($histMax < $v);
            }
            close(G);

            open(G, "> $termDir/$asm.qc.readdepth");
            print G "\n[Read Depth Histogram]\n";
            print G "d    < 3Kbp    < 10Kbp   < 1Mbp    < inf\n";
            for (my $v=0; $v<=$histMax; $v++) {
                printf(G "%-4d %-10d %-10d %-10d %-10d\n", $v, int($H1[$v]), int($H2[$v]), int($H3[$v]), int($H4[$v]));
            }
        }

        #  Remove our temporary files.

        unlink "$termDir/$asm.posmap.frgscf.histogram1";
        unlink "$termDir/$asm.posmap.frgscf.histogram2";
        unlink "$termDir/$asm.posmap.frgscf.histogram3";
        unlink "$termDir/$asm.posmap.frgscf.histogram4";
    }


    ########################################
    #
    #  Generate statistics.
    #
    ########################################

    if (! -e "$termDir/$asm.qc") {
        my $qcOptions;

        if ( -e "$wrk/$asm.frg" ) {
            link "$wrk/$asm.frg", "$termDir/$asm.frg";
            $qcOptions = "-metrics";
        }
        if ( -e "$wrk/$asm.catmap" && !-e "$termDir/$asm.catmap" )  {
            link "$wrk/$asm.catmap", "$termDir/$asm.catmap";
        }
        if ( -e "$wrk/$asm.seq.features" && !-e "$termDir/$asm.seq.features" )  {
            link "$wrk/$asm.seq.features", "$termDir/$asm.seq.features";
        }
        if (runCommand("$termDir", "$perl $bin/caqc.pl -euid $qcOptions $termDir/$asm.asm")) {
            rename "$termDir/$asm.qc", "$termDir/$asm.qc.FAILED";
        }

        summarizeConsensusStatistics("$wrk/5-consensus");
        summarizeConsensusStatistics("$wrk/8-consensus");

        open(F, ">> $termDir/$asm.qc") or caFailure("failed to append to '$termDir/$asm.qc'", undef);

        if (-e "$wrk/5-consensus/consensus.stats.summary") {
            print F "\n[Unitig Consensus]\n";
            open(G, "<  $wrk/5-consensus/consensus.stats.summary") or caFailure("failed to open '$wrk/5-consensus/consensus.stats.summary'", undef);
            while (<G>) {
                print F $_;
            }
            close(G);
        }

        if (-e "$wrk/8-consensus/consensus.stats.summary") {
            print F "\n[Contig Consensus]\n";
            open(G, "<  $wrk/8-consensus/consensus.stats.summary") or caFailure("failed to open '$wrk/8-consensus/consensus.stats.summary'", undef);
            while (<G>) {
                print F $_;
            }
            close(G);
        }

        if (-e "$termDir/$asm.qc.readdepth") {
            open(G, "< $termDir/$asm.qc.readdepth") or caFailure("failed to open '$termDir/$asm.qc.readdepth'", undef);
            while (<G>) {
                print F $_;
            }
            close(G);
        }

        close(F);

        unlink "$wrk/5-consensus/consensus.stats.summary";
        unlink "$wrk/8-consensus/consensus.stats.summary";
        unlink "$termDir/$asm.qc.readdepth";
    }


    ########################################
    #
    #  Mercy merQC
    #
    ########################################


    if ((getGlobal("merQC") > 0) &&
        (! -e "$termDir/$asm.merQC") &&
        (merylVersion() eq "Mighty")) {

        system("mkdir $termDir/mercy") if (! -e "$termDir/mercy");

        my $ms      = getGlobal("merQCmerSize");
        my $mem     = getGlobal("merQCmemory");
        my $verbose = "";

        if (! -e "$termDir/mercy/$asm-ms$ms-frgFull.mcidx") {
            $cmd  = "$bin/meryl -B -C -m $ms -threads 4 -memory $mem $verbose ";
            $cmd .= "-s $wrk/$asm.gkpStore:untrim ";
            $cmd .= "-o $termDir/mercy/$asm-ms$ms-frgFull";
            if (runCommand("$termDir/mercy", $cmd)) {
                print STDERR "Failed.\n";
                unlink "$termDir/mercy/$asm-ms$ms-frgFull.mcidx";
                unlink "$termDir/mercy/$asm-ms$ms-frgFull.mcdat";
            }
        }
        if (! -e "$termDir/mercy/$asm-ms$ms-frgTrim.mcidx") {
            $cmd  = "$bin/meryl -B -C -m $ms -threads 4 -memory $mem $verbose ";
            $cmd .= "-s $wrk/$asm.gkpStore ";
            $cmd .= "-o $termDir/mercy/$asm-ms$ms-frgTrim";
            if (runCommand("$termDir/mercy", $cmd)) {
                print STDERR "Failed.\n";
                unlink "$termDir/mercy/$asm-ms$ms-frgTrim.mcidx";
                unlink "$termDir/mercy/$asm-ms$ms-frgTrim.mcdat";
            }
        }

        #  XXX This can likely be optimized -- by feeding
        #  asmOutputcontigsFasta directly to meryl.  It'd be harder
        #  (but great) if only one pass through the asm file could be
        #  made.  Easier then if we write all three files at the same
        #  time.

        if (! -e "$termDir/mercy/$asm.ctgNorm.fasta") {
            link "$termDir/$asm.ctg.fasta", "$termDir/mercy/$asm.ctgNorm.fasta";
        }
        if (! -e "$termDir/mercy/$asm.ctgDreg.fasta") {
            link "$termDir/$asm.deg.fasta", "$termDir/mercy/$asm.ctgDreg.fasta";
        }
        if (! -e "$termDir/mercy/$asm.ctgAll.fasta") {
            system "cat $termDir/$asm.{ctg,deg}.fasta > $termDir/mercy/$asm.ctgAll.fasta";
        }

        if ((! -e "$termDir/mercy/$asm-ms$ms-ctgNorm.mcidx") &&
            (-e "$termDir/mercy/$asm.ctgNorm.fasta")) {
            $cmd  = "$bin/meryl -B -C -m $ms -threads 4 -segments 4 $verbose ";
            $cmd .= "-s $termDir/mercy/$asm.ctgNorm.fasta ";
            $cmd .= "-o $termDir/mercy/$asm-ms$ms-ctgNorm";
            if (runCommand("$termDir/mercy", $cmd)) {
                print STDERR "Failed.\n";
                unlink "$termDir/mercy/$asm-ms$ms-ctgNorm.mcidx";
                unlink "$termDir/mercy/$asm-ms$ms-ctgNorm.mcdat";
            }
        }
        if ((! -e "$termDir/mercy/$asm-ms$ms-ctgDreg.mcidx") &&
            (-e "$termDir/mercy/$asm.ctgDreg.fasta")) {
            $cmd  = "$bin/meryl -B -C -m $ms -threads 4 -segments 4 $verbose ";
            $cmd .= "-s $termDir/mercy/$asm.ctgDreg.fasta ";
            $cmd .= "-o $termDir/mercy/$asm-ms$ms-ctgDreg";
            if (runCommand("$termDir/mercy", $cmd)) {
                print STDERR "Failed.\n";
                unlink "$termDir/mercy/$asm-ms$ms-ctgDreg.mcidx";
                unlink "$termDir/mercy/$asm-ms$ms-ctgDreg.mcdat";
            }
        }
        if ((! -e "$termDir/mercy/$asm-ms$ms-ctgAll.mcidx") &&
            (-e "$termDir/mercy/$asm.ctgAll.fasta")) {
            $cmd  = "$bin/meryl -B -C -m $ms -threads 4 -segments 4 $verbose ";
            $cmd .= "-s $termDir/mercy/$asm.ctgAll.fasta ";
            $cmd .= "-o $termDir/mercy/$asm-ms$ms-ctgAll";
            if (runCommand("$termDir/mercy", $cmd)) {
                print STDERR "Failed.\n";
                unlink "$termDir/mercy/$asm-ms$ms-ctgAll.mcidx";
                unlink "$termDir/mercy/$asm-ms$ms-ctgAll.mcdat";
            }
        }

        if (! -e "$termDir/$asm-ms$ms.merQC") {
            $cmd  = "$bin/mercy ";
            $cmd .= "-af $termDir/mercy/$asm-ms$ms-frgFull "  if (-e "$termDir/mercy/$asm-ms$ms-frgFull.mcidx");
            $cmd .= "-tf $termDir/mercy/$asm-ms$ms-frgTrim "  if (-e "$termDir/mercy/$asm-ms$ms-frgTrim.mcidx");
            $cmd .= "-co $termDir/mercy/$asm-ms$ms-ctgNorm "  if (-e "$termDir/mercy/$asm-ms$ms-ctgNorm.mcidx");
            $cmd .= "-dc $termDir/mercy/$asm-ms$ms-ctgDreg "  if (-e "$termDir/mercy/$asm-ms$ms-ctgDreg.mcidx");
            $cmd .= "-ac $termDir/mercy/$asm-ms$ms-ctgAll "   if (-e "$termDir/mercy/$asm-ms$ms-ctgAll.mcidx");
            $cmd .= "> $termDir/$asm-ms$ms.merQC";
            if (runCommand("$termDir/mercy", $cmd)) {
                print STDERR "Failed.\n";
                rename "$termDir/$asm-ms$ms.merQC", "$termDir/$asm-ms$ms.merQC.FAILED";
            }
        }
    }


    ########################################
    #
    #  AGP and ACE file generation
    #
    ########################################


    if (getGlobal("createAGP") > 0) {
        if (! -e "$termDir/$asm.agp") {
            if (runCommand($termDir, "$perl $bin/asmToAGP.pl < $termDir/$asm.asm > $termDir/$asm.agp")) {
                rename "$termDir/$asm.agp", "$termDir/$asm.agp.FAILED";
            }
        }
    }

    if (getGlobal("createACE") > 0) {
        if (! -e "$termDir/$asm.ace.bz2") {
            if (! -e "$termDir/$asm.frg") {
                if (runCommand($termDir, "$bin/gatekeeper -dumpfrg -allreads $wrk/$asm.gkpStore > $termDir/$asm.frg 2> $termDir/gatekeeper.err")) {
                    caFailure("gatekeeper failed to dump fragments for ACE generation", "$termDir/gatekeeper.err");
                }
                unlink "$termDir/gatekeeper.err";
            }
            if (runCommand($termDir, "$perl $bin/ca2ace.pl $termDir/$asm.asm")) {
                rename "$termDir/$asm.ace.bz2", "$termDir/$asm.ace.FAILED.bz2";
            }
        }
    }

    ########################################
    #
    #  FASTQ dump
    #
    ########################################

    dumpInfo();

    if ((getGlobal("dumpFASTQ") > 0) &&
        (-e "$wrk/$asm.gkpStore.info")) {

        open(I, "< $wrk/$asm.gkpStore.info");
        while (<I>) {
            my ($libIID, $bgnIID, $endIID, $active, $deleted, $mated, $totLen, $clrLen, $libName) = split '\s+', $_;

            if ($libIID > 0) {
                my $prefix = sprintf "$wrk/9-terminator/$asm.lib.%03d.$libName", $libIID;

                if (! -e "$prefix.unmated.fastq") {
                    $cmd  = "$bin/gatekeeper ";
                    $cmd .= " -randomsubset $libIID 1.0 ";
                    $cmd .= " -dumpfastq $prefix ";
                    $cmd .= " $wrk/$asm.gkpStore ";
                    $cmd .= "> $prefix.err 2>&1";

                    runCommand($termDir, $cmd);
                }
            }
        }
        close(I);
    }


    ########################################
    #
    #  Link into the work directory
    #
    ########################################

    unlink "$wrk/$asm.asm";
    unlink "$wrk/$asm.qc";

    link "$termDir/$asm.asm", "$wrk/$asm.asm";
    link "$termDir/$asm.qc",  "$wrk/$asm.qc";

    return(0);
}

################################################################################
################################################################################
################################################################################


#  Assembly all done, toggle the unitigs and re-run CGW and subsequent steps of the assembly.

sub toggler () {
    my $toggledDir = "10-toggledAsm";
    my $ecrEdits = "";

    return if (-d "$wrk/$toggledDir/$asm.asm");
    return if (getGlobal("doToggle") == 0);

    my $minLength = getGlobal("toggleUnitigLength");
    my $numInstances = getGlobal("toggleNumInstances");
    my $maxDistance = getGlobal("toggleMaxDistance");

    system("mkdir $wrk/$toggledDir") if (! -d "$wrk/$toggledDir");

    #  A simple link to the ovlStore suffices.
    #
    if (! -e "$wrk/$toggledDir/$asm.ovlStore") {
        system("ln -s ../$asm.ovlStore $wrk/$toggledDir/$asm.ovlStore");
    }

    #  The gatekeeper store must be paritally copied, so that we can first undo
    #  clear range changes made by any previous ECR, and allow clear range changes
    #  to be made by future ECR.
    #
    if (! -e "$wrk/$toggledDir/$asm.gkpStore") {
        system("mkdir $wrk/$toggledDir/$asm.gkpStore");

        system("cd $wrk/$toggledDir/$asm.gkpStore && ln -s ../..//$asm.gkpStore/[qsu]?? .");
        system("cp $wrk/$asm.gkpStore/[filp]?? $wrk/$toggledDir/$asm.gkpStore/");
        system("cp $wrk/$asm.gkpStore/clr-*    $wrk/$toggledDir/$asm.gkpStore/");

        $cmd  = "$bin/gatekeeper";
        $cmd .= " --revertclear OBTCHIMERA $wrk/$toggledDir/$asm.gkpStore";
        $cmd .= " > $wrk/$toggledDir/$asm.gkpStore.resetClearRange.err 2>&1";

        if (runCommand("$wrk/$toggledDir", $cmd)) {
            caFailure("failed to get pre-ECR clear-ranges for toggling", "$wrk/$toggledDir/$asm.gkpStore.resetClearRange.err");
        }
    }

    #  The tigStore needs only a partial copy, and links suffice.
    #
    if (! -e "$wrk/$toggledDir/$asm.tigStore") {
        system("mkdir $wrk/$toggledDir/$asm.tigStore") ;

        system("cd $wrk/$toggledDir/$asm.tigStore && ln -s ../../$asm.tigStore/*v00[12345]* .");
    }

    system("ln -s ../4-unitigger $wrk/$toggledDir") if (! -e "$wrk/$toggledDir/4-unitigger");
    system("ln -s ../5-consensus $wrk/$toggledDir") if (! -e "$wrk/$toggledDir/5-consensus");

    #  Update the tigStore, flipping repeat untigs to unique unitigs.
    #
    if (! -e "$wrk/$toggledDir/toggled.success") {
        $cmd  = "$bin/markUniqueUnique ";
        $cmd .= " -a $wrk/9-terminator/$asm.asm ";
        $cmd .= " -t $wrk/$toggledDir/$asm.tigStore 5 ";
        $cmd .= " -l $minLength ";
        $cmd .= " -n $numInstances ";
        $cmd .= " -d $maxDistance ";
        $cmd .= " > $wrk/$toggledDir/toggle.err 2>&1";

        if (runCommand("$wrk/$toggledDir", $cmd)) {
            caFailure("failed to toggle unitigs ", "$wrk/$toggledDir/toggle.err");
        }

        touch("$wrk/$toggledDir/toggled.success");
    }

    my $numToggles = `tail -n 1 $wrk/$toggledDir/toggle.err | awk '{print \$2}'`;

    if ($numToggles == 0) {
        print "No toggling occured. Finished.\n";
        return;
    }

    my $oldwrk = $wrk;

    $wrk = "$wrk/$toggledDir";

    scaffolder();
    postScaffolderConsensus();
    terminate();
    cleaner();

    $wrk = $oldwrk;
}
################################################################################
################################################################################
################################################################################


#  Assembly all done, remove some of the crud.

sub cleaner () {
    my $cleanType = getGlobal("cleanup");
    my $cleanValu = 0;

    print STDERR "The Cleaner has arrived.  Doing '$cleanType'.\n";

    $cleanValu = 0  if ($cleanType =~ m/none/);
    $cleanValu = 1  if ($cleanType =~ m/light/);
    $cleanValu = 2  if ($cleanType =~ m/heavy/);
    $cleanValu = 3  if ($cleanType =~ m/aggressive/);


    if ($cleanValu >= 1) {
        #
        #  Remove some of the more useless output files,
        #  and many of the stores and whatnot that can be recreated.
        #
        rmrf("$asm.obtStore");
        rmrf("0-mercounts/*blocks", "0-mercounts/*sequence");
        rmrf("0-overlaptrim-overlap/overlap*out");
        rmrf("1-overlapper/overlap*out");
        rmrf("4-unitigger/$asm.fge", "4-unitigger/$asm.fgv");
        rmrf("7*/rezlog");
    }


    if ($cleanValu >= 2) {
        #
        #
        #
    }


    if ($cleanValu >= 3) {
        #
        #  Nuke everything except 9-terminator.  Be paranoid about doing it.
        #
        rmrf("0-mercounts");
        rmrf("0-overlaptrim");
        rmrf("0-overlaptrim-overlap");
        rmrf("1-overlapper");
        rmrf("2-frgcorr");
        rmrf("3-ovlcorr");
        rmrf("4-unitigger");
        rmrf("5-consensus");
        rmrf("7-[0-9]-CGW");
        rmrf("7-[0-9]-ECR");
        rmrf("7-CGW");
        rmrf("8-consensus");
        rmrf("$asm.SeqStore");
        rmrf("$asm.asm");
        rmrf("$asm.frg");
        rmrf("$asm.gkpStore");
        rmrf("$asm.obtStore");
        rmrf("$asm.ovlStore");
        rmrf("$asm.qc");
    }


    if ($cleanType =~ m/compress/) {
        #  Compress *.err (usually tiny)
        #  Compress overlaps (*ovb)
        #  Compress checkpoints (*ckp.*[0-9])
    }
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

    } elsif (($arg =~ /\.frg$|frg\.gz$|frg\.bz2$|frg\.xz$/i) && (-e $arg)) {
        $arg = "$ENV{'PWD'}/$arg" if ($arg !~ m!^/!);
        push @fragFiles, $arg;
        $commandLineOptions .= " \"$arg\"";

    } elsif (($arg =~ /\.sff$|sff\.gz$|sff\.bz2$|sff\.xz$/i) && (-e $arg)) {
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

preoverlap(@fragFiles);
merTrim();       #  merTrim() MUST be before overlapTrim().
overlapTrim();
createOverlapJobs("normal");
checkOverlap("normal");
createOverlapStore();
overlapCorrection();
classifyMates();
unitigger();
postUnitiggerConsensus();
scaffolder();
postScaffolderConsensus();
terminate();
toggler();
cleaner();

exit(0);
