package ca3g::Defaults;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(getCommandLineOptions addCommandLineOption writeLog caFailure printHelp setParametersFromFile setParametersFromCommandLine setParameters getGlobal setGlobal setDefaults);

use strict;
use Carp;
use Sys::Hostname;

my %global;
my %synops;
my %synnam;

my $cLineOpts = "";
my $specLog   = "";

sub getCommandLineOptions () {
    return($cLineOpts);
}

sub addCommandLineOption ($) {
    if ($cLineOpts =~ m/\s$/) {
        $cLineOpts .= "$_[0]";
    } else {
        $cLineOpts .= " $_[0]";
    }
}


sub writeLog ($) {
    my $wrk = shift @_;

    my $time = time();
    my $host = hostname();
    my $pid  = $$;

    open(F, "> $wrk/runCA-logs/${time}_${host}_${pid}_ca3g");
    print F $specLog;
    close(F);
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




sub printHelp ($) {
    my $bin = shift @_;  #  Can't include ca3g::Execution without a loop.

    if (getGlobal("version")) {
        system("$bin/gatekeeper    --version");
        system("$bin/overlapInCore --version");
        system("$bin/bogart        --version");
        system("$bin/utgcns        --version");
        exit(0);
    }

    if (getGlobal("options")) {
        foreach my $k (sort vals %synnam) {
            my $o = substr("$k                                    ", 0, 35);
            my $d = substr($global{$k}   . "                      ", 0, 20);
            my $u = $synops{$k};

            if (!defined($global{$k})) {
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

        addCommandLineOption("\"$var=$val\"");
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


sub setParameters ($) {
    my $bin = shift @_;  #  Can't include ca3g::Execution without a loop.

    #
    #  Update obsolete usages.
    #

    setGlobal("unitigger", "unitigger")  if (getGlobal("unitigger") eq "utg");

    #
    #  Fiddle with filenames to make them absolute paths.
    #

    makeAbsolute("pathMap");

    makeAbsolute("ovlFrequentMers");

    #
    #  Adjust case on some of them
    #

    fixCase("doChimeraDetection");
    fixCase("obtOverlapper");
    fixCase("ovlOverlapper");
    fixCase("unitigger");
    fixCase("stopBefore");
    fixCase("stopAfter");
    fixCase("consensus");
    fixCase("cleanup");

    #
    #  Check for invalid usage
    #

    if ((getGlobal("doChimeraDetection") ne "off") &&
        (getGlobal("doChimeraDetection") ne "normal") &&
        (getGlobal("doChimeraDetection") ne "aggressive")) {
        caFailure("invalid doChimeraDetection specified (" . getGlobal("doChimeraDetection") . "); must be 'off', 'normal', or 'aggressive'", undef);
    }
    if ((getGlobal("obtOverlapper") ne "mhap") &&
        (getGlobal("obtOverlapper") ne "ovl")) {
        caFailure("invalid obtOverlapper specified (" . getGlobal("obtOverlapper") . "); must be 'mer' or 'ovl' (or DEVEL_ONLY 'ovm')", undef);
    }
    if ((getGlobal("ovlOverlapper") ne "mhap") &&
        (getGlobal("ovlOverlapper") ne "ovl")) {
        caFailure("invalid ovlOverlapper specified (" . getGlobal("ovlOverlapper") . "); must be 'mer' or 'ovl' (or DEVEL_ONLY 'ovm')", undef);
    }
    if ((getGlobal("unitigger") ne "unitigger") &&
        (getGlobal("unitigger") ne "bogart")) {
        caFailure("invalid unitigger specified (" . getGlobal("unitigger") . "); must be 'unitigger' or 'bogart'", undef);
    }
    if ((getGlobal("consensus") ne "utgcns") &&
        (getGlobal("consensus") ne "seqan") &&
        (getGlobal("consensus") ne "pbdagcon") &&
        (getGlobal("consensus") ne "pbutgcns")) {
        caFailure("invalid consensus specified (" . getGlobal("consensus") . "); must be 'cns' or 'seqan' or 'pbdagcon' or 'pbutgcns'", undef);
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
                          "unitig");

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

        my @stopAfter = ("gatekeeper",
                         "meryl",
                         "overlap-configure",
                         "overlap",
                         "unitig",
                         "utgcns");

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

    #
    #  PIck a nice looking set of binaries, and check them.
    #

    caFailure("can't find 'gatekeeperCreate' program in $bin.  Possibly incomplete installation", undef) if (! -x "$bin/gatekeeperCreate");
    caFailure("can't find 'meryl' program in $bin.  Possibly incomplete installation", undef)            if (! -x "$bin/meryl");
    caFailure("can't find 'overlapInCore' program in $bin.  Possibly incomplete installation", undef)    if (! -x "$bin/overlapInCore");
    caFailure("can't find 'bogart' program in $bin.  Possibly incomplete installation", undef)           if (! -x "$bin/bogart");
    caFailure("can't find 'utgcns' program in $bin.  Possibly incomplete installation", undef)           if (! -x "$bin/utgcns");
}






#  Return the second argument, unless the first argument is found in
#  %global, in which case return that.
#
sub getGlobal ($) {
    my $var = shift @_;

    $var =~ tr/A-Z/a-z/;

    caFailure("script error -- global '$var' has no defined value", undef) if (!exists($global{$var}));

    return($global{$var});
}


sub setGlobal (%$$) {
    my $var = shift @_;
    my $val = shift @_;
    my $set = shift @_;

    $var =~ tr/A-Z/a-z/;

    #  If no value, set the field to undefined, the default for many of the options.

    $val = undef  if ($val eq "");

    #  Handle special cases.

    if (($var eq "gridEngine") && ($val eq "SGE")) {
        setGlobal("gridEngineSubmitCommand",      "qsub");
        setGlobal("gridEngineHoldOption",         "-hold_jid \"WAIT_TAG\"");
        setGlobal("gridEngineHoldOptionNoArray",  undef);
        setGlobal("gridEngineSyncOption",         "-sync y");
        setGlobal("gridEngineNameOption",         "-cwd -N");
        setGlobal("gridEngineArrayOption",        "-t ARRAY_JOBS");
        setGlobal("gridEngineArrayName",          "ARRAY_NAME");
        setGlobal("gridEngineOutputOption",       "-j y -o");
        setGlobal("gridEnginePropagateCommand",   "qalter -hold_jid \"WAIT_TAG\"");
        setGlobal("gridEngineNameToJobIDCommand", undef);
        setGlobal("gridEngineNameToJobIDCommandNoArray", undef);
        setGlobal("gridEngineTaskID",             "SGE_TASK_ID");
        setGlobal("gridEngineArraySubmitID",      "\\\$TASK_ID");
        setGlobal("gridEngineJobID",              "JOB_ID");
    }

    if (($var eq "gridEngine") && ($val eq "PBS")) {
        setGlobal("gridEngineSubmitCommand",      "qsub");
        setGlobal("gridEngineHoldOption",         "-W depend=afterany:\"WAIT_TAG\"");
        setGlobal("gridEngineHoldOptionNoArray",  undef);
        setGlobal("gridEngineSyncOption",         "");
        setGlobal("gridEngineNameOption",         "-d `pwd` -N");
        setGlobal("gridEngineArrayOption",        "-t ARRAY_JOBS");
        setGlobal("gridEngineArrayName",          "ARRAY_NAME\[ARRAY_JOBS\]");
        setGlobal("gridEngineOutputOption",       "-j oe -o");
        setGlobal("gridEnginePropagateCommand",   "qalter -W depend=afterany:\"WAIT_TAG\"");
        setGlobal("gridEngineNameToJobIDCommand", undef);
        setGlobal("gridEngineNameToJobIDCommandNoArray", undef);
        setGlobal("gridEngineTaskID",             "PBS_TASKNUM");
        setGlobal("gridEngineArraySubmitID",      "\\\$PBS_TASKNUM");
        setGlobal("gridEngineJobID",              "PBS_JOBID");
    }

    if (($var eq "gridEngine") && ($val eq "LSF")) {
        setGlobal("gridEngineSubmitCommand",      "bsub");
        setGlobal("gridEngineHoldOption",         "-w \"numended\(\"WAIT_TAG\", \*\)\"");
        setGlobal("gridEngineHoldOptionNoArray",  "-w \"done\(\"WAIT_TAG\"\)\"");
        setGlobal("gridEngineSyncOption",         "-K");
        setGlobal("gridEngineNameOption",         "-J");
        setGlobal("gridEngineArrayOption",        "");
        setGlobal("gridEngineArrayName",          "ARRAY_NAME\[ARRAY_JOBS\]");
        setGlobal("gridEngineOutputOption",       "-o");
        setGlobal("gridEnginePropagateCommand",   "bmodify -w \"done\(\"WAIT_TAG\"\)\"");
        setGlobal("gridEngineNameToJobIDCommand", "bjobs -A -J \"WAIT_TAG\" | grep -v JOBID");
        setGlobal("gridEngineNameToJobIDCommandNoArray", "bjobs -J \"WAIT_TAG\" | grep -v JOBID");
        setGlobal("gridEngineTaskID",             "LSB_JOBINDEX");
        setGlobal("gridEngineArraySubmitID",      "%I");
        setGlobal("gridEngineJobID",              "LSB_JOBID");
    }

    #  Update obsolete usage.


    #  Update aliases.

    $var = "doOverlapBasedTrimming"    if ($var eq "doOBT");
    $var = "doExtendClearRanges"       if ($var eq "doECR");

    #  If "help" exists, we're parsing command line options, and will catch this failure in
    #  printHelp().  Otherwise, this is an internal error, and we should bomb now.
    #
    if ((!defined($set)) && (!exists($global{$var}))) {
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

    $global{"utgErrorRate"}                = 0.030;
    $synops{"utgErrorRate"}                = "Overlaps at or below this error rate are used to construct unitigs (BOG and UTG)";

    $global{"utgGraphErrorRate"}           = 0.030;
    $synops{"utgGraphErrorRate"}           = "Overlaps at or below this error rate are used to construct unitigs (BOGART)";

    $global{"utgMergeErrorRate"}           = 0.045;
    $synops{"utgMergeErrorRate"}           = "Overlaps at or below this error rate are used to construct unitigs (BOGART)";

    $global{"cnsErrorRate"}                = 0.06;
    $synops{"cnsErrorRate"}                = "Consensus expects alignments at about this error rate";

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

    $global{"gridEngine"}                        = "SGE";

    $global{"gridEngineSubmitCommand"}           = "qsub";

    $global{"gridEngineHoldOption"}              = "-hold_jid \"WAIT_TAG\"";         # for lsf: -w "done("WAIT_TAG")"

    $global{"gridEngineHoldOptionNoArray"}       = undef;

    $global{"gridEngineSyncOption"}              = "-sync y";                        # for lsf: -K

    $global{"gridEngineNameOption"}              = "-cwd -N";                        # for lsf: -J

    $global{"gridEngineArrayOption"}             = "-t ARRAY_JOBS";                  # for lsf: empty ("")

    $global{"gridEngineArrayName"}               = "ARRAY_NAME";                     # for lsf: ARRAY_NAME[ARRAY_JOBS]

    $global{"gridEngineOutputOption"}            = "-j y -o";                        # for lsf: -o

    $global{"gridEnginePropagateCommand"}        = "qalter -hold_jid \"WAIT_TAG\"";  # for lsf: bmodify -w "done(WAIT_TAG)"

    $global{"gridEngineNameToJobIDCommand"}      = undef;                            # for lsf: bjobs -J "WAIT_TAG" | grep -v JOBID

    $global{"gridEngineNameToJobIDCommandNoArray"} = undef;

    $global{"gridEngineTaskID"}                  = "SGE_TASK_ID";                    # for lsf: LSB_JOBINDEX

    $global{"gridEngineArraySubmitID"}           = "\\\$TASK_ID";                    # for lsf: %I

    $global{"gridEngineJobID"}                   = "JOB_ID";                         # for lsf: LSB_JOBID

    #####  Sun Grid Engine

    $global{"useGrid"}                     = 0;
    $synops{"useGrid"}                     = "Enable SGE globally";

    $global{"useGridScript"}               = 0;
    $synops{"useGridScript"}               = "Enable SGE for runCA (and unitigger, scaffolder, other sequential phases)";

    $global{"useGridOVL"}                  = 1;
    $synops{"useGridOVL"}                  = "Enable SGE for overlap computations";

    $global{"useGridOVS"}                  = 0;
    $synops{"useGridOVS"}                  = "Enable OverlapStore Build on Grid";

    $global{"useGridFEC"}                  = 0;
    $synops{"useGridFEC"}                  = "Enable SGE for the fragment error correction";

    $global{"useGridOEC"}                  = 0;
    $synops{"useGridOEC"}                  = "Enable SGE for the overlap error correction";

    $global{"useGridCNS"}                  = 1;
    $synops{"useGridCNS"}                  = "Enable SGE for consensus";

    $global{"gridOptions"}                 = undef;
    $synops{"gridOptions"}                 = "SGE options applied to all SGE jobs";

    $global{"gridOptionsJobName"}          = undef;
    $synops{"gridOptionsJobName"}          = "SGE jobs name suffix";

    $global{"gridOptionsScript"}           = undef;
    $synops{"gridOptionsScript"}           = "SGE options applied to runCA jobs (and unitigger, scaffolder, other sequential phases)";

    $global{"gridOptionsOVL"}              = undef;
    $synops{"gridOptionsOVL"}              = "SGE options applied to overlap computation jobs";

    $global{"gridOptionsCNS"}              = undef;
    $synops{"gridOptionsCNS"}              = "SGE options applied to consensus jobs";

    $global{"gridOptionsFEC"}              = undef;
    $synops{"gridOptionsCEC"}              = "SGE options applied to fragment error correction jobs";

    $global{"gridOptionsOEC"}              = undef;
    $synops{"gridOptionsOEC"}              = "SGE options applied to overlap error correction jobs";

    #$global{"sgePropagateHold"}            = undef;
    #synnam{"sgePropagateHold"}            = "sgePropagateHold";
    #$synops{"sgePropagateHold"}            = undef;  #  Internal option

    #####  Gatekeeper (was 'Preoverlap')


    #####  Overlap Based Trimming

    $global{"doOverlapBasedTrimming"}      = 1;
    $synops{"doOverlapBasedTrimming"}      = "Enable the Overlap Based Trimming module (doOBT and doOverlapTrimming are aliases)";

    $global{"doDeDuplication"}             = 1;
    $synops{"doDeDuplication"}             = "Enable the OBT duplication detection and cleaning module for 454 reads, enabled automatically";

    $global{"doChimeraDetection"}          = "normal";
    $synops{"doChimeraDetection"}          = "Enable the OBT chimera detection and cleaning module; 'off', 'normal' or 'aggressive'";

    #####  Overlapper

    $global{"obtOverlapper"}               = "ovl";
    $synops{"obtOverlapper"}               = "Which overlap algorithm to use for OBT overlaps";

    $global{"ovlOverlapper"}               = "ovl";
    $synops{"ovlOverlapper"}               = "Which overlap algorithm to use for OVL (unitigger) overlaps";

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
    $synops{"ovlMerThreshold"}             = "K-mer frequency threshold; mers more frequent than this count are ignored";

    $global{"ovlMerDistinct"}              = undef;
    $synops{"ovlMerDistinct"}              = "K-mer frequency threshold; the least frequent fraction of distinct mers can seed overlaps";

    $global{"ovlMerTotal"}                 = undef;
    $synops{"ovlMerTotal"}                 = "K-mer frequency threshold; the least frequent fraction of all mers can seed overlaps";

    $global{"ovlFrequentMers"}             = undef;
    $synops{"ovlFrequentMers"}             = "Do not seed overlaps with these kmers (fasta format)";


    $global{"ovlHashLibrary"}               = "0";
    $synops{"ovlHashLibrary"}               = "For ovl overlaps, only load hash fragments from specified lib, 0 means all";

    $global{"ovlRefLibrary"}                = "0";
    $synops{"ovlRefLibrary"}                = "For ovl overlaps, only load ref fragments from specified lib, 0 means all";

    $global{"ovlCheckLibrary"}              = 1;
    $synops{"ovlCheckLibrary"}              = "Check that all libraries are used during ovl overlaps";


    $global{"merCompression"}              = 1;
    $synops{"merCompression"}              = "K-mer size";

    $global{"saveOverlaps"}                = 0;
    $synops{"saveOverlaps"}                = "Save intermediate overlap files";

    ##### Overlap Store

    $global{"ovlStoreMemory"}              = 2048;
    $synops{"ovlStoreMemory"}              = "How much memory (MB) to use when constructing overlap stores";

    $global{"ovlStoreMethod"}              = "sequential";
    $synops{"ovlStoreMethod"}              = "Use the 'sequential' or 'parallel' algorithm for constructing an overlap store";

    $global{"ovlStoreSlices"}              = 128;
    $synops{"ovlStoreSlices"}              = "How many pieces to split the sorting into, for the parallel store build";

    $global{"osbConcurrency"}              = 12;
    $synops{"osbConcurrency"}              = "How many bucketizing jobs to run concurrently, for the parallel store build";

    $global{"ossConcurrency"}              = 4;
    $synops{"ossConcurrency"}              = "How many sorting jobs to run concurrently, for the parallel store build";

    #####  Mers

    $global{"merylMemory"}                 = 4096;
    $synops{"merylMemory"}                 = "Amount of memory, in MB, to use for mer counting (conflicts with merylSegments)";

    $global{"merylSegments"}               = undef;
    $synops{"merylSegments"}               = "Number of segments to compute (overrides merylMemory)";

    $global{"merylThreads"}                = 1;
    $synops{"merylThreads"}                = "Number of threads to use for mer counting";

    #####  Fragment/Overlap Error Correction

    $global{"enableOEA"}                   = 1;
    $synops{"enableOEA"}                   = "Do overlap error adjustment - comprises two steps: read error detection (RED) and overlap error adjustment (OEA)";

    $global{"useGridRED"}                  = 0;
    $synops{"useGridRED"}                  = "Use grid engine for read error detection computes";

    $global{"redBatchSize"}                = 200000;
    $synops{"redBatchSize"}                = "Number of fragments per fragment error detection batch, directly affects memory usage";

    $global{"redThreads"}                  = 2;
    $synops{"redThreads"}                  = "Number of threads to use while computing fragment errors";

    $global{"redConcurrency"}              = 1;
    $synops{"redConcurrency"}              = "If not SGE, number of fragment error detection processes to run at the same time";

    $global{"useGridOEA"}                  = 0;
    $synops{"useGridOEA"}                  = "Use grid engine for overlap error adjustment computes";

    $global{"oeaBatchSize"}                = 200000;
    $synops{"oeaBatchSize"}                = "Number of fragments per overlap error correction batch";

    $global{"oeaConcurrency"}              = 4;
    $synops{"oeaConcurrency"}              = "If not SGE, number of overlap error correction processes to run at the same time";

    #####  Unitigger & BOG & bogart Options

    $global{"unitigger"}                   = "bogart";
    $synops{"unitigger"}                   = "Which unitig algorithm to use; utg or bogart (defalut)";

    $global{"utgGenomeSize"}               = 0;
    $synops{"utgGenomeSize"}               = "An estimate of the size of the genome; decides if unitigs are unique or repeats";

    $global{"utgBubblePopping"}            = 1;
    $synops{"utgBubblePopping"}            = "Smooth polymorphic regions";

    $global{"utgRecalibrateGAR"}           = 1;
    $synops{"utgRecalibrateGAR"}           = "Use an experimental algorithm to decide unique/repeat";

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

    #####  Consensus Options

    $global{"cnsPartitions"}               = 128;
    $synops{"cnsPartitions"}               = "Partition consensus into N jobs";

    $global{"cnsMinFrags"}                 = 75000;
    $synops{"cnsMinFrags"}                 = "Don't make a consensus partition with fewer than N fragments";

    $global{"cnsConcurrency"}              = 2;
    $synops{"cnsConcurrency"}              = "If not SGE, number of consensus jobs to run at the same time";

    $global{"cnsMaxCoverage"}              = 0;
    $synops{"cnsMaxCoverage"}              = "Limit unitig consensus to to at most this coverage";

    $global{"consensus"}                   = "utgcns";
    $synops{"consensus"}                   = "Which consensus algorithm to use; currently only 'cns' is supported";

    #####  Terminator Options

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

    #####  Ugly, command line options passed to printHelp()

    $global{"help"}                        = "";
    $synops{"help"}                        = undef;

    $global{"version"}                     = 0;
    $synops{"version"}                     = undef;

    $global{"options"}                     = 0;
    $synops{"options"}                     = undef;

    #  Convert all the keys to lowercase, and remember the case-sensitive version

    foreach my $k (keys %global) {
        (my $l = $k) =~ tr/A-Z/a-z/;

        if (! exists($synnam{$l})) {
            $synnam{$l} = $k;

            if (!exists($global{$l})) {
                $global{$l} = $global{$k};
                delete $global{$k};
            }

            #print "$k -> $l\n";
        }
    }

    #  If this is set, it breaks the consensus.sh and overlap.sh scripts.  Good grief!  Why
    #  are you running runCA in a task array!?

    if (exists($ENV{getGlobal("gridEngineTaskID")})) {
        undef $ENV{getGlobal("gridEngineTaskID")};
        print STDERR "ENV: ", getGlobal("gridEngineTaskID"), " needs to be unset, done.\n";
    }
}


1;
