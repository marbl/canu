package ca3g::Defaults;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(getCommandLineOptions addCommandLineOption writeLog caFailure printHelp setParametersFromFile setParametersFromCommandLine checkParameters getGlobal setGlobal setErrorRate setDefaults);

use strict;
use Carp;
use Sys::Hostname;

my %global;
my %synops;
my %synnam;

my $cLineOpts = "";
my $specLog   = "";



#  Return the second argument, unless the first argument is found in
#  %global, in which case return that.
#
sub getGlobal ($) {
    my $var = shift @_;

    $var =~ tr/A-Z/a-z/;

    caFailure("script error -- parameter '$var' is not known", undef) if (!exists($global{$var}));

    return($global{$var});
}



sub setGlobal ($$) {
    my $var = shift @_;
    my $val = shift @_;

    $var =~ tr/A-Z/a-z/;
    $val = undef  if ($val eq "");  #  Set to undefined, the default for many of the options.

    caFailure("script error -- paramter '$var' is not known", undef) if (!exists($global{$var}));

    $global{$var} = $val;
}



sub setGlobalIfUndef ($$) {
    my $var = shift @_;
    my $val = shift @_;

    return  if (defined($global{$var}));

    $var =~ tr/A-Z/a-z/;
    $val = undef  if ($val eq "");  #  Set to undefined, the default for many of the options.

    $global{$var} = $val;
}



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
    print STDERR "Stack trace:\n";
    print STDERR "\n";
    carp;
    print STDERR "\n";

    if (-e $log) {
        print STDERR "Last few lines of the relevant log file ($log):\n";
        print STDERR "\n";
        system("tail -n 50 $log");
    }
    print STDERR "\n";
    print STDERR "ca3g failed with '$msg'.\n";

    exit(1);
}



sub printHelp ($) {
    my $bin = shift @_;  #  Can't include ca3g::Execution without a loop.

    if (getGlobal("version")) {
        system("$bin/gatekeeperCreate --version");
        system("$bin/overlapInCore    --version");
        system("$bin/bogart           --version");
        system("$bin/utgcns           --version");
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
        print "usage: ca3g.pl [run | correct | trim | assemble] \\\n";
        print "               -p <assembly-prefix> \\\n";
        print "               -d <assembly-directory> \\\n";
        print "               -s <assembly-specifications-file> \\\n";
        print "               genomeSize=Ng \\\n";
        print "               errorRate=0.X \\\n";
        print "               [other-options] \\\n";
        print "               [read-type] *fastq\n";
        print "\n";
        print "    fully automatic modes:\n";
        print "      run      - generate an assembly, automagically applying the best correction, trimming\n";
        print "                 and assembly parameters\n";
        print "\n";
        print "    semi-automatic modes:\n";
        print "      correct  - generate corrected reads\n";
        print "      trim     - generate trimmed reads\n";
        print "      assemble - generate an assembly\n";
        print "\n";
        print "    a full list of options can be printed with '-options'.  All parameters except mode of\n";
        print "    operation, -p, -d and -s can be supplied in the sepc file.\n";
        print "\n";
        print "    reads can be either FASTA or FASTQ format, uncompressed, or compressed with gz, bz2 or xz:\n";
        print "      -pacbio-raw         \n";
        print "      -pacbio-corrected   \n";
        print "      -nanopore-raw       \n";
        print "      -nanopore-corrected \n";
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



sub checkParameters ($) {
    my $bin = shift @_;  #  Can't include ca3g::Execution without a loop.

    #
    #  PIck a nice looking set of binaries, and check them.
    #

    caFailure("can't find 'gatekeeperCreate' program in $bin.  Possibly incomplete installation", undef) if (! -x "$bin/gatekeeperCreate");
    caFailure("can't find 'meryl' program in $bin.  Possibly incomplete installation", undef)            if (! -x "$bin/meryl");
    caFailure("can't find 'overlapInCore' program in $bin.  Possibly incomplete installation", undef)    if (! -x "$bin/overlapInCore");
    caFailure("can't find 'bogart' program in $bin.  Possibly incomplete installation", undef)           if (! -x "$bin/bogart");
    caFailure("can't find 'utgcns' program in $bin.  Possibly incomplete installation", undef)           if (! -x "$bin/utgcns");

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

    fixCase("overlapper");
    fixCase("unitigger");
    fixCase("stopBefore");
    fixCase("stopAfter");
    fixCase("consensus");

    #
    #  Check for invalid usage
    #

    #if ((getGlobal("doChimeraDetection") ne "off") &&
    #    (getGlobal("doChimeraDetection") ne "normal") &&
    #    (getGlobal("doChimeraDetection") ne "aggressive")) {
    #    caFailure("invalid doChimeraDetection specified (" . getGlobal("doChimeraDetection") . "); must be 'off', 'normal', or 'aggressive'", undef);
    #}
    if ((getGlobal("overlapper") ne "mhap") &&
        (getGlobal("overlapper") ne "ovl")) {
        caFailure("invalid 'overlapper' specified (" . getGlobal("overlapper") . "); must be 'mhap' or 'ovl'", undef);
    }
    if ((getGlobal("unitigger") ne "unitigger") &&
        (getGlobal("unitigger") ne "bogart")) {
        caFailure("invalid 'unitigger' specified (" . getGlobal("unitigger") . "); must be 'unitigger' or 'bogart'", undef);
    }
    if ((getGlobal("consensus") ne "utgcns") &&
        (getGlobal("consensus") ne "falcon") &&
        (getGlobal("consensus") ne "falconpipe") &&
        (getGlobal("consensus") ne "pbdagcon") &&
        (getGlobal("consensus") ne "pbutgcns")) {
        caFailure("invalid 'consensus' specified (" . getGlobal("consensus") . "); must be 'utgcns' or 'falcon' or 'falconpipe' or 'pbdagcon' or 'pbutgcns'", undef);
    }
    #if ((getGlobal("cleanup") ne "none") &&
    #    (getGlobal("cleanup") ne "light") &&
    #    (getGlobal("cleanup") ne "heavy") &&
    #    (getGlobal("cleanup") ne "aggressive")) {
    #    caFailure("invalid cleaup specified (" . getGlobal("cleanup") . "); must be 'none', 'light', 'heavy' or 'aggressive'", undef);
    #}

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

    if (! defined(getGlobal("errorRate"))) {
        caFailure("ERROR: 'errorRate' is not set", undef);
    }

    if (! defined(getGlobal("genomeSize"))) {
        caFailure("ERROR: 'genomeSize' is not set", undef);
    }

    setGlobal("genomeSize", $1 * 1000)        if (getGlobal("genomeSize") =~ m/(\d+.*\d+)k/i);
    setGlobal("genomeSize", $1 * 1000000)     if (getGlobal("genomeSize") =~ m/(\d+.*\d+)m/i);
    setGlobal("genomeSize", $1 * 1000000000)  if (getGlobal("genomeSize") =~ m/(\d+.*\d+)g/i);

    #
    #  Finish grid configuration.  If any of these are set, they were set by the user.
    #

    #  Handle special cases.

    if (uc(getGlobal("gridEngine")) eq "SGE") {
        setGlobalIfUndef("gridEngineSubmitCommand",              "qsub");
        setGlobalIfUndef("gridEngineHoldOption",                 "-hold_jid \"WAIT_TAG\"");
        setGlobalIfUndef("gridEngineHoldOptionNoArray",          undef);
        setGlobalIfUndef("gridEngineSyncOption",                 "-sync y");
        setGlobalIfUndef("gridEngineNameOption",                 "-cwd -N");
        setGlobalIfUndef("gridEngineArrayOption",                "-t ARRAY_JOBS");
        setGlobalIfUndef("gridEngineArrayName",                  "ARRAY_NAME");
        setGlobalIfUndef("gridEngineOutputOption",               "-j y -o");
        setGlobalIfUndef("gridEnginePropagateCommand",           "qalter -hold_jid \"WAIT_TAG\"");
        setGlobalIfUndef("gridEngineThreadsOption",              undef);  #"-pe threads THREADS");
        setGlobalIfUndef("gridEngineMemoryOption",               "-l mem=MEMORY");
        setGlobalIfUndef("gridEngineNameToJobIDCommand",         undef);
        setGlobalIfUndef("gridEngineNameToJobIDCommandNoArray",  undef);
        setGlobalIfUndef("gridEngineTaskID",                     "SGE_TASK_ID");
        setGlobalIfUndef("gridEngineArraySubmitID",              "\\\$TASK_ID");
        setGlobalIfUndef("gridEngineJobID",                      "JOB_ID");

        #  Try to figure out the name of the threaded job execution environment.
        #  It's the one with allocation_rule of $pe_slots.

        my @env = `qconf -spl`;  chomp @env;
        my $thr = undef;
        my $nth = 0;

        foreach my $env (@env) {
            open(F, "qconf -sp $env |");
            while (<F>) {
                if (m/allocation_rule.*pe_slots/) {
                    $thr = $env;
                    $nth++;
                }
            }
            close(F);
        }

        if ($nth == 1) {
            print STDERR "-- Detected Grid Engine environment '$thr'.\n";

            setGlobal("gridEngineThreadsOption", "-pe $thr THREADS");

        } else {
            print STDERR "-- WARNING:  Couldn't determine the SGE parallel environment to run multi-threaded codes.\n";
            print STDERR "--           Set 'gridEngineThreadsOption' manually; example: '-pe threaded THREADS'.\n";
        }

        #  Try to figure out the name of the memory resource.

        my $mem = undef;
        my $nmm = 0;

        open(F, "qconf -sc |");
        while (<F>) {
            my @vals = split '\s+', $_;

            next  if ($vals[5] ne "YES");        #  Not a consumable resource.
            next  if ($vals[2] ne "MEMORY");     #  Not a memory resource.
            next  if ($vals[0] =~ m/swap/);      #  Don't care about swap.
            next  if ($vals[0] =~ m/virtual/);   #  Don't care about vm space.

            $mem .= " " if (defined($mem));
            $mem  = $vals[0];
            $nmm++;
        }
        close(F);

        if      ($nmm == 1) {
            print STDERR "-- Detected Grid Engine consumable '$mem'.\n";

            setGlobal("gridEngineMemoryOption", "-l $mem=MEMORY");

        } elsif ($nmm > 1) {
            print STDERR "-- WARNING:  Couldn't determine the SGE resource to request memory.\n";
            print STDERR "--           Found $nmm choices: $mem\n";
            print STDERR "--           Set 'gridEngineMemoryOption' manually; example: '-l mem=MEMORY'.\n";

        } else {
            print STDERR "-- WARNING:  Couldn't determine the SGE resource to request memory.\n";
            print STDERR "--           Set 'gridEngineMemoryOption' manually; example: '-l mem=MEMORY'.\n";
        }
    }

    if (uc(getGlobal("gridEngine")) eq "PBS") {
        setGlobalIfUndef("gridEngineSubmitCommand",              "qsub");
        setGlobalIfUndef("gridEngineHoldOption",                 "-W depend=afterany:\"WAIT_TAG\"");
        setGlobalIfUndef("gridEngineHoldOptionNoArray",          undef);
        setGlobalIfUndef("gridEngineSyncOption",                 "");
        setGlobalIfUndef("gridEngineNameOption",                 "-d `pwd` -N");
        setGlobalIfUndef("gridEngineArrayOption",                "-t ARRAY_JOBS");
        setGlobalIfUndef("gridEngineArrayName",                  "ARRAY_NAME\[ARRAY_JOBS\]");
        setGlobalIfUndef("gridEngineOutputOption",               "-j oe -o");
        setGlobalIfUndef("gridEnginePropagateCommand",           "qalter -W depend=afterany:\"WAIT_TAG\"");
        setGlobalIfUndef("gridEngineThreadsOption",              undef);
        setGlobalIfUndef("gridEngineMemoryOption",               undef);
        setGlobalIfUndef("gridEngineNameToJobIDCommand",         undef);
        setGlobalIfUndef("gridEngineNameToJobIDCommandNoArray",  undef);
        setGlobalIfUndef("gridEngineTaskID",                     "PBS_TASKNUM");
        setGlobalIfUndef("gridEngineArraySubmitID",              "\\\$PBS_TASKNUM");
        setGlobalIfUndef("gridEngineJobID",                      "PBS_JOBID");
    }

    if (uc(getGlobal("gridEngine")) eq "LSF") {
        setGlobalIfUndef("gridEngineSubmitCommand",              "bsub");
        setGlobalIfUndef("gridEngineHoldOption",                 "-w \"numended\(\"WAIT_TAG\", \*\)\"");
        setGlobalIfUndef("gridEngineHoldOptionNoArray",          "-w \"done\(\"WAIT_TAG\"\)\"");
        setGlobalIfUndef("gridEngineSyncOption",                 "-K");
        setGlobalIfUndef("gridEngineNameOption",                 "-J");
        setGlobalIfUndef("gridEngineArrayOption",                "");
        setGlobalIfUndef("gridEngineArrayName",                  "ARRAY_NAME\[ARRAY_JOBS\]");
        setGlobalIfUndef("gridEngineOutputOption",               "-o");
        setGlobalIfUndef("gridEnginePropagateCommand",           "bmodify -w \"done\(\"WAIT_TAG\"\)\"");
        setGlobalIfUndef("gridEngineThreadsOption",              undef);
        setGlobalIfUndef("gridEngineMemoryOption",               undef);
        setGlobalIfUndef("gridEngineNameToJobIDCommand",         "bjobs -A -J \"WAIT_TAG\" | grep -v JOBID");
        setGlobalIfUndef("gridEngineNameToJobIDCommandNoArray",  "bjobs -J \"WAIT_TAG\" | grep -v JOBID");
        setGlobalIfUndef("gridEngineTaskID",                     "LSB_JOBINDEX");
        setGlobalIfUndef("gridEngineArraySubmitID",              "%I");
        setGlobalIfUndef("gridEngineJobID",                      "LSB_JOBID");
    }

    #
    #  Set default error rates based on the per-read error rate.
    #

    setGlobalIfUndef("ovlErrorRate",      3.0 * getGlobal("errorRate"));
    setGlobalIfUndef("obtErrorRate",      3.0 * getGlobal("errorRate"));
    setGlobalIfUndef("utgErrorRate",      3.0 * getGlobal("errorRate"));
    setGlobalIfUndef("utgGraphErrorRate", 3.0 * getGlobal("errorRate"));
    setGlobalIfUndef("utgMergeErrorRate", 3.0 * getGlobal("errorRate"));
    setGlobalIfUndef("cnsErrorRate",      3.0 * getGlobal("errorRate"));

    #
    #  Report.
    #

    print STDERR "-- \n";
    print STDERR "-- genomeSize        = ", getGlobal("genomeSize"), "\n";
    print STDERR "-- errorRate         = ", getGlobal("errorRate"), "\n";
    print STDERR "-- \n";
    print STDERR "-- ovlErrorRate      = ", getGlobal("ovlErrorRate"), "\n";
    print STDERR "-- obtErrorRate      = ", getGlobal("obtErrorRate"), "\n";
    print STDERR "-- utgErrorRate      = ", getGlobal("utgErrorRate"), "\n";
    print STDERR "-- utgGraphErrorRate = ", getGlobal("utgGraphErrorRate"), "\n";
    print STDERR "-- utgMergeErrorRate = ", getGlobal("utgMergeErrorRate"), "\n";
    print STDERR "-- cnsErrorRate      = ", getGlobal("cnsErrorRate"), "\n";
    print STDERR "-- \n";
}


sub setExecDefaults ($$$$$$) {
    my $tag         = shift @_;
    my $name        = shift @_;
    my $usegrid     = shift @_;
    my $memory      = shift @_;
    my $threads     = shift @_;
    my $concurrent  = shift @_;

    $global{"useGrid${tag}"}       = $usegrid;
    $synops{"useGrid${tag}"}       = "Use grid engine for $name computes";

    $global{"gridOptions${tag}"}   = undef;
    $synops{"gridOptions${tag}"}   = "SGE options applied to $name jobs";

    $global{"${tag}Memory"}        = $memory;
    $synops{"${tag}Memory"}        = "Amount of memory, in gigabytes, to use for $name jobs";

    $global{"${tag}Threads"}       = $threads;
    $synops{"${tag}Threads"}       = "Number of threads to use for $name jobs";

    $global{"${tag}Concurrency"}   = $concurrent;
    $synops{"${tag}Concurrency"}   = "If not SGE, number of $name jobs to run at the same time; default is n_proc / n_threads";
}



sub setErrorRate ($) {
    my $er = shift @_;

    setGlobal("errorRate",          $er);

    setGlobal("ovlErrorRate",       $er * 3);
    setGlobal("obtErrorRate",       $er * 3);

    setGlobal("utgGraphErrorRate",  $er * 3);
    setGlobal("utgBubbleErrorRate", $er * 3 + 0.5 * $er);
    setGlobal("utgMergeErrorRate",  $er * 3 - 0.5 * $er);
    setGlobal("utgRepeatErrorRate", $er * 3);

    setGlobal("cnsErrorRate",       $er);
}


sub setDefaults () {

    #####  General Configuration Options (aka miscellany)

    $global{"showNext"}                    = undef;
    $synops{"showNext"}                    = "Don't run any commands, just report what would run";

    $global{"pathMap"}                     = undef;
    $synops{"pathMap"}                     = "File with a hostname to binary directory map";

    $global{"shell"}                       = "/bin/sh";
    $synops{"shell"}                       = "Command interpreter to use; sh-compatible (e.g., bash), NOT C-shell (csh or tcsh)";

    #####  Cleanup options

    $global{"saveOverlaps"}                = 0;
    $synops{"saveOverlaps"}                = "Save intermediate overlap files, almost never a good idea";

    #####  Error Rates

    $global{"errorRate"}                   = undef;
    $synops{"errorRate"}                   = "The expected error rate in the input reads";

    $global{"ovlErrorRate"}                = undef;
    $synops{"ovlErrorRate"}                = "Overlaps above this error rate are not computed";

    $global{"obtErrorRate"}                = undef;
    $synops{"obtErrorRate"}                = "Overlaps at or below this error rate are used to trim reads";

    #$global{"utgErrorRate"}                = undef;
    #$synops{"utgErrorRate"}                = "Overlaps at or below this error rate are used to construct unitigs (BOG and UTG)";

    $global{"utgGraphErrorRate"}           = undef;
    $synops{"utgGraphErrorRate"}           = "Overlaps at or below this error rate are used to construct unitigs (BOGART)";

    $global{"utgBubbleErrorRate"}          = undef;
    $synops{"utgBubbleErrorRate"}          = "Overlaps at or below this error rate are used to construct unitigs (BOGART)";

    $global{"utgMergeErrorRate"}           = undef;
    $synops{"utgMergeErrorRate"}           = "Overlaps at or below this error rate are used to construct unitigs (BOGART)";

    $global{"utgRepeatErrorRate"}          = undef;
    $synops{"utgRepeatErrorRate"}          = "Overlaps at or below this error rate are used to construct unitigs (BOGART)";

    $global{"cnsErrorRate"}                = undef;
    $synops{"cnsErrorRate"}                = "Consensus expects alignments at about this error rate";

    #####  Minimums

    $global{"minReadLength"}               = 500;
    $synops{"minReadLength"}               = "Reads shorter than this length are not loaded into the assembler";

    $global{"minOverlapLength"}            = 40;
    $synops{"minOverlapLength"}            = "Overlaps shorter than this length are not computed";

    #####  Stopping conditions

    $global{"stopBefore"}                  = undef;
    $synops{"stopBefore"}                  = "Tell ca3g when to halt execution";

    $global{"stopAfter"}                   = undef;
    $synops{"stopAfter"}                   = "Tell ca3g when to halt execution";

    #####  Grid Engine configuration, internal parameters

    $global{"gridEngine"}                           = undef;
    $global{"gridEngineSubmitCommand"}              = undef;
    $global{"gridEngineHoldOption"}                 = undef;
    $global{"gridEngineHoldOptionNoArray"}          = undef;
    $global{"gridEngineSyncOption"}                 = undef;
    $global{"gridEngineNameOption"}                 = undef;
    $global{"gridEngineArrayOption"}                = undef;
    $global{"gridEngineArrayName"}                  = undef;
    $global{"gridEngineOutputOption"}               = undef;
    $global{"gridEnginePropagateCommand"}           = undef;
    $global{"gridEngineThreadsOption"}              = undef;
    $global{"gridEngineMemoryOption"}               = undef;
    $global{"gridEngineNameToJobIDCommand"}         = undef;
    $global{"gridEngineNameToJobIDCommandNoArray"}  = undef;
    $global{"gridEngineTaskID"}                     = undef;
    $global{"gridEngineArraySubmitID"}              = undef;
    $global{"gridEngineJobID"}                      = undef;

    #  Try to decide which grid engine we have.

    if (defined($ENV{'SGE_ROOT'})) {
        print STDERR "-- Detected Sun Grid Engine in '$ENV{'SGE_ROOT'}/$ENV{'SGE_CELL'}'.\n";
        $global{"gridEngine"} = "SGE";
    } else {
        print STDERR "-- No grid engine detected.\n";
    }


    #####  Grid Engine Pipeline

    $global{"useGrid"}                     = 0;
    $synops{"useGrid"}                     = "Enable SGE; if unset, no grid will be used";

    #####  Grid Engine configuration, for each step of the pipeline

    $global{"gridOptions"}                 = undef;
    $synops{"gridOptions"}                 = "SGE options applied to all SGE jobs";

    $global{"gridOptionsJobName"}          = undef;
    $synops{"gridOptionsJobName"}          = "SGE jobs name suffix";

    #$global{"gridEngineCPUs"}              = undef;
    #$synops{"gridEngineCPUs"}              = "Number of CPUs available per grid host";

    #$global{"gridEngineMemory"}            = undef;
    #$synops{"gridEngineMemory"}            = "Amount of memory, in gigabytes, available per grid host";

    #####  Grid Engine configuration and parameters, for each step of the pipeline

    setExecDefaults("CNS",    "unitig consensus",           1,  4, 1, undef);  #  Params are useGrid, memory, threads and concurrency
    setExecDefaults("COR",    "read correction",            1,  4, 1, undef);  #  Default concurrency is n_cpu / n_threads
    setExecDefaults("RED",    "read error detection",       1,  4, 4, undef);
    setExecDefaults("OEA",    "overlap error adjustment",   1,  4, 4, undef);
    setExecDefaults("OVL",    "overlaps",                   1,  4, 4, undef);
    setExecDefaults("MHAP",   "mhap overlaps",              1, 16, 8, undef);
    setExecDefaults("OVB",    "overlap store bucketizing",  1,  2, 1, undef);
    setExecDefaults("OVS",    "overlap store sorting",      1,  4, 1, undef);
    setExecDefaults("MASTER", "master script",              0, 16, 1, undef);  #  Broken; bogart blows the limits

    #####  Overlapper, ovl

    $global{"overlapper"}                  = "ovl";
    $synops{"overlapper"}                  = "Which overlap algorithm to use for computing overlaps for assembly";

    $global{"ovlHashBlockLength"}          = 100000000;
    $synops{"ovlHashBlockLength"}          = "Amount of sequence (bp) to load into the overlap hash table";

    $global{"ovlRefBlockSize"}             = 2000000;
    $synops{"ovlRefBlockSize"}             = "Number of reads to search against the hash table per batch";

    $global{"ovlRefBlockLength"}           = 0;
    $synops{"ovlRefBlockLength"}           = "Amount of sequence (bp) to search against the hash table per batch";

    $global{"ovlHashBits"}                 = 22;
    $synops{"ovlHashBits"}                 = "Width of the kmer hash.  Width 22=1gb, 23=2gb, 24=4gb, 25=8gb.  Plus 10b per ovlHashBlockLength";

    $global{"ovlHashLoad"}                 = 0.75;
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

    #####  Overlapper, mhap

    $global{"mhapBlockSize"}               = 20000;
    $synops{"mhapBlockSize"}               = "Number of reads ....";

    $global{"mhapMerSize"}                 = 16;
    $synops{"mhapMerSize"}                 = "K-mer size for seeds in mhap";

    $global{"mhapReAlign"}                 = 0;
    $synops{"mhapReAlign"}                 = "Compute actual alignments from mhap overlaps; 'raw' from mhap output, 'final' from overlap store";

    $global{"mhapSensitivity"}             = "normal";
    $synops{"mhapSensitivity"}             = "Coarse sensitivity level: 'normal' or 'high'";

    #  PROBLEM: want to define mhap and ovl parameters independently, but then need to duplicate
    #  all the kmer stuff below for mhap.

    ##### Overlap Store

    $global{"ovlStoreMemory"}              = 2048;
    $synops{"ovlStoreMemory"}              = "How much memory (MB) to use when constructing overlap stores";

    $global{"ovlStoreMethod"}              = "sequential";
    $synops{"ovlStoreMethod"}              = "Use the 'sequential' or 'parallel' algorithm for constructing an overlap store";

    $global{"ovlStoreSlices"}              = 128;
    $synops{"ovlStoreSlices"}              = "How many pieces to split the sorting into, for the parallel store build";

    #####  Mers

    $global{"merylMemory"}                 = 4;
    $synops{"merylMemory"}                 = "Amount of memory, in gigabytes, to use for mer counting (conflicts with merylSegments)";

    $global{"merylSegments"}               = undef;
    $synops{"merylSegments"}               = "Number of segments to compute (overrides merylMemory)";

    $global{"merylThreads"}                = 1;
    $synops{"merylThreads"}                = "Number of threads to use for mer counting";

    #####  Overlap Based Trimming

    $global{"trimReadsOverlap"}            = 1;
    $synops{"trimReadsOverlap"}            = "Minimum overlap between evidence to make contiguous trim";

    $global{"trimReadsCoverage"}           = 1;
    $synops{"trimReadsCoverage"}           = "Minimum depth of evidence to retain bases";

    #$global{"splitReads..."}               = 1;
    #$synops{"splitReads..."}               = "";

    #####  Fragment/Overlap Error Correction

    $global{"enableOEA"}                   = 1;
    $synops{"enableOEA"}                   = "Do overlap error adjustment - comprises two steps: read error detection (RED) and overlap error adjustment (OEA)";

    $global{"redBatchSize"}                = 200000;
    $synops{"redBatchSize"}                = "Number of reads per fragment error detection batch, directly affects memory usage";

    $global{"oeaBatchSize"}                = 200000;
    $synops{"oeaBatchSize"}                = "Number of reads per overlap error correction batch";

    #####  Unitigger & BOG & bogart Options

    $global{"unitigger"}                   = "bogart";
    $synops{"unitigger"}                   = "Which unitig algorithm to use; utg or bogart (defalut)";

    $global{"genomeSize"}                  = undef;
    $synops{"genomeSize"}                  = "An estimate of the size of the genome";

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
    $synops{"cnsMinFrags"}                 = "Don't make a consensus partition with fewer than N reads";

    $global{"cnsMaxCoverage"}              = 0;
    $synops{"cnsMaxCoverage"}              = "Limit unitig consensus to to at most this coverage";

    $global{"consensus"}                   = "utgcns";
    $synops{"consensus"}                   = "Which consensus algorithm to use; currently only 'cns' is supported";

    $global{"falcon"}                      = undef;
    $global{"falcon"}                      = "/home/walenzb/ca3g/ca3g-build/src/falcon_sense/falcon_sense.Linux-amd64.bin" if (-e "/home/walenzb/ca3g/ca3g-build/src/falcon_sense/falcon_sense.Linux-amd64.bin");
    $global{"falcon"}                      = "/nbacc/scratch/bri/ca3g-build/src/falcon_sense/falcon_sense.Linux-amd64.bin" if (-e "/nbacc/scratch/bri/ca3g-build/src/falcon_sense/falcon_sense.Linux-amd64.bin");
    $global{"falcon"}                      = "/work/software/falcon/install/fc_env/bin/fc_consensus.py"                    if (-e "/work/software/falcon/install/fc_env/bin/fc_consensus.py");
    $synops{"falcon"}                      = "Path to fc_consensus.py";

    $global{"falconThreads"}               = 4;
    $synops{"falconThreads"}               = "Number of compute threads to use for falcon consensus";

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
    #  are you running this in a task array!?

    if (exists($ENV{getGlobal("gridEngineTaskID")})) {
        undef $ENV{getGlobal("gridEngineTaskID")};
        print STDERR "ENV: ", getGlobal("gridEngineTaskID"), " needs to be unset, done.\n";
    }
}

1;
