package ca3g::Defaults;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(getCommandLineOptions addCommandLineOption writeLog caExit caFailure getNumberOfCPUs getPhysicalMemorySize diskSpace printHelp setParametersFromFile setParametersFromCommandLine checkParameters getGlobal setGlobal showErrorRates setErrorRate setDefaults);

use strict;
use Carp qw(cluck);
use Sys::Hostname;

use Filesys::Df;  #  for diskSpace()

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

    caFailure("parameter '$var' is not known", undef) if (!exists($global{$var}));

    return($global{$var});
}


sub setGlobalSpecialization ($@) {
    my $val = shift @_;

    foreach my $var (@_) {
        #print STDERR "set specialization $var = $val\n";
        $global{$var} = $val;
    }

    return(1);
}

sub setGlobal ($$) {
    my $var = shift @_;
    my $val = shift @_;
    my $set = 0;

    $var =~ tr/A-Z/a-z/;
    $val = undef  if ($val eq "");  #  Set to undefined, the default for many of the options.

    #  Translate from generic to specialized var

    foreach my $alg ("ovl", "mhap") {
        foreach my $opt ("usegrid", "gridoptions") {
            $set += setGlobalSpecialization($val, ("${opt}cor${alg}", "${opt}obt${alg}", "${opt}utg${alg}"))  if ($var eq "${opt}${alg}");
        }

        foreach my $opt ("memory", "threads", "concurrency") {
            $set += setGlobalSpecialization($val, ("cor${alg}${opt}", "obt${alg}${opt}", "utg${alg}${opt}"))  if ($var eq "${alg}${opt}");
        }
    }

    foreach my $opt ("overlapper") {
        $set += setGlobalSpecialization($val, ("cor${opt}", "obt${opt}", "utg${opt}"))  if ($var eq "${opt}");
    }

    #  e.g., corOvlHashBlockLength
    foreach my $opt ("ovlerrorrate", "ovlhashblocklength", "ovlrefblocksize", "ovlrefblocklength", "ovlhashbits", "ovlhashload", "ovlmersize", "ovlmerthreshold", "ovlmerdistinct", "ovlmertotal", "ovlfrequentmers") {
        $set += setGlobalSpecialization($val, ("cor${opt}", "obt${opt}", "utg${opt}"))  if ($var eq "${opt}");
    }

    #  e.g., corMhapBlockSize
    foreach my $opt ("mhapblocksize", "mhapmersize", "mhaprealign", "mhapsensitivity") {
        $set += setGlobalSpecialization($val, ("cor${opt}", "obt${opt}", "utg${opt}"))  if ($var eq "${opt}");
    }

    return  if ($set > 0);

    if ($var eq "errorrate") {
        setErrorRate($val);
        return;
    }

    caFailure("paramter '$var' is not known", undef) if (!exists($global{$var}));

    $global{$var} = $val;
}



sub setGlobalIfUndef ($$) {
    my $var = shift @_;
    my $val = shift @_;

    $var =~ tr/A-Z/a-z/;
    $val = undef  if ($val eq "");  #  Set to undefined, the default for many of the options.

    return  if (defined($global{$var}));

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



#  Use caExit() for transient errors, like not opening files, processes that die, etc.
sub caExit ($$) {
    my  $msg   = shift @_;
    my  $log   = shift @_;

    print STDERR "================================================================================\n";
    print STDERR "Don't panic, but a mostly harmless error occurred and ca3g failed.\n";
    print STDERR "\n";

    #  Really should pass in $wrk
    if (defined($log)) {
        my  $df = diskSpace($log);

        print STDERR "Disk space available:  $df GB\n";
        print STDERR "\n";
    }

    if (-e $log) {
        print STDERR "Last 50 lines of the relevant log file ($log):\n";
        print STDERR "\n";
        system("tail -n 50 $log");
        print STDERR "\n";
    }

    print STDERR "ca3g failed with '$msg'.\n";
    print STDERR "\n";

    exit(1);
}


#  Use caFailure() for errors that definitely will require code changes to fix.
sub caFailure ($$) {
    my  $msg   = shift @_;
    my  $log   = shift @_;

    print STDERR "================================================================================\n";
    print STDERR "Please panic.  ca3g failed, and it shouldn't have.\n";
    print STDERR "\n";
    print STDERR "Stack trace:\n";
    print STDERR "\n";
    cluck;
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


#
#  Host management - these really belong in 'Execution.pm' (or 'Utilities.pm') but can't go there
#  (Execution.pm) and be used here too.
#

sub getNumberOfCPUs () {
    my $os   = $^O;
    my $ncpu = 1;

    #  See http://stackoverflow.com/questions/6481005/obtain-the-number-of-cpus-cores-in-linux

    if ($os eq "freebsd") {
        $ncpu = int(`/sbin/sysctl -n hw.ncpu`);
    }

    if ($os eq "darwin") {
        $ncpu = int(`/usr/bin/getconf _NPROCESSORS_ONLN`);
    }

    if ($os eq "linux") {
        $ncpu = int(`getconf _NPROCESSORS_ONLN`);
    }

    return($ncpu);
}


sub getPhysicalMemorySize () {
    my $os     = $^O;
    my $memory = 1;

    if ($os eq "freebsd") {
        $memory = `/sbin/sysctl -n hw.physmem` / 1024 / 1024 / 1024;
    }

    if ($os eq "darwin") {
        $memory = `/usr/sbin/sysctl -n hw.memsize` / 1024 / 1024 / 1024;
    }

    if ($os eq "linux") {
        open(F, "< /proc/meminfo");        #  Way to go, Linux!  Make it easy on us!
        while (<F>) {
            if (m/MemTotal:\s+(\d+)/) {
                $memory = $1 / 1024 / 1024;
            }
        }
        close(F);
    }

    return(int($memory + 0.5));  #  Poor man's rounding
}


sub dirname ($) {
    my $d = shift @_;

    return($d)  if (-d $d);

    my @d = split '/', $d;
    pop @d;

    $d = join('/', @d);

    return($d);
}


sub diskSpace ($) {
    my $wrk   = dirname($_[0]);
    my $df    = df($wrk, 1024);

    my $total = int(10 * $df->{blocks} / 1048576) / 10;
    my $used  = int(10 * $df->{used}   / 1048576) / 10;
    my $free  = int(10 * $df->{bfree}  / 1048576) / 10;
    my $avail = int(10 * $df->{bavail} / 1048576) / 10;

    #print STDERR "Disk space: total $total GB, used $used GB, free $free GB, available $avail GB\n";

    return (wantarray) ? ($total, $used, $free, $avail) : $avail;
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
        foreach my $k (sort values %synnam) {
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
        print "\n";
        print "usage: ca3g.pl [run | correct | trim | assemble] \\\n";
        print "               -p <assembly-prefix> \\\n";
        print "               -d <assembly-directory> \\\n";
        print "               -s <assembly-specifications-file> \\\n";
        print "               genomeSize=Ng \\\n";
        print "               errorRate=0.X \\\n";
        print "               [other-options] \\\n";
        print "               [read-type *fastq]\n";
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

    open(F, "< $specFile") or caExit("can't open '$specFile' for reading: $!", undef);

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

    caExit("can't find 'gatekeeperCreate' program in $bin.  Possibly incomplete installation", undef) if (! -x "$bin/gatekeeperCreate");
    caExit("can't find 'meryl' program in $bin.  Possibly incomplete installation", undef)            if (! -x "$bin/meryl");
    caExit("can't find 'overlapInCore' program in $bin.  Possibly incomplete installation", undef)    if (! -x "$bin/overlapInCore");
    caExit("can't find 'bogart' program in $bin.  Possibly incomplete installation", undef)           if (! -x "$bin/bogart");
    caExit("can't find 'utgcns' program in $bin.  Possibly incomplete installation", undef)           if (! -x "$bin/utgcns");

    #
    #  Update obsolete usages.
    #

    setGlobal("unitigger", "unitigger")  if (getGlobal("unitigger") eq "utg");

    #
    #  Fiddle with filenames to make them absolute paths.
    #

    makeAbsolute("pathMap");

    makeAbsolute("corOvlFrequentMers");
    makeAbsolute("obtOvlFrequentMers");
    makeAbsolute("utgOvlFrequentMers");

    #
    #  Adjust case on some of them
    #

    fixCase("corOverlapper");
    fixCase("obtOverlapper");
    fixCase("utgOverlapper");

    fixCase("corConsensus");
    fixCase("cnsConsensus");

    fixCase("corFilter");

    fixCase("unitigger");
    fixCase("stopBefore");
    fixCase("stopAfter");

    #
    #  Check for invalid usage
    #

    foreach my $tag ("cor", "obt", "utg") {
        if ((getGlobal("${tag}Overlapper") ne "mhap") &&
            (getGlobal("${tag}Overlapper") ne "ovl")) {
            caExit("invalid '${tag}Overlapper' specified (" . getGlobal("${tag}Overlapper") . "); must be 'mhap' or 'ovl'", undef);
        }
    }

    if ((getGlobal("unitigger") ne "unitigger") &&
        (getGlobal("unitigger") ne "bogart")) {
        caExit("invalid 'unitigger' specified (" . getGlobal("unitigger") . "); must be 'unitigger' or 'bogart'", undef);
    }

    foreach my $tag ("cor", "cns") {
        if ((getGlobal("${tag}consensus") ne "utgcns") &&
            (getGlobal("${tag}consensus") ne "falcon") &&
            (getGlobal("${tag}consensus") ne "falconpipe") &&
            (getGlobal("${tag}consensus") ne "pbdagcon") &&
            (getGlobal("${tag}consensus") ne "pbutgcns")) {
            caExit("invalid 'consensus' specified (" . getGlobal("${tag}consensus") . "); must be 'utgcns' or 'falcon' or 'falconpipe' or 'pbdagcon' or 'pbutgcns'", undef);
        }
    }

    #if ((getGlobal("cleanup") ne "none") &&
    #    (getGlobal("cleanup") ne "light") &&
    #    (getGlobal("cleanup") ne "heavy") &&
    #    (getGlobal("cleanup") ne "aggressive")) {
    #    caExit("invalid cleaup specified (" . getGlobal("cleanup") . "); must be 'none', 'light', 'heavy' or 'aggressive'", undef);
    #}

    if ((getGlobal("corFilter") ne "quick") &&
        (getGlobal("corFilter") ne "expensive")) {
        caExit("invalid 'corFilter' specified (" . getGlobal("corFilter") . "); must be 'quick' or 'expensive'", undef);
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

        caExit($failureString, undef) if ($ok == 0);
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

        caExit($failureString, undef) if ($ok == 0);
    }

    caExit("required parameter 'errorRate' is not set", undef)   if (! defined(getGlobal("errorRate")));
    caExit("required parameter 'genomeSize' is not set", undef)  if (! defined(getGlobal("genomeSize")));

    setGlobal("genomeSize", $1 * 1000)        if (getGlobal("genomeSize") =~ m/(\d+.*\d+)k/i);
    setGlobal("genomeSize", $1 * 1000000)     if (getGlobal("genomeSize") =~ m/(\d+.*\d+)m/i);
    setGlobal("genomeSize", $1 * 1000000000)  if (getGlobal("genomeSize") =~ m/(\d+.*\d+)g/i);

    #
    #  Java?  Need JRE 1.8.
    #

    if ((getGlobal("corOverlapper") eq "mhap") ||
        (getGlobal("obtOverlapper") eq "mhap") ||
        (getGlobal("utgOverlapper") eq "mhap")) {
        my $java       = getGlobal("java");
        my $versionStr = "unknown";
        my $version    = 0;

        open(F, "$java -showversion 2>&1 |");
        while (<F>) {
            #  First word is either "java" or "openjdk" or ...
            if (m/^.*\s+version\s+\"(\d+.\d+)(.*)\"$/) {
                $versionStr = "$1$2";
                $version    =  $1;
            }
        }
        close(F);

        print STDERR "-- Detected Java(TM) Runtime Environment '$versionStr' (from '$java').\n";

        caExit("mhap overlapper requires java version at least 1.8.0; you have $versionStr", undef)  if ($version < 1.8);
    }

    #
    #  Falcon?  Need to find it.
    #

    if ((getGlobal("corConsensus") eq "falcon") ||
        (getGlobal("corConsensus") eq "falconpipe")) {
        my $falcon = getGlobal("falconSense");

        caExit("didn't find falcon program with option falconSense='$falcon'", undef)  if (! -e $falcon);
    }


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
        setGlobalIfUndef("gridEngineTaskID",                     "\$SGE_TASK_ID");
        setGlobalIfUndef("gridEngineArraySubmitID",              "\\\$TASK_ID");
        setGlobalIfUndef("gridEngineJobID",                      "JOB_ID");

        #  Try to figure out the name of the threaded job execution environment.
        #  It's the one with allocation_rule of $pe_slots.

        my @env = `qconf -spl`;  chomp @env;
        my $thr = undef;
        my $nth = 0;

        foreach my $env (@env) {
            my $ar = 0;
            my $cs = 0;
            my $jf = 0;

            open(F, "qconf -sp $env |");
            while (<F>) {
                $ar = 1  if (m/allocation_rule.*pe_slots/);
                $cs = 1  if (m/control_slaves.*FALSE/);
                $jf = 1  if (m/job_is_first_task.*TRUE/);
            }
            close(F);

            if (($ar == 1) && ($cs == 1) && ($jf == 1)) {
                $thr = $env;
                $nth++;
            }
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
        setGlobalIfUndef("gridEngineTaskID",                     "\$PBS_TASKNUM");
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
        setGlobalIfUndef("gridEngineTaskID",                     "\$LSB_JOBINDEX");
        setGlobalIfUndef("gridEngineArraySubmitID",              "%I");
        setGlobalIfUndef("gridEngineJobID",                      "LSB_JOBID");
    }

    #
    #  Set default error rates based on the per-read error rate.
    #

    setGlobalIfUndef("corOvlErrorRate",      3.0 * getGlobal("errorRate"));
    setGlobalIfUndef("obtOvlErrorRate",      3.0 * getGlobal("errorRate"));
    setGlobalIfUndef("utgOvlErrorRate",      3.0 * getGlobal("errorRate"));

    setGlobalIfUndef("utgGraphErrorRate",    3.0 * getGlobal("errorRate"));
    setGlobalIfUndef("utgBubbleErrorRate",   3.0 * getGlobal("errorRate") + 0.5 * getGlobal("errorRate"));
    setGlobalIfUndef("utgMergeErrorRate",    3.0 * getGlobal("errorRate") - 0.5 * getGlobal("errorRate"));
    setGlobalIfUndef("utgRepeatErrorRate",   3.0 * getGlobal("errorRate"));

    setGlobalIfUndef("cnsErrorRate",         3.0 * getGlobal("errorRate"));
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



sub showErrorRates ($) {
    my $prefix = shift @_;

    print STDERR "${prefix}\n";
    print STDERR "${prefix}genomeSize          -- ", getGlobal("genomeSize"), "\n";
    print STDERR "${prefix}errorRate           -- ", getGlobal("errorRate"), "\n";
    print STDERR "${prefix}\n";
    print STDERR "${prefix}corOvlErrorRate     -- ", getGlobal("corOvlErrorRate"), "\n";
    print STDERR "${prefix}obtOvlErrorRate     -- ", getGlobal("obtOvlErrorRate"), "\n";
    print STDERR "${prefix}utgOvlErrorRate     -- ", getGlobal("utgOvlErrorRate"), "\n";
    print STDERR "${prefix}\n";
    print STDERR "${prefix}obtErrorRate        -- ", getGlobal("obtErrorRate"), "\n";
    print STDERR "${prefix}\n";
    print STDERR "${prefix}utgGraphErrorRate   -- ", getGlobal("utgGraphErrorRate"), "\n";
    print STDERR "${prefix}utgBubbleErrorRate  -- ", getGlobal("utgBubbleErrorRate"), "\n";
    print STDERR "${prefix}utgMergeErrorRate   -- ", getGlobal("utgMergeErrorRate"), "\n";
    print STDERR "${prefix}utgRepeatErrorRate  -- ", getGlobal("utgRepeatErrorRate"), "\n";
    print STDERR "${prefix}\n";
    print STDERR "${prefix}cnsErrorRate        -- ", getGlobal("cnsErrorRate"), "\n";
    print STDERR "${prefix}\n";
}



#  Defaults are set for yeast:
#    trimming   errorRate = 0.009  obtOvlErrorRate = 0.06  obtErrorRate = 0.035
#    assembly   errorRate = 0.009  utgOvlErrorRate = 0.06  bogart 0.035
#
sub setErrorRate ($@) {
    my $er      = shift @_;
    my $verbose = shift @_;

    print STDERR "-- Set errorRate to $er (verbose='$verbose')\n"  if (defined($verbose));

    #  Can NOT call setGlobal() for this, because it calls setErrorRate()!.
    $global{"errorrate"} = $er;

    setGlobal("corOvlErrorRate",    $er * 6);  #  Not used, except for realigning
    setGlobal("obtOvlErrorRate",    $er * 6);  #  Generally must be smaller than utgGraphErrorRate
    setGlobal("utgOvlErrorRate",    $er * 6);

    setGlobal("obtErrorRate",       $er * 3);

    setGlobal("utgGraphErrorRate",  $er * 3);
    setGlobal("utgBubbleErrorRate", $er * 3 + 0.5 * $er);  #  Not tested!
    setGlobal("utgMergeErrorRate",  $er * 3 - 0.5 * $er);
    setGlobal("utgRepeatErrorRate", $er * 3);

    setGlobal("cnsErrorRate",       $er);

    showErrorRates("-- ")  if (defined($verbose));
}



sub setOverlapDefaults ($$$) {
    my $tag     = shift @_;  #  If 'cor', some parameters are loosened for raw pacbio reads
    my $name    = shift @_;
    my $default = shift @_;  #  Sets ${tag}Overlapper

    #  Which overlapper to use.

    $global{"${tag}Overlapper"}                  = $default;
    $synops{"${tag}Overlapper"}                  = "Which overlap algorithm to use for $name";

    #  OverlapInCore parameters.

    $global{"${tag}OvlHashBlockLength"}          = ($tag eq "cor") ? 2500000 : 100000000;
    $synops{"${tag}OvlHashBlockLength"}          = "Amount of sequence (bp) to load into the overlap hash table";

    $global{"${tag}OvlRefBlockSize"}             = ($tag eq "cor") ? 20000 : 2000000;
    $synops{"${tag}OvlRefBlockSize"}             = "Number of reads to search against the hash table per batch";

    $global{"${tag}OvlRefBlockLength"}           = 0;
    $synops{"${tag}OvlRefBlockLength"}           = "Amount of sequence (bp) to search against the hash table per batch";

    $global{"${tag}OvlHashBits"}                 = ($tag eq "cor") ? 18 : 22;
    $synops{"${tag}OvlHashBits"}                 = "Width of the kmer hash.  Width 22=1gb, 23=2gb, 24=4gb, 25=8gb.  Plus 10b per ${tag}OvlHashBlockLength";

    $global{"${tag}OvlHashLoad"}                 = 0.75;
    $synops{"${tag}OvlHashLoad"}                 = "Maximum hash table load.  If set too high, table lookups are inefficent; if too low, search overhead dominates run time";

    $global{"${tag}OvlMerSize"}                  = ($tag eq "cor") ? 19 : 22;
    $synops{"${tag}OvlMerSize"}                  = "K-mer size for seeds in overlaps";

    $global{"${tag}OvlMerThreshold"}             = "auto";
    $synops{"${tag}OvlMerThreshold"}             = "K-mer frequency threshold; mers more frequent than this count are ignored";

    $global{"${tag}OvlMerDistinct"}              = undef;
    $synops{"${tag}OvlMerDistinct"}              = "K-mer frequency threshold; the least frequent fraction of distinct mers can seed overlaps";

    $global{"${tag}OvlMerTotal"}                 = undef;
    $synops{"${tag}OvlMerTotal"}                 = "K-mer frequency threshold; the least frequent fraction of all mers can seed overlaps";

    $global{"${tag}OvlFrequentMers"}             = undef;
    $synops{"${tag}OvlFrequentMers"}             = "Do not seed overlaps with these kmers (fasta format)";

    #  Mhap parameters.

    $global{"${tag}MhapBlockSize"}            = 20000;
    $synops{"${tag}MhapBlockSize"}            = "Number of reads per block; one block is loaded into memory per job";

    $global{"${tag}MhapMerSize"}              = ($tag eq "cor") ? 16 : 22;
    $synops{"${tag}MhapMerSize"}              = "K-mer size for seeds in mhap";

    $global{"${tag}MhapReAlign"}              = undef;
    $synops{"${tag}MhapReAlign"}              = "Compute actual alignments from mhap overlaps; 'raw' from mhap output, 'final' from overlap store; uses either obtErrorRate or ovlErrorRate, depending on which overlaps are computed";

    $global{"${tag}MhapSensitivity"}          = "normal";
    $synops{"${tag}MhapSensitivity"}          = "Coarse sensitivity level: 'normal' or 'high'";
}



sub setDefaults () {

    #####  General Configuration Options (aka miscellany)

    $global{"ca3gIteration"}               = 0;  #  See documentation in Execution.pm
    $global{"ca3gIterationMax"}            = 2;

    $global{"showNext"}                    = undef;
    $synops{"showNext"}                    = "Don't run any commands, just report what would run";

    $global{"pathMap"}                     = undef;
    $synops{"pathMap"}                     = "File with a hostname to binary directory map";

    $global{"shell"}                       = "/bin/sh";
    $synops{"shell"}                       = "Command interpreter to use; sh-compatible (e.g., bash), NOT C-shell (csh or tcsh)";

    $global{"java"}                        = "java";
    #$global{"java"}                        = "/usr/bin/java"                               if (-e "/usr/bin/java");
    #$global{"java"}                        = "/usr/local/bin/java"                         if (-e "/usr/local/bin/java");
    $global{"java"}                        = "/nbacc/local/packages/jdk1.8.0_25/bin/java"  if (-e "/nbacc/local/packages/jdk1.8.0_25/bin/java");
    $synops{"java"}                        = "Java interpreter to use; at least version 1.8";

    #####  Cleanup options

    $global{"saveOverlaps"}                = 0;
    $synops{"saveOverlaps"}                = "Save intermediate overlap files, almost never a good idea";

    $global{"saveMerCounts"}               = 0;
    $synops{"saveMerCounts"}               = "Save full mer counting results, sometimes useful";

    #####  Error Rates

    $global{"errorRate"}                   = undef;
    $synops{"errorRate"}                   = "The expected error rate in the input reads";

    $global{"corOvlErrorRate"}             = undef;
    $synops{"corOvlErrorRate"}             = "Overlaps above this error rate are not computed";

    $global{"obtOvlErrorRate"}             = undef;
    $synops{"obtOvlErrorRate"}             = "Overlaps at or below this error rate are used to trim reads";

    $global{"utgOvlErrorRate"}             = undef;
    $synops{"utgOvlErrorRate"}             = "Overlaps at or below this error rate are used to trim reads";

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

    $global{"minReadLength"}               = 1000;
    $synops{"minReadLength"}               = "Reads shorter than this length are not loaded into the assembler";

    $global{"minOverlapLength"}            = 500;
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
        print STDERR "-- No grid engine detected, grid disabled.\n";
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

    setExecDefaults("cns",    "unitig consensus",                       1,  4, 1, undef);  #  Params are useGrid, memory, threads and concurrency
    setExecDefaults("cor",    "read correction",                        1,  4, 4, undef);  #    Default concurrency is n_cpu / n_threads
    setExecDefaults("red",    "read error detection",                   1,  4, 4, undef);  #    Default gridOptions are undef
    setExecDefaults("oea",    "overlap error adjustment",               1,  4, 1, undef);

    setExecDefaults("corovl",  "overlaps for correction",               1,  4, 4, undef);
    setExecDefaults("obtovl",  "overlaps for trimming",                 1,  4, 4, undef);
    setExecDefaults("utgovl",  "overlaps for unitig construction",      1,  4, 4, undef);

    setExecDefaults("cormhap", "mhap overlaps for correction",          1, 16, 8, undef);
    setExecDefaults("obtmhap", "mhap overlaps for trimming",            1, 16, 8, undef);
    setExecDefaults("utgmhap", "mhap overlaps for unitig construction", 1, 16, 8, undef);

    setExecDefaults("ovb",    "overlap store bucketizing",              1,  2, 1, undef);
    setExecDefaults("ovs",    "overlap store sorting",                  1,  4, 1, undef);
    setExecDefaults("master", "master script",                          0, 16, 1, undef);  #  Broken; bogart blows the limits

    #####  Overlapper

    setOverlapDefaults("cor", "correction",             "mhap");  #  Overlaps computed for correction
    setOverlapDefaults("obt", "overlap based trimming", "ovl");   #  Overlaps computed for trimming
    setOverlapDefaults("utg", "unitig construction",    "ovl");   #  Overlaps computed for unitigging

    ##### Overlap Store

    $global{"ovlStoreMemory"}              = 2048;
    $synops{"ovlStoreMemory"}              = "How much memory (MB) to use when constructing overlap stores";

    $global{"ovlStoreMethod"}              = "sequential";
    $synops{"ovlStoreMethod"}              = "Use the 'sequential' or 'parallel' algorithm for constructing an overlap store";

    $global{"ovlStoreSlices"}              = 128;
    $synops{"ovlStoreSlices"}              = "How many pieces to split the sorting into, for the parallel store build";

    #####  Mers

    $global{"merylMemory"}                 = undef;
    $synops{"merylMemory"}                 = "Amount of memory, in gigabytes, to use for mer counting";

    $global{"merylThreads"}                = undef;
    $synops{"merylThreads"}                = "Number of threads to use for mer counting";

    #####  Overlap Based Trimming

    $global{"obtErrorRate"}                = undef;
    $synops{"obtErrorRate"}                = "Stringency of overlaps to use for trimming";

    $global{"trimReadsOverlap"}            = 1;
    $synops{"trimReadsOverlap"}            = "Minimum overlap between evidence to make contiguous trim";

    $global{"trimReadsCoverage"}           = 1;
    $synops{"trimReadsCoverage"}           = "Minimum depth of evidence to retain bases";

    #$global{"splitReads..."}               = 1;
    #$synops{"splitReads..."}               = "";

    #####  Fragment/Overlap Error Correction

    $global{"enableOEA"}                   = 1;
    $synops{"enableOEA"}                   = "Do overlap error adjustment - comprises two steps: read error detection (RED) and overlap error adjustment (OEA)";

    $global{"redBatchSize"}                = 0;
    $synops{"redBatchSize"}                = "Number of reads per fragment error detection batch";

    $global{"redBatchLength"}              = 0;
    $synops{"redBatchLength"}              = "Number of bases per fragment error detection batch";

    $global{"oeaBatchSize"}                = 0;
    $synops{"oeaBatchSize"}                = "Number of reads per overlap error correction batch";

    $global{"oeaBatchLength"}              = 0;
    $synops{"oeaBatchLength"}              = "Number of bases per overlap error correction batch";

    #####  Unitigger & BOG & bogart Options

    $global{"unitigger"}                   = "bogart";
    $synops{"unitigger"}                   = "Which unitig algorithm to use; utg or bogart (defalut)";

    $global{"genomeSize"}                  = undef;
    $synops{"genomeSize"}                  = "An estimate of the size of the genome";

    $global{"utgBubblePopping"}            = 1;
    $synops{"utgBubblePopping"}            = "Smooth polymorphic regions";

    $global{"utgRecalibrateGAR"}           = 1;
    $synops{"utgRecalibrateGAR"}           = "Use an experimental algorithm to decide unique/repeat";

    $global{"batOptions"}                  = undef;
    $synops{"batOptions"}                  = "Advanced options to bogart";

    $global{"batMemory"}                   = undef;
    $synops{"batMemory"}                   = "Approximate maximum memory usage for loading overlaps, in gigabytes, default is unlimited";

    $global{"batThreads"}                  = undef;
    $synops{"batThreads"}                  = "Number of threads to use in the Merge/Split/Join phase; default is whatever OpenMP wants";

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

    $global{"cnsPartitionMin"}             = 75000;
    $synops{"cnsPartitionMin"}             = "Don't make a consensus partition with fewer than N reads";

    $global{"cnsMaxCoverage"}              = 0;
    $synops{"cnsMaxCoverage"}              = "Limit unitig consensus to at most this coverage";

    $global{"cnsConsensus"}                = "utgcns";
    $synops{"cnsConsensus"}                = "Which consensus algorithm to use; only 'utgcns' is supported";

    #####  Correction Options

    $global{"corPartitions"}               = 128;
    $synops{"corPartitions"}               = "Partition read correction into N jobs";

    $global{"corPartitionMin"}             = 25000;
    $synops{"corPartitionMin"}             = "Don't make a read correction partition with fewer than N reads";

    $global{"corMinEvidenceLength"}        = undef;
    $synops{"corMinEvidenceLength"}        = "Limit read correction to only overlaps longer than this; default: unlimited";

    $global{"corMaxEvidenceErate"}         = undef;
    $synops{"corMaxEvidenceErate"}         = "Limit read correction to only overlaps at or below this fraction error; default: unlimited";

    $global{"corMaxEvidenceCoverageGlobal"}= "1.0x";
    $synops{"corMaxEvidenceCoverageGlobal"}= "Limit reads used for correction to supporting at most this coverage; default: 1.0 * estimated coverage";

    $global{"corMaxEvidenceCoverageLocal"} = "1.5x";
    $synops{"corMaxEvidenceCoverageLocal"} = "Limit reads being corrected to at most this much evidence coverage; default: 1.5 * estimated coverage";

    $global{"corOutCoverage"}              = 40;
    $synops{"corOutCoverage"}              = "Only correct the longest reads up to this coverage; default 40";

    $global{"corFilter"}                   = "expensive";
    $synops{"corFilter"}                   = "Method to filter short reads from correction; 'quick' or 'expensive'";
    
    $global{"corConsensus"}                = "falconpipe";
    $synops{"corConsensus"}                = "Which consensus algorithm to use; only 'falcon' and 'falconpipe' are supported";

    $global{"falconSense"}                 = undef;
    $global{"falconSense"}                 = "/home/walenzb/ca3g/src/falcon_sense/falcon_sense.Linux-amd64.bin"                 if (-e "/home/walenzb/ca3g/src/falcon_sense/falcon_sense.Linux-amd64.bin");
    $global{"falconSense"}                 = "/nbacc/scratch/bri/ca3g/ca3g-build/src/falcon_sense/falcon_sense.Linux-amd64.bin" if (-e "/nbacc/scratch/bri/ca3g/ca3g-build/src/falcon_sense/falcon_sense.Linux-amd64.bin");
    $global{"falconSense"}                 = "/work/software/falcon/install/fc_env/bin/fc_consensus.py"                         if (-e "/work/software/falcon/install/fc_env/bin/fc_consensus.py");
    $synops{"falconSense"}                 = "Path to fc_consensus.py or falcon_sense.bin";

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

        next  if ($k eq "version");
        next  if ($k eq "options");
        next  if ($k eq "help");

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

    #  Finally, set the global default error rate

    setErrorRate(0.01);
}

1;
