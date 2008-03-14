use strict;


# Submit batch jobs in groups. This function allows one to limit 
# the maximum number of jobs submitted at once to the grid. The jobs are 
# split into multiple submissions, dependent on each other.
#
# NOTE: inefficient because all jobs in a previous 
# submission must finish for the next submission to be scheduled.
#
sub submitBatchJobs($$$$) {
   my $id = shift @_;
   my $sgeParam = shift @_;
   my $jobs = shift @_;
   my $jobsSubmitted = 0;
   my $max = 1;
   my $min = 1;
   my $SGE = "";
   my $numThreads = shift @_;

   my $maxSize = getGlobal("maxGridJobSize");
   if (defined($maxSize)) {
      $maxSize = floor($maxSize / $numThreads);
   }
   else {
      $maxSize = $jobs;
   }

   while ($jobsSubmitted < $jobs) {
      my $max = $jobsSubmitted + $maxSize;
      my $min = $jobsSubmitted + 1;

      if ($max > $jobs) {
         $max = $jobs;
      }

       $SGE = $sgeParam;
       $SGE =~ s/NAME/"$id\_$asm\_$max"/e;
       if ($jobsSubmitted != 0) {
          $SGE =~ s/MINMAX/"$min-$max -hold_jid \"$id\_$asm\_$jobsSubmitted\" "/e;
       }
       else {
          $SGE =~ s/MINMAX/"$min-$max"/e;
       }

       if (runningOnGrid()) {
          system($SGE) and caFailure("Failed to submit overlap jobs.");
       } else {
          pleaseExecute($SGE);
       }
       $jobsSubmitted = $max;
   }

   return "$id\_$asm\_$jobs";
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
        open(F, "< $pathMap") or caFailure("Failed to open pathMap '$pathMap'.\n");
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
        open(PM, "< $pathMap") or caFailure("Failed to open pathMap '$pathMap'.\n");
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
    caFailure("ERROR: $var has no defined value!\n") if (!exists($global{$var}));
    return($global{$var});
}

sub setGlobal ($$) {
    my $var = shift @_;
    my $val = shift @_;
    #  If no value, set the field to undefined, the default for many of the options.
    if ($val eq "") {
        $val = undef;
    }
    #  Special case -- merSize sets both obtMerSize and ovlMerSize.
    if ($var eq "merSize") {
        setGlobal("obtMerSize", $val);
        setGlobal("ovlMerSize", $val);
        return;
    }
    #  Special case -- overlapper sets both obtOverlapper and ovlOverlapper.
    if ($var eq "overlapper") {
        setGlobal("obtOverlapper", $val);
        setGlobal("ovlOverlapper", $val);
        return;
    }
    caFailure("ERROR: $var is not a valid option.\n") if (!exists($global{$var}));
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

    #####  Error Rates

    $global{"ovlErrorRate"}                = 0.06;
    $global{"utgErrorRate"}                = 15;
    $global{"cnsErrorRate"}                = 0.06;
    $global{"cgwErrorRate"}                = 0.10;

    #####  Stopping conditions

    $global{"stopAfter"}                   = undef;

    #####  General Configuration Options (aka miscellany)

    $global{"doBackupFragStore"}           = 1;

    $global{"fakeUIDs"}                    = 0;
    $global{"uidServer"}                   = undef;

    $global{"pathMap"}                     = undef;

    #####  Sun Grid Engine

    $global{"useGrid"}                     = 0;
    $global{"scriptOnGrid"}                = 0;

    $global{"ovlOnGrid"}                   = 1;
    $global{"frgCorrOnGrid"}               = 0;
    $global{"ovlCorrOnGrid"}               = 0;
    $global{"cnsOnGrid"}                   = 1;

    $global{"maxGridJobSize"}		   = undef;

    $global{"sge"}                         = undef;
    $global{"sgeScript"}                   = undef;
    $global{"sgeOverlap"}                  = undef;
    $global{"sgeConsensus"}                = undef;
    $global{"sgeFragmentCorrection"}       = undef;
    $global{"sgeOverlapCorrection"}        = undef;

    #####  Preoverlap

    $global{"gkpFixInsertSizes"}           = 1;

    #####  Overlap Based Trimming

    $global{"doOverlapTrimming"}           = 1;
    $global{"vectorIntersect"}             = undef;
    $global{"vectorTrimmer"}               = "ca";
    $global{"figaroFlags"}                 = "-T 30 -M 100 -E 500 -V f";

    #####  Overlapper

    $global{"obtOverlapper"}               = "ovl";
    $global{"ovlOverlapper"}               = "ovl";
    $global{"ovlStoreMemory"}              = 1024;

    $global{"ovlThreads"}                  = 2;
    $global{"ovlStart"}                    = 1;
    $global{"ovlHashBlockSize"}            = 200000;
    $global{"ovlRefBlockSize"}             = 2000000;
    $global{"ovlMemory"}                   = "2GB";

    $global{"ovlMerSize"}                  = 22;
    $global{"ovlMerThreshold"}             = 500;

    $global{"obtMerSize"}                  = 22;
    $global{"obtMerThreshold"}             = 1000;

    $global{"merCompression"}              = 1;
    $global{"merOverlapperThreads"}        = 2;
    $global{"merOverlapperSeedBatchSize"}  = 100000;
    $global{"merOverlapperExtendBatchSize"}= 75000;

    $global{"merOverlapperSeedConcurrency"}  = 1;
    $global{"merOverlapperExtendConcurrency"}= 1;

    $global{"umdOverlapperFlags"}          = "-use-uncleaned-reads -trim-error-rate 0.03 -max-minimizer-cutoff 150";

    #####  Mers

    $global{"merylMemory"}                 = 800;
    $global{"merylThreads"}                = 1;

    #####  Fragment/Overlap Error Correction

    $global{"frgCorrBatchSize"}            = 200000;
    $global{"doFragmentCorrection"}        = 1;
    $global{"frgCorrThreads"}              = 2;
    $global{"frgCorrConcurrency"}          = 1;
    $global{"ovlCorrBatchSize"}            = 200000;
    $global{"ovlCorrConcurrency"}          = 4;

    #####  Unitigger & BOG Options

    $global{"unitigger"}                   = "utg";

    $global{"utgGenomeSize"}               = undef;
    $global{"utgEdges"}                    = undef;
    $global{"utgFragments"}                = undef;
    $global{"utgBubblePopping"}            = 1;
    $global{"utgRecalibrateGAR"}           = 1;

    $global{"bogPromiscuous"}              = 0;
    $global{"bogEjectUnhappyContain"}      = 0;
    $global{"bogBadMateDepth"}             = 7;

    #####  Scaffolder Options

    $global{"cgwOutputIntermediate"}       = 0;
    $global{"cgwPurgeCheckpoints"}         = 1;
    $global{"cgwDemoteRBP"}                = 1;

    $global{"astatLowBound"}		   = 1;
    $global{"astatHighBound"}		   = 5;

    $global{"stoneLevel"}                  = 2;

    $global{"computeInsertSize"}           = 0;
    $global{"cgwDistanceSampleSize"}       = 100;

    $global{"doResolveSurrogates"}         = 1;

    $global{"doExtendClearRanges"}         = 2;
    $global{"extendClearRangesStepSize"}   = undef;

    #####  Consensus Options

    $global{"cnsPartitions"}               = 128;
    $global{"cnsMinFrags"}                 = 75000;
    $global{"cnsConcurrency"}              = 2;
    $global{"consensus"}                   = "cns";

    #####  Terminator Options

    $global{"createAGP"}                   = 0;
    $global{"createACE"}                   = 0;
    $global{"createPosMap"}                = 1;

    $global{"merQC"}                       = 0;
    $global{"merQCmemory"}                 = 1024;
    $global{"merQCmerSize"}                = 22;

    $global{"cleanup"}                     = "none";

    #####  Ugly, command line options passed to printHelp()

    $global{"help"}                        = 0;
    $global{"version"}			   = 0;
    $global{"fields"}			   = 0;

    $global{"specFile"}                    = undef;
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

sub setParameters ($@) {
    my $specFile = shift @_;
    my @specOpts = @_;

    if (exists($ENV{'AS_OVL_ERROR_RATE'})) {
        setGlobal("ovlErrorRate", $ENV{'AS_OVL_ERROR_RATE'});
        print STDERR "ENV: ovlErrorRate $ENV{'AS_OVL_ERROR_RATE'}\n";
    }
    if (exists($ENV{'AS_CGW_ERROR_RATE'})) {
        setGlobal("cgwErrorRate", $ENV{'AS_CGW_ERROR_RATE'});
        print STDERR "cgwErrorRate $ENV{'AS_CGW_ERROR_RATE'}\n";
    }
    if (exists($ENV{'AS_CNS_ERROR_RATE'})) {
        setGlobal("cnsErrorRate", $ENV{'AS_CNS_ERROR_RATE'});
        print STDERR "cnsErrorRate $ENV{'AS_CNS_ERROR_RATE'}\n";
    }

    #  If the user didn't give us a specFile, see if there is a
    #  system-wide one defined (set by your localDefaults()).
    #
    if( !defined($specFile) ) {
    	$specFile = getGlobal("specFile");
    }

    if (defined($specFile)) {
        my $bin = "$FindBin::RealBin/spec";

        if (-e $specFile && ! -d $specFile) {
            open(F, "< $specFile") or caFailure("Couldn't open '$specFile'\n");
        } elsif (-e "$bin/$specFile") {
            open(F, "< $bin/$specFile") or caFailure("Couldn't open '$bin/$specFile'\n");
        } elsif (-e "$bin/$specFile.specFile") {
            open(F, "< $bin/$specFile.specFile") or caFailure("Couldn't open '$bin/$specFile.specFile'\n");
        } else {
            caFailure("You gave me a specFile, but I couldn't find '$specFile' or '$bin/$specFile' or '$bin/$specFile.specFile'!\n");
        }
        while (<F>) {
            chomp;
            next if (m/^\s*\#/);
            next if (m/^\s*$/);
            if (m/\s*(\w*)\s*=(.*)/) {
                my ($var, $val) = ($1, $2);
                print STDERR $_,"\n"; # echo the spec file
                $var =~ s/^\s+//; $var =~ s/\s+$//;
                $val =~ s/^\s+//; $val =~ s/\s+$//;
                undef $val if ($val eq "undef");
                setGlobal($var, $val);
            } else {
                print STDERR "WARNING!  Invalid specFile line '$_'\n";
            }
        }
        close(F);
    }

    foreach my $s (@specOpts) {
        if ($s =~ m/\s*(\w*)\s*=(.*)/) {
            my ($var, $val) = ($1, $2);
            $var =~ s/^\s+//; $var =~ s/\s+$//;
            $val =~ s/^\s+//; $val =~ s/\s+$//;
            setGlobal($var, $val);
        } else {
            print STDERR "WARNING!  Misformed specOption '$s'\n";
        }
    }

    #  Fiddle with filenames to make them absolute paths.
    #
    makeAbsolute("vectorIntersect");
    makeAbsolute("pathMap");

    #  PIck a nice looking set of binaries, and check them.
    #
    {
        my $bin = getBinDirectory();

        caFailure("Can't find 'gatekeeper' program in $bin.  Possibly incomplete installation.\n") if (! -x "$bin/gatekeeper");
        caFailure("Can't find 'meryl' program in $bin.  Possibly incomplete installation.\n") if (! -x "$bin/meryl");
        caFailure("Can't find 'overlap' program in $bin.  Possibly incomplete installation.\n") if (! -x "$bin/overlap");
        caFailure("Can't find 'unitigger' program in $bin.  Possibly incomplete installation.\n") if (! -x "$bin/unitigger");
        caFailure("Can't find 'cgw' program in $bin.  Possibly incomplete installation.\n") if (! -x "$bin/cgw");
        caFailure("Can't find 'consensus' program in $bin.  Possibly incomplete installation.\n") if (! -x "$bin/consensus");
        caFailure("Can't find 'terminator' program in $bin.  Possibly incomplete installation.\n") if (! -x "$bin/terminator");
    }

    #  Set the globally accessible error rates.  Adjust them if they
    #  look strange.
    #
    #  We must have:     ovl <= cns <= cgw
    #  We usually have:  ovl == cns <= cgw
    #
    my $ovlER = getGlobal("ovlErrorRate");
    my $cgwER = getGlobal("cgwErrorRate");
    my $cnsER = getGlobal("cnsErrorRate");

    if (($ovlER < 0.0) || (0.25 < $ovlER)) {
        caFailure("ovlErrorRate is $ovlER, this MUST be between 0.0 and 0.25.\n");
    }
    if (($cgwER < 0.0) || (0.25 < $cgwER)) {
        caFailure("cgwErrorRate is $cgwER, this MUST be between 0.0 and 0.25.\n");
    }
    if (($cnsER < 0.0) || (0.25 < $cnsER)) {
        caFailure("cnsErrorRate is $cnsER, this MUST be between 0.0 and 0.25.\n");
    }
    if ($ovlER > $cnsER) {
        caFailure("ovlErrorRate is $ovlER, this MUST be <= cnsErrorRate ($cnsER)\n");
    }
    if ($ovlER > $cgwER) {
        caFailure("ovlErrorRate is $ovlER, this MUST be <= cgwErrorRate ($cgwER)\n");
    }
    if ($cnsER > $cgwER) {
        caFailure("cnsErrorRate is $cnsER, this MUST be <= cgwErrorRate ($cgwER)\n");
    }
    $ENV{'AS_OVL_ERROR_RATE'} = $ovlER;
    $ENV{'AS_CGW_ERROR_RATE'} = $cgwER;
    $ENV{'AS_CNS_ERROR_RATE'} = $cnsER;
}


sub printHelp () {

    if ( getGlobal("version")) {
	my @full_pName = split('/',$0);
	my $pName = $full_pName[$#full_pName]; 
	print "$pName $MY_VERSION\n";
	exit(0);
    }

    if (getGlobal("help")) {
        print $HELPTEXT;
	exit(0);
    }

    if ( getGlobal("fields") ) {
        foreach my $k (sort keys %global) {
            if (defined(getGlobal($k))) {
                print substr("$k                             ", 0, 30) . getGlobal($k) . "\n";
            } else {
                print substr("$k                             ", 0, 30) . "<not defined>\n";
            }
        }
        exit(0);
    }
}



sub checkDirectories () {

    #  Check that we were supplied a work directory, and that it
    #  exists, or we can create it.
    #
    die "ERROR: I need a directory to run the assembly in (-d option).\n" if (!defined($wrk));

    system("mkdir -p $wrk") if (! -d $wrk);
    chmod 0755, "$wrk";

    caFailure("ERROR: Directory '$wrk' doesn't exist (-d option) and couldn't be created.\n") if (! -d $wrk);
}


sub findFirstCheckpoint ($) {
    my $dir      = shift @_;
    my $firstckp = 0;

    $dir = "$wrk/$dir" if (! -d $dir);

    open(F, "ls -1 $dir/*ckp* |");
    while (<F>) {
        chomp;

        if (m/ckp.(\d+)$/) {
            $firstckp = $1 if ($1 < $firstckp);
        } else {
            print STDERR "WARNING: Didn't match $_\n";
        }
    }
    close(F);

    return($firstckp);
}

sub findLastCheckpoint ($) {
    my $dir     = shift @_;
    my $lastckp = 0;

    $dir = "$wrk/$dir" if (! -d $dir);

    open(F, "ls -1 $dir/*ckp* |");
    while (<F>) {
        chomp;

        if (m/ckp.(\d+)$/) {
            $lastckp = $1 if ($1 > $lastckp);
        } else {
            print STDERR "WARNING: Didn't match $_\n";
        }
    }
    close(F);

    return($lastckp);
}

sub findNumScaffoldsInCheckpoint ($$) {
    my $dir     = shift @_;
    my $lastckp = shift @_;
    my $bin     = getBinDirectory();

    open(F, "cd $wrk/$dir && $bin/getNumScaffolds ../$asm.gkpStore $asm $lastckp 2> /dev/null |");
    my $numscaf = <F>;  chomp $numscaf;
    close(F);
    $numscaf = int($numscaf);

    caFailure("findNumScaffoldsInCheckpoint($dir, $lastckp) found no scaffolds?!") if ($numscaf == 0);
    print STDERR "Found $numscaf scaffolds in $dir checkpoint number $lastckp.\n";
    return($numscaf);
}


sub getNumberOfFragsInStore ($$) {
    my $wrk = shift @_;
    my $asm = shift @_;
    my $bin = getBinDirectory();

    return(0) if (! -e "$wrk/$asm.gkpStore/frg");

    open(F, "$bin/gatekeeper -L $wrk/$asm.gkpStore 2> /dev/null |") or caFailure("Failed to run gatekeeper to get the number of frags in the store.");
    $_ = <F>;    chomp $_;
    close(F);

    $numFrags = $1 if (m/^Last frag in store is iid = (\d+)$/);
    caFailure("No frags in the store?\n") if ($numFrags == 0);
    return($numFrags);
}


#  Decide if we have the CA meryl or the Mighty one.
#
sub merylVersion () {
    my $bin = getBinDirectory();
    my $ver = "unknown";

    open(F, "$bin/meryl -V |");
    while (<F>) {
        $ver = "CA"     if (m/CA/);
        $ver = "Mighty" if (m/Mighty/);
    }
    close(F);
    return($ver);
}



sub backupFragStore ($) {
    my $backupName = shift @_;

    return if (getGlobal("doBackupFragStore") == 0);

    if (-e "$wrk/$asm.gkpStore/frg.$backupName") {

        print STDERR "Found a backup for $backupName!  Restoring!\n";

        unlink "$wrk/$asm.gkpStore/frg";
        if (system("cp -p $wrk/$asm.gkpStore/frg.$backupName $wrk/$asm.gkpStore/frg")) {
            unlink "$wrk/$asm.gkpStore/frg";
            caFailure("Failed to restore gkpStore from backup.\n");
        }
    }
    if (! -e "$wrk/$asm.gkpStore/frg.$backupName") {

        print STDERR "Backing up the gkpStore to $backupName.\n";

        if (system("cp -p $wrk/$asm.gkpStore/frg $wrk/$asm.gkpStore/frg.$backupName")) {
            unlink "$wrk/$asm.gkpStore/frg.$backupName";
            caFailure("Failed to backup gkpStore.\n");
        }
    }
}



sub stopAfter ($) {
    my $stopAfter = shift @_;
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
    my $cmd;

    open(F, "> $script") or caFailure("Failed to open '$script' for writing\n");
    print F "#!/bin/sh\n";
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

    print F "hostname\n";
    print F "echo \$bin\n";

    print F "/usr/bin/env perl \$bin/runCA $commandLineOptions\n";
    close(F);

    my $sge       = getGlobal("sge");
    my $sgeScript = getGlobal("sgeScript");

    if (!defined($waitTag) || ($waitTag eq "")) {
        $cmd = "qsub $sge $sgeScript -cwd -r y -N runCA_${asm} -j y -o $output $script";
    } else {
        $cmd = "qsub $sge $sgeScript -cwd -r y -N runCA_${asm} -j y -o $output -hold_jid \"$waitTag\" $script";
    }

    system("chmod +x $script");
    system($cmd) and caFailure("Failed to submit script.\n");

    exit(0);
}


sub caFailure ($) {
    my  $msg = shift @_;
    localFailure($msg);
    die $msg;
}



#  Create an empty file.  Much faster than system("touch ...").
#
sub touch ($) {
    open(F, "> $_[0]") or caFailure("Failed to touch '$_[0]'\n");
    #print F "$wrk\n";
    #print F "process id: " . getppid() . "\n";
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

1;
