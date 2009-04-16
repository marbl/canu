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

    #####  General Configuration Options (aka miscellany)

    $global{"pathMap"}                     = undef;
    $synops{"pathMap"}                     = "File with a hostname to binary directory map";

    #####  Error Rates

    $global{"ovlErrorRate"}                = 0.06;
    $synops{"ovlErrorRate"}                = "Overlaps above this error rate are not computed";

    $global{"utgErrorRate"}                = 0.015;
    $synops{"utgErrorRate"}                = "Overlaps above this error rate are not used to construct unitigs";

    $global{"cnsErrorRate"}                = 0.06;
    $synops{"cnsErrorRate"}                = "Consensus expects alignments at about this error rate";

    $global{"cgwErrorRate"}                = 0.10;
    $synops{"cgwErrorRate"}                = "Unitigs/Contigs are not merged if they align above this error rate";

    #####  Stopping conditions

    $global{"stopAfter"}                   = undef;
    $synops{"stopAfter"}                   = "Tell runCA when to halt execution";

    #####  Sun Grid Engine

    $global{"useGrid"}                     = 0;
    $synops{"useGrid"}                     = "Enable SGE globally";

    $global{"scriptOnGrid"}                = 0;
    $synops{"scriptOnGrid"}                = "Enable SGE for runCA (and unitigger, scaffolder, other sequential phases)";

    $global{"ovlOnGrid"}                   = 1;
    $synops{"ovlOnGrid"}                   = "Enable SGE for overlap computations";

    $global{"frgCorrOnGrid"}               = 0;
    $synops{"frgCorrOnGrid"}               = "Enable SGE for the fragment error correction";

    $global{"ovlCorrOnGrid"}               = 0;
    $synops{"ovlCorrOnGrid"}               = "Enable SGE for the overlap error correction";

    $global{"cnsOnGrid"}                   = 1;
    $synops{"cnsOnGrid"}                   = "Enable SGE for consensus";

    $global{"maxGridJobSize"}              = undef;
    $synops{"maxGridJobSize"}              = "";

    $global{"sge"}                         = undef;
    $synops{"sge"}                         = "SGE options applied to all SGE jobs";

    $global{"sgeScript"}                   = undef;
    $synops{"sgeScript"}                   = "SGE options applied to runCA jobs (and unitigger, scaffolder, other sequential phases)";

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

    #####  Vector Trimming

    $global{"vectorIntersect"}             = undef;
    $synops{"vectorIntersect"}             = "File of vector clear ranges";

    $global{"vectorTrimmer"}               = "ca";
    $synops{"vectorTrimmer"}               = "Use the CA default vector trimmer, or figaro";

    $global{"figaroFlags"}                 = "-T 30 -M 100 -E 500 -V f";
    $synops{"figaroFlags"}                 = "Options to the figaro vector trimmer";

    #####  Overlap Based Trimming

    $global{"perfectTrimming"}             = undef;  #  SECRET!
    $synops{"perfectTrimming"}             = undef;  #  SECRET!

    $global{"doOverlapTrimming"}           = 1;
    $synops{"doOverlapTrimming"}           = "Enable the Overlap Based Trimming module";

    #####  Overlapper

    $global{"obtOverlapper"}               = "ovl";
    $synops{"obtOverlapper"}               = "Which overlap algorithm to use for OBT overlaps";

    $global{"ovlOverlapper"}               = "ovl";
    $synops{"ovlOverlapper"}               = "Which overlap algorithm to use for OVL (unitigger) overlaps";

    $global{"ovlStoreMemory"}              = 1024;
    $synops{"ovlStoreMemory"}              = "How much memory (MB) to use when constructing overlap stores";

    $global{"ovlThreads"}                  = 2;
    $synops{"ovlThreads"}                  = "Number of threads to use when computing overlaps";

    $global{"ovlStart"}                    = 1;
    $synops{"ovlStart"}                    = "Starting fragment for overlaps (EXPERT!)";

    $global{"ovlHashBlockSize"}            = 200000;
    $synops{"ovlHashBlockSize"}            = "Number of fragments to load into the in-core overlap hash table";

    $global{"ovlRefBlockSize"}             = 2000000;
    $synops{"ovlRefBlockSize"}             = "Number of fragments to search against the hash table per batch";

    $global{"ovlMemory"}                   = "2GB";
    $synops{"ovlMemory"}                   = "Amount of memory to use for overlaps";

    $global{"ovlMerSize"}                  = 22;
    $synops{"ovlMerSize"}                  = "K-mer size for seeds in overlaps";

    $global{"ovlMerThreshold"}             = "auto";
    $synops{"ovlMerThreshold"}             = "K-mer frequency threshold; mers more frequent than this are ignored";

    $global{"obtMerSize"}                  = 22;
    $synops{"obtMerSize"}                  = "K-mer size";

    $global{"obtMerThreshold"}             = "auto";
    $synops{"obtMerThreshold"}             = "K-mer frequency threshold; mers more frequent than this are ignored";

    $global{"merCompression"}              = 1;
    $synops{"merCompression"}              = "K-mer size";

    $global{"merOverlapperThreads"}        = 2;
    $synops{"merOverlapperThreads"}        = "Number of threads to use in the mer overlapper";

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

    $global{"merOverlapperThreads"}        = 2;
    $synops{"merOverlapperThreads"}        = "Number of threads to use for both mer overlapper seed finding and extension jobs";

    $global{"umdOverlapperFlags"}          = "-use-uncleaned-reads -trim-error-rate 0.03 -max-minimizer-cutoff 150";
    $synops{"umdOverlapperFlags"}          = "Options for the UMD overlapper";

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

    #####  Unitigger & BOG Options

    $global{"unitigger"}                   = "utg";
    $synops{"unitigger"}                   = "Which unitig algorithm to use; utg or bog (Best Overlap Graph)";

    $global{"utgGenomeSize"}               = undef;
    $synops{"utgGenomeSize"}               = "An estimate of the size of the genome; decides if unitigs are unique or repeats";

    $global{"utgBubblePopping"}            = 1;
    $synops{"utgBubblePopping"}            = "Smooth polymorphic regions";

    $global{"utgRecalibrateGAR"}           = 1;
    $synops{"utgRecalibrateGAR"}           = "Use an experimental algorithm to decide unique/repeat";

    $global{"bogPromiscuous"}              = 0;
    $synops{"bogPromiscuous"}              = "EXPERT!";

    $global{"bogEjectUnhappyContain"}      = 0;
    $synops{"bogEjectUnhappyContain"}      = "EXPERT!";

    $global{"bogBadMateDepth"}             = 7;
    $synops{"bogBadMateDepth"}             = "EXPERT!";

    #####  Scaffolder Options

    $global{"cgwOutputIntermediate"}       = 0;
    $synops{"cgwOutputIntermediate"}       = "Output .cgw files for intermediate scaffolding (advanced)";

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

    $global{"computeInsertSize"}           = 0;
    $synops{"computeInsertSize"}           = "Compute a scratch scaffolding to estimate insert sizes";

    $global{"cgwDistanceSampleSize"}       = 100;
    $synops{"cgwDistanceSampleSize"}       = "Require N mates to reestimate insert sizes";

    $global{"doResolveSurrogates"}         = 1;
    $synops{"doResolveSurrogates"}         = "Place fragments in surrogates in the final assembly";

    $global{"doExtendClearRanges"}         = 2;
    $synops{"doExtendClearRanges"}         = "Enable the clear range extension heuristic";

    $global{"extendClearRangesStepSize"}   = undef;
    $synops{"extendClearRangesStepSize"}   = "Batch N scaffolds per ECR run";

    #####  Consensus Options

    $global{"cnsPartitions"}               = 128;
    $synops{"cnsPartitions"}               = "Partition consensus into N jobs";

    $global{"cnsMinFrags"}                 = 75000;
    $synops{"cnsMinFrags"}                 = "Don't make a consensus partition with fewer than N fragments";

    $global{"cnsConcurrency"}              = 2;
    $synops{"cnsConcurrency"}              = "If not SGE, number of consensus jobs to run at the same time";

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

    #####  Ugly, command line options passed to printHelp()

    $global{"help"}                        = 0;
    $synops{"help"}                        = undef;

    $global{"version"}                     = 0;
    $synops{"version"}                     = undef;

    $global{"help"}                        = 0;
    $synops{"help"}                        = undef;

    $global{"specFile"}                    = undef;
    $synops{"specFile"}                    = undef;
    
    #### Closure Options

    $global{"closureEdges"}               = undef;
    $synops{"closureEdges"}               = undef;

    $global{"closureOverlaps"}             = 0;
    $synops{"closureOverlaps"}             = undef;
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

            if (m/\s*(\w*)\s*=([^#]*)#*.*$/) {
                my ($var, $val) = ($1, $2);
                print STDERR $_,"\n"; # echo the spec file
                $var =~ s/^\s+//; $var =~ s/\s+$//;
                $val =~ s/^\s+//; $val =~ s/\s+$//;
                undef $val if ($val eq "undef");
                setGlobal($var, $val);
            } else {
                my $xx = $_;
                $xx = "$ENV{'PWD'}/$xx" if ($xx !~ m!^/!);
                if (-e $xx) {
                    push @fragFiles, $xx;
                } else {
                    print STDERR "WARNING!  Invalid specFile line '$_'\n";
                }
            }
        }
        close(F);
    }

    return(@fragFiles);
}


sub setParametersFromCommandLine(@) {
    my @specOpts = @_;

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
}


sub setParameters () {

    #  Fiddle with filenames to make them absolute paths.
    #
    makeAbsolute("vectorIntersect");
    makeAbsolute("pathMap");

    #  Adjust case on some of them
    #
    fixCase("obtOverlapper");
    fixCase("ovlOverlapper");
    fixCase("unitigger");
    fixCase("vectorTrimmer");
    #fixCase("stopAfter");
    fixCase("consensus");
    fixCase("cleanup");

    if ((getGlobal("obtOverlapper") ne "mer") && (getGlobal("obtOverlapper") ne "ovl")) {
        caFailure("Invalid obtOverlapper specified (" . getGlobal("obtOverlapper") . "); must be 'mer' or 'ovl'\n");
    }
    if ((getGlobal("ovlOverlapper") ne "mer") && (getGlobal("ovlOverlapper") ne "ovl")) {
        caFailure("Invalid ovlOverlapper specified (" . getGlobal("ovlOverlapper") . "); must be 'mer' or 'ovl'\n");
    }
    if ((getGlobal("unitigger") ne "utg") && (getGlobal("unitigger") ne "bog")) {
        caFailure("Invalid unitigger specified (" . getGlobal("unitigger") . "); must be 'utg' or 'bog'\n");
    }
    if ((getGlobal("vectorTrimmer") ne "ca") && (getGlobal("vectorTrimmer") ne "figaro")) {
        caFailure("Invalid vectorTrimmer specified (" . getGlobal("vectorTrimmer") . "); must be 'ca' or 'figaro'\n");
    }
    if ((getGlobal("consensus") ne "cns") && (getGlobal("consensus") ne "seqan")) {
        caFailure("Invalid consensus specified (" . getGlobal("consensus") . "); must be 'cns' or 'seqan'\n");
    }
    if ((getGlobal("cleanup") ne "none") &&
        (getGlobal("cleanup") ne "light") &&
        (getGlobal("cleanup") ne "heavy") &&
        (getGlobal("cleanup") ne "aggressive")) {
        caFailure("Invalid cleaup specified (" . getGlobal("cleanup") . "); must be 'none', 'light', 'heavy' or 'aggressive'\n");
    }

    #  PIck a nice looking set of binaries, and check them.
    #
    {
        my $bin = getBinDirectory();

        caFailure("Can't find 'gatekeeper' program in $bin.  Possibly incomplete installation.\n") if (! -x "$bin/gatekeeper");
        caFailure("Can't find 'meryl' program in $bin.  Possibly incomplete installation.\n")      if (! -x "$bin/meryl");
        caFailure("Can't find 'overlap' program in $bin.  Possibly incomplete installation.\n")    if (! -x "$bin/overlap");
        caFailure("Can't find 'unitigger' program in $bin.  Possibly incomplete installation.\n")  if (! -x "$bin/unitigger");
        caFailure("Can't find 'cgw' program in $bin.  Possibly incomplete installation.\n")        if (! -x "$bin/cgw");
        caFailure("Can't find 'consensus' program in $bin.  Possibly incomplete installation.\n")  if (! -x "$bin/consensus");
        caFailure("Can't find 'terminator' program in $bin.  Possibly incomplete installation.\n") if (! -x "$bin/terminator");

        if ((getGlobal("obtOverlapper") eq "mer") || (getGlobal("ovlOverlapper") eq "mer")) {
            caFailure("Can't find 'overmerry' program in $bin.  Possibly incomplete installation.\n") if (! -x "$bin/overmerry");
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
        caFailure("ovlErrorRate is $ovlER, this MUST be between 0.00 and 0.25.\n");
    }
    if (($utgER < 0.0) || (0.25 < $utgER)) {
        caFailure("utgErrorRate is $utgER, this MUST be between 0.00 and 0.25.\n");
    }
    if (($cgwER < 0.0) || (0.25 < $cgwER)) {
        caFailure("cgwErrorRate is $cgwER, this MUST be between 0.00 and 0.25.\n");
    }
    if (($cnsER < 0.0) || (0.25 < $cnsER)) {
        caFailure("cnsErrorRate is $cnsER, this MUST be between 0.00 and 0.25.\n");
    }
    if ($utgER > $ovlER) {
        caFailure("utgErrorRate is $utgER, this MUST be <= ovlErrorRate ($ovlER)\n");
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

    if (getGlobal("version")) {
        my $bin = getBinDirectory();

        system("$bin/gatekeeper   --version");
        system("$bin/overlap      --version");
        system("$bin/unitigger    --version");
        system("$bin/buildUnitigs --version");
        system("$bin/cgw          --version");
        system("$bin/consensus    --version");
        system("$bin/terminator   --version");

        exit(0);
    }

    if (getGlobal("help")) {
        foreach my $k (sort keys %global) {
            my $o = substr("$k                             ", 0, 35);
            my $d = substr(getGlobal($k) . "               ", 0, 20);
            my $u = $synops{$k};

            if (!defined(getGlobal($k))) {
                $d = substr("<unset>" . "               ", 0, 16);
            }

            print "$o$d($u)\n";
        }
        exit(0);
    }

    if (getGlobal("help")) {
        print "usage: runCA -d <dir> -p <prefix> [options] <frg>\n";
        print "  -d <dir>          Use <dir> as the working directory.  Required\n";
        print "  -p <prefix>       Use <prefix> as the output prefix.  Required\n";
        print "\n";
        print "  -s <specFile>     Read options from the specifications file <specfile>.\n";
        print "                      <specfile> can also be one of the following key words:\n";
        print "                      [no]OBT - run with[out] OBT\n";
        print "                      noVec   - run with OBT but without Vector\n";
        print "\n";
        print "  -version          Version information\n";
        print "  -help             Describe specFile options, and show default values\n";
        print "\n";
        print "  <frg>             CA formatted fragment file\n";
        print "\n";
        print "Complete documentation at http://wgs-assembler.sourceforge.net/\n";
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

    $ENV{'AS_RUNCA_DIRECTORY'} = $wrk;

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

    $dir = "$wrk/$dir" if (-d "$wrk/$dir");

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

    open(F, "$bin/gatekeeper -lastfragiid $wrk/$asm.gkpStore 2> /dev/null |") or caFailure("Failed to run gatekeeper to get the number of frags in the store.");
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



sub removeFragStoreBackup ($) {
    my $backupName = shift @_;

    unlink "$wrk/$asm.gkpStore/frg.$backupName";
}

sub restoreFragStoreBackup ($) {
    my $backupName = shift @_;

    if (-e "$wrk/$asm.gkpStore/frg.$backupName") {
        print STDERR "Restoring the gkpStore backup from $backupName.\n";
        unlink "$wrk/$asm.gkpStore/frg.FAILED";
        rename "$wrk/$asm.gkpStore/frg", "$wrk/$asm.gkpStore/frg.$backupName.FAILED";
        rename "$wrk/$asm.gkpStore/frg.$backupName", "$wrk/$asm.gkpStore/frg";
    }
}

sub backupFragStore ($) {
    my $backupName = shift @_;

    return if (-e "$wrk/$asm.gkpStore/frg.$backupName");

    if (system("cp -p $wrk/$asm.gkpStore/frg $wrk/$asm.gkpStore/frg.$backupName")) {
        unlink "$wrk/$asm.gkpStore/frg.$backupName";
        caFailure("Failed to backup gkpStore.\n");
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

    system("chmod +x $script");

    my $sge         = getGlobal("sge");
    my $sgeScript   = getGlobal("sgeScript");
    my $sgePropHold = getGlobal("sgePropagateHold");

    $waitTag = "-hold_jid \"$waitTag\"" if (defined($waitTag));

    my $qcmd = "qsub $sge $sgeScript -cwd -N \"runCA_${asm}\" -j y -o $output $waitTag $script";
    print STDERR "$qcmd\n";
    system($qcmd) and caFailure("Failed to submit script.\n");

    if (defined($sgePropHold)) {
        my $acmd = "qalter -hold_jid \"runCA_${asm}\" \"$sgePropHold\"";
        print STDERR "$acmd\n";
        system($acmd) and print STDERR "WARNING: Failed to reset hold_jid trigger on '$sgePropHold'.\n";
    }

    exit(0);
}


sub caFailure ($) {
    my  $msg = shift @_;
    localFailure($msg);
    die $msg;
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
    open(F, "> $_[0]") or caFailure("Failed to touch '$_[0]'\n");
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

sub setupFilesForClosure() {
    makeAbsolute("closureEdges");

    my $closureEdges = getGlobal("closureEdges");

    if (defined($closureEdges)) {
       #default to only allowing overlaps between closure reads and other reads. No overlaps within the set of closure reads  
       setGlobal("closureOverlaps", 2);
       if (-e "$wrk/closureEdges") {
          return;
       }
       system("cp $closureEdges $wrk/$asm.gkpStore.closureEdges");
    }
}
1;
