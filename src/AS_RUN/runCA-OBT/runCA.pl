

my $specFile = undef;
my @specOpts;
my @fragFiles;

my $isContinuation = 0;
my @cgbFiles;
my $cgiFile;
my $scaffoldDir;

my @steps = (
             'Pre-overlap',               'preoverlap(@fragFiles)',
             'OverlapTrim',               'overlapTrim()',
             'CreateOverlapJobs',         'createOverlapJobs("normal")',
             'CheckOverlap',              'checkOverlap("normal")',
             'CreateOverlapStore',        'createOverlapStore()',
             'OverlapCorrection',         'overlapCorrection()',
             'Unitigger',                 '@cgbFiles = unitigger(@cgbFiles)',
             'PostUnitiggerConsensus',    'postUnitiggerConsensus(@cgbFiles)',
             'Scaffolder',                'scaffolder($cgiFile)',
             'PostScaffolderConsensus',   'postScaffolderConsensus($scaffoldDir)',
             'Terminator',                'terminate($scaffoldDir)',
             'Cleaner',                   'cleaner()',
             );


setDefaults();
localDefaults();

# echo invocation to aid debugging
my $oldSep = $,;
$, = ' ';
print STDERR $0,@ARGV,"\n\n";
$, = $oldSep;

#  Stash the original options, quoted, for later use.  We need to use
#  these when we resubmit ourself to SGE.
foreach my $a (@ARGV) {
    $invocation .= " $a";
    $commandLineOptions .= " \"$a\" ";
}


while (scalar(@ARGV)) {
    my $arg = shift @ARGV;

    my $found = 0;
    ($found, @ARGV) = localOption($arg, @ARGV);

    if ($found == 1 ) {
    	next;
    }

    if (!defined($arg)) {
        last;
    }

    if      ($arg =~ m/^-d/) {
        $wrk = shift @ARGV;
        $wrk = "$ENV{'PWD'}/$wrk" if ($wrk !~ m!^/!);
    } elsif ($arg =~ m/^-p/) {
        $asm = shift @ARGV;
    } elsif ($arg =~ m/^-s/) {
        $specFile = shift @ARGV;
    } elsif ($arg =~ m/^-h/) {
        setGlobal("help", 1);
    } elsif ($arg =~ m/^-v/ or $arg =~ m/^-V/) {
        setGlobal("version", 1);
    } elsif ($arg =~ m/^-f/) {
        setGlobal("fields", 1);
    } elsif (($arg =~ /\.frg$|frg\.gz$|frg\.bz2$/) && (-e $arg)) {
        $arg = "$ENV{'PWD'}/$arg" if ($arg !~ m!^/!);
        push @fragFiles, $arg;
    } elsif (($arg =~ /\.sff$|sff\.gz$|sff\.bz2$/) && (-e $arg)) {
        $arg = "$ENV{'PWD'}/$arg" if ($arg !~ m!^/!);
        push @fragFiles, $arg;
    } elsif (($arg =~ /\.cgb$/) && (-e $arg)) {
        $isContinuation = 1;
        $arg = "$ENV{'PWD'}/$arg" if ($arg !~ m!^/!);
        push @cgbFiles, $arg;
    } elsif (($arg =~ /\.cgi$/) && (-e $arg)) {
        $isContinuation = 1;
        $cgiFile = $arg;
        $cgiFile = "$ENV{'PWD'}/$cgiFile" if ($cgiFile !~ m!^/!);
    } elsif (-d $arg) {
        $isContinuation = 1;
        $scaffoldDir  = $arg;
        $scaffoldDir  = "$ENV{'PWD'}/$scaffoldDir" if ($scaffoldDir !~ m!^/!);
    } elsif ($arg =~ m/=/) {
        push @specOpts, $arg;
    } else {
        die "$0: Unknown argument '$arg'\n";
    }
}

setGlobal("help",1) unless $asm;

@fragFiles = setParametersFromFile($specFile, @fragFiles);

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
    die "No frg files given on the command line.\nTry $0 -h for help.\n\n";
}

checkDirectories();

localSetup(scalar(@steps) / 2);

#  If this is a continuation, we don't want to do obt or fragment
#  error correction, or a bunch of other stuff.  We could surround
#  those steps below with if's, but the whole design of this script is
#  that each piece checks if it is done or not.  So, we disable those
#  pieces.
#
if ($isContinuation) {
    setGlobal("doOverlapTrimming", 0);
    setGlobal("doFragmentCorrection", 0);

    #  If given cgb files, we don't need to do anything more

    #  If given cgi files, we need to tell unitigger and consensus that
    #  we're done.
    if (defined($cgiFile)) {
        system("mkdir $wrk/4-unitigger") if (! -d "$wrk/4-unitigger");
        touch("$wrk/4-unitigger/unitigger.success");

        system("mkdir $wrk/5-consensus") if (! -d "$wrk/5-consensus");
        touch("$wrk/5-consensus/jobsCreated.success");
        touch ("$wrk/5-consensus/consensus.success");
    }

    #  If given a scaffold directory, tell unitigger, consensus and
    #  scaffolder that they are done.
    if (defined($scaffoldDir)) {
        system("mkdir $wrk/4-unitigger") if (! -d "$wrk/4-unitigger");
        touch("$wrk/4-unitigger/unitigger.success");

        system("mkdir $wrk/5-consensus") if (! -d "$wrk/5-consensus");
        touch("$wrk/5-consensus/jobsCreated.success");
        touch ("$wrk/5-consensus/consensus.success");

        system("mkdir $wrk/7-CGW") if (! -d "$wrk/7-CGW");
        touch ("$wrk/7-CGW/cgw.success");
    }
}


#  If not already on the grid, see if we should be on the grid.
#  N.B. the arg MUST BE undef.
#
submitScript(undef) if (!runningOnGrid());


#  Begin


while (scalar(@steps) > 0) {
    my $stepName = shift @steps;
    my $stepCmd  = shift @steps;

    localStart($stepName);
    eval($stepCmd);
    if ($@) {
        chomp $@;
        die "step $stepName failed with '$@'\n" ;
    }
    localFinish($stepName);
}

localFinalize();
exit(0);
