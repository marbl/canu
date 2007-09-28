

my $specFile = undef;
my @specOpts;
my @fragFiles;

my $isContinuation = 0;
my @cgbFiles;
my $cgiFile;
my $scaffoldDir;

setDefaults();

#  Stash the original options, quoted, for later use.
foreach my $a (@ARGV) {
    $commandLineOptions .= " \"$a\" ";
}

while (scalar(@ARGV)) {
    my $arg = shift @ARGV;

    if      ($arg =~ m/^-d/) {
        $wrk = shift @ARGV;
        $wrk = "$ENV{'PWD'}/$wrk" if ($wrk !~ m!^/!);
    } elsif ($arg =~ m/^-p/) {
        $asm = shift @ARGV;
    } elsif ($arg =~ m/^-s/) {
        $specFile = shift @ARGV;
    } elsif ($arg =~ m/^-h/) {
        setGlobal("help", 1);
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
        print STDERR "$0: Unknown argument '$arg'\n";
    }
}

setParameters($specFile, @specOpts);
printHelp();
checkDirectories();

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

    #  If given a scaffold directory, tell unitigeger, consensus and
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

#  Begin

#  If not already on the grid, see if we should be on the grid.
#
submitScript("") if (!runningOnGrid());

preoverlap(@fragFiles);

overlapTrim();

createOverlapJobs("normal");
checkOverlap("normal");

createOverlapStore();

createFragmentCorrectionJobs();
mergeFragmentCorrection();
createOverlapCorrectionJobs();
applyOverlapCorrection();

@cgbFiles = unitigger(@cgbFiles);

postUnitiggerConsensus(@cgbFiles);

scaffolder($cgiFile);

postScaffolderConsensus($scaffoldDir);

terminate($scaffoldDir);

exit(0);

