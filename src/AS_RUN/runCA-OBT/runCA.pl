

my $specFile = undef;
my @specOpts;
my @fragFiles;

my @cgbFiles;
my $cgiFile;
my $scaffoldDir;

setDefaults();

#  At some pain, we stash the original options for later use.  We need
#  to use these when we resubmit ourself to SGE.
#
#  We can't simply dump all of @ARGV into here, because we need to
#  fix up relative paths.
#
$commandLineOptions = "";

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
        $specFile = shift @ARGV;
        $commandLineOptions .= " -s \"$specFile\"";

    } elsif ($arg eq "-version") {
        setGlobal("version", 1);

    } elsif ($arg eq "-options") {
        setGlobal("options", 1);

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
        setGlobal("help", 1);
    }
}

print STDERR "commandLineOptions=$commandLineOptions\n";

setGlobal("help", 1) if (!defined($asm));
setGlobal("help", 1) if (!defined($wrk));

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

#setup closure stuff
setupFilesForClosure();

#  If not already on the grid, see if we should be on the grid.
#  N.B. the arg MUST BE undef.
#
submitScript(undef) if (!runningOnGrid());

#  Begin

preoverlap(@fragFiles);
overlapTrim();
createOverlapJobs("normal");
checkOverlap("normal");
createOverlapStore();
overlapCorrection();
@cgbFiles = unitigger(@cgbFiles);
postUnitiggerConsensus(@cgbFiles);
scaffolder($cgiFile);
postScaffolderConsensus($scaffoldDir);
terminate($scaffoldDir);
cleaner();
toggler();

exit(0);
