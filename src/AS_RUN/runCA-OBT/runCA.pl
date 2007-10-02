

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
    $invocation .= " $a";
    $commandLineOptions .= " \"$a\" ";
}

while (scalar(@ARGV)) {
    my $arg = shift @ARGV;

    if ($arg =~ m/^-alias/ && $JCVI == 1 ) {
        setGlobal("alias",shift @ARGV);
    } elsif ($arg =~ m/^-maxCopy/ && $JCVI == 1 ) {
        setGlobal("copy",getGlobal("maxCopy"));
    } elsif ($arg =~ m/^-minCopy/ && $JCVI == 1 ) {
        setGlobal("copy",getGlobal("minCopy"));
    } elsif ($arg =~ m/^-medCopy/ && $JCVI == 1 ) {
        setGlobal("copy",getGlobal("medCopy"));
    } elsif ($arg =~ m/^-noCopy/ && $JCVI == 1 ) {
        setGlobal("copy",getGlobal("noCopy"));
    } elsif ($arg =~ m/^-(no)?notify/ && $JCVI == 1 ) {
        setGlobal("notify",1) if ( $1 ne 'no' );
        setGlobal("notify",0) if ( $1 eq 'no' );
    } elsif ($arg =~ m/^-test/ && $JCVI == 1 ) {    
        setGlobal("test",1);
    } elsif ($arg =~ m/^-e/) {
        setGlobal("utgErrorRate",shift @ARGV);
    } elsif ($arg =~ m/^-g/) {
        setGlobal("utgGenomeSize",shift @ARGV);
    } elsif ($arg =~ m/^-j/) {
        setGlobal("astatLowBound",shift @ARGV);
    } elsif ($arg =~ m/^-k/) {    
        setGlobal("astatHighBound",shift @ARGV);
    } elsif ($arg =~ m/^-(no)?ubs/) {
        setGlobal("utgBubblePopping",1) if ( $1 ne 'no' );
        setGlobal("utgBubblePopping",0) if ( $1 eq 'no' );
    } elsif ($arg =~ m/^-d/) {
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
        print STDERR "$0: Unknown argument '$arg'\n";
    }
}

setParameters($specFile, @specOpts);
printHelp();

my $retVal;
my @steps = 
(
	'Pre-overlap',			'$retVal = preoverlap(@fragFiles)',
	'OverlapTrim',			'$retVal = overlapTrim()',
	'CreateOverlapJobs',		'$retVal = createOverlapJobs("normal")',
	'CheckOverlap',		'$retVal = checkOverlap("normal")',
	'CreateOverlapStore',		'$retVal = createOverlapStore()',
	'CreateFragmentCorrections',	'$retVal = createFragmentCorrectionJobs()',
	'MergeFragmentCorrection',	'$retVal = mergeFragmentCorrection()',
	'CreateOverlapCorrections',	'$retVal = createOverlapCorrectionJobs()',
	'ApplyOverlapCorrection',	'$retVal = applyOverlapCorrection()',
	'Unitigger',			'($retVal,@cgbFiles) = unitigger(@cgbFiles)',
	'PostUnitiggerConsensus',	'$retVal = postUnitiggerConsensus(@cgbFiles)',
	'Scaffolder',			'$retVal = scaffolder($cgiFile)',
	'PostScaffolderConsensus',	'$retVal = postScaffolderConsensus($scaffoldDir)',
	'Terminator',			'$retVal = terminate($scaffoldDir)'
);

if ( $JCVI == 1 && !runningOnGrid()) {
    $request_id = asdbInit();
    print "Your Assembly Console request id is: $request_id\n";
}

checkDirectories();

print "Your work directory is: '$wrk'\n";

createInvocationScript() && init_prop_file( undef, $asm,(scalar @steps)/2)
	if ( $JCVI == 1 && !runningOnGrid());

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
my $finished = 0;

for( my $index = 0 ; $index < scalar @steps ; $index+=2) {
	my $stepName = $steps[$index];
	my $stepCmd = $steps[$index+1];
	if ( $JCVI == 1 && ! -e "$wrk/log/$stepName.started") {
	    touch("$wrk/log/$stepName.started");
	    start($stepName);
	}
	print "Executing '$stepName'\n";
	eval($stepCmd);
	print "Finished executing '$stepName'\nReturn value: $retVal\n";
	if ($retVal == -1) {
	    failure($stepName);
	} elsif ( $JCVI == 1 && !-e "$wrk/log/$stepName.finished") {
	    touch("$wrk/log/$stepName.finished");
	    finish($stepName);
	    $finished = 1 if ( $index + 2 >= scalar @steps );
	}
}

copy() if ( $JCVI == 1 && $finished == 1);

exit(0);

