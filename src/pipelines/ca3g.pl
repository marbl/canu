#!/usr/bin/perl

use strict;

use FindBin;
use Cwd;

use File::Path qw(make_path remove_tree);

use Carp;
#use FileHandle;

use POSIX "ceil";
use POSIX "floor";
use POSIX "sys_wait_h";  #  waitpid()

use lib "$FindBin::RealBin";
use lib "$FindBin::RealBin/ca3g/lib/perl5";

use ca3g::Defaults;
use ca3g::Execution;

use ca3g::Gatekeeper;
use ca3g::Meryl;
use ca3g::OverlapInCore;
use ca3g::OverlapStore;
use ca3g::OverlapErrorAdjustment;
use ca3g::Unitig;


my $bin = undef;  #  Path to binaries, set once in main.
my $cmd = undef;  #  Temporary string passed to system().
my $wrk = undef;  #  Path to our assembly directory.
my $asm = undef;  #  Name of our assembly.

my %global;       #  Global parameters
my %synops;       #  Global parameters - description
my %synnam;       #  Global parameters - case sensitive name

my $specLog            = "";

my @specFiles;
my @specOpts;
my @inputFiles;



#require "./runCA-jobSubmission.pl";
#require "./runCA-defaults.pl";
#require "./runCA-gatekeeper.pl";






#  Initialize our defaults.

setDefaults();

#  Check for the presence of a -options switch BEFORE we do any work.
#  This lets us print the default values of options.

foreach my $arg (@ARGV) {
    if ($arg eq "-options") {
        setGlobal("options", 1);
        printHelp($bin);
    }
}

#  At some pain, we stash the original options for later use.  We need
#  to use these when we resubmit ourself to the grid.  We can't simply dump
#  all of @ARGV into here, because we need to fix up relative paths.

while (scalar(@ARGV)) {
    my $arg = shift @ARGV;

    if      ($arg =~ m/^-d/) {
        $wrk = shift @ARGV;
        $wrk = "$ENV{'PWD'}/$wrk" if ($wrk !~ m!^/!);
        addCommandLineOption("-d \"$wrk\"");

    } elsif ($arg eq "-p") {
        $asm = shift @ARGV;
        addCommandLineOption("-p \"$asm\"");

    } elsif ($arg eq "-s") {
        push @specFiles, shift @ARGV;

    } elsif ($arg eq "-version") {
        setGlobal("version", 1);

    } elsif ($arg eq "-options") {
        #  Do nothing.  Handled above, but we still need to process it here.
        #setGlobal("options", 1);

    } elsif (-e $arg) {
        $arg = "$ENV{'PWD'}/$arg" if ($arg !~ m!^/!);
        push @inputFiles, $arg;
        addCommandLineOption("\"$arg\"");

    } elsif ($arg =~ m/=/) {
        push @specOpts, $arg;
        addCommandLineOption("\"$arg\"");

    } else {
        setGlobal("help",
                  getGlobal("help") . "File not found or invalid command line option '$arg'\n");
    }
}


setGlobal("help", getGlobal("help") . "ERROR:  Assembly name prefix not supplied with -p.\n") if (!defined($asm));
setGlobal("help", getGlobal("help") . "ERROR:  Directory not supplied with -d.\n")            if (!defined($wrk));


$bin = getBinDirectory();

@inputFiles = setParametersFromFile("$bin/spec/runCA.default.specFile", @inputFiles)   if (-e "$bin/spec/runCA.default.specFile");
@inputFiles = setParametersFromFile("$ENV{'HOME'}/.runCA",              @inputFiles)   if (-e "$ENV{'HOME'}/.runCA");


#  For each of the specfiles on the command line, find the actual file and make it an absolute path.
#  These can be in the current directory (e.g., 'my.spec'), or in the installed directory ('$bin/spec').

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

    addCommandLineOption("-s \"$specFile\"");
}

#  For each of the spec files, parse it, setting parameters and remembering any input files discovered.

foreach my $specFile (@specFiles) {
    @inputFiles = setParametersFromFile($specFile, @inputFiles);
}

#  Set parameters from the command line.

setParametersFromCommandLine(@specOpts);

#  Finish setting parameters.

setParameters($bin);

#  If anything complained, global{'help'} will be defined, and we'll print help (and the error) and
#  stop.

printHelp($bin);

#  Fail immediately if we run the script on the grid, and the gkpStore
#  directory doesn't exist and we have no input files.  Without this
#  check we'd fail only after being scheduled on the grid.

if ((! -d "$wrk/$asm.gkpStore") &&
    (scalar(@inputFiles) == 0)) {
    caFailure("no input files specified, and store not already created", undef);
}

#  Check that we were supplied a work directory, and that it
#  exists, or we can create it.

caFailure("no run directory (-d option) specified", undef)  if (!defined($wrk));

make_path("$wrk")             if (! -d "$wrk/runCA-logs");
make_path("$wrk/runCA-logs")  if (! -d "$wrk/runCA-logs");

caFailure("run directory (-d option) '$wrk' doesn't exist and couldn't be created", undef)  if (! -d $wrk);

#  This environment variable tells the binaries to log their execution in runCA-logs/

$ENV{'AS_RUNCA_DIRECTORY'} = $wrk;

#  Report the parameters used.

writeLog($wrk);

#  Submit ourself for grid execution?  If not grid enabled, or already running on the grid, this
#  call just returns.  The arg MUST be undef.

submitScript($wrk, $asm, undef);

#  Begin

gatekeeper($wrk, $asm, @inputFiles);

meryl($wrk, $asm);

overlapConfigure($wrk, $asm, "normal");
#overlap($wrk, $asm, "normal");
overlapCheck($wrk, $asm, "normal", 0);
overlapCheck($wrk, $asm, "normal", 1);

createOverlapStore($wrk, $asm, "sequential");
#createOverlapStore($wrk, $asm, "parallel");

#if (0) {
#    initialTrim($wrk, $asm);
#    finalTrim($wrk, $asm);
#    chimeraDetection($wrk, $asm);
#}

#readErrorDetection($wrk, $asm);

overlapErrorAdjustment($wrk, $asm);

unitig($wrk, $asm);

#utgcns();
#make_consensus();
#pbdagcon();
#falcon_sense();

#output();

exit(0);




