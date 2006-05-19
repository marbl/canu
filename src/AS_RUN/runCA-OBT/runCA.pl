#!/usr/bin/perl

use strict;
use Config;  #  for @signame
use FindBin;

use vars qw($bin $gin $wrk $asm);

use vars qw($numFrags);

#  Set some not reasonable defaults.
$wrk = undef;
$asm = "asm";

my %global;
my $specFile = undef;
my @specOpts;
my @fragFiles;

my $isContinuation = 0;
my @cgbFiles;
my $cgiFile;
my $scaffoldDir;

setDefaults();

while (scalar(@ARGV)) {
    my $arg = shift @ARGV;

    if      ($arg =~ m/^-d/) {
        $wrk = shift @ARGV;
        if ($wrk !~ m!^/!) {
            $wrk = "$ENV{'PWD'}/$wrk";
            print STDERR "WARNING:  Constructing full path to work directory: $wrk\n";
        }
    } elsif ($arg =~ m/^-p/) {
        $asm = shift @ARGV;
    } elsif ($arg =~ m/^-s/) {
        $specFile = shift @ARGV;
    } elsif ($arg =~ m/^-h/) {
        setGlobal("help", 1);
    } elsif (($arg =~ /\.frg$|frg\.gz$|frg\.bz2$/) && (-e $arg)) {
        push @fragFiles, $arg;
    } elsif (($arg =~ /\.cgb$/) && (-e $arg)) {
        $isContinuation = 1;
        $arg  = "$ENV{'PWD'}/$arg" if ($arg !~ m/^\//);
        push @cgbFiles, $arg;
    } elsif (($arg =~ /\.cgi$/) && (-e $arg)) {
        $isContinuation = 1;
        $cgiFile        = $arg;
        $cgiFile        = "$ENV{'PWD'}/$cgiFile" if ($cgiFile !~ m/^\//);
    } elsif (-d $arg) {
        $isContinuation = 1;
        $scaffoldDir  = $arg;
        $scaffoldDir  = "$ENV{'PWD'}/$scaffoldDir" if ($scaffoldDir !~ m/^\//);
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

    print STDERR "IS CONTINUE!\n";

    #  If given cgb files, we don't need to do anything more

    #  If given cgi files, we need to tell unitigger and consensus that
    #  we're done.
    if (defined($cgiFile)) {
        print STDERR "IS CONTINUE SCAFFOLDER!\n";

        system("mkdir $wrk/4-unitigger") if (! -d "$wrk/4-unitigger");
        touch("$wrk/4-unitigger/unitigger.success");

        system("mkdir $wrk/5-consensus") if (! -d "$wrk/5-consensus");
        touch("$wrk/5-consensus/jobsCreated.success");
        touch ("$wrk/5-consensus/consensus.success");
    }

    #  If given a scaffold directory, tell unitigeger, consensus and
    #  scaffolder that they are done.
    if (defined($scaffoldDir)) {
        print STDERR "IS CONTINUE CONSENSUS!\n";
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

createPostUnitiggerConsensusJobs(@cgbFiles);
checkPostUnitiggerConsensus(@cgbFiles);

scaffolder($cgiFile);

createConsensusJobs($scaffoldDir);

terminate($scaffoldDir);

exit(0);

