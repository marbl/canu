use strict;


#  Decide what host we are on, and the the bin directory
#  appropriately
#
sub setBinDirectory ($$) {
    my $host    = shift @_;
    my $mach    = shift @_;
    my $binRoot;
    my $binDir;

    #  See if we're being forced to a host/maching combination.
    #
    $host = `uname`     if (!defined($host) || ($host eq ""));
    $mach = `uname -m`  if (!defined($mach) || ($mach eq ""));
    chomp $host;
    chomp $mach;


    #  Test if we are in an installed wgs-assembler tree.  If so, set the
    #  root to the root of the assembler tree.
    #
    if (-e "$FindBin::Bin/gatekeeper") {
        $binRoot = "$FindBin::Bin";
        my @t = split '/', $binRoot;
        $#t -= 2;
        $binRoot = join '/', @t;
    }

    #  Let the user override this root, for no good reason.  We'll warn, though.
    #
    my $t = getGlobal("binRoot");
    if (defined($binRoot) && (defined($t)) && ($binRoot ne $t)) {
        print STDERR "WARNING: I appear to be installed in $binRoot, but you\n";
        print STDERR "         specified a bin directory of $t in the spec file.\n";
        $binRoot = $t;
    }

    if (!defined($binRoot)) {
        $binRoot = $t;
    }

    if (!defined($binRoot)) {
        print STDERR "ERROR: I'm not installed in the assembler tree, and you\n";
        print STDERR "       didn't specify a binDir in the spec file.\n";
        exit(1);
    }



    if      (($host eq "Linux") && ($mach eq "i686")) {
        #  Linux, on Intel 686
        $binDir       = "$binRoot/Linux/bin";
    } elsif (($host eq "Linux") && ($mach eq "x86_64")) {
        #  Linux on Opteron
        $binDir       = "$binRoot/Linux64/bin";
    } elsif (($host eq "FreeBSD") && ($mach eq "i386")) {
        #  FreeBSD on Intel
        $binDir       = "$binRoot/FreeBSD/bin";
    } elsif ($host eq "Darwin") {
        #  Darwin (UNTESTED)
        $binDir       = "$binRoot/Darwin/bin";
    } elsif (($mach eq "alpha") || ($mach eq "OSF1")) {
        #  OSF1 on alpha
        $binDir       = "$binRoot/OSF1/bin";
    } else {
        print STDERR "Unknown host '$host' and machine '$mach' combination.\n";
        exit(1);
    }

    if (! -d $binRoot) {
        die "ERROR: Failed to find the binary directory $binDir.\n";
    }

    return($binDir);
}


#  Return the second argument, unless the first argument is found in
#  %global, in which case return that.
#
sub getGlobal ($) {
    my $var = shift @_;
    die "ERROR: $var has no defined value!\n" if (!exists($global{$var}));
    return($global{$var});
}

sub setGlobal ($$) {
    my $var = shift @_;
    my $val = shift @_;
    die "ERROR: $var is not a valid option.\n" if (!exists($global{$var}));
    $global{$var} = $val;
}

sub setDefaults {
    $global{"binRoot"}                     = undef;
    $global{"cnsPartitions"}               = 128;
    $global{"cnsMinFrags"}                 = 75000;
    $global{"cnsConcurrency"}              = 2;
    $global{"cnsOnGrid"}                   = 1;
    $global{"delayInterleavedMerging"}     = 0;
    $global{"doBackupFragStore"}           = 1;
    $global{"doExtendClearRanges"}         = 2;
    $global{"doFragmentCorrection"}        = 1;
    $global{"doOverlapTrimming"}           = 1;
    $global{"doResolveSurrogates"}         = 1;
    $global{"doUpdateDistanceRecords"}     = 1;
    $global{"fakeUIDs"}                    = 0;
    $global{"frgCorrBatchSize"}            = 175000;
    $global{"frgCorrOnGrid"}               = 0;
    $global{"frgCorrThreads"}              = 2;
    $global{"frgCorrConcurrency"}          = 1;
    $global{"gridHost"}                    = "Linux";
    $global{"gridMachine"}                 = "i686";
    $global{"immutableFrags"}              = undef;
    $global{"localHost"}                   = undef;
    $global{"localMachine"}                = undef;
    $global{"ovlCorrBatchSize"}            = 175000;
    $global{"ovlCorrOnGrid"}               = 0;
    $global{"ovlCorrConcurrency"}          = 4;
    $global{"ovlHashBlockSize"}            = 40000;
    $global{"ovlMemory"}                   = "1GB";
    $global{"ovlRefBlockSize"}             = 2000000;
    $global{"ovlSortMemory"}               = 1000;
    $global{"ovlStoreMemory"}              = 512;
    $global{"ovlThreads"}                  = 2;
    $global{"ovlOnGrid"}                   = 1;
    $global{"processStats"}                = undef;
    $global{"scratch"}                     = "/scratch";
    $global{"stoneLevel"}                  = 2;
    $global{"uidServer"}                   = undef;
    $global{"updateDistanceType"}          = "pre";
    $global{"utgEdges"}                    = undef;
    $global{"utgErrorRate"}                = 15;
    $global{"utgFragments"}                = undef;
    $global{"utgBubblePopping"}            = 1;
    $global{"utgGenomeSize"}               = undef;
    $global{"useGrid"}                     = 0;
    $global{"useBogUnitig"}                = 0;
    $global{"vectorIntersect"}             = undef;
}

sub setParameters {
    my $specFile = shift @_;
    my @specOpts = @_;

    if (defined($specFile)) {
        open(F, "< $specFile") or die "Failed to open '$specFile'\n";
        while (<F>) {
            chomp;
            next if (m/^\s*\#/);
            next if (m/^\s*$/);
            if (m/\s*(\w*)\s*=(.*)/) {
                my ($var, $val) = ($1, $2);
                $var =~ s/^\s+//; $var =~ s/\s+$//;
                $val =~ s/^\s+//; $val =~ s/\s+$//;
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

    #  If we have been given a configuration in the spec file, use that instead.
    #
    my $binhost = getGlobal("localHost");
    my $binmach = getGlobal("localMachine");
    my $ginhost = getGlobal("gridHost");
    my $ginmach = getGlobal("gridMachine");

    $bin        = setBinDirectory($binhost, $binmach);
    $gin        = setBinDirectory($binhost, $binmach)  if (getGlobal("useGrid") == 0);  #  bin for the Grid
    $gin        = setBinDirectory($ginhost, $ginmach)  if (getGlobal("useGrid") == 1);  #  bin for the Grid

    die "Can't find local bin/gatekeeper in $bin\n" if (! -e "$bin/gatekeeper");
    die "Can't find grid bin/gatekeeper in $gin\n"  if (! -e "$gin/gatekeeper");
}


sub checkDirectories {

    #  Check that we were supplied a work directory, and that it
    #  exists, or we can create it.
    #
    die "ERROR: I need a directory to run the assembly in (-d option).\n" if (!defined($wrk));
    system("mkdir -p $wrk") if (! -d $wrk);
    die "ERROR: Directory '$wrk' doesn't exist (-d option) and couldn't be created.\n" if (! -d $wrk);

    #  Check that we have scratch space, or try to make one in the
    #  work directory.

    #  See if we can use the supplied scratch space
    #
    my $scratch = getGlobal("scratch");
    system("mkdir -p $scratch") if (! -d $scratch);

    #  If not created, warn, and try to make one in the work directory.
    #
    if (! -d $scratch) {
        print STDERR "WARNING: Scratch directory '$scratch' doesn't exist and couldn't be created; trying '$wrk/scratch' instead.\n";
        $scratch = "$wrk/scratch";
        system("mkdir -p $scratch");
        setGlobal("scratch", $scratch);
    }

    #  If still not created, die.
    #
    if (! -d $scratch) {
        die "ERROR:  Scratch directory '$scratch' doesn't exist, and couldn't be created!\n";
    }
}


sub findLastCheckpoint ($) {
    my $dir     = shift @_;
    my $lastckp;

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

    open(F, "cd $wrk/$dir && $bin/getNumScaffolds $asm $lastckp 2> /dev/null |");
    my $numscaf = <F>;  chomp $numscaf;
    close(F);
    $numscaf = int($numscaf);

    die "findNumSCaffoldsInCheckpoint($dir, $lastckp) found no scaffolds?!" if ($numscaf == 0);
    print STDERR "Found $numscaf in $dir checkpoint $lastckp\n";
    return($numscaf);
}


sub getNumberOfFragsInStore {
    my $bin = shift @_;
    my $wrk = shift @_;
    my $asm = shift @_;

    return(0) if (! -e "$wrk/$asm.frgStore/db.frg");

    open(F, "$bin/lastfraginstore $wrk/$asm.frgStore 2> /dev/null |") or die;
    $_ = <F>;    chomp $_;
    close(F);

    $numFrags = $1 if (m/^Last frag in store is iid = (\d+)$/);
    die "No frags in the store?\n" if ($numFrags == 0);
    return($numFrags);
}


#  Create an empty file.  Much faster than system("touch ...").
#
sub touch {
    open(F, "> $_[0]") or die "Failed to touch '$_[0]'\n";
    close(F);
}


sub pleaseExecute {
    my $file = shift @_;

    print STDERR "Please execute:\n";
    print STDERR "  $file\n";
    print STDERR "to submit jobs to the grid, then restart this script when all\n";
    print STDERR "jobs finish.  I'll make sure all jobs finished properly.\n";
}


#  Utility to run a command and check the exit status, report time used.
#
sub runCommand {
    my $cmd = shift @_;

    my $t = localtime();
    my $d = time();
    print STDERR "----------------------------------------START $t\n$cmd\n";

    my $rc = 0xffff & system($cmd);

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

    my $error = "ERROR: $cmd\n        failed with ";

    if ($rc == 0xff00) {
        $error .= "$!\n";
    } elsif ($rc > 0x80) {
        $rc >>= 8;
        $error .= "exit status $rc\n";
    } else {
        if ($rc & 0x80) {
            $rc &= ~0x80;
            $error .= "coredump from ";
        }
        if (defined($signame[$rc])) {
            $error .= "signal $signame[$rc]\n";
        } else {
            $error .= "signal $rc\n";
        }
    }

    print STDERR $error;

    return(1);
}

1;
