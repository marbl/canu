use strict;

############################################################
#
#  Utility functions
#
############################################################


#  Decide what host we are on, and the the bin directory
#  appropriately
#
sub setBinDirectory {
    my $host    = `uname`;     chomp $host;
    my $mach    = `uname -m`;  chomp $mach;
    my $binRoot;
    my $binDir;

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
    my $t = getGlobal("binRoot", undef);
    if (defined($binRoot) && (defined($t))) {
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

    #  Finally, we might be being forced to configure for a host that
    #  isn't the one we are currently executing on -- this lets us set
    #  the host of the grid.
    #
    $host    = shift @_ if (scalar(@_));
    $mach    = shift @_ if (scalar(@_));

    if      (($host eq "Linux") && ($mach eq "i686")) {
        #  Linux, on Intel 686
        $binDir       = "$binRoot/Linux/bin";
        $scratch      = "/scratch";
        $processStats = "/home/bwalenz/linux/bin/processstats";
    } elsif (($host eq "Linux") && ($mach eq "x86_64")) {
        #  Linux on Opteron
        $binDir       = "$binRoot/Linux64/bin";
        $scratch      = "/scratch";
        $processStats = "time";
    } elsif (($host eq "FreeBSD") && ($mach eq "i386")) {
        #  FreeBSD on Intel
        $binDir       = "$binRoot/FreeBSD/bin";
        $scratch      = "/scratch";
        $processStats = "time";
    } elsif ($host eq "Darwin") {
        #  Darwin (UNTESTED)
        $binDir       = "$binRoot/Darwin/bin";
        $scratch      = "/scratch";
        $processStats = "time";
    } elsif (($mach eq "alpha") || ($mach eq "OSF1")) {
        #  OSF1 on alpha
        $binDir       = "$binRoot/OSF1/bin";
        $scratch      = "/bioinfo/assembly/walenz/scratch";
        $processStats = "time";
    } else {
        print STDERR "Unknown host '$host' and machine '$mach' combination.\n";
        exit(1);
    }

    if (! -d $binRoot) {
        die "ERROR: Failed to find the binary directory $binDir.\n";
    }

    return($binDir);
}



my %global;


#  Return the second argument, unless the first argument is found in
#  %global, in which case return that.
#
sub getGlobal {
    my $var = shift @_;
    my $val = shift @_;
    $val = $global{$var} if (defined($global{$var}));
    #print STDERR "spec: return '$val' -> '$var'\n";
    return($val);
}


sub setParameters {
    my $specFile = shift @_;

    if (defined($specFile)) {
        open(F, "< $specFile") or die "Failed to open '$specFile'\n";
        while (<F>) {
            my ($var, $val) = split '=', $_;
            $var =~ s/^\s+//; $var =~ s/\s+$//;
            $val =~ s/^\s+//; $val =~ s/\s+$//;
            if (($var ne "") && ($val ne "")) {
                $global{$var} = $val;
                print STDERR "spec: '$var' <- '$val'\n";
            }
        }
        close(F);
    }


    ############################################################
    #
    #  Job Parameters
    #

    #  If we have been given a configuration in the spec file, use that instead.
    #
    my $binhost = getGlobal("localHost",    undef);
    my $binmach = getGlobal("localMachine", undef);
    my $ginhost = getGlobal("gridHost",     "Linux");
    my $ginmach = getGlobal("gridMachine",  "i686");

    $useGrid = getGlobal("useGrid", 1);
    $bin     = setBinDirectory($binhost, $binmach);
    $gin     = setBinDirectory($binhost, $binmach)  if ($useGrid == 0);  #  bin for the Grid
    $gin     = setBinDirectory($ginhost, $ginmach)  if ($useGrid == 1);  #  bin for the Grid



    #print STDERR "bin = $bin\n";
    #print STDERR "gin = $gin\n";


    ############################################################
    #
    #  Overlap based trim parameters
    #
    $doOverlapTrimming  = getGlobal("doOverlapTrimming", 1);
    $vectorIntersect    = getGlobal("vectorIntersect", "");


    ############################################################
    #
    #  Overlapper parameters
    #
    #  ovlThreads -- number of threads per overlapper job.  On the VI
    #  grid, hosts are hyper-threaded dual-Xeons.  Previous trivial tests
    #  indicated that using two threads per CPU gave about 20% better
    #  performance than one thread per CPU.
    #
    #  ovlHashBlockSize -- number of fragments to build the in-core hash
    #  table with.  200,000 fragments of mean size ~550bp ran in 4GB at
    #  Celera.  VI fragments are ~800bp.  A host at VI has 2GB, and can
    #  run two processes, so figure on ~600MB per job.  The script for
    #  running dog used 40000.  30000 was used here to prevent any chance
    #  of paging.
    #
    #  ovlRefBlockSize -- to better utilize CPU and to make jobs shorter, we
    #  can segment the number of fragments we run by the overlapper in one
    #  run.  1,000,000 (with 40,000 ovlHaskBlockSize) is reported to give
    #  about one hour of run time.
    #
    $ovlThreads        = getGlobal("ovlThreads", 2);
    $ovlHashBlockSize  = getGlobal("ovlHashBlockSize", 40000);
    $ovlRefBlockSize   = getGlobal("ovlRefBlockSize", 2000000);
    $ovlMemory         = getGlobal("ovlMemory", "1GB");

    $ovlThreads        = getGlobal("ovlThreads", 2);
    $ovlHashBlockSize  = getGlobal("ovlHashBlockSize", 150000);
    $ovlRefBlockSize   = getGlobal("ovlRefBlockSize", 5000000);
    $ovlMemory         = getGlobal("ovlMemory", "2GB");

    $ovlStoreMemory    = getGlobal("ovlStoreMemory", 16384);

    ############################################################
    #
    #  Fragment correction parameters
    #
    #  Using the store directly (recent commit to the tree, compile time
    #  option), testing shows that, for a human assembly using 22-mer
    #  and standard overlaps:
    #
    #      10,000 frags per batch needs   132 MB
    #      50,000 frags per batch needs   650 MB
    #     100,000 frags per batch needs  1300 MB
    #     200,000 frags per batch needs  2500 MB
    #     500,000 frags per batch needs  6300 MB
    #   1,000,000 frags per batch needs 13000 MB
    #   2,000,000 frags per batch needs 23000 MB
    #   2,500,000 frags per batch needs 30000 MB (died)
    #
    #  3 million fragments work on a 32GB box (usually), assuming
    #  correct-frags doesn't use a temporary internal store.
    #
    #  This will do 12 batches of roughly the same size.
    #
    $frgCorrBatchSize  = getGlobal("frgCorrBatchSize", int(26742484 / 12));
    $frgCorrBatchSize  = getGlobal("frgCorrBatchSize", 175000);
    $frgCorrThreads    = getGlobal("frgCorrThreads", 2);
    $frgCorrOnGrid     = getGlobal("frgCorrOnGrid", 0);


    ############################################################
    #
    #  Overlap correction parameters
    #
    #  Don't know anything about sizes here....
    #
    $ovlCorrBatchSize    = getGlobal("ovlCorrBatchSize", int(26742484 / 12));
    $ovlCorrBatchSize    = getGlobal("ovlCorrBatchSize", 175000);
    $ovlCorrOnGrid       = getGlobal("ovlCorrOnGrid", 0);


    ############################################################
    #
    #  CGW parameters
    #
    $stoneLevel          = getGlobal("throwStones", 2);


    ############################################################
    #
    #  Extend Clear Ranges -- don't throw stones if you set this.
    #
    #  N.B. NOT IMPLEMENTED
    #
    $doExtendClearRanges = getGlobal("doExtendClearRanges", 0);
}







sub findLastCheckpoint {
    my $dir     = shift @_;
    my $lastckp = "000";

    open(F, "ls -1 $wrk/$dir/*ckp* |");
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





#  Create an empty file.  Much faster than system("touch ...").
#
sub touch {
    open(F, "> $_[0]") or die "Failed to touch '$_[0]'\n";
    close(F);
}


#  Utility to run a command and check the exit status, grabbed from
#  ESTmapper util
#
sub runCommand {
    my $cmd = shift @_;

    print STDERR "----------------------------------------START\n$cmd\n";

    my $rc = 0xffff & system($cmd);

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


sub getNumberOfFragsInStore {
    my $bin = shift @_;
    my $wrk = shift @_;
    my $asm = shift @_;

    return(0) if (! -e "$wrk/$asm.frgStore");

    open(F, "$bin/lastfraginstore $wrk/$asm.frgStore 2> /dev/null |") or die;
    $_ = <F>;    chomp $_;
    close(F);

    $numFrags = $1 if (m/^Last frag in store is iid = (\d+)$/);
    die "No frags in the store?\n" if ($numFrags == 0);
    return($numFrags);
}


sub pleaseExecute {
    my $file = shift @_;

    print STDERR "Please execute:\n";
    print STDERR "  $file\n";
    print STDERR "to submit jobs to the grid, then restart this script when all\n";
    print STDERR "jobs finish.  I'll make sure all jobs finished properly.\n";
}



sub copyStoresLocally {
    my $bin = shift @_;
    my $wrk = shift @_;
    my $asm = shift @_;

    if ((-e "$scratch/$asm.frgStore") && (-e "$scratch/$asm.ovlStore") && (! -e "$scratch/$asm.copying")) {
        print STDERR "Found local frgStore and local ovlStore!\n";
        return;
    }

    #  Someone already copying the data.  Wait for completion.
    #
    if (-e "$scratch/$asm.copying") {
        print STDERR "Waiting for someone to copy local frgStore and local ovlStore!\n";
        while (-e "$scratch/$asm.copying") {
            sleep 10;
        }
        print STDERR "Found local frgStore and local ovlStore!\n";
        return;
    }

    #  Nobody copying, start the copies.
    #
    open(F, "> $scratch/$asm.copying");
    close(F);

    open(F, "df -k $scratch | ") or die "Failed.\n";
    $_ = <F>;
    $_ = <F>;
    my @df = split '\s+', $_;
    my $freeSpace = $df[3];
    close(F);

    print STDERR "Found $freeSpace free KB\n";

    if ($freeSpace > 24743977 + 12815696 + 2048) {
        print STDERR "Copying.\n";
        system("date");
        if (runCommand("cp -rp /bioinfo/assembly/projects/HuRef31/$asm.frgStore $scratch/$asm.frgStore")) {
            print STDERR "Failed copy frgStore!\n";
            system("rm -rf $scratch/$asm.frgStore");
        }
        if (runComand("cp -rp /bioinfo/assembly/projects/HuRef31/$asm.ovlStore $scratch/$asm.ovlStore")) {
            print STDERR "Failed copy frgStore!\n";
            system("rm -rf $scratch/$asm.ovlStore");
        }
        system("date");
    }

    if (! -d "$scratch/$asm.frgStore") {
        print STDERR "Symlink frgStore.\n";
        if (runCommand("ln -s /bioinfo/assembly/projects/Huref31/$asm.frgStore $scratch/$asm.frgStore")) {
            print STDERR "Failed ln -s frgStore!\n";
        }
    }
    if (! -d "$scratch/$asm.ovlStore") {
        print STDERR "Symlink gkpStore.\n";
        if (runCommand("ln -s /bioinfo/assembly/projects/Huref31/$asm.ovlStore $scratch/$asm.ovlStore")) {
            print STDERR "Failed ln -s ovlStore!\n";
        }
    }

    unlink "$scratch/$asm.copying";
}


#  Attempt to check what the process limits are.  Print out and make
#  the user look, but don't check.
#
#  We use /bin/sh to run subcommands on the farm, and system() in here
#  (which appears to be using /bin/sh, or nothing at all, depending on
#  if you use shell features or not).  Supposedly, 'limit' works in
#  sh, but I see differences between Tru64 and Linux, so we just make
#  the user take care of fixing limits.
#
sub checkProcessLimits {
    open (F, "> /tmp/ch$$.tcsh");
    print F "#!/bin/tcsh\n";
    print F "limit\n";
    close(F);
    my @limits;
    my $datasize = 0;
    open(F, "tcsh /tmp/ch$$.tcsh |");
    while (<F>) {
        push @limits, "****    $_";
        $datasize = $1        if (m/datasize\s+(\d+)\s+kbytes/);
        $datasize = 999999999 if (m/datasize\s+unlimited/);
    }
    close(F);
    unlink "/tmp/ch$$.tcsh";

    if ($datasize < 3000000) {
        print STDERR "****\n";
        print STDERR "****    CHECK YOUR PROCESS LIMITS!  Ensure 'datasize' is big.\n";
        print STDERR "****    Fix with (tcsh) 'unlimit' or (bash) 'ulimit -a'\n";
        print STDERR "****\n";
        print @limits;
        print STDERR "****\n";
        exit;
    }
}



############################################################
#
#  End of utility functions
#
############################################################

1;
