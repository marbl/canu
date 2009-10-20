use strict;

sub createOverlapStore {

    goto alldone if (-d "$wrk/$asm.ovlStore");

    if (runCommand($wrk, "find $wrk/1-overlapper \\( -name \\*ovb.gz -or -name \\*ovb \\) -print > $wrk/$asm.ovlStore.list")) {
        caFailure("failed to generate a list of all the overlap files", undef);
    }

    my $bin = getBinDirectory();
    my $cmd;
    $cmd  = "$bin/overlapStore ";
    $cmd .= " -c $wrk/$asm.ovlStore.BUILDING ";
    $cmd .= " -g $wrk/$asm.gkpStore ";
    
    if (defined(getGlobal("closureOverlaps"))){
      $cmd .= " -i " . getGlobal("closureOverlaps");
    }
    
    $cmd .= " -M " . getGlobal("ovlStoreMemory");
    $cmd .= " -L $wrk/$asm.ovlStore.list ";
    $cmd .= " > $wrk/$asm.ovlStore.err 2>&1";

    if (runCommand($wrk, $cmd)) {
        caFailure("failed to create the overlap store", "$wrk/$asm.ovlStore.err");
    }

    rename "$wrk/$asm.ovlStore.BUILDING", "$wrk/$asm.ovlStore";

    if (getGlobal("saveOverlaps") == 0) {
        open(F, "< $wrk/$asm.ovlStore.list");
        while (<F>) {
            chomp;
            unlink $_;
        }
        close(F);
    }

    rmrf("$wrk/$asm.ovlStore.list");
    rmrf("$wrk/$asm.ovlStore.err");

  alldone:
    stopAfter("overlapper");
}

1;
