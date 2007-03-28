use strict;

sub createOverlapStore {
    my $ovlStoreMemory    = getGlobal("ovlStoreMemory");

    goto alldone if (-d "$wrk/$asm.ovlStore");

    if (! -e "$wrk/1-overlapper/all-overlaps.ovllist") {
        if (runCommand("$wrk/1-overlapper",
                       "find $wrk/1-overlapper/ \\( -name \\*ovb -o -name \\*ovb.bz2 \\) -print > $wrk/1-overlapper/all-overlaps.ovllist")) {
            rename "$wrk/1-overlapper/all-overlaps.ovllist", "$wrk/1-overlapper/all-overlaps.ovllist.FAILED";
            die "Failed to generate a list of all the overlap files.\n";
        }
    }

    my $cmd;
    $cmd  = "$bin/overlapStore ";
    $cmd .= "-c $wrk/$asm.ovlStore ";
    $cmd .= "-M $ovlStoreMemory ";
    $cmd .= "-m $numFrags ";
    $cmd .= "-L $wrk/1-overlapper/all-overlaps.ovllist ";
    $cmd .= "> $wrk/1-overlapper/grow-olap-store.err 2>&1";

    if (runCommand("$wrk", $cmd)) {
        rename "$wrk/$asm.ovlStore", "$wrk/$asm.ovlStore.FAILED";
        die "Failed to create the overlap store.\n";
    }

  alldone:
    stopAfter("overlapper");
}

1;
