use strict;

#  Create the overlap store
#
#  This could be split into smaller batches, use option -A instead of -cf

sub createOverlapStore {
    my $ovlStoreMemory    = getGlobal("ovlStoreMemory");

    return if (-d "$wrk/$asm.ovlStore");

    if (! -e "$wrk/1-overlapper/all-overlaps.ovllist") {
        if (runCommand("find $wrk/1-overlapper/ -name \\*ovl -print > $wrk/1-overlapper/all-overlaps.ovllist")) {
            rename "$wrk/1-overlapper/all-overlaps.ovllist", "$wrk/1-overlapper/all-overlaps.ovllist.FAILED";
            die "Failed to generate a list of all the overlap files.\n";
        }
    }

    my $cmd;
    $cmd  = "$bin/grow-olap-store ";
    $cmd .= "-cfS ";
    $cmd .= "-M $ovlStoreMemory ";
    $cmd .= "-o $wrk/$asm.ovlStore ";
    $cmd .= "-L $wrk/1-overlapper/all-overlaps.ovllist ";
    $cmd .= "> $wrk/1-overlapper/grow-olap-store.out ";
    $cmd .= "2> $wrk/1-overlapper/grow-olap-store.err ";
    if (runCommand($cmd)) {
        rename "$wrk/$asm.ovlStore", "$wrk/$asm.ovlStore.FAILED";
        die "Failed to grow the overlap store.\n";
    }
}

1;
