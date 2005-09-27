use strict;

#  Create the overlap store
#
#  This could be split into smaller batches, use option -A instead of -cf

sub createOverlapStore {

    if (! -e "$wrk/1-overlapper/all-overlaps.ovllist") {
        if (runCommand("find $wrk/1-overlapper/ -name \\*ovl -print > $wrk/1-overlapper/all-overlaps.ovllist")) {
            print STDERR "Failed to generate a list of all the overlap files.\n";
            rename "$wrk/1-overlapper/all-overlaps.ovllist", "$wrk/1-overlapper/all-overlaps.ovllist.FAILED";
            exit(1);
        }
    }

    if (! -e "$wrk/$asm.ovlStore") {
        print STDERR "Starting -- make overlap store\n";

        my $cmd;
        $cmd  = "$bin/grow-olap-store ";
        $cmd .= "-cfS ";
        $cmd .= "-M $ovlStoreMemory ";
        $cmd .= "-o $wrk/$asm.ovlStore ";
        $cmd .= "-L $wrk/1-overlapper/all-overlaps.ovllist ";
        $cmd .= "> $wrk/1-overlapper/grow-olap-store.out ";
        $cmd .= "2> $wrk/1-overlapper/grow-olap-store.err ";
        if (runCommand($cmd)) {
            print STDERR "Failed to grow the overlap store.\n";
            rename "$wrk/$asm.ovlStore", "$wrk/$asm.ovlStore.FAILED";
            exit(1);
        }
    }
}

1;
