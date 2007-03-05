use strict;


sub preoverlap {
    my @fragFiles = @_;

    $numFrags = getNumberOfFragsInStore($bin, $wrk, $asm);

    #  Return if there are fragments in the store, and die if there
    #  are no fragments and no source files.
    #
    goto alldone if ($numFrags > 0);

    die "ERROR: No fragment files specified, and stores not already created.\n" if (scalar(@fragFiles) == 0);

    system("mkdir $wrk/0-preoverlap") if (! -d "$wrk/0-preoverlap");

    if ((! -d "$wrk/$asm.gkpStore") || (! -e "$wrk/$asm.gkpStore/frg")) {

        #  Make sure all the inputs are here
        #
        my $failedFiles = 0;
        my $gkpInput = "";
        foreach my $frg (@fragFiles) {
            if (! -e $frg) {
                print STDERR "MISSING: $frg\n";
                $failedFiles++;
            }
            $gkpInput .= " $frg";
        }
        die if ($failedFiles);

        my $cmd;
        $cmd  = "$bin/gatekeeper -e 10000000 ";
        $cmd .= "-Q -T " if (getGlobal("doOverlapTrimming"));
        $cmd .= "-o $wrk/$asm.gkpStore ";
        $cmd .= "$gkpInput ";
        $cmd .= "> $wrk/0-preoverlap/gatekeeper.out ";
        $cmd .= "2> $wrk/0-preoverlap/gatekeeper.err";

        if (runCommand("$wrk/0-preoverlap", $cmd)) {
            print STDERR "Failed.\n";
            rename "$wrk/0-preoverlap/$asm.inp", "$wrk/0-preoverlap/$asm.inp.FAILED";
            rename "$wrk/$asm.gkpStore", "$wrk/$asm.gkpStore.FAILED";
            exit(1);
        }
    }

    $numFrags = getNumberOfFragsInStore($bin, $wrk, $asm);

alldone:
    stopAfter("initialStoreBuilding");
}

1;
