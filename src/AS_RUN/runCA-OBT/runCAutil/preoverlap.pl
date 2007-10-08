use strict;


sub preoverlap {
    my @fragFiles = @_;

    $numFrags = getNumberOfFragsInStore($bin, $wrk, $asm);

    #  Return if there are fragments in the store, and die if there
    #  are no fragments and no source files.
    #
    if ($numFrags > 0) {
        goto stopafter;
    }

    caFailure("ERROR: No fragment files specified, and stores not already created.\n")
    	if (scalar(@fragFiles) == 0);

    system("mkdir $wrk/0-preoverlap") if (! -d "$wrk/0-preoverlap");

    if ((! -d "$wrk/$asm.gkpStore") ||
        (! -e "$wrk/$asm.gkpStore/frg")) {

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
        caFailure("Files supplied on command line not found.\n") if ($failedFiles);

        my $cmd;
        $cmd  = "$bin/gatekeeper -o $wrk/$asm.gkpStore ";
        $cmd .= " -T " if (getGlobal("doOverlapTrimming"));
        $cmd .= " -D " if (getGlobal("gkpBelieveInputStdDev"));
        $cmd .= " -E $wrk/0-preoverlap/gatekeeper.errors ";
        $cmd .= "$gkpInput ";
        $cmd .= "> $wrk/0-preoverlap/gatekeeper.err 2>&1";
        if (runCommand("$wrk/0-preoverlap", $cmd)) {
            print STDERR "Failed.\n";
            rename "$wrk/0-preoverlap/$asm.inp", "$wrk/0-preoverlap/$asm.inp.FAILED";
            rename "$wrk/$asm.gkpStore", "$wrk/$asm.gkpStore.FAILED";
            caFailure("gatekeeper failed.  Check your input files!\n");
        }
    }

    my $vi = getGlobal("vectorIntersect");
    if ((defined($vi)) && (! -e "$wrk/0-preoverlap/$asm.vectorClearLoaded")) {
        if (runCommand("$wrk/0-preoverlap", "$bin/gatekeeper -a -v $vi -o $wrk/$asm.gkpStore > $wrk/0-preoverlap/$asm.vectorClearLoaded.err 2>&1")) {
            caFailure("Failed.\n");
        }
        touch("$wrk/0-preoverlap/$asm.vectorClearLoaded");
    }
        
    $numFrags = getNumberOfFragsInStore($bin, $wrk, $asm);

    if (runCommand($wrk, "$bin/gatekeeper -dumpfrg $wrk/$asm.gkpStore 2> $wrk/gatekeeper.err | grep -v 'No source' > $wrk/$asm.frg")) {
        unlink "$wrk/$asm.frg";
    }    

 stopafter:
    stopAfter("initialStoreBuilding");    
}

1;
