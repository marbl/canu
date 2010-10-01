use strict;

#  Check that the overlapper jobs properly executed.  If not,
#  complain, but don't help the user fix things.


sub checkOverlapper ($) {
    my $isTrim = shift @_;

    my $outDir = "1-overlapper";
    my $ovlOpt = "";

    if ($isTrim eq "trim") {
        $outDir = "0-overlaptrim-overlap";
        $ovlOpt = "-G";
    }

    my $failedJobs = 0;
    my $failureMessage = "";

    open(B, "< $wrk/$outDir/ovlbat") or caFailure("failed to open '$wrk/$outDir/ovlbat'", undef);
    open(J, "< $wrk/$outDir/ovljob") or caFailure("failed to open '$wrk/$outDir/ovljob'", undef);

    while (!eof(B) && !eof(J)) {
        my $b = <B>;  chomp $b;
        my $j = <J>;  chomp $j;

        if ((! -e "$wrk/$outDir/$b/$j.ovb.gz") &&
            (! -e "$wrk/$outDir/$b/$j.ovb")) {
            $failureMessage .= "ERROR:  Overlap job $wrk/$outDir/$b/$j FAILED.\n";
            $failedJobs++;
        }
    }

    if (!eof(B) || !eof(J)) {
        print STDERR "Partitioning error; '$wrk/$outDir/ovlbat' and '$wrk/$outDir/ovljob' have extra lines.\n";
    }

    #  FAILUREHELPME
    #
    $failureMessage .= "\n$failedJobs overlapper jobs failed";
    caFailure($failureMessage, undef) if ($failedJobs);
}


sub checkMerOverlapper ($) {
    my $isTrim = shift @_;

    my $outDir = "1-overlapper";

    if ($isTrim eq "trim") {
        $outDir = "0-overlaptrim-overlap";
    }

    my $batchSize  = getGlobal("merOverlapperExtendBatchSize");
    my $jobs       = int($numFrags / $batchSize) + (($numFrags % $batchSize == 0) ? 0 : 1);
    my $failedJobs = 0;

    for (my $i=1; $i<=$jobs; $i++) {
        my $job = substr("0000" . $i, -4);

        if ((! -e "$wrk/$outDir/olaps/$job.ovb.gz") &&
            (! -e "$wrk/$outDir/olaps/$job.ovb")) {
            print STDERR "$wrk/$outDir/olaps/$job failed.\n";
            $failedJobs++;
        }
    }
    
    caFailure("$failedJobs overlapper jobs failed", undef) if ($failedJobs);
}


sub checkOverlap {
    my $isTrim = shift @_;

    caFailure("overlap checker needs to know if trimming or assembling", undef) if (!defined($isTrim));

    if ($isTrim eq "trim") {
        return if (-d "$wrk/$asm.obtStore");
        if      (getGlobal("obtOverlapper") eq "ovl") {
            checkOverlapper($isTrim);
        } elsif (getGlobal("obtOverlapper") eq "mer") {
            checkMerOverlapper($isTrim);
        } elsif (getGlobal("obtOverlapper") eq "umd") {
            caError("checkOverlap() wanted to check umd overlapper for obt?\n");
        } else {
            caError("checkOverlap() unknown obt overlapper?\n");
        }
    } else {
        return if (-d "$wrk/$asm.ovlStore");
        if      (getGlobal("ovlOverlapper") eq "ovl") {
            checkOverlapper($isTrim);
        } elsif (getGlobal("ovlOverlapper") eq "mer") {
            checkMerOverlapper($isTrim);
        } elsif (getGlobal("ovlOverlapper") eq "umd") {
            #  Nop.
        } else {
            caError("checkOverlap() unknown ovl overlapper?\n");
        }
    }
}

1;

