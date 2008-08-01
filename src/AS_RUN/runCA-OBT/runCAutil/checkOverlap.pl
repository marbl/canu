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

    open(F, "< $wrk/$outDir/ovljobs.dat") or caFailure("Failed to open '$wrk/$outDir/ovljobs.dat'\n");
    $_ = <F>;
    my @bat = split '\s+', $_;
    $_ = <F>;
    my @job = split '\s+', $_;
    close(F);

    my $jobIndex   = 1;
    my $failedJobs = 0;

    while (scalar(@bat) > 0) {
        my $batchName = shift @bat;
        my $jobName   = shift @job;

        if (! -e "$wrk/$outDir/$batchName/$jobName.ovb.gz") {
            print STDERR "$wrk/$outDir/$batchName/$jobName failed, job index $jobIndex.\n";
            $failedJobs++;
        }

        $jobIndex++;
    }

    caFailure("$failedJobs overlapper jobs failed.\n") if ($failedJobs);
}


sub checkMerOverlapper ($) {
    my $isTrim = shift @_;

    my $outDir = "1-overlapper";

    if ($isTrim eq "trim") {
        $outDir = "0-overlaptrim-overlap";
    }

    my $batchSize  = getGlobal("ovlCorrBatchSize");
    my $jobs       = int($numFrags / ($batchSize-1)) + 1;
    my $failedJobs = 0;

    for (my $i=1; $i<=$jobs; $i++) {
        my $job = substr("0000" . $i, -4);

        if (! -e "$wrk/$outDir/olaps/$job.success") {
            print STDERR "$wrk/$outDir/olaps/$job failed.\n";
            $failedJobs++;
        }
    }
    
    caFailure("$failedJobs overlapper jobs failed.\n") if ($failedJobs);
}


sub checkOverlap {
    my $isTrim = shift @_;

    caFailure("checkOverlap()-- I need to know if I'm trimming or assembling!\n") if (!defined($isTrim));

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

