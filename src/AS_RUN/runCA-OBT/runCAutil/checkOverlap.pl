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

    open(F, "> $wrk/$outDir/overlap-restart.sh");
    print F "#!/bin/sh\n\n";

    while (scalar(@bat) > 0) {
        my $batchName = shift @bat;
        my $jobName   = shift @job;

        if (! -e "$wrk/$outDir/$batchName/$jobName.success") {
            print STDERR "$wrk/$outDir/$batchName/$jobName failed, job index $jobIndex.\n";

            if (getGlobal("useGrid") && getGlobal("ovlOnGrid")) {
                my $sge        = getGlobal("sge");
                my $sgeOverlap = getGlobal("sgeOverlap");

                print F "qsub $sge $sgeOverlap -r y -N ovl_${asm} \\\n";
                print F "  -t $jobIndex \\\n";
                print F "  -j y -o $wrk/$outDir/overlap.\\\$TASK_ID.out \\\n";
                print F "  -e $wrk/$outDir/overlap.\\\$TASK_ID.err \\\n";
                print F "  $wrk/$outDir/overlap.sh\n";
            } else {
                my $out = substr("0000" . $jobIndex, -4);
                print F "$wrk/$outDir/overlap.sh $jobIndex > $wrk/$outDir/$out.out 2>&1\n";
            }

            $failedJobs++;
        }

        $jobIndex++;
    }

    close(F);

    if ($failedJobs) {
        caFailure("$failedJobs failed.  See $wrk/$outDir/overlap-restart.sh for resubmission commands.\n");
    }
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

        if (! -e "$wrk/$outDir/olaps/$asm.$job.success") {
            print STDERR "$wrk/$outDir/olaps/$job failed.\n";
            $failedJobs++;
        }
    }
    if ($failedJobs) {
        caFailure("$failedJobs failed.  See $wrk/$outDir/overlap-restart.sh for resubmission commands.\n");
    }
}


sub checkOverlap {
    my $isTrim = shift @_;

    return if (-d "$wrk/$asm.ovlStore");

    caFailure("checkOverlap()-- I need to know if I'm trimming or assembling!\n") if (!defined($isTrim));

    if ((getGlobal("merOverlap") eq "both") ||
        ((getGlobal("merOverlap") eq "obt") && ($isTrim eq "trim")) ||
        ((getGlobal("merOverlap") eq "ovl") && ($isTrim ne "trim"))) {
        checkMerOverlapper($isTrim);
    } else {
        checkOverlapper($isTrim);
    }
}

1;

