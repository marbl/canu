use strict;

#  Check that the overlapper jobs properly executed.  If not,
#  complain, but don't help the user fix things.

sub checkOverlap {
    my $isTrim = shift @_;

    if (!defined($isTrim)) {
        die "checkOverlap()-- I need to know if I'm trimming or assembling!\n";
    }

    my $outDir = "1-overlapper";
    my $ovlOpt = "";

    if ($isTrim eq "trim") {
        $outDir = "0-overlaptrim-overlap";
        $ovlOpt = "-G";
    }

    open(F, "< $wrk/$outDir/ovljobs.dat") or die "Failed to open '$wrk/$outDir/ovljobs.dat'\n";
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

            if ($useGrid) {
                print F "qsub -p 0 -r y -N ovl_${asm} \\\n";
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
        print STDERR "$failedJobs failed.  See $wrk/$outDir/overlap-restart.sh for resubmission.\n";
        exit(1);
    }
}

1;

