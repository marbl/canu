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

    while (scalar(@bat) > 0) {
        my $batchName = shift @bat;
        my $jobName   = shift @job;

        if (! -e "$wrk/$outDir/$batchName/$jobName.success") {
            print STDERR "$wrk/$outDir/$batchName/$jobName failed, job index $jobIndex.\n";
            $failedJobs++;
        }

        $jobIndex++;
    }

    if ($failedJobs) {
        print STDERR "$failedJobs failed.  Good luck.\n\n";
        print STDERR "(Yes, I should give you some help, I guess.)\n";
        exit(1);
    }
}

1;

