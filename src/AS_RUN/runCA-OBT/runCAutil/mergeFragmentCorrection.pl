use strict;

#  Check and merge the fragment correction results

sub mergeFragmentCorrection {

    return if (getGlobal("doFragmentCorrection") == 0);

    my $batchSize  = getGlobal("frgCorrBatchSize");
    my $jobs       = int($numFrags / ($batchSize-1)) + 1;
    my $failedJobs = 0;

    if (! -e "$wrk/2-frgcorr/$asm.corr") {

        open(F, "> $wrk/2-frgcorr/cat-corrects.corrlist");
        for (my $i=1; $i<=$jobs; $i++) {
            my $jobid = substr("0000" . $i, -4);

            if (! -e "$wrk/2-frgcorr/$jobid.success") {
                print STDERR "Fragment correction job $jobid failed.\n";
                $failedJobs++;
            }

            print F "$wrk/2-frgcorr/$jobid.frgcorr\n";
        }
        close(F);

        caFailure("$failedJobs failed.  Good luck.\n") if ($failedJobs);
            
        my $bin = getBinDirectory();
        my $cmd;
        $cmd  = "$bin/cat-corrects ";
        $cmd .= "-L $wrk/2-frgcorr/cat-corrects.corrlist ";
        $cmd .= "-o $wrk/2-frgcorr/$asm.corr ";
        $cmd .= "> $wrk/2-frgcorr/cat-corrects.err 2>&1";

        if (runCommand("$wrk/2-frgcorr", $cmd)) {
            rename "$wrk/2-frgcorr/$asm.corr", "$wrk/2-frgcorr/$asm.corr.FAILED";
            caFailure("Failed to concatenate the fragment corrections.\n");
        }
    }
}

1;
