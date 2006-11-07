use strict;

#  Check and merge the fragment correction results

sub mergeFragmentCorrection {

    return if (getGlobal("doFragmentCorrection") == 0);

    my $frgCorrBatchSize  = getGlobal("frgCorrBatchSize");

    if (! -e "$wrk/2-frgcorr/$asm.corr") {
        my $failedJobs = 0;

        open(F, "> $wrk/2-frgcorr/cat-corrects.corrlist");

        my $jobs = int($numFrags / ($frgCorrBatchSize-1)) + 1;
        for (my $i=1; $i<=$jobs; $i++) {
            my $frgBeg = $i * $frgCorrBatchSize - $frgCorrBatchSize + 1;
            my $frgEnd = $i * $frgCorrBatchSize;
            if ($frgEnd > $numFrags) {
                $frgEnd = $numFrags - 1;
            }
            $frgBeg = substr("00000000$frgBeg", -8);
            $frgEnd = substr("00000000$frgEnd", -8);

            if (! -e "$wrk/2-frgcorr/$asm-$frgBeg-$frgEnd.success") {
                print STDERR "Fragment correction job $i ($wrk/2-frgcorr/$asm-$frgBeg-$frgEnd) failed.\n";
                $failedJobs++;
            }

            print F "$wrk/2-frgcorr/$asm-$frgBeg-$frgEnd.frgcorr\n";
        }

        close(F);

        die "$failedJobs failed.  Good luck.\n" if ($failedJobs);
            
        my $cmd;
        $cmd  = "$bin/cat-corrects ";
        $cmd .= "-L $wrk/2-frgcorr/cat-corrects.corrlist ";
        $cmd .= "-o $wrk/2-frgcorr/$asm.corr ";
        $cmd .= "> $wrk/2-frgcorr/cat-corrects.err 2>&1";
        if (runCommand("$wrk/2-frgcorr", $cmd)) {
            rename "$wrk/2-frgcorr/$asm.corr", "$wrk/2-frgcorr/$asm.corr.FAILED";
            die "Failed to concatenate the fragment corrections.\n";
        }
    }
}

1;
