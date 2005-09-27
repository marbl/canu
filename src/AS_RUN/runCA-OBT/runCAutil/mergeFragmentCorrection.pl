use strict;

#  Check and merge the fragment correction results

sub mergeFragmentCorrection {

    if (! -e "$wrk/2-frgcorr/$asm.corr") {
        my $failedJobs = 0;

        my $frgBeg = 1;
        my $frgEnd = 0;

        while ($frgBeg < $numFrags) {
            $frgEnd = $frgBeg + $frgCorrBatchSize - 1;
            $frgEnd = $numFrags if ($frgEnd > $numFrags);

            my $jobName   = substr("0000000000" . $frgBeg, -8);
            my $batchName = substr($jobName, 0, 2);

            if (! -e "$wrk/2-frgcorr/$batchName/$jobName.success") {
                print STDERR "$wrk/2-frgcorr/$batchName/$jobName failed.\n";
                $failedJobs++;
            }

            $frgBeg = $frgEnd + 1;
        }

        if ($failedJobs) {
            print STDERR "$failedJobs failed.  Good luck.\n";
            exit(1);
        }

        ########################################

        if (! -e "$wrk/2-frgcorr/all-corrections.corrlist") {
            if (runCommand("find $wrk/2-frgcorr -name \\*.frgcorr -print | sort > $wrk/2-frgcorr/all-corrections.corrlist")) {
                print STDERR "Failed to generate a list of all the fragment correction files.\n";
                rename "$wrk/2-frgcorr/all-corrections.corrlist", "$wrk/2-frgcorr/all-corrections.corrlist.FAILED";
                exit(1);
            }
        }

        ########################################

        my $cmd;
        $cmd  = "$bin/cat-corrects ";
        $cmd .= "-L $wrk/2-frgcorr/all-corrections.corrlist ";
        $cmd .= "-o $wrk/2-frgcorr/$asm.corr ";
        $cmd .= "> $wrk/2-frgcorr/cat-corrects.out ";
        $cmd .= "2> $wrk/2-frgcorr/cat-corrects.err";
        if (runCommand($cmd)) {
            print STDERR "Failed to concatenate the fragment corrections.\n";
            rename "$wrk/2-frgcorr/$asm.corr", "$wrk/2-frgcorr/$asm.corr.FAILED";
            exit(1);
        }
    }
}

1;
