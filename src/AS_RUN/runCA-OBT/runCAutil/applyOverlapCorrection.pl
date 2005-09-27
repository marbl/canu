use strict;

#  Check and apply the overlap correction results.

sub applyOverlapCorrection {

    if (! -e "$wrk/3-ovlcorr/$asm.eratesupdated") {
        my $failedJobs = 0;

        my $frgBeg = 1;
        my $frgEnd = 0;

        while ($frgBeg < $numFrags) {
            $frgEnd = $frgBeg + $ovlCorrBatchSize - 1;
            $frgEnd = $numFrags if ($frgEnd > $numFrags);

            my $jobName   = substr("0000000000" . $frgBeg, -8);

            if (! -e "$wrk/3-ovlcorr/$jobName.success") {
                print STDERR "$wrk/3-ovlcorr/$jobName failed.\n";
                $failedJobs++;
            }

            $frgBeg = $frgEnd + 1;
        }

        if ($failedJobs) {
            print STDERR "$failedJobs failed.  Good luck.\n";
            exit(1);
        }

        ########################################

        if (! -e "$wrk/3-ovlcorr/all-corrections.eratelist") {
            if (runCommand("find $wrk/3-ovlcorr -name \\*.erate -print | sort > $wrk/3-ovlcorr/all-corrections.eratelist")) {
                print STDERR "Failed to generate a list of all the overlap correction files.\n";
                rename "$wrk/3-ovlcorr/all-corrections.eratelist", "$wrk/3-ovlcorr/all-corrections.eratelist.FAILED";
                exit(1);
            }
        }

        ########################################

        if (! -e "$wrk/3-ovlcorr/$asm.erates") {
            my $cmd;
            $cmd  = "$bin/cat-erates ";
            $cmd .= "-L $wrk/3-ovlcorr/all-corrections.eratelist ";
            $cmd .= "-o $wrk/3-ovlcorr/$asm.erates ";
            $cmd .= "> $wrk/3-ovlcorr/cat-erates.out ";
            $cmd .= "2> $wrk/3-ovlcorr/cat-erates.err";
            if (runCommand($cmd)) {
                print STDERR "Failed to concatenate the fragment corrections.\n";
                exit(1);
            }
        }

        ########################################

        my $cmd;
        $cmd  = "$bin/update-erates ";
        $cmd .= "$wrk/$asm.ovlStore ";
        $cmd .= "$wrk/3-ovlcorr/$asm.erates";
        $cmd .= "> $wrk/3-ovlcorr/update-erates.out";
        $cmd .= "2> $wrk/3-ovlcorr/update-erates.err";
        if (runCommand($cmd)) {
            print STDERR "Failed to apply the overlap corrections.\n";
            exit(1);
        }

        touch("$wrk/3-ovlcorr/$asm.eratesupdated");
    }
}

1;
