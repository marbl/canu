use strict;
#use scheduler;

#  Use the fragment correction results to update the overlap store.

sub createOverlapCorrectionJobs {
    my $ovlCorrBatchSize    = getGlobal("ovlCorrBatchSize");
    my $ovlCorrOnGrid       = getGlobal("ovlCorrOnGrid");
    my $scratch             = getGlobal("scratch");

    system("mkdir $wrk/3-ovlcorr") if (! -e "$wrk/3-ovlcorr");

    return if (-e "$wrk/3-ovlcorr/jobsCreated.success");

    if ($ovlCorrOnGrid) {
        open(SUB, "> $wrk/3-ovlcorr/submit.sh") or die;
        print SUB "#!/bin/sh\n";
    }

    my $frgBeg = 1;
    my $frgEnd = 0;

    while ($frgBeg < $numFrags) {
        $frgEnd = $frgBeg + $ovlCorrBatchSize - 1;
        $frgEnd = $numFrags if ($frgEnd > $numFrags);

        my $jobName   = substr("0000000000" . $frgBeg, -8);


        if ($ovlCorrOnGrid == 0) {
            #  Run right here, right now

            open(F, "> $wrk/3-ovlcorr/$jobName.sh") or die;
            print F "#!/bin/sh\n";
            #print F "$processStats \\\n";
            print F "$bin/correct-olaps \\\n";
            print F "  -S $wrk/$asm.ovlStore \\\n";
            print F "  -e $wrk/3-ovlcorr/$jobName.erate \\\n";
            print F "  $wrk/$asm.frgStore \\\n";
            print F "  $wrk/2-frgcorr/$asm.corr \\\n";
            print F "  $frgBeg $frgEnd \\\n";
            print F "&& \\\n";
            print F "touch $wrk/3-ovlcorr/$jobName.success\n";
            close(F);

            #runCommand("sh $wrk/3-ovlcorr/$jobName.sh") and die;
            &scheduler::schedulerSubmit("sh $wrk/3-ovlcorr/$jobName.sh");
        } else {
            #  Run on the grid

            open(F, "> $wrk/3-ovlcorr/$jobName.sh") or die;
            print F "#!/bin/sh\n";
            #print F "$processStats \\\n";
            print F "$gin/correct-olaps \\\n";
            print F "  -S $wrk/$asm.ovlStore \\\n";
            print F "  -e $scratch/ovlcorr-$jobName.erate \\\n";
            print F "  $wrk/$asm.frgStore \\\n";
            print F "  $wrk/2-frgcorr/$asm.corr \\\n";
            print F "  $frgBeg $frgEnd \\\n";
            print F "&&  \\\n";
            print F "mv $scratch/ovlcorr-$jobName.erate \\\n";
            print F "   $wrk/3-ovlcorr/$jobName.erate \\\n";
            print F "&& \\\n";
            print F "touch $wrk/3-ovlcorr/$jobName.success\n";
            close(F);

            chmod 0755, "$wrk/3-ovlcorr/$jobName.sh";

            print SUB "qsub ";
            print SUB "-p 0 ";      #  Priority
            print SUB "-r y ";      #  Rerunnable
            print SUB "-N ovc_${asm}_$jobName ";
            print SUB "-o $wrk/3-ovlcorr/$jobName.out ";
            print SUB "-e $wrk/3-ovlcorr/$jobName.err ";
            print SUB "$wrk/3-ovlcorr/$jobName.sh\n";
        }

        $frgBeg = $frgEnd + 1;
    }

    touch("$wrk/3-ovlcorr/jobsCreated.success");

    if ($ovlCorrOnGrid) {
        close(SUB);
        pleaseExecute("$wrk/3-ovlcorr/submit.sh");
        exit(0);
    } else {
        &scheduler::schedulerSetNumberOfProcesses(4);
        &scheduler::schedulerFinish();
    }
}

1;
