use strict;
#use scheduler;

#  Use the fragment correction results to update the overlap store.

sub createOverlapCorrectionJobs {

    return if (getGlobal("doFragmentCorrection") == 0);

    my $ovlCorrBatchSize    = getGlobal("ovlCorrBatchSize");
    my $ovlCorrOnGrid       = getGlobal("ovlCorrOnGrid") && getGlobal("useGrid");
    my $scratch             = getGlobal("scratch");

    system("mkdir $wrk/3-ovlcorr") if (! -d "$wrk/3-ovlcorr");

    return if (-e "$wrk/3-ovlcorr/jobsCreated.success");

    if ($ovlCorrOnGrid) {
        open(SUB, "> $wrk/3-ovlcorr/submit.sh") or die;
        print SUB "#!/bin/sh\n\n";
    }

    my $frgBeg = 1;
    my $frgEnd = 0;

    while ($frgBeg < $numFrags) {
        $frgEnd = $frgBeg + $ovlCorrBatchSize - 1;
        $frgEnd = $numFrags if ($frgEnd > $numFrags);

        my $jobName   = substr("0000000000" . $frgBeg, -8);

        if ($ovlCorrOnGrid) {
            #  Run on the grid

            open(F, "> $wrk/3-ovlcorr/$jobName.sh") or die;
            print F "#!/bin/sh\n\n";
            print F "jid=\$\$\n";
            #print F "$processStats \\\n";
            print F "$gin/correct-olaps \\\n";
            print F "  -S $wrk/$asm.ovlStore \\\n";
            print F "  -e $scratch/ovlcorr-$jobName.\$jid.erate \\\n";
            print F "  $wrk/$asm.frgStore \\\n";
            print F "  $wrk/2-frgcorr/$asm.corr \\\n";
            print F "  $frgBeg $frgEnd \\\n";
            print F " > $wrk/3-ovlcorr/$jobName.err 2>&1 \\\n";
            print F "&&  \\\n";
            print F "mv $scratch/ovlcorr-$jobName.\$jid.erate \\\n";
            print F "   $wrk/3-ovlcorr/$jobName.erate \\\n";
            print F "&& \\\n";
            print F "touch $wrk/3-ovlcorr/$jobName.success\n";
            close(F);

            chmod 0755, "$wrk/3-ovlcorr/$jobName.sh";

            my $sge                  = getGlobal("sge");
            my $sgeOverlapCorrection = getGlobal("sgeOverlapCorrection");

            print SUB "qsub $sge $sgeOverlapCorrection ";
            print SUB " -r y -N ovc_${asm}_$jobName ";
            print SUB " -j y -o $wrk/3-ovlcorr/$jobName.grid.err ";
            print SUB "$wrk/3-ovlcorr/$jobName.sh\n";
        } else {
            #  Run right here, right now

            open(F, "> $wrk/3-ovlcorr/$jobName.sh") or die;
            print F "#!/bin/sh\n\n";
            #print F "$processStats \\\n";
            print F "$bin/correct-olaps \\\n";
            print F "  -S $wrk/$asm.ovlStore \\\n";
            print F "  -e $wrk/3-ovlcorr/$jobName.erate \\\n";
            print F "  $wrk/$asm.frgStore \\\n";
            print F "  $wrk/2-frgcorr/$asm.corr \\\n";
            print F "  $frgBeg $frgEnd \\\n";
            print F " > $wrk/3-ovlcorr/$jobName.err 2>&1 \\\n";
            print F "&& \\\n";
            print F "touch $wrk/3-ovlcorr/$jobName.success\n";
            close(F);

            &scheduler::schedulerSubmit("sh $wrk/3-ovlcorr/$jobName.sh > /dev/null 2>&1");
        }

        $frgBeg = $frgEnd + 1;
    }

    if ($ovlCorrOnGrid) {
        close(SUB);
        pleaseExecute("$wrk/3-ovlcorr/submit.sh");
        exit(0);
    } else {
        &scheduler::schedulerSetNumberOfProcesses($global{"ovlCorrConcurrency"});
        &scheduler::schedulerFinish();
    }

    touch("$wrk/3-ovlcorr/jobsCreated.success");
}

1;
