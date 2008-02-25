use strict;

sub createFragmentCorrectionJobs {
    my $batchSize   = getGlobal("frgCorrBatchSize");
    my $numThreads  = getGlobal("frgCorrThreads");
    my $scratch     = getGlobal("scratch");

    return if (getGlobal("doFragmentCorrection") == 0);
    return if (-e "$wrk/2-frgcorr/jobsCreated.success");
    system("mkdir $wrk/2-frgcorr") if (! -e "$wrk/2-frgcorr");

    #  Figure out how many jobs there are
    my $jobs = int($numFrags / ($batchSize-1)) + 1;

    open(F, "> $wrk/2-frgcorr/correct.sh") or caFailure("Failed to write to '$wrk/2-frgcorr/correct.sh'\n");
    print F "#!/bin/sh\n\n";
    print F "jobid=\$SGE_TASK_ID\n";
    print F "if [ x\$jobid = x -o x\$jobid = xundefined ]; then\n";
    print F "  jobid=\$1\n";
    print F "fi\n";
    print F "if [ x\$jobid = x ]; then\n";
    print F "  echo Error: I need SGE_TASK_ID set, or a job index on the command line.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "\n";
    print F "jobid=`printf %04d \$jobid`\n";
    print F "minid=`expr \$jobid \\* $batchSize - $batchSize + 1`\n";
    print F "maxid=`expr \$jobid \\* $batchSize`\n";
    print F "runid=\$\$\n";
    print F "\n";
    print F "if [ \$maxid -gt $numFrags ] ; then\n";
    print F "  maxid=$numFrags\n";
    print F "fi\n";
    print F "if [ \$minid -gt \$maxid ] ; then\n";
    print F "  echo Job partitioning error -- minid=\$minid maxid=\$maxid.\n";
    print F "  exit\n";
    print F "fi\n";
    print F "\n";
    print F "AS_OVL_ERROR_RATE=", getGlobal("ovlErrorRate"), "\n";
    print F "AS_CNS_ERROR_RATE=", getGlobal("cnsErrorRate"), "\n";
    print F "AS_CGW_ERROR_RATE=", getGlobal("cgwErrorRate"), "\n";
    print F "export AS_OVL_ERROR_RATE AS_CNS_ERROR_RATE AS_CGW_ERROR_RATE\n";
    print F "\n";
    print F "if [ -e $wrk/2-frgcorr/\$jobid.success ] ; then\n";
    print F "  echo Job previously completed successfully.\n";
    print F "  exit\n";
    print F "fi\n";

    print F getBinDirectoryShellCode();

    print F "\$bin/correct-frags \\\n";
    print F "  -t $numThreads \\\n";
    print F "  -S $wrk/$asm.ovlStore \\\n";
    print F "  -o $wrk/2-frgcorr/\$jobid.frgcorr \\\n";
    print F "  $wrk/$asm.gkpStore \\\n";
    print F "  \$minid \$maxid \\\n";
    print F "&& \\\n";
    print F "touch $wrk/2-frgcorr/\$jobid.success\n";

    close(F);

    chmod 0755, "$wrk/2-frgcorr/correct.sh";

    if (getGlobal("frgCorrOnGrid") && getGlobal("useGrid")) {
        #  Run the correction job on the grid.

        my $sge                   = getGlobal("sge");
        my $sgeFragmentCorrection = getGlobal("sgeFragmentCorrection");

        my $SGE;
        $SGE  = "qsub $sge $sgeFragmentCorrection -r y -N NAME ";
        $SGE .= "-t MINMAX ";
        $SGE .= " -j y -o $wrk/2-frgcorr/\\\$TASK_ID.err ";
        $SGE .= "$wrk/2-frgcorr/correct.sh\n";

        my $waitTag = submitBatchJobs("frg", $SGE, $jobs, $numThreads);

        if (runningOnGrid()) {
            touch("$wrk/2-frgcorr/jobsCreated.success");
            submitScript("$waitTag");
            exit(0);
        } else {
            touch("$wrk/2-frgcorr/jobsCreated.success");
            exit(0);
        }
    } else {
        #  Run the correction job right here, right now.

        for (my $i=1; $i<=$jobs; $i++) {
            my $out = substr("0000" . $i, -4);
            &scheduler::schedulerSubmit("sh $wrk/2-frgcorr/correct.sh $i > $wrk/2-frgcorr/$out.err 2>&1");
        }

        &scheduler::schedulerSetNumberOfProcesses($global{"frgCorrConcurrency"});
        &scheduler::schedulerFinish();
    }

    touch("$wrk/2-frgcorr/jobsCreated.success");
    
}

1;
