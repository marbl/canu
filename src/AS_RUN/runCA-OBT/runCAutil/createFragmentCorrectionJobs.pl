use strict;

sub createFragmentCorrectionJobs {
    my $frgCorrBatchSize  = getGlobal("frgCorrBatchSize");
    my $frgCorrThreads    = getGlobal("frgCorrThreads");
    my $scratch           = getGlobal("scratch");

    return if (getGlobal("doFragmentCorrection") == 0);
    return if (-e "$wrk/2-frgcorr/jobsCreated.success");
    system("mkdir $wrk/2-frgcorr") if (! -e "$wrk/2-frgcorr");

    #  Figure out how many jobs there are
    my $jobs = int($numFrags / ($frgCorrBatchSize-1)) + 1;

    open(F, "> $wrk/2-frgcorr/correct.sh") or die;
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
    print F "if [ \$jobid -gt $jobs ] ; then\n";
    print F "  exit\n";
    print F "fi\n";
    print F "\n";
    print F "frgBeg=`expr \$jobid \\* $frgCorrBatchSize - $frgCorrBatchSize + 1`\n";
    print F "frgEnd=`expr \$jobid \\* $frgCorrBatchSize`\n";
    print F "if [ \$frgEnd -ge $numFrags ] ; then\n";
    print F "  frgEnd=`expr $numFrags - 1`\n";
    print F "fi\n";
    print F "frgBeg=`printf %08d \$frgBeg`\n";
    print F "frgEnd=`printf %08d \$frgEnd`\n";
    print F "\n";
    print F "if [ ! -e $wrk/2-frgcorr/$asm-\$frgBeg-\$frgEnd.success ] ; then\n";
    if (getGlobal('frgCorrOnGrid')) {
        print F "  $gin/correct-frags \\\n";
    } else {
        print F "  $bin/correct-frags \\\n";
    }
    print F "    -S $wrk/$asm.ovlStore \\\n";
    print F "    -o $scratch/$asm-\$frgBeg-\$frgEnd.frgcorr \\\n";
    print F "    $wrk/$asm.frgStore \\\n";
    print F "    \$frgBeg \$frgEnd \\\n";
    print F "   > $wrk/2-frgcorr/$asm-\$frgBeg-\$frgEnd.err 2>&1 \\\n";
    print F "  &&  \\\n";
    print F "  cp -p $scratch/$asm-\$frgBeg-\$frgEnd.frgcorr \\\n";
    print F "        $wrk/2-frgcorr/$asm-\$frgBeg-\$frgEnd.frgcorr \\\n";
    print F "  && \\\n";
    print F "  rm -f $scratch/$asm-\$frgBeg-\$frgEnd.frgcorr \\\n";
    print F "  && \\\n";
    print F "  touch $wrk/2-frgcorr/$asm-\$frgBeg-\$frgEnd.success\n";
    print F "fi\n";
    close(F);

    chmod 0755, "$wrk/2-frgcorr/correct.sh";

    if (getGlobal("frgCorrOnGrid") && getGlobal("useGrid")) {
        #  Run the correction job on the grid.

        my $sge                   = getGlobal("sge");
        my $sgeFragmentCorrection = getGlobal("sgeFragmentCorrection");

        my $SGE;
        $SGE  = "qsub $sge $sgeFragmentCorrection -r y -N frg_$asm ";
        $SGE .= "-t 1-$jobs ";
        $SGE .= " -j y -o /dev/null ";
        $SGE .= "$wrk/2-frgcorr/correct.sh\n";

        if (runningOnGrid()) {
            touch("$wrk/2-frgcorr/jobsCreated.success");
            system($SGE) and die "Failed to submit fragment correction jobs.\n";
            submitScript("frg_$asm");
            exit(0);
        } else {
            pleaseExecute($SGE);
            touch("$wrk/2-frgcorr/jobsCreated.success");
            exit(0);
        }
    } else {
        #  Run the correction job right here, right now.

        for (my $i=1; $i<=$jobs; $i++) {
            &scheduler::schedulerSubmit("sh $wrk/2-frgcorr/correct.sh $i > /dev/null 2>&1");
        }

        &scheduler::schedulerSetNumberOfProcesses($global{"frgCorrConcurrency"});
        &scheduler::schedulerFinish();
    }

    touch("$wrk/2-frgcorr/jobsCreated.success");
}

1;
