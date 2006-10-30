use strict;

sub createOverlapCorrectionJobs {
    my $ovlCorrBatchSize    = getGlobal("ovlCorrBatchSize");
    my $scratch             = getGlobal("scratch");

    return if (getGlobal("doFragmentCorrection") == 0);
    return if (-e "$wrk/3-ovlcorr/jobsCreated.success");
    system("mkdir $wrk/3-ovlcorr") if (! -d "$wrk/3-ovlcorr");

    #  Figure out how many jobs there are
    my $jobs = int($numFrags / ($ovlCorrBatchSize-1)) + 1;

    open(F, "> $wrk/3-ovlcorr/correct.sh") or die;
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
    print F "frgBeg=`expr \$jobid \\* $ovlCorrBatchSize - $ovlCorrBatchSize + 1`\n";
    print F "frgEnd=`expr \$jobid \\* $ovlCorrBatchSize`\n";
    print F "if [ \$frgEnd -ge $numFrags ] ; then\n";
    print F "  frgEnd=`expr $numFrags - 1`\n";
    print F "fi\n";
    print F "frgBeg=`printf %08d \$frgBeg`\n";
    print F "frgEnd=`printf %08d \$frgEnd`\n";
    print F "\n";
    print F "if [ ! -e $wrk/3-ovlcorr/$asm-\$frgBeg-\$frgEnd.success ] ; then\n";
    print F "  $gin/correct-olaps \\\n";
    print F "    -S $wrk/$asm.ovlStore \\\n";
    print F "    -e $scratch/$asm-\$frgBeg-\$frgEnd.erate \\\n";
    print F "    $wrk/$asm.frgStore \\\n";
    print F "    $wrk/2-frgcorr/$asm.corr \\\n";
    print F "    \$frgBeg \$frgEnd \\\n";
    print F "   > $wrk/3-ovlcorr/$asm-\$frgBeg-\$frgEnd.err 2>&1 \\\n";
    print F "  &&  \\\n";
    print F "  cp -p $scratch/$asm-\$frgBeg-\$frgEnd.erate \\\n";
    print F "        $wrk/3-ovlcorr/$asm-\$frgBeg-\$frgEnd.erate \\\n";
    print F "  &&  \\\n";
    print F "  rm -f $scratch/$asm-\$frgBeg-\$frgEnd.erate \\\n";
    print F "  && \\\n";
    print F "  touch $wrk/3-ovlcorr/$asm-\$frgBeg-\$frgEnd.success\n";
    print F "fi\n";
    close(F);

    chmod 0755, "$wrk/3-ovlcorr/correct.sh";

    if (getGlobal("ovlCorrOnGrid") && getGlobal("useGrid")) {
        #  Run the correction job on the grid.

        my $sge                   = getGlobal("sge");
        my $sgeOverlapCorrection  = getGlobal("sgeOverlapCorrection");

        my $SGE;
        $SGE  = "qsub $sge $sgeOverlapCorrection -r y -N ovlcorr_$asm ";
        $SGE .= "-t 1-$jobs ";
        $SGE .= " -j y -o /dev/null ";
        $SGE .= "$wrk/3-ovlcorr/correct.sh\n";

        if (runningOnGrid()) {
            touch("$wrk/3-ovlcorr/jobsCreated.success");
            system($SGE) and die "Failed to submit overlap correction jobs.\n";
            submitScript("ovlcorr_$asm");
            exit(0);
        } else {
            pleaseExecute($SGE);
            touch("$wrk/3-ovlcorr/jobsCreated.success");
            exit(0);
        }
    } else {
        #  Run the correction job right here, right now.

        for (my $i=1; $i<=$jobs; $i++) {
            &scheduler::schedulerSubmit("sh $wrk/3-ovlcorr/correct.sh $i > /dev/null 2>&1");
        }

        &scheduler::schedulerSetNumberOfProcesses($global{"ovlCorrConcurrency"});
        &scheduler::schedulerFinish();
    }

    touch("$wrk/3-ovlcorr/jobsCreated.success");
}

1;
