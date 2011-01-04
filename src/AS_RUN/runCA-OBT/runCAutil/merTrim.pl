use strict;



sub findMBTFailures ($) {
    my $mbtJobs  = shift @_;
    my $failures = 0;

    for (my $i=1; $i<=$mbtJobs; $i++) {
        my $jobid = substr("0000" . $i, -4);
        if (-e "$wrk/0-mertrim/$asm.merTrim.\$jobid.log.WORKING") {
            $failures++;
        }
    }

    return $failures;
}

sub findMBTSuccess ($) {
    my $mbtJobs  = shift @_;
    my $successes= 0;

    for (my $i=1; $i<=$mbtJobs; $i++) {
        my $jobid = substr("0000" . $i, -4);
        if (-e "$wrk/0-mertrim/$asm.merTrim.\$jobid.log") {
            $successes++;
        }
    }

    return($successes == $mbtJobs);
}




sub merTrim {

    #return if (getGlobal("doMerBasedTrimming") == 0);
    return if (getGlobal("doOverlapBasedTrimming") == 0);
    return if (getGlobal("ovlOverlapper") eq "umd");

    #  Skip mer based trimming if it is done, or if the ovlStore already exists.
    #
    goto alldone if (-e "$wrk/0-mertrim/mertrim.success");
    goto alldone if (-d "$wrk/$asm.ovlStore");

    system("mkdir $wrk/0-mertrim")         if (! -d "$wrk/0-mertrim");

    my $bin = getBinDirectory();

    #  Decide if any libraries request mer based trimming.  There is a simpler method to get
    #  this (see unitigger.pl), but for unity with overlapTrim.pl, we count the number
    #  of libraries that want MBT.
    #
    my $mbtNeeded = 0;

    open(F, "$bin/gatekeeper -nouid -dumplibraries $wrk/$asm.gkpStore |");
    while (<F>) {
        $mbtNeeded++ if (m/doMerBasedTrimming.*=.*1/);
    }
    close(F);

    if ($mbtNeeded == 0) {
        touch("$wrk/0-mertrim/mertrim.success");
        goto alldone;
    }

    #
    #  Run mer trim on the grid.
    #
    #  Well, we'd LIKE to run it on the grid, but can't.  merTrim updates the gkpStore
    #  after each fragment it processes, so we can't (reliably) have multiple processes
    #  running at the same time.  We could make merTrim write a log of changes, then apply
    #  then when all processes are done.  I'm not sure there will be any benefit -- we'll still
    #  have one giant process that needs to read/write the whole gkpStore.
    #

    meryl();

    my $mbtBatchSize = getGlobal("mbtBatchSize");
    my $mbtJobs      = int($numFrags / $mbtBatchSize) + (($numFrags % $mbtBatchSize == 0) ? 0 : 1);

    my $merSize      = getGlobal("obtMerSize");
    my $merComp      = 0;  # getGlobal("merCompression");

    if (! -e "$wrk/0-mertrim/mertrim.sh") {
        open(F, "> $wrk/0-mertrim/mertrim.sh") or caFailure("can't open '$wrk/0-mertrim/mertrim.sh'", undef);
        print F "#!" . getGlobal("shell") . "\n";
        print F "\n";
        print F "perl='/usr/bin/env perl'\n";
        print F "\n";
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
        print F "minid=`expr \$jobid \\* $mbtBatchSize - $mbtBatchSize + 1`\n";
        print F "maxid=`expr \$jobid \\* $mbtBatchSize`\n";
        print F "\n";
        print F "if [ \$maxid -gt $numFrags ] ; then\n";
        print F "  maxid=$numFrags\n";
        print F "fi\n";
        print F "\n";
        print F "if [ \$minid -gt \$maxid ] ; then\n";
        print F "  echo Job partitioning error -- minid=\$minid maxid=\$maxid.\n";
        print F "  exit\n";
        print F "fi\n";
        print F "\n";
        print F "if [ -e $wrk/0-mertrim/$asm.merTrim.\$jobid.log ]; then\n";
        print F "  echo Job previously completed successfully.\n";
        print F "  exit\n";
        print F "fi\n";
        print F "\n";

        print F "AS_OVL_ERROR_RATE=", getGlobal("ovlErrorRate"), "\n";
        print F "AS_CNS_ERROR_RATE=", getGlobal("cnsErrorRate"), "\n";
        print F "AS_CGW_ERROR_RATE=", getGlobal("cgwErrorRate"), "\n";
        print F "export AS_OVL_ERROR_RATE AS_CNS_ERROR_RATE AS_CGW_ERROR_RATE\n";

        print F getBinDirectoryShellCode();

        print F "\$bin/merTrim \\\n";
        print F " -b  \$minid \\\n";
        print F " -e  \$maxid \\\n";
        print F " -g  $wrk/$asm.gkpStore \\\n";
        print F " -m  $merSize \\\n";
        print F " -c  $merComp \\\n";
        print F " -mc $wrk/0-mercounts/$asm-C-ms$merSize-cm$merComp \\\n";
        print F " -l  $wrk/0-mertrim/$asm.merTrim.\$jobid.log.WORKING \\\n";
        print F " >   $wrk/0-mertrim/$asm.merTrim.\$jobid.err 2>&1 \\\n";
        print F "&& \\\n";
        print F "mv $wrk/0-mertrim/$asm.merTrim.\$jobid.log.WORKING $wrk/0-mertrim/$asm.merTrim.\$jobid.log\n";
        print F "\n";
        print F "exit 0\n";
        close(F);

        system("chmod +x $wrk/0-mertrim/mertrim.sh");
    }

    stopBefore("initialTrim", undef);

    #  Don't try to rerun failures.
    #
    #  FAILUREHELPME
    #
    if (findMBTFailures($mbtJobs) > 0) {
        caFailure("merTrim failed.  See *.err in $wrk/0-mertrim", undef);
    }

    #  Submit to the grid (or tell the user to do it), or just run
    #  things here
    #
    if (findMBTSuccess($mbtJobs) == 0) {
        if (getGlobal("useGrid") && getGlobal("mbtOnGrid")) {
            my $sge        = getGlobal("sge");
            my $sgeName    = getGlobal("sgeName");
            my $sgeMerTrim = getGlobal("sgeMerTrim");

            $sgeName = "_$sgeName" if (defined($sgeName));

            my $SGE;
            $SGE  = "qsub $sge $sgeMerTrim -cwd -N mbt_$asm$sgeName \\\n";
            $SGE .= "  -t 1-$mbtJobs \\\n";
            $SGE .= "  -j y -o $wrk/0-mertrim/$asm.merTrim.\$TASK_ID.err \\\n";
            $SGE .= "  $wrk/0-mertrim/mertrim.sh\n";

            submitBatchJobs($SGE, "mbt_$asm$sgeName");
            exit(0);
        } else {
            for (my $i=1; $i<=$mbtJobs; $i++) {
                my $out = substr("0000" . $i, -4);
                &scheduler::schedulerSubmit("$wrk/0-mertrim/mertrim.sh $i > $wrk/0-mertrim/$asm.merTrim.err 2>&1");
            }
            
            &scheduler::schedulerSetNumberOfProcesses(getGlobal("mbtConcurrency"));
            &scheduler::schedulerFinish();
        }
    }

    #  Make sure everything finished ok.
    #
    #  FAILUREHELPME
    #
    if (findMBTFailures($mbtJobs) > 0) {
        caFailure("merTrim failed.  See *.err in $wrk/0-mertrim", undef);
    }


    if (runCommand($wrk, "find $wrk/0-mertrim -name \\*.merTrim.\\*.log -print > $wrk/0-mertrim/$asm.merTrim.list")) {
        caFailure("failed to generate a list of all the merTrim results", undef);
    }

    {
        my $cmd;

        $cmd  = "$bin/merTrimApply \\\n";
        $cmd .= " -g $wrk/$asm.gkpStore \\\n";
        $cmd .= " -L $wrk/0-mertrim/$asm.merTrim.list \\\n";
        $cmd .= " -l $wrk/0-mertrim/$asm.merTrimLog \\\n";
        $cmd .= " > $wrk/0-mertrim/$asm.merTrim.err 2>&1";

        if (runCommand($wrk, $cmd)) {
            rename "$wrk/0-mertrim/$asm.merTrimLog", "$wrk/0-mertrim/$asm.merTrimLog.FAILED";
            caError("merTrimApply failed", "$wrk/0-mertrim/$asm.merTrim.err");
        }
    }

    touch("$wrk/0-mertrim/mertrim.success");

  alldone:
    #stopAfter("overlapBasedTrimming");
    #stopAfter("OBT");
}

1;
