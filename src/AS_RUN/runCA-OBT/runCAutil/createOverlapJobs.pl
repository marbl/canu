use strict;


sub createOverlapJobs($) {
    my $isTrim = shift @_;

    return if (-d "$wrk/$asm.ovlStore");

    caFailure("overlapper detected no fragments", undef) if ($numFrags == 0);
    caFailure("overlapper needs to know if trimming or assembling", undef) if (!defined($isTrim));

    my $ovlThreads        = getGlobal("ovlThreads");
    my $ovlMemory         = getGlobal("ovlMemory");

    my $outDir  = "1-overlapper";
    my $ovlOpt  = "";
    my $merSize = getGlobal("ovlMerSize");
    my $merComp = getGlobal("merCompression");

    if ($isTrim eq "trim") {
        $outDir  = "0-overlaptrim-overlap";
        $ovlOpt  = "-G";
        $merSize = getGlobal("obtMerSize");
    }

    system("mkdir $wrk/$outDir") if (! -d "$wrk/$outDir");

    return if (-e "$wrk/$outDir/overlap.sh");

    #  umd overlapper here
    #
    if (getGlobal("ovlOverlapper") eq "umd") {
        #  For Sergey:
        #
        #  UMDoverlapper() needs to dump the gkpstore, run UMD, build
        #  the ovlStore and update gkpStore with new clear ranges.
        #  The explicit call to UMDoverlapper in main() can then go away.
        #  OBT is smart enough to disable itself if umd is enabled.
        #
        UMDoverlapper();
        return;
    }

    #  mer overlapper here
    #
    if ((($isTrim eq "trim") && (getGlobal("obtOverlapper") eq "mer")) ||
        (($isTrim ne "trim") && (getGlobal("ovlOverlapper") eq "mer"))) {
        merOverlapper($isTrim);
        return;
    }

    #  To prevent infinite loops -- stop now if the overlap script
    #  exists.  This will unfortunately make restarting from transient
    #  failures non-trivial.
    #
    #  FAILUREHELPME
    #
    caFailure("overlapper failed\nmanual restart needed to prevent infinite loops\nremove file '$wrk/$outDir/overlap.sh'", undef) if (-e "$wrk/$outDir/overlap.sh");

    meryl();

    #  We make a giant job array for this -- we need to know hashBeg,
    #  hashEnd, refBeg and refEnd -- from that we compute batchName
    #  and jobName.
    #
    #  ovlopts.pl returns the batch name ($batchName), the job name
    #  ($jobName) and options to pass to overlap (-h $hashBeg-$hashEnd
    #  -r $refBeg-$refEnd).  From those, we can construct the command
    #  to run.
    #
    open(F, "> $wrk/$outDir/overlap.sh") or caFailure("can't open '$wrk/$outDir/overlap.sh'", undef);
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
    print F "bat=`head -n \$jobid $wrk/$outDir/ovlbat | tail -n 1`\n";
    print F "job=`head -n \$jobid $wrk/$outDir/ovljob | tail -n 1`\n";
    print F "opt=`head -n \$jobid $wrk/$outDir/ovlopt | tail -n 1`\n";
    print F "jid=\$\$\n";
    print F "\n";
    print F "if [ ! -d $wrk/$outDir/\$bat ]; then\n";
    print F "  mkdir $wrk/$outDir/\$bat\n";
    print F "fi\n";
    print F "\n";
    print F "if [ -e $wrk/$outDir/\$bat/\$job.ovb.gz ]; then\n";
    print F "  echo Job previously completed successfully.\n";
    print F "  exit\n";
    print F "fi\n";
    print F "\n";
    print F "if [ x\$bat = x ]; then\n";
    print F "  echo Error: Job index out of range.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "\n";
    print F "AS_OVL_ERROR_RATE=", getGlobal("ovlErrorRate"), "\n";
    print F "AS_CNS_ERROR_RATE=", getGlobal("cnsErrorRate"), "\n";
    print F "AS_CGW_ERROR_RATE=", getGlobal("cgwErrorRate"), "\n";
    print F "export AS_OVL_ERROR_RATE AS_CNS_ERROR_RATE AS_CGW_ERROR_RATE\n";

    print F getBinDirectoryShellCode();

    print F "\$bin/overlap $ovlOpt -M $ovlMemory -t $ovlThreads \\\n";
    print F "  \$opt \\\n";
    print F "  -k $merSize \\\n";
    print F "  -k $wrk/0-mercounts/$asm.nmers.obt.fasta \\\n" if ($isTrim eq "trim");
    print F "  -k $wrk/0-mercounts/$asm.nmers.ovl.fasta \\\n" if ($isTrim ne "trim");
    print F "  -o $wrk/$outDir/\$bat/\$job.ovb.WORKING.gz \\\n";
    print F "  $wrk/$asm.gkpStore \\\n";
    print F "&& \\\n";
    print F "mv $wrk/$outDir/\$bat/\$job.ovb.WORKING.gz $wrk/$outDir/\$bat/\$job.ovb.gz\n";
    print F "\n";
    print F "exit 0\n";
    close(F);

    system("chmod +x $wrk/$outDir/overlap.sh");

    #  We segment the hash into $numFrags / $ovlHashBlockSize pieces,
    #  and the stream into $numFrags / $ovlRefBlockSize pieces.  Put
    #  all runs for the same hash into a subdirectory.

    my ($hashBeg, $hashEnd, $refBeg, $refEnd) = (1, 0, 1, 0);

    my $ovlHashBlockLength = getGlobal("ovlHashBlockLength");
    my $ovlHashBlockSize   = getGlobal("ovlHashBlockSize");
    my $ovlRefBlockSize    = getGlobal("ovlRefBlockSize");

    if (defined($ovlHashBlockLength)) {
        my $bin = getBinDirectory();

        print STDERR "Partitioning overlap jobs by fragment length.  $ovlHashBlockLength bp per block.\n";
        open(FL, "$bin/gatekeeper -nouid -dumpfragments -tabular $wrk/$asm.gkpStore |");
        $_ = <FL>;  #  Header
    }

    #  Saved for output to ovlopts.pl
    my @bat;
    my @job;
    my @opt;

    #  Number of jobs per batch directory
    #
    my $batchMax  = 1000;
    my $batchSize = 0;

    my $batchName = "001";
    my $jobName   = "000001";

    #  hashBeg and hashEnd are INCLUSIVE ranges, not C-style ranges.

    while ($hashBeg < $numFrags) {
        my $maxNumFrags = 0;
        my $maxLength   = 0;

        if (defined($ovlHashBlockLength)) {
            die "Partitioning error hashBeg=$hashBeg hashEnd=$hashEnd\n" if ($hashEnd != $hashBeg - 1);

            my $len = 0;

            while (<FL>) {
                my @v = split '\s+', $_;

                #  Even deleted fragments contribute to the length (by one byte, the terminating zero)
                $len += $v[11] - $v[10] if ($v[6] == 0);
                $len += 1;
                $hashEnd = $v[1];

                last if ($len >= $ovlHashBlockLength);
            }

            if (eof(FL)) {
                if ($hashEnd != $numFrags) {
                    print STDERR "WARNING:  End of initialTrimLog before end of fragments.  hashEnd=$hashEnd numFrags=$numFrags\n";
                    $hashEnd = $numFrags;
                }
            }

            $maxNumFrags = ($maxNumFrags > $hashEnd - $hashBeg + 1) ? ($maxNumFrags) : ($hashEnd - $hashBeg + 1);
            $maxLength   = ($maxLength   > $len)                    ? ($maxLength)   : ($len);

            print STDERR "Batch $batchName/$jobName from $hashBeg to $hashEnd (", $hashEnd - $hashBeg + 1, " fragments) with length $len\n";
        } else {
            $hashEnd = $hashBeg + $ovlHashBlockSize - 1;
            $hashEnd = $numFrags if ($hashEnd > $numFrags);
        }

        $refBeg = 0;
        $refEnd = 0;

        while ($refBeg < $hashEnd) {
            $refEnd = $refBeg + $ovlRefBlockSize - 1;
            $refEnd = $numFrags if ($refEnd > $numFrags);

            push @bat, "$batchName";
            push @job, "$jobName";

            if (defined($ovlHashBlockLength)) {
                push @opt, "-h $hashBeg-$hashEnd  -r $refBeg-$refEnd --hashstrings $maxNumFrags --hashdatalen $maxLength";
            } else {
                push @opt, "-h $hashBeg-$hashEnd  -r $refBeg-$refEnd";
            }

            $refBeg = $refEnd + 1;

            $batchSize++;
            if ($batchSize >= $batchMax) {
                $batchSize = 0;
                $batchName++;
            }

            $jobName++;
        }

        $hashBeg = $hashEnd + 1;
    }

    open(BAT, "> $wrk/$outDir/ovlbat") or caFailure("failed to open '$wrk/$outDir/ovlbat'", undef);
    open(JOB, "> $wrk/$outDir/ovljob") or caFailure("failed to open '$wrk/$outDir/ovljob'", undef);
    open(OPT, "> $wrk/$outDir/ovlopt") or caFailure("failed to open '$wrk/$outDir/ovlopt'", undef);

    foreach my $b (@bat) { print BAT "$b\n"; }
    foreach my $b (@job) { print JOB "$b\n"; }
    foreach my $b (@opt) { print OPT "$b\n"; }

    close(BAT);
    close(JOB);
    close(OPT);

    my $jobs = scalar(@opt);

    print STDERR "Created $jobs overlap jobs.  Last batch '$batchName', last job '$jobName'.\n";

    #  Submit to the grid (or tell the user to do it), or just run
    #  things here
    #
    if (getGlobal("useGrid") && getGlobal("ovlOnGrid")) {
        my $sge        = getGlobal("sge");
        my $sgeName    = getGlobal("sgeName");
        my $sgeOverlap = getGlobal("sgeOverlap");

        $sgeName = "_$sgeName" if (defined($sgeName));

        my $SGE;
        $SGE  = "qsub $sge $sgeOverlap -cwd -N ovl_$asm$sgeName \\\n";
        $SGE .= "  -t 1-$jobs \\\n";
        $SGE .= "  -j y -o $wrk/$outDir/\\\$TASK_ID.out \\\n";
        $SGE .= "  $wrk/$outDir/overlap.sh\n";

	submitBatchJobs($SGE, "ovl_$asm$sgeName");
        exit(0);
    } else {
        for (my $i=1; $i<=$jobs; $i++) {
            my $out = substr("000000" . $i, -6);
            &scheduler::schedulerSubmit("$wrk/$outDir/overlap.sh $i > $wrk/$outDir/$out.out 2>&1");
        }

        &scheduler::schedulerSetNumberOfProcesses(getGlobal("ovlConcurrency"));
        &scheduler::schedulerFinish();
    }
}

1;
