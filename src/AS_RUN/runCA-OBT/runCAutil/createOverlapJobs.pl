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
    print F "bat=`\$perl $wrk/$outDir/ovlopts.pl bat \$jobid`\n";
    print F "job=`\$perl $wrk/$outDir/ovlopts.pl job \$jobid`\n";
    print F "opt=`\$perl $wrk/$outDir/ovlopts.pl opt \$jobid`\n";
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

    my ($hashBeg, $hashEnd, $refBeg, $refEnd) = (getGlobal("ovlStart"), 0, 1, 0);

    my $ovlHashBlockSize  = getGlobal("ovlHashBlockSize");
    my $ovlRefBlockSize   = getGlobal("ovlRefBlockSize");

    #  Saved for output to ovlopts.pl
    my @bat;
    my @job;
    my @opt;

    #  Number of jobs per batch directory
    #
    my $batchMax  = 200;
    my $batchSize = 0;
    my $batch     = 1;

    my $batchName = substr("0000000000" . $batch, -10);

    while ($hashBeg < $numFrags) {
        $hashEnd = $hashBeg + $ovlHashBlockSize - 1;
        $hashEnd = $numFrags if ($hashEnd > $numFrags);
        $refBeg = 0;
        $refEnd = 0;

        while ($refBeg < $hashEnd) {
            $refEnd = $refBeg + $ovlRefBlockSize - 1;
            $refEnd = $numFrags if ($refEnd > $numFrags);

            #print STDERR "hash: $hashBeg-$hashEnd ref: $refBeg-$refEnd\n";

            my $jobName;
            $jobName .= "h" . substr("0000000000" . $hashBeg, -10);
            $jobName .= "r" . substr("0000000000" . $refBeg, -10);

            push @bat, "$batchName";
            push @job, "$jobName";
            push @opt, "-h $hashBeg-$hashEnd  -r $refBeg-$refEnd";

            $refBeg = $refEnd + 1;

            $batchSize++;
            if ($batchSize >= $batchMax) {
                $batch++;
                $batchName = substr("0000000000" . $batch, -10);
                $batchSize = 0;
            }
        }

        $hashBeg = $hashEnd + 1;
    }

    open(SUB, "> $wrk/$outDir/ovlopts.pl") or caFailure("failed to open '$wrk/$outDir/ovlopts.pl'", undef);
    print SUB "#!/usr/bin/env perl\n";
    print SUB "use strict;\n";
    print SUB "my \@bat = (\n";  foreach my $b (@bat) { print SUB "\"$b\",\n"; }  print SUB ");\n";
    print SUB "my \@job = (\n";  foreach my $b (@job) { print SUB "\"$b\",\n"; }  print SUB ");\n";
    print SUB "my \@opt = (\n";  foreach my $b (@opt) { print SUB "\"$b\",\n"; }  print SUB ");\n";
    print SUB "my \$idx = int(\$ARGV[1]) - 1;\n";
    print SUB "if      (\$ARGV[0] eq \"bat\") {\n";
    print SUB "    print \"\$bat[\$idx]\";\n";
    print SUB "} elsif (\$ARGV[0] eq \"job\") {\n";
    print SUB "    print \"\$job[\$idx]\";\n";
    print SUB "} elsif (\$ARGV[0] eq \"opt\") {\n";
    print SUB "    print \"\$opt[\$idx]\";\n";
    print SUB "} else {\n";
    print SUB "    print STDOUT \"Got '\$ARGV[0]' and don't know what to do!\\n\";\n";
    print SUB "    print STDERR \"Got '\$ARGV[0]' and don't know what to do!\\n\";\n";
    print SUB "    die;\n";
    print SUB "}\n";
    print SUB "exit(0);\n";
    close(SUB);

    open(SUB, "> $wrk/$outDir/ovljobs.dat") or caFailure("failed to open '$wrk/$outDir/ovljobs.dat'", undef);
    foreach my $b (@bat) { print SUB "$b "; }  print SUB "\n";
    foreach my $b (@job) { print SUB "$b "; }  print SUB "\n";
    close(SUB);


    my $jobs = scalar(@opt);

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
        $SGE .= "  -j y -o $wrk/$outDir/overlap.\\\$TASK_ID.out \\\n";
        $SGE .= "  $wrk/$outDir/overlap.sh\n";

	submitBatchJobs($SGE, "ovl_$asm$sgeName");
        exit(0);
    } else {
        for (my $i=1; $i<=$jobs; $i++) {
            my $out = substr("0000" . $i, -4);
            &scheduler::schedulerSubmit("$wrk/$outDir/overlap.sh $i > $wrk/$outDir/overlap.$out.out 2>&1");
        }

        &scheduler::schedulerSetNumberOfProcesses(getGlobal("ovlConcurrency"));
        &scheduler::schedulerFinish();
    }
}

1;
