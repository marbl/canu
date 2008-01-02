use strict;

sub createPostUnitiggerConsensusJobs (@) {
    my @cgbFiles  = @_;

    if (! -e "$wrk/5-consensus/$asm.partitioned") {

        #  Then, build a partition information file, and do the partitioning.
        #
        open(G, "> $wrk/5-consensus/$asm.partFile") or caFailure("Failed to write '$wrk/5-consensus/$asm.partFile'\n");
        foreach my $f (@cgbFiles) {
            if ($f =~ m/^.*(\d\d\d).cgb$/) {
                my $part = $1;
                open(F, "grep ^mid: $f |") or caFailure("Failed to grep '^mid: $f'\n");
                while (<F>) {
                    print G "$part $1\n" if (m/^mid:(\d+)$/);                        
                }
                close(F);
            } else {
                caFailure("CGB file didn't match ###.cgb!\n");
            }
        }
        close(G);

        my $cmd;
        $cmd  = "$bin/gatekeeper -P ";
        $cmd .= "$wrk/5-consensus/$asm.partFile ";
        $cmd .= "$wrk/$asm.gkpStore ";
        $cmd .= "> $wrk/5-consensus/$asm.partitioned.err 2>&1";
        if (runCommand("$wrk/5-consensus", $cmd)) {
            rename "$wrk/5-consensus/$asm.partFile", "$wrk/5-consensus/$asm.partFile.FAILED";
            caFailure("Failed to partition the fragStore.\n");
        }

        touch "$wrk/5-consensus/$asm.partitioned";
    }

    ########################################
    #
    #  Build consensus jobs for the grid -- this is very similar to that in createConsensusJobs.pl
    #
    #  Create a set of shell scripts to run consensus, one per cgb
    #  batch.  The last batch is not used, the small tests BPW has
    #  tried always as an empty file there.
    #
    my $jobP;
    my $jobs = 0;

    open(F, "> $wrk/5-consensus/consensus.cgi.input") or caFailure("Failed to open '$wrk/5-consensus/consensus.cgi.input'\n");
    foreach my $f (@cgbFiles) {
        print F "$f\n";

        if ($f =~ m/^.*(\d\d\d).cgb/) {
            $jobP .= "$1\t";
            $jobs++;
        } else {
            print STDERR "WARNING: didn't match $f for CGB filename!\n";
        }
    }
    close(F);

    $jobP = join ' ', sort { $a <=> $b } split '\s+', $jobP;

    open(F, "> $wrk/5-consensus/consensus.sh") or caFailure("Can't open '$wrk/5-consensus/consensus.sh'\n");
    print F "#!/bin/sh\n";
    print F "\n";
    print F "jobid=\$SGE_TASK_ID\n";
    print F "if [ x\$jobid = x -o x\$jobid = xundefined ]; then\n";
    print F "  jobid=\$1\n";
    print F "fi\n";
    print F "if [ x\$jobid = x ]; then\n";
    print F "  echo Error: I need SGE_TASK_ID set, or a job index on the command line.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "jobp=`echo $jobP | cut -d' ' -f \$jobid`\n";
    print F "cgbfile=`head -\$jobid < $wrk/5-consensus/consensus.cgi.input | tail -1`\n";
    print F "\n";
    print F "if [ -e $wrk/5-consensus/${asm}_\$jobp.success ] ; then\n";
    print F "  exit 0\n";
    print F "fi\n";
    print F "\n";
    print F "AS_OVL_ERROR_RATE=", getGlobal("ovlErrorRate"), "\n";
    print F "AS_CNS_ERROR_RATE=", getGlobal("cnsErrorRate"), "\n";
    print F "AS_CGW_ERROR_RATE=", getGlobal("cgwErrorRate"), "\n";
    print F "export AS_OVL_ERROR_RATE AS_CNS_ERROR_RATE AS_CGW_ERROR_RATE\n";
    print F "\n";
    print F "echo \\\n";
    print F "$gin/consensus \\\n";
    print F "  -G -U \\\n";
    print F "  -m -S \$jobp \\\n";
    print F "  -o $wrk/5-consensus/${asm}_\$jobp.cgi \\\n";
    print F "  $wrk/$asm.gkpStore \\\n";
    print F "  \$cgbfile \\\n";
    print F " \\> $wrk/5-consensus/${asm}_\$jobp.err 2\\>\\&1\n";
    print F "\n";
    print F "$gin/consensus \\\n";
    print F "  -G -U \\\n";
    print F "  -m -S \$jobp \\\n";
    print F "  -o $wrk/5-consensus/${asm}_\$jobp.cgi \\\n";
    print F "  $wrk/$asm.gkpStore \\\n";
    print F "  \$cgbfile \\\n";
    print F " >> $wrk/5-consensus/${asm}_\$jobp.err 2>&1 \\\n";
    print F "&& \\\n";
    print F "touch $wrk/5-consensus/${asm}_\$jobp.success\n";
    close(F);

    chmod 0755, "$wrk/5-consensus/consensus.sh";

    if (getGlobal("useGrid") && getGlobal("cnsOnGrid")) {
        my $sge          = getGlobal("sge");
        my $sgeConsensus = getGlobal("sgeConsensus");

        my $SGE;
        $SGE  = "qsub $sge $sgeConsensus -r y -N NAME ";
        $SGE .= "-t MINMAX ";
        $SGE .= "-j y -o /dev/null ";
        $SGE .= "$wrk/5-consensus/consensus.sh\n";

	my $numThreads = 1;
        my $waitTag = submitBatchJobs("cns1", $SGE, $jobs, $numThreads);

        if (runningOnGrid()) {
            touch("$wrk/5-consensus/jobsCreated.success");
            submitScript("$waitTag");
            exit(0);
        } else {
            touch("$wrk/5-consensus/jobsCreated.success");
            exit(0);
        }
    } else {
        for (my $i=1; $i<=$jobs; $i++) {
            &scheduler::schedulerSubmit("sh $wrk/5-consensus/consensus.sh $i > /dev/null 2>&1");
        }

        &scheduler::schedulerSetNumberOfProcesses(getGlobal("cnsConcurrency"));
        &scheduler::schedulerFinish();

        touch("$wrk/5-consensus/jobsCreated.success");
    }
}



sub postUnitiggerConsensus (@) {
    my @cgbFiles  = @_;

    system("mkdir $wrk/5-consensus") if (! -d "$wrk/5-consensus");

    goto alldone if (-e "$wrk/5-consensus/consensus.success");

    #
    #  Create and/or run consensus jobs
    #

    createPostUnitiggerConsensusJobs(@cgbFiles) if (! -e "$wrk/5-consensus/jobsCreated.success");

    #
    #  Check that consensus finished properly
    #

    my $failedJobs = 0;
    my @cgbIndices;

    foreach my $f (@cgbFiles) {
        if ($f =~ m/^.*(\d\d\d).cgb$/) {
            push @cgbIndices, $1;
        } else {
            caFailure("Didn't match '$f' for CGB filename!\n");
        }
    }

    foreach my $f (@cgbIndices) {
        if ((! -e "$wrk/5-consensus/${asm}_$f.success") ||
            (! -e "$wrk/5-consensus/${asm}_$f.cgi")) {
            print STDERR "$wrk/5-consensus/$f failed -- no .success or no .cgi!\n";
            $failedJobs++;
        }
    }

    caFailure("$failedJobs consensusAfterUnitigger jobs failed.  Good luck.\n") if ($failedJobs);

    #  All jobs finished.  Remove the partitioning from the gatekeeper
    #  store.  The gatekeeper store is currently (5 Mar 2007) tolerant
    #  of someone asking for a partition that isn't there -- it'll
    #  fallback to the complete store.  So, if you happen to want to
    #  run consensus again, it'll still work, just a little slower.
    #
    open(F, "ls $wrk/$asm.gkpStore |");
    while (<F>) {
        chomp;
        if (m/^\S\S\S\.\d\d\d$/) {
            unlink $_;
        }
    }
    close(F);


    #
    #  Consolidate all the output
    #

    open(G, "> $wrk/5-consensus/$asm.cgi") or caFailure("Failed to write '$wrk/5-consensus/$asm.cgi'\n");
    foreach my $fid (@cgbIndices) {
        open(F, "< $wrk/5-consensus/${asm}_$fid.cgi") or caFailure("Failed to open '$wrk/5-consensus/${asm}_$fid.cgi'\n");
        while (<F>) {
            print G $_;
        }
        close(F);
    }
    close(G);

    touch ("$wrk/5-consensus/consensus.success");

  alldone:
    stopAfter("consensusAfterUnitigger");
}

1;
