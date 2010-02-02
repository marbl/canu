use strict;

#  Prepare for consensus on the grid
#    Partition the contigs
#    Repartition the frag store

sub createPostScaffolderConsensusJobs () {
    my $consensusType = getGlobal("consensus");

    return if (-e "$wrk/8-consensus/consensus.sh");

    caFailure("contig consensus didn't find '$wrk/$asm.tigStore'", undef)    if (! -d "$wrk/$asm.tigStore");

    my $tigVersion = findLastCheckpoint("$wrk/7-CGW");
    caFailure("contig consensus didn't find any checkpoints in '$wrk/7-CGW'", undef) if (!defined($tigVersion));

    ########################################
    #
    #  Partition the gkpStore for consensus.
    #
    if (! -e "$wrk/8-consensus/$asm.partitioned") {
        my $bin = getBinDirectory();
        my $cmd;
        $cmd  = "$bin/gatekeeper -P $wrk/7-CGW/$asm.partitioning $wrk/$asm.gkpStore ";
        $cmd .= "> $wrk/8-consensus/$asm.partitioned.err 2>&1";

        caFailure("gatekeeper partitioning failed", "$wrk/8-consensus/$asm.partitioned.err") if (runCommand("$wrk/8-consensus", $cmd));
        touch("$wrk/8-consensus/$asm.partitioned");
    }

    ########################################
    #
    #  Build consensus jobs for the grid -- this is very similar to that in createPostUnitiggerConsensus.pl
    #
    my $jobs = 0;

    open(F, "< $wrk/7-CGW/$asm.partitionInfo") or caFailure("can't open '$wrk/7-CGW/$asm.partitionInfo'", undef);
    while (<F>) {
        if (m/Partition\s+(\d+)\s+has\s+(\d+)\s+contigs\sand\s+(\d+)\s+fragments./) {
            $jobs = $1;
        }
    }
    close(F);

    open(F, "> $wrk/8-consensus/consensus.sh") or caFailure("can't open '$wrk/8-consensus/consensus.sh'", undef);
    print F "#!" . getGlobal("shell") . "\n";
    print F "\n";
    print F "jobid=\$SGE_TASK_ID\n";
    print F "if [ x\$jobid = x -o x\$jobid = xundefined ]; then\n";
    print F "  jobid=\$1\n";
    print F "fi\n";
    print F "if [ x\$jobid = x ]; then\n";
    print F "  echo Error: I need SGE_TASK_ID set, or a job index on the command line.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "if [ \$jobid -gt $jobs ]; then\n";
    print F "  echo Error: Only $jobs partitions, you asked for \$jobid.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "\n";
    print F "jobid=`printf %03d \$jobid`\n";
    print F "\n";
    print F "if [ -e $wrk/8-consensus/${asm}_\$jobid.success ] ; then\n";
    print F "  exit 0\n";
    print F "fi\n";
    print F "\n";
    print F "AS_OVL_ERROR_RATE=", getGlobal("ovlErrorRate"), "\n";
    print F "AS_CNS_ERROR_RATE=", getGlobal("cnsErrorRate"), "\n";
    print F "AS_CGW_ERROR_RATE=", getGlobal("cgwErrorRate"), "\n";
    print F "export AS_OVL_ERROR_RATE AS_CNS_ERROR_RATE AS_CGW_ERROR_RATE\n";

    print F getBinDirectoryShellCode();

    if ($consensusType eq "cns") {
       print F "\$bin/ctgcns \\\n";
       print F "  -g $wrk/$asm.gkpStore \\\n";
       print F "  -t $wrk/$asm.tigStore $tigVersion \$jobid \\\n";
       print F "  -P ", getGlobal("cnsPhasing"), "\\\n";
       print F " > $wrk/8-consensus/${asm}_\$jobid.err 2>&1 \\\n";
       print F "&& \\\n";
       print F "touch $wrk/8-consensus/${asm}_\$jobid.success\n";
    } elsif ($consensusType eq "seqan") {
       print F "\$bin/SeqAn_CNS \\\n";
       print F "  -G $wrk/$asm.gkpStore \\\n";
       print F "  -u $wrk/$asm.SeqStore \\\n";
       print F "  -V $tigVersion \\\n";
       print F "  -p \$jobid \\\n";
       print F "  -S \$jobid \\\n";
       #print F "  -c $cgwDir/$asm.cgw_contigs.\$jobid \\\n";
       print F "  -s \$bin/graph_consensus \\\n";
       print F "  -w $wrk/8-consensus/ \\\n";
       print F "  -o $wrk/8-consensus/$asm.cns_contigs.\$jobid \\\n";
       print F " > $wrk/8-consensus/$asm.cns_contigs.\$jobid.err 2>&1 \\\n";
       print F "&& \\\n";
       print F "touch $wrk/8-consensus/$asm.cns_contigs.\$jobid.success\n";
    } else {
       caFailure("unknown consensus type $consensusType; must be 'cns' or 'seqan'", undef);
    }
    print F "exit 0\n";
    close(F);

    chmod 0755, "$wrk/8-consensus/consensus.sh";

    if (getGlobal("cnsOnGrid") && getGlobal("useGrid")) {
        my $sge          = getGlobal("sge");
        my $sgeName      = getGlobal("sgeName");
        my $sgeConsensus = getGlobal("sgeConsensus");

        $sgeName = "_$sgeName" if (defined($sgeName));

        my $SGE;
        $SGE  = "qsub $sge $sgeConsensus -cwd -N ctg_$asm$sgeName ";
        $SGE .= "-t 1-$jobs ";
        $SGE .= "-j y -o /dev/null ";
        $SGE .= "$wrk/8-consensus/consensus.sh\n";

        submitBatchJobs($SGE, "ctg_$asm$sgeName");
        exit(0);
    } else {
        for (my $i=1; $i<=$jobs; $i++) {
            &scheduler::schedulerSubmit("$wrk/8-consensus/consensus.sh $i > /dev/null 2>&1");
        }

        &scheduler::schedulerSetNumberOfProcesses(getGlobal("cnsConcurrency"));
        &scheduler::schedulerFinish();
    }
}


sub postScaffolderConsensus () {

    system("mkdir $wrk/8-consensus") if (! -d "$wrk/8-consensus");

    goto alldone if (-e "$wrk/8-consensus/consensus.success");

    createPostScaffolderConsensusJobs();

    #
    #  Check that consensus finished properly
    #
    my $failedJobs = 0;

    open(F, "< $wrk/7-CGW/$asm.partitionInfo") or caFailure("can't open '$wrk/7-CGW/$asm.partitionInfo'", undef);
    while (<F>) {
        if (m/Partition\s+(\d+)\s+has\s+(\d+)\s+contigs\sand\s+(\d+)\s+fragments./) {
            my $id = substr("000" . $1, -3);

            if (! -e "$wrk/8-consensus/${asm}_$id.success") {
                print STDERR "$wrk/8-consensus/${asm}_$id failed -- no .success.\n";
                $failedJobs++;
            }
        }
    }
    close(F);

    #  FAILUREHELPME
    #
    caFailure("$failedJobs consensusAfterScaffolder jobs failed; remove 8-consensus/consensus.sh to try again", undef) if ($failedJobs);

    #  All jobs finished.  Remove the partitioning from the gatekeeper
    #  store.  The gatekeeper store is currently (5 Mar 2007) tolerant
    #  of someone asking for a partition that isn't there -- it'll
    #  fallback to the complete store.  So, if you happen to want to
    #  run consensus again, it'll still work, just a little slower.
    #
    #  (This block appears in both createPostUnitiggerConsensus.pl and createConsensusJobs.pl)
    #
    system("rm -f $wrk/$asm.gkpStore/frg.[0-9][0-9][0-9]");
    system("rm -f $wrk/$asm.gkpStore/hps.[0-9][0-9][0-9]");
    system("rm -f $wrk/$asm.gkpStore/qlt.[0-9][0-9][0-9]");
    system("rm -f $wrk/$asm.gkpStore/src.[0-9][0-9][0-9]");

    touch("$wrk/8-consensus/consensus.success");

  alldone:
    stopAfter("ctgcns");
    stopAfter("consensusAfterScaffolder");
}

1;
