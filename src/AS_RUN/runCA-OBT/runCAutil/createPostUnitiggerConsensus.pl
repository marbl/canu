use strict;

sub createPostUnitiggerConsensusJobs (@) {
    my @cgbFiles  = @_;
    my $consensusType = getGlobal("consensus");

    return if (-e "$wrk/5-consensus/consensus.sh");

    if (! -e "$wrk/5-consensus/$asm.partitioned") {

        #  Then, build a partition information file, and do the partitioning.
        #
        open(G, "> $wrk/5-consensus/$asm.partFile") or caFailure("failed to write '$wrk/5-consensus/$asm.partFile'", undef);
        foreach my $f (@cgbFiles) {
            if ($f =~ m/^.*(\d\d\d).cgb$/) {
                my $part = $1;
                open(F, "grep ^mid: $f |") or caFailure("failed to grep '^mid: $f'", undef);
                while (<F>) {
                    print G "$part $1\n" if (m/^mid:(\d+)$/);
                }
                close(F);
            } else {
                caFailure("unitigger file '$f' didn't match ###.cgb", undef);
            }
        }
        close(G);

        my $bin = getBinDirectory();
        my $cmd;
        $cmd  = "$bin/gatekeeper -P ";
        $cmd .= "$wrk/5-consensus/$asm.partFile ";
        $cmd .= "$wrk/$asm.gkpStore ";
        $cmd .= "> $wrk/5-consensus/$asm.partitioned.err 2>&1";
        if (runCommand("$wrk/5-consensus", $cmd)) {
            rename "$wrk/5-consensus/$asm.partFile", "$wrk/5-consensus/$asm.partFile.FAILED";
            caFailure("failed to partition the fragStore", "$wrk/5-consensus/$asm.partitioned.err");
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

    open(F, "> $wrk/5-consensus/consensus.cgi.input") or caFailure("failed to open '$wrk/5-consensus/consensus.cgi.input'", undef);
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

    open(F, "> $wrk/5-consensus/consensus.sh") or caFailure("can't open '$wrk/5-consensus/consensus.sh'", undef);
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
    print F "jobp=`echo $jobP | cut -d' ' -f \$jobid`\n";
    print F "cgbfile=`head -n \$jobid < $wrk/5-consensus/consensus.cgi.input | tail -n 1`\n";
    print F "\n";
    print F "if [ -e $wrk/5-consensus/${asm}_\$jobp.success ] ; then\n";
    print F "  exit 0\n";
    print F "fi\n";
    print F "\n";
    print F "AS_OVL_ERROR_RATE=", getGlobal("ovlErrorRate"), "\n";
    print F "AS_CNS_ERROR_RATE=", getGlobal("cnsErrorRate"), "\n";
    print F "AS_CGW_ERROR_RATE=", getGlobal("cgwErrorRate"), "\n";
    print F "export AS_OVL_ERROR_RATE AS_CNS_ERROR_RATE AS_CGW_ERROR_RATE\n";

    print F getBinDirectoryShellCode();

    if ($consensusType eq "cns") {
       print F "\$bin/consensus -U -m \\\n";
       #  Too expensive in general, let the fixUnitigs handle anything that fails here.
       #print F "  -O $wrk/$asm.ovlStore \\\n" if (getGlobal('unitigger') eq "bog");
       print F "  -S \$jobp \\\n";
       print F "  -o $wrk/5-consensus/${asm}_\$jobp.cgi \\\n";
       print F "  $wrk/$asm.gkpStore \\\n";
       print F "  \$cgbfile \\\n";
       print F " > $wrk/5-consensus/${asm}_\$jobp.err 2>&1\n";
       print F "\n";
       print F "if [ -e $wrk/5-consensus/${asm}_\$jobp.cgi_tmp ] ; then\n";
       print F "  echo Yikes!  Consensus crashed.\n";
       print F "  exit 1\n";
       print F "fi\n";
       print F "\n";
       print F "\n";
       print F "#  Attempt to autofix problems.\n";
       print F "if [ -e $wrk/5-consensus/${asm}_\$jobp.cgi.failed ] ; then\n";
       print F "  mv $wrk/5-consensus/${asm}_\$jobp.cgi.failed \\\n";
       print F "     $wrk/5-consensus/${asm}_\$jobp.autofix.orig\n";
       print F "\n";
       print F "  #  Consensus will remove this if successful.\n";
       print F "  touch $wrk/5-consensus/${asm}_\$jobp.autofix.cgi.failed\n";
       print F "\n";
       print F "  \$bin/fixUnitigs -O $wrk/$asm.ovlStore \\\n";
       print F "    < $wrk/5-consensus/${asm}_\$jobp.autofix.orig \\\n";
       print F "    > $wrk/5-consensus/${asm}_\$jobp.autofix \\\n";
       print F "   2> $wrk/5-consensus/${asm}_\$jobp.autofix.log \\\n";
       print F "  && \\\n";
       print F "  \$bin/consensus -U -m \\\n";
       print F "    -D verbosemultialign \\\n";
       print F "    -O $wrk/$asm.ovlStore \\\n";
       print F "    -S \$jobp \\\n";
       print F "    -o $wrk/5-consensus/${asm}_\$jobp.autofix.cgi \\\n";
       print F "    $wrk/$asm.gkpStore \\\n";
       print F "    $wrk/5-consensus/${asm}_\$jobp.autofix \\\n";
       print F "   > $wrk/5-consensus/${asm}_\$jobp.autofix.err 2>&1\n";
       print F "  \n";
       print F "fi\n";
       print F "\n";
       print F "\n";
       print F "\n";
       print F "if [ ! -e $wrk/5-consensus/${asm}_\$jobp.cgi.failed -a \\\n";
       print F "     ! -e $wrk/5-consensus/${asm}_\$jobp.autofix.cgi.failed ] ; then\n";
       print F "  touch $wrk/5-consensus/${asm}_\$jobp.success\n";
       print F "fi\n";
       print F "\n";
    } elsif ($consensusType eq "seqan") {
       print F "\$bin/SeqAn_CNS \\\n";
       print F "  -G $wrk/$asm.gkpStore \\\n";
       print F "  -c \$cgbfile \\\n";
       print F "  -s \$bin/graph_consensus \\\n";
       print F "  -w $wrk/5-consensus/ \\\n";
       print F "  -o $wrk/5-consensus/${asm}_\$jobp.cgi \\\n";
       print F " > $wrk/5-consensus/${asm}_\$jobp.err 2>&1 \\\n";
       print F "&& \\\n";
       print F "touch $wrk/5-consensus/${asm}_\$jobp.success\n";
    } else {
       caFailure("unknown consensus type $consensusType; should be 'cns' or 'seqan'", undef);
    }
    close(F);

    chmod 0755, "$wrk/5-consensus/consensus.sh";

    if (getGlobal("useGrid") && getGlobal("cnsOnGrid")) {
        my $sge          = getGlobal("sge");
        my $sgeConsensus = getGlobal("sgeConsensus");

        my $SGE;
        $SGE  = "qsub $sge $sgeConsensus -cwd -N utg_$asm ";
        $SGE .= "-t 1-$jobs ";
        $SGE .= "-j y -o /dev/null ";
        $SGE .= "$wrk/5-consensus/consensus.sh\n";

        submitBatchJobs($SGE, "utg_$asm");
        exit(0);
    } else {
        for (my $i=1; $i<=$jobs; $i++) {
            &scheduler::schedulerSubmit("$wrk/5-consensus/consensus.sh $i > /dev/null 2>&1");
        }

        &scheduler::schedulerSetNumberOfProcesses(getGlobal("cnsConcurrency"));
        &scheduler::schedulerFinish();
    }
}



sub postUnitiggerConsensus (@) {
    my @cgbFiles  = @_;

    system("mkdir $wrk/5-consensus") if (! -d "$wrk/5-consensus");

    goto alldone if (-e "$wrk/5-consensus/consensus.success");

    createPostUnitiggerConsensusJobs(@cgbFiles);

    #
    #  Check that consensus finished properly
    #

    my $failedJobs = 0;

    foreach my $f (@cgbFiles) {
        if ($f =~ m/^.*(\d\d\d).cgb$/) {
            if ((! -e "$wrk/5-consensus/${asm}_$1.success") ||
                (! -e "$wrk/5-consensus/${asm}_$1.cgi")) {
                print STDERR "$wrk/5-consensus/${asm}_$1 failed -- no .success or no .cgi!\n";
                $failedJobs++;
            }
        } else {
            caFailure("unitigger file '$f' didn't match ###.cgb", undef);
        }
    }

    #  FAILUREHELPME

    caFailure("$failedJobs consensusAfterUnitigger jobs failed", undef) if ($failedJobs);

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

    touch("$wrk/5-consensus/consensus.success");

  alldone:
    stopAfter("consensusAfterUnitigger");
}

1;
