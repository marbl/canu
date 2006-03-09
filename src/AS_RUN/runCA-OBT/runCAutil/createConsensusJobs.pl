use strict;

#  Prepare for consensus on the grid
#    Partition the contigs
#    Repartition the frag store

sub createConsensusJobs {
    my $pstats            = getGlobal("processStats");

    return if (-e "$wrk/8-consensus/jobsCreated.success");

    my $partitionSize = int($numFrags / getGlobal("cnsPartitions"));
    $partitionSize = getGlobal("cnsMinFrags") if ($partitionSize < getGlobal("cnsMinFrags"));

    system("mkdir $wrk/8-consensus") if (! -d "$wrk/8-consensus");

    my $lastckpt = findLastCheckpoint("7-CGW");

    if (! -e "$wrk/8-consensus/partitionSDB1.success") {
        my $cmd;
        $cmd  = "cd $wrk/8-consensus && ";
        $cmd .= "$bin/PartitionSDB1 $wrk/7-CGW/$asm.SeqStore $lastckpt $partitionSize $wrk/7-CGW/$asm.cgw_contigs ";
        $cmd .= "> $wrk/8-consensus/partitionSDB1.err 2>&1";

        die "Failed.\n" if (runCommand($cmd));
        touch("$wrk/8-consensus/partitionSDB1.success");
    }

    if (! -e "$wrk/8-consensus/partitionSDB2.success") {
        my $cmd;
        $cmd  = "cd $wrk/8-consensus && ";
        $cmd .= "$bin/PartitionSDB2 $wrk/7-CGW/$asm.SeqStore $lastckpt $wrk/8-consensus/UnitigPartition.txt ";
        $cmd .= "> $wrk/8-consensus/partitionSDB2.err 2>&1";

        die "Failed.\n" if (runCommand($cmd));
        touch("$wrk/8-consensus/partitionSDB2.success");
    }
    
    if (! -e "$wrk/8-consensus/partitionFragStore.success") {
        my $cmd;
        $cmd  = "cd $wrk/8-consensus && ";
        $cmd .= "$bin/partitionFragStore $wrk/8-consensus/FragPartition.txt $wrk/$asm.frgStore $wrk/$asm.frgStore_cns2part ";
        $cmd .= "> $wrk/8-consensus/partitionFragStore.out 2>&1";

        die "Failed.\n" if (runCommand($cmd));
        touch("$wrk/8-consensus/partitionFragStore.success");
    }

    ########################################
    #
    #  Build consensus jobs for the grid -- this is very similar to that in createPostUnitiggerConsensus.pl
    #
    my $jobP;
    my $jobs = 0;

    open(CGW, "ls $wrk/7-CGW/$asm.cgw_contigs.* |") or die;
    while (<CGW>) {
        if (m/cgw_contigs.(\d+)/) {
            $jobP .= "$1\t";
            $jobs++;
        } else {
            print STDERR "Didn't match cgw_contigs.# in $_\n";
        }
    }
    close(CGW);

    $jobP = join ' ', sort { $a <=> $b } split '\s+', $jobP;

    open(F, "> $wrk/8-consensus/consensus.sh") or die "Can't open '$wrk/8-consensus/consensus.sh'\n";
    print F "#!/bin/sh\n";
    print F "\n";
    print F "jobid=\$SGE_TASK_ID\n";
    print F "if [ x\$jobid = x ]; then\n";
    print F "  jobid=\$1\n";
    print F "fi\n";
    print F "if [ x\$jobid = x ]; then\n";
    print F "  echo Error: I need SGE_TASK_ID set, or a job index on the command line.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "jobp=`echo $jobP | cut -d' ' -f \$jobid`\n";
    print F "\n";
    print F "if [ -e $wrk/8-consensus/$asm.cns_contigs.\$jobp.success ] ; then\n";
    print F "  exit 0\n";
    print F "fi\n";
    print F "\n";
    print F "echo \\\n";
    print F "$gin/consensus \\\n";
    print F "  -P \\\n";
    print F "  -s $wrk/7-CGW/$asm.SeqStore \\\n";
    print F "  -V $lastckpt \\\n";
    print F "  -p \$jobp \\\n";
    print F "  -S \$jobp \\\n";
    print F "  -m \\\n";
    print F "  -o $wrk/8-consensus/$asm.cns_contigs.\$jobp \\\n";
    print F "  $wrk/$asm.frgStore_cns2part \\\n";
    print F "  $wrk/7-CGW/$asm.cgw_contigs.\$jobp \\\n";
    print F " \\> $wrk/8-consensus/$asm.cns_contigs.\$jobp.err 2\\>\\&1\n";
    print F "\n";
    print F "$pstats \\\n" if (defined($pstats));
    print F "$gin/consensus \\\n";
    print F "  -P \\\n";
    print F "  -s $wrk/7-CGW/$asm.SeqStore \\\n";
    print F "  -V $lastckpt \\\n";
    print F "  -p \$jobp \\\n";
    print F "  -S \$jobp \\\n";
    print F "  -m \\\n";
    print F "  -o $wrk/8-consensus/$asm.cns_contigs.\$jobp \\\n";
    print F "  $wrk/$asm.frgStore_cns2part \\\n";
    print F "  $wrk/7-CGW/$asm.cgw_contigs.\$jobp \\\n";
    print F " > $wrk/8-consensus/$asm.cns_contigs.\$jobp.err 2>&1 \\\n";
    print F "&& \\\n";
    print F "touch $wrk/8-consensus/$asm.cns_contigs.\$jobp.success\n";
    print F "exit 0\n";
    close(F);

    chmod 0755, "$wrk/8-consensus/consensus.sh";

    if ($useGrid) {
        my $cmd;
        $cmd  = "qsub -p 0 -r y -N cns2_${asm} ";
        $cmd .= "-t 1-$jobs ";
        $cmd .= "-j y -o /dev/null ";
        $cmd .= "$wrk/8-consensus/consensus.sh\n";
        pleaseExecute($cmd);
        touch("$wrk/8-consensus/jobsCreated.success");
        exit(0);
    } else {
        for (my $i=1; $i<=$jobs; $i++) {
            &scheduler::schedulerSubmit("sh $wrk/8-consensus/consensus.sh $i");

            #$ENV{'SGE_TASK_ID'} = $i;
            #if (runCommand("$wrk/8-consensus/consensus.sh")) {
            #    print STDERR "Failed job $i\n";
            #    exit(1);
            #}
            #delete $ENV{'SGE_TASK_ID'};
        }

        &scheduler::schedulerSetNumberOfProcesses(getGlobal("cnsConcurrency"));
        &scheduler::schedulerFinish();

        touch("$wrk/8-consensus/jobsCreated.success");
    }
}

1;
