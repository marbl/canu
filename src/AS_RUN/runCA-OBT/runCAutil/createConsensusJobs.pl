use strict;

#  Prepare for consensus on the grid
#    Partition the contigs
#    Repartition the frag store


sub createConsensusJobs {

    system("mkdir $wrk/8-consensus") if (! -d "$wrk/8-consensus");

    my $lastckpt = findLastCheckpoint("7-CGW");


    ########################################

    #  I guess this is the parition size, for human, it was 400000
    my $partitionSize = 400000;

    if (! -e "$wrk/8-consensus/partitionSDB1.success") {
        my $cmd;
        $cmd  = "cd $wrk/8-consensus && ";
        $cmd .= "$bin/PartitionSDB1 $wrk/7-CGW/$asm.SeqStore $lastckpt $partitionSize $wrk/7-CGW/$asm.cgw_contigs ";
        $cmd .= "> $wrk/8-consensus/partitionSDB1.out ";
        $cmd .= "2> $wrk/8-consensus/partitionSDB1.err";

        if (runCommand($cmd)) {
            print STDERR "Failed.\n";
            exit(1);
        }

        touch("$wrk/8-consensus/partitionSDB1.success");
    }

    if (! -e "$wrk/8-consensus/partitionSDB2.success") {
        my $cmd;
        $cmd  = "cd $wrk/8-consensus && ";
        $cmd .= "$bin/PartitionSDB2 $wrk/7-CGW/$asm.SeqStore $lastckpt $wrk/8-consensus/UnitigPartition.txt ";
        $cmd .= "> $wrk/8-consensus/partitionSDB2.out ";
        $cmd .= "2> $wrk/8-consensus/partitionSDB2.err";

        if (runCommand($cmd)) {
            print STDERR "Failed.\n";
            exit(1);
        }

        touch("$wrk/8-consensus/partitionSDB2.success");
    }

    
    if (! -e "$wrk/8-consensus/partitionFragStore.success") {
        my $cmd;
        $cmd  = "cd $wrk/8-consensus && ";
        $cmd .= "$bin/partitionFragStore $wrk/8-consensus/FragPartition.txt $wrk/$asm.frgStore $wrk/$asm.frgStore_cns2part ";
        $cmd .= "> $wrk/8-consensus/partitionFragStore.out ";
        $cmd .= "2> $wrk/8-consensus/partitionFragStore.err";

        if (runCommand($cmd)) {
            print STDERR "Failed.\n";
            exit(1);
        }

        touch("$wrk/8-consensus/partitionFragStore.success");
    }

    ########################################

    #  Build consensus jobs for the grid
    #
    if (! -e "$wrk/8-consensus/jobsCreated.success") {

        if ($useGrid) {
            open(SUB, "> $wrk/8-consensus/submit.sh") or die;
            print SUB "#!/bin/sh\n";
        }

        open(CGW, "ls $wrk/7-CGW/$asm.cgw_contigs.* |") or die;
        while (<CGW>) {
            chomp;

            if (m/cgw_contigs.(\d+)/) {
                my $jobName   = substr("0000000000" . $1, -8);

                open(F, "> $wrk/8-consensus/$jobName.sh") or die;
                print F "#!/bin/sh\n";
                print F "$processStats \\\n";
                print F "$gin/consensus \\\n";
                print F "  -P \\\n";
                print F "  -s $wrk/7-CGW/$asm.SeqStore \\\n";
                print F "  -V $lastckpt \\\n";
                print F "  -p $1 \\\n";
                print F "  -S $1 \\\n";
                print F "  -m \\\n";
                print F "  -o $wrk/8-consensus/$asm.cns_contigs.$1 \\\n";
                print F "  $wrk/$asm.frgStore_cns2part \\\n";
                print F "  $wrk/7-CGW/$asm.cgw_contigs.$1 \\\n";
                print F "&& \\\n";
                print F "touch $wrk/8-consensus/$jobName.success\n";
                close(F);

                chmod 0755, "$wrk/8-consensus/$jobName.sh";

                if ($useGrid) {
                    print SUB "qsub ";
                    print SUB "-p 0 ";  #  Priority
                    print SUB "-r y ";  #  Rerunnable
                    print SUB "-N cns2_${asm}_$jobName ";
                    print SUB "-o $wrk/8-consensus/$jobName.out ";
                    print SUB "-e $wrk/8-consensus/$jobName.err ";
                    print SUB "$wrk/8-consensus/$jobName.sh\n";
                } else {
                    if (runCommand("sh $wrk/8-consensus/$jobName.sh")) {
                        print STDERR "Failed.\n";
                        exit(1);
                    }
                }
            } else {
                print STDERR "Didn't match cgw_contigs.# in $_\n";
            }
        }
        close(CGW);

        if ($useGrid) {
            close(SUB);
        }

        touch("$wrk/8-consensus/jobsCreated.success");

        if ($useGrid) {
            pleaseExecute("$wrk/8-consensus/submit.sh");
            exit(0);
        }
    }
}

1;
