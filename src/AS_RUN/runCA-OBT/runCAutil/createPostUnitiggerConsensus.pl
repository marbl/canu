 use strict;

#  Create the post-unitigger consensus jobs.

sub createPostUnitiggerConsensusJobs {

    system("mkdir $wrk/5-consensus") if (! -d "$wrk/5-consensus");

    #  This is $AS_ROOT/scripts/make_partitionFile, but much safer.
    #  And just a little more obvious, sigh.

    if (! -e "$wrk/5-consensus/$asm.partFile") {
        print STDERR "Starting c -- partition fragStore\n";

        open(F, "ls $wrk/4-unitigger/*.cgb |") or die;
        while (<F>) {
            chomp;

            if (m/^.*(\d\d\d).cgb$/) {
                if (runCommand("grep mid: $_ | sed 's/mid:/$1 /' >> $wrk/5-consensus/$asm.partFile")) {
                    print STDERR "Failed to grep mid: from CGB output in $_.\n";
                    exit(1);
                }
            } else {
                print STDERR "WARNING: didn't match $_ for CGB filename!\n";
            }
        }
        close(F);

        my $cmd;
        $cmd  = "$bin/partitionFragStore ";
        $cmd .= "$wrk/5-consensus/$asm.partFile ";
        $cmd .= "$wrk/$asm.frgStore ";
        $cmd .= "$wrk/$asm.frgStore_cns1part ";
        $cmd .= "> $wrk/5-consensus/partitionfragstore.out ";
        $cmd .= "2> $wrk/5-consensus/partitionfragstore.err";

        if (runCommand($cmd)) {
            print STDERR "Failed to partition the fragStore.\n";
            rename "$wrk/5-consensus/$asm.partFile", "$wrk/5-consensus/$asm.partFile.FAILED";
            exit(1);
        }
    }

    #  Nasty bits of awk here.  It creates a set of shell scripts to
    #  run consensus, one per cgb batch.  The last batch is not used,
    #  the small tests BPW has tried always as an empty file there.

    #  Consensus also uses the prefix of the input cgb to name the
    #  output cgi, so we link those in.

    if (! -e "$wrk/5-consensus/jobsCreated.success") {
        print STDERR "Starting d -- create post-unitigger consensus\n";

        if (runCommand("ln -f $wrk/4-unitigger/*.cgb $wrk/5-consensus")) {
            print STDERR "Linking CGB input into consensus directory failed.\n";
            exit(1);
        }

        if ($useGrid) {
            open(SUB, "> $wrk/5-consensus/submit.sh") or die;
            print SUB "#!/bin/sh\n";
        }

        open(CGB, "ls $wrk/5-consensus/*.cgb |") or die;
        while (<CGB>) {
            chomp;

            if (m/^.*(\d\d\d).cgb$/) {
                my $jobName   = $1;

                open(F, "> $wrk/5-consensus/$jobName.sh") or die;
                print F "#!/bin/sh\n";
                print F "$processStats \\\n";
                print F "$gin/consensus -P -S $1 -m -U -z $wrk/$asm.frgStore_cns1part $wrk/5-consensus/${asm}_$1.cgb \\\n";
                print F "&& \\\n";
                print F "touch $wrk/5-consensus/$jobName.success\n";
                close(F);

                chmod 0755, "$wrk/5-consensus/$jobName.sh";

                if ($useGrid) {
                    print SUB "qsub ";
                    print SUB "-p 0 ";  #  Priority
                    print SUB "-r y ";  #  Rerunnable
                    print SUB "-N cns1_${asm}_$1 ";
                    print SUB "-o $wrk/5-consensus/$jobName.out ";
                    print SUB "-e $wrk/5-consensus/$jobName.err ";
                    print SUB "$wrk/5-consensus/$jobName.sh\n";
                } else {
                    if (runCommand("sh $wrk/5-consensus/$jobName.sh")) {
                        print STDERR "Failed.\n";
                        exit(1);
                    }
                }

            } else {
                print STDERR "WARNING: didn't match $_ for CGB filename!\n";
            }
        }
        close(CGB);

        if ($useGrid) {
            close(SUB);
        }

        touch("$wrk/5-consensus/jobsCreated.success");

        if ($useGrid) {
            pleaseExecute("$wrk/5-consensus/submit.sh");
            exit(0);
        }
    }
}

1;
