use strict;



sub overlapCorrection {
    my $cleanup = 1;

    return if (getGlobal("doFragmentCorrection") == 0);

    return if (-e "$wrk/3-overlapcorrection/$asm.erates.updated");
    return if (-e "$wrk/$asm.ovlStore/corrected");

    system("mkdir $wrk/3-overlapcorrection") if (! -e "$wrk/3-overlapcorrection");

    if ((getGlobal("ovlOverlapper") eq "ovl") && (! -e "$wrk/3-overlapcorrection/frgcorr.sh")) {
        my $batchSize   = getGlobal("frgCorrBatchSize");
        my $numThreads  = getGlobal("frgCorrThreads");
        my $jobs        = int($numFrags / $batchSize) + (($numFrags % $batchSize == 0) ? 0 : 1);

        open(F, "> $wrk/3-overlapcorrection/frgcorr.sh") or caFailure("failed to write to '$wrk/3-overlapcorrection/frgcorr.sh'", undef);
        print F "#!" . getGlobal("shell") . "\n\n";
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
        print F "minid=`expr \$jobid \\* $batchSize - $batchSize + 1`\n";
        print F "maxid=`expr \$jobid \\* $batchSize`\n";
        print F "runid=\$\$\n";
        print F "\n";
        print F "if [ \$maxid -gt $numFrags ] ; then\n";
        print F "  maxid=$numFrags\n";
        print F "fi\n";
        print F "if [ \$minid -gt \$maxid ] ; then\n";
        print F "  echo Job partitioning error -- minid=\$minid maxid=\$maxid.\n";
        print F "  exit\n";
        print F "fi\n";
        print F "\n";
        print F "AS_OVL_ERROR_RATE=", getGlobal("ovlErrorRate"), "\n";
        print F "AS_CNS_ERROR_RATE=", getGlobal("cnsErrorRate"), "\n";
        print F "AS_CGW_ERROR_RATE=", getGlobal("cgwErrorRate"), "\n";
        print F "export AS_OVL_ERROR_RATE AS_CNS_ERROR_RATE AS_CGW_ERROR_RATE\n";
        print F "\n";
        print F "if [ -e $wrk/3-overlapcorrection/\$jobid.frgcorr ] ; then\n";
        print F "  echo Job previously completed successfully.\n";
        print F "  exit\n";
        print F "fi\n";

        print F getBinDirectoryShellCode();

        print F "\$bin/correct-frags \\\n";
        print F "  -t $numThreads \\\n";
        print F "  -S $wrk/$asm.ovlStore \\\n";
        print F "  -o $wrk/3-overlapcorrection/\$jobid.frgcorr.WORKING \\\n";
        print F "  $wrk/$asm.gkpStore \\\n";
        print F "  \$minid \$maxid \\\n";
        print F "&& \\\n";
        print F "mv $wrk/3-overlapcorrection/\$jobid.frgcorr.WORKING $wrk/3-overlapcorrection/\$jobid.frgcorr\n";

        close(F);

        chmod 0755, "$wrk/3-overlapcorrection/frgcorr.sh";

        if (getGlobal("frgCorrOnGrid") && getGlobal("useGrid")) {
            #  Run the correction job on the grid.

            my $sge                   = getGlobal("sge");
            my $sgeName               = getGlobal("sgeName");
            my $sgeFragmentCorrection = getGlobal("sgeFragmentCorrection");

            $sgeName = "_$sgeName" if (defined($sgeName));

            my $SGE;
            $SGE  = "qsub $sge $sgeFragmentCorrection -cwd -N frg_$asm$sgeName ";
            $SGE .= "-t 1-$jobs ";
            $SGE .= " -j y -o $wrk/3-overlapcorrection/\\\$TASK_ID.err ";
            $SGE .= "$wrk/3-overlapcorrection/frgcorr.sh\n";

            submitBatchJobs($SGE, "frg_$asm$sgeName");
            exit(0);
        } else {
            #  Run the correction job right here, right now.

            for (my $i=1; $i<=$jobs; $i++) {
                my $out = substr("0000" . $i, -4);
                &scheduler::schedulerSubmit("$wrk/3-overlapcorrection/frgcorr.sh $i > $wrk/3-overlapcorrection/$out.err 2>&1");
            }

            &scheduler::schedulerSetNumberOfProcesses($global{"frgCorrConcurrency"});
            &scheduler::schedulerFinish();
        }
    }

    #
    #  MERGE CORRECTION
    #

    if (! -e "$wrk/3-overlapcorrection/$asm.frgcorr") {
        my $batchSize  = (getGlobal("ovlOverlapper") eq "mer") ? getGlobal("merOverlapperExtendBatchSize") : getGlobal("frgCorrBatchSize");
        my $jobs       = int($numFrags / $batchSize) + (($numFrags % $batchSize == 0) ? 0 : 1);
        my $failedJobs = 0;

        open(F, "> $wrk/3-overlapcorrection/cat-corrects.frgcorrlist");
        for (my $i=1; $i<=$jobs; $i++) {
            my $jobid = substr("0000" . $i, -4);

            if (! -e "$wrk/3-overlapcorrection/$jobid.frgcorr") {
                print STDERR "Fragment correction job $jobid failed.\n";
                $failedJobs++;
            }

            print F "$wrk/3-overlapcorrection/$jobid.frgcorr\n";
        }
        close(F);

        #  FAILUREHELPME

        if ($failedJobs) {
            if (getGlobal("ovlOverlapper") eq "ovl") {
                caFailure("$failedJobs overlap jobs failed; remove $wrk/3-overlapcorrection/frgcorr.sh to try again", undef);
            } else {
                caFailure("$failedJobs overlap jobs failed due to mer overlap seed extension", undef);
            }
        }

        my $bin = getBinDirectory();
        my $cmd;
        $cmd  = "$bin/cat-corrects ";
        $cmd .= "-L $wrk/3-overlapcorrection/cat-corrects.frgcorrlist ";
        $cmd .= "-o $wrk/3-overlapcorrection/$asm.frgcorr ";
        $cmd .= "> $wrk/3-overlapcorrection/cat-corrects.err 2>&1";

        if (runCommand("$wrk/3-overlapcorrection", $cmd)) {
            rename "$wrk/3-overlapcorrection/$asm.frgcorr", "$wrk/3-overlapcorrection/$asm.frgcorr.FAILED";
            caFailure("failed to concatenate the fragment corrections", "$wrk/3-overlapcorrection/cat-corrects.err");
        }

        if ($cleanup) {
            open(F, "< $wrk/3-overlapcorrection/cat-corrects.frgcorrlist");
            while (<F>) {
                if (m/^(.*)\/([0-9]*).frgcorr/) {
                    #unlink "$1/$2.frgcorr";
                    #unlink "$1/$2.err";
                    my $sge = int($2);
                    #unlink "$1/$sge.err";
                }
            }
            close(F);
            #unlink "$wrk/3-overlapcorrection/cat-corrects.frgcorrlist";
            #unlink "$wrk/3-overlapcorrection/cat-corrects.err";
        }
    }

    #
    #  CREATE OVERLAP CORRECTION
    #

    if (! -e "$wrk/3-overlapcorrection/ovlcorr.sh") {
        my $batchSize  = getGlobal("ovlCorrBatchSize");
        my $jobs       = int($numFrags / $batchSize) + (($numFrags % $batchSize == 0) ? 0 : 1);

        open(F, "> $wrk/3-overlapcorrection/ovlcorr.sh") or caFailure("failed to write '$wrk/3-overlapcorrection/ovlcorr.sh'", undef);
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
        print F "jobid=`printf %04d \$jobid`\n";
        print F "frgBeg=`expr \$jobid \\* $batchSize - $batchSize + 1`\n";
        print F "frgEnd=`expr \$jobid \\* $batchSize`\n";
        print F "if [ \$frgEnd -ge $numFrags ] ; then\n";
        print F "  frgEnd=$numFrags\n";
        print F "fi\n";
        print F "frgBeg=`printf %08d \$frgBeg`\n";
        print F "frgEnd=`printf %08d \$frgEnd`\n";

        print F getBinDirectoryShellCode();

        print F "if [ ! -e $wrk/3-overlapcorrection/\$jobid.erate ] ; then\n";
        print F "  \$bin/correct-olaps \\\n";
        print F "    -S $wrk/$asm.ovlStore \\\n";
        print F "    -e $wrk/3-overlapcorrection/\$jobid.erate.WORKING \\\n";
        print F "    $wrk/$asm.gkpStore \\\n";
        print F "    $wrk/3-overlapcorrection/$asm.frgcorr \\\n";
        print F "    \$frgBeg \$frgEnd \\\n";
        print F "  &&  \\\n";
        print F "  mv $wrk/3-overlapcorrection/\$jobid.erate.WORKING $wrk/3-overlapcorrection/\$jobid.erate\n";
        print F "fi\n";
        close(F);

        chmod 0755, "$wrk/3-overlapcorrection/ovlcorr.sh";

        if (getGlobal("ovlCorrOnGrid") && getGlobal("useGrid")) {
            #  Run the correction job on the grid.

            my $sge                   = getGlobal("sge");
            my $sgeName               = getGlobal("sgeName");
            my $sgeOverlapCorrection  = getGlobal("sgeOverlapCorrection");

            $sgeName = "_$sgeName" if (defined($sgeName));

            my $SGE;
            $SGE  = "qsub $sge $sgeOverlapCorrection -cwd -N ovc_$asm$sgeName ";
            $SGE .= "-t 1-$jobs ";
            $SGE .= " -j y -o $wrk/3-overlapcorrection/\\\$TASK_ID.err ";
            $SGE .= "$wrk/3-overlapcorrection/ovlcorr.sh\n";

            submitBatchJobs($SGE, "ovc_$asm$sgeName");
            exit(0);
        } else {
            #  Run the correction job right here, right now.

            for (my $i=1; $i<=$jobs; $i++) {
                my $out = substr("0000" . $i, -4);
                &scheduler::schedulerSubmit("$wrk/3-overlapcorrection/ovlcorr.sh $i > $wrk/3-overlapcorrection/$out.err 2>&1");
            }

            &scheduler::schedulerSetNumberOfProcesses($global{"ovlCorrConcurrency"});
            &scheduler::schedulerFinish();
        }
    }

    #
    #  APPLY OVERLAP CORRECTION
    #

    if (! -e "$wrk/3-overlapcorrection/$asm.erates.updated") {
        my $batchSize   = getGlobal("ovlCorrBatchSize");
        my $bin         = getBinDirectory();
        my $failedJobs  = 0;
        my $jobs        = int($numFrags / $batchSize) + (($numFrags % $batchSize == 0) ? 0 : 1);
        my $cmd;

        open(F, "> $wrk/3-overlapcorrection/cat-erates.eratelist");
        for (my $i=1; $i<=$jobs; $i++) {
            my $jobid = substr("0000" . $i, -4);

            if (! -e "$wrk/3-overlapcorrection/$jobid.erate") {
                print STDERR "Overlap correction job $i ($wrk/3-overlapcorrection/$jobid) failed.\n";
                $failedJobs++;
            }

            print F "$wrk/3-overlapcorrection/$jobid.erate\n";
        }
        close(F);

        #  FAILUREHELPME

        if ($failedJobs) {
            caFailure("$failedJobs overlap correction jobs failed; remove $wrk/3-overlapcorrection/ovlcorr.sh (or run by hand) to try again", undef);
        }

        #unlink "$wrk/3-overlapcorrection/$asm.frgcorr" if ($cleanup);

        $cmd  = "$bin/cat-erates ";
        $cmd .= "-L $wrk/3-overlapcorrection/cat-erates.eratelist ";
        $cmd .= "-o $wrk/3-overlapcorrection/$asm.erates ";
        $cmd .= "> $wrk/3-overlapcorrection/cat-erates.err 2>&1";
        if (runCommand("$wrk/3-overlapcorrection", $cmd)) {
            rename "$wrk/3-overlapcorrection/$asm.erates", "$wrk/3-overlapcorrection/$asm.erates.FAILED";
            caFailure("failed to concatenate the overlap erate corrections", "$wrk/3-overlapcorrection/cat-erates.err");
        }

        $cmd  = "$bin/overlapStore ";
        $cmd .= " -u $wrk/$asm.ovlStore ";
        $cmd .= " $wrk/3-overlapcorrection/$asm.erates";
        $cmd .= "> $wrk/3-overlapcorrection/overlapStore-update-erates.err 2>&1";
        if (runCommand("$wrk/3-overlapcorrection", $cmd)) {
            caFailure("failed to apply the overlap corrections", "$wrk/3-overlapcorrection/overlapStore-update-erates.err");
        }

        touch("$wrk/3-overlapcorrection/$asm.erates.updated");
        touch("$wrk/$asm.ovlStore/corrected");

        if ($cleanup) {
            open(F, "< $wrk/3-overlapcorrection/cat-erates.eratelist");
            while (<F>) {
                if (m/^(.*)\/([0-9]*).erate/) {
                    #unlink "$1/$2.erate";
                    #unlink "$1/$2.err";
                    my $sge = int($2);
                    #unlink "$1/$sge.err";
                }
            }
            close(F);

            #unlink "$wrk/3-overlapcorrection/overlapStore-update-erates.err";
            #unlink "$wrk/3-overlapcorrection/$asm.erates";

            #unlink "$wrk/3-overlapcorrection/cat-erates.err";
            #unlink "$wrk/3-overlapcorrection/cat-erates.eratelist";

            #unlink "$wrk/3-overlapcorrection/frgcorr.sh";
            #unlink "$wrk/3-overlapcorrection/ovlcorr.sh";
        }
    }
}

1;
