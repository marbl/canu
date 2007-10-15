use strict;

sub merOverlapper($) {
    my $isTrim = shift @_;

    return if (-d "$wrk/$asm.ovlStore");

    caFailure("createOverlapJobs()-- Help!  I have no frags!\n") if ($numFrags == 0);
    caFailure("createOverlapJobs()-- I need to know if I'm trimming or assembling!\n") if (!defined($isTrim));

    my $ovlThreads        = getGlobal("ovlThreads");
    my $ovlMemory         = getGlobal("ovlMemory");
    my $scratch           = getGlobal("scratch");

    my $outDir  = "1-overlapper";
    my $ovlOpt  = "";
    my $merSize = getGlobal("merSizeOvl");
    my $merComp = getGlobal("merCompression");

    if ($isTrim eq "trim") {
        $outDir  = "0-overlaptrim-overlap";
        $ovlOpt  = "-G";
        $merSize = getGlobal("merSizeObt");
    }

    system("mkdir $wrk/$outDir") if (! -d "$wrk/$outDir");

    return if (-e "$wrk/$outDir/jobsCreated.success");

    my $cmd;

    #
    #  Currently, overmerry runs as one gigantic process, threaded.  A
    #  microbe M063 with 21,150 sanger reads and a full run of flx,
    #  450,403 reads, takes 2600 seconds using 8 threads.  It used
    #  1.6GB memory.
    #
    #  To make it run in multiple segments, we need to have meryl counts
    #  around, AND overmerry needs to know how to use those.
    #
    if (! -e "$wrk/$outDir/$asm.ovm") {
        $cmd  = "$bin/overmerry";
        $cmd .= " -t " . getGlobal("merOverlapThreads");
        $cmd .= " -g $wrk/$asm.gkpStore";
        $cmd .= " -m $merSize";
        $cmd .= " -c $merComp";
        $cmd .= " -o $wrk/$outDir/$asm.ovm";
        $cmd .= " > $wrk/$outDir/$asm.ovm.err 2>&1";
        if (runCommand("$wrk/$outDir", $cmd)) {
            rename "$wrk/$outDir/$asm.ovm", "$wrk/$outDir/$asm.ovm.FAILED";
            caFailure("Failed.\n");
        }
    }

    #  Same microbe, 26 seconds.
    #
    if (! -e "$wrk/$outDir/$asm.merStore") {
        $cmd  = "$bin/overlapStore";
        $cmd .= " -c $wrk/$outDir/$asm.merStore";
        $cmd .= " -M " . getGlobal("ovlStoreMemory");
        $cmd .= " -m $numFrags";
        $cmd .= " $wrk/$outDir/$asm.ovm";
        $cmd .= " > $wrk/$outDir/$asm.merStore.err 2>&1";
        if (runCommand("$wrk/$outDir", $cmd)) {
            rename "$wrk/$outDir/$asm.merStore", "$wrk/$outDir/$asm.merStore.FAILED";
            caFailure("Failed.\n");
        }
    }

    #
    #  2378 seconds.
    #

    if (! -e "$wrk/$outDir/jobsCreated.success") {
        my $batchSize = getGlobal("ovlCorrBatchSize");
        my $jobs      = int($numFrags / ($batchSize-1)) + 1;

        open(F, "> $wrk/$outDir/olap-from-seeds.sh") or caFailure("Can't open '$wrk/$outDir/olap-from-seeds.sh'\n");

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

        print F "if [ ! -d $wrk/$outDir/olaps ]; then\n";
        print F "  mkdir $wrk/$outDir/olaps\n";
        print F "fi\n";
        print F "\n";

        print F "if [ -e $wrk/$outDir/olaps/\$jobid.success ]; then\n";
        print F "  echo Job previously completed successfully.\n";
        print F "  exit\n";
        print F "fi\n";
        print F "\n";

        print F "echo \\\n";
        print F "$bin/olap-from-seeds \\\n";
        print F " -a -b \\\n";
        print F " -t 1 \\\n";
        print F " -S $wrk/$outDir/$asm.merStore \\\n";
        print F " -c $wrk/2-frgcorr/\$jobid.frgcorr \\\n"       if ($isTrim ne "trim");
        print F " -o $wrk/$outDir/olaps/$asm.\$jobid.ovb \\\n"  if ($isTrim ne "trim");
        print F " -o $wrk/$outDir/olaps/$asm.\$jobid.ovr \\\n"  if ($isTrim eq "trim");
        print F " $wrk/$asm.gkpStore \\\n";
        print F " \$minid \$maxid \\\n";
        print F " \\> $wrk/$outDir/olaps/$asm.\$jobid.ovb.err 2\\>\\&1\n";
        print F "\n";

        print F "$bin/olap-from-seeds \\\n";
        print F " -a -b \\\n";
        print F " -t 1 \\\n";
        print F " -S $wrk/$outDir/$asm.merStore \\\n";
        print F " -c $wrk/2-frgcorr/\$jobid.frgcorr \\\n"       if ($isTrim ne "trim");
        print F " -o $wrk/$outDir/olaps/$asm.\$jobid.ovb \\\n"  if ($isTrim ne "trim");
        print F " -o $wrk/$outDir/olaps/$asm.\$jobid.ovr \\\n"  if ($isTrim eq "trim");
        print F " $wrk/$asm.gkpStore \\\n";
        print F " \$minid \$maxid \\\n";
        print F " > $wrk/$outDir/olaps/$asm.\$jobid.ovb.err 2>&1 \\\n";

        if ($isTrim eq "trim") {
            print F " &&  \\\n";
            print F "$gin/acceptableOBToverlap \\\n";
            print F " < $wrk/$outDir/$asm.ovr \\\n";
            print F " > $wrk/$outDir/$asm.ovb \\\n";
        } else {
            print F "&& \\\n";
            print F "touch $wrk/2-frgcorr/\$jobid.success \\\n";
        }

        print F "&& \\\n";
        print F "touch $wrk/$outDir/olaps/\$jobid.success\n";

        close(F);


        #  Make the 2-frgcorr directory (to hold the corrections
        #  output) and claim that fragment correction is all done.
        #  after this, the rest of the fragment/overlap correction
        #  pipeline Just Works.
        #
        if ($isTrim ne "trim") {
            system("mkdir $wrk/2-frgcorr") if (! -d "$wrk/2-frgcorr");
            touch("$wrk/2-frgcorr/jobsCreated.success");
        }


        #  Submit to the grid (or tell the user to do it), or just run
        #  things here
        #
        if (getGlobal("useGrid") && getGlobal("ovlOnGrid")) {
            my $sge        = getGlobal("sge");
            my $sgeOverlap = getGlobal("sgeOverlap");

            my $SGE;
            $SGE  = "qsub $sge $sgeOverlap -r y -N NAME \\\n";
            $SGE .= "  -t MINMAX \\\n";
            $SGE .= "  -j y -o $wrk/$outDir/olap-from-seeds.\\\$TASK_ID.out \\\n";
            $SGE .= "  $wrk/$outDir/olap-from-seeds.sh\n";

            my $waitTag = submitBatchJobs("ovl", $SGE, $jobs, $ovlThreads);

            if (runningOnGrid()) {
                touch("$wrk/$outDir/jobsCreated.success");
                submitScript("$waitTag");
                exit(0);
            } else {
                touch("$wrk/$outDir/jobsCreated.success");
                exit(0);
            }
        } else {
            my $failures = 0;
            for (my $i=1; $i<=$jobs; $i++) {
                my $out = substr("0000" . $i, -4);
                if (runCommand("$wrk/$outDir", "$wrk/$outDir/overlap.sh $i > $wrk/$outDir/olap-from-seeds.$out.out 2>&1")) {
                    $failures++;
                }
            }
            if ($failures == 0) {
                touch("$wrk/$outDir/jobsCreated.success");
            }
        }
    }
}
