package ca3g::OverlapErrorAdjustment;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(readErrorDetectionConfigure readErrorDetectionCheck overlapErrorAdjustmentConfigure overlapErrorAdjustmentCheck updateOverlapStore);

use strict;

use File::Path qw(make_path remove_tree);

use ca3g::Defaults;
use ca3g::Execution;
use ca3g::Gatekeeper;




sub concatOutput ($@) {
    my $outName     = shift @_;
    my @successJobs =       @_;

    open(O, "> $outName");
    binmode(O);

    foreach my $f (@successJobs) {
        open(F, "< $f");
        binmode(F);

        my $buf;
        my $len = sysread(F, $buf, 1024 * 1024);

        while ($len > 0) {
            syswrite(O, $buf, $len);

            $len = sysread(F, $buf, 1024 * 1024);
        }

        close(F);
    }

    close(O);

    foreach my $f (@successJobs) {
        unlink $f;
    }
}




sub readErrorDetectionConfigure ($$) {
    my $wrk   = shift @_;
    my $asm   = shift @_;

    return  if (getGlobal("enableOEA") == 0);
    return  if (skipStage($wrk, $asm, "readErrorDetectionConfigure") == 1);
    return  if (-e "$wrk/3-overlapErrorAdjustment/red.red");

    return if (-e "$wrk/$asm.ovlStore/adjustedErates");
    return if (-d "$wrk/$asm.tigStore");

    make_path("$wrk/3-overlapErrorAdjustment")  if (! -d "$wrk/3-overlapErrorAdjustment");

    my $batchSize   = getGlobal("redBatchSize");
    my $numThreads  = getGlobal("redThreads");

    my $numReads    = getNumberOfReadsInStore($wrk, $asm);
    my $numJobs     = int($numReads / $batchSize) + (($numReads % $batchSize == 0) ? 0 : 1);

    open(F, "> $wrk/3-overlapErrorAdjustment/red.sh") or caWarn("can't open '$wrk/3-overlapErrorAdjustment/red.sh' for writing: $!", undef);

    print F "#!" . getGlobal("shell") . "\n\n";
    print F "\n";
    print F "jobid=" . getGlobal("gridEngineTaskID") . "\n";
    print F "if [ x\$jobid = x -o x\$jobid = xundefined -o x\$jobid = x0 ]; then\n";
    print F "  jobid=\$1\n";
    print F "fi\n";
    print F "if [ x\$jobid = x ]; then\n";
    print F "  echo Error: I need " . getGlobal("gridEngineTaskID") . " set, or a job index on the command line.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "\n";
    print F "jobid=`printf %04d \$jobid`\n";
    print F "minid=`expr \$jobid \\* $batchSize - $batchSize + 1`\n";
    print F "maxid=`expr \$jobid \\* $batchSize`\n";
    print F "runid=\$\$\n";
    print F "\n";
    print F "if [ \$maxid -gt $numReads ] ; then\n";
    print F "  maxid=$numReads\n";
    print F "fi\n";
     print F "if [ \$minid -gt \$maxid ] ; then\n";
    print F "  echo Job partitioning error -- minid=\$minid maxid=\$maxid.\n";
    print F "  exit\n";
    print F "fi\n";
    print F "\n";
    print F "if [ -e $wrk/3-overlapErrorAdjustment/\$jobid.red ] ; then\n";
    print F "  echo Job previously completed successfully.\n";
    print F "  exit\n";
    print F "fi\n";

    print F getBinDirectoryShellCode();

    print F "if [ ! -e $wrk/3-overlapErrorAdjustment/\$jobid.red ] ; then\n";
    print F "  \$bin/findErrors \\\n";
    print F "    -G $wrk/$asm.gkpStore \\\n";
    print F "    -O $wrk/$asm.ovlStore \\\n";
    print F "    -R \$minid \$maxid \\\n";
    print F "    -e " . getGlobal("ovlErrorRate") . " -l " . getGlobal("minOverlapLength") . " \\\n";
    print F "    -o $wrk/3-overlapErrorAdjustment/\$jobid.red.WORKING \\\n";
    print F "    -t $numThreads \\\n";
    print F "  && \\\n";
    print F "  mv $wrk/3-overlapErrorAdjustment/\$jobid.red.WORKING $wrk/3-overlapErrorAdjustment/\$jobid.red\n";
    print F "fi\n";

    close(F);

    chmod 0755, "$wrk/3-overlapErrorAdjustment/red.sh";

    emitStage($wrk, $asm, "readErrorDetectionConfigure");
}





sub readErrorDetectionCheck ($$$) {
    my $wrk     = shift @_;
    my $asm     = shift @_;
    my $attempt = shift @_;

    return  if (getGlobal("enableOEA") == 0);
    return  if (skipStage($wrk, $asm, "readErrorDetectionCheck", $attempt) == 1);
    return  if (-e "$wrk/3-overlapErrorAdjustment/red.red");

    return if (-e "$wrk/$asm.ovlStore/adjustedErates");
    return if (-d "$wrk/$asm.tigStore");

    my $batchSize   = getGlobal("redBatchSize");
    my $numReads    = getNumberOfReadsInStore($wrk, $asm);
    my $numJobs     = int($numReads / $batchSize) + (($numReads % $batchSize == 0) ? 0 : 1);

    my $currentJobID   = 1;
    my @successJobs;
    my @failedJobs;
    my $failureMessage = "";

    for (my $currentJobID=1; $currentJobID <= $numJobs; $currentJobID++) {
        my $ji = substr("0000" . $currentJobID, -4);
        my $jn = "$wrk/3-overlapErrorAdjustment/$ji.red";

        if (! -e $jn) {
            $failureMessage .= "  job $jn FAILED.\n";
            push @failedJobs, $currentJobID;
        } else {
            push @successJobs, $jn;
        }
    }

    #  No failed jobs?  Success!

    #  I didn't wan't to concat all the corrections, but it is _vastly_ easier to do so, compared to
    #  hacking correctOverlaps to handle multiple corrections files.  Plus, it is now really just a
    #  concat; before, the files needed to be parsed to strip off a header.

    if (scalar(@failedJobs) == 0) {
        concatOutput("$wrk/3-overlapErrorAdjustment/red.red", @successJobs);
        emitStage($wrk, $asm, "readErrorDetectionCheck");
        return;
    }

    #  If not the first attempt, report the jobs that failed, and that we're recomputing.

    if ($attempt > 1) {
        print STDERR "\n";
        print STDERR scalar(@failedJobs), " read error detection jobs failed:\n";
        print STDERR $failureMessage;
        print STDERR "\n";
    }


    #  If too many attempts, give up.

    if ($attempt > 2) {
        caExit("failed to detect errors in reads.  Made " . ($attempt-1) . " attempts, jobs still failed", undef);
    }

    #  Otherwise, run some jobs.

    print STDERR "readErrorDetectionCheck() -- attempt $attempt begins with ", scalar(@successJobs), " finished, and ", scalar(@failedJobs), " to compute.\n";

    emitStage($wrk, $asm, "readErrorDetectionCheck", $attempt);

    submitOrRunParallelJob($wrk, $asm, "red", "$wrk/3-overlapErrorAdjustment", "red", @failedJobs);
}





sub overlapErrorAdjustmentConfigure ($$) {
    my $wrk   = shift @_;
    my $asm   = shift @_;

    return  if (getGlobal("enableOEA") == 0);
    return  if (skipStage($wrk, $asm, "overlapErrorAdjustmentConfigure") == 1);
    return  if (-e "$wrk/3-overlapErrorAdjustment/oea.sh");

    return if (-e "$wrk/$asm.ovlStore/adjustedErates");
    return if (-d "$wrk/$asm.tigStore");

    my $batchSize   = getGlobal("oeaBatchSize");
    my $numReads    = getNumberOfReadsInStore($wrk, $asm);
    my $numJobs     = int($numReads / $batchSize) + (($numReads % $batchSize == 0) ? 0 : 1);

    open(F, "> $wrk/3-overlapErrorAdjustment/oea.sh") or caExit("can't open '$wrk/3-overlapErrorAdjustment/oea.sh' for writing: $!", undef);

    print F "#!" . getGlobal("shell") . "\n\n";
    print F "\n";
    print F "jobid=" . getGlobal("gridEngineTaskID") . "\n";
    print F "if [ x\$jobid = x -o x\$jobid = xundefined -o x\$jobid = x0 ]; then\n";
    print F "  jobid=\$1\n";
    print F "fi\n";
    print F "if [ x\$jobid = x ]; then\n";
    print F "  echo Error: I need " . getGlobal("gridEngineTaskID") . " set, or a job index on the command line.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "\n";
    print F "if [ \$jobid -gt $numJobs ] ; then\n";
    print F "  exit\n";
    print F "fi\n";
    print F "\n";
    print F "jobid=`printf %04d \$jobid`\n";
    print F "minid=`expr \$jobid \\* $batchSize - $batchSize + 1`\n";
    print F "maxid=`expr \$jobid \\* $batchSize`\n";
    print F "if [ \$maxid -ge $numReads ] ; then\n";
    print F "  maxid=$numReads\n";
    print F "fi\n";
    print F "minid=`printf %08d \$minid`\n";
    print F "maxid=`printf %08d \$maxid`\n";

    print F getBinDirectoryShellCode();

    print F "if [ ! -e $wrk/3-overlapErrorAdjustment/\$jobid.oea ] ; then\n";
    print F "  \$bin/correctOverlaps \\\n";
    print F "    -G $wrk/$asm.gkpStore \\\n";
    print F "    -O $wrk/$asm.ovlStore \\\n";
    print F "    -R \$minid \$maxid \\\n";
    print F "    -e " . getGlobal("ovlErrorRate") . " -l " . getGlobal("minOverlapLength") . " \\\n";
    print F "    -c $wrk/3-overlapErrorAdjustment/red.red \\\n";
    print F "    -o $wrk/3-overlapErrorAdjustment/\$jobid.oea.WORKING \\\n";
    print F "  && \\\n";
    print F "  mv $wrk/3-overlapErrorAdjustment/\$jobid.oea.WORKING $wrk/3-overlapErrorAdjustment/\$jobid.oea\n";
    print F "fi\n";

    close(F);

    chmod 0755, "$wrk/3-overlapErrorAdjustment/oea.sh";

    emitStage($wrk, $asm, "overlapErrorAdjustmentConfigure");
}





sub overlapErrorAdjustmentCheck ($$$) {
    my $wrk     = shift @_;
    my $asm     = shift @_;
    my $attempt = shift @_;

    return  if (getGlobal("enableOEA") == 0);
    return  if (skipStage($wrk, $asm, "overlapErrorAdjustmentCheck", $attempt) == 1);
    return  if (-e "$wrk/3-overlapErrorAdjustment/oea.files");

    my $batchSize   = getGlobal("oeaBatchSize");
    my $failedJobs  = 0;
    my $numReads    = getNumberOfReadsInStore($wrk, $asm);
    my $numJobs     = int($numReads / $batchSize) + (($numReads % $batchSize == 0) ? 0 : 1);

    my $currentJobID   = 1;
    my @successJobs;
    my @failedJobs;
    my $failureMessage = "";

    for (my $currentJobID=1; $currentJobID <= $numJobs; $currentJobID++) {
        my $ji = substr("0000" . $currentJobID, -4);
        my $jn = "$wrk/3-overlapErrorAdjustment/$ji.oea";

        if (! -e $jn) {
            $failureMessage .= "  job $jn FAILED.\n";
            push @failedJobs, $currentJobID;
        } else {
            push @successJobs, $jn;
        }
    }

    #  No failed jobs?  Success!

    if (scalar(@failedJobs) == 0) {
        open(L, "> $wrk/3-overlapErrorAdjustment/oea.files") or caExit("can't open '$wrk/3-overlapErrorAdjustment/oea.files' for writing: $!", undef);
        foreach my $f (@successJobs) {
            print L "$f\n";
        }
        close(L);
        return;
    }

    #  If not the first attempt, report the jobs that failed, and that we're recomputing.

    if ($attempt > 1) {
        print STDERR "\n";
        print STDERR scalar(@failedJobs), " overlap error adjustment jobs failed:\n";
        print STDERR $failureMessage;
        print STDERR "\n";
    }

    #  If too many attempts, give up.

    if ($attempt > 2) {
        caExit("failed to adjust overlap error rates.  Made " . ($attempt-1) . " attempts, jobs still failed", undef);
    }

    #  Otherwise, run some jobs.

    print STDERR "overlapErrorAdjustmentCheck() -- attempt $attempt begins with ", scalar(@successJobs), " finished, and ", scalar(@failedJobs), " to compute.\n";

    emitStage($wrk, $asm, "overlapErrorAdjustmentCheck", $attempt);

    submitOrRunParallelJob($wrk, $asm, "oea", "$wrk/3-overlapErrorAdjustment", "oea", @failedJobs);
}




sub updateOverlapStore ($$) {
    my $wrk  = shift @_;
    my $asm  = shift @_;
    my $bin  = getBinDirectory();
    my $cmd;

    return  if (getGlobal("enableOEA") == 0);
    return  if (skipStage($wrk, $asm, "updateOverlapStore") == 1);
    return  if (-e "$wrk/$asm.ovlStore/erates");

    return if (-e "$wrk/$asm.ovlStore/adjustedErates");
    return if (-d "$wrk/$asm.tigStore");

    caExit("didn't find '$wrk/3-overlapErrorAdjustment/oea.files' to add to store, yet overlapper finished", undef)  if (! -e "$wrk/3-overlapErrorAdjustment/oea.files");

    $cmd  = "$bin/ovStoreBuild \\\n";
    $cmd .= "  -g $wrk/$asm.gkpStore \\\n";
    $cmd .= "  -o $wrk/$asm.ovlStore \\\n";
    $cmd .= "  -erates \\\n";
    $cmd .= "  -L $wrk/3-overlapErrorAdjustment/oea.files \\\n";
    $cmd .= "> $wrk/3-overlapErrorAdjustment/oea.apply.err 2>&1\n";

    if (runCommand("$wrk/3-overlapErrorAdjustment", $cmd)) {
        unlink "$wrk/$asm.ovlStore/erates";
        caExit("failed to add error rates to overlap store", "$wrk/3-overlapErrorAdjustment/oea.apply.err");
    }

    emitStage($wrk, $asm, "updateOverlapStore");
}
