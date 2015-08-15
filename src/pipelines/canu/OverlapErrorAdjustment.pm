package canu::OverlapErrorAdjustment;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(readErrorDetectionConfigure readErrorDetectionCheck overlapErrorAdjustmentConfigure overlapErrorAdjustmentCheck updateOverlapStore);

use strict;

use File::Path qw(make_path remove_tree);

use canu::Defaults;
use canu::Execution;
use canu::Gatekeeper;

#  Hardcoded to use utgOvlErrorRate



sub concatOutput ($@) {
    my $outName     = shift @_;
    my @successJobs =       @_;

    open(O, "> $outName");
    binmode(O);

    foreach my $f (@successJobs) {
        print STDERR "OverlapErrorAdjustment::concatOutput() $f -->\n";

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

    print STDERR "OverlapErrorAdjustment::concatOutput() --> $outName\n";

    #foreach my $f (@successJobs) {
    #    unlink $f;
    #}
}




sub readErrorDetectionConfigure ($$) {
    my $wrk   = shift @_;
    my $asm   = shift @_;
    my $bin   = getBinDirectory();
    my $path  = "$wrk/3-overlapErrorAdjustment";

    return  if (getGlobal("enableOEA") == 0);
    return  if (skipStage($wrk, $asm, "readErrorDetectionConfigure") == 1);
    return  if (-e "$path/red.red");

    return if (-e "$wrk/$asm.ovlStore/adjustedEvalues");
    return if (-d "$wrk/$asm.tigStore");

    make_path("$path")  if (! -d "$path");

    #  RED uses 13 bytes/base plus 12 bytes/overlap + space for evidence reads.

    my @readLengths;
    my @numOlaps;

    print STDERR "$bin/gatekeeperDumpMetaData -G $wrk/$asm.gkpStore -reads\n";
    open(F, "$bin/gatekeeperDumpMetaData -G $wrk/$asm.gkpStore -reads |");
    while (<F>) {
        my @v = split '\s+', $_;
        $readLengths[$v[0]] = $v[2];
    }
    close(F);

    print STDERR "$bin/ovStoreDump -G $wrk/$asm.gkpStore -O $wrk/$asm.ovlStore -d -counts\n";
    open(F, "$bin/ovStoreDump -G $wrk/$asm.gkpStore -O $wrk/$asm.ovlStore -d -counts |");
    while (<F>) {
        my @v = split '\s+', $_;
        $numOlaps[$v[0]] = $v[1];
    }
    close(F);

    #  Make an array of partitions, putting as many reads into each as will fit in the desired memory.

    my @bgn;
    my @end;
    my $nj = 0;

    my $maxID    = getNumberOfReadsInStore($wrk, $asm);
    my $maxMem   = getGlobal("redMemory") * 1024 * 1024 * 1024;
    my $maxReads = getGlobal("redBatchSize");
    my $maxBases = getGlobal("redBatchLength");

    my $reads    = 0;
    my $bases    = 0;
    my $olaps    = 0;

    my $coverage = getExpectedCoverage($wrk, $asm);

    push @bgn, 1;

    for (my $id = 1; $id < $maxID; $id++) {
        $reads += 1;
        $bases += $readLengths[$id];
        $olaps += $numOlaps[$id];

        #  Guess how much extra memory used for overlapping reads.  Small genomes tend to load every read in the store,
        #  large genomes ... load repeats + 2 * coverage * bases in reads (times 2 for overlaps off of each end)

        my $memory = (13 * $bases) + (12 * $olaps) + (2 * $bases * $coverage);

        if ((($maxMem   > 0) && ($memory >= $maxMem)) ||
            (($maxReads > 0) && ($reads  >= $maxReads)) ||
            (($maxBases > 0) && ($bases  >= $maxBases)) ||
            (($id == $maxID - 1))) {
            push @end, $id;

            printf(STDERR "RED job %3u from read %9u to read %9u using %7.3f GB for %7u reads, %7.3f GB for %9u olaps and %7.3f GB for evidence\n",
                   $nj + 1, $bgn[$nj], $end[$nj], $reads,
                   13 * $bases / 1024 / 1024 / 1024, $bases,
                   12 * $olaps / 1024 / 1024 / 1024, $olaps,
                   2 * $bases * $coverage / 1024 / 1024 / 1024);

            $nj++;

            $reads = 0;
            $bases = 0;
            $olaps = 0;

            push @bgn, $id + 1;  #  RED expects inclusive ranges.
        }
    }

    #  Dump a script.

    my $batchSize   = getGlobal("redBatchSize");
    my $numThreads  = getGlobal("redThreads");

    my $numReads    = getNumberOfReadsInStore($wrk, $asm);

    open(F, "> $path/red.sh") or caWarn("can't open '$path/red.sh' for writing: $!", undef);

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

    for (my $jj=1; $jj <= $nj; $jj++) {
        print F "if [ \$jobid = $jj ] ; then\n";
        print F "  minid=$bgn[$jj-1]\n";
        print F "  maxid=$end[$jj-1]\n";
        print F "fi\n";
    }

    print F "jobid=`printf %04d \$jobid`\n";
    print F "\n";
    print F "if [ -e $path/\$jobid.red ] ; then\n";
    print F "  echo Job previously completed successfully.\n";
    print F "  exit\n";
    print F "fi\n";

    print F getBinDirectoryShellCode();

    print F "if [ ! -e $path/\$jobid.red ] ; then\n";
    print F "  \$bin/findErrors \\\n";
    print F "    -G $wrk/$asm.gkpStore \\\n";
    print F "    -O $wrk/$asm.ovlStore \\\n";
    print F "    -R \$minid \$maxid \\\n";
    print F "    -e " . getGlobal("utgOvlErrorRate") . " -l " . getGlobal("minOverlapLength") . " \\\n";
    print F "    -o $path/\$jobid.red.WORKING \\\n";
    print F "    -t $numThreads \\\n";
    print F "  && \\\n";
    print F "  mv $path/\$jobid.red.WORKING $path/\$jobid.red\n";
    print F "fi\n";

    close(F);

    chmod 0755, "$path/red.sh";

    emitStage($wrk, $asm, "readErrorDetectionConfigure");
}





sub readErrorDetectionCheck ($$$) {
    my $wrk     = shift @_;
    my $asm     = shift @_;
    my $attempt = shift @_;
    my $path    = "$wrk/3-overlapErrorAdjustment";

    return  if (getGlobal("enableOEA") == 0);
    return  if (skipStage($wrk, $asm, "readErrorDetectionCheck", $attempt) == 1);
    return  if (-e "$path/red.red");

    return if (-e "$wrk/$asm.ovlStore/adjustedEvalues");
    return if (-d "$wrk/$asm.tigStore");

    my $numReads    = getNumberOfReadsInStore($wrk, $asm);
    my $numJobs     = 0;

    my $currentJobID   = 1;
    my @successJobs;
    my @failedJobs;
    my $failureMessage = "";

    #  Need to read script to find number of jobs!

    open(A, "< $path/red.sh") or caExit("can't open '$path/red.sh' for reading: $!\n", undef);

    while (<A>) {
        if (m/if.*jobid\s+=\s+(\d+)\s+.*then/) {
            $numJobs = ($numJobs < $1) ? $1 : $numJobs;
        }
    }

    close(A);

    for (my $currentJobID=1; $currentJobID <= $numJobs; $currentJobID++) {
        my $ji = substr("0000" . $currentJobID, -4);
        my $jn = "$path/$ji.red";

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
        concatOutput("$path/red.red", @successJobs);
        setGlobal("canuIteration", 0);
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

    submitOrRunParallelJob($wrk, $asm, "red", "$path", "red", @failedJobs);
}





sub overlapErrorAdjustmentConfigure ($$) {
    my $wrk   = shift @_;
    my $asm   = shift @_;
    my $bin   = getBinDirectory();
    my $path  = "$wrk/3-overlapErrorAdjustment";

    return  if (getGlobal("enableOEA") == 0);
    return  if (skipStage($wrk, $asm, "overlapErrorAdjustmentConfigure") == 1);
    return  if (-e "$path/oea.sh");

    return if (-e "$wrk/$asm.ovlStore/adjustedEvalues");
    return if (-d "$wrk/$asm.tigStore");

    #  OEA uses 1 byte/base + 8 bytes/adjustment + 28 bytes/overlap.  We don't know the number of adjustments, but that's
    #  basically error rate.  No adjustment is output for mismatches.

    my @readLengths;
    my @numOlaps;

    print STDERR "$bin/gatekeeperDumpMetaData -G $wrk/$asm.gkpStore -reads\n";
    open(F, "$bin/gatekeeperDumpMetaData -G $wrk/$asm.gkpStore -reads |");
    while (<F>) {
        my @v = split '\s+', $_;
        $readLengths[$v[0]] = $v[2];
    }
    close(F);

    print STDERR "$bin/ovStoreDump -G $wrk/$asm.gkpStore -O $wrk/$asm.ovlStore -d -counts\n";
    open(F, "$bin/ovStoreDump -G $wrk/$asm.gkpStore -O $wrk/$asm.ovlStore -d -counts |");
    while (<F>) {
        my @v = split '\s+', $_;
        $numOlaps[$v[0]] = $v[1];
    }
    close(F);

    #  Make an array of partitions, putting as many reads into each as will fit in the desired memory.

    my @bgn;
    my @end;
    my $nj = 0;

    my $maxID    = getNumberOfReadsInStore($wrk, $asm);
    my $maxMem   = getGlobal("oeaMemory") * 1024 * 1024 * 1024;
    my $maxReads = getGlobal("oeaBatchSize");
    my $maxBases = getGlobal("oeaBatchLength");

    my $reads    = 0;
    my $bases    = 0;
    my $olaps    = 0;

    my $coverage = getExpectedCoverage($wrk, $asm);

    push @bgn, 1;

    for (my $id = 1; $id < $maxID; $id++) {
        $reads += 1;
        $bases += $readLengths[$id];
        $olaps += $numOlaps[$id];

        my $memory = (1 * $bases) + (28 * $olaps) + (8 * $bases * getGlobal("errorRate"));

        if ((($maxMem   > 0) && ($memory >= $maxMem)) ||
            (($maxReads > 0) && ($reads  >= $maxReads)) ||
            (($maxBases > 0) && ($bases  >= $maxBases)) ||
            (($id == $maxID - 1))) {
            push @end, $id;

            printf(STDERR "OEA job %3u from read %9u to read %9u using %7.3f GB for reads, %7.3f GB for olaps and %7.3f GB for adjustments\n",
                   $nj + 1, $bgn[$nj], $end[$nj],
                   1 * $bases / 1024 / 1024 / 1024,
                   12 * $olaps / 1024 / 1024 / 1024,
                   8 * $bases * getGlobal("errorRate") / 1024 / 1024 / 1024);

            $nj++;

            $reads = 0;
            $bases = 0;
            $olaps = 0;

            push @bgn, $id + 1;  #  OEA expects inclusive ranges.
        }
    }

    #  Dump a script

    open(F, "> $path/oea.sh") or caExit("can't open '$path/oea.sh' for writing: $!", undef);

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

    for (my $jj=1; $jj <= $nj; $jj++) {
        print F "if [ \$jobid = $jj ] ; then\n";
        print F "  minid=$bgn[$jj-1]\n";
        print F "  maxid=$end[$jj-1]\n";
        print F "fi\n";
    }

    print F "jobid=`printf %04d \$jobid`\n";
    print F "\n";
    print F "if [ -e $path/\$jobid.oea ] ; then\n";
    print F "  echo Job previously completed successfully.\n";
    print F "  exit\n";
    print F "fi\n";

    print F getBinDirectoryShellCode();

    print F "if [ ! -e $path/\$jobid.oea ] ; then\n";
    print F "  \$bin/correctOverlaps \\\n";
    print F "    -G $wrk/$asm.gkpStore \\\n";
    print F "    -O $wrk/$asm.ovlStore \\\n";
    print F "    -R \$minid \$maxid \\\n";
    print F "    -e " . getGlobal("utgOvlErrorRate") . " -l " . getGlobal("minOverlapLength") . " \\\n";
    print F "    -c $path/red.red \\\n";
    print F "    -o $path/\$jobid.oea.WORKING \\\n";
    print F "  && \\\n";
    print F "  mv $path/\$jobid.oea.WORKING $path/\$jobid.oea\n";
    print F "fi\n";

    close(F);

    chmod 0755, "$path/oea.sh";

    emitStage($wrk, $asm, "overlapErrorAdjustmentConfigure");
}





sub overlapErrorAdjustmentCheck ($$$) {
    my $wrk     = shift @_;
    my $asm     = shift @_;
    my $attempt = shift @_;
    my $path    = "$wrk/3-overlapErrorAdjustment";

    return  if (getGlobal("enableOEA") == 0);
    return  if (skipStage($wrk, $asm, "overlapErrorAdjustmentCheck", $attempt) == 1);
    return  if (-e "$path/oea.files");

    my $batchSize   = getGlobal("oeaBatchSize");
    my $failedJobs  = 0;
    my $numReads    = getNumberOfReadsInStore($wrk, $asm);
    my $numJobs     = 0;  #int($numReads / $batchSize) + (($numReads % $batchSize == 0) ? 0 : 1);

    #  Need to read script to find number of jobs!

    open(A, "< $path/oea.sh") or caExit("can't open '$path/oea.sh' for reading: $!\n", undef);

    while (<A>) {
        if (m/if.*jobid\s+=\s+(\d+)\s+.*then/) {
            $numJobs = ($numJobs < $1) ? $1 : $numJobs;
        }
    }

    close(A);

    my $currentJobID   = 1;
    my @successJobs;
    my @failedJobs;
    my $failureMessage = "";

    for (my $currentJobID=1; $currentJobID <= $numJobs; $currentJobID++) {
        my $ji = substr("0000" . $currentJobID, -4);
        my $jn = "$path/$ji.oea";

        if (! -e $jn) {
            $failureMessage .= "  job $jn FAILED.\n";
            push @failedJobs, $currentJobID;
        } else {
            push @successJobs, $jn;
        }
    }

    #  No failed jobs?  Success!

    if (scalar(@failedJobs) == 0) {
        open(L, "> $path/oea.files") or caExit("can't open '$path/oea.files' for writing: $!", undef);
        foreach my $f (@successJobs) {
            print L "$f\n";
        }
        close(L);
        setGlobal("canuIteration", 0);
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

    submitOrRunParallelJob($wrk, $asm, "oea", "$path", "oea", @failedJobs);
}




sub updateOverlapStore ($$) {
    my $wrk  = shift @_;
    my $asm  = shift @_;
    my $bin  = getBinDirectory();
    my $cmd;
    my $path = "$wrk/3-overlapErrorAdjustment";

    return  if (getGlobal("enableOEA") == 0);
    return  if (skipStage($wrk, $asm, "updateOverlapStore") == 1);
    return  if (-e "$wrk/$asm.ovlStore/evalues");

    return if (-e "$wrk/$asm.ovlStore/adjustedEvalues");
    return if (-d "$wrk/$asm.tigStore");

    caExit("didn't find '$path/oea.files' to add to store, yet overlapper finished", undef)  if (! -e "$path/oea.files");

    $cmd  = "$bin/ovStoreBuild \\\n";
    $cmd .= "  -g $wrk/$asm.gkpStore \\\n";
    $cmd .= "  -o $wrk/$asm.ovlStore \\\n";
    $cmd .= "  -evalues \\\n";
    $cmd .= "  -L $path/oea.files \\\n";
    $cmd .= "> $path/oea.apply.err 2>&1\n";

    if (runCommand("$path", $cmd)) {
        unlink "$wrk/$asm.ovlStore/evalues";
        caExit("failed to add error rates to overlap store", "$path/oea.apply.err");
    }

    emitStage($wrk, $asm, "updateOverlapStore");
}
