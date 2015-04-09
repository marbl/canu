package ca3g::OverlapErrorAdjustment;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(overlapErrorAdjustment);

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

    return  if (-e "$wrk/3-overlapErrorAdjustment/red.red");

    my $batchSize   = getGlobal("redBatchSize");
    my $numThreads  = getGlobal("redThreads");

    my $numReads    = getNumberOfReadsInStore($wrk, $asm);
    my $numJobs     = int($numReads / $batchSize) + (($numReads % $batchSize == 0) ? 0 : 1);

    open(F, "> $wrk/3-overlapErrorAdjustment/red.sh") or caFailure("failed to write to '$wrk/3-overlapErrorAdjustment/red.sh'", undef);

    print F "#!" . getGlobal("shell") . "\n\n";
    print F "\n";
    print F "jobid=\$" . getGlobal("gridEngineTaskID") . "\n";
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
    print F "    -e " . getGlobal("ovlErrorRate") . " -l " . getGlobal("ovlMinLen") . " \\\n";
    print F "    -o $wrk/3-overlapErrorAdjustment/\$jobid.red.WORKING \\\n";
    print F "    -t $numThreads \\\n";
    print F "  && \\\n";
    print F "  mv $wrk/3-overlapErrorAdjustment/\$jobid.red.WORKING $wrk/3-overlapErrorAdjustment/\$jobid.red\n";
    print F "fi\n";

    close(F);

    chmod 0755, "$wrk/3-overlapErrorAdjustment/red.sh";

    submitOrRunParallelJob($wrk, $asm, "oea", "$wrk/3-overlapErrorAdjustment", "red", getGlobal("redConcurrency"), "1-$numJobs");
}





sub readErrorDetectionCheck ($$$) {
    my $wrk     = shift @_;
    my $asm     = shift @_;
    my $attempt = shift @_;

    return  if (-e "$wrk/3-overlapErrorAdjustment/red.red");

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

    #  I didn't wan't to concat all the corrections, but it is _vastly_ easier to do so, compared to
    #  hacking correctOverlaps to handle multiple corrections files.  Plus, it is now really just a
    #  concat; before, the files needed to be parsed to strip off a header.

    if (scalar(@failedJobs) == 0) {
        concatOutput("$wrk/3-overlapErrorAdjustment/red.red", @successJobs);
        return;
    }

    print STDERR "\n";
    print STDERR scalar(@failedJobs), " read error detection jobs failed:\n";
    print STDERR $failureMessage;
    print STDERR "\n";

    if ($attempt < 1) {
        submitOrRunParallelJob($wrk, $asm, "oea", "$wrk/3-overlapErrorAdjustment", "red", getGlobal("redConcurrency"), @failedJobs);
    } else {
        caFailure("failed to detect errors in reads.  Made $attempt attempts, jobs still failed", undef);
    }
}





sub overlapErrorAdjustmentConfigure ($$) {
    my $wrk   = shift @_;
    my $asm   = shift @_;

    my $batchSize   = getGlobal("oeaBatchSize");
    my $numReads    = getNumberOfReadsInStore($wrk, $asm);
    my $numJobs     = int($numReads / $batchSize) + (($numReads % $batchSize == 0) ? 0 : 1);

    open(F, "> $wrk/3-overlapErrorAdjustment/oea.sh") or caFailure("failed to write '$wrk/3-overlapErrorAdjustment/oea.sh'", undef);

    print F "#!" . getGlobal("shell") . "\n\n";
    print F "\n";
    print F "jobid=\$" . getGlobal("gridEngineTaskID") . "\n";
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
    print F "    -e " . getGlobal("ovlErrorRate") . " -l " . getGlobal("ovlMinLen") . " \\\n";
    print F "    -c $wrk/3-overlapErrorAdjustment/red.red \\\n";
    print F "    -o $wrk/3-overlapErrorAdjustment/\$jobid.oea.WORKING \\\n";
    print F "  && \\\n";
    print F "  mv $wrk/3-overlapErrorAdjustment/\$jobid.oea.WORKING $wrk/3-overlapErrorAdjustment/\$jobid.oea\n";
    print F "fi\n";

    close(F);

    chmod 0755, "$wrk/3-overlapErrorAdjustment/oea.sh";

    submitOrRunParallelJob($wrk, $asm, "oea", "$wrk/3-overlapErrorAdjustment", "oea", getGlobal("oeaConcurrency"), "1-$numJobs");
}





sub overlapErrorAdjustmentCheck ($$$) {
    my $wrk     = shift @_;
    my $asm     = shift @_;
    my $attempt = shift @_;

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

    if (scalar(@failedJobs) == 0) {
        open(L, "> $wrk/3-overlapErrorAdjustment/oea.files") or caFailure("failed to open '$wrk/3-overlapErrorAdjustment/oea.files'", undef);
        foreach my $f (@successJobs) {
            print L "$f\n";
        }
        close(L);
        return;
    }

    #if (scalar(@failedJobs) == 0) {
    #    concatOutput("$wrk/3-overlapErrorAdjustment/oea.oea", @successJobs);
    #    return;
    #}

    print STDERR "\n";
    print STDERR scalar(@failedJobs), " overlap error adjustment jobs failed:\n";
    print STDERR $failureMessage;
    print STDERR "\n";

    if ($attempt < 1) {
        submitOrRunParallelJob($wrk, $asm, "oea", "$wrk/3-overlapErrorAdjustment", "oea", getGlobal("oeaConcurrency"), @failedJobs);
    } else {
        caFailure("failed to adjust overlap error rates.  Made $attempt attempts, jobs still failed", undef);
    }
}






sub updateOverlapStore ($$) {
    my $wrk  = shift @_;
    my $asm  = shift @_;
    my $bin  = getBinDirectory();
    my $cmd;

    return  if (-e "$wrk/$asm.ovlStore/erates");

    caFailure("didn't find '$wrk/3-overlapErrorAdjustment/oea.files' to add to store", undef)  if (! -e "$wrk/3-overlapErrorAdjustment/oea.files");

    $cmd  = "$bin/ovStoreBuild \\\n";
    $cmd .= "  -g $wrk/$asm.gkpStore \\\n";
    $cmd .= "  -o $wrk/$asm.ovlStore \\\n";
    $cmd .= "  -erates \\\n";
    $cmd .= "  -L $wrk/3-overlapErrorAdjustment/oea.files \\\n";
    $cmd .= "> $wrk/3-overlapErrorAdjustment/oea.apply.err 2>&1\n";

    if (runCommand("$wrk/3-overlapErrorAdjustment", $cmd)) {
        unlink "$wrk/$asm.ovlStore/erates";
        caFailure("failed to add error rates to overlap store", "$wrk/3-overlapErrorAdjustment/oea.apply.err");
    }
}





#  overlap error adjustment is two steps:
#    read error deteciton
#    recompute each overlap with read errors fixed
#
#  The original version (mostly untouched since Celera) needed to merge the outputs of each step,
#  and the last step would essentially rewrite the ovlStore to update error rates.  It was ugly.
#
sub overlapErrorAdjustment ($$) {
    my $wrk  = shift @_;
    my $asm  = shift @_;
    my $bin  = getBinDirectory();
    my $cmd;

    return if (getGlobal("enableOEA") == 0);

    return if (-e "$wrk/$asm.ovlStore/adjustedErates");
    return if (-d "$wrk/$asm.tigStore");

    make_path("$wrk/3-overlapErrorAdjustment")  if (! -d "$wrk/3-overlapErrorAdjustment");

    readErrorDetectionConfigure($wrk, $asm);
    readErrorDetectionCheck($wrk, $asm, 0);
    readErrorDetectionCheck($wrk, $asm, 1);

    overlapErrorAdjustmentConfigure($wrk, $asm);
    overlapErrorAdjustmentCheck($wrk, $asm, 0);
    overlapErrorAdjustmentCheck($wrk, $asm, 1);

    updateOverlapStore($wrk, $asm);
}

