
###############################################################################
 #
 #  This file is part of canu, a software program that assembles whole-genome
 #  sequencing reads into contigs.
 #
 #  This software is based on:
 #    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 #    the 'kmer package' (http://kmer.sourceforge.net)
 #  both originally distributed by Applera Corporation under the GNU General
 #  Public License, version 2.
 #
 #  Canu branched from Celera Assembler at its revision 4587.
 #  Canu branched from the kmer project at its revision 1994.
 #
 #  This file is derived from:
 #
 #    src/pipelines/ca3g/OverlapErrorAdjustment.pm
 #
 #  Modifications by:
 #
 #    Brian P. Walenz from 2015-FEB-27 to 2015-SEP-21
 #      are Copyright 2015 Battelle National Biodefense Institute, and
 #      are subject to the BSD 3-Clause License
 #
 #    Brian P. Walenz beginning on 2015-NOV-03
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #    Sergey Koren beginning on 2016-MAR-27
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

package canu::OverlapErrorAdjustment;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(readErrorDetectionConfigure readErrorDetectionCheck overlapErrorAdjustmentConfigure overlapErrorAdjustmentCheck updateOverlapStore);

use strict;

use File::Path qw(make_path remove_tree);

use canu::Defaults;
use canu::Execution;
use canu::Gatekeeper;
use canu::HTML;

#  Hardcoded to use utgOvlErrorRate



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
    my $WRK     = shift @_;           #  Root work directory (the -d option to canu)
    my $wrk     = "$WRK/unitigging";  #  Local work directory
    my $asm     = shift @_;
    my $bin     = getBinDirectory();
    my $path    = "$wrk/3-overlapErrorAdjustment";

    return         if (getGlobal("enableOEA") == 0);

    goto allDone   if (skipStage($wrk, $asm, "readErrorDetectionConfigure") == 1);
    goto allDone   if (-e "$path/red.red");

    goto allDone   if (-e "$wrk/$asm.ovlStore/adjustedEvalues");
    goto allDone   if (-d "$wrk/$asm.tigStore");

    make_path("$path")  if (! -d "$path");

    #  RED uses 13 bytes/base plus 12 bytes/overlap + space for evidence reads.

    my @readLengths;
    my @numOlaps;

    #print STDERR "$bin/gatekeeperDumpMetaData -G $wrk/$asm.gkpStore -reads\n";
    open(F, "$bin/gatekeeperDumpMetaData -G $wrk/$asm.gkpStore -reads |");
    while (<F>) {
        my @v = split '\s+', $_;
        $readLengths[$v[0]] = $v[2];
    }
    close(F);

    #print STDERR "$bin/ovStoreDump -G $wrk/$asm.gkpStore -O $wrk/$asm.ovlStore -d -counts\n";
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

    #getAllowedResources("", "red");

    my $maxID    = getNumberOfReadsInStore($wrk, $asm);
    my $maxMem   = getGlobal("redMemory") * 1024 * 1024 * 1024;
    my $maxReads = getGlobal("redBatchSize");
    my $maxBases = getGlobal("redBatchLength");

    print STDERR "\n";
    print STDERR "Configure RED for ", getGlobal("redMemory"), "gb memory with batches of at most ", ($maxReads > 0) ? $maxReads : "(unlimited)", " reads and ", ($maxBases > 0) ? $maxBases : "(unlimited)", " bases.\n";
    print STDERR "\n";

    my $reads    = 0;
    my $bases    = 0;
    my $olaps    = 0;

    my $coverage = getExpectedCoverage($wrk, $asm);

    push @bgn, 1;

    for (my $id = 1; $id <= $maxID; $id++) {
        $reads += 1;
        $bases += $readLengths[$id];
        $olaps += $numOlaps[$id];

        #  Guess how much extra memory used for overlapping reads.  Small genomes tend to load every read in the store,
        #  large genomes ... load repeats + 2 * coverage * bases in reads (times 2 for overlaps off of each end)

        my $memory = (13 * $bases) + (12 * $olaps) + (2 * $bases * $coverage);

        if ((($maxMem   > 0) && ($memory >= $maxMem * 0.75)) ||    #  Allow 25% slop (10% is probably sufficient)
            (($maxReads > 0) && ($reads  >= $maxReads))      ||
            (($maxBases > 0) && ($bases  >= $maxBases))      ||
            (($id == $maxID - 1))) {
            push @end, $id;

            printf(STDERR "RED job %3u from read %9u to read %9u - %7.3f GB for %7u reads - %7.3f GB for %9u olaps - %7.3f GB for evidence\n",
                   $nj + 1, $bgn[$nj], $end[$nj],
                   $memory / 1024 / 1024 / 1024, $reads,
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

    open(F, "> $path/red.sh") or caExit("can't open '$path/red.sh' for writing: $!", undef);

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

  finishStage:
    emitStage($WRK, $asm, "readErrorDetectionConfigure");
    buildHTML($WRK, $asm, "utg");

  allDone:
}





sub readErrorDetectionCheck ($$) {
    my $WRK     = shift @_;           #  Root work directory (the -d option to canu)
    my $wrk     = "$WRK/unitigging";  #  Local work directory
    my $asm     = shift @_;
    my $attempt = getGlobal("canuIteration");
    my $path    = "$wrk/3-overlapErrorAdjustment";

    return         if (getGlobal("enableOEA") == 0);
    goto allDone   if (skipStage($wrk, $asm, "readErrorDetectionCheck", $attempt) == 1);
    goto allDone   if (-e "$path/red.red");

    goto allDone   if (-e "$wrk/$asm.ovlStore/adjustedEvalues");
    goto allDone   if (-d "$wrk/$asm.tigStore");

    #  Figure out if all the tasks finished correctly.

    my @successJobs;
    my @failedJobs;
    my $failureMessage = "";

    open(A, "< $path/red.sh") or caExit("can't open '$path/red.sh' for reading: $!\n", undef);
    while (<A>) {
        if (m/if.*jobid\s+=\s+(\d+)\s+.*then/) {
            my $ji = substr("0000" . $1, -4);
            my $jn = "$path/$ji.red";

            if (! -e $jn) {
                $failureMessage .= "--   job $jn FAILED.\n";
                push @failedJobs, $1;
            } else {
                push @successJobs, $jn;
            }
        }
    }
    close(A);

    #  Failed jobs, retry.

    if (scalar(@failedJobs) > 0) {

        #  If not the first attempt, report the jobs that failed, and that we're recomputing.

        if ($attempt > 1) {
            print STDERR "--\n";
            print STDERR "-- ", scalar(@failedJobs), " read error detection jobs failed:\n";
            print STDERR $failureMessage;
            print STDERR "--\n";
        }

        #  If too many attempts, give up.

        if ($attempt > getGlobal("canuIterationMax")) {
            caExit("failed to detect errors in reads.  Made " . ($attempt-1) . " attempts, jobs still failed", undef);
        }

        #  Otherwise, run some jobs.

        print STDERR "-- read error detection attempt $attempt begins with ", scalar(@successJobs), " finished, and ", scalar(@failedJobs), " to compute.\n";

        emitStage($WRK, $asm, "readErrorDetectionCheck", $attempt);
        buildHTML($WRK, $asm, "utg");

        submitOrRunParallelJob($WRK, $asm, "red", "$path", "red", @failedJobs);
        return;
    }


  finishStage:
    print STDERR "-- Found ", scalar(@successJobs), " read error detection output files.\n";

    #  I didn't wan't to concat all the corrections, but it is _vastly_ easier to do so, compared to
    #  hacking correctOverlaps to handle multiple corrections files.  Plus, it is now really just a
    #  concat; before, the files needed to be parsed to strip off a header.

    concatOutput("$path/red.red", @successJobs);

    setGlobal("canuIteration", 0);
    emitStage($WRK, $asm, "readErrorDetectionCheck");
    buildHTML($WRK, $asm, "utg");
    stopAfter("red");

  allDone:
}





sub overlapErrorAdjustmentConfigure ($$) {
    my $WRK     = shift @_;           #  Root work directory (the -d option to canu)
    my $wrk     = "$WRK/unitigging";  #  Local work directory
    my $asm     = shift @_;
    my $bin     = getBinDirectory();
    my $path    = "$wrk/3-overlapErrorAdjustment";

    return         if (getGlobal("enableOEA") == 0);

    goto allDone   if (skipStage($wrk, $asm, "overlapErrorAdjustmentConfigure") == 1);
    goto allDone   if (-e "$path/oea.sh");

    goto allDone   if (-e "$wrk/$asm.ovlStore/adjustedEvalues");
    goto allDone   if (-d "$wrk/$asm.tigStore");

    #  OEA uses 1 byte/base + 8 bytes/adjustment + 28 bytes/overlap.  We don't know the number of adjustments, but that's
    #  basically error rate.  No adjustment is output for mismatches.

    my @readLengths;
    my @numOlaps;

    #print STDERR "$bin/gatekeeperDumpMetaData -G $wrk/$asm.gkpStore -reads\n";
    open(F, "$bin/gatekeeperDumpMetaData -G $wrk/$asm.gkpStore -reads |");
    while (<F>) {
        my @v = split '\s+', $_;
        $readLengths[$v[0]] = $v[2];
    }
    close(F);

    #print STDERR "$bin/ovStoreDump -G $wrk/$asm.gkpStore -O $wrk/$asm.ovlStore -d -counts\n";
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

    #getAllowedResources("", "oea");

    my $maxID    = getNumberOfReadsInStore($wrk, $asm);
    my $maxMem   = getGlobal("oeaMemory") * 1024 * 1024 * 1024;
    my $maxReads = getGlobal("oeaBatchSize");
    my $maxBases = getGlobal("oeaBatchLength");

    print STDERR "\n";
    print STDERR "Configure OEA for ", getGlobal("oeaMemory"), "gb memory with batches of at most ", ($maxReads > 0) ? $maxReads : "(unlimited)", " reads and ", ($maxBases > 0) ? $maxBases : "(unlimited)", " bases.\n";
    print STDERR "\n";

    my $reads    = 0;
    my $bases    = 0;
    my $olaps    = 0;

    my $coverage = getExpectedCoverage($wrk, $asm);
    my $corrSize = (-s "$path/red.red");

    push @bgn, 1;

    for (my $id = 1; $id <= $maxID; $id++) {
        $reads += 1;
        $bases += $readLengths[$id];
        $olaps += $numOlaps[$id];

        #  Hacked to attempt to estimate adjustment size better.  Olaps should only require 12 bytes each.

        my $memBases  = (1   * $bases);              #  Corrected reads for this batch
        my $memAdj1   = (8   * $corrSize) * 0.33;    #  Overestimate of the size of the indel adjustments needed (total size includes mismatches)
        my $memReads  = (32  * $reads);              #  Read data in the batch
        my $memOlaps  = (32  * $olaps);              #  Loaded overlaps
        my $memSeq    = (4   * 2097152);             #  two char arrays of 2*maxReadLen
        my $memAdj2   = (16  * 2097152);             #  two Adjust_t arrays of maxReadLen
        my $memWA     = (32  * 1048576);             #  Work area (16mb) and edit array (16mb)
        my $memMisc   = (256 * 1048576);             #  Work area (16mb) and edit array (16mb) and (192mb) slop

        my $memory = $memBases + $memAdj1 + $memReads + $memOlaps + $memSeq + $memAdj2 + $memWA + $memMisc;

        if ((($maxMem   > 0) && ($memory >= $maxMem * 0.75)) ||
            (($maxReads > 0) && ($reads  >= $maxReads))      ||
            (($maxBases > 0) && ($bases  >= $maxBases))      ||
            (($id == $maxID))) {
            push @end, $id;

            printf(STDERR "OEA job %3u from read %9u to read %9u - %4.1f bases + %4.1f adjusts + %4.1f reads + %4.1f olaps + %4.1f fseq/rseq + %4.1f fadj/radj + %4.1f work + %4.1f misc = %5.1f MB\n",
                   $nj + 1, $bgn[$nj], $end[$nj],
                   $memBases / 1024 / 1024,
                   $memAdj1  / 1024 / 1024,
                   $memReads / 1024 / 1024,
                   $memOlaps / 1024 / 1024,
                   $memSeq   / 1024 / 1024,
                   $memAdj2  / 1024 / 1024,
                   $memWA    / 1024 / 1024,
                   $memMisc  / 1024 / 1024,
                   $memory   / 1024 / 1024);

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
    print F "jobid=\$" . getGlobal("gridEngineTaskID") . "\n";
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

  finishStage:
    emitStage($WRK, $asm, "overlapErrorAdjustmentConfigure");
    buildHTML($WRK, $asm, "utg");

  allDone:
}





sub overlapErrorAdjustmentCheck ($$) {
    my $WRK     = shift @_;           #  Root work directory (the -d option to canu)
    my $wrk     = "$WRK/unitigging";  #  Local work directory
    my $asm     = shift @_;
    my $attempt = getGlobal("canuIteration");
    my $path    = "$wrk/3-overlapErrorAdjustment";

    return         if (getGlobal("enableOEA") == 0);

    goto allDone   if (skipStage($wrk, $asm, "overlapErrorAdjustmentCheck", $attempt) == 1);
    goto allDone   if (-e "$path/oea.files");

    #  Figure out if all the tasks finished correctly.

    my $batchSize   = getGlobal("oeaBatchSize");
    my $failedJobs  = 0;
    my $numReads    = getNumberOfReadsInStore($wrk, $asm);
    my $numJobs     = 0;  #int($numReads / $batchSize) + (($numReads % $batchSize == 0) ? 0 : 1);

    #  Need to read script to find number of jobs!

    my @successJobs;
    my @failedJobs;
    my $failureMessage = "";

    open(A, "< $path/oea.sh") or caExit("can't open '$path/oea.sh' for reading: $!\n", undef);
    while (<A>) {
        if (m/if.*jobid\s+=\s+(\d+)\s+.*then/) {
            my $ji = substr("0000" . $1, -4);
            my $jn = "$path/$ji.oea";

            if (! -e $jn) {
                $failureMessage .= "--   job $jn FAILED.\n";
                push @failedJobs, $1;
            } else {
                push @successJobs, $jn;
            }
        }
    }
    close(A);

    #  Failed jobs, retry.

    if (scalar(@failedJobs) > 0) {

        #  If not the first attempt, report the jobs that failed, and that we're recomputing.

        if ($attempt > 1) {
            print STDERR "--\n";
            print STDERR "-- ", scalar(@failedJobs), " overlap error adjustment jobs failed:\n";
            print STDERR $failureMessage;
            print STDERR "--\n";
        }

        #  If too many attempts, give up.

        if ($attempt > getGlobal("canuIterationMax")) {
            caExit("failed to adjust overlap error rates.  Made " . ($attempt-1) . " attempts, jobs still failed", undef);
        }

        #  Otherwise, run some jobs.

        print STDERR "-- overlap error adjustment attempt $attempt begins with ", scalar(@successJobs), " finished, and ", scalar(@failedJobs), " to compute.\n";

        emitStage($WRK, $asm, "overlapErrorAdjustmentCheck", $attempt);
        buildHTML($WRK, $asm, "utg");

        submitOrRunParallelJob($WRK, $asm, "oea", "$path", "oea", @failedJobs);
        return;
    }

  finishStage:
    print STDERR "-- Found ", scalar(@successJobs), " overlap error adjustment output files.\n";

    open(L, "> $path/oea.files") or caExit("can't open '$path/oea.files' for writing: $!", undef);
    foreach my $f (@successJobs) {
        print L "$f\n";
    }
    close(L);

    setGlobal("canuIteration", 0);
    emitStage($WRK, $asm, "overlapErrorAdjustmentCheck");
    buildHTML($WRK, $asm, "utg");
    stopAfter("oea");

  allDone:
}




sub updateOverlapStore ($$) {
    my $WRK     = shift @_;           #  Root work directory (the -d option to canu)
    my $wrk     = "$WRK/unitigging";  #  Local work directory
    my $asm     = shift @_;
    my $bin     = getBinDirectory();
    my $cmd;
    my $path    = "$wrk/3-overlapErrorAdjustment";

    return         if (getGlobal("enableOEA") == 0);

    goto allDone   if (skipStage($wrk, $asm, "updateOverlapStore") == 1);
    goto allDone   if (-e "$wrk/$asm.ovlStore/evalues");

    goto allDone   if (-e "$wrk/$asm.ovlStore/adjustedEvalues");
    goto allDone   if (-d "$wrk/$asm.tigStore");

    caExit("didn't find '$path/oea.files' to add to store, yet overlapper finished", undef)  if (! -e "$path/oea.files");

    $cmd  = "$bin/ovStoreBuild \\\n";
    $cmd .= "  -G $wrk/$asm.gkpStore \\\n";
    $cmd .= "  -O $wrk/$asm.ovlStore \\\n";
    $cmd .= "  -evalues \\\n";
    $cmd .= "  -L $path/oea.files \\\n";
    $cmd .= "> $path/oea.apply.err 2>&1";

    if (runCommand("$path", $cmd)) {
        unlink "$wrk/$asm.ovlStore/evalues";
        caExit("failed to add error rates to overlap store", "$path/oea.apply.err");
    }

  finishStage:
    emitStage($WRK, $asm, "updateOverlapStore");
    buildHTML($WRK, $asm, "utg");

  allDone:
}
