
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

use File::Path 2.08 qw(make_path remove_tree);

use canu::Defaults;
use canu::Execution;
use canu::Gatekeeper;
use canu::Report;
use canu::HTML;
use canu::Grid_Cloud;

#  Hardcoded to use utgOvlErrorRate



sub readErrorDetectionConfigure ($) {
    my $asm     = shift @_;
    my $bin     = getBinDirectory();
    my $path    = "unitigging/3-overlapErrorAdjustment";

    return         if (getGlobal("enableOEA") == 0);

    goto allDone   if (fileExists("$path/red.sh"));    #  Script exists
    goto allDone   if (fileExists("$path/red.red"));   #  Result exists

    goto allDone   if (skipStage($asm, "readErrorDetectionConfigure") == 1);

    goto allDone   if (fileExists("unitigging/$asm.ovlStore/evalues"));   #  Stage entrely finished
    goto allDone   if (-d "unitigging/$asm.ctgStore");                    #  Assembly finished

    make_path("$path")  if (! -d "$path");

    #  RED uses 13 bytes/base plus 12 bytes/overlap + space for evidence reads.

    my @readLengths;
    my @numOlaps;

    #print STDERR "$bin/gatekeeperDumpMetaData -G unitigging/$asm.gkpStore -reads\n";
    open(F, "$bin/gatekeeperDumpMetaData -G unitigging/$asm.gkpStore -reads |");
    while (<F>) {
        my @v = split '\s+', $_;
        $readLengths[$v[0]] = $v[2];
    }
    close(F);

    #  NEEDS OPTIMIZE - only need counts here, not the whole store

    fetchStore("unitigging/$asm.ovlStore");

    #print STDERR "$bin/ovStoreDump -G unitigging/$asm.gkpStore -O unitigging/$asm.ovlStore -d -counts\n";
    open(F, "$bin/ovStoreDump -G unitigging/$asm.gkpStore -O unitigging/$asm.ovlStore -d -counts |");
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

    my $maxID    = getNumberOfReadsInStore("unitigging", $asm);
    my $maxMem   = getGlobal("redMemory") * 1024 * 1024 * 1024;
    my $maxReads = getGlobal("redBatchSize");
    my $maxBases = getGlobal("redBatchLength");

    print STDERR "\n";
    print STDERR "Configure RED for ", getGlobal("redMemory"), "gb memory with batches of at most ", ($maxReads > 0) ? $maxReads : "(unlimited)", " reads and ", ($maxBases > 0) ? $maxBases : "(unlimited)", " bases.\n";
    print STDERR "\n";

    my $reads    = 0;
    my $bases    = 0;
    my $olaps    = 0;

    my $coverage = getExpectedCoverage("unitigging", $asm);

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
            (($id == $maxID))) {
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

    my $numReads    = getNumberOfReadsInStore("unitigging", $asm);

    open(F, "> $path/red.sh") or caExit("can't open '$path/red.sh' for writing: $!", undef);

    print F "#!" . getGlobal("shell") . "\n\n";
    print F "\n";
    print F getBinDirectoryShellCode();
    print F "\n";
    print F setWorkDirectoryShellCode($path);
    print F fetchStoreShellCode("unitigging/$asm.gkpStore", $path, "");
    print F fetchStoreShellCode("unitigging/$asm.ovlStore", $path, "");
    print F "\n";
    print F getJobIDShellCode();
    print F "\n";

    for (my $jj=1; $jj <= $nj; $jj++) {
        print F "if [ \$jobid = $jj ] ; then\n";
        print F "  minid=$bgn[$jj-1]\n";
        print F "  maxid=$end[$jj-1]\n";
        print F "fi\n";
    }

    print F "jobid=`printf %05d \$jobid`\n";
    print F "\n";
    print F "if [ -e ./\$jobid.red ] ; then\n";
    print F "  echo Job previously completed successfully.\n";
    print F "  exit\n";
    print F "fi\n";
    print F "\n";
    print F "\$bin/findErrors \\\n";
    print F "  -G ../$asm.gkpStore \\\n";
    print F "  -O ../$asm.ovlStore \\\n";
    print F "  -R \$minid \$maxid \\\n";
    print F "  -e " . getGlobal("utgOvlErrorRate") . " -l " . getGlobal("minOverlapLength") . " \\\n";
    print F "  -o ./\$jobid.red.WORKING \\\n";
    print F "  -t $numThreads \\\n";
    print F "&& \\\n";
    print F "mv ./\$jobid.red.WORKING ./\$jobid.red\n";
    print F "\n";
    print F stashFileShellCode("$path", "\$jobid.red", "");
    print F "\n";

    close(F);

    makeExecutable("$path/red.sh");
    stashFile("$path/red.sh");

  finishStage:
    emitStage($asm, "readErrorDetectionConfigure");
    buildHTML($asm, "utg");

  allDone:
}





sub readErrorDetectionCheck ($) {
    my $asm     = shift @_;
    my $attempt = getGlobal("canuIteration");
    my $path    = "unitigging/3-overlapErrorAdjustment";

    return         if (getGlobal("enableOEA") == 0);

    goto allDone   if (fileExists("$path/red.red"));       #  Output exists

    goto allDone   if (skipStage($asm, "readErrorDetectionCheck", $attempt) == 1);

    goto allDone   if (fileExists("unitigging/$asm.ovlStore/evalues"));   #  Stage entrely finished
    goto allDone   if (-d "unitigging/$asm.ctgStore");                    #  Assembly finished

    fetchFile("$path/red.sh");

    #  Figure out if all the tasks finished correctly.

    my @successJobs;
    my @failedJobs;
    my $failureMessage = "";

    open(A, "< $path/red.sh") or caExit("can't open '$path/red.sh' for reading: $!", undef);
    while (<A>) {
        if (m/if.*jobid\s+=\s+(\d+)\s+.*then/) {
            my $ji = substr("00000" . $1, -5);
            my $jn = "unitigging/3-overlapErrorAdjustment/$ji.red";

            if (! fileExists($jn)) {
                $failureMessage .= "--   job $ji.red FAILED.\n";
                push @failedJobs, $1;
            } else {
                push @successJobs, $jn;
            }
        }
    }
    close(A);

    #  Failed jobs, retry.

    if (scalar(@failedJobs) > 0) {

        #  If too many attempts, give up.

        if ($attempt >= getGlobal("canuIterationMax")) {
            print STDERR "--\n";
            print STDERR "-- Read error detection jobs failed, tried $attempt times, giving up.\n";
            print STDERR $failureMessage;
            print STDERR "--\n";
            caExit(undef, undef);
        }

        if ($attempt > 0) {
            print STDERR "--\n";
            print STDERR "-- Read error detection jobs failed, retry.\n";
            print STDERR $failureMessage;
            print STDERR "--\n";
        }

        #  Otherwise, run some jobs.

        emitStage($asm, "readErrorDetectionCheck", $attempt);
        buildHTML($asm, "utg");

        submitOrRunParallelJob($asm, "red", $path, "red", @failedJobs);
        return;
    }


  finishStage:
    print STDERR "-- Found ", scalar(@successJobs), " read error detection output files.\n";

    #  I didn't wan't to concat all the corrections, but it is _vastly_ easier to do so, compared to
    #  hacking correctOverlaps to handle multiple corrections files.  Plus, it is now really just a
    #  concat; before, the files needed to be parsed to strip off a header.

    open(O, "> $path/red.red") or caExit("can't open '$path/red.red' for writing: $!", undef);
    binmode(O);

    foreach my $f (@successJobs) {
        fetchFile($f);

        open(F, "< $f") or caExit("can't open '$f' for reading: $!", undef);
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

    stashFile("$path/red.red");

    foreach my $f (@successJobs) {
        unlink $f;
    }

    emitStage($asm, "readErrorDetectionCheck");
    buildHTML($asm, "utg");

  allDone:
}





sub overlapErrorAdjustmentConfigure ($) {
    my $asm     = shift @_;
    my $bin     = getBinDirectory();
    my $path    = "unitigging/3-overlapErrorAdjustment";

    return         if (getGlobal("enableOEA") == 0);

    goto allDone   if (fileExists("$path/oea.sh"));   #  Script exists

    goto allDone   if (skipStage($asm, "overlapErrorAdjustmentConfigure") == 1);

    goto allDone   if (fileExists("unitigging/$asm.ovlStore/evalues"));   #  Stage entrely finished
    goto allDone   if (-d "unitigging/$asm.ctgStore");                    #  Assembly finished

    #  OEA uses 1 byte/base + 8 bytes/adjustment + 28 bytes/overlap.  We don't know the number of adjustments, but that's
    #  basically error rate.  No adjustment is output for mismatches.

    my @readLengths;
    my @numOlaps;

    #print STDERR "$bin/gatekeeperDumpMetaData -G unitigging/$asm.gkpStore -reads\n";
    open(F, "$bin/gatekeeperDumpMetaData -G unitigging/$asm.gkpStore -reads |");
    while (<F>) {
        my @v = split '\s+', $_;
        $readLengths[$v[0]] = $v[2];
    }
    close(F);

    #print STDERR "$bin/ovStoreDump -G unitigging/$asm.gkpStore -O unitigging/$asm.ovlStore -d -counts\n";
    open(F, "$bin/ovStoreDump -G unitigging/$asm.gkpStore -O unitigging/$asm.ovlStore -d -counts |");
    while (<F>) {
        my @v = split '\s+', $_;
        $numOlaps[$v[0]] = $v[1];
    }
    close(F);

    #  Make an array of partitions, putting as many reads into each as will fit in the desired memory.

  tryOEAagain:
    my @bgn;   undef @bgn;
    my @end;   undef @end;
    my @log;   undef @log;

    my $nj = 0;

    my $maxID    = getNumberOfReadsInStore("unitigging", $asm);
    my $maxMem   = getGlobal("oeaMemory") * 1024 * 1024 * 1024;
    my $maxReads = getGlobal("oeaBatchSize");
    my $maxBases = getGlobal("oeaBatchLength");

    print STDERR "\n";
    print STDERR "Configure OEA for ", getGlobal("oeaMemory"), "gb memory with batches of at most ", ($maxReads > 0) ? $maxReads : "(unlimited)", " reads and ", ($maxBases > 0) ? $maxBases : "(unlimited)", " bases.\n";

    my $reads    = 0;
    my $bases    = 0;
    my $olaps    = 0;

    fetchFile("$path/red.red");

    my $coverage     = getExpectedCoverage("unitigging", $asm);
    my $corrSize     = (-s "$path/red.red");

    my $smallJobs    = 0;
    my $smallJobSize = 1024;

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

            $smallJobs++   if ($end[$nj] - $bgn[$nj] < $smallJobSize);

            push @log, sprintf("OEA job %3u from read %9u to read %9u - %4.1f bases + %4.1f adjusts + %4.1f reads + %4.1f olaps + %4.1f fseq/rseq + %4.1f fadj/radj + %4.1f work + %4.1f misc = %5.1f MB\n",
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

    #  If too many small jobs, increase memory and try again.  We'll allow any size jobs as long as
    #  there are 8 or less, but then demand there are at most 2 small jobs.

    if (($nj > 8) && ($smallJobs >= 2)) {
        my $curMem =            getGlobal("oeaMemory");
        my $newMem = int(1000 * getGlobal("oeaMemory") * 1.25) / 1000;

        print STDERR "  FAILED - configured $nj jobs, but $smallJobs jobs process $smallJobSize reads or less each.  Increasing memory from $curMem GB to $newMem GB.\n";

        setGlobal("oeaMemory", $newMem);

        goto tryOEAagain;
    }

    #  Report.

    print STDERR "Configured $nj jobs.\n";
    print STDERR "\n";

    foreach my $l (@log) {
        print STDERR $l;
    }

    print STDERR "\n";

    #  Dump a script

    open(F, "> $path/oea.sh") or caExit("can't open '$path/oea.sh' for writing: $!", undef);

    print F "#!" . getGlobal("shell") . "\n\n";
    print F "\n";
    print F getBinDirectoryShellCode();
    print F "\n";
    print F setWorkDirectoryShellCode($path);
    print F fetchStoreShellCode("unitigging/$asm.gkpStore", $path, "");
    print F fetchStoreShellCode("unitigging/$asm.ovlStore", $path, "");
    print F "\n";
    print F getJobIDShellCode();
    print F "\n";

    for (my $jj=1; $jj <= $nj; $jj++) {
        print F "if [ \$jobid = $jj ] ; then\n";
        print F "  minid=$bgn[$jj-1]\n";
        print F "  maxid=$end[$jj-1]\n";
        print F "fi\n";
    }

    print F "jobid=`printf %05d \$jobid`\n";
    print F "\n";
    print F "if [ -e ./\$jobid.oea ] ; then\n";
    print F "  echo Job previously completed successfully.\n";
    print F "  exit\n";
    print F "fi\n";
    print F "\n";
    print F fetchFileShellCode("unitigging/3-overlapErrorAdjustment", "red.red", "");
    print F "\n";
    print F "\$bin/correctOverlaps \\\n";
    print F "  -G ../$asm.gkpStore \\\n";
    print F "  -O ../$asm.ovlStore \\\n";
    print F "  -R \$minid \$maxid \\\n";
    print F "  -e " . getGlobal("utgOvlErrorRate") . " -l " . getGlobal("minOverlapLength") . " \\\n";
    print F "  -c ./red.red \\\n";
    print F "  -o ./\$jobid.oea.WORKING \\\n";
    print F "&& \\\n";
    print F "mv ./\$jobid.oea.WORKING ./\$jobid.oea\n";
    print F "\n";
    print F stashFileShellCode("$path", "\$jobid.oea", "");
    print F "\n";

    close(F);

    makeExecutable("$path/oea.sh");
    stashFile("$path/oea.sh");

  finishStage:
    emitStage($asm, "overlapErrorAdjustmentConfigure");
    buildHTML($asm, "utg");

  allDone:
}





sub overlapErrorAdjustmentCheck ($) {
    my $asm     = shift @_;
    my $attempt = getGlobal("canuIteration");
    my $path    = "unitigging/3-overlapErrorAdjustment";

    return         if (getGlobal("enableOEA") == 0);

    goto allDone   if (fileExists("$path/oea.files"));   #  Output exists

    goto allDone   if (skipStage($asm, "overlapErrorAdjustmentCheck", $attempt) == 1);

    goto allDone   if (fileExists("unitigging/$asm.ovlStore/evalues"));   #  Stage entrely finished
    goto allDone   if (-d "unitigging/$asm.ctgStore");                    #  Assembly finished

    #  Figure out if all the tasks finished correctly.

    my $batchSize   = getGlobal("oeaBatchSize");
    my $failedJobs  = 0;
    my $numReads    = getNumberOfReadsInStore("unitigging", $asm);
    my $numJobs     = 0;  #int($numReads / $batchSize) + (($numReads % $batchSize == 0) ? 0 : 1);

    #  Need to read script to find number of jobs!

    my @successJobs;
    my @failedJobs;
    my $failureMessage = "";

    fetchFile("$path/oea.sh");

    open(A, "< $path/oea.sh") or caExit("can't open '$path/oea.sh' for reading: $!", undef);
    while (<A>) {
        if (m/if.*jobid\s+=\s+(\d+)\s+.*then/) {
            my $ji = substr("00000" . $1, -5);

            if (! fileExists("unitigging/3-overlapErrorAdjustment/$ji.oea")) {
                $failureMessage .= "--   job $ji.oea FAILED.\n";
                push @failedJobs, $1;
            } else {
                push @successJobs, "./$ji.oea";
            }
        }
    }
    close(A);

    #  Failed jobs, retry.

    if (scalar(@failedJobs) > 0) {

        #  If too many attempts, give up.

        if ($attempt >= getGlobal("canuIterationMax")) {
            print STDERR "--\n";
            print STDERR "-- Overlap error adjustment jobs failed, tried $attempt times, giving up.\n";
            print STDERR $failureMessage;
            print STDERR "--\n";
            caExit(undef, undef);
        }

        if ($attempt > 0) {
            print STDERR "--\n";
            print STDERR "-- Overlap error adjustment jobs failed, retry.\n";
            print STDERR $failureMessage;
            print STDERR "--\n";
        }

        #  Otherwise, run some jobs.

        emitStage($asm, "overlapErrorAdjustmentCheck", $attempt);
        buildHTML($asm, "utg");

        submitOrRunParallelJob($asm, "oea", $path, "oea", @failedJobs);
        return;
    }

  finishStage:
    print STDERR "-- Found ", scalar(@successJobs), " overlap error adjustment output files.\n";

    open(L, "> $path/oea.files") or caExit("can't open '$path/oea.files' for writing: $!", undef);
    foreach my $f (@successJobs) {
        print L "$f\n";
    }
    close(L);

    stashFile("$path/oea.files");

    emitStage($asm, "overlapErrorAdjustmentCheck");
    buildHTML($asm, "utg");

  allDone:
}




sub updateOverlapStore ($) {
    my $asm     = shift @_;
    my $bin     = getBinDirectory();
    my $cmd;
    my $path    = "unitigging/3-overlapErrorAdjustment";

    return         if (getGlobal("enableOEA") == 0);

    goto allDone   if (skipStage($asm, "updateOverlapStore") == 1);

    goto allDone   if (fileExists("unitigging/$asm.ovlStore/evalues"));   #  Stage entrely finished
    goto allDone   if (-d "unitigging/$asm.ctgStore");                    #  Assembly finished

    fetchFile("unitigging/3-overlapErrorAdjustment/oea.files");

    caExit("didn't find '$path/oea.files' to add to store, yet jobs claim to be finished", undef)  if (! -e "$path/oea.files");

    open(F, "< $path/oea.files");
    while (<F>) {
        chomp;
        fetchFile("$path/$_");
    }
    close(F);

    fetchStore("unitigging/$asm.ovlStore");

    $cmd  = "$bin/ovStoreBuild \\\n";
    $cmd .= "  -G ../$asm.gkpStore \\\n";
    $cmd .= "  -O ../$asm.ovlStore \\\n";
    $cmd .= "  -evalues \\\n";
    $cmd .= "  -L ./oea.files \\\n";
    $cmd .= "> ./oea.apply.err 2>&1";

    if (runCommand($path, $cmd)) {
        unlink "unitigging/$asm.ovlStore/evalues";
        caExit("failed to add error rates to overlap store", "$path/oea.apply.err");
    }

    stashFile("unitigging/$asm.ovlStore/evalues");

    my $report = "-- No report available.\n";

    #open(F, "< $path/oea.apply.stats") or caExit("Failed to open error rate adjustment statistics in '$path/oea.apply.stats': $!", undef);
    #while (<F>) {
    #}
    #close(F);

    addToReport("adjustments", $report);

  finishStage:
    emitStage($asm, "updateOverlapStore");
    buildHTML($asm, "utg");

  allDone:
}
