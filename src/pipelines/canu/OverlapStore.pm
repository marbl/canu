
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
 #    src/pipelines/ca3g/OverlapStore.pm
 #
 #  Modifications by:
 #
 #    Brian P. Walenz from 2015-FEB-27 to 2015-SEP-21
 #      are Copyright 2015 Battelle National Biodefense Institute, and
 #      are subject to the BSD 3-Clause License
 #
 #    Brian P. Walenz beginning on 2015-OCT-10
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #    Sergey Koren beginning on 2015-DEC-08
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

package canu::OverlapStore;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(createOverlapStore);

use strict;
use File::Basename;   #  dirname
use File::Path 2.08 qw(make_path remove_tree);

use POSIX qw(ceil);
use canu::Defaults;
use canu::Execution;
use canu::Gatekeeper;
use canu::Report;
use canu::Grid_Cloud;


#  Parallel documentation: Each overlap job is converted into a single bucket of overlaps.  Within
#  each bucket, the overlaps are distributed into many slices, one per sort job.  The sort jobs then
#  load the same slice from each bucket.


#  NOT FILTERING overlaps by error rate when building the parallel store.


sub createOverlapStoreSequential ($$$) {
    my $base    = shift @_;
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $bin     = getBinDirectory();
    my $cmd;

    #  Fetch inputs.  If you're not cloud-based, this does nothing.  Really.  Trust me.

    fetchFile("$base/1-overlapper/ovljob.files");

    open(F, "< $base/1-overlapper/ovljob.files") or caExit("failed to open overlapper output list in '$base/1-overlapper/ovljob.files'", undef);
    while (<F>) {
        chomp;

        if (m/^(.*).ovb$/) {
            fetchFile("$base/$1.ovb");
            fetchFile("$base/$1.counts");
        } else {
            caExit("didn't recognize ovljob.files line '$_'", undef);
        }
    }
    close(F);

    #  This is running in the canu process itself.  Execution.pm has special case code
    #  to submit canu to grids using the maximum of 4gb and this memory limit.

    $cmd  = "$bin/ovStoreBuild \\\n";
    $cmd .= " -O ./$asm.ovlStore.BUILDING \\\n";
    $cmd .= " -G ./$asm.gkpStore \\\n";
    $cmd .= " -C ./$asm.ovlStore.config \\\n";
    $cmd .= " > ./$asm.ovlStore.err 2>&1";

    if (runCommand($base, $cmd)) {
        caExit("failed to create the overlap store", "$base/$asm.ovlStore.err");
    }

    unlink("$base/$asm.ovlStore.err");

    rename("$base/$asm.ovlStore.BUILDING", "$base/$asm.ovlStore");

    stashStore("$base/$asm.ovlStore");
}



sub createOverlapStoreParallel ($$$$$$) {
    my $base       = shift @_;
    my $asm        = shift @_;
    my $tag        = shift @_;
    my $numBuckets = shift @_;
    my $numSlices  = shift @_;
    my $sortMemory = shift @_;
    my $bin        = getBinDirectory();
    my $cmd;
    my $path       = "$base/$asm.ovlStore.BUILDING";

    goto allDone   if (skipStage($asm, "$tag-overlapStoreConfigure") == 1);
    goto allDone   if ((-e "$path/scripts/1-bucketize.sh") &&
                       (-e "$path/scripts/2-sort.sh"));
    goto allDone   if (-d "$base/$asm.ovlStore");


    my $path       = "$base/$asm.ovlStore.BUILDING";

    make_path("$path/scripts");
    make_path("$path/logs");

    #  Parallel jobs for bucketizing.

    if (! -e "$path/scripts/1-bucketize.sh") {
        open(F, "> $path/scripts/1-bucketize.sh") or die;
        print F "#!" . getGlobal("shell") . "\n";
        print F "\n";
        print F getBinDirectoryShellCode();
        print F "\n";
        print F setWorkDirectoryShellCode($base);
        print F "\n";
        print F getJobIDShellCode();
        print F "\n";
        print F getLimitShellCode();
        print F "\n";
        print F "#  This script should be executed from $path/, but the binary needs\n";
        print F "#  to run from $base/ (all the paths in the config are relative to there).\n";
        print F "\n";
        print F "cd ..\n";
        print F "\n";
        print F "\$bin/ovStoreBucketizer \\\n";
        print F "  -O ./$asm.ovlStore.BUILDING \\\n";
        print F "  -G ./$asm.gkpStore \\\n";
        print F "  -C ./$asm.ovlStore.config \\\n";
        print F "  -b \$jobid \n";
        #print F "  -e " . getGlobal("") . " \\\n"  if (defined(getGlobal("")));
        close(F);
    }

    #  Parallel jobs for sorting each bucket

    if (! -e "$path/scripts/2-sort.sh") {
        open(F, "> $path/scripts/2-sort.sh") or die;
        print F "#!" . getGlobal("shell") . "\n";
        print F "\n";
        print F getBinDirectoryShellCode();
        print F "\n";
        print F setWorkDirectoryShellCode($base);
        print F "\n";
        print F getJobIDShellCode();
        print F "\n";
        print F getLimitShellCode();
        print F "\n";
        print F "#  This script should be executed from $path/, but the binary needs\n";
        print F "#  to run from $base/ (all the paths in the config are relative to there).\n";
        print F "\n";
        print F "cd ..\n";
        print F "\n";
        print F "\$bin/ovStoreSorter \\\n";
        #print F "  -deletelate \\\n";  #  Choices -deleteearly -deletelate or nothing
        print F "  -O ./$asm.ovlStore.BUILDING \\\n";
        print F "  -G ./$asm.gkpStore \\\n";
        print F "  -C ./$asm.ovlStore.config \\\n";
        print F "  -s \$jobid \\\n";
        print F "  -M $sortMemory \n";
        close(F);
    }

    #  A final job to merge the indices.

    makeExecutable("$path/scripts/1-bucketize.sh");
    makeExecutable("$path/scripts/2-sort.sh");

    stashFile("$path/scripts/1-bucketize.sh");
    stashFile("$path/scripts/2-sort.sh");

  finishStage:
    emitStage($asm, "$tag-overlapStoreConfigure");

  allDone:
    stopAfter("overlapStoreConfigure");
}



sub overlapStoreBucketizerCheck ($$$$$) {
    my $base       = shift @_;
    my $asm        = shift @_;
    my $tag        = shift @_;
    my $numBuckets = shift @_;
    my $numSlices  = shift @_;
    my $attempt    = getGlobal("canuIteration");
    my $path       = "$base/$asm.ovlStore.BUILDING";

    goto allDone   if (skipStage($asm, "$tag-overlapStoreBucketizerCheck", $attempt) == 1);
    goto allDone   if (-d "$base/$asm.ovlStore");
    goto allDone   if (-e "$path/1-bucketize.success");

    fetchFile("scripts/1-bucketize/1-bucketize.sh");

    #  Figure out if all the tasks finished correctly.

    my $currentJobID   = 1;
    my $bucketID       = "0001";

    my @successJobs;
    my @failedJobs;
    my $failureMessage = "";

    #  Two ways to check for completeness, either 'sliceSizes' exists, or the 'bucket' directory
    #  exists.  The compute is done in a 'create' directory, which is renamed to 'bucket' just
    #  before the job completes.

    for (my $bb=1; $bb<=$numBuckets; $bb++) {
        if (! -e "$path/bucket$bucketID") {
            $failureMessage .= "--   job $path/bucket$bucketID FAILED.\n";
            push @failedJobs, $currentJobID;
        } else {
            push @successJobs, $currentJobID;
        }

        $currentJobID++;
        $bucketID++;
    }

    close(F);

    #  Failed jobs, retry.

    if (scalar(@failedJobs) > 0) {

        #  If too many attempts, give up.

        if ($attempt >= getGlobal("canuIterationMax")) {
            print STDERR "--\n";
            print STDERR "-- Overlap store bucketizer jobs failed, tried $attempt times, giving up.\n";
            print STDERR $failureMessage;
            print STDERR "--\n";
            caExit(undef, undef);
        }

        if ($attempt > 0) {
            print STDERR "--\n";
            print STDERR "-- Overlap store bucketizer jobs failed, retry.\n";
            print STDERR $failureMessage;
            print STDERR "--\n";
        }

        #  Otherwise, run some jobs.

        emitStage($asm, "$tag-overlapStoreBucketizerCheck", $attempt);

        submitOrRunParallelJob($asm, "ovB", $path, "scripts/1-bucketize", @failedJobs);
        return;
    }

  finishStage:
    print STDERR "-- Overlap store bucketizer finished.\n";

    touch("$path/1-bucketize.success");

    emitStage($asm, "$tag-overlapStoreBucketizerCheck");

  allDone:
}





sub overlapStoreSorterCheck ($$$$$) {
    my $base       = shift @_;
    my $asm        = shift @_;
    my $tag        = shift @_;
    my $numBuckets = shift @_;
    my $numSlices  = shift @_;
    my $attempt    = getGlobal("canuIteration");
    my $path       = "$base/$asm.ovlStore.BUILDING";

    goto allDone   if (skipStage($asm, "$tag-overlapStoreSorterCheck", $attempt) == 1);
    goto allDone   if (-d "$base/$asm.ovlStore");
    goto allDone   if (-e "$path/2-sorter.success");

    fetchFile("scripts/1-bucketize/2-sort.sh");

    #  Figure out if all the tasks finished correctly.

    my $currentJobID   = 1;
    my $sliceID        = "0001";

    my @successJobs;
    my @failedJobs;
    my $failureMessage = "";

    #  A valid result has three files:
    #    $path/$sliceID<001>
    #    $path/$sliceID.index
    #    $path/$sliceID.info
    #
    #  A crashed result has one file, if it crashes before output
    #    $path/$sliceID.ovs  (PROBABLY NOT TRUE AS OF MARCH 2018)
    #
    #  On out of disk, the .info is missing.  It's the last thing created.
    #
    for (my $ss=1; $ss<=$numSlices; $ss++) {
        if ((! -e "$path/$sliceID<001>") ||
            (! -e "$path/$sliceID.info") ||
            (  -e "$path/$sliceID.ovs")) {
            $failureMessage .= "--   job $path/$sliceID FAILED.\n";
            unlink "$path/$sliceID.ovs";
            push @failedJobs, $currentJobID;
        } else {
            push @successJobs, $currentJobID;
        }

        $currentJobID++;
        $sliceID++;
    }

    close(F);

    #  Failed jobs, retry.

    if (scalar(@failedJobs) > 0) {

        #  If too many attempts, give up.


        if ($attempt >= getGlobal("canuIterationMax")) {
            print STDERR "--\n";
            print STDERR "-- Overlap store sorting jobs failed, tried $attempt times, giving up.\n";
            print STDERR $failureMessage;
            print STDERR "--\n";
            caExit(undef, undef);
        }

        if ($attempt > 0) {
            print STDERR "--\n";
            print STDERR "-- Overlap store sorting jobs failed, retry.\n";
            print STDERR $failureMessage;
            print STDERR "--\n";
        }

        #  Otherwise, run some jobs.

        emitStage($asm, "$tag-overlapStoreSorterCheck", $attempt);

        submitOrRunParallelJob($asm, "ovS", $path, "scripts/2-sort", @failedJobs);
        return;
    }

  finishStage:
    print STDERR "-- Overlap store sorter finished.\n";

    touch("$path/2-sorter.success");

    emitStage($asm, "$tag-overlapStoreSorterCheck");

  allDone:
}




sub checkOverlapStore ($$) {
    my $base    = shift @_;
    my $asm     = shift @_;

    my $bin   = getBinDirectory();
    my $cmd;

    $cmd  = "$bin/ovStoreDump \\\n";
    $cmd .= " -G ./$asm.gkpStore \\\n";
    $cmd .= " -O ./$asm.ovlStore \\\n";
    $cmd .= " -counts \\\n";
    $cmd .= " > ./$asm.ovlStore/counts.dat 2> ./$asm.ovlStore/counts.err";

    print STDERR "-- Checking store.\n";

    if (runCommand($base, $cmd)) {
        caExit("failed to dump counts of overlaps; invalid store?", "$base/$asm.ovlStore/counts.err");
    }

    my $totOvl   = 0;
    my $nulReads = 0;
    my $ovlReads = 0;

    open(F, "< ./$base/$asm.ovlStore/counts.dat") or die;
    while (<F>) {
        my @v = split '\s+', $_;

        $nulReads += 1       if ($v[1] < 1);
        $ovlReads += 1       if ($v[1] > 0);
        $totOvl   += $v[1];
    }
    close(F);

    print STDERR "--\n";
    print STDERR "-- Overlap store '$base/$asm.ovlStore' successfully constructed.\n";
    print STDERR "-- Found $totOvl overlaps for $ovlReads reads; $nulReads reads have no overlaps.\n";
    print STDERR "--\n";

    unlink "./$base/$asm.ovlStore/counts.dat";
    unlink "./$base/$asm.ovlStore/counts.err";
}



#  Delete the inputs and directories.
#
#    Directories - Viciously remove the whole thing (after all files are deleted, so we
#                  can get the sizes).
#    Files       - Sum the size, remove the file, and try to remove the directory.  In
#                  particular, we don't want to remove_tree() this directory - there could
#                  be other stuff in it - only remove if empty.
#
#  Ideally, every directory we have in our list should be empty after we delete the files in the
#  list.  But they won't be.  Usually because there are empty directories in there too.  Maybe
#  some stray files we didn't track.  Regardless, just blow them away.
#
#  Previous (to July 2017) versions tried to gently rmdir things, but it was ugly and didn't
#  quite work.
#
sub deleteOverlapIntermediateFiles ($$) {
    my $base    = shift @_;
    my $asm     = shift @_;

    my %directories;
    my $bytes = 0;
    my $files = 0;

    foreach my $file ("$base/1-overlapper/ovljob.files",
                      "$base/1-overlapper/ovljob.more.files",
                      "$base/1-overlapper/mhap.files",
                      "$base/1-overlapper/mmap.files",
                      "$base/1-overlapper/precompute.files") {
        next  if (! -e $file);

        open(F, "< $file") or caExit("can't open '$file' for reading: $!\n", undef);
        while (<F>) {
            chomp;

            if (-d "$base/$_") {
                $directories{$_}++;

            } elsif (-e "$base/$_") {
                $bytes += -s "$base/$_";
                $files += 1;

                unlink "$base/$_";
                rmdir dirname("$base/$_");  #  Try to rmdir the directory the file is in.  If empty, yay!
            }
        }
        close(F);
    }

    foreach my $dir (keys %directories) {
        remove_tree("$base/$dir");
    }

    unlink "$base/1-overlapper/ovljob.files";
    unlink "$base/1-overlapper/ovljob.more.files";
    unlink "$base/1-overlapper/mhap.files";
    unlink "$base/1-overlapper/mmap.files";
    unlink "$base/1-overlapper/precompute.files";

    print STDERR "--\n";
    print STDERR "-- Purged ", int(1000 * $bytes / 1024 / 1024 / 1024) / 1000, " GB in $files overlap output files.\n";
}



sub createOverlapStore ($$$) {
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $seq     = shift @_;
    my $cmd;
    my $bin     = getBinDirectory();

    my $base;

    $base = "correction"  if ($tag eq "cor");
    $base = "trimming"    if ($tag eq "obt");
    $base = "unitigging"  if ($tag eq "utg");

    goto allDone   if (skipStage($asm, "$tag-createOverlapStore") == 1);
    goto allDone   if ((-d "$base/$asm.ovlStore") || (fileExists("$base/$asm.ovlStore.tar")));
    goto allDone   if ((-d "$base/$asm.ctgStore") || (fileExists("$base/$asm.ctgStore.tar")));

    #  Did we _really_ complete?

    caExit("overlapper claims to be finished, but no job list found in '$base/1-overlapper/ovljob.files'", undef)  if (! fileExists("$base/1-overlapper/ovljob.files"));

    #  Figure out how to build the store.

    if (! -e "$base/$asm.ovlStore.config") {
        $cmd  = "$bin/ovStoreConfig \\\n";
        $cmd .= " -G ./$asm.gkpStore \\\n";
        $cmd .= " -M " . getGlobal("ovsMemory") . " \\\n";    #  User supplied memory limit, reset below
        $cmd .= " -L ./1-overlapper/ovljob.files \\\n";
        $cmd .= " -create ./$asm.ovlStore.config \\\n";
        $cmd .= " > ./$asm.ovlStore.config.txt \\\n";
        $cmd .= "2> ./$asm.ovlStore.config.err\n";

        if (runCommand($base, $cmd)) {
            caExit("failed to configure the overlap store", "$base/$asm.ovlStore.config.err");
        }

        unlink "$base/$asm.ovlStore.config.err";
    }

    if (! -e "$base/$asm.ovlStore.config.txt") {
        $cmd  = "$bin/ovStoreConfig -describe ./$asm.ovlStore.config \\\n";
        $cmd .= " > ./$asm.ovlStore.config.txt \\\n";
        $cmd .= "2> ./$asm.ovlStore.config.err\n";

        if (runCommand($base, $cmd)) {
            caExit("failed to dump overlap store configuration", "$base/$asm.ovlStore.config.err");
        }

        unlink "$base/$asm.ovlStore.config.err";
    }

    my $numBuckets = 0;
    my $numSlices  = 0;
    my $sortMemory = 0;

    open(F, "< $base/$asm.ovlStore.config.txt") or caExit("can't open '$base/$asm.ovlStore.config.txt' for reading: $!\n", undef);
    while (<F>) {
        $numBuckets = $1  if (m/numBuckets\s+(\d+)/);
        $numSlices  = $1  if (m/numSlices\s+(\d+)/);
        $sortMemory = $1  if (m/sortMemory\s+(\d+)\s+GB/);
    }
    close(F);

    printf STDERR "--\n";
    printf STDERR "-- Creating overlap store $base/$asm.ovlStore using:\n";
    printf STDERR "--   %4d bucket%s\n", $numBuckets, ($numBuckets == 1) ? "" : "s";
    printf STDERR "--   %4d slice%s\n",  $numSlices,  ($numSlices  == 1) ? "" : "s";
    printf STDERR "--        using at most %d GB memory each\n", $sortMemory;

    setGlobal("ovsMemory", $sortMemory);  #  Actual memory usage of sort jobs (rounded up).

    #  If only one slice, do it all in core, otherwise, use the big gun and run it in parallel.

    if ($numSlices == 1) {
        createOverlapStoreSequential($base, $asm, $tag);
    } else {
        createOverlapStoreParallel ($base, $asm, $tag, $numBuckets, $numSlices, $sortMemory);
        overlapStoreBucketizerCheck($base, $asm, $tag, $numBuckets, $numSlices)   foreach (1..getGlobal("canuIterationMax") + 1);
        overlapStoreSorterCheck    ($base, $asm, $tag, $numBuckets, $numSlices)   foreach (1..getGlobal("canuIterationMax") + 1);

        $cmd  = "$bin/ovStoreIndexer \\\n";
        $cmd .= "  -O ./$asm.ovlStore.BUILDING \\\n";
        $cmd .= "  -G ./$asm.gkpStore \\\n";
        $cmd .= "  -C ./$asm.ovlStore.config \\\n";
        #$cmd .= "  -delete \\\n";
        $cmd .= "> ./$asm.ovlStore.BUILDING.index.err 2>&1";

        if (runCommand("$base", $cmd)) {
            caExit("failed to build index for overlap store", "$base/$asm.ovlStore.BUILDING.index.err");
        }

        unlink "$base/$asm.ovlStore.BUILDING.index.err";
        rename "$base/$asm.ovlStore.BUILDING", "$base/$asm.ovlStore";
    }

    checkOverlapStore($base, $asm);

    goto finishStage  if (getGlobal("saveOverlaps") eq "1");

    deleteOverlapIntermediateFiles($base, $asm);

    #  Now all done!

  finishStage:
    if ($tag eq "utg") {
        $cmd  = "$bin/ovStoreStats \\\n";
        $cmd .= " -G ./$asm.gkpStore \\\n";
        $cmd .= " -O ./$asm.ovlStore \\\n";
        $cmd .= " -C " . getExpectedCoverage("utg", $asm). " \\\n";
        $cmd .= " -o ./$asm.ovlStore \\\n";
        $cmd .= " > ./$asm.ovlStore.summary.err 2>&1";

        if (runCommand($base, $cmd)) {
            #  Do nothing if it fails.  This is just logging.
        }

        unlink "$base/$asm.ovlStore.summary.err";

        if (! -e "$base/$asm.ovlStore.summary") {
            print STDERR "--\n";
            print STDERR "-- WARNING: failed to generate statistics for the overlap store; no summary will appear in report.\n";
            print STDERR "--\n";
            print STDERR "----------------------------------------\n";
        } else {
            print STDERR "--\n";
            print STDERR "-- Overlap store '$base/$asm.ovlStore' contains:\n";
            print STDERR "--\n";

            my $report;

            open(F, "< $base/$asm.ovlStore.summary") or caExit("Failed to open overlap store statistics in '$base/$asm.ovlStore': $!", undef);
            while (<F>) {
                $report .= "--   $_";
            }
            close(F);

            addToReport("overlaps", $report);   #  Also shows it.
        }
    }

    emitStage($asm, "$tag-createOverlapStore");

  allDone:
    stopAfter("overlapStore");
}
