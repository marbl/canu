
###############################################################################
 #
 #  This file is part of canu, a software program that assembles whole-genome
 #  sequencing reads into contigs.
 #
 #  This software is based on:
 #    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 #    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 #
 #  Except as indicated otherwise, this is a 'United States Government Work',
 #  and is released in the public domain.
 #
 #  File 'README.licenses' in the root directory of this distribution
 #  contains full conditions and disclaimers.
 ##

package canu::OverlapStore;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(createOverlapStore);

use strict;
use warnings "all";
no  warnings "uninitialized";

use File::Basename;   #  dirname
use File::Path 2.08 qw(make_path remove_tree);
use POSIX qw(ceil);

use canu::Defaults;
use canu::Execution;

use canu::SequenceStore;
use canu::Report;

use canu::Grid_Cloud;


#  Parallel documentation: Each overlap job is converted into a single bucket of overlaps.  Within
#  each bucket, the overlaps are distributed into many slices, one per sort job.  The sort jobs then
#  load the same slice from each bucket.


#  NOT FILTERING overlaps by error rate when building the parallel store.


#
#  Fetch overlap outputs.  If you're not cloud-based, this does nothing.  Really.  Trust me.
#
sub fetchOverlapData ($$@) {
    my $asm    = shift @_;   #  Not actually used.
    my $base   = shift @_;   #  What stage we're in, 'correction', 'trimming', etc.
    my $getD   = shift @_;   #  If set, only fetch count data.

    fetchFile("$base/1-overlapper/ovljob.files");

    open(F, "< $base/1-overlapper/ovljob.files") or caExit("failed to open overlapper output list in '$base/1-overlapper/ovljob.files'", undef);
    while (<F>) {
        my $dir;
        my $num;

        if (m/^(.*)\/([0-9]*).ovb$/) {
            $dir = $1;
            $num = $2;
        } else {
            caExit("didn't recognize ovljob.files line '$_'", undef);
        }

        make_path("$base/$dir");                                  #  Make the output directory.
        fetchFile("$base/$dir/$num.oc");                          #  Fetch The Counts.   https://www.youtube.com/watch?v=vC0uvUuXVh8
        fetchFile("$base/$dir/$num.ovb")   if (defined($getD));   #  Fetch the overlap data if told to.
    }
    close(F);
}



sub fetchOverlapDataShellCode ($$) {
    my $asm    = shift @_;   #  Not actually used.
    my $base   = shift @_;   #  What stage we're in, 'correction', 'trimming', etc.
    my $string;

    fetchFile("$base/1-overlapper/ovljob.files");

    open(I, "< $base/1-overlapper/ovljob.files") or caExit("failed to open overlapper output list in '$base/1-overlapper/ovljob.files'", undef);
    while (<I>) {
        my $dir;
        my $num;

        if (m/^(.*)\/([0-9]*).ovb$/) {
            $dir = $1;
            $num = $2;
        } else {
            caExit("didn't recognize ovljob.files line '$_'", undef);
        }

        $string .= fetchFileShellCode("$base", "$dir/$num.oc",  "");
        $string .= fetchFileShellCode("$base", "$dir/$num.ovb", "");
    }
    close(I);

    return($string);
}



sub createOverlapStoreSequential ($$$) {
    my $base       = shift @_;
    my $asm        = shift @_;
    my $tag        = shift @_;
    my $bin        = getBinDirectory();
    my $cmd;

    #  Fetch all the data.

    #fetchOverlapData($asm, $base, 1);

    #  This is running in the canu process itself.  Execution.pm has special case code
    #  to submit canu to grids using the maximum of 4gb and this memory limit.

    if (! -e "./$base/$asm.ovlStore.sh") {
        open(F, "> ./$base/$asm.ovlStore.sh") or caExit("can't open './$base/$asm.ovlStore.sh' for writing: $!\n", undef);

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
        print F fetchSeqStoreShellCode($asm, $base, "");
        print F "\n";
        print F fetchFileShellCode($base, "$asm.ovlStore.config", "");
        print F "\n";
        print F fetchOverlapDataShellCode($asm, $base);
        print F "\n";
        print F "$bin/ovStoreBuild \\\n";
        print F " -O  ./$asm.ovlStore.BUILDING \\\n";
        print F" -S ../$asm.seqStore \\\n";
        print F " -C  ./$asm.ovlStore.config \\\n";
        print F " > ./$asm.ovlStore.err 2>&1 \\\n";
        print F "&& \\\n";
        print F "mv ./$asm.ovlStore.BUILDING ./$asm.ovlStore\n";
        print F "\n";
        print F stashOvlStoreShellCode($asm, $base);
        print F "\n";
        close(F);
    }

    makeExecutable("./$base/$asm.ovlStore.sh");
    stashFile("./$base/$asm.ovlStore.sh");
}



sub overlapStoreCheck ($$$) {
    my $base       = shift @_;
    my $asm        = shift @_;
    my $tag        = shift @_;
    my $attempt    = getGlobal("canuIteration");

    goto allDone   if ((-d "$base/$asm.ovlStore") || (fileExists("$base/$asm.ovlStore.tar.gz")));

    #  Since there is only one job, if we get here, we're not done.  Any other 'check' function
    #  shows how to process multiple jobs.  This only checks for the existence of the final outputs.
    #  (meryl, unitig, overlapStoreSequential are the same)

    #  If too many attempts, give up.

    if ($attempt >= getGlobal("canuIterationMax")) {
        print STDERR "--\n";
        print STDERR "-- Overlap store failed, tried $attempt times, giving up.\n";
        print STDERR "--\n";
        caExit(undef, "$base/$asm.ovlStore.err")   if (-e "$base/$asm.ovlStore.err");
        caExit(undef, undef);
    }

    if ($attempt > 0) {
        print STDERR "--\n";
        print STDERR "-- Overlap store failed, retry.\n";
        print STDERR "--\n";
    }

    #  Otherwise, run some jobs.

    generateReport($asm);

    submitOrRunParallelJob($asm, "ovS", $base, "$asm.ovlStore", (1));
    return;

  finishStage:
    print STDERR "-- Overlap store finished.\n";

    generateReport($asm);
    resetIteration("$tag-overlapStoreCheck");

  allDone:
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

    goto allDone   if ((fileExists("$path/scripts/1-bucketize.sh")) &&
                       (fileExists("$path/scripts/2-sort.sh")));
    goto allDone   if ((-d "$base/$asm.ovlStore") || (fileExists("$base/$asm.ovlStore.tar.gz")));

    make_path("$path/scripts");
    make_path("$path/logs");

    #  Parallel jobs for bucketizing.

    fetchFile("$path/scripts/1-bucketize.sh");
    fetchFile("$path/scripts/2-sort.sh");

    #  The -f option to bucketizer forces it to overwrite a partial result found
    #  in any 'create####' directory.
    #
    #  Both the binary and the script will stop if the output directory exists.

    if (! fileExists("$path/scripts/1-bucketize.sh")) {
        open(F, "> $path/scripts/1-bucketize.sh") or caExit("can't open '$path/scripts/1-bucketize.sh' for writing: $!\n", undef);
        print F "#!" . getGlobal("shell") . "\n";
        print F "\n";
        print F getBinDirectoryShellCode();
        print F "\n";
        print F setWorkDirectoryShellCode("$base/$asm.ovlStore.BUILDING");
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
        print F "jobname=`printf %04d \$jobid`\n";
        print F "\n";
        print F "if [ -e ./$asm.ovlStore.BUILDING/bucket\$jobname ] ; then\n";
        print F "  echo \"Bucketizing job finished; directory './$asm.ovlStore.BUILDING/bucket\$jobname' exists.\"\n";
        print F "  exit\n";
        print F "fi\n";
        print F "\n";
        print F fetchSeqStoreShellCode($asm, $base, "");
        print F "\n";
        print F fetchFileShellCode($base, "$asm.ovlStore.config", "");
        print F "\n";

        if (defined(getGlobal("objectStore"))) {
            print F "#\n";
            print F "#  Fetch all the inputs.\n";
            print F "#\n";
            print F "\n";
            print F "inputs=`\$bin/ovStoreConfig -describe ./$asm.ovlStore.config -listinputs \$jobid`\n";
            print F "\n";
            print F "for file in \$inputs ; do\n";                 #  The oc file doesn't appear
            print F fetchFileShellCode($base, "\$file", "  ");     #  to be needed for bucketizing.
            print F "done\n";
            print F "\n";
        }

        print F "#\n";
        print F "#  Bucketize!\n";
        print F "#\n";
        print F "\n";
        print F "\$bin/ovStoreBucketizer \\\n";
        print F "  -delete \\\n"    if ((getGlobal("purgeOverlaps") eq "aggressive") || (getGlobal("purgeOverlaps") eq "dangerous"));
        print F "  -O  ./$asm.ovlStore.BUILDING \\\n";
        print F "  -S ../$asm.seqStore \\\n";
        print F "  -C  ./$asm.ovlStore.config \\\n";
        print F "  -f \\\n";
        print F "  -b \$jobid \n";
        print F "\n";

        if (defined(getGlobal("objectStore"))) {
            print F "#\n";
            print F "#  Stash the outputs.\n";
            print F "#\n";
            print F "\n";
            print F "cd ./$asm.ovlStore.BUILDING/bucket\$jobname\n";
            print F "\n";
            print F stashFilesShellCode("$base/$asm.ovlStore.BUILDING/bucket\$jobname", "slice????", "", "purge");
            print F "\n";
            print F stashFileShellCode("$base/$asm.ovlStore.BUILDING/bucket\$jobname", "sliceSizes", "");
            print F "\n";
            print F "cd -\n";
            print F "\n";
        }

        close(F);
    }

    #  The -f option forces sorting to run even if the job appears to be still running.
    #
    #  Just after starting, ####.started is created.  This is removed on success.
    #  When writing, ####<001> exists.
    #  When finished, ####.started doesn't exist, and ####.info does.

    if (! fileExists("$path/scripts/2-sort.sh")) {
        open(F, "> $path/scripts/2-sort.sh") or caExit("can't open '$path/scripts/2-sort.sh' for writing: $!\n", undef);
        print F "#!" . getGlobal("shell") . "\n";
        print F "\n";
        print F getBinDirectoryShellCode();
        print F "\n";
        print F setWorkDirectoryShellCode("$base/$asm.ovlStore.BUILDING");
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
        print F "jobname=`printf %04d \$jobid`\n";
        print F "\n";
        print F "if [ -e ./$asm.ovlStore.BUILDING/\$jobname.info -a ! -e ./$asm.ovlStore.BUILDING/\$jobname.started ] ; then\n";
        print F "  echo \"Sorting job finished; info file './$asm.ovlStore.BUILDING/\$jobname.info' exists.\"\n";
        print F "  exit\n";
        print F "fi\n";

        if (defined(getGlobal("objectStore"))) {
            print F "#\n";
            print F "#  Fetch all the input slices, and each sliceSizes file.\n";
            print F "#\n";
            print F "\n";
            print F fetchSeqStoreShellCode($asm, $base, "");
            print F "\n";
            print F fetchFileShellCode($base, "$asm.ovlStore.config", "");
            print F "\n";
            print F "inputs=`\$bin/ovStoreConfig -describe ./$asm.ovlStore.config -listslices \$jobid`\n";
            print F "\n";
            print F "for file in \$inputs ; do\n";                                        #  $inputs is relative to
            print F fetchFileShellCode($base, "./$asm.ovlStore.BUILDING/\$file", "  ");   #  ovStore.BUILDING!
            print F "done\n";
            print F "\n";
        }

        print F "#\n";
        print F "#  Sort!\n";
        print F "#\n";
        print F "\n";
        print F "\$bin/ovStoreSorter \\\n";
        print F "  -deleteearly \\\n"    if ((getGlobal("purgeOverlaps") eq "dangerous"));
        print F "  -deletelate \\\n"     if ((getGlobal("purgeOverlaps") eq "normal") || (getGlobal("purgeOverlaps") eq "aggressive"));
        print F "  -O  ./$asm.ovlStore.BUILDING \\\n";
        print F "  -S ../$asm.seqStore \\\n";
        print F "  -C  ./$asm.ovlStore.config \\\n";
        print F "  -f \\\n";
        print F "  -s \$jobid \\\n";
        print F "  -M $sortMemory \n";
        print F "\n";

        if (defined(getGlobal("objectStore"))) {
            print F "#\n";
            print F "#  Stash the outputs.\n";
            print F "#\n";
            print F "\n";
            print F "jobid=`printf %04d \$jobid`\n";
            print F "\n";
            print F "cd ./$asm.ovlStore.BUILDING\n";
            print F "\n";
            print F stashFilesShellCode("$base/$asm.ovlStore.BUILDING", "\$jobid\\<???\\>", "", "do-not-purge");
            print F "\n";
            print F "tar -cf - \\\n";
            print F "  ./\$jobid.info \\\n";
            print F "  ./\$jobid.index \\\n";
            print F "  ./\$jobid*.statistics \\\n";
            print F "| \\\n";
            print F "gzip -1vc > ./\$jobid.statistics.tar.gz\n";
            print F "\n";
            print F stashFileShellCode("$base/$asm.ovlStore.BUILDING", "\$jobid.statistics.tar.gz", "");
            print F "\n";
            print F "cd ..\n";
            print F "\n";
        }

        close(F);
    }

    #  The final job to merge all the indices is done in createOverlapStore() below.

    makeExecutable("$path/scripts/1-bucketize.sh");
    makeExecutable("$path/scripts/2-sort.sh");

    stashFile("$path/scripts/1-bucketize.sh");
    stashFile("$path/scripts/2-sort.sh");

  finishStage:
    generateReport($asm);
    resetIteration("$tag-overlapStoreConfigure");

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

    goto allDone   if (fileExists("$path/1-bucketize.success"));
    goto allDone   if ((-d "$base/$asm.ovlStore") || (fileExists("$base/$asm.ovlStore.tar.gz")));

    make_path("$path/scripts");
    make_path("$path/logs");

    fetchFile("$path/scripts/1-bucketize.sh");

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
        if (! fileExists("$path/bucket$bucketID/sliceSizes")) {
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

        generateReport($asm);

        submitOrRunParallelJob($asm, "ovB", $path, "scripts/1-bucketize", @failedJobs);
        return;
    }

  finishStage:
    print STDERR "-- Overlap store bucketizer finished.\n";

    #  Note that we're done.

    touch("$path/1-bucketize.success");
    stashFile("$path/1-bucketize.success");

    #  All overlap inputs are bucketized, so delete them.

    if ((getGlobal("purgeOverlaps") eq "aggressive") ||
        (getGlobal("purgeOverlaps") eq "dangerous")) {
        deleteOverlapIntermediateFiles($base, $asm);
    }

    #  And show a report.

    generateReport($asm);
    resetIteration("$tag-overlapStoreBucketizerCheck");

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

    goto allDone   if (fileExists("$path/2-sorter.success"));
    goto allDone   if ((-d "$base/$asm.ovlStore") || (fileExists("$base/$asm.ovlStore.tar.gz")));

    make_path("$path/scripts");
    make_path("$path/logs");

    fetchFile("$path/scripts/2-sort.sh");

    #  Figure out if all the tasks finished correctly.

    my $currentJobID   = 1;
    my $sliceID        = "0001";

    my @successJobs;
    my @failedJobs;
    my $failureMessage = "";

    #  A sort job is done if either:
    #    the info exists            (for normal storage)
    #    a statistics.tar.gz exists (for object storage)

    for (my $ss=1; $ss<=$numSlices; $ss++) {
        if ((! fileExists("$path/$sliceID.info")) &&
            (! fileExists("$path/$sliceID.statistics.tar.gz"))) {
            $failureMessage .= "--   job $path/$sliceID FAILED.\n";
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

        generateReport($asm);

        submitOrRunParallelJob($asm, "ovS", $path, "scripts/2-sort", @failedJobs);
        return;
    }

  finishStage:
    print STDERR "-- Overlap store sorter finished.\n";

    touch("$path/2-sorter.success");
    stashFile("$path/2-sorter.success");

    generateReport($asm);
    resetIteration("$tag-overlapStoreSorterCheck");

  allDone:
}



sub overlapStoreIndexerCheck ($$$$$) {
    my $base       = shift @_;
    my $asm        = shift @_;
    my $tag        = shift @_;
    my $numBuckets = shift @_;
    my $numSlices  = shift @_;

    my $bin   = getBinDirectory();
    my $cmd;

    #  If running with showNext, ovStoreIndexer will be run (manually
    #  possibly) leaving an index file behind, but not renaming the store.
    #  In that case, jump straight to finishStage and rename it.

    goto allDone      if ((-d "$base/$asm.ovlStore") || (fileExists("$base/$asm.ovlStore.tar.gz")));

    goto finishStage  if (-e  "$base/$asm.ovlStore.BUILDING/index");

    #  Fetch the stats and index data.  If not using an object store, the
    #  fetch does nothing, and since there is no file, the gzip/tar are
    #  skipped.

    for (my $ss=1; $ss<=$numSlices; $ss++) {
        my $slice = substr("0000" . $ss, -4);

        fetchFile("$base/$asm.ovlStore.BUILDING/$slice.statistics.tar.gz");

        if (-e "$base/$asm.ovlStore.BUILDING/$slice.statistics.tar.gz") {
            runCommandSilently("$base/$asm.ovlStore.BUILDING", "gzip -dc $slice.statistics.tar.gz | tar -xf -", 1);
            unlink("$base/$asm.ovlStore.BUILDING/$slice.statistics.tar.gz");
        }
    }

    $cmd  = "$bin/ovStoreIndexer \\\n";
    $cmd .= "  -O  ./$asm.ovlStore.BUILDING \\\n";
    $cmd .= "  -S ../$asm.seqStore \\\n";
    $cmd .= "  -C  ./$asm.ovlStore.config \\\n";
    $cmd .= "  -delete \\\n";
    $cmd .= "> ./$asm.ovlStore.BUILDING.index.err 2>&1";

    if (runCommand("$base", $cmd)) {
        caExit("failed to build index for overlap store", "$base/$asm.ovlStore.BUILDING.index.err");
    }

    unlink "$base/$asm.ovlStore.BUILDING.index.err";

  finishStage:
    print STDERR "-- Overlap store indexer finished.\n";
 
    rename "$base/$asm.ovlStore.BUILDING", "$base/$asm.ovlStore";

    renameStashedFile("$base/$asm.ovlStore.BUILDING", "$base/$asm.ovlStore");
    stashOvlStore($asm, $base);

    #resetIteration("$tag-overlapStoreIndexerCheck");

  allDone:
}



sub checkOverlapStore ($$) {
    my $base    = shift @_;
    my $asm     = shift @_;

    my $bin   = getBinDirectory();
    my $cmd;

    $cmd  = "$bin/ovStoreDump \\\n";
    $cmd .= " -S ../$asm.seqStore \\\n";
    $cmd .= " -O  ./$asm.ovlStore \\\n";
    $cmd .= " -counts \\\n";
    $cmd .= " > ./$asm.ovlStore/counts.dat 2> ./$asm.ovlStore/counts.err";

    print STDERR "-- Checking store.\n";

    if (runCommand($base, $cmd)) {
        caExit("failed to dump counts of overlaps; invalid store?", "$base/$asm.ovlStore/counts.err");
    }

    my $totOvl   = 0;
    my $nulReads = 0;
    my $ovlReads = 0;

    open(F, "< ./$base/$asm.ovlStore/counts.dat") or caExit("can't open './$base/$asm.ovlStore/counts.dat' for reading: $!\n", undef);
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
sub deleteOverlapIntermediateFiles ($$@) {
    my $base    = shift @_;
    my $asm     = shift @_;
    my @what    = @_;        #  types of data to delete
    my @files;               #  list of files to read

    #  If a list of stuff to delete is supplied, delete only that stuff.  Otherwise,
    #  delete everything.

    if (scalar(@what) == 0) {
        #print STDERR "SET what to default list\n";
        @what = ( "precompute", "mhap", "mmap", "ovljob", "ovljob.more" );
    }

    foreach my $file (@what) {
        #print STDERR "ADD $base/1-overlapper/$file.files ?\n";
        push @files, "$base/1-overlapper/$file.files"   if (-e "$base/1-overlapper/$file.files");
    }

    #  Now open each file and delete what ever files are listed in there.

    my %directories;
    my $bytes = 0;
    my $files = 0;

    foreach my $file (@files) {
        #print STDERR "UNLINK from $file\n";

        open(F, "< $file") or caExit("can't open '$file' for reading: $!\n", undef);
        while (<F>) {
            chomp;

            if (-d "$base/$_") {
                $directories{$_}++;

            } elsif (-e "$base/$_") {
                $bytes += -s "$base/$_";
                $files += 1;

                #print STDERR "UNLINK '$base/$_'\n";
                unlink "$base/$_";
                rmdir dirname("$base/$_");  #  Try to rmdir the directory the file is in.  If empty, yay!
            }
        }
        close(F);
    }

    foreach my $dir (keys %directories) {
        #print STDERR "RMTREE '$base/$dir'\n";
        remove_tree("$base/$dir");
    }

#    foreach my $file (@files) {
#        #print STDERR "UNLINK '$file'  SOURCE\n";
#        unlink $file;
#    }

    print STDERR "--\n";
    print STDERR "-- Purged ", int(1000 * $bytes / 1024 / 1024 / 1024) / 1000, " GB in $files overlap output files.\n";
}



#
#  Scan the overlap store for old-format file names.  If any are found,
#  create new-format names as symlinks to the original files.
#
sub updateOverlapStoreFiles ($$$) {
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $base    = shift @_;

    if (-e "$base/$asm.ovlStore") {
        my %oldToNew;
        my @newExists;

        #  Scan the store for any old-format files.
        #  Remember if new-format links exist too.
        open(L, "ls $base/$asm.ovlStore/ | ");
        while (<L>) {
            chomp;

            if (m/^(\d\d\d\d)<(\d\d\d)>$/) {
                my $old = "$1<$2>";
                my $new = "$1-$2";

                $oldToNew{$old} = $new;
            }

            if (m/^(\d\d\d\d)-(\d\d\d)$/) {
                my $old = "$1<$2>";
                my $new = "$1-$2";

                push @newExists, $old;
            }
        }
        close(L);

        #  Forget about old-format names with new-format links already present.
        foreach my $nn (@newExists) {
            delete $oldToNew{$nn};
        }

        #  If there are any old-format names left, make links for them, logging
        #  whatever we do.
        if (scalar(%oldToNew) > 0) {
            print STDERR "--\n";
            print STDERR "-- UPDATING $base/$asm.ovlStore data files to new name format:\n";

            foreach my $oo (keys %oldToNew) {
                print STDERR "--   '$oo' -> '$oldToNew{$oo}'\n";
                system("ln -s '$oo' '$base/$asm.ovlStore/$oldToNew{$oo}'");
            }

            print STDERR "--\n";
        }
    }
}



sub createOverlapStore ($$) {
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $cmd;
    my $bin     = getBinDirectory();

    my $base;

    $base = "correction"  if ($tag eq "cor");
    $base = "trimming"    if ($tag eq "obt");
    $base = "unitigging"  if ($tag eq "utg");

    #  A compatibility fix for renaming ovlStore data files from 0000<000> to 0000-000.
    #  This change occurred around July 7 2020.

    updateOverlapStoreFiles($asm, $tag, $base);

    #  Now, if the directory exists, we have a store.

    goto allDone   if ((-d "$base/$asm.ovlStore") || (fileExists("$base/$asm.ovlStore.tar.gz")));
    goto allDone   if ((-d "$base/$asm.ctgStore") || (fileExists("$base/$asm.ctgStore.tar.gz")));

    #  Did we _really_ complete?

    caExit("overlapper claims to be finished, but no job list found in '$base/1-overlapper/ovljob.files'", undef)  if (! fileExists("$base/1-overlapper/ovljob.files"));

    #  Apparently so.  Delete the mhap precompute data if we're being aggressive.

    if ((getGlobal("purgeOverlaps") eq "aggressive") ||
        (getGlobal("purgeOverlaps") eq "dangerous")) {
        deleteOverlapIntermediateFiles($base, $asm, "precompute");
    }

    #  Figure out how to build the store.

    fetchFile("$base/1-overlapper/ovljob.files");
    fetchFile("$base/$asm.ovlStore.config");
    fetchFile("$base/$asm.ovlStore.config.txt");

    if (! -e "$base/$asm.ovlStore.config") {
        $cmd  = "$bin/ovStoreConfig \\\n";
        $cmd .= " -S ../$asm.seqStore \\\n";
        $cmd .= " -M " . getGlobal("ovsMemory") . " \\\n";    #  User supplied memory limit, reset below
        $cmd .= " -L ./1-overlapper/ovljob.files \\\n";
        $cmd .= " -create ./$asm.ovlStore.config \\\n";
        $cmd .= " > ./$asm.ovlStore.config.txt \\\n";
        $cmd .= "2> ./$asm.ovlStore.config.err";

        fetchOverlapData($asm, $base);   #  Fetch overlap data, if not here.

        if (runCommand($base, $cmd)) {
            caExit("failed to configure the overlap store", "$base/$asm.ovlStore.config.err");
        }

        unlink "$base/$asm.ovlStore.config.err";

        stashFile("$base/$asm.ovlStore.config");
        stashFile("$base/$asm.ovlStore.config.txt");
    }

    if (! -e "$base/$asm.ovlStore.config.txt") {
        fetchFile("$base/$asm.ovlStore.config");

        $cmd  = "$bin/ovStoreConfig -describe ./$asm.ovlStore.config \\\n";
        $cmd .= " > ./$asm.ovlStore.config.txt \\\n";
        $cmd .= "2> ./$asm.ovlStore.config.err";

        if (runCommand($base, $cmd)) {
            caExit("failed to dump overlap store configuration", "$base/$asm.ovlStore.config.err");
        }

        unlink "$base/$asm.ovlStore.config.err";

        stashFile("$base/$asm.ovlStore.config.txt");
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

    setGlobal("ovsMemory", $sortMemory + 2);  #  Actual memory usage of sort jobs (rounded up).

    #  If only one slice, do it all in core, otherwise, use the big gun and run it in parallel.

    if ($numSlices == 1) {
        createOverlapStoreSequential($base, $asm, $tag);
        overlapStoreCheck           ($base, $asm, $tag)   foreach (1..getGlobal("canuIterationMax") + 1);
    }

    else {
        createOverlapStoreParallel ($base, $asm, $tag, $numBuckets, $numSlices, $sortMemory);
        overlapStoreBucketizerCheck($base, $asm, $tag, $numBuckets, $numSlices)   foreach (1..getGlobal("canuIterationMax") + 1);
        overlapStoreSorterCheck    ($base, $asm, $tag, $numBuckets, $numSlices)   foreach (1..getGlobal("canuIterationMax") + 1);
        overlapStoreIndexerCheck   ($base, $asm, $tag, $numBuckets, $numSlices);
    }

  finishStage:
    checkOverlapStore($base, $asm);

    #  The store is fully created and tested, so delete the overlapper outputs.
    if ((getGlobal("purgeOverlaps") eq "normal") ||
        (getGlobal("purgeOverlaps") eq "aggressive") ||
        (getGlobal("purgeOverlaps") eq "dangerous")) {
        deleteOverlapIntermediateFiles($base, $asm);
    }

    if ($tag eq "utg") {
        $cmd  = "$bin/ovStoreStats \\\n";
        $cmd .= " -C " . getExpectedCoverage($asm, "utg"). " \\\n";
        $cmd .= " -S ../$asm.seqStore \\\n";
        $cmd .= " -O  ./$asm.ovlStore \\\n";
        $cmd .= " -o  ./$asm.ovlStore \\\n";
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

    generateReport($asm);
    resetIteration("$tag-createOverlapStore");

  allDone:
    stopAfter("overlapStore");
}
