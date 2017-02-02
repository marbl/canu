
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

use POSIX qw(ceil);
use canu::Defaults;
use canu::Execution;
use canu::HTML;


#  Parallel documentation: Each overlap job is converted into a single bucket of overlaps.  Within
#  each bucket, the overlaps are distributed into many slices, one per sort job.  The sort jobs then
#  load the same slice from each bucket.


#  NOT FILTERING overlaps by error rate when building the parallel store.


sub createOverlapStoreSequential ($$$$) {
    my $base    = shift @_;
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $files   = shift @_;
    my $bin     = getBinDirectory();
    my $cmd;

    #  This is running in the canu process itself.  Execution.pm has special case code
    #  to submit canu to grids using the maximum of 4gb and this memory limit.

    my $memSize = getGlobal("ovsMemory");

    #  The parallel store build will unlimit 'max user processes'.  The sequential method usually
    #  runs out of open file handles first (meaning it has never run out of processes yet).

    $cmd  = "$bin/ovStoreBuild \\\n";
    $cmd .= " -O ./$asm.ovlStore.BUILDING \\\n";
    $cmd .= " -G ./$asm.gkpStore \\\n";
    $cmd .= " -M $memSize \\\n";
    $cmd .= " -L ./1-overlapper/ovljob.files \\\n";
    $cmd .= " > ./$asm.ovlStore.err 2>&1";

    if (runCommand($base, $cmd)) {
        caExit("failed to create the overlap store", "$base/$asm.ovlStore.err");
    }

    unlink "$base/$asm.ovlStore.err";

    rename "$base/$asm.ovlStore.BUILDING", "$base/$asm.ovlStore";
}




#  Count the number of inputs.  We don't expect any to be missing (they were just checked
#  by overlapCheck()) but feel silly not checking again.
#
#  Note that these files are rooted in '$base' (because that's where we run the overlap store
#  building) but canu.pl itself is rooted in the same directory as '$base', so we need to
#  add in '$base'.

sub countOverlapStoreInputs ($) {
    my $base      = shift @_;
    my $numInputs = 0;

    open(F, "< $base/1-overlapper/ovljob.files") or caExit("Failed to open overlap store input file '$base/1-overlapper/ovljob.files': $0", undef);
    while (<F>) {
        chomp;
        caExit("overlapper output '$_' not found", undef)   if (! -e "$base/$_");
        $numInputs++;
    }
    close(F);

    return($numInputs);
}




sub getNumOlapsAndSlices ($$) {
    my $base = shift @_;
    my $asm  = shift @_;

    my $numOlaps   = undef;
    my $numSlices  = undef;
    my $memLimit   = undef;

    open(F, "< $base/$asm.ovlStore.BUILDING/config.err") or caExit("can't open '$base/$asm.ovlStore.BUILDING/config.err' for reading: $!\n", undef);
    while (<F>) {
        if (m/Will sort (\d+.\d+) million overlaps per bucket, using (\d+) buckets (\d+.\d+) GB per bucket./) {
            $numOlaps  = $1;
            $numSlices = $2;
            $memLimit  = ceil($3);
        }
    }
    close(F);

    if (!defined($numOlaps) || !defined($numSlices) || !defined($memLimit)) {
        caExit("Failed to find any overlaps ($numOlaps) or slices ($numSlices) or memory limit ($memLimit)", undef);
    }

    return($numOlaps, $numSlices, $memLimit);
}



sub overlapStoreConfigure ($$$$) {
    my $base    = shift @_;
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $files   = shift @_;
    my $bin     = getBinDirectory();
    my $cmd;

    goto allDone   if (skipStage($asm, "$tag-overlapStoreConfigure") == 1);
    goto allDone   if (-d "$base/$asm.ovlStore");

    my $numInputs  = countOverlapStoreInputs($base);

    #  Create an output directory, and populate it with more directories and scripts

    system("mkdir -p $base/$asm.ovlStore.BUILDING")                   if (! -d "$base/$asm.ovlStore.BUILDING");
    system("mkdir -p $base/$asm.ovlStore.BUILDING/scripts")           if (! -d "$base/$asm.ovlStore.BUILDING/scripts");
    system("mkdir -p $base/$asm.ovlStore.BUILDING/logs")              if (! -d "$base/$asm.ovlStore.BUILDING/logs");

    #  Run the normal store build, but just to get the partitioning.  ovStoreBuild internally
    #  writes to config.WORKING, then renames when it is finished.  No need for the script
    #  to be overly careful about incomplete files.

    if (! -e "$base/$asm.ovlStore.BUILDING/scripts/0-config.sh") {
        open(F, "> $base/$asm.ovlStore.BUILDING/scripts/0-config.sh") or die;
        print F "#!" . getGlobal("shell") . "\n";
        print F "\n";
        print F getLimitShellCode("processes");
        print F getLimitShellCode("files");
        print F "\n";
        print F getBinDirectoryShellCode();
        print F "\n";
        print F "\$bin/ovStoreBuild \\\n";
        print F " -G ./$asm.gkpStore \\\n";
        print F " -O ./$asm.ovlStore \\\n";  #  NOT created!
        print F " -M " . getGlobal("ovsMemory") . " \\\n";
        print F " -config ./$asm.ovlStore.BUILDING/config \\\n";
        print F " -L ./1-overlapper/ovljob.files \\\n";
        close(F);
    }
    system("chmod +x $base/$asm.ovlStore.BUILDING/scripts/0-config.sh");

    if (! -e "$base/$asm.ovlStore.BUILDING/config") {
        $cmd  = "./$asm.ovlStore.BUILDING/scripts/0-config.sh \\\n";
        $cmd .= "> ./$asm.ovlStore.BUILDING/config.err 2>&1\n";

        if (runCommand($base, $cmd)) {
            caExit("failed to generate configuration for building overlap store", "$base/$asm.ovlStore.BUILDING/config.err");
        }
    }

    #  Parse the output to find the number of jobs we need to sort and the memory
    #  ovs store memory is left as a range (e.g. 4-16) so building can scale itself to (hopefully) fit both into memory and into max system open files
    my ($numOlaps, $numSlices, $memLimit) = getNumOlapsAndSlices($base, $asm);
    setGlobal("ovsMemory", ceil($memLimit * 1.1)); # request more memory to avoid memory issues on machines with no swap

    #  Parallel jobs for bucketizing.  This should really be part of overlap computation itself.

    #getAllowedResources("", "ovb");

    if (! -e "$base/$asm.ovlStore.BUILDING/scripts/1-bucketize.sh") {
        open(F, "> $base/$asm.ovlStore.BUILDING/scripts/1-bucketize.sh") or die;
        print F "#!" . getGlobal("shell") . "\n";
        print F "\n";
        print F getJobIDShellCode();
        print F "\n";
        print F "bn=`printf %04d \$jobid`\n";
        print F "jn=\"undefined\"\n";
        print F "\n";
        print F "#  This script runs in $base/$asm.ovlStore.BUILDING/, but the overlap file list\n";
        print F "#  is relative to $base/, so we need to add a few dots to make things work.\n";
        print F "\n";

        my $tstid = 1;

        open(I, "< $files") or die "Failed to open '$files': $0\n";

        while (<I>) {
            chomp;

            print F "if [ \"\$jobid\" -eq \"$tstid\" ] ; then jn=\"../$_\"; fi\n";
            $tstid++;
        }

        close(I);

        print F "\n";
        print F "if [ \$jn = \"undefined\" ] ; then\n";
        print F "  echo \"Job out of range.\"\n";
        print F "  exit\n";
        print F "fi\n";
        print F "\n";
        print F "if [ -e \"./bucket\$bn/sliceSizes\" ] ; then\n";
        print F "  echo \"Bucket bucket\$bn finished successfully.\"\n";
        print F "  exit\n";
        print F "fi\n";
        print F "\n";
        print F "if [ -e \"./create\$bn\" ] ; then\n";
        print F "  echo \"Removing incomplete bucket create\$bn\"\n";
        print F "  rm -rf \"./create\$bn\"\n";
        print F "fi\n";
        print F "\n";
        print F getLimitShellCode("processes");
        print F getLimitShellCode("files");
        print F "\n";
        print F getBinDirectoryShellCode();
        print F "\n";
        print F "\$bin/ovStoreBucketizer \\\n";
        print F "  -O . \\\n";
        print F "  -G ../$asm.gkpStore \\\n";
        print F "  -C ./config \\\n";
        #print F "  -e " . getGlobal("") . " \\\n"  if (defined(getGlobal("")));
        print F "  -job \$jobid \\\n";
        print F "  -i   \$jn\n";
        print F "\n";
        print F "if [ \$? = 0 ] ; then\n";
        print F "  echo Success.\n";
        print F "  exit 0\n";
        print F "else\n";
        print F "  echo Failure.\n";
        print F "  exit 1\n";
        print F "fi\n";
        close(F);
    }

    #  Parallel jobs for sorting each bucket

    #getAllowedResources("", "ovs");

    if (! -e "$base/$asm.ovlStore.BUILDING/scripts/2-sort.sh") {
        open(F, "> $base/$asm.ovlStore.BUILDING/scripts/2-sort.sh") or die;
        print F "#!" . getGlobal("shell") . "\n";
        print F "\n";
        print F getJobIDShellCode();
        print F "\n";
        print F getLimitShellCode("processes");
        print F getLimitShellCode("files");
        print F "\n";
        print F getBinDirectoryShellCode();
        print F "\n";
        print F "\$bin/ovStoreSorter \\\n";
        print F "  -deletelate \\\n";  #  Choices -deleteearly -deletelate or nothing
        print F "  -M $memLimit \\\n";
        print F "  -O . \\\n";
        print F "  -G ../$asm.gkpStore \\\n";
        print F "  -F $numSlices \\\n";
        print F "  -job \$jobid $numInputs\n";
        print F "\n";
        print F "if [ \$? = 0 ] ; then\n";
        print F "  echo Success.\n";
        print F "  exit 0\n";
        print F "else\n";
        print F "  echo Failure.\n";
        print F "  exit 1\n";
        print F "fi\n";
        close(F);
    }

    #  A final job to merge the indices.

    if (! -e "$base/$asm.ovlStore.BUILDING/scripts/3-index.sh") {
        open(F, "> $base/$asm.ovlStore.BUILDING/scripts/3-index.sh") or die;
        print F "#!" . getGlobal("shell") . "\n";
        print F "\n";
        print F getBinDirectoryShellCode();
        print F "\n";
        print F "\$bin/ovStoreIndexer \\\n";
        #print F "  -nodelete \\\n";  #  Choices -nodelete or nothing
        print F "  -O . \\\n";
        print F "  -F $numSlices\n";
        print F "\n";
        print F "if [ \$? = 0 ] ; then\n";
        print F "  echo Success.\n";
        print F "  exit 0\n";
        print F "else\n";
        print F "  echo Failure.\n";
        print F "  exit 1\n";
        print F "fi\n";
        close(F);
    }

    system("chmod +x $base/$asm.ovlStore.BUILDING/scripts/1-bucketize.sh");
    system("chmod +x $base/$asm.ovlStore.BUILDING/scripts/2-sort.sh");
    system("chmod +x $base/$asm.ovlStore.BUILDING/scripts/3-index.sh");

  finishStage:
    emitStage($asm, "$tag-overlapStoreConfigure");
    buildHTML($asm, $tag);

  allDone:
    stopAfter("overlapStoreConfigure");
}



sub overlapStoreBucketizerCheck ($$$$) {
    my $base    = shift @_;
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $files   = shift @_;
    my $attempt = getGlobal("canuIteration");

    goto allDone   if (skipStage($asm, "$tag-overlapStoreBucketizerCheck", $attempt) == 1);
    goto allDone   if (-d "$base/$asm.ovlStore");
    goto allDone   if (-e "$base/$asm.ovlStore.BUILDING/1-bucketize.success");

    #  Figure out if all the tasks finished correctly.

    my $numInputs      = countOverlapStoreInputs($base);
    my $currentJobID   = 1;
    my @successJobs;
    my @failedJobs;
    my $failureMessage = "";

    my $bucketID       = "0001";

    #  Two ways to check for completeness, either 'sliceSizes' exists, or the 'bucket' directory
    #  exists.  The compute is done in a 'create' directory, which is renamed to 'bucket' just
    #  before the job completes.

    open(F, "< $files") or caExit("can't open '$files' for reading: $!", undef);

    while (<F>) {
        chomp;

        if (! -e "$base/$asm.ovlStore.BUILDING/bucket$bucketID") {
            $failureMessage .= "--   job $base/$asm.ovlStore.BUILDING/bucket$bucketID FAILED.\n";
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

        #  If not the first attempt, report the jobs that failed, and that we're recomputing.

        if ($attempt > 1) {
            print STDERR "--\n";
            print STDERR "-- ", scalar(@failedJobs), " overlap store bucketizer jobs failed:\n";
            print STDERR $failureMessage;
            print STDERR "--\n";
        }

        #  If too many attempts, give up.

        if ($attempt > getGlobal("canuIterationMax")) {
            caExit("failed to overlapStoreBucketize.  Made " . ($attempt-1) . " attempts, jobs still failed", undef);
        }

        #  Otherwise, run some jobs.

        print STDERR "-- overlap store bucketizer attempt $attempt begins with ", scalar(@successJobs), " finished, and ", scalar(@failedJobs), " to compute.\n";

        emitStage($asm, "$tag-overlapStoreBucketizerCheck", $attempt);
        buildHTML($asm, $tag);

        submitOrRunParallelJob($asm, "ovB", "$base/$asm.ovlStore.BUILDING", "scripts/1-bucketize", @failedJobs);
        return;
    }

  finishStage:
    print STDERR "-- Overlap store bucketizer finished.\n";

    touch("$base/$asm.ovlStore.BUILDING/1-bucketize.success");

    setGlobal("canuIteration", 1);
    emitStage($asm, "$tag-overlapStoreBucketizerCheck");
    buildHTML($asm, $tag);

  allDone:
}





sub overlapStoreSorterCheck ($$$$) {
    my $base    = shift @_;
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $files   = shift @_;
    my $attempt = getGlobal("canuIteration");

    goto allDone   if (skipStage($asm, "$tag-overlapStoreSorterCheck", $attempt) == 1);
    goto allDone   if (-d "$base/$asm.ovlStore");
    goto allDone   if (-e "$base/$asm.ovlStore.BUILDING/2-sorter.success");

    #  Figure out if all the tasks finished correctly.

    my ($numOlaps, $numSlices) = getNumOlapsAndSlices($base, $asm);

    my $currentJobID   = 1;
    my @successJobs;
    my @failedJobs;
    my $failureMessage = "";

    my $sortID       = "0001";

    open(F, "< $files") or caExit("can't open '$files' for reading: $!", undef);

    #  A valid result has three files:
    #    $base/$asm.ovlStore.BUILDING/$sortID
    #    $base/$asm.ovlStore.BUILDING/$sortID.index
    #    $base/$asm.ovlStore.BUILDING/$sortID.info
    #
    #  A crashed result has one file, if it crashes before output
    #    $base/$asm.ovlStore.BUILDING/$sortID.ovs
    #
    #  On out of disk, the .info is missing.  It's the last thing created.
    #
    while ($currentJobID <= $numSlices) {

        if ((! -e "$base/$asm.ovlStore.BUILDING/$sortID") ||
            (! -e "$base/$asm.ovlStore.BUILDING/$sortID.info") ||
            (  -e "$base/$asm.ovlStore.BUILDING/$sortID.ovs")) {
            $failureMessage .= "--   job $base/$asm.ovlStore.BUILDING/$sortID FAILED.\n";
            unlink "$base/$asm.ovlStore.BUILDING/$sortID.ovs";
            push @failedJobs, $currentJobID;
        } else {
            push @successJobs, $currentJobID;
        }

        $currentJobID++;
        $sortID++;
    }

    close(F);

    #  Failed jobs, retry.

    if (scalar(@failedJobs) > 0) {

        #  If not the first attempt, report the jobs that failed, and that we're recomputing.

        if ($attempt > 1) {
            print STDERR "--\n";
            print STDERR "-- ", scalar(@failedJobs), " overlap store sorter jobs failed:\n";
            print STDERR $failureMessage;
            print STDERR "--\n";
        }

        #  If too many attempts, give up.

        if ($attempt > getGlobal("canuIterationMax")) {
            caExit("failed to overlapStoreSorter.  Made " . ($attempt-1) . " attempts, jobs still failed", undef);
        }

        #  Otherwise, run some jobs.

        print STDERR "-- overlap store sorter attempt $attempt begins with ", scalar(@successJobs), " finished, and ", scalar(@failedJobs), " to compute.\n";

        emitStage($asm, "$tag-overlapStoreSorterCheck", $attempt);
        buildHTML($asm, $tag);

        submitOrRunParallelJob($asm, "ovS", "$base/$asm.ovlStore.BUILDING", "scripts/2-sort", @failedJobs);
        return;
    }

  finishStage:
    print STDERR "-- Overlap store sorter finished.\n";

    touch("$base/$asm.ovlStore.BUILDING/2-sorter.success");

    setGlobal("canuIteration", 1);
    emitStage($asm, "$tag-overlapStoreSorterCheck");
    buildHTML($asm, $tag);

  allDone:
}




sub createOverlapStoreParallel ($$$$) {
    my $base    = shift @_;
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $files   = shift @_;

    overlapStoreConfigure      ($base, $asm, $tag, $files);
    overlapStoreBucketizerCheck($base, $asm, $tag, $files)   foreach (1..getGlobal("canuIterationMax") + 1);
    overlapStoreSorterCheck    ($base, $asm, $tag, $files)   foreach (1..getGlobal("canuIterationMax") + 1);

    if (runCommand("$base/$asm.ovlStore.BUILDING", "./scripts/3-index.sh > ./logs/3-index.err 2>&1")) {
        caExit("failed to build index for overlap store", "$base/$asm.ovlStore.BUILDING/logs/3-index.err");
    }

    rename "$base/$asm.ovlStore.BUILDING", "$base/$asm.ovlStore";
}


sub generateOverlapStoreStats ($$) {
    my $base    = shift @_;
    my $asm     = shift @_;

    my $bin   = getBinDirectory();
    my $cmd;

    $cmd  = "$bin/ovStoreStats \\\n";
    $cmd .= " -G ./$asm.gkpStore \\\n";
    $cmd .= " -O ./$asm.ovlStore \\\n";
    $cmd .= " -o ./$asm.ovlStore \\\n";
    $cmd .= " > ./$asm.ovlStore.summary.err 2>&1";

    if (runCommand($base, $cmd)) {
        print STDERR "--\n";
        print STDERR "-- WARNING: failed to generate statistics for the overlap store; no summary will appear in HTML output.\n";
        print STDERR "--\n";
        print STDERR "----------------------------------------\n";
    }
}


sub createOverlapStore ($$$) {
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $seq     = shift @_;

    my $base;

    $base = "correction"  if ($tag eq "cor");
    $base = "trimming"    if ($tag eq "obt");
    $base = "unitigging"  if ($tag eq "utg");

    my $path  = "$base/1-overlapper";

    goto allDone   if (skipStage($asm, "$tag-createOverlapStore") == 1);
    goto allDone   if (-d "$base/$asm.ovlStore");
    goto allDone   if (-d "$base/$asm.ctgStore");

    #  Did we _really_ complete?

    caExit("overlapper claims to be finished, but no job list found in '$path/ovljob.files'", undef)  if (! -e "$path/ovljob.files");

    #  Then just build the store!  Simple!

    createOverlapStoreSequential($base, $asm, $tag, "$path/ovljob.files")  if ($seq eq "sequential");
    createOverlapStoreParallel  ($base, $asm, $tag, "$path/ovljob.files")  if ($seq eq "parallel");

    print STDERR "--\n";
    print STDERR "-- Overlap store '$base/$asm.ovlStore' successfully constructed.\n";

    goto finishStage  if (getGlobal("saveOverlaps") eq "1");

    #  Delete the inputs and directories.  Some contortions are needed to get directory deletes in order.
    #  In particular, mhap's queries directory needs to be deleted after it's subdirectories are.

    my %directories;
    my $bytes = 0;
    my $files = 0;

    foreach my $file ("$path/ovljob.files", "$path/ovljob.more.files", "$path/mhap.files", "$path/mmap.files", "$path/precompute.files") {
        next  if (! -e $file);

        open(F, "< $file") or caExit("can't open '$file' for reading: $!\n", undef);
        while (<F>) {
            chomp;

            #  Decide what to do.  Directories - register for later deletion.  Files - sum size and
            #  delete.

            if (-d "$base/$_") {
                $directories{"$base/$_"}++;

            } elsif (-e "$base/$_") {
                $bytes += -s "$base/$_";
                $files += 1;

                unlink "$base/$_";
            }

            #  If the path isn't a directory register the directory it is in for deletion.

            $directories{dirname("$base/$_")}++   if (! -d "$base/$_");
        }
        close(F);

        unlink $file;
    }

    #  Ideally, every directory we have in our list should be empty.  But they won't be.  So, loop until we fail to delete anything.

    my $dirs = 0;
    my $deleted = 1;

    while ($deleted > 0) {
        $deleted = 0;

        foreach my $dir (keys %directories) {
            if (-d $dir) {
                rmdir $dir;
                $dirs++;
            }

            if (! -d $dir) {  #  If really removed, remove it from our list.
                delete $directories{$dir};
                $deleted++;
            }
        }
    }

    print STDERR "--\n";
    print STDERR "-- Purged ", int(1000 * $bytes / 1024 / 1024 / 1024) / 1000, " GB in $files overlap output files and $dirs directories.\n";

    #  Now all done!

  finishStage:
    if ($tag eq "utg") {
        generateOverlapStoreStats($base, $asm);
    }

    if (-e "$base/$asm.ovlStore.summary") {
        print STDERR "--\n";
        print STDERR "-- Overlap store '$base/$asm.ovlStore' contains:\n";
        print STDERR "--\n";

        open(F, "< $base/$asm.ovlStore.summary") or caExit("Failed to open overlap store statistics in '$base/$asm.ovlStore': $!", undef);
        while (<F>) {
            print STDERR "--   $_";
        }
        close(F);

    } else {
        print STDERR "-- Overlap store '$base/$asm.ovlStore' statistics not available (skipped in correction and trimming stages).\n";
    }

    emitStage($asm, "$tag-createOverlapStore");
    buildHTML($asm, $tag);

  allDone:
    stopAfter("overlapStore");
}
