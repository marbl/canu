
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

use POSIX qw(ceil);
use canu::Defaults;
use canu::Execution;
use canu::HTML;


#  Parallel documentation: Each overlap job is converted into a single bucket of overlaps.  Within
#  each bucket, the overlaps are distributed into many slices, one per sort job.  The sort jobs then
#  load the same slice from each bucket.


#  NOT FILTERING overlaps by error rate when building the parallel store.


sub createOverlapStoreSequential ($$$$) {
    my $WRK     = shift @_;  #  Root work directory (the -d option to canu)
    my $wrk     = $WRK;      #  Local work directory
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $files   = shift @_;
    my $bin     = getBinDirectory();
    my $cmd;

    $wrk = "$wrk/correction"  if ($tag eq "cor");
    $wrk = "$wrk/trimming"    if ($tag eq "obt");
    $wrk = "$wrk/unitigging"  if ($tag eq "utg");

    #  This is running in the canu process itself.  Execution.pm has special case code
    #  to submit canu to grids using the maximum of 4gb and this memory limit.

    my $memSize = getGlobal("ovsMemory");

    #  The parallel store build will unlimit 'max user processes'.  The sequential method usually
    #  runs out of open file handles first (meaning it has never run out of processes yet).

    $cmd  = "$bin/ovStoreBuild \\\n";
    $cmd .= " -O $wrk/$asm.ovlStore.BUILDING \\\n";
    $cmd .= " -G $wrk/$asm.gkpStore \\\n";
    $cmd .= " -M $memSize \\\n";
    $cmd .= " -L $files \\\n";
    $cmd .= " > $wrk/$asm.ovlStore.err 2>&1";

    if (runCommand($wrk, $cmd)) {
        caExit("failed to create the overlap store", "$wrk/$asm.ovlStore.err");
    }

    unlink "$wrk/$asm.ovlStore.err";

    rename "$wrk/$asm.ovlStore.BUILDING", "$wrk/$asm.ovlStore";
}




#  Count the number of inputs.  We don't expect any to be missing (they were just checked
#  by overlapCheck()) but feel silly not checking again.

sub countOverlapStoreInputs ($) {
    my $inputs    = shift @_;
    my $numInputs = 0;

    open(F, "< $inputs") or die "Failed to open overlap store input file '$inputs': $0\n";
    while (<F>) {
        chomp;
        die "overlapper output '$_' not found\n"  if (! -e $_);
        $numInputs++;
    }
    close(F);

    return($numInputs);
}




sub getNumOlapsAndSlices ($$) {
    my $wrk = shift @_;
    my $asm = shift @_;

    my $numOlaps   = undef;
    my $numSlices  = undef;
    my $memLimit   = undef;

    open(F, "< $wrk/$asm.ovlStore.BUILDING/config.err") or caExit("can't open '$wrk/$asm.ovlStore.BUILDING/config.err' for reading: $!\n", undef);
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
    my $WRK     = shift @_;  #  Root work directory (the -d option to canu)
    my $wrk     = $WRK;      #  Local work directory
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $files   = shift @_;
    my $bin     = getBinDirectory();
    my $cmd;

    $wrk = "$wrk/correction"  if ($tag eq "cor");
    $wrk = "$wrk/trimming"    if ($tag eq "obt");
    $wrk = "$wrk/unitigging"  if ($tag eq "utg");

    goto allDone   if (skipStage($WRK, $asm, "$tag-overlapStoreConfigure") == 1);
    goto allDone   if (-d "$wrk/$asm.ovlStore");

    my $numInputs  = countOverlapStoreInputs($files);

    #  Create an output directory, and populate it with more directories and scripts

    system("mkdir -p $wrk/$asm.ovlStore.BUILDING")                   if (! -d "$wrk/$asm.ovlStore.BUILDING");
    system("mkdir -p $wrk/$asm.ovlStore.BUILDING/scripts")           if (! -d "$wrk/$asm.ovlStore.BUILDING/scripts");
    system("mkdir -p $wrk/$asm.ovlStore.BUILDING/logs")              if (! -d "$wrk/$asm.ovlStore.BUILDING/logs");

    #  Run the normal store build, but just to get the partitioning.

    if (! -e "$wrk/$asm.ovlStore.BUILDING/config") {
        $cmd  = "$bin/ovStoreBuild \\\n";
        $cmd .= " -G $wrk/$asm.gkpStore \\\n";
        $cmd .= " -O $wrk/$asm.ovlStore \\\n";  #  NOT created!
        $cmd .= " -M " . getGlobal("ovsMemory") . " \\\n";
        $cmd .= " -config $wrk/$asm.ovlStore.BUILDING/config \\\n";
        $cmd .= " -L $files \\\n";
        $cmd .= "> $wrk/$asm.ovlStore.BUILDING/config.err 2>&1\n";

        if (runCommand($wrk, $cmd)) {
            caExit("failed to generate configuration for building overlap store", "$wrk/$asm.ovlStore.BUILDING/config.err");
        }
    }

    #  Parse the output to find the number of jobs we need to sort and the memory
    #  ovs store memory is left as a range (e.g. 4-16) so building can scale itself to (hopefully) fit both into memory and into max system open files
    my ($numOlaps, $numSlices, $memLimit) = getNumOlapsAndSlices($wrk, $asm);
    setGlobal("ovsMemory", ceil($memLimit * 1.1)); # request more memory to avoid memory issues on machines with no swap

    #  Parallel jobs for bucketizing.  This should really be part of overlap computation itself.

    #getAllowedResources("", "ovb");

    if (! -e "$wrk/$asm.ovlStore.BUILDING/scripts/1-bucketize.sh") {
        open(F, "> $wrk/$asm.ovlStore.BUILDING/scripts/1-bucketize.sh") or die;
        print F "#!" . getGlobal("shell") . "\n";
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
        print F "bn=`printf %04d \$jobid`\n";
        print F "jn=\"undefined\"\n";
        print F "\n";

        my $tstid = 1;

        open(I, "< $files") or die "Failed to open '$files': $0\n";

        while (<I>) {
            chomp;

            print F "if [ \"\$jobid\" -eq \"$tstid\" ] ; then jn=\"$_\"; fi\n";
            $tstid++;
        }

        close(I);

        print F "\n";
        print F "if [ \$jn = \"undefined\" ] ; then\n";
        print F "  echo \"Job out of range.\"\n";
        print F "  exit\n";
        print F "fi\n";
        print F "\n";
        print F "if [ -e \"$wrk/$asm.ovlStore.BUILDING/bucket\$bn/sliceSizes\" ] ; then\n";
        print F "  echo \"Bucket $wrk/$asm.ovlStore.BUILDING/bucket\$bn finished successfully.\"\n";
        print F "  exit\n";
        print F "fi\n";
        print F "\n";
        print F "if [ -e \"$wrk/$asm.ovlStore.BUILDING/create\$bn\" ] ; then\n";
        print F "  echo \"Removing incomplete bucket $wrk/$asm.ovlStore.BUILDING/create\$bn\"\n";
        print F "  rm -rf \"$wrk/$asm.ovlStore.BUILDING/create\$bn\"\n";
        print F "fi\n";
        print F "\n";
        print F "max=`ulimit -Hu`\n";
        print F "bef=`ulimit -Su`\n";
        print F "if [ \$bef -lt \$max ] ; then\n";
        print F "  ulimit -Su \$max\n";
        print F "  aft=`ulimit -Su`\n";
        print F "  echo \"Changed max processes per user from \$bef to \$aft (max \$max).\"\n";
        print F "  echo \"\"\n";
        print F "else\n";
        print F "  echo \"Max processes per user limited to \$bef, no increase possible.\"\n";
        print F "  echo \"\"\n";
        print F "fi\n";
        print F "\n";
        print F getBinDirectoryShellCode();
        print F "\n";
        print F "\$bin/ovStoreBucketizer \\\n";
        print F "  -O $wrk/$asm.ovlStore.BUILDING \\\n";
        print F "  -G $wrk/$asm.gkpStore \\\n";
        print F "  -C $wrk/$asm.ovlStore.BUILDING/config \\\n";
        #print F "  -e " . getGlobal("") . " \\\n"  if (defined(getGlobal("")));
        print F "  -job \$jobid \\\n";
        print F "  -i   \$jn\n";
        close(F);
    }

    #  Parallel jobs for sorting each bucket

    #getAllowedResources("", "ovs");

    if (! -e "$wrk/$asm.ovlStore.BUILDING/scripts/2-sort.sh") {
        open(F, "> $wrk/$asm.ovlStore.BUILDING/scripts/2-sort.sh") or die;
        print F "#!" . getGlobal("shell") . "\n";
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
        print F "max=`ulimit -Hu`\n";
        print F "bef=`ulimit -Su`\n";
        print F "if [ \$bef -lt \$max ] ; then\n";
        print F "  ulimit -Su \$max\n";
        print F "  aft=`ulimit -Su`\n";
        print F "  echo \"Changed max processes per user from \$bef to \$aft (max \$max).\"\n";
        print F "  echo \"\"\n";
        print F "else\n";
        print F "  echo \"Max processes per user limited to \$bef, no increase possible.\"\n";
        print F "  echo \"\"\n";
        print F "fi\n";
        print F "\n";
        print F getBinDirectoryShellCode();
        print F "\n";
        print F "\$bin/ovStoreSorter \\\n";
        print F "  -deletelate \\\n";  #  Choices -deleteearly -deletelate or nothing
        print F "  -M $memLimit \\\n";
        print F "  -O $wrk/$asm.ovlStore.BUILDING \\\n";
        print F "  -G $wrk/$asm.gkpStore \\\n";
        print F "  -F $numSlices \\\n";
        print F "  -job \$jobid $numInputs\n";
        print F "\n";
        print F "if [ \$? = 0 ] ; then\n";
        print F "  echo Success.\n";
        print F "else\n";
        print F "  echo Failure.\n";
        print F "fi\n";
        close(F);
    }

    #  A final job to merge the indices.

    if (! -e "$wrk/$asm.ovlStore.BUILDING/scripts/3-index.sh") {
        open(F, "> $wrk/$asm.ovlStore.BUILDING/scripts/3-index.sh") or die;
        print F "#!" . getGlobal("shell") . "\n";
        print F "\n";
        print F getBinDirectoryShellCode();
        print F "\n";
        print F "\$bin/ovStoreIndexer \\\n";
        #print F "  -nodelete \\\n";  #  Choices -nodelete or nothing
        print F "  -O $wrk/$asm.ovlStore.BUILDING \\\n";
        print F "  -F $numSlices\n";
        print F "\n";
        print F "if [ \$? = 0 ] ; then\n";
        print F "  echo Success.\n";
        print F "else\n";
        print F "  echo Failure.\n";
        print F "fi\n";
        close(F);
    }

    system("chmod +x $wrk/$asm.ovlStore.BUILDING/scripts/1-bucketize.sh");
    system("chmod +x $wrk/$asm.ovlStore.BUILDING/scripts/2-sort.sh");
    system("chmod +x $wrk/$asm.ovlStore.BUILDING/scripts/3-index.sh");

  finishStage:
    emitStage($WRK, $asm, "$tag-overlapStoreConfigure");
    buildHTML($WRK, $asm, $tag);
    stopAfter("overlapStoreConfigure");

  allDone:
}



sub overlapStoreBucketizerCheck ($$$$) {
    my $WRK     = shift @_;  #  Root work directory (the -d option to canu)
    my $wrk     = $WRK;      #  Local work directory
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $files   = shift @_;
    my $attempt = getGlobal("canuIteration");

    $wrk = "$wrk/correction"  if ($tag eq "cor");
    $wrk = "$wrk/trimming"    if ($tag eq "obt");
    $wrk = "$wrk/unitigging"  if ($tag eq "utg");

    goto allDone   if (skipStage($WRK, $asm, "$tag-overlapStoreBucketizerCheck", $attempt) == 1);
    goto allDone   if (-d "$wrk/$asm.ovlStore");
    goto allDone   if (-e "$wrk/$asm.ovlStore.BUILDING/1-bucketize.success");

    #  Figure out if all the tasks finished correctly.

    my $numInputs      = countOverlapStoreInputs($files);
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

        if (! -e "$wrk/$asm.ovlStore.BUILDING/bucket$bucketID") {
            $failureMessage .= "--   job $wrk/$asm.ovlStore.BUILDING/bucket$bucketID FAILED.\n";
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

        emitStage($WRK, $asm, "$tag-overlapStoreBucketizerCheck", $attempt);
        buildHTML($WRK, $asm, $tag);

        submitOrRunParallelJob($WRK, $asm, "ovB", "$wrk/$asm.ovlStore.BUILDING", "scripts/1-bucketize", @failedJobs);
        return;
    }

  finishStage:
    print STDERR "-- Overlap store bucketizer finished.\n";

    touch("$wrk/$asm.ovlStore.BUILDING/1-bucketize.success");

    setGlobal("canuIteration", 0);
    emitStage($WRK, $asm, "$tag-overlapStoreBucketizerCheck");
    buildHTML($WRK, $asm, $tag);
    stopAfter("overlapBucketizer");

  allDone:
}





sub overlapStoreSorterCheck ($$$$) {
    my $WRK     = shift @_;  #  Root work directory (the -d option to canu)
    my $wrk     = $WRK;      #  Local work directory
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $files   = shift @_;
    my $attempt = getGlobal("canuIteration");

    $wrk = "$wrk/correction"  if ($tag eq "cor");
    $wrk = "$wrk/trimming"    if ($tag eq "obt");
    $wrk = "$wrk/unitigging"  if ($tag eq "utg");

    goto allDone   if (skipStage($WRK, $asm, "$tag-overlapStoreSorterCheck", $attempt) == 1);
    goto allDone   if (-d "$wrk/$asm.ovlStore");
    goto allDone   if (-e "$wrk/$asm.ovlStore.BUILDING/2-sorter.success");

    #  Figure out if all the tasks finished correctly.

    my ($numOlaps, $numSlices) = getNumOlapsAndSlices($wrk, $asm);

    my $currentJobID   = 1;
    my @successJobs;
    my @failedJobs;
    my $failureMessage = "";

    my $sortID       = "0001";

    open(F, "< $files") or caExit("can't open '$files' for reading: $!", undef);

    #  A valid result has three files:
    #    $wrk/$asm.ovlStore.BUILDING/$sortID
    #    $wrk/$asm.ovlStore.BUILDING/$sortID.index
    #    $wrk/$asm.ovlStore.BUILDING/$sortID.info
    #
    #  A crashed result has one file, if it crashes before output
    #    $wrk/$asm.ovlStore.BUILDING/$sortID.ovs
    #
    #  On out of disk, the .info is missing.  It's the last thing created.
    #
    while ($currentJobID <= $numSlices) {

        if ((! -e "$wrk/$asm.ovlStore.BUILDING/$sortID") ||
            (! -e "$wrk/$asm.ovlStore.BUILDING/$sortID.info") ||
            (  -e "$wrk/$asm.ovlStore.BUILDING/$sortID.ovs")) {
            $failureMessage .= "--   job $wrk/$asm.ovlStore.BUILDING/$sortID FAILED.\n";
            unlink "$wrk/$asm.ovlStore.BUILDING/$sortID.ovs";
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

        emitStage($WRK, $asm, "$tag-overlapStoreSorterCheck", $attempt);
        buildHTML($WRK, $asm, $tag);

        submitOrRunParallelJob($WRK, $asm, "ovS", "$wrk/$asm.ovlStore.BUILDING", "scripts/2-sort", @failedJobs);
        return;
    }

  finishStage:
    print STDERR "-- Overlap store sorter finished.\n";

    touch("$wrk/$asm.ovlStore.BUILDING/2-sorter.success");

    setGlobal("canuIteration", 0);
    emitStage($WRK, $asm, "$tag-overlapStoreSorterCheck");
    buildHTML($WRK, $asm, $tag);
    stopAfter("overlapSorter");

  allDone:
}




sub createOverlapStoreParallel ($$$$) {
    my $WRK     = shift @_;  #  Root work directory (the -d option to canu)
    my $wrk     = $WRK;      #  Local work directory
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $files   = shift @_;

    $wrk = "$wrk/correction"  if ($tag eq "cor");
    $wrk = "$wrk/trimming"    if ($tag eq "obt");
    $wrk = "$wrk/unitigging"  if ($tag eq "utg");

    overlapStoreConfigure($WRK, $asm, $tag, $files);
    overlapStoreBucketizerCheck($WRK, $asm, $tag, $files)  foreach (1..getGlobal("canuIterationMax") + 1);
    overlapStoreSorterCheck($WRK, $asm, $tag, $files)      foreach (1..getGlobal("canuIterationMax") + 1);

    if (runCommand("$wrk/$asm.ovlStore.BUILDING", "$wrk/$asm.ovlStore.BUILDING/scripts/3-index.sh > $wrk/$asm.ovlStore.BUILDING/logs/3-index.err 2>&1")) {
        caExit("failed to build index for overlap store", "$wrk/$asm.ovlStore.BUILDING/logs/3-index.err");
    }

    rename "$wrk/$asm.ovlStore.BUILDING", "$wrk/$asm.ovlStore";
}


sub generateOverlapStoreStats ($$) {
    my $WRK     = shift @_;  #  Root work directory (the -d option to canu)
    my $wrk     = $WRK;      #  Local work directory
    my $asm     = shift @_;

    my $bin   = getBinDirectory();
    my $cmd;

    $cmd  = "$bin/ovStoreStats \\\n";
    $cmd .= " -G $wrk/$asm.gkpStore \\\n";
    $cmd .= " -O $wrk/$asm.ovlStore \\\n";
    $cmd .= " -o $wrk/$asm.ovlStore \\\n";
    $cmd .= " > $wrk/$asm.ovlStore.summary.err 2>&1";

    if (runCommand($wrk, $cmd)) {
        print STDERR "--\n";
        print STDERR "-- WARNING: failed to generate statistics for the overlap store; no summary will appear in HTML output.\n";
        print STDERR "--\n";
        print STDERR "----------------------------------------\n";
    }
}


sub createOverlapStore ($$$$) {
    my $WRK     = shift @_;  #  Root work directory (the -d option to canu)
    my $wrk     = $WRK;      #  Local work directory
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $seq     = shift @_;

    $wrk = "$wrk/correction"  if ($tag eq "cor");
    $wrk = "$wrk/trimming"    if ($tag eq "obt");
    $wrk = "$wrk/unitigging"  if ($tag eq "utg");

    my $path  = "$wrk/1-overlapper";

    goto allDone   if (skipStage($WRK, $asm, "$tag-createOverlapStore") == 1);
    goto allDone   if (-d "$wrk/$asm.ovlStore");
    goto allDone   if (-d "$wrk/$asm.tigStore");

    #  Did we _really_ complete?

    caExit("overlapper claims to be finished, but no job list found in '$path/ovljob.files'", undef)  if (! -e "$path/ovljob.files");

    #  Then just build the store!  Simple!

    createOverlapStoreSequential($WRK, $asm, $tag, "$path/ovljob.files")  if ($seq eq "sequential");
    createOverlapStoreParallel  ($WRK, $asm, $tag, "$path/ovljob.files")  if ($seq eq "parallel");

    print STDERR "--\n";
    print STDERR "-- Overlap store '$wrk/$asm.ovlStore' successfully constructed.\n";

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

            if (-d $_) {
                $directories{$_}++;

            } elsif (-e $_) {
                $bytes += -s $_;
                $files += 1;

                unlink $_;
            }

            #  If the file isn't a directory -- the file itself could not exist above -- register
            #  the directory it is in for deletion.

            if (! -d $_) {
                my @components = split '/', $_;
                pop @components;
                my $dir = join '/', @components;

                $directories{$dir}++;
            }

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
    generateOverlapStoreStats($wrk, $asm);

    if (-e "$wrk/$asm.ovlStore.summary") {
        print STDERR "--\n";
        print STDERR "-- Overlap store '$wrk/$asm.ovlStore' contains:\n";
        print STDERR "--\n";

        open(F, "< $wrk/$asm.ovlStore.summary") or caExit("Failed to open overlap store statistics in '$wrk/$asm.ovlStore': $!", undef);
        while (<F>) {
            print STDERR "--   $_";
        }
        close(F);

    } else {
        print STDERR "-- Overlap store '$wrk/$asm.ovlStore' statistics not available.\n";
    }

    emitStage($WRK, $asm, "$tag-createOverlapStore");
    buildHTML($WRK, $asm, $tag);
    stopAfter("overlapStore");

  allDone:
}
