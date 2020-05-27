
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

package canu::OverlapMMap;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(mmapConfigure mmapPrecomputeCheck mmapCheck);

use strict;
use warnings "all";
no  warnings "uninitialized";

use File::Path 2.08 qw(make_path remove_tree);

use canu::Defaults;
use canu::Execution;

use canu::SequenceStore;
use canu::Report;

use canu::Grid_Cloud;

#  Map long reads to long reads with minimap.

sub setParameters ($$) {
    my $asm  = shift @_;
    my $base = shift @_;

    my ($numRaw,
        $numRawTri,
        $numCor,
        $numCorTri,
        $numPacBio,
        $numNanopore,
        $numHiFi) = getSequenceStoreStats($asm);

    my $parameters;

    if    (($numRaw > 0) && ($numCor == 0) && ($numPacBio > 0)) {
       $parameters = "-x ava-pb";
    }

    elsif (($numRaw > 0) && ($numCor == 0) && ($numNanopore > 0)) {
       $parameters = "-x ava-ont";
    }

    elsif (($numCor > 0) && ($numHiFi > 0 || $numPacBio > 0)) {
       $parameters = "-x ava-pb"; # -Hk21 -w14"; #tuned to find 1000bp 5% error
    }

    elsif (($numCor > 0) && ($numNanopore > 0)) {
       $parameters = "-x ava-ont"; # -k17 -w11"; #tuned to find 1000bp 15% error
    }

    else {
       caFailure("no known read types found in $base/$asm.seqStore/libraries.txt", undef)
    }

    return($parameters);
}



sub mmapConfigure ($$$) {
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $typ     = shift @_;
    my $bin     = getBinDirectory();
    my $minimap = getGlobal("minimap");

    my $base;                #  e.g., $base/1-overlapper/mhap.sh
    my $path;                #  e.g., $path/mhap.sh

    $base = "correction"  if ($tag eq "cor");
    $base = "trimming"    if ($tag eq "obt");
    $base = "unitigging"  if ($tag eq "utg");

    $path = "$base/1-overlapper";

    caFailure("invalid type '$typ'", undef)  if (($typ ne "partial") && ($typ ne "normal"));

    goto allDone   if (fileExists("$path/ovljob.files"));
    goto allDone   if ((-d "$base/$asm.ovlStore") || (fileExists("$base/$asm.ovlStore.tar.gz")));

    if (fileExists("$path/queries.tar")) {
        print STDERR "--\n";
        print STDERR "-- OVERLAPPER (mmap) (correction) complete, not rewriting scripts.\n"  if ($tag eq "cor");
        print STDERR "-- OVERLAPPER (mmap) (trimming) complete, not rewriting scripts.\n"    if ($tag eq "obt");
        print STDERR "-- OVERLAPPER (mmap) (assembly) complete, not rewriting scripts.\n"    if ($tag eq "utg");
        print STDERR "--\n";

        goto allDone;
    }

    my $parameters = setParameters($asm, $base);

    print STDERR "--\n";
    print STDERR "-- OVERLAPPER (mmap) (correction) with $parameters\n"  if ($tag eq "cor");
    print STDERR "-- OVERLAPPER (mmap) (trimming) with $parameters\n"    if ($tag eq "obt");
    print STDERR "-- OVERLAPPER (mmap) (assembly) with $parameters\n"    if ($tag eq "utg");
    print STDERR "--\n";

    make_path($path) if (! -d $path);

    #  Constants.

    my $numReads      = getNumberOfReadsInStore($asm, "all");   #  Need to iterate over all read IDs!
    my $memorySize    = getGlobal("${tag}mmapMemory");
    my $blockPerGb    = getGlobal("${tag}MMapBlockSize");
    my $blockSize = int($blockPerGb * $memorySize);

    print STDERR "-- Given $memorySize GB, can fit $blockSize reads per block.\n";

    #  Divide the reads into blocks of ovlHashBlockSize.  Each one of these blocks is used as the
    #  table in mmap.  Several of these blocks are used as the queries.

    my @blocks;    #  Range of reads to extract for this block
    my @blockBgn;  #  First read in the block
    my @blockLen;  #  Number of reads in the block

    my @hashes;    #  One for each job, the block that is the hash table
    my @skipSelf;  #  One for each job, jobs that would search block N against hash block N need to be handled special
    my @convert;   #  One for each job, flags to the mmap-ovl conversion program

    push @blocks,   "no zeroth block, makes the loop where this is used easier";
    push @blockBgn, "no zeroth block";
    push @blockLen, "no zeroth block";
    push @hashes,   "no zeroth job";
    push @skipSelf, "no zeroth job";

    for (my $bgn=1; $bgn < $numReads; $bgn += $blockSize) {
        my $end = $bgn + $blockSize - 1;
        $end = $numReads  if ($end > $numReads);

        #print STDERR "BLOCK ", scalar(@blocks), " reads from $bgn through $end\n";

        push @blocks, "-r $bgn-$end";

        push @blockBgn, $bgn;
        push @blockLen, $end - $bgn + 1;
    }

    #  Each mmap job will process one block against a set of other blocks.  We'll pick, arbitrarily,
    #  to use num_blocks/4 for that size, unless it is too small.

    my $numBlocks = scalar(@blocks);
    my $qryStride = ($numBlocks < 16) ? (2) : int($numBlocks / 4);

    print STDERR "-- For $numBlocks blocks, set stride to $qryStride blocks.\n";
    print STDERR "-- Logging partitioning to '$path/partitioning.log'.\n";

    open(L, "> $path/partitioning.log") or caExit("can't open '$path/partitioning.log' for writing: $!\n", undef);

    #  Make queries.  Each hask block needs to search against all blocks less than or equal to it.
    #  Each job will search at most $qryStride blocks at once.  So, queries could be:
    #  1:  1 vs 1,2,3  (with self-allowed, and quert block 1 implicitly included)
    #  2:  1 vs 4,5    (with no-self allowed)
    #  3:  2 vs 2,3,4
    #  4:  2 vs 5
    #  5:  3 vs 3,4,5
    #  6:  4 vs 4,5
    #  7:  5 vs 5

    make_path("$path/queries");

    for (my $bid=1; $bid < $numBlocks; $bid++) {

        #  Note that we never do qbgn = bid; the self-self overlap is special cased.

        for (my $qbgn = $bid; $qbgn < $numBlocks; $qbgn += $qryStride) {
            my $andSelf = "";

            if ($bid == $qbgn) {
                #  The hash block bid is in the range and we need to compute hash-to-hash
                #  overlaps, and exclude the block from the queries.
                push @skipSelf, "true";
                $qbgn++;
                $andSelf = " (and self)";
            } else {
                #  The hash block bid isn't in the range of query blocks, don't allow hash-to-hash
                #  overlaps.
                push @skipSelf, "";
            }

            my $qend = $qbgn + $qryStride - 1;                 #  Block bid searches reads in dat files from
            $qend = $numBlocks-1   if ($qend >= $numBlocks);   #  qbgn to qend (inclusive).


            my $job = substr("000000" . scalar(@hashes), -6);  #  Unique ID for this compute

            #  Make a place to save queries.  If this is the last-block-special-case, make a directory,
            #  but don't link in any files.  Without the directory, we'd need even more special case
            #  code down in mmap.sh to exclude the -q option for this last block.

            make_path("$path/queries/$job");

            if ($qbgn < $numBlocks) {
                print L "Job ", scalar(@hashes), " computes block $bid vs blocks $qbgn-$qend$andSelf,\n";

                for (my $qid=$qbgn; $qid <= $qend; $qid++) {
                    my $qry = substr("000000" . $qid, -6);             #  Name for the query block

                    symlink("../../blocks/$qry.fasta", "$path/queries/$job/$qry.fasta");
                }

            } else {
                print L "Job ", scalar(@hashes), " computes block $bid vs itself.\n";
                $qbgn = $bid;  #  Otherwise, the @convert -q value is bogus
            }

            #  This is easy, the ID of the hash.

            push @hashes, substr("000000" . $bid, -6);  #  One new job for block bid with qend-qbgn query files in it
        }
    }

    close(L);

    #  Create a script to generate precomputed blocks, including extracting the reads from seqStore.

    #OPTIMIZE
    #OPTIMIZE  Probably a big optimization for cloud assemblies, the block fasta inputs can be
    #OPTIMIZE  computed ahead of time, stashed, and then fetched to do the actual precompute.
    #OPTIMIZE

    open(F, "> $path/precompute.sh") or caFailure("can't open '$path/precompute.sh' for writing: $!", undef);

    print F "#!" . getGlobal("shell") . "\n";
    print F "\n";
    print F getBinDirectoryShellCode();
    print F "\n";
    print F setWorkDirectoryShellCode($path);
    print F fetchSeqStoreShellCode($asm, $path, "");
    print F "\n";
    print F getJobIDShellCode();
    print F "\n";
    for (my $ii=1; $ii < scalar(@blocks); $ii++) {
        print F "if [ \$jobid -eq $ii ] ; then\n";
        print F "  rge=\"$blocks[$ii]\"\n";
        print F "  job=\"", substr("000000" . $ii, -6), "\"\n";
        print F "fi\n";
        print F "\n";
    }
    print F "\n";
    print F "if [ x\$job = x ] ; then\n";
    print F "  echo Job partitioning error.  jobid \$jobid is invalid.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "\n";
    print F "if [ ! -d ./blocks ]; then\n";
    print F "  mkdir -p ./blocks\n";
    print F "fi\n";
    print F "\n";
    print F fileExistsShellCode("exists", $path, "blocks/\$job.mmi");
    print F "if [ \$exists = true ] ; then\n";
    print F "  echo Job previously completed successfully.\n";
    print F "  exit\n";
    print F "fi\n";
    print F "\n";
    print F "\$bin/sqStoreDumpFASTQ \\\n";
    print F "  -S ../../$asm.seqStore \\\n";
    print F "  \$rge \\\n";
    print F "  -nolibname \\\n";
    print F "  -noreadname \\\n";
    print F "  -fasta \\\n";
    print F "  -o ./blocks/\$job.input \\\n";
    print F "&& \\\n";
    print F "mv -f ./blocks/\$job.input.fasta ./blocks/\$job.fasta\n";
    print F "if [ ! -e ./blocks/\$job.fasta ] ; then\n";
    print F "  echo Failed to extract fasta.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "\n";
    print F "\n";
    print F "echo \"\"\n";
    print F "echo Starting minimap precompute.\n";
    print F "echo \"\"\n";
    print F "\n";
    print F "$minimap \\\n";
    print F "  $parameters -t ", getGlobal("${tag}mmapThreads"), " \\\n";
    print F "  -d ./blocks/\$job.input.mmi\\\n";
    print F "  ./blocks/\$job.fasta \\\n";
    print F "&& \\\n";
    print F "mv -f ./blocks/\$job.input.mmi ./blocks/\$job.mmi\n";
    print F "\n";
    print F "if [ ! -e ./blocks/\$job.mmi ] ; then\n";
    print F "  echo minimap failed.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "\n";
    print F stashFileShellCode("$base/1-overlapper/blocks", "\$job.mmi", "");
    print F "\n";
    print F "exit 0\n";

    close(F);

    #  Create a script to run mmap.

    open(F, "> $path/mmap.sh") or caFailure("can't open '$path/mmap.sh' for writing: $!", undef);

    print F "#!" . getGlobal("shell") . "\n";
    print F "\n";
    print F getBinDirectoryShellCode();
    print F "\n";
    print F setWorkDirectoryShellCode($path);
    print F fetchSeqStoreShellCode($asm, $path, "");
    print F "\n";
    print F getJobIDShellCode();
    print F "\n";
    for (my $ii=1; $ii < scalar(@hashes); $ii++) {
        print F "if [ \$jobid -eq $ii ] ; then\n";
        print F "  blk=\"$hashes[$ii]\"\n";
        print F "  slf=\"$skipSelf[$ii]\"\n";
        print F "  qry=\"", substr("000000" . $ii, -6), "\"\n";
        print F "fi\n";
        print F "\n";
    }

    print F "\n";
    print F "if [ x\$qry = x ]; then\n";
    print F "  echo Error: Job index out of range.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "\n";
    print F "if [ -e ./results/\$qry.ovb ]; then\n";
    print F "  echo Job previously completed successfully.\n";
    print F "  exit\n";
    print F "fi\n";
    print F "\n";

    print F fetchFileShellCode("$path", "queries.tar", "");
    print F "\n";
    print F "if [ -e ./queries.tar -a ! -d ./queries ] ; then\n";
    print F "  tar -xf ./queries.tar\n";
    print F "fi\n";
    print F "\n";

    print F "if [ ! -d ./results ]; then\n";
    print F "  mkdir -p ./results\n";
    print F "fi\n";
    print F "\n";

    print F fetchFileShellCode("$path", "blocks/\$blk.mmi", "");

    print F "for ii in `ls ./queries/\$qry` ; do\n";
    print F "  echo Fetch blocks/\$ii\n";
    print F    fetchFileShellCode("$path", "blocks/\$ii", "  ");
    print F "done\n";
    print F "\n";

    #  Begin comparison, we loop through query and compare current block to it, if we need to do
    #  self first compare to self, otherwise initialize as empty

    print F "if [ x\$slf = x ]; then\n";
    print F "   >  ./results/\$qry.mmap.WORKING\n";
    print F "else\n";
    print F "  $minimap \\\n";
    print F "    $parameters -t ", getGlobal("${tag}mmapThreads"), " \\\n";
    print F "    ./blocks/\$blk.mmi \\\n";
    print F "    ./blocks/\$blk.fasta \\\n";
    print F "  > ./results/\$qry.mmap.WORKING \n";
    print F " \n";
    print F "fi\n";
    print F "\n";

    print F "for file in `ls queries/\$qry/*.fasta`; do\n";
    print F "  $minimap \\\n";
    print F "    $parameters -t ", getGlobal("${tag}mmapThreads"), " \\\n";
    print F "    ./blocks/\$blk.mmi \\\n";
    print F "    \$file \\\n";
    print F "  >> ./results/\$qry.mmap.WORKING \n";
    print F "done\n";
    print F "\n";
    print F "mv  ./results/\$qry.mmap.WORKING  ./results/\$qry.mmap\n";
    print F "\n";

    print F "if [   -e ./results/\$qry.mmap -a \\\n";
    print F "     ! -e ./results/\$qry.ovb ] ; then\n";
    print F "  \$bin/mmapConvert \\\n";
    print F "    -S ../../$asm.seqStore \\\n";
    print F "    -o ./results/\$qry.mmap.ovb.WORKING \\\n";
    print F "    -e " . getGlobal("${tag}OvlErrorRate");
    print F "    -partial \\\n"  if ($typ eq "partial");
    print F "    -len "  , getGlobal("minOverlapLength"),  " \\\n";
    print F "    ./results/\$qry.mmap \\\n";
    print F "  && \\\n";
    print F "  mv ./results/\$qry.mmap.ovb.WORKING ./results/\$qry.mmap.ovb\n";
    print F "fi\n";
    print F "\n";

    if (getGlobal('purgeOverlaps') ne "never") {
        print F "if [   -e ./results/\$qry.mmap -a \\\n";
        print F "       -e ./results/\$qry.mmap.ovb ] ; then\n";
        print F "  rm -f ./results/\$qry.mmap\n";
        print F "fi\n";
        print F "\n";
    }

    print F "if [ -e ./results/\$qry.mmap.ovb ] ; then\n";
    if (getGlobal("${tag}ReAlign") eq "1") {
        print F "  \$bin/overlapPair \\\n";
        print F "    -S ../../$asm.seqStore \\\n";
        print F "    -O ./results/\$qry.mmap.ovb \\\n";
        print F "    -o ./results/\$qry.ovb \\\n";
        print F "    -partial \\\n"  if ($typ eq "partial");
        print F "    -erate ", getGlobal("corOvlErrorRate"), " \\\n"  if ($tag eq "cor");   #  Explicitly using proper name for grepability.
        print F "    -erate ", getGlobal("obtOvlErrorRate"), " \\\n"  if ($tag eq "obt");
        print F "    -erate ", getGlobal("utgOvlErrorRate"), " \\\n"  if ($tag eq "utg");
        print F "    -memory " . getGlobal("${tag}mmapMemory") . " \\\n";
        print F "    -t " . getGlobal("${tag}mmapThreads") . " \n";
    } else {
        print F "  mv -f ./results/\$qry.mmap.ovb ./results/\$qry.ovb\n";
        print F "  mv -f ./results/\$qry.mmap.oc  ./results/\$qry.oc\n";
    }
    print F "fi\n";

    print F stashFileShellCode("$path", "results/\$qry.ovb", "");
    print F stashFileShellCode("$path", "results/\$qry.oc",  "");
    print F "\n";

    print F "\n";
    print F "exit 0\n";

    close(F);

    if (-e "$path/precompute.sh") {
        my $numJobs = 0;
        open(F, "< $path/precompute.sh") or caFailure("can't open '$path/precompute.sh' for reading: $!", undef);
        while (<F>) {
            $numJobs++   if (m/^\s+job=/);
        }
        close(F);

        print STDERR "-- Configured $numJobs mmap precompute jobs.\n";
    }

    if (-e "$path/mmap.sh") {
        my $numJobs = 0;
        open(F, "< $path/mmap.sh") or caFailure("can't open '$path/mmap.sh' for reading: $!", undef);
        while (<F>) {
            $numJobs++  if (m/^\s+qry=/);
        }
        close(F);

        print STDERR "-- Configured $numJobs mmap overlap jobs.\n";
    }

    #  Tar up the queries directory and save scripts.  The queries.tar is
    #  useful for cloud support, but we also use it to decide if this step
    #  has finished.

    runCommandSilently($path, "tar -cf queries.tar queries", 1);

    makeExecutable("$path/precompute.sh");
    makeExecutable("$path/mmap.sh");

    stashFile("$path/precompute.sh");
    stashFile("$path/mmap.sh");
    stashFile("$path/queries.tar");

  finishStage:
    generateReport($asm);
    resetIteration("$tag-mmapConfigure");

  allDone:
    stopAfter("overlapConfigure");
}


sub mmapPrecomputeCheck ($$$) {
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $typ     = shift @_;
    my $attempt = getGlobal("canuIteration");

    my $base;                #  e.g., $base/1-overlapper/mhap.sh
    my $path;                #  e.g., $path/mhap.sh

    $base = "correction"  if ($tag eq "cor");
    $base = "trimming"    if ($tag eq "obt");
    $base = "unitigging"  if ($tag eq "utg");

    $path = "$base/1-overlapper";

    goto allDone   if (fileExists("$path/precompute.files"));
    goto allDone   if (fileExists("$path/mmap.files"));
    goto allDone   if ((-d "$base/$asm.ovlStore") || (fileExists("$base/$asm.ovlStore.tar.gz")));

    fetchFile("$path/precompute.sh");

    #  Figure out if all the tasks finished correctly.

    my $currentJobID   = 1;
    my @successJobs;
    my @failedJobs;
    my $failureMessage = "";

    open(F, "< $path/precompute.sh") or caFailure("can't open '$path/precompute.sh' for reading: $!", undef);
    while (<F>) {
        if (m/^\s+job=\"(\d+)\"$/) {
            if (fileExists("$path/blocks/$1.mmi")) {
                push @successJobs, "1-overlapper/blocks/$1.mmi\\n";
            } else {
                $failureMessage .= "--   job $path/blocks/$1.fasta FAILED.\n";
                push @failedJobs, $currentJobID;
            }

            $currentJobID++;
        }
    }
    close(F);

    #  Failed jobs, retry.

    if (scalar(@failedJobs) > 0) {

        #  If too many attempts, give up.

        if ($attempt >= getGlobal("canuIterationMax")) {
            print STDERR "--\n";
            print STDERR "-- MiniMap precompute jobs failed, tried $attempt times, giving up.\n";
            print STDERR $failureMessage;
            print STDERR "--\n";
            caExit(undef, undef);
        }

        if ($attempt > 0) {
            print STDERR "--\n";
            print STDERR "-- MiniMap precompute jobs failed, retry.\n";
            print STDERR $failureMessage;
            print STDERR "--\n";
        }

        #  Otherwise, run some jobs.

        generateReport($asm);

        submitOrRunParallelJob($asm, "${tag}mmap", $path, "precompute", @failedJobs);
        return;
    }

  finishStage:
    print STDERR "-- All ", scalar(@successJobs), " mmap precompute jobs finished successfully.\n";

    open(L, "> $path/precompute.files") or caExit("failed to open '$path/precompute.files'", undef);
    print L @successJobs;
    close(L);

    stashFile("$path/precompute.files");

    generateReport($asm);
    resetIteration("$tag-mmapPrecomputeCheck");

  allDone:
}



sub mmapCheck ($$$) {
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $typ     = shift @_;
    my $attempt = getGlobal("canuIteration");

    my $base;                #  e.g., $base/1-overlapper/mhap.sh
    my $path;                #  e.g., $path/mhap.sh

    $base = "correction"  if ($tag eq "cor");
    $base = "trimming"    if ($tag eq "obt");
    $base = "unitigging"  if ($tag eq "utg");

    $path = "$base/1-overlapper";

    goto allDone   if (fileExists("$path/mmap.files"));
    goto allDone   if ((-d "$base/$asm.ovlStore") || (fileExists("$base/$asm.ovlStore.tar.gz")));

    fetchFile("$path/mmap.sh");

    #  Figure out if all the tasks finished correctly.

    my $currentJobID   = 1;
    my @mmapJobs;
    my @successJobs;
    my @miscJobs;
    my @failedJobs;
    my $failureMessage = "";

    open(F, "< $path/mmap.sh") or caExit("failed to open '$path/mmap.sh'", undef);
    while (<F>) {
        if (m/^\s+qry=\"(\d+)\"$/) {
            if      (fileExists("$path/results/$1.ovb.gz")) {
                push @mmapJobs,    "1-overlapper/results/$1.mmap\n";
                push @successJobs, "1-overlapper/results/$1.ovb.gz\n";
                push @miscJobs,    "1-overlapper/results/$1.stats\n";
                push @miscJobs,    "1-overlapper/results/$1.oc\n";

            } elsif (fileExists("$path/results/$1.ovb")) {
                push @mmapJobs,    "1-overlapper/results/$1.mmap\n";
                push @successJobs, "1-overlapper/results/$1.ovb\n";
                push @miscJobs,    "1-overlapper/results/$1.stats\n";
                push @miscJobs,    "1-overlapper/results/$1.oc\n";

            } elsif (fileExists("$path/results/$1.ovb.bz2")) {
                push @mmapJobs,    "1-overlapper/results/$1.mmap\n";
                push @successJobs, "1-overlapper/results/$1.ovb.bz2\n";
                push @miscJobs,    "1-overlapper/results/$1.stats\n";
                push @miscJobs,    "1-overlapper/results/$1.oc\n";

            } elsif (fileExists("$path/results/$1.ovb.xz")) {
                push @mmapJobs,    "1-overlapper/results/$1.mmap\n";
                push @successJobs, "1-overlapper/results/$1.ovb.xz\n";
                push @miscJobs,    "1-overlapper/results/$1.stats\n";
                push @miscJobs,    "1-overlapper/results/$1.oc\n";

            } else {
                $failureMessage .= "--   job $path/results/$1.ovb FAILED.\n";
                push @failedJobs, $currentJobID;
            }

            $currentJobID++;
        }
    }
    close(F);

    #  Also find the queries symlinks so we can remove those.  And the query directories, because
    #  the last directory can be empty, and so we'd never see it at all if only finding files.

    open(F, "cd $base && find 1-overlapper/queries -print |");
    while (<F>) {
        push @mmapJobs, $_;
    }
    close(F);

    #  Failed jobs, retry.

    if (scalar(@failedJobs) > 0) {

        #  If too many attempts, give up.

        if ($attempt >= getGlobal("canuIterationMax")) {
            print STDERR "--\n";
            print STDERR "-- MiniMap overlap jobs failed, tried $attempt times, giving up.\n";
            print STDERR $failureMessage;
            print STDERR "--\n";
            caExit(undef, undef);
        }

        if ($attempt > 0) {
            print STDERR "--\n";
            print STDERR "-- MiniMap overlap jobs failed, retry.\n";
            print STDERR $failureMessage;
            print STDERR "--\n";
        }

        #  Otherwise, run some jobs.

        generateReport($asm);

        submitOrRunParallelJob($asm, "${tag}mmap", $path, "mmap", @failedJobs);
        return;
    }

  finishStage:
    print STDERR "-- Found ", scalar(@successJobs), " mmap overlap output files.\n";

    open(L, "> $path/mmap.files") or caExit("failed to open '$path/mmap.files'", undef);
    print L @mmapJobs;
    close(L);

    open(L, "> $path/ovljob.files") or caExit("failed to open '$path/ovljob.files'", undef);
    print L @successJobs;
    close(L);

    open(L, "> $path/ovljob.more.files") or caExit("failed to open '$path/ovljob.more.files'", undef);
    print L @miscJobs;
    close(L);

    stashFile("$path/mmap.files");
    stashFile("$path/ovljob.files");
    stashFile("$path/ovljob.more.files");

    generateReport($asm);
    resetIteration("$tag-mmapCheck");

  allDone:
    stopAfter("overlap");
}

