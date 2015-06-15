package ca3g::OverlapMhap;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(mhapConfigure mhapPrecomputeCheck mhapCheck);

use strict;

use File::Path qw(make_path remove_tree);

use ca3g::Defaults;
use ca3g::Execution;
use ca3g::Gatekeeper;

#  Map long reads to long reads with mhap.

#  Problems:
#  - mhap .dat output can't be verified.  if the job is killed, a partial .dat is output.  It needs to write to
#    .dat.WORKING, then let the script rename it when the job finishes successfully.  There is no output file name option.



sub mhapConfigure ($$$$) {
    my $WRK  = shift @_;
    my $wrk  = $WRK;
    my $asm  = shift @_;
    my $tag  = shift @_;
    my $typ  = shift @_;
    my $bin  = getBinDirectory();

    $wrk = "$wrk/correction"  if ($tag eq "cor");
    $wrk = "$wrk/trimming"    if ($tag eq "obt");

    my $path = "$wrk/1-overlapper";

    caFailure("invalid type '$typ'", undef)  if (($typ ne "partial") && ($typ ne "normal"));

    return  if (skipStage($WRK, $asm, "$tag-mhapConfigure") == 1);
    return  if (-e "$path/ovljob.files");
    return  if (-e "$path/mhap.sh");

    make_path("$path") if (! -d "$path");

    #  Constants.

    my $merSize       = getGlobal("${tag}MhapMerSize");

    my $taskID        = getGlobal("gridEngineTaskID");
    my $submitTaskID  = getGlobal("gridEngineArraySubmitID");

    my $numReads      = getNumberOfReadsInStore($wrk, $asm);
    my $blockSize     = getGlobal("${tag}MhapBlockSize");

    #  Divide the reads into blocks of ovlHashBlockSize.  Each one of these blocks is used as the
    #  table in mhap.  Several of these blocks are used as the queries.

    my @blocks;    #  Range of reads to extract for this block
    my @blockBgn;  #  First read in the block
    my @blockLen;  #  Number of reads in the block

    my @hashes;    #  One for each job, the block that is the hash table
    my @skipSelf;  #  One for each job, jobs that would search block N against hash block N need to be handled special
    my @convert;   #  One for each job, flags to the mhap-ovl conversion program

    push @blocks,   "no zeroth block, makes the loop where this is used easier";
    push @blockBgn, "no zeroth block";
    push @blockLen, "no zeroth block";
    push @hashes,   "no zeroth job";
    push @skipSelf, "no zeroth job";
    push @convert,  "no zeroth job";

    for (my $bgn=1; $bgn < $numReads; $bgn += $blockSize) {
        my $end = $bgn + $blockSize - 1;
        $end = $numReads  if ($end > $numReads);

        #print STDERR "BLOCK ", scalar(@blocks), " reads from $bgn through $end\n";

        push @blocks, "-b $bgn -e $end";

        push @blockBgn, $bgn;
        push @blockLen, $end - $bgn + 1;
    }

    #  Each mhap job will process one block against a set of other blocks.  We'll pick, arbitrarily,
    #  to use num_blocks/4 for that size, unless it is too small.

    my $numBlocks = scalar(@blocks);
    my $qryStride = ($numBlocks < 16) ? (2) : int($numBlocks / 4);

    print STDERR "For $numBlocks blocks, set stride to $qryStride blocks.\n";

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

            #  mhap's read labeling is dumb.  It dumps all reads into one array, labeling the hash
            #  reads 0 to N, and the query reads N+1 to N+1+M.  This makes it impossible to tell if
            #  we are comparing a read against itself (e.g., when comparing hash block 1 against
            #  query blocks 1,2,3,4 -- compared to hash block 1 against query blocks 5,6,7,8).
            #
            #  So, there is a giant special case to compare the hash block against itself enabled by
            #  default.  This needs to be disabled for the second example above.
            #
            #  This makes for ANOTHER special case, that of the last block.  There are no query
            #  blocks, just the hash block, and we need to omit the flag...do we?

            my $andSelf = "";

            if ($bid == $qbgn) {
                #  The hash block bid is in the range and we need to compute hash-to-hash
                #  overlaps, and exclude the block from the queries.
                push @skipSelf, "";
                $qbgn++;
                $andSelf = " (and self)";
            } else {
                #  The hash block bid isn't in the range of query blocks, don't allow hash-to-hash
                #  overlaps.
                push @skipSelf, "--no-self";
            }

            my $qend = $qbgn + $qryStride - 1;                 #  Block bid searches reads in dat files from
            $qend = $numBlocks-1   if ($qend >= $numBlocks);   #  qbgn to qend (inclusive).


            my $job = substr("000000" . scalar(@hashes), -6);  #  Unique ID for this compute

            #  Make a place to save queries.  If this is the last-block-special-case, make a directory,
            #  but don't link in any files.  Without the directory, we'd need even more special case
            #  code down in mhap.sh to exclude the -q option for this last block.

            make_path("$path/queries/$job");

            if ($qbgn < $numBlocks) {
                print STDERR "JOB ", scalar(@hashes), " BLOCK $bid vs BLOCKS $qbgn-$qend$andSelf\n";

                for (my $qid=$qbgn; $qid <= $qend; $qid++) {
                    my $qry = substr("000000" . $qid, -6);             #  Name for the query block

                    symlink("$path/blocks/$qry.dat", "$path/queries/$job/$qry.dat");
                }

            } else {
                print STDERR "JOB ", scalar(@hashes), " BLOCK $bid vs (self)\n";
                $qbgn = $bid;  #  Otherwise, the @convert -q value is bogus
            }

            #  This is easy, the ID of the hash.

            push @hashes, substr("000000" . $bid, -6);  #  One new job for block bid with qend-qbgn query files in it

            #  Annoyingly, if we're against 'self', then the conversion needs to know that the query IDs
            #  aren't offset by the number of hash reads.

            if ($andSelf eq "") {
                push @convert, "-h $blockBgn[$bid] $blockLen[$bid] -q $blockBgn[$qbgn]";
            } else {
                push @convert, "-h $blockBgn[$bid] 0 -q $blockBgn[$bid]";
            }


            #print STDERR " -- $hashes[scalar(@hashes)-1] $convert[scalar(@convert)-1]\n";
            #print STDERR "\n";
        }
    }


    #  The ignore file is created in Meryl.pm


    #  The seed length is the shortest read such that all reads longer than this sum to 50x genome size.
    #  genomeSize must be set (ca3g should be failing early if it isn't).

    my $seedLength = 500;

    {
        my @readLengths;

        open(F, "$bin/gatekeeperDumpMetaData -reads -G $wrk/$asm.gkpStore 2> /dev/null |") or caExit("failed to get read lengths from store", undef);
        while (<F>) {
            my @v = split '\s+', $_;
            push @readLengths, $v[2];
        }
        close(F);

        @readLengths = sort { $b <=> $a } @readLengths;

        my $readLengthSum = 0;
        my $targetSum     = getGlobal("corOutCoverage") * getGlobal("genomeSize");

        foreach my $l (@readLengths) {
            $readLengthSum += $l;

            if ($readLengthSum > $targetSum) {
                $seedLength = $l;
                last;
            }
        }

        print STDERR "Computed seed length $seedLength from desired output coverage ", getGlobal("corOutCoverage"), " and genome size ", getGlobal("genomeSize"), "\n";
    }

    #  Mhap parameters - filterThreshold needs to be a string, else it is printed as 5e-06.

    my $numHashes       = (getGlobal("${tag}MhapSensitivity") eq "normal") ? "512"        : "768";
    my $minNumMatches   = (getGlobal("${tag}MhapSensitivity") eq "normal") ?   "3"        :   "2";
    my $threshold       = (getGlobal("${tag}MhapSensitivity") eq "normal") ?   "0.04"     :   "0.04";
    my $filterThreshold = (getGlobal("${tag}MhapSensitivity") eq "normal") ?   "0.000005" :   "0.000005";  #  Also set in Meryl.pm

    #  Create a script to generate precomputed blocks, including extracting the reads from gkpStore.

    my $javaPath = getGlobal("java");

    open(F, "> $path/precompute.sh") or caFailure("can't open '$path/precompute.sh' for writing: $!", undef);

    print F "#!" . getGlobal("shell") . "\n";
    print F "\n";
    print F "jobid=$taskID\n";
    print F "if [ x\$jobid = x -o x\$jobid = xundefined -o x\$jobid = x0 ]; then\n";
    print F "  jobid=\$1\n";
    print F "fi\n";
    print F "if [ x\$jobid = x ]; then\n";
    print F "  echo Error: I need $taskID set, or a job index on the command line.\n";
    print F "  exit 1\n";
    print F "fi\n";
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
    print F "if [ ! -d $path/blocks ]; then\n";
    print F "  mkdir $path/blocks\n";
    print F "fi\n";
    print F "\n";
    print F "if [ -e $path/blocks/\$job.dat ]; then\n";
    print F "  echo Job previously completed successfully.\n";
    print F "  exit\n";
    print F "fi\n";
    print F "\n";
    print F "#  If the fasta exists, our job failed, and we should try again.\n";
    print F "if [ -e \"$path/blocks/\$job.fasta\" ] ; then\n";
    print F "  rm -f $path/blocks/\$job.dat\n";
    print F "fi\n";
    print F "\n";
    print F getBinDirectoryShellCode();
    print F "\n";
    print F "\$bin/gatekeeperDumpFASTQ \\\n";
    print F "  -G $wrk/$asm.gkpStore \\\n";
    print F "  \$rge \\\n";
    print F "  -nolibname \\\n";
    print F "  -fasta \\\n";
    print F "  -o $path/blocks/\$job \\\n";
    print F "|| \\\n";
    print F "mv -f $path/blocks/\$job.fasta $path/blocks/\$job.fasta.FAILED\n";
    print F "\n";
    print F "\n";
    print F "if [ ! -e \"$path/blocks/\$job.fasta\" ] ; then\n";
    print F "  echo Failed to extract fasta.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "\n";
    print F "echo Starting mhap precompute.\n";
    print F "\n";
    print F "#  So mhap writes its output in the correct spot.\n";
    print F "cd $path/blocks\n";
    print F "\n";
    print F "$javaPath -server -Xmx", getGlobal("${tag}mhapMemory"), "g \\\n";
    print F "  -jar \$bin/mhap-1.5b1.jar \\\n";
    print F "  --weighted -k $merSize \\\n";
    print F "  --num-hashes $numHashes \\\n";
    print F "  --num-min-matches $minNumMatches \\\n";
    print F "  --threshold $threshold \\\n";
    print F "  --filter-threshold $filterThreshold \\\n";
    print F "  --min-store-length " . ($seedLength-1) . " \\\n";
    print F "  --num-threads ", getGlobal("${tag}mhapThreads"), " \\\n";
    print F "  -f $wrk/0-mercounts/$asm.ms$merSize.frequentMers.ignore \\\n"   if (-e "$wrk/0-mercounts/$asm.ms$merSize.frequentMers.ignore");
    print F "  -p $path/blocks/\$job.fasta \\\n";
    print F "  -q $path/blocks \\\n";
    print F "|| \\\n";
    print F "mv -f $path/blocks/\$job.dat $path/blocks/\$job.dat.FAILED\n";
    print F "\n";
    print F "if [ ! -e \"$path/blocks/\$job.dat\" ] ; then\n";
    print F "  echo Mhap failed.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "\n";
    print F "#  Clean up, remove the fasta input\n";
    print F "rm -f $path/blocks/\$job.fasta\n";
    print F "\n";
    print F "exit 0\n";

    close(F);

    #  Create a script to run mhap.

    open(F, "> $path/mhap.sh") or caFailure("can't open '$path/mhap.sh' for writing: $!", undef);

    print F "#!" . getGlobal("shell") . "\n";
    print F "\n";
    print F "jobid=$taskID\n";
    print F "if [ x\$jobid = x -o x\$jobid = xundefined -o x\$jobid = x0 ]; then\n";
    print F "  jobid=\$1\n";
    print F "fi\n";
    print F "if [ x\$jobid = x ]; then\n";
    print F "  echo Error: I need $taskID set, or a job index on the command line.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "\n";
    for (my $ii=1; $ii < scalar(@hashes); $ii++) {
        print F "if [ \$jobid -eq $ii ] ; then\n";
        print F "  blk=\"$hashes[$ii]\"\n";
        print F "  slf=\"$skipSelf[$ii]\"\n";
        print F "  cvt=\"$convert[$ii]\"\n";
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
    print F "if [ -e $path/blocks/\$qry.ovb.gz ]; then\n";
    print F "  echo Job previously completed successfully.\n";
    print F "  exit\n";
    print F "fi\n";
    print F "\n";
    print F "if [ ! -d $path/results ]; then\n";
    print F "  mkdir $path/results\n";
    print F "fi\n";
    print F "\n";

    print F "echo Running block \$blk in query \$qry\n";

    print F "\n";
    print F getBinDirectoryShellCode();
    print F "\n";
    print F "if [ ! -e \"$path/results/\$qry.mhap\" ] ; then\n";
    print F "  $javaPath -server -Xmx", getGlobal("${tag}mhapMemory"), "g \\\n";
    print F "    -jar \$bin/mhap-1.5b1.jar \\\n";
    print F "    --weighted -k $merSize \\\n";
    print F "    --num-hashes $numHashes \\\n";
    print F "    --num-min-matches $minNumMatches \\\n";
    print F "    --threshold $threshold \\\n";
    print F "    --filter-threshold $filterThreshold \\\n";
    print F "    --min-store-length " . ($seedLength-1) . " \\\n";
    print F "    --num-threads ", getGlobal("${tag}mhapThreads"), " \\\n";
    print F "    -f $wrk/0-mercounts/$asm.ms$merSize.frequentMers.ignore \\\n"   if (-e "$wrk/0-mercounts/$asm.ms$merSize.frequentMers.ignore");
    print F "    -s $path/blocks/\$blk.dat \$slf \\\n";
    print F "    -q $path/queries/\$qry \\\n";
    print F "  > $path/results/\$qry.mhap.WORKING \\\n";
    print F "  && \\\n";
    print F "  mv -f $path/results/\$qry.mhap.WORKING $path/results/\$qry.mhap\n";
    print F "fi\n";

    print F "\n";

    print F "if [   -e \"$path/results/\$qry.mhap\" -a \\\n";
    print F "     ! -e \"$path/results/\$qry.ovb.gz\" ] ; then\n";
    print F "  \$bin/mhapConvert \\\n";
    print F "    \$cvt \\\n";
    print F "    -o $path/results/\$qry.mhap.ovb.gz \\\n";
    print F "    $path/results/\$qry.mhap\n";
    print F "fi\n";

    print F "\n";

    if (getGlobal("${tag}mhapReAlign") eq "raw") {
        print F "if [ -e \"$path/results/\$qry.mhap.ovb.gz\" ] ; then\n";
        print F "  \$bin/overlapPair \\\n";
        print F "    -G $wrk/$asm.gkpStore \\\n";
        print F "    -O $path/results/\$qry.mhap.ovb.gz \\\n";
        print F "    -o $path/results/\$qry.ovb.gz \\\n";
        print F "    -partial \\\n"  if ($typ eq "partial");
        print F "    -erate ", getGlobal("obtOvlErrorRate"), " \\\n"  if ($typ eq "partial");
        print F "    -erate ", getGlobal("utgOvlErrorRate"), " \\\n"  if ($typ eq "normal");
        print F "    -memory 10 \\\n";
        print F "    -t 12\n";
        print F "fi\n";
    } else {
        print F "mv -f \"$path/results/\$qry.mhap.ovb.gz\" \"$path/results/\$qry.ovb.gz\"\n";
    }

    print F "\n";
    print F "\n";
    #print F "rm -rf $path/queries/\$qry\n";
    print F "\n";
    print F "exit 0\n";

    close(F);

    emitStage($WRK, $asm, "$tag-mhapConfigure");
}




sub mhapPrecomputeCheck ($$$$$) {
    my $WRK     = shift @_;
    my $wrk     = $WRK;
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $typ     = shift @_;
    my $attempt = shift @_;

    $wrk = "$wrk/correction"  if ($tag eq "cor");
    $wrk = "$wrk/trimming"    if ($tag eq "obt");

    my $path    = "$wrk/1-overlapper";

    return  if (skipStage($WRK, $asm, "$tag-mhapPrecomputeCheck", $attempt) == 1);
    return  if (-e "$path/ovljob.files");

    my $currentJobID   = 1;
    my @successJobs;
    my @failedJobs;
    my $failureMessage = "";

    open(F, "< $path/precompute.sh") or caFailure("can't open '$path/precompute.sh' for reading: $!", undef);
    while (<F>) {
        if (m/^\s+job=\"(\d+)\"$/) {
            if (-e "$path/blocks/$1.dat") {
                push @successJobs, $1;
            } else {
                $failureMessage .= "   job $path/blocks/$1.dat FAILED.\n";
                push @failedJobs, $currentJobID;
            }

            $currentJobID++;
        }
    }
    close(F);

    #  No failed jobs?  Success!

    if (scalar(@failedJobs) == 0) {
        emitStage($WRK, $asm, "$tag-mhapPrecomputeCheck");
        return;
    }

    #  If not the first attempt, report the jobs that failed, and that we're recomputing.

    if ($attempt > 1) {
        print STDERR "\n";
        print STDERR scalar(@failedJobs), " mhap precompute jobs failed:\n";
        print STDERR $failureMessage;
        print STDERR "\n";
    }

    #  If too many attempts, give up.

    if ($attempt > 2) {
        caExit("failed to precompute mhap indices.  Made " . ($attempt-1) . " attempts, jobs still failed", undef);
    }

    #  Otherwise, run some jobs.

    print STDERR "mhapPrecomputeCheck() -- attempt $attempt begins with ", scalar(@successJobs), " finished, and ", scalar(@failedJobs), " to compute.\n";

    emitStage($WRK, $asm, "$tag-mhapPrecomputeCheck", $attempt);

    submitOrRunParallelJob($WRK, $asm, "${tag}mhap", $path, "precompute", @failedJobs);
}






sub mhapCheck ($$$$$) {
    my $WRK     = shift @_;
    my $wrk     = $WRK;
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $typ     = shift @_;
    my $attempt = shift @_;

    $wrk = "$wrk/correction"  if ($tag eq "cor");
    $wrk = "$wrk/trimming"    if ($tag eq "obt");

    my $path    = "$wrk/1-overlapper";

    return  if (skipStage($WRK, $asm, "$tag-mhapCheck", $attempt) == 1);
    return  if (-e "$path/ovljob.files");

    my $currentJobID   = 1;
    my @successJobs;
    my @failedJobs;
    my $failureMessage = "";

    open(F, "< $path/mhap.sh") or caExit("failed to open '$path/mhap.sh'", undef);
    while (<F>) {
        if (m/^\s+qry=\"(\d+)\"$/) {
            if      (-e "$path/results/$1.ovb.gz") {
                push @successJobs, "$path/results/$1.ovb.gz\n";

            } elsif (-e "$path/results/$1.ovb") {
                push @successJobs, "$path/results/$1.ovb\n";

            } elsif (-e "$path/results/$1.ovb.bz2") {
                push @successJobs, "$path/results/$1.ovb.bz2\n";

            } elsif (-e "$path/results/$1.ovb.xz") {
                push @successJobs, "$path/results/$1.ovb.xz\n";

            } else {
                $failureMessage .= "   job $path/results/$1.ovb FAILED.\n";
                push @failedJobs, $currentJobID;
            }

            $currentJobID++;
        }
    }
    close(F);

    #  No failed jobs?  Success!

    if (scalar(@failedJobs) == 0) {
        open(L, "> $path/ovljob.files") or caExit("failed to open '$path/ovljob.files'", undef);
        print L @successJobs;
        close(L);
        emitStage($WRK, $asm, "$tag-mhapCheck");
        return;
    }

    #  If not the first attempt, report the jobs that failed, and that we're recomputing.

    if ($attempt > 1) {
        print STDERR "\n";
        print STDERR scalar(@failedJobs), " mhap jobs failed:\n";
        print STDERR $failureMessage;
        print STDERR "\n";
    }

    #  If too many attempts, give up.

    if ($attempt > 2) {
        caExit("failed to compute mhap overlaps.  Made " . ($attempt-1) . " attempts, jobs still failed", undef);
    }

    #  Otherwise, run some jobs.

    print STDERR "mhapCheck() -- attempt $attempt begins with ", scalar(@successJobs), " finished, and ", scalar(@failedJobs), " to compute.\n";

    emitStage($WRK, $asm, "$tag-mhapCheck", $attempt);

    submitOrRunParallelJob($WRK, $asm, "${tag}mhap", $path, "mhap", @failedJobs);
}

