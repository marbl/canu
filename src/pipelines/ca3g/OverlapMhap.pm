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


my $javaPath = "java";

sub mhapConfigure ($$$) {
    my $wrk          = shift @_;
    my $asm          = shift @_;
    my $bin          = getBinDirectory();
    my $type         = shift @_;

    my $path    = "$wrk/1-overlapper";

    caFailure("invalid type '$type'", undef)  if (($type ne "partial") && ($type ne "normal"));

    return  if (-e "$path/ovljob.files");
    return  if (-e "$path/mhap.sh");

    make_path("$path") if (! -d "$path");

    #  Constants.

    my $merSize       = getGlobal("ovlMerSize");
    my $ovlThreads    = getGlobal("ovlThreads");

    my $taskID        = getGlobal("gridEngineTaskID");
    my $submitTaskID  = getGlobal("gridEngineArraySubmitID");

    my $numReads      = getNumberOfReadsInStore($wrk, $asm);
    my $blockSize     = getGlobal("mhapBlockSize");

    #  Divide the reads into blocks of ovlHashBlockSize.  Each one of these blocks is used as the
    #  table in mhap.  Several of these blocks are used as the queries.

    my @blocks;
    my @blockBgn;
    my @blockLen;

    my @queries;
    my @convert;

    push @blocks,   "no zeroth block, makes the loop where this is used easier";
    push @blockBgn, "no zeroth block";
    push @blockLen, "no zeroth block";
    push @queries,  "no zeroth job";
    push @convert,  "no zeroth job";

    for (my $bgn=1; $bgn < $numReads; $bgn += $blockSize) {
        my $end = $bgn + $blockSize - 1;
        $end = $numReads  if ($end > $numReads);

        print STDERR "BLOCK ", scalar(@blocks), " reads from $bgn through $end\n";

        push @blocks, "-b $bgn -e $end";

        push @blockBgn, $bgn;
        push @blockLen, $end - $bgn + 1;
    }

    #  Each mhap job will process one block against a set of other blocks.  We'll pick, arbitrarily,
    #  to use num_blocks/4 for that size, unless it is too small.

    my $numBlocks = scalar(@blocks);
    my $qryStride = ($numBlocks < 16) ? (2) : ($numBlocks / 4);

    #  Make queries.  Each block needs to search against all blocks less than or equal to it.
    #  Each job will search at most $qryStride blocks at once.  So, queries could be:
    #  1:  1 vs 1
    #  2:  2 vs 1,2
    #  3:  3 vs 1,2,3
    #  4:  4 vs 1,2,3
    #  5:  4 vs 4
    #  6:  5 vs 1,2,3
    #  7:  5 vs 4,5

    make_path("$path/queries");

    for (my $bid=1; $bid < $numBlocks; $bid++) {
        for (my $qbgn = 1; $qbgn <= $bid; $qbgn += $qryStride) {
            my $qend = $qbgn + $qryStride - 1;  #  Block bid searches reads in dat files from
            $qend = $bid  if ($qend > $bid);    #  qbgn to qend.

            my $job = substr("000000" . scalar(@queries), -6);  #  Unique ID for this compute

            make_path("$path/queries/$job");  #  Make a place for it to work.

            for (my $qid=$qbgn; $qid <= $qend; $qid++) {
                my $qry = substr("000000" . $qid, -6);             #  Name for the query block

                print STDERR "JOB ", scalar(@queries), " BLOCK $bid vs BLOCK $qid ($path/blocks/$qry.dat -> $path/queries/$job/$qry.dat)\n";

                symlink("$path/blocks/$qry.dat", "$path/queries/$job/$qry.dat");
            }

            push @queries, substr("000000" . $bid, -6);  #  One new job for block bid with qend-qbgn query files in it
            push @convert, "-h $blockBgn[$bid] $blockLen[$bid] -q $blockBgn[$qbgn]";

            #print STDERR " -- $queries[scalar(@queries)-1] $convert[scalar(@convert)-1]\n";
            print STDERR "\n";
        }
    }

    #  Create the ignore file.
    #
    #    TTTTGTTTTTTTTTTT        0.0000044602    589     132055862
    #
    #  The fraction is just the $3/$4.  I assume this is used with "--filter-threshold 0.000005".

    if (! -e "$wrk/0-mercounts/$asm.ms$merSize.frequentMers.fasta") {
        runCommand("$wrk/0-mercounts", "meryl -Dt -n 10 -s $wrk/0-mercounts/$asm.ms$merSize > $wrk/0-mercounts/$asm.ms$merSize.frequentMers.fasta");
    }

    if (! -e "$wrk/0-mercounts/$asm.ms$merSize.frequentMers.mhap_ignore") {
        my $totalMers = 0;

        open(F, "< $wrk/0-mercounts/$asm.ms$merSize.histogram.err") or die "Failed to open '$wrk/0-mercounts/$asm.ms$merSize.histogram.err' for reading: $!\n";
        while (<F>) {
            if (m/Found\s+(\d+)\s+mers./) {
                $totalMers = $1;
            }
        }
        close(F);


        open(F, "< $wrk/0-mercounts/$asm.ms$merSize.frequentMers.fasta") or die "Failed to open '$wrk/0-mercounts/$asm.ms$merSize.frequentMers.fasta' for reading: $!\n";
        open(O, "> $wrk/0-mercounts/$asm.ms$merSize.frequentMers.mhap_ignore") or die "Failed to open '$wrk/0-mercounts/$asm.ms$merSize.frequentMers.mhap_ignore' for writing: $!\n";

        while (!eof(F)) {
            my $h = <F>;
            my $m = <F>;  chomp $m;

            if ($h =~ m/^>(\d+)/) {
                printf(O "%s\t%.16f\n", $m, $1 / $totalMers);
            }
        }

        close(O);
        close(F);
    }


    #  The seed length is the shortest read such that all reads longer than this sum to 50x genome size.

    my @readLengths;

    open(F, "$bin/gatekeeperDumpMetaData -reads -G $wrk/$asm.gkpStore |") or caFailure("failed to get read lengths from store", undef);
    while (<F>) {
        my @v = split '\s+', $_;
        push @readLengths, $v[2];
    }
    close(F);

    @readLengths = sort { $b <=> $a } @readLengths;

    my $seedLength    = 500;
    my $readLengthSum = 0;
    my $targetSum     = 50 * getGlobal("genomeSize");

    foreach my $l (@readLengths) {
        $readLengthSum += $l;

        if ($readLengthSum > $targetSum) {
            $seedLength = $l;
            last;
        }
    }

    undef @readLengths;

    print STDERR "Computed seed length $seedLength from genome size ", getGlobal("genomeSize"), "\n";



    #  Create a script to generate precomputed blocks, including extracting the reads from gkpStore.

    open(F, "> $path/precompute.sh") or caFailure("can't open '$path/precompute.sh'", undef);
    print F "#!" . getGlobal("shell") . "\n";
    print F "\n";
    print F "perl='/usr/bin/env perl'\n";
    print F "\n";
    print F "jobid=\$$taskID\n";
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
    print F "echo Starting mhap.\n";
    print F "\n";
    print F "#  So mhap writes its output in the correct spot.\n";
    print F "cd $path/blocks\n";
    print F "\n";
    print F "$javaPath -XX:+UseG1GC -server -Xmx10g \\\n";
    print F "  -jar \$bin/mhap-0.1-ob.jar FastAlignMain \\\n";
    print F "  -k $merSize \\\n";
    print F "  --num-hashes 512 \\\n";
    print F "  --num-min-matches 3 \\\n";
    print F "  --threshold 0.04 \\\n";
    print F "  --min-store-length " . ($seedLength-1) . " \\\n";
    print F "  --num-threads $ovlThreads \\\n";
    print F "  -f $wrk/$asm.ignore \\\n"      if ( -e "$wrk/$asm.ignore");
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

    #  Create a script to run mhap.

    open(F, "> $path/mhap.sh") or caFailure("can't open '$path/mhap.sh'", undef);
    print F "#!" . getGlobal("shell") . "\n";
    print F "\n";
    print F "perl='/usr/bin/env perl'\n";
    print F "\n";
    print F "jobid=\$$taskID\n";
    print F "if [ x\$jobid = x -o x\$jobid = xundefined -o x\$jobid = x0 ]; then\n";
    print F "  jobid=\$1\n";
    print F "fi\n";
    print F "if [ x\$jobid = x ]; then\n";
    print F "  echo Error: I need $taskID set, or a job index on the command line.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "\n";
    for (my $ii=1; $ii < scalar(@queries); $ii++) {
        print F "if [ \$jobid -eq $ii ] ; then\n";
        print F "  blk=\"$queries[$ii]\"\n";
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
    print F "if [ -e $path/blocks/\$qry.mhap ]; then\n";
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
    print F "$javaPath -server -Xmx10g \\\n";
    print F "  -jar \$bin/mhap-0.1-ob.jar FastAlignMain \\\n";
    print F "  --weighted -k 16 --num-hashes 512 --num-min-matches 3 --threshold 0.04 --filter-threshold 0.000005 --min-store-length 16 --num-threads 12 \\\n";
    print F "  --no-self \\\n";
    print F "  -f $wrk/$asm.ignore \\\n"      if ( -e "$wrk/$asm.ignore");
    print F "  -s $path/blocks/\$blk.dat \\\n";
    print F "  -q $path/queries/\$qry \\\n";
    print F "> $path/results/\$qry.mhap.WORKING \\\n";
    print F "&& \\\n";
    print F "mv -f $path/results/\$qry.mhap.WORKING $path/results/\$qry.mhap\n";
    print F "\n";

    print F "\n";
    print F "if [ -e \"$path/results/\$qry.mhap\" ] ; then\n";
    print F "  \$bin/mhapConvert \\\n";
    print F "    \$cvt \\\n";
    print F "    -o $path/results/\$qry.ovb.gz \\\n";
    print F "    $path/results/\$qry.mhap\n";
    print F "fi\n";
    print F "\n";
    print F "\n";
    #print F "rm -rf $path/queries/\$qry\n";
    print F "\n";
    print F "exit 0\n";
}




sub mhapPrecomputeCheck ($$$$) {
    my $wrk          = shift @_;
    my $asm          = shift @_;
    my $type         = shift @_;
    my $attempt      = shift @_;

    my $path    = "$wrk/1-overlapper";
    my $script  = "precompute";
    my $jobType = "ovl";

    return  if (-e "$path/ovljob.files");

    my $currentJobID   = 1;
    my @successJobs;
    my @failedJobs;
    my $failureMessage = "";

    open(F, "< $path/precompute.sh") or caFailure("failed to open '$path/precompute.sh'", undef);
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
        
    if (scalar(@failedJobs) == 0) {
        #open(L, "> $path/ovljob.files") or caFailure("failed to open '$path/ovljob.files'", undef);
        #print L @successJobs;
        #close(L);
        return;
    }

    if ($attempt > 0) {
        print STDERR "\n";
        print STDERR scalar(@failedJobs), " mhap precompute jobs failed:\n";
        print STDERR $failureMessage;
        print STDERR "\n";
    }

    print STDERR "mhapPrecomputeCheck() -- attempt $attempt begins with ", scalar(@successJobs), " finished, and ", scalar(@failedJobs), " to compute.\n";

    if ($attempt < 1) {
        submitOrRunParallelJob($wrk, $asm, $jobType, $path, $script, getGlobal("ovlConcurrency"), @failedJobs);
    } else {
        caFailure("failed to precompute mhap indices.  Made $attempt attempts, jobs still failed.\n", undef);
    }
}






sub mhapCheck ($$$$) {
    my $wrk          = shift @_;
    my $asm          = shift @_;
    my $type         = shift @_;
    my $attempt      = shift @_;

    my $path    = "$wrk/1-overlapper";
    my $script  = "mhap";
    my $jobType = "ovl";

    return  if (-e "$path/ovljob.files");

    my $currentJobID   = 1;
    my @successJobs;
    my @failedJobs;
    my $failureMessage = "";

    open(F, "< $path/mhap.sh") or caFailure("failed to open '$path/mhap.sh'", undef);
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
        
    if (scalar(@failedJobs) == 0) {
        open(L, "> $path/ovljob.files") or caFailure("failed to open '$path/ovljob.files'", undef);
        print L @successJobs;
        close(L);
        return;
    }

    if ($attempt > 0) {
        print STDERR "\n";
        print STDERR scalar(@failedJobs), " mhap jobs failed:\n";
        print STDERR $failureMessage;
        print STDERR "\n";
    }

    print STDERR "mhapCheck() -- attempt $attempt begins with ", scalar(@successJobs), " finished, and ", scalar(@failedJobs), " to compute.\n";

    if ($attempt < 1) {
        submitOrRunParallelJob($wrk, $asm, $jobType, $path, $script, getGlobal("ovlConcurrency"), @failedJobs);
    } else {
        caFailure("failed to compute mhap overlaps.  Made $attempt attempts, jobs still failed.\n", undef);
    }
}

