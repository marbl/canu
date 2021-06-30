
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

package canu::OverlapMhap;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(mhapConfigure mhapPrecomputeCheck mhapCheck);

use strict;
use warnings "all";
no  warnings "uninitialized";
use POSIX qw(floor);

use File::Path 2.08 qw(make_path remove_tree);

use canu::Defaults;
use canu::Execution;

use canu::SequenceStore;
use canu::Report;

use canu::Grid_Cloud;

#  Map long reads to long reads with mhap.

#  Problems:
#  - mhap .dat output can't be verified.  if the job is killed, a partial .dat is output.  It needs to write to
#    .dat.WORKING, then let the script rename it when the job finishes successfully.  There is no output file name option.



sub mhapConfigure ($$$) {
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $typ     = shift @_;
    my $bin     = getBinDirectory();

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
        print STDERR "-- OVERLAPPER (mhap) (correction) complete, not rewriting scripts.\n"  if ($tag eq "cor");
        print STDERR "-- OVERLAPPER (mhap) (trimming) complete, not rewriting scripts.\n"    if ($tag eq "obt");
        print STDERR "-- OVERLAPPER (mhap) (assembly) complete, not rewriting scripts.\n"    if ($tag eq "utg");
        print STDERR "--\n";

        goto allDone;
    }

    print STDERR "--\n";
    print STDERR "-- OVERLAPPER (mhap) (correction)\n"  if ($tag eq "cor");
    print STDERR "-- OVERLAPPER (mhap) (trimming)\n"    if ($tag eq "obt");
    print STDERR "-- OVERLAPPER (mhap) (assembly)\n"    if ($tag eq "utg");
    print STDERR "--\n";

    make_path($path) if (! -d $path);

    #  Mhap parameters

    my ($numHashes, $minNumMatches, $threshold, $ordSketch, $ordSketchMer);

    if (!defined(getGlobal("${tag}MhapSensitivity"))) {
        my $cov = getExpectedCoverage($asm, $tag);

        setGlobal("${tag}MhapSensitivity", "low");                          #  Yup, super inefficient.  The code is
        setGlobal("${tag}MhapSensitivity", "normal")   if ($cov <  60);     #  compact and clear and runs once.
        setGlobal("${tag}MhapSensitivity", "high")     if ($cov <= 30);     #  Live with it.

        print STDERR "-- Set ${tag}MhapSensitivity=", getGlobal("${tag}MhapSensitivity"), " based on read coverage of $cov.\n";
    }

    if      (getGlobal("${tag}MhapSensitivity") eq "low") {
        $numHashes     =  256;
        $minNumMatches =    3;
        $threshold     =    0.80;
        $ordSketch     = 1000;
        $ordSketchMer  = getGlobal("${tag}MhapOrderedMerSize") + 2;

    } elsif (getGlobal("${tag}MhapSensitivity") eq "normal") {
        $numHashes     =  512;
        $minNumMatches =    3;
        $threshold     =    0.78;
        $ordSketch     = 1000;
        $ordSketchMer  = getGlobal("${tag}MhapOrderedMerSize") + 2;

    } elsif (getGlobal("${tag}MhapSensitivity") eq "high") {
        $numHashes     =  768;
        $minNumMatches =    2;
        $threshold     =    0.73;
        $ordSketch     = 1536;
        $ordSketchMer  = getGlobal("${tag}MhapOrderedMerSize");

    } else {
        caFailure("invalid ${tag}MhapSensitivity=" . getGlobal("${tag}MhapSensitivity"), undef);
    }

    # quick guess parameter adjustment for corrected reads, hack for now and should better take error rate into account
    if (($tag eq "obt") || ($tag eq "utg")) {
       $numHashes      = "128";
       $minNumMatches  = 5;
       $ordSketch      = 1000;
       $threshold      = 1-getGlobal("${tag}OvlErrorRate");
    }

    my $filterThreshold = getGlobal("${tag}MhapFilterThreshold");

    #  Constants.

    my $merSize       = getGlobal("${tag}MhapMerSize");

    my $numReads      = getNumberOfReadsInStore($asm, "all");   #  Need to iterate over all read IDs!
    my $memorySize    = getGlobal("${tag}mhapMemory") * 0.9;
    my $blockPerGb    = getGlobal("${tag}MhapBlockSize") / (($numHashes < 768) ? 1 : 2);

    my $javaPath      = getGlobal("java");
    my $javaOpt       = "-d64" if (defined(getGlobal("javaUse64Bit")) && getGlobal("javaUse64Bit") == 1);
    my $javaMemory    = int(getGlobal("${tag}mhapMemory") * 1024 * 0.9 + 0.5);

    print STDERR "--\n";
    print STDERR "-- PARAMETERS: hashes=$numHashes, minMatches=$minNumMatches, threshold=$threshold\n";
    print STDERR "--\n";

    my $blockSize = int($blockPerGb * $memorySize);

    print STDERR "-- Given $memorySize GB, can fit $blockSize reads per block.\n";

    #  Divide the reads into blocks of ovlHashBlockSize.  Each one of these blocks is used as the
    #  table in mhap.  Several of these blocks are used as the queries.

    my @blocks;    #  Range of reads to extract for this block
    my @blockBgn;  #  First read in the block
    my @blockLen;  #  Number of reads in the block

    my @hashes;    #  One for each job, the block that is the hash table
    my @skipSelf;  #  One for each job, jobs that would search block N against hash block N need to be handled special

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

    #  Each mhap job will process one block against a set of other blocks.  We'll pick, arbitrarily,
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
                print L "Job ", scalar(@hashes), " computes block $bid vs blocks $qbgn-$qend$andSelf,\n";

                for (my $qid=$qbgn; $qid <= $qend; $qid++) {
                    my $qry = substr("000000" . $qid, -6);             #  Name for the query block

                    symlink("../../blocks/$qry.dat", "$path/queries/$job/$qry.dat");
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

    my $cygA       = ($^O ne "cygwin") ? "" : "\$(cygpath -w ";   #  Some baloney to convert nice POSIX paths into
    my $cygB       = ($^O ne "cygwin") ? "" : ")";                #  whatever abomination Windows expects.

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
    print F fileExistsShellCode("exists", "$path", "blocks/\$job.dat");
    print F "if [ \$exists = true ] ; then\n";
    print F "  echo Job previously completed successfully.\n";
    print F "  exit\n";
    print F "fi\n";
    print F "\n";
    print F "#  Grab the kmers.ignore.gz file.\n";
    print F "mkdir -p ../0-mercounts/\n";
    print F "cd       ../0-mercounts/\n";
    print F fetchFileShellCode("$base/0-mercounts", "$asm.ms$merSize.ignore.gz", "");
    print F "cd -\n";
    print F "\n";
    print F "#  Extract the sequences we want to compute on.\n";
    print F "\$bin/sqStoreDumpFASTQ \\\n";
    print F "  -S ../../$asm.seqStore \\\n";
    print F "  \$rge \\\n";
    print F "  -nolibname \\\n";
    print F "  -noreadname \\\n";
    print F "  -fasta \\\n";
    print F "  -o ./blocks/\$job.input \\\n";
    print F "|| \\\n";
    print F "mv -f ./blocks/\$job.input.fasta ./blocks/\$job.input.fasta.FAILED\n";
    print F "\n";
    print F "if [ ! -e ./blocks/\$job.input.fasta ] ; then\n";
    print F "  echo Failed to extract fasta.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "\n";
    print F "echo \"\"\n";
    print F "echo Starting mhap precompute.\n";
    print F "echo \"\"\n";
    print F "\n";
    print F "#  Remove any previous result.\n";
    print F "rm -f ./blocks/\$job.input.dat\n";
    print F "\n";
    print F "#  So mhap writes its output in the correct spot.\n";
    print F "cd ./blocks\n";
    print F "\n";
    print F "$javaPath $javaOpt -XX:ParallelGCThreads=",  getGlobal("${tag}mhapThreads"), " -server -Xms", $javaMemory, "m -Xmx", $javaMemory, "m \\\n";
    print F "  -jar $cygA \$bin/../share/java/classes/mhap-" . getGlobal("${tag}MhapVersion") . ".jar $cygB \\\n";
    print F "  --repeat-weight 0.9 --repeat-idf-scale 10 -k $merSize \\\n";
    print F "  --supress-noise 2 \\\n"  if (defined(getGlobal("${tag}MhapFilterUnique")) && getGlobal("${tag}MhapFilterUnique") == 1);
    print F "  --no-tf \\\n"            if (defined(getGlobal("${tag}MhapNoTf")) && getGlobal("${tag}MhapNoTf") == 1);
    print F "  --store-full-id \\\n";
    print F "  --num-hashes $numHashes \\\n";
    print F "  --num-min-matches $minNumMatches \\\n";
    print F "  --ordered-sketch-size $ordSketch \\\n";
    print F "  --ordered-kmer-size $ordSketchMer \\\n";
    print F "  --threshold $threshold \\\n";
    print F "  --filter-threshold $filterThreshold \\\n";
    print F "  --min-olap-length ", getGlobal("minOverlapLength"), " \\\n";
    print F "  --num-threads ", getGlobal("${tag}mhapThreads"), " \\\n";
    print F " " . getGlobal("${tag}MhapOptions")         . " \\\n"   if (defined(getGlobal("${tag}MhapOptions")));
    print F "  -f $cygA ../../0-mercounts/$asm.ms$merSize.ignore.gz $cygB \\\n"   if (fileExists("$base/0-mercounts/$asm.ms$merSize.ignore.gz"));
    print F "  -p $cygA ./\$job.input.fasta $cygB \\\n";
    print F "  -q $cygA . $cygB \\\n";
    print F "&& \\\n";
    print F "mv -f ./\$job.input.dat ./\$job.dat\n";
    print F "\n";
    print F "if [ ! -e ./\$job.dat ] ; then\n";
    print F "  echo Mhap failed.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "\n";
    print F "#  Clean up, remove the fasta input\n";
    print F "rm -f ./\$job.input.fasta\n";
    print F "\n";
    print F stashFileShellCode("$path/blocks", "\$job.dat", "");
    print F "\n";
    print F "exit 0\n";

    close(F);

    #  Create a script to run mhap.

    open(F, "> $path/mhap.sh") or caFailure("can't open '$path/mhap.sh' for writing: $!", undef);

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

    my $stageDir = getGlobal("stageDirectory");
    my $mhapPipe = getGlobal("${tag}mhapPipe");

    print F "outPath=\"$stageDir\"\n"          if ( defined($stageDir) && ($mhapPipe == 0));
    print F "outPath=\"./results\"\n"          if (!defined($stageDir) && ($mhapPipe == 0));
    print F "\n";
    print F "if [ ! -d ./results ]; then  mkdir -p ./results;  fi\n";
    print F "if [ ! -d ./blocks  ]; then  mkdir -p ./blocks;   fi\n";
    print F "if [ ! -d \$outPath ]; then  mkdir -p \$outPath;  fi\n"   if ($mhapPipe == 0);
    print F "\n";

    print F fetchFileShellCode("$path", "blocks/\$blk.dat", "");

    print F "for ii in `ls ./queries/\$qry` ; do\n";
    print F "  echo Fetch blocks/\$ii\n";
    print F    fetchFileShellCode("$path", "blocks/\$ii", "  ");
    print F "done\n";
    print F "\n";
    print F "echo \"\"\n";
    print F "echo Running block \$blk in query \$qry\n";
    print F "echo \"\"\n";
    print F "\n";
    print F "if [ ! -e ./results/\$qry.mhap.ovb ] ; then\n";

    if ($mhapPipe) {
        print F "  #  Make a fifo so we can check return status on both\n";
        print F "  #  mhap and mhapConvert, and still pipe results so we\n";
        print F "  #  stop filling up disks.\n";
        print F "  rm -f  \$qry-pipe\n";
        print F "  mkfifo \$qry-pipe\n";
        print F "\n";
        print F "  #  Start up the consumer.\n";
        print F "  \$bin/mhapConvert \\\n";
        print F "    -S ../../$asm.seqStore \\\n";
        print F "    -o ./results/\$qry.mhap.ovb.WORKING \\\n";
        print F "    -minlength ", getGlobal("minOverlapLength"), " \\\n";
        print F "    \$qry-pipe \\\n";
        print F "  && \\\n";
        print F "  touch ./results/\$qry.mcvt.success &\n";
        print F "\n";
    }

    print F "  #  Start up the producer.\n";
    print F "  $javaPath $javaOpt -XX:ParallelGCThreads=",  getGlobal("${tag}mhapThreads"), " -server -Xms", $javaMemory, "m -Xmx", $javaMemory, "m \\\n";
    print F "    -jar $cygA \$bin/../share/java/classes/mhap-" . getGlobal("${tag}MhapVersion") . ".jar $cygB \\\n";
    print F "    --repeat-weight 0.9 --repeat-idf-scale 10 -k $merSize \\\n";
    print F "    --supress-noise 2 \\\n"  if (defined(getGlobal("${tag}MhapFilterUnique")) && getGlobal("${tag}MhapFilterUnique") == 1);
    print F "    --no-tf \\\n"            if (defined(getGlobal("${tag}MhapNoTf")) && getGlobal("${tag}MhapNoTf") == 1);
    print F "    --store-full-id \\\n";
    print F "    --num-hashes $numHashes \\\n";
    print F "    --num-min-matches $minNumMatches \\\n";
    print F "    --threshold $threshold \\\n";
    print F "    --filter-threshold $filterThreshold \\\n";
    print F "    --ordered-sketch-size $ordSketch \\\n";
    print F "    --ordered-kmer-size $ordSketchMer \\\n";
    print F "    --min-olap-length ", getGlobal("minOverlapLength"), " \\\n";
    print F "    --num-threads ", getGlobal("${tag}mhapThreads"), " \\\n";
    print F " " . getGlobal("${tag}MhapOptions")         . " \\\n"   if (defined(getGlobal("${tag}MhapOptions")));
    print F "    -s $cygA ./blocks/\$blk.dat \$slf $cygB \\\n";
    print F "    -q $cygA queries/\$qry $cygB \\\n";

    if ($mhapPipe) {
        print F "  > \$qry-pipe \\\n";
        print F "  && \\\n";
        print F "  touch ./results/\$qry.mhap.success\n";
        print F "\n";
        print F "  #  Wait for mhapConvert to finish.\n";
        print F "  wait\n";
        print F "\n";
        print F "  #  Now that they're done, check status.\n";
        print F "  if [ -e ./results/\$qry.mhap.success -a -e ./results/\$qry.mcvt.success ] ; then\n";
        print F "    mv ./results/\$qry.mhap.ovb.WORKING ./results/\$qry.mhap.ovb\n";
        print F "    rm -f ./results/\$qry.mhap.success\n";
        print F "    rm -f ./results/\$qry.mcvt.success\n";
        print F "  fi\n";
        print F "\n";
        print F "  #  And destroy the pipe.\n";
        print F "  rm -f  \$qry-pipe\n";
        print F "fi\n";
        print F "\n";
    } else {
        print F "  > \$outPath/\$qry.mhap.WORKING \\\n";
        print F "  && \\\n";
        print F "  mv -f \$outPath/\$qry.mhap.WORKING \$outPath/\$qry.mhap\n";
        print F "fi\n";
        print F "\n";
        print F "if [   -e \$outPath/\$qry.mhap -a \\\n";
        print F "     ! -e ./results/\$qry.ovb ] ; then\n";
        print F "  \$bin/mhapConvert \\\n";
        print F "    -S ../../$asm.seqStore \\\n";
        print F "    -o ./results/\$qry.mhap.ovb.WORKING \\\n";
        print F "    -minlength ", getGlobal("minOverlapLength"), " \\\n";
        print F "    \$outPath/\$qry.mhap \\\n";
        print F "  && \\\n";
        print F "  mv ./results/\$qry.mhap.ovb.WORKING ./results/\$qry.mhap.ovb\n";
        print F "fi\n";
        print F "\n";

        if (getGlobal('purgeOverlaps') ne "never") {
            print F "if [   -e \$outPath/\$qry.mhap -a \\\n";
            print F "       -e ./results/\$qry.mhap.ovb ] ; then\n";
            print F "  rm -f \$outPath/\$qry.mhap\n";
            print F "fi\n";
            print F "\n";
        }
    }

    if (getGlobal("${tag}ReAlign") eq "1") {
        print F "if [ -e ./results/\$qry.mhap.ovb ] ; then\n";
        print F "  \$bin/overlapPair \\\n";
        print F "    -S ../../$asm.seqStore \\\n";
        print F "    -O ./results/\$qry.mhap.ovb \\\n";
        print F "    -o ./results/\$qry.WORKING.ovb \\\n";
        print F "    -partial \\\n"  if ($typ eq "partial");
        print F "    -len "  , getGlobal("minOverlapLength"),  " \\\n";
        print F "    -erate ", getGlobal("corOvlErrorRate"), " \\\n"  if ($tag eq "cor");   #  Explicitly using proper name for grepability.
        print F "    -erate ", getGlobal("obtOvlErrorRate"), " \\\n"  if ($tag eq "obt");
        print F "    -erate ", getGlobal("utgOvlErrorRate"), " \\\n"  if ($tag eq "utg");
        print F "    -memory " . getGlobal("${tag}mhapMemory") . " \\\n";
        print F "    -t " . getGlobal("${tag}mhapThreads") . " \\\n";
        print F "  && \\\n";
        print F "  mv -f ./results/\$qry.WORKING.ovb ./results/\$qry.ovb\n";
        print F "fi\n";
    } else {
        print F "if [ -e ./results/\$qry.mhap.ovb ] ; then\n";
        print F "  mv -f ./results/\$qry.mhap.ovb ./results/\$qry.ovb\n";
        print F "fi\n";
    }

    print F "\n";
    print F stashFileShellCode($path, "results/\$qry.ovb", "");
    print F stashFileShellCode($path, "results/\$qry.oc",  "");
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

        print STDERR "-- Configured $numJobs mhap precompute jobs.\n";
    }

    if (-e "$path/mhap.sh") {
        my $numJobs = 0;
        open(F, "< $path/mhap.sh") or caFailure("can't open '$path/mhap.sh' for reading: $!", undef);
        while (<F>) {
            $numJobs++  if (m/^\s+qry=/);
        }
        close(F);

        print STDERR "-- Configured $numJobs mhap overlap jobs.\n";
    }

    #  Tar up the queries directory and save scripts.  The queries.tar is
    #  useful for cloud support, but we also use it to decide if this step
    #  has finished.

    runCommandSilently($path, "tar -cf queries.tar queries", 1);

    makeExecutable("$path/precompute.sh");
    makeExecutable("$path/mhap.sh");

    stashFile("$path/precompute.sh");
    stashFile("$path/mhap.sh");
    stashFile("$path/queries.tar");

  finishStage:
    generateReport($asm);
    resetIteration("$tag-mhapConfigure");

  allDone:
    stopAfter("overlapConfigure");
}




sub mhapPrecomputeCheck ($$$) {
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
    goto allDone   if (fileExists("$path/mhap.files"));
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
            if (fileExists("$path/blocks/$1.dat")) {
                push @successJobs, "1-overlapper/blocks/$1.dat\n";
            } else {
                $failureMessage .= "--   job $path/blocks/$1.dat FAILED.\n";
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
            print STDERR "-- Mhap precompute jobs failed, tried $attempt times, giving up.\n";
            print STDERR $failureMessage;
            print STDERR "--\n";
            caExit(undef, undef);
        }

        if ($attempt > 0) {
            print STDERR "--\n";
            print STDERR "-- Mhap precompute jobs failed, retry.\n";
            print STDERR $failureMessage;
            print STDERR "--\n";
        }

        #  Otherwise, run some jobs.

        generateReport($asm);

        submitOrRunParallelJob($asm, "${tag}mhap", $path, "precompute", @failedJobs);
        return;
    }

  finishStage:
    print STDERR "-- All ", scalar(@successJobs), " mhap precompute jobs finished successfully.\n";

    open(L, "> $path/precompute.files") or caExit("failed to open '$path/precompute.files'", undef);
    print L @successJobs;
    close(L);

    stashFile("$path/precompute.files");

    generateReport($asm);
    resetIteration("$tag-mhapPrecomputeCheck");

  allDone:
}






sub mhapCheck ($$$) {
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

    goto allDone   if (fileExists("$path/mhap.files"));
    goto allDone   if ((-d "$base/$asm.ovlStore") || (fileExists("$base/$asm.ovlStore.tar.gz")));

    fetchFile("$path/mhap.sh");

    #  Figure out if all the tasks finished correctly.

    my $currentJobID   = 1;
    my @mhapJobs;
    my @successJobs;
    my @miscJobs;
    my @failedJobs;
    my $failureMessage = "";

    open(F, "< $path/mhap.sh") or caExit("failed to open '$path/mhap.sh'", undef);
    while (<F>) {
        if (m/^\s+qry=\"(\d+)\"$/) {
            if      (fileExists("$path/results/$1.ovb.gz")) {
                push @mhapJobs,    "1-overlapper/results/$1.mhap\n";
                push @successJobs, "1-overlapper/results/$1.ovb.gz\n";
                push @miscJobs,    "1-overlapper/results/$1.stats\n";
                push @miscJobs,    "1-overlapper/results/$1.oc\n";

            } elsif (fileExists("$path/results/$1.ovb")) {
                push @mhapJobs,    "1-overlapper/results/$1.mhap\n";
                push @successJobs, "1-overlapper/results/$1.ovb\n";
                push @miscJobs,    "1-overlapper/results/$1.stats\n";
                push @miscJobs,    "1-overlapper/results/$1.oc\n";

            } elsif (fileExists("$path/results/$1.ovb.bz2")) {
                push @mhapJobs,    "1-overlapper/results/$1.mhap\n";
                push @successJobs, "1-overlapper/results/$1.ovb.bz2\n";
                push @miscJobs,    "1-overlapper/results/$1.stats\n";
                push @miscJobs,    "1-overlapper/results/$1.oc\n";

            } elsif (fileExists("$path/results/$1.ovb.xz")) {
                push @mhapJobs,    "1-overlapper/results/$1.mhap\n";
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
    #  (In cloud mode, the 'queries' directory probably doesn't exist when this is executed.)

    if (-e "$base/1-overlapper/queries") {
        open(F, "cd $base && find 1-overlapper/queries -print |");
        while (<F>) {
            push @mhapJobs, $_;
        }
        close(F);
    }

    #  Failed jobs, retry.

    if (scalar(@failedJobs) > 0) {

        #  If too many attempts, give up.

        if ($attempt >= getGlobal("canuIterationMax")) {
            print STDERR "--\n";
            print STDERR "-- Mhap overlap jobs failed, tried $attempt times, giving up.\n";
            print STDERR $failureMessage;
            print STDERR "--\n";
            caExit(undef, undef);
        }

        if ($attempt > 0) {
            print STDERR "--\n";
            print STDERR "-- Mhap overlap jobs failed, retry.\n";
            print STDERR $failureMessage;
            print STDERR "--\n";
        }

        #  Otherwise, run some jobs.

        generateReport($asm);

        submitOrRunParallelJob($asm, "${tag}mhap", $path, "mhap", @failedJobs);
        return;
    }

  finishStage:
    print STDERR "-- Found ", scalar(@successJobs), " mhap overlap output files.\n";

    open(L, "> $path/mhap.files") or caExit("failed to open '$path/mhap.files'", undef);
    print L @mhapJobs;
    close(L);

    open(L, "> $path/ovljob.files") or caExit("failed to open '$path/ovljob.files'", undef);
    print L @successJobs;
    close(L);

    open(L, "> $path/ovljob.more.files") or caExit("failed to open '$path/ovljob.more.files'", undef);
    print L @miscJobs;
    close(L);

    stashFile("$path/mhap.files");
    stashFile("$path/ovljob.files");
    stashFile("$path/ovljob.more.files");

    generateReport($asm);
    resetIteration("$tag-mhapCheck");

  allDone:
    stopAfter("overlap");
}

