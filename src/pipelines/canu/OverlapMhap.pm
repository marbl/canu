
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
 #    src/pipelines/ca3g/OverlapMhap.pm
 #
 #  Modifications by:
 #
 #    Brian P. Walenz from 2015-MAR-27 to 2015-SEP-21
 #      are Copyright 2015 Battelle National Biodefense Institute, and
 #      are subject to the BSD 3-Clause License
 #
 #    Brian P. Walenz beginning on 2015-NOV-03
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #    Sergey Koren beginning on 2015-NOV-20
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

package canu::OverlapMhap;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(mhapConfigure mhapPrecomputeCheck mhapCheck);

use strict;
use POSIX qw(floor);

use File::Path 2.08 qw(make_path remove_tree);

use canu::Defaults;
use canu::Execution;
use canu::Gatekeeper;
use canu::HTML;
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

    goto allDone   if (skipStage($asm, "$tag-mhapConfigure") == 1);
    goto allDone   if (fileExists("$path/precompute.sh")) && (fileExists("$path/mhap.sh"));
    goto allDone   if (fileExists("$path/ovljob.files"));
    goto allDone   if (-e "$base/$asm.ovlStore");
    goto allDone   if (fileExists("$base/$asm.ovlStore.tar"));

    print STDERR "--\n";
    print STDERR "-- OVERLAPPER (mhap) (correction)\n"  if ($tag eq "cor");
    print STDERR "-- OVERLAPPER (mhap) (trimming)\n"    if ($tag eq "obt");
    print STDERR "-- OVERLAPPER (mhap) (assembly)\n"    if ($tag eq "utg");
    print STDERR "--\n";

    make_path($path) if (! -d $path);

    #  Mhap parameters

    my ($numHashes, $minNumMatches, $threshold, $ordSketch, $ordSketchMer);

    if (!defined(getGlobal("${tag}MhapSensitivity"))) {
        my $cov = getExpectedCoverage($base, $asm);

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
        $ordSketch     = 1536;
        $ordSketchMer  = getGlobal("${tag}MhapOrderedMerSize");

    } elsif (getGlobal("${tag}MhapSensitivity") eq "high") {
        $numHashes     =  768;
        $minNumMatches =    2;
        $threshold     =    0.73;
        $ordSketch     = 1536;
        $ordSketchMer  = getGlobal("${tag}MhapOrderedMerSize");

    } else {
        caFailure("invalid ${tag}MhapSensitivity=" . getGlobal("${tag}MhapSensitivity"), undef);
    }

    # due to systematic bias in nanopore data, adjust threshold up by 5%
    my $numNanoporeRaw = 0;
    open(L, "< $base/$asm.gkpStore/libraries.txt") or caExit("can't open '$base/$asm.gkpStore/libraries.txt' for reading: $!", undef);
    while (<L>) {
        $numNanoporeRaw++         if (m/nanopore-raw/);
    }
    close(L);

    $threshold += 0.05   if ($numNanoporeRaw > 0);

    my $filterThreshold = getGlobal("${tag}MhapFilterThreshold");

    #  Constants.

    my $merSize       = getGlobal("${tag}MhapMerSize");

    my $numReads      = getNumberOfReadsInStore($base, $asm);
    my $memorySize    = getGlobal("${tag}mhapMemory");
    my $blockPerGb    = getGlobal("${tag}MhapBlockSize");
    if ($numHashes >= 768) {
       $blockPerGb = int($blockPerGb / 2);
    }

    # quick guess parameter adjustment for corrected reads, hack for now and should better take error rate into account
    if (($tag eq "obt") || ($tag eq "utg")) {
       $numHashes      = "128";
       $minNumMatches  = 5;
       $ordSketch      = 1000;
       $threshold      = 1-getGlobal("${tag}OvlErrorRate");
    }

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

            #  Annoyingly, if we're against 'self', then the conversion needs to know that the query IDs
            #  aren't offset by the number of hash reads.

            if ($andSelf eq "") {
                push @convert, "-h $blockBgn[$bid] $blockLen[$bid] -q $blockBgn[$qbgn]";
            } else {
                push @convert, "-h $blockBgn[$bid] 0 -q $blockBgn[$bid]";
            }
        }
    }

    close(L);

    #  Tar up the queries directory.  Only useful for cloud support.

    runCommandSilently($path, "tar -cf queries.tar queries", 1);
    stashFile("$path/queries.tar");

    #  Create a script to generate precomputed blocks, including extracting the reads from gkpStore.

    #OPTIMIZE
    #OPTIMIZE  Probably a big optimization for cloud assemblies, the block fasta inputs can be
    #OPTIMIZE  computed ahead of time, stashed, and then fetched to do the actual precompute.
    #OPTIMIZE

    my $javaPath   = getGlobal("java");
    my $javaMemory = int(getGlobal("${tag}mhapMemory") * 1024 + 0.5);

    my $cygA       = ($^O ne "cygwin") ? "" : "\$(cygpath -w ";   #  Some baloney to convert nice POSIX paths into
    my $cygB       = ($^O ne "cygwin") ? "" : ")";                #  whatever abomination Windows expects.

    open(F, "> $path/precompute.sh") or caFailure("can't open '$path/precompute.sh' for writing: $!", undef);

    print F "#!" . getGlobal("shell") . "\n";
    print F "\n";
    print F getBinDirectoryShellCode();
    print F "\n";
    print F setWorkDirectoryShellCode($path);
    print F fetchStoreShellCode("$base/$asm.gkpStore", "$base/1-overlapper", "");
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
    print F fileExistsShellCode("./blocks/\$job.dat");
    print F "  echo Job previously completed successfully.\n";
    print F "  exit\n";
    print F "fi\n";
    print F "\n";
    print F "#  Remove any previous result.\n";
    print F "rm -f ./blocks/\$job.input.dat\n";
    print F "\n";
    #print F fetchFileShellCode("./blocks/\$job.input.fasta");
    print F "\n";
    print F "\$bin/gatekeeperDumpFASTQ \\\n";
    print F "  -G ../$asm.gkpStore \\\n";
    print F "  \$rge \\\n";
    print F "  -nolibname \\\n";
    print F "  -fasta \\\n";
    print F "  -o ./blocks/\$job.input \\\n";
    print F "|| \\\n";
    print F "mv -f ./blocks/\$job.input.fasta ./blocks/\$job.input.fasta.FAILED\n";
    print F "\n";
    print F "\n";
    print F "if [ ! -e ./blocks/\$job.input.fasta ] ; then\n";
    print F "  echo Failed to extract fasta.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "\n";
    print F fetchFileShellCode("$base/0-mercounts", "$asm.ms$merSize.frequentMers.ignore.gz", "");
    print F "\n";
    print F "echo \"\"\n";
    print F "echo Starting mhap precompute.\n";
    print F "echo \"\"\n";
    print F "\n";
    print F "#  So mhap writes its output in the correct spot.\n";
    print F "cd ./blocks\n";
    print F "\n";
    print F "$javaPath -d64 -server -Xmx", $javaMemory, "m \\\n";
    print F "  -jar $cygA \$bin/mhap-" . getGlobal("${tag}MhapVersion") . ".jar $cygB \\\n";
    print F "  --repeat-weight 0.9 --repeat-idf-scale 10 -k $merSize \\\n";
    print F "  --supress-noise 2 \\\n"  if (defined(getGlobal("${tag}MhapFilterUnique")) && getGlobal("${tag}MhapFilterUnique") == 1);
    print F "  --no-tf \\\n"            if (defined(getGlobal("${tag}MhapNoTf")) && getGlobal("${tag}MhapNoTf") == 1);
    print F "  --num-hashes $numHashes \\\n";
    print F "  --num-min-matches $minNumMatches \\\n";
    print F "  --ordered-sketch-size $ordSketch \\\n";
    print F "  --ordered-kmer-size $ordSketchMer \\\n";
    print F "  --threshold $threshold \\\n";
    print F "  --filter-threshold $filterThreshold \\\n";
    print F "  --min-olap-length ", getGlobal("minOverlapLength"), " \\\n";
    print F "  --num-threads ", getGlobal("${tag}mhapThreads"), " \\\n";
    print F " " . getGlobal("${tag}MhapOptions")         . " \\\n"   if (defined(getGlobal("${tag}MhapOptions")));
    print F "  -f $cygA ../../0-mercounts/$asm.ms$merSize.frequentMers.ignore.gz $cygB \\\n"   if (-e "$base/0-mercounts/$asm.ms$merSize.frequentMers.ignore.gz");
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
    print F stashFileShellCode("$base/1-overlapper/blocks", "\$job.dat", "");
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
    print F fetchStoreShellCode("$base/$asm.gkpStore", "$base/1-overlapper", "");
    print F "\n";
    print F getJobIDShellCode();
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
    print F "if [ ! -d ./blocks ] ; then\n";
    print F "  mkdir -p ./blocks\n";
    print F "fi\n";

    print F fetchFileShellCode("$path", "blocks/\$blk.dat", "");

    print F "for ii in `ls ./queries/\$qry` ; do\n";
    print F "  echo Fetch blocks/\$ii\n";
    print F    fetchFileShellCode("$path", "blocks/\$ii", "  ");
    print F "done\n";
    print F "\n";

    print F fetchFileShellCode("$base/0-mercounts", "$asm.ms$merSize.frequentMers.ignore.gz", "");
    print F "\n";

    print F "echo \"\"\n";
    print F "echo Running block \$blk in query \$qry\n";
    print F "echo \"\"\n";
    print F "\n";
    print F "if [ ! -e ./results/\$qry.mhap ] ; then\n";
    print F "  $javaPath -d64 -server -Xmx", $javaMemory, "m \\\n";
    print F "    -jar $cygA \$bin/mhap-" . getGlobal("${tag}MhapVersion") . ".jar $cygB \\\n";
    print F "    --repeat-weight 0.9 --repeat-idf-scale 10 -k $merSize \\\n";
    print F "    --supress-noise 2 \\\n"  if (defined(getGlobal("${tag}MhapFilterUnique")) && getGlobal("${tag}MhapFilterUnique") == 1);
    print F "    --no-tf \\\n"            if (defined(getGlobal("${tag}MhapNoTf")) && getGlobal("${tag}MhapNoTf") == 1);
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
    print F "  > ./results/\$qry.mhap.WORKING \\\n";
    print F "  && \\\n";
    print F "  mv -f ./results/\$qry.mhap.WORKING ./results/\$qry.mhap\n";
    print F "fi\n";
    print F "\n";

    print F "if [   -e ./results/\$qry.mhap -a \\\n";
    print F "     ! -e ./results/\$qry.ovb ] ; then\n";
    print F "  \$bin/mhapConvert \\\n";
    print F "    -G ../$asm.gkpStore \\\n";
    print F "    \$cvt \\\n";
    print F "    -o ./results/\$qry.mhap.ovb.WORKING \\\n";
    print F "    ./results/\$qry.mhap \\\n";
    print F "  && \\\n";
    print F "  mv ./results/\$qry.mhap.ovb.WORKING ./results/\$qry.mhap.ovb\n";
    print F "fi\n";
    print F "\n";

    if (getGlobal('saveOverlaps') eq "0") {
        print F "if [   -e ./results/\$qry.mhap -a \\\n";
        print F "       -e ./results/\$qry.mhap.ovb ] ; then\n";
        print F "  rm -f ./results/\$qry.mhap\n";
        print F "fi\n";
        print F "\n";
    }

    if (getGlobal("${tag}ReAlign") eq "1") {
        print F "if [ -e ./results/\$qry.mhap.ovb ] ; then\n";
        print F "  \$bin/overlapPair \\\n";
        print F "    -G ../$asm.gkpStore \\\n";
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
    print F stashFileShellCode("$path", "results/\$qry.ovb",    "");
    print F stashFileShellCode("$path", "results/\$qry.counts", "");
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

    makeExecutable("$path/precompute.sh");
    makeExecutable("$path/mhap.sh");

    stashFile("$path/precompute.sh");
    stashFile("$path/mhap.sh");

  finishStage:
    emitStage($asm, "$tag-mhapConfigure");
    buildHTML($asm, $tag);

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

    goto allDone   if (skipStage($asm, "$tag-mhapPrecomputeCheck", $attempt) == 1);
    goto allDone   if (fileExists("$path/precompute.files"));
    goto allDone   if (-e "$base/$asm.ovlStore");
    goto allDone   if (fileExists("$base/$asm.ovlStore.tar"));

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
                $failureMessage .= "--   job 1-overlapper/blocks/$1.dat FAILED.\n";
                push @failedJobs, $currentJobID;
            }

            $currentJobID++;
        }
    }
    close(F);

    #  Failed jobs, retry.

    if (scalar(@failedJobs) > 0) {

        #  If not the first attempt, report the jobs that failed, and that we're recomputing.

        if ($attempt > 1) {
            print STDERR "--\n";
            print STDERR "-- ", scalar(@failedJobs), " mhap precompute jobs failed:\n";
            print STDERR $failureMessage;
            print STDERR "--\n";
        }

        #  If too many attempts, give up.

        if ($attempt > getGlobal("canuIterationMax")) {
            caExit("failed to precompute mhap indices.  Made " . ($attempt-1) . " attempts, jobs still failed", undef);
        }

        #  Otherwise, run some jobs.

        emitStage($asm, "$tag-mhapPrecomputeCheck", $attempt);
        buildHTML($asm, $tag);

        submitOrRunParallelJob($asm, "${tag}mhap", $path, "precompute", @failedJobs);
        return;
    }

  finishStage:
    print STDERR "-- All ", scalar(@successJobs), " mhap precompute jobs finished successfully.\n";

    open(L, "> $path/precompute.files") or caExit("failed to open '$path/precompute.files'", undef);
    print L @successJobs;
    close(L);

    stashFile("$path/precompute.files");

    emitStage($asm, "$tag-mhapPrecomputeCheck");
    buildHTML($asm, $tag);

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

    goto allDone   if (skipStage($asm, "$tag-mhapCheck", $attempt) == 1);
    goto allDone   if (fileExists("$path/mhap.files"));
    goto allDone   if (-e "$base/$asm.ovlStore");
    goto allDone   if (fileExists("$base/$asm.ovlStore.tar"));

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
                push @miscJobs,    "1-overlapper/results/$1.counts\n";

            } elsif (fileExists("$path/results/$1.ovb")) {
                push @mhapJobs,    "1-overlapper/results/$1.mhap\n";
                push @successJobs, "1-overlapper/results/$1.ovb\n";
                push @miscJobs,    "1-overlapper/results/$1.stats\n";
                push @miscJobs,    "1-overlapper/results/$1.counts\n";

            } elsif (fileExists("$path/results/$1.ovb.bz2")) {
                push @mhapJobs,    "1-overlapper/results/$1.mhap\n";
                push @successJobs, "1-overlapper/results/$1.ovb.bz2\n";
                push @miscJobs,    "1-overlapper/results/$1.stats\n";
                push @miscJobs,    "1-overlapper/results/$1.counts\n";

            } elsif (fileExists("$path/results/$1.ovb.xz")) {
                push @mhapJobs,    "1-overlapper/results/$1.mhap\n";
                push @successJobs, "1-overlapper/results/$1.ovb.xz\n";
                push @miscJobs,    "1-overlapper/results/$1.stats\n";
                push @miscJobs,    "1-overlapper/results/$1.counts\n";

            } else {
                $failureMessage .= "--   job 1-overlapper/results/$1.ovb FAILED.\n";
                push @failedJobs, $currentJobID;
            }

            $currentJobID++;
        }
    }
    close(F);

    #  Also find the queries symlinks so we can remove those.  And the query directories, because
    #  the last directory can be empty, and so we'd never see it at all if only finding files.

    open(F, "find $path/queries -print |");
    while (<F>) {
        push @mhapJobs, $_;
    }
    close(F);

    #  Failed jobs, retry.

    if (scalar(@failedJobs) > 0) {

        #  If not the first attempt, report the jobs that failed, and that we're recomputing.

        if ($attempt > 1) {
            print STDERR "--\n";
            print STDERR "-- ", scalar(@failedJobs), " mhap jobs failed:\n";
            print STDERR $failureMessage;
            print STDERR "--\n";
        }

        #  If too many attempts, give up.

        if ($attempt > getGlobal("canuIterationMax")) {
            caExit("failed to compute mhap overlaps.  Made " . ($attempt-1) . " attempts, jobs still failed", undef);
        }

        #  Otherwise, run some jobs.

        emitStage($asm, "$tag-mhapCheck", $attempt);
        buildHTML($asm, $tag);
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

    emitStage($asm, "$tag-mhapCheck");
    buildHTML($asm, $tag);

  allDone:
    stopAfter("overlap");
}

