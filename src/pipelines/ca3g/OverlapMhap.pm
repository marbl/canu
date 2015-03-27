package ca3g::OverlapMhap;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(mhapConfigure mhapCheck);

use strict;

use File::Path qw(make_path remove_tree);

use ca3g::Defaults;
use ca3g::Execution;
use ca3g::Gatekeeper;

#  Map long reads to long reads with mhap.

my $javaPath = "java";

sub mhapConfigure ($$$) {
    my $wrk          = shift @_;
    my $asm          = shift @_;
    my $type         = shift @_;

    caFailure("invalid type '$type'", undef)  if (($type ne "partial") && ($type ne "normal"));

    #return  if  (-d "$wrk/$asm.tigStore");

    #return  if ((-d "$wrk/$asm.obtStore") && ($type eq "partial"));
    #return  if ((-d "$wrk/$asm.ovlStore") && ($type eq "normal"));

    #return  if (-e "$wrk/1-overlapper/ovlopt");

    make_path("$wrk/1-overlapper") if (! -d "$wrk/1-overlapper");

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

    push @blocks, "no first block, makes the loop where this is used easier";

    for (my $bgn=1; $bgn < $numReads; $bgn += $blockSize) {
        my $end = $bgn + $blockSize - 1;
        $end = $numReads  if ($end > $numReads);

        push @blocks, "-b $bgn -e $end";
    }

    my $numBlocks = scalar(@blocks);
    my $qryStride = 0;

    #  Each mhap job will process one block against a set of other blocks.  We'll pick, arbitrarily,
    #  to use num_blocks/4 for that size, unless it is too small.

    if      ($numBlocks < 16) {
        $qryStride = 4;
    } else {
        $qryStride = $numBlocks / 4;
    }

    #  For each block, make a work directory of inputs to search.  This is just symlinks to other blocks.

    make_path("$wrk/1-overlapper/queries");

    my $numJobs = 1;

    for (my $bid=0; $bid < $numBlocks; $bid++) {
        for (my $qbgn = 1; $qbgn < $numBlocks; $qbgn += $qryStride) {
            my $qend = $qbgn + $qryStride;
            $qend = $numBlocks  if ($qend > $numBlocks);

            my $jobDir  = substr("000000" . $numJobs, -6);
            my $localID = "0";

            make_path("$wrk/1-overlapper/queries/$jobDir");

            symlink("$wrk/1-overlapper/\$bat/\$job.dat",
                    "$wrk/1-overlapper/queries/$jobDir/");

            $numJobs++;
            $localID++;
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

    print STDERR "NOT REALLY COMPUTING SEED LENGTH!\n";
    my $seedLength = 500;

    #  Create a script to generate precomputed blocks, including extracting the reads from gkpStore.

    open(F, "> $wrk/1-overlapper/precompute.sh") or caFailure("can't open '$wrk/1-overlapper/precompute.sh'", undef);
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
        print F "  bat=\"001\"\n";
        print F "  job=`printf %06d \$ii`\n";
        print F "fi\n";
        print F "\n";
    }
    print F "\n";
    print F "if [ x\$bat = x ]; then\n";
    print F "  echo Error: Job index out of range.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "\n";
    print F "if [ ! -d $wrk/1-overlapper/\$bat ]; then\n";
    print F "  mkdir $wrk/1-overlapper/\$bat\n";
    print F "fi\n";
    print F "\n";
    print F "if [ -e $wrk/1-overlapper/\$bat/\$job.ovb.gz ]; then\n";
    print F "  echo Job previously completed successfully.\n";
    print F "  exit\n";
    print F "fi\n";
    print F "\n";
    print F getBinDirectoryShellCode();
    print F "\n";
    print F "\$bin/gatekeeperDumpFASTQ -G $wrk/$asm.gkpStore \$rge -nolibname -fasta -o $wrk/1-overlapper/\$bat/\$job \n";
    print F "\n";
    print F "echo Starting mhap.\n";
    print F "\n";
    print F "#  So mhap writes its output in the correct spot.\n";
    print F "cd $wrk/1-overlapper/\$bat\n";
    print F "\n";
    print F "$javaPath -XX:+UseG1GC -server -Xmx10g \\\n";
    print F "  -jar \$bin/java/mhap-0.1-ob.jar FastAlignMain \\\n";
    print F "  -k $merSize \\\n";
    print F "  --num-hashes 512 \\\n";
    print F "  --num-min-matches 3 \\\n";
    print F "  --threshold 0.04 \\\n";
    print F "  --min-store-length " . ($seedLength-1) . " \\\n";
    print F "  --num-threads $ovlThreads \\\n";
    print F "  -f $wrk/$asm.ignore \\\n"      if ( -e "$wrk/$asm.ignore");
    print F "  -p $wrk/1-overlapper/\$bat/\$job.fasta \\\n";
    print F "  -q $wrk/1-overlapper/ \\\n";
    print F ">    $wrk/1-overlapper/\$bat/\$job.hash.err 2>&1\n";
    print F "\n";
    print F "#  Clean up, remove the fasta input\n";
    print F "rm -f $wrk/1-overlapper/\$bat/\$job.fasta\n";
    print F "\n";
    print F "exit 0\n";

    #  Create a script to run mhap.

    open(F, "> $wrk/1-overlapper/mhap.sh") or caFailure("can't open '$wrk/1-overlapper/mhap.sh'", undef);
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
    for (my $ii=0; $ii<scalar(@blocks); $ii++) {
        print F "if [ \$jobid -eq 1 ] ; then\n";
        print F "  hashrange=\"$blocks[$ii]\"\n";
        print F "  bat=\"001\"\n";
        print F "  job=`printf %06d \$ii`\n";
        print F "fi\n";
        print F "\n";
        print F "\n";
    }

    print F "\n";
    print F "if [ ! -d $wrk/1-overlapper/\$bat ]; then\n";
    print F "  mkdir $wrk/1-overlapper/\$bat\n";
    print F "fi\n";
    print F "\n";
    print F "if [ -e $wrk/1-overlapper/\$bat/\$job.ovb ]; then\n";
    print F "  echo Job previously completed successfully.\n";
    print F "  exit\n";
    print F "fi\n";
    print F "\n";
    print F "if [ x\$bat = x ]; then\n";
    print F "  echo Error: Job index out of range.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "\n";
    print F "options=(\$opt)\n";
    
    #  screws up indenting
    #print F "hashParams=(\$(echo \${options[1]}  | tr \"-\" \"\\n\")) \n";
    #print F "hashStart=\${hashParams[0]} \n";
    #print F "hashStart=\$((\$hashStart - 1))\n";
    #print F "hashEnd=\${hashParams[1]} \n";
    #print F "params=(\$(echo \${options[3]}  | tr \"-\" \"\\n\")) \n";

    #print F "numSA=" , ($numIndicies-1) , "\n";


    #print F "start=\${params[0]} \n";
    #print F "start=\$((\$start - 1))\n";

    #print F "end=\${params[1]} \n";

    #print F "total=\$((\$end - \$start))\n"; 

    #print F "refOffset=\$((\$start - \$hashEnd))\n";

    #print F "echo \"Running partition \$job with options \$opt start \$start end \$end total \$total  and \"\n";

    #print F "\n";

    #print F " rm -rf $wrk/1-overlapper/stream_\$jobid\n";
    #print F " mkdir -p $wrk/1-overlapper/stream_\$jobid\n";


    #print F " startIndex=`perl -w -e \"use POSIX; print floor(\$start/$ovlRefBlockSize)\"`\n";
    #print F " endIndex=`perl -w -e \"use POSIX; print ceil(\$end/$ovlRefBlockSize)\"`\n";


    #print F " for i in \$(seq \$((\$startIndex+1)) \$endIndex); do\n";
    #print F "    leadingI=`printf %06d \$i`\n";
    #print F "    ln -s $wrk/1-overlapper/correct_reads_part\$i.dat $wrk/1-overlapper/stream_\$jobid/correct_reads_part\$leadingI.dat\n";
    #print F " done\n";

    #  Comparing block X vs blocks Y-Z.
    #
    #  If the last block, there isn't another block to compare against.

    #print F " if [ \$sa -eq \$numSA ]; then\n";
    #print F "    $javaCmd              -s $wrk/1-overlapper/correct_reads_part\$sa.dat 2> $wrk/1-overlapper/\$jobid.err | $convertCmd | \$bin/convertOverlap -ovl > $wrk/1-overlapper/\$bat/\$job.ovb\n";
    #print F " else\n";
    #print F "    if [ \$start -eq \$hashEnd ]; then\n";
    #print F "       $javaCmd           -s $wrk/1-overlapper/correct_reads_part\$sa.dat -q $wrk/1-overlapper/stream_\$jobid 2> $wrk/1-overlapper/\$jobid.err | $convertCmd | \$bin/convertOverlap -ovl > $wrk/1-overlapper/\$bat/\$job.ovb\n";
    #print F "    else\n";
    #print F "       $javaCmd --no-self -s $wrk/1-overlapper/correct_reads_part\$sa.dat -q $wrk/1-overlapper/stream_\$jobid 2> $wrk/1-overlapper/\$jobid.err | $convertCmd | \n";
    #print F "    fi\n";
    #print F " fi\n";



    print F "$javaPath -server -Xmx10g \\\n";
    print F "  -jar \$bin/java/mhap-0.1-ob.jar FastAlignMain \\\n";
    print F "  --weighted -k 16 --num-hashes 512 --num-min-matches 3 --threshold 0.04 --filter-threshold 0.000005 --min-store-length 16 --num-threads 12 \\\n";
    print F "  --no-self \\\n";
    print F "  -f $wrk/$asm.ignore \\\n"      if ( -e "$wrk/$asm.ignore");
    print F "  -s $wrk/1-overlapper/correct_reads_part\$sa.dat \\\n";
    print F "  -q $wrk/1-overlapper/stream_\$jobid \\\n";   #  ONLY if $sa -ne $numSA ???
    print F "| \\\n";

    #  hash start is the first read id in the hash table, minus one
    #  ref offset is 'first_read_in_reference - last_read_in_hash_table'

    #my $convertCmd = "awk -v REF_OFFSET=\$refOffset -v OFFSET=\$hashStart '{if (\$5 == 0 && \$9 == 0) { ORI=\"N\"; } if (\$5 == 0 && \$9 == 1) { ORI=\"I\"; } if (\$8 <= \$12 && (\$7-\$6) / \$8 > 0.9) print \$1+OFFSET+REF_OFFSET\"\\t\"\$2+OFFSET\"\\t\"ORI\"\\t\"(-1*\$10)\"\\t\"(\$12-\$11)\"\\t\"\$3/5\"\\t\"\$3/5; else if (\$8 > \$12 && (\$11 - \$10) / \$12 > 0.9) print \$1+OFFSET+REF_OFFSET\"\\t\"\$2+OFFSET\"\\t\"ORI\"\\t\"\$6\"\\t\"(-1*(\$8-\$7))\"\\t\"\$3/5\"\\t\"\$3/5; }'";

    print F "| \\\n";
    print F "\$bin/convertOverlap -ovl > $wrk/1-overlapper/\$bat/\$job.ovb\n";


    print F "\n";
    print F "rm -rf $wrk/1-overlapper/stream_\$jobid\n";

    print F "\n";
    print F "exit 0\n";
}




sub mhapCheck ($$$$) {
    my $wrk          = shift @_;
    my $asm          = shift @_;
    my $type         = shift @_;
    my $attempt      = shift @_;

    my $path    = "$wrk/1-overlapper";
    my $script  = "overlap";
    my $jobType = "ovl";

    return  if ((-d "$wrk/$asm.obtStore") && ($type eq "partial"));
    return  if ((-d "$wrk/$asm.ovlStore") && ($type eq "normal"));
    return  if  (-d "$wrk/$asm.tigStore");

}

