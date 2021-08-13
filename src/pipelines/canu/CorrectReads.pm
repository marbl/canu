
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

package canu::CorrectReads;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(setupCorrectionParameters buildCorrectionLayoutsConfigure buildCorrectionLayoutsCheck filterCorrectionLayouts generateCorrectedReadsConfigure generateCorrectedReadsCheck loadCorrectedReads dumpCorrectedReads);

use strict;
use warnings "all";
no  warnings "uninitialized";

use File::Path 2.08 qw(make_path remove_tree);

use canu::Defaults;
use canu::Execution;

use canu::Configure;
use canu::SequenceStore;
use canu::Report;
use canu::Output;

use canu::Grid_Cloud;


#  Returns a coverage:
#    If $cov not defined, default to desired output coverage * 1.0.
#    Otherwise, if defined but ends in an 'x', that's desired output coverage * whatever
#    Otherwise, the coverage is as defined.
#
sub getCorCov ($$) {
    my $asm     = shift @_;
    my $typ     = shift @_;
    my $cov     = getGlobal("corMaxEvidenceCoverage$typ");

    my $exp = getExpectedCoverage($asm, "cor");
    my $des = getGlobal("corOutCoverage");

    if (!defined($cov)) {
        $cov = $des;
    } elsif ($cov =~ m/(.*)x/) {
        $cov = int($des * $1);
    }

    return($cov);
}



sub setupCorrectionParameters ($) {
    my $asm     = shift @_;
    my $bin     = getBinDirectory();
    my $cmd;

    my $base    = "correction";
    my $path    = "correction/2-correction";

    make_path("$path")  if (! -d "$path");

    #  Set the minimum coverage for a corrected read based on coverage in input reads.

    if (!defined(getGlobal("corMinCoverage"))) {
        my $cov = getExpectedCoverage($asm, "cor");

        setGlobal("corMinCoverage", 4);
        setGlobal("corMinCoverage", 4)   if ($cov <  60);
        setGlobal("corMinCoverage", 0)   if ($cov <= 20);

        print STDERR "-- Set corMinCoverage=", getGlobal("corMinCoverage"), " based on read coverage of $cov.\n";
    }
}



sub buildCorrectionLayoutsConfigure ($) {
    my $asm     = shift @_;
    my $bin     = getBinDirectory();
    my $cmd;

    my $base    = "correction";
    my $path    = "correction/2-correction";

    goto allDone       if (-d "$base/$asm.corStore");                             #  Jobs all finished
    goto allDone       if (fileExists("$base/$asm.corStore/seqDB.v001.dat"));
    goto allDone       if (fileExists("$base/$asm.corStore/seqDB.v001.tig"));

    goto finishStage   if ((-e "$base/$asm.corStore.WORKING/seqDB.v001.dat") &&   #  Job ran manually, showNext
                           (-e "$base/$asm.corStore.WORKING/seqDB.v001.tig"));

    #  Make layouts for each corrected read.
    #
    #  The global filter is estimated from data saved in ovlStore.  To
    #  compute it exactly, use filterCorrectionOverlaps, with parameter
    #  -scores <outfile>.globalScores, and pass that to
    #  generateCorrectionLayouts using the -scores parameter.
    #
    #  This was the original way to apply the global filter, but once scores
    #  were saved in the overlap store, it has not been tested.
    #
    #  (the pipeline code for this was removed somewhere after
    #  81fe8a145ee9006d35503a9e8fba259cbde88a51)

    fetchOvlStore($asm, $base);

    print STDERR "-- Computing correction layouts.\n";
    print STDERR "--   Local  filter coverage   ", getCorCov($asm, "Local"),  "\n";
    print STDERR "--   Global filter coverage   ", getCorCov($asm, "Global"), "\n";

    $cmd  = "$bin/generateCorrectionLayouts \\\n";
    $cmd .= "  -S ../$asm.seqStore \\\n";
    $cmd .= "  -O  ./$asm.ovlStore \\\n";
    $cmd .= "  -C  ./$asm.corStore.WORKING \\\n";
    $cmd .= "  -eL " . getGlobal("corMinEvidenceLength") . " \\\n"  if (defined(getGlobal("corMinEvidenceLength")));
    $cmd .= "  -eE " . getGlobal("corMaxEvidenceErate")  . " \\\n"  if (defined(getGlobal("corMaxEvidenceErate")));
    $cmd .= "  -eC " . getCorCov($asm, "Local") . " \\\n";
    $cmd .= "  -xC " . getCorCov($asm, "Global") . " \\\n";
    $cmd .= "> ./$asm.corStore.err 2>&1";

    if (runCommand($base, $cmd)) {
        caExit("failed to generate correction layouts", "$base/$asm.corStore.err");
    }

    unlink "$base/$asm.corStore.err";

  finishStage:
    rename "$base/$asm.corStore.WORKING", "$base/$asm.corStore";

    stashFile("$base/$asm.corStore/seqDB.v001.dat");
    stashFile("$base/$asm.corStore/seqDB.v001.tig");

    generateReport($asm);
    resetIteration("cor-buildCorrectionLayoutsConfigure");

  allDone:
}



sub buildCorrectionLayoutsCheck ($) {
    my $asm     = shift @_;
    my $bin     = getBinDirectory();
    my $cmd;

    my $base    = "correction";
    my $path    = "correction/2-correction";

    goto allDone   if (fileExists("$base/$asm.corStore"));                      #  Jobs all finished

    #  Eventually, we'll run generateCorrectionLayouts on the grid.  Then we'll need to load the new
    #  tigs into the corStore here.

  finishStage:
    generateReport($asm);
    resetIteration("cor-buildCorrectionLayoutsCheck");

  allDone:
}



sub filterCorrectionLayouts ($) {
    my $asm     = shift @_;
    my $bin     = getBinDirectory();
    my $cmd;

    my $base    = "correction";
    my $path    = "correction/2-correction";

    goto allDone       if (fileExists("$path/$asm.readsToCorrect"));                #  Jobs all finished

    goto finishStage   if ((-e "$path/$asm.readsToCorrect.WORKING.stats") &&        #  Job ran manually, showNext
                           (-e "$path/$asm.readsToCorrect.WORKING.log"));

    #  Analyze the corStore to decide what reads we want to correct.

    fetchOvlStore($asm, $base);

    print STDERR "-- Computing correction layouts.\n";

    $cmd  = "$bin/filterCorrectionLayouts \\\n";
    $cmd .= "  -S  ../../$asm.seqStore \\\n";
    $cmd .= "  -C     ../$asm.corStore \\\n";
    $cmd .= "  -R      ./$asm.readsToCorrect.WORKING \\\n";
    $cmd .= "  -cc " . getGlobal("corMinCoverage") . " \\\n";
    $cmd .= "  -cl " . getGlobal("minReadLength")  . " \\\n";
    $cmd .= "  -g  " . getGlobal("genomeSize")     . " \\\n";
    $cmd .= "  -c  " . getGlobal("corOutCoverage") . " \\\n";
    $cmd .= "> ./$asm.readsToCorrect.err 2>&1";

    if (runCommand($path, $cmd)) {
        caExit("failed to generate list of reads to correct", "$path/$asm.readsToCorrect.err");
    }

  finishStage:
    rename "$path/$asm.readsToCorrect.WORKING",       "$path/$asm.readsToCorrect";
    rename "$path/$asm.readsToCorrect.WORKING.stats", "$path/$asm.readsToCorrect.stats";
    rename "$path/$asm.readsToCorrect.WORKING.log",   "$path/$asm.readsToCorrect.log";

    stashFile("$path/$asm.readsToCorrect");
    stashFile("$path/$asm.readsToCorrect.stats");
    stashFile("$path/$asm.readsToCorrect.log");

    my $report = getFromReport("corLayout");

    open(F, "< $path/$asm.readsToCorrect.stats") or caExit("can't open '$path/$asm.readsToCorrect.stats' for reading: $!", undef);
    while (<F>) {
        $report .= "--   $_";
    }
    close(F);

    addToReport("corLayout", $report);

    generateReport($asm);
    resetIteration("cor-filterCorrectionLayouts");

  allDone:
}



sub generateCorrectedReadsConfigure ($) {
    my $asm     = shift @_;
    my $bin     = getBinDirectory();
    my $cmd;

    my $base    = "correction";
    my $path    = "correction/2-correction";

    #  We cannot do this; we still need to figure out if corMemory is set
    #  appropriately.
    #goto allDone   if (fileExists("$path/correctReads.sh"));

    make_path("$path/results")  if (! -d "$path/results");

    #  Set corMemory to the maximum of:
    #    minmem:  genome size based minimum size
    #    estmem:  falconsense estimated memory requirement + 2 GB
    #    cormem:  user supplied corMemory
    #
    #  The genome size based minimum is to give more memory to larger genomes
    #  that (should) have more reads so falconsense can fit more evidence
    #  sequence in memory.

    my $minmem   = 0;
    my $estmem   = 0;
    my $maxmem   = getGlobal("corMemory");

    #  Figure out the genome size based minimum memory.

    if      (getGlobal("genomeSize") < adjustGenomeSize("40m")) {
        $minmem =  6;
    } elsif (getGlobal("genomeSize") < adjustGenomeSize("1g")) {
        $minmem = 12;
    } else {
        $minmem = 16;
    }

    #  Find the largest estimated memory needed to compute a corrected read.
    #
    #  This does not include overhead for opening stores or for loading reads
    #  into memory.

    fetchFile("$path/$asm.readsToCorrect.stats");

    if (! -e "$path/$asm.readsToCorrect.stats") {
        caExit("failed to find `$path/$asm.readsToCorrect.stats` to read maximum estimated memory needed for correction", undef);
    }

    open(F, "< $path/$asm.readsToCorrect.stats") or caExit("can't open '$path/$asm.readsToCorrect.stats' for reading: $!", undef);
    while (<F>) {
        if (m/Maximum\s+Memory\s+(\d+)/) {
            $estmem = ($estmem < $1) ? $1 : $estmem;
        }
    }
    close(F);

    $estmem = int($estmem / 1024.0 / 1024.0 / 1024.0 * 1000) / 1000;   #  truncate trailing digits past X.XXX.

    #  Pick the maximum.  The estimated memory is padded by 2gb to give some
    #  extra room in case estimates are incorrect.

    my $reqmem;

    $reqmem =                         $minmem;
    $reqmem = ($estmem+2 > $reqmem) ? $estmem+2 : $reqmem;
    $reqmem = ($maxmem   > $reqmem) ? $maxmem   : $reqmem;

    #  If there was a corMemory request (maxmem), and it is smaller than what
    #  we want, fail.  On the otherhand, if the corMemory request is larger
    #  than what we estimate, use the original request.

    if (($maxmem > 0) && ($maxmem < $reqmem)) {
        caExit("not enough memory for correction; increase corMemory to at least $reqmem", undef); 
    }

    #  Reset corMemory to whatever our request is.  Then call getAllowedResources() to
    #  set up other stuff, like corConcurrency, based on memory and threads.

    {
        print STDERR "--\n";  
        print STDERR "-- Correction jobs estimated to need at most $estmem GB for computation.\n"; 
        print STDERR "-- Correction jobs will request $reqmem GB each.\n";

        setGlobal("corMemory", $reqmem);

        my ($err, $all) = getAllowedResources("", "cor", "", "", 0);

        print STDERR "--\n";
        print STDERR $all;
        print STDERR "--\n";
    }

    #  If both the batch output and the correction script exist, we're done here.

    if ((-e "$path/correctReadsPartition.batches") &&
        (-e "$path/correctReads.sh")) {
        goto allDone;
    }


    #  Generate a script to compute the partitioning of correction jobs, given
    #  the memory allowed per process and memory needed for a correction.

    if (! -e "$path/correctReadsPartition.batches") {
        my $par = getGlobal("corPartitions");
        my $rds = getGlobal("corPartitionMin");

        print STDERR "--\n";  
        print STDERR "-- Configuring correction jobs:\n";  
        print STDERR "--   Reads estimated to need at most $estmem GB for computation.\n"; 
        print STDERR "--   Jobs will request $reqmem GB each.\n";
 
        open(F, "> $path/correctReadsPartition.sh") or caExit("can't open '$path/correctReadsPartition.sh' for writing: $!", undef);

        print F "#!" . getGlobal("shell") . "\n";
        print F "\n";
        print F getBinDirectoryShellCode();
        print F "\n";
        print F "\$bin/falconsense \\\n";
        print F "  -partition $reqmem $estmem $par $rds \\\n";
        print F "  -S ../../$asm.seqStore \\\n";
        print F "  -C ../$asm.corStore \\\n";
        print F "  -R ./$asm.readsToCorrect \\\n"   if ( fileExists("$path/$asm.readsToCorrect"));
        print F "  -t  " .        getGlobal("corThreads")       . " \\\n";
        print F "  -cc " .        getGlobal("corMinCoverage")   . " \\\n";
        print F "  -cl " .        getGlobal("minReadLength")    . " \\\n";
        print F "  -oi " . (1.0 - getGlobal("corErrorRate"))    . " \\\n";
        print F "  -ol " .        getGlobal("minOverlapLength") . " \\\n";
        print F "  -p ./correctReadsPartition.WORKING \\\n";
        print F "&& \\\n";
        print F "mv ./correctReadsPartition.WORKING.batches ./correctReadsPartition.batches \\\n";
        print F "&& \\\n";
        print F "exit 0\n";
        print F "\n";
        print F "exit 1\n";

        close(F);

        makeExecutable("$path/correctReadsPartition.sh");
        stashFile("$path/correctReadsPartition.sh");

        if (runCommand($path, "./correctReadsPartition.sh > ./correctReadsPartition.err 2>&1")) {
            caExit("failed to partition reads for correction", "$path/correctReadsPartition.err");
        }

        unlink("$path/correctReadsPartition.err");

        stashFile("$path/correctReadsPartition.batches");
    }

    #  Generate a script for computing corrected reads, using the batches file
    #  as a template.

    if (! -e "$path/correctReads.sh") {
        open(F, "> $path/correctReads.sh") or caExit("can't open '$path/correctReads.sh' for writing: $!", undef);

        print F "#!" . getGlobal("shell") . "\n";
        print F "\n";
        print F getBinDirectoryShellCode();
        print F "\n";
        print F setWorkDirectoryShellCode($path);
        print F "\n";
        print F getJobIDShellCode();
        print F "\n";
        print F "bgnid=0\n";
        print F "endid=0\n";
        print F "\n";

        my $nJobs = 0;

        open(B, "< $path/correctReadsPartition.batches") or caExit("can't open '$path/correctReadsPartition.batches' for reading: $!", undef);
        $_ = <B>;    #  Skip header line 1
        $_ = <B>;    #  Skip header line 2
        while (<B>) {
            s/^\s+//;
            s/\s+$//;

            my ($jobID, $bgnID, $endID, $nReads, $reqmem) = split '\s+', $_;

            print  F "if [ \$jobid -eq $jobID ] ; then\n";
            printf F "  jobid=%04d\n", $jobID;   #  Parsed in Check() below.
            print  F "  bgnid=$bgnID\n";
            print  F "  endid=$endID\n";
            print  F "fi\n";

            $nJobs = $jobID;
        }
        close(B);

        print F "\n";
        print F "if [ \$bgnid -eq 0 ]; then\n";
        print F "  echo Error: Invalid job \$jobid requested, must be between 1 and $nJobs.\n";
        print F "  exit 1\n";
        print F "fi\n";
        print F "\n";

        print F "\n";
        print F "if [ -e ./results/\$jobid.cns ] ; then\n";
        print F "  echo Job finished successfully.\n";
        print F "  exit 0\n";
        print F "fi\n";
        print F "\n";
        print F "if [ ! -d ./results ] ; then\n";
        print F "  mkdir -p ./results\n";
        print F "fi\n";
        print F "\n";

        print F fetchSeqStoreShellCode($asm, $path, "");
        print F "\n";
        print F fetchTigStoreShellCode("correction/2-correction", $asm, "corStore", "001", "");
        print F "\n";
        print F fetchFileShellCode($path, "$asm.readsToCorrect", "");
        print F "\n";

        print F "seqStore=../../$asm.seqStore\n";
        print F "\n";

        my $stageDir = getGlobal("stageDirectory");
        my $storeDir = "$stageDir/$asm-\$jobid.seqStore";    #  jobid is expanded during the job, in the shell

        if (defined($stageDir)) {
            print F "#  Try to make the seqStore directory in the stage location.\n";
            print F "#  If that fails, fallback to using the original store.\n";
            print F "\n";
            print F "mkdir -p $storeDir\n";
            print F "\n";
            print F "if [ -d $storeDir ] ; then\n";
            print F "  seqStore=$storeDir\n";
            print F "  cp -p ../../$asm.seqStore/info      $storeDir/info\n";
            print F "  cp -p ../../$asm.seqStore/libraries $storeDir/libraries\n";
            print F "  cp -p ../../$asm.seqStore/reads*    $storeDir/\n";
            print F "  cp -p ../../$asm.seqStore/blobs.*   $storeDir/\n";
            print F "fi\n";
            print F "\n";
        }

        print F "\$bin/falconsense \\\n";
        print F "  -S \$seqStore \\\n";
        print F "  -C ../$asm.corStore \\\n";
        print F "  -R ./$asm.readsToCorrect \\\n"   if ( fileExists("$path/$asm.readsToCorrect"));
        print F "  -r \$bgnid-\$endid \\\n";
        print F "  -t  " .        getGlobal("corThreads")       . " \\\n";
        print F "  -cc " .        getGlobal("corMinCoverage")   . " \\\n";
        print F "  -cl " .        getGlobal("minReadLength")    . " \\\n";
        print F "  -oi " . (1.0 - getGlobal("corErrorRate"))    . " \\\n";
        print F "  -ol " .        getGlobal("minOverlapLength") . " \\\n";
        print F "  -p ./results/\$jobid.WORKING \\\n";
        print F "  -cns \\\n";
        print F "  > ./results/\$jobid.err 2>&1 \\\n";
        print F "&& \\\n";
        print F "mv ./results/\$jobid.WORKING.cns ./results/\$jobid.cns \\\n";
        print F "\n";

        print F stashFileShellCode("$path", "results/\$jobid.cns", "");

        if (defined($stageDir)) {
            print F "rm -rf $storeDir\n";
            print F "rmdir  $stageDir\n";
            print F "\n";
        }

        close(F);

        makeExecutable("$path/correctReads.sh");
        stashFile("$path/correctReads.sh");
    }

  finishStage:
    generateReport($asm);
    resetIteration("cor-generateCorrectedReadsConfigure");

  allDone:
    stopAfter("correctionConfigure");
}



sub generateCorrectedReadsCheck ($) {
    my $asm     = shift @_;
    my $attempt = getGlobal("canuIteration");
    my $bin     = getBinDirectory();

    my $path    = "correction/2-correction";

    #  Compute the size of seqStore for staging

    setGlobal("corStageSpace", getSizeOfSequenceStore($asm));

    #  Figure out if all the tasks finished correctly.

    my @successJobs;
    my @failedJobs;
    my $failureMessage = "";

    fetchFile("$path/correctReads.sh");

    open(F, "< $path/correctReads.sh") or caExit("can't open '$path/correctReads.sh' for reading: $!", undef);

    while (<F>) {
        if (m/^\s+jobid=(\d+)$/) {
            if (fileExists("$path/results/$1.cns")) {
                push @successJobs, "2-correction/results/$1.cns\n";

            } else {
                $failureMessage .= "--   job 2-correction/results/$1.cns FAILED.\n";
                push @failedJobs, $1;
            }
        }
    }

    close(F);

    #  Failed jobs, retry.

    if (scalar(@failedJobs) > 0) {

        #  If too many attempts, give up.

        if ($attempt >= getGlobal("canuIterationMax")) {
            print STDERR "--\n";
            print STDERR "-- Read correction jobs failed, tried $attempt times, giving up.\n";
            print STDERR $failureMessage;
            print STDERR "--\n";
            caExit(undef, undef);
        }

        if ($attempt > 0) {
            print STDERR "--\n";
            print STDERR "-- Read correction jobs failed, retry.\n";
            print STDERR $failureMessage;
            print STDERR "--\n";
        }

        #  Otherwise, run some jobs.

        generateReport($asm);

        submitOrRunParallelJob($asm, "cor", $path, "correctReads", @failedJobs);
        return;
    }

  finishStage:
    print STDERR "-- Found ", scalar(@successJobs), " read correction output files.\n";

    open(L, "> $path/corjob.files") or caExit("failed to open '$path/corjob.files'", undef);
    print L @successJobs;
    close(L);

    stashFile("$path/corjob.files");

    generateReport($asm);
    resetIteration("cor-generateCorrectedReadsCheck");

  allDone:
}




sub loadCorrectedReads ($) {
    my $asm     = shift @_;
    my $bin     = getBinDirectory();
    my $cmd;

    my $base    = "correction";
    my $path    = "correction/2-correction";

    goto allDone   if (getNumberOfBasesInStore($asm, "obt") > 0);

    print STDERR "--\n";
    print STDERR "-- Loading corrected reads into corStore and seqStore.\n";

    #  Grab the correction outputs.

    fetchFile("$path/corjob.files");

    open(F, "< $path/corjob.files") or caExit("failed to open '$path/corjob.files' for reading: $!", undef);
    while (<F>) {
        chomp;
        fetchFile("correction/$_");
    }
    close(F);

    #  Grab the stores we're going to load into.

    fetchFile("$base/$asm.corStore/seqDB.v001.dat");
    fetchFile("$base/$asm.corStore/seqDB.v001.tig");

    #  Load the results into the stores.

    $cmd  = "$bin/loadCorrectedReads \\\n";
    $cmd .= "  -S ../$asm.seqStore \\\n";
    $cmd .= "  -C ./$asm.corStore \\\n";
    $cmd .= "  -L ./2-correction/corjob.files \\\n";
    $cmd .= ">  ./$asm.loadCorrectedReads.log \\\n";
    $cmd .= "2> ./$asm.loadCorrectedReads.err";

    if (runCommand($base, $cmd)) {
        caExit("failed to load corrected reads into store", "$base/$asm.loadCorrectedReads.err");
    }

    unlink("$base/$asm.loadCorrectedReads.err");

    #  Summarize and save updated stores.

    generateReadLengthHistogram("obt", $asm);
    stashSeqStore($asm);

    stashFile("$base/$asm.corStore/seqDB.v002.dat");
    stashFile("$base/$asm.corStore/seqDB.v002.tig");

    #  Now that all outputs are (re)written, cleanup the job outputs.
    #  (unless there are no corrected reads, then leave things alone for debugging)

    my $Ncns     = 0;
    my $Nerr     = 0;
    my $Nlog     = 0;

    if (getNumberOfReadsInStore($asm, "obt") == 0) {
        print STDERR "--\n";
        print STDERR "-- No corrected reads generated; correctReads output saved.\n";
    }
    elsif (getGlobal("saveReadCorrections") != 1) {
        print STDERR "--\n";
        print STDERR "-- Purging correctReads output after loading into stores.\n";

        open(F, "< $path/corjob.files") or caExit("can't open '$path/corjob.files' for reading: $!", undef);
        while (<F>) {
            chomp;

            if (m/^(.*)\/results\/0*(\d+).cns$/) {
                my $ID6 = substr("00000" . $2, -6);
                my $ID4 = substr("000"   . $2, -4);
                my $ID0 = $2;

                if (-e "correction/$1/results/$ID4.cns")      { $Ncns++;  unlink "correction/$1/results/$ID4.cns";      }
                if (-e "correction/$1/results/$ID4.err")      { $Nlog++;  unlink "correction/$1/results/$ID4.err";      }

                if (-e "correction/$1/correctReads.$ID6.out") { $Nlog++;  unlink "correction/$1/correctReads.$ID6.out"; }
                if (-e "correction/$1/correctReads.$ID0.out") { $Nlog++;  unlink "correction/$1/correctReads.$ID0.out"; }

            } else {
                caExit("unknown correctReads job name '$_'\n", undef);
            }
        }
        close(F);

        print STDERR "-- Purged $Ncns .cns outputs.\n"                  if ($Ncns > 0);
        print STDERR "-- Purged $Nerr .err outputs.\n"                  if ($Nerr > 0);
        print STDERR "-- Purged $Nlog .out job log outputs.\n"          if ($Nlog > 0);
    }
    else {
        print STDERR "--\n";
        print STDERR "-- Purging correctReads output disabled by saveReadCorrections=true.\n"  if (getGlobal("saveReadCorrections") == 1);
    }

    #  And purge the usually massive overlap store.

    if      (getNumberOfReadsInStore($asm, "obt") > 0) {
        print STDERR "--\n";
        print STDERR "-- No corrected reads generated, overlaps used for correction saved.\n";
    }
    elsif (getGlobal("saveOverlaps") eq "0") {
        print STDERR "--\n";
        print STDERR "-- Purging overlaps used for correction.\n";
        remove_tree("correction/$asm.ovlStore")
    }
    else {
        print STDERR "--\n";
        print STDERR "-- Overlaps used for correction saved.\n";
    }

  finishStage:
    generateReport($asm);
    resetIteration("cor-loadCorrectedReads");

  allDone:
}




sub dumpCorrectedReads ($) {
    my $asm     = shift @_;
    my $bin     = getBinDirectory();
    my $cmd;

    goto allDone   if (fileExists("$asm.correctedReads.fasta.gz"));
    goto allDone   if (getGlobal("saveReads") == 0);

    $cmd  = "$bin/sqStoreDumpFASTQ \\\n";
    $cmd .= "  -corrected \\\n";
    $cmd .= "  -S ./$asm.seqStore \\\n";
    $cmd .= "  -o ./$asm.correctedReads.gz \\\n";
    $cmd .= "  -fasta \\\n";
    $cmd .= "  -nolibname \\\n";
    $cmd .= "> $asm.correctedReads.fasta.err 2>&1";

    if (runCommand(".", $cmd)) {
        caExit("failed to output corrected reads", "./$asm.correctedReads.fasta.err");
    }

    unlink "./$asm.correctedReads.fasta.err";

    stashFile("$asm.correctedReads.fasta.gz");

    print STDERR "--\n";
    print STDERR "-- Corrected reads saved in '$asm.correctedReads.fasta.gz'.\n";

  finishStage:
    generateReport($asm);
    resetIteration("cor-dumpCorrectedReads");

  allDone:
    stopAfter("correction");
}
