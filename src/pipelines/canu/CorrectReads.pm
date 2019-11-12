
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
 #    src/pipelines/ca3g/CorrectReads.pm
 #
 #  Modifications by:
 #
 #    Brian P. Walenz from 2015-APR-09 to 2015-SEP-03
 #      are Copyright 2015 Battelle National Biodefense Institute, and
 #      are subject to the BSD 3-Clause License
 #
 #    Brian P. Walenz beginning on 2015-OCT-19
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #    Sergey Koren beginning on 2015-NOV-17
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
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

    my $exp = getExpectedCoverage("cor", $asm);
    my $des = getGlobal("corOutCoverage");

    if (!defined($cov)) {
        $cov = $des;
    } elsif ($cov =~ m/(.*)x/) {
        $cov = int($des * $1);
    }

    return($cov);
}


#  Query seqStore to find the read types involved.  Return an error rate that is appropriate for
#  aligning reads of that type to each other.
sub getCorIdentity ($) {
    my $asm     = shift @_;
    my $bin     = getBinDirectory();
    my $erate   = getGlobal("corErrorRate");

    if (defined($erate)) {
        print STDERR "-- Using overlaps no worse than $erate fraction error for correcting reads (from corErrorRate parameter).\n";
        return(1 - $erate);
    }

    my $numPacBioRaw         = 0;
    my $numPacBioCorrected   = 0;
    my $numNanoporeRaw       = 0;
    my $numNanoporeCorrected = 0;

    open(L, "< ./$asm.seqStore/libraries.txt") or caExit("can't open './$asm.seqStore/libraries.txt' for reading: $!", undef);
    while (<L>) {
        $numPacBioRaw++           if (m/pacbio-raw/);
        $numPacBioCorrected++     if (m/pacbio-corrected/);
        $numNanoporeRaw++         if (m/nanopore-raw/);
        $numNanoporeCorrected++   if (m/nanopore-corrected/);
    }
    close(L);

    $erate = 0.10;                              #  Default; user is stupid and forced correction of corrected reads.
    $erate = 0.30   if ($numPacBioRaw   > 0);
    $erate = 0.50   if ($numNanoporeRaw > 0);

    print STDERR "-- Found $numPacBioRaw raw and $numPacBioCorrected corrected PacBio libraries.\n";
    print STDERR "-- Found $numNanoporeRaw raw and $numNanoporeCorrected corrected Nanopore libraries.\n";
    print STDERR "-- Using overlaps no worse than $erate fraction error for correcting reads.\n";

    return(1 - $erate);
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
        my $cov = getExpectedCoverage("cor", $asm);

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

    goto allDone   if (-d "$base/$asm.corStore");                         #  Jobs all finished
    goto allDone   if (fileExists("$base/$asm.corStore/seqDB.v001.dat"));
    goto allDone   if (fileExists("$base/$asm.corStore/seqDB.v001.tig"));

    #  The global filter can be estimated from data saved in ovlStore.  This code will compute it exactly.
    #
    #  IT HAS NOT BEEN UPDATED OR TESTED.

    fetchFile("$path/$asm.globalScores");

    my $computeGlobalScores = 0;

    if ($computeGlobalScores) {
        if (! fileExists("$path/$asm.globalScores")) {
            print STDERR "-- Computing global filter scores '$path/$asm.globalScores'.\n";

            fetchOvlStore($asm, $base);

            $cmd  = "$bin/filterCorrectionOverlaps \\\n";
            $cmd .= "  -estimate -nolog \\\n";
            $cmd .= "  -S ../../$asm.seqStore \\\n";
            $cmd .= "  -O    ../$asm.ovlStore \\\n";
            $cmd .= "  -scores ./$asm.globalScores.WORKING \\\n";
            $cmd .= "  -c " . getCorCov($asm, "Global") . " \\\n";
            $cmd .= "  -l " . getGlobal("corMinEvidenceLength") . " \\\n"  if (defined(getGlobal("corMinEvidenceLength")));
            $cmd .= "  -e " . getGlobal("corMaxEvidenceErate")  . " \\\n"  if (defined(getGlobal("corMaxEvidenceErate")));
            $cmd .= "> ./$asm.globalScores.err 2>&1";

            if (runCommand($path, $cmd)) {
                caExit("failed to globally filter overlaps for correction", "$path/$asm.globalScores.err");
            }

            rename "$path/$asm.globalScores.WORKING",       "$path/$asm.globalScores";
            rename "$path/$asm.globalScores.WORKING.stats", "$path/$asm.globalScores.stats";
            rename "$path/$asm.globalScores.WORKING.log",   "$path/$asm.globalScores.log";
            unlink "$path/$asm.globalScores.err";

            stashFile("$path/$asm.globalScores");

            my $report = getFromReport("corFilter");

            open(F, "< $path/$asm.globalScores.stats") or caExit("can't open '$path/$asm.globalScores.stats' for reading: $!", undef);
            while(<F>) {
                $report .= "--  $_";
            }
            close(F);

            addToReport("corFilter", $report);

        } else {
            print STDERR "-- Global filter scores found in '$path/$asm.globalScores'.\n";
        }
    } else {
        print STDERR "-- Global filter scores will be estimated.\n";
    }

    #  Make layouts for each corrected read.

    fetchOvlStore($asm, $base);

    print STDERR "-- Computing correction layouts.\n";

    $cmd  = "$bin/generateCorrectionLayouts \\\n";
    $cmd .= "  -S ../$asm.seqStore \\\n";
    $cmd .= "  -O  ./$asm.ovlStore \\\n";
    $cmd .= "  -C  ./$asm.corStore.WORKING \\\n";
    $cmd .= "  -scores 2-correction/$asm.globalScores \\\n"         if (-e "$path/$asm.globalScores");
    $cmd .= "  -eL " . getGlobal("corMinEvidenceLength") . " \\\n"  if (defined(getGlobal("corMinEvidenceLength")));
    $cmd .= "  -eE " . getGlobal("corMaxEvidenceErate")  . " \\\n"  if (defined(getGlobal("corMaxEvidenceErate")));
    $cmd .= "  -eC " . getCorCov($asm, "Local") . " \\\n";
    $cmd .= "> ./$asm.corStore.err 2>&1";

    if (runCommand($base, $cmd)) {
        caExit("failed to generate correction layouts", "$base/$asm.corStore.err");
    }

    rename "$base/$asm.corStore.WORKING", "$base/$asm.corStore";
    unlink "$base/$asm.corStore.err";

    stashFile("$base/$asm.corStore/seqDB.v001.dat");
    stashFile("$base/$asm.corStore/seqDB.v001.tig");

  finishStage:
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

    goto allDone   if (fileExists("$path/$asm.readsToCorrect"));                #  Jobs all finished

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

  finishStage:
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

    goto allDone   if (fileExists("$path/correctReads.sh"));

    make_path("$path/results")  if (! -d "$path/results");

    #  Figure out the maximum memory needed to correct a read, then configure correction
    #  jobs so they fit in memory.

    my $cnsmem   = 1 * 1024 * 1024 * 1024;
    my $mem      = getGlobal("corMemory");
    my $par      = getGlobal("corPartitions");
    my $rds      = getGlobal("corPartitionMin");
    my $remain   = 0;

    fetchFile("$path/$asm.readsToCorrect.stats");

    if (! -e "$path/$asm.readsToCorrect.stats") {
        caExit("failed to find `$path/$asm.readsToCorrect.stats` to read maximum estimated memory needed for correction", undef);
    }

    #  Read the (vast overestimate) of memory needed to compute consensus for each read,
    #  and return the maximum.  This does not include overhead for opening stores or
    #  for loading reads into memory.

    open(F, "< $path/$asm.readsToCorrect.stats") or caExit("can't open '$path/$asm.readsToCorrect.stats' for reading: $!", undef);
    while (<F>) {
        if (m/Maximum\s+Memory\s+(\d+)/) {
            $cnsmem = ($cnsmem < $1) ? $1 : $cnsmem;
        }
    }
    close(F);

    $cnsmem /= 1024.0;
    $cnsmem /= 1024.0;
    $cnsmem /= 1024.0;

    $cnsmem = int($cnsmem * 1000) / 1000;
    $remain = int(($mem - $cnsmem) * 1000) / 1000;

    #  Generate a script to compute the partitioning of correction jobs, given
    #  the memory allowed per process and memory needed for a correction.

    open(F, "> $path/correctReadsPartition.sh") or caExit("can't open '$path/correctReadsPartition.sh' for writing: $!", undef);

    print F "#!" . getGlobal("shell") . "\n";
    print F "set -eu\n";
    print F "\n";
    print F getBinDirectoryShellCode();
    print F "\n";
    print F "\$bin/falconsense \\\n";
    print F "  -partition $mem $cnsmem $par $rds \\\n";
    print F "  -S ../../$asm.seqStore \\\n";
    print F "  -C ../$asm.corStore \\\n";
    print F "  -R ./$asm.readsToCorrect \\\n"                if ( fileExists("$path/$asm.readsToCorrect"));
    print F "  -t  " . getGlobal("corThreads") . " \\\n";
    print F "  -cc " . getGlobal("corMinCoverage") . " \\\n";
    print F "  -cl " . getGlobal("minReadLength") . " \\\n";
    print F "  -oi " . getCorIdentity($asm) . " \\\n";
    print F "  -ol " . getGlobal("minOverlapLength") . " \\\n";
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

    print STDERR "--\n";
    print STDERR "-- Configuring correction jobs:\n";
    print STDERR "--   Jobs limited to $mem GB per job (via option corMemory).\n";
    print STDERR "--   Reads estimated to need at most $cnsmem GB for computation.\n";
    print STDERR "--   Leaving $remain GB memory for read data.\n";

    if ($remain < 1.0) {
        #print STDERR "--\n";
        #print STDERR "-- ERROR: Not enough memory for correction.  Increase corMemory.\n";
        caExit("not enough memory for correction; increase corMemory", undef);
    }

    if (runCommand($path, "./correctReadsPartition.sh > ./correctReadsPartition.err 2>&1")) {
        caExit("failed to partition reads for correction", "$path/correctReadsPartition.err");
    }

    stashFile("$path/correctReadsPartition.batches");
    unlink("$path/correctReadsPartition.err");

    #  Generate a script for computing corrected reads, using the batches file
    #  as a template.

    open(F, "> $path/correctReads.sh") or caExit("can't open '$path/correctReads.sh' for writing: $!", undef);

    print F "#!" . getGlobal("shell") . "\n";
    print F "set -eu\n";
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

        my ($jobID, $bgnID, $endID, $nReads, $mem) = split '\s+', $_;

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
    print F "if [ -e \"./results/\$jobid.cns\" ] ; then\n";
    print F "  echo Job finished successfully.\n";
    print F "  exit 0\n";
    print F "fi\n";
    print F "\n";
    print F "if [ ! -d \"./results\" ] ; then\n";
    print F "  mkdir -p \"./results\"\n";
    print F "fi\n";
    print F "\n";

    print F fetchSeqStoreShellCode($asm, $path, "");
    print F "\n";
    print F fetchTigStoreShellCode("correction/2-correction", $asm, "corStore", "001", "");
    print F "\n";
    print F fetchFileShellCode($path, "$asm.readsToCorrect", "");
    print F "\n";

    print F "seqStore=\"../../$asm.seqStore\"\n";
    print F "\n";

    my $stageDir = getGlobal("stageDirectory");

    if (defined($stageDir)) {
        print F "if [ ! -d $stageDir ] ; then\n";
        print F "  mkdir -p $stageDir\n";
        print F "fi\n";
        print F "\n";
        print F "mkdir -p $stageDir/$asm.seqStore\n";
        print F "\n";
        print F "echo Start copy at `date`\n";
        print F "cp -p \$seqStore/info      $stageDir/$asm.seqStore/info\n";
        print F "cp -p \$seqStore/libraries $stageDir/$asm.seqStore/libraries\n";
        print F "cp -p \$seqStore/reads     $stageDir/$asm.seqStore/reads\n";
        print F "cp -p \$seqStore/blobs.*   $stageDir/$asm.seqStore/\n";
        print F "echo Finished   at `date`\n";
        print F "\n";
        print F "seqStore=\"$stageDir/$asm.seqStore\"\n";
        print F "\n";
    }

    print F "\n";
    print F "\$bin/falconsense \\\n";
    print F "  -S \$seqStore \\\n";
    print F "  -C ../$asm.corStore \\\n";
    print F "  -R ./$asm.readsToCorrect \\\n"                if ( fileExists("$path/$asm.readsToCorrect"));
    print F "  -r \$bgnid-\$endid \\\n";
    print F "  -t  " . getGlobal("corThreads") . " \\\n";
    print F "  -cc " . getGlobal("corMinCoverage") . " \\\n";
    print F "  -cl " . getGlobal("minReadLength") . " \\\n";
    print F "  -oi " . getCorIdentity($asm) . " \\\n";
    print F "  -ol " . getGlobal("minOverlapLength") . " \\\n";
    print F "  -p ./results/\$jobid.WORKING \\\n";
    print F "  -cns \\\n";
    print F "  > ./results/\$jobid.err 2>&1 \\\n";
    print F "&& \\\n";
    print F "mv ./results/\$jobid.WORKING.cns ./results/\$jobid.cns \\\n";
    print F "\n";

    print F stashFileShellCode("$path", "results/\$jobid.cns", "");

    if (defined($stageDir)) {
        print F "rm -rf $stageDir/$asm.seqStore\n";   #  Prevent accidents of 'rm -rf /' if stageDir = "/".
        print F "rmdir  $stageDir\n";
        print F "\n";
    }

    print F "\n";
    print F "exit 0\n";

    close(F);

    makeExecutable("$path/correctReads.sh");
    stashFile("$path/correctReads.sh");


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

    #  Save updated stores.

    stashSeqStore($asm);

    stashFile("$base/$asm.corStore/seqDB.v002.dat");
    stashFile("$base/$asm.corStore/seqDB.v002.tig");

    #  Report reads.

    addToReport("obtSeqStore", generateReadLengthHistogram("obt", $asm));

    #  Now that all outputs are (re)written, cleanup the job outputs.

    my $Ncns     = 0;
    my $Nerr     = 0;
    my $Nlog     = 0;

    if (getGlobal("saveReadCorrections") != 1) {
        print STDERR "--\n";
        print STDERR "-- Purging correctReads output after loading into stores.\n";

        open(F, "< $path/corjob.files") or caExit("can't open '$path/corjob.files' for reading: $!", undef);
        while (<F>) {
            chomp;

            if (m/^(.*)\/results\/0*(\d+).cns$/) {
                my $ID6 = substr("00000" . $2, -6);
                my $ID4 = substr("000"   . $2, -4);
                my $ID0 = $2;

                if (-e "correction/$1/results/$ID4.cns") {
                    $Ncns++;
                    unlink "correction/$1/results/$ID4.cns";
                }
                if (-e "correction/$1/results/$ID4.err") {
                    $Nlog++;
                    unlink "correction/$1/results/$ID4.err";
                }

                if (-e "correction/$1/correctReads.$ID6.out") {
                    $Nlog++;
                    unlink "correction/$1/correctReads.$ID6.out";
                }
                if (-e "correction/$1/correctReads.$ID0.out") {
                    $Nlog++;
                    unlink "correction/$1/correctReads.$ID0.out";
                }

            } else {
                caExit("unknown correctReads job name '$_'\n", undef);
            }
        }
        close(F);

        print STDERR "-- Purged $Ncns .cns outputs.\n"                  if ($Ncns > 0);
        print STDERR "-- Purged $Nerr .err outputs.\n"                  if ($Nerr > 0);
        print STDERR "-- Purged $Nlog .out job log outputs.\n"          if ($Nlog > 0);
    } else {
        print STDERR "--\n";
        print STDERR "-- Purging correctReads output disabled by saveReadCorrections=true.\n"  if (getGlobal("saveReadCorrections") == 1);
    }

    #  And purge the usually massive overlap store.

    if (getGlobal("saveOverlaps") eq "0") {
        print STDERR "--\n";
        print STDERR "-- Purging overlaps used for correction.\n";

        remove_tree("correction/$asm.ovlStore")
    } else {
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
    goto allDone   if (fileExists("$asm.correctedReads.fastq.gz"));
    goto allDone   if (getGlobal("saveReads") == 0);

    #  We need to skip this entire function if corrected reads were not computed.
    #  Otherwise, we incorrectly declare that no corrected reads were generated
    #  and halt the assembly.
    #
    #  In cloud mode, we dump reads right after loading them, and to load them,
    #  the 2-correction directory must exist, so we use that to decide if correction
    #  was attempted.

    return         if (! -d "correction/2-correction");

    #  If no corrected reads exist, don't bother trying to dump them.

    if (getNumberOfReadsInStore($asm, "obt") > 0) {
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
    }

    #  If the corrected reads file exists, report so.
    #  Otherwise, report no corrected reads, and generate fake outputs so we terminate.

    my $out;

    $out = "$asm.correctedReads.fasta.gz"   if (fileExists("$asm.correctedReads.fasta.gz"));
    $out = "$asm.correctedReads.fastq.gz"   if (fileExists("$asm.correctedReads.fastq.gz"));

    if (defined($out)) {
        print STDERR "--\n";
        print STDERR "-- Corrected reads saved in '$out'.\n";
    } else {
        print STDERR "--\n";
        print STDERR "-- Yikes!  No corrected reads generated!\n";
        print STDERR "-- Can't proceed!\n";
        print STDERR "--\n";
        print STDERR "-- Generating empty outputs.\n";

        runCommandSilently(".", "gzip -1vc < /dev/null > $asm.correctedReads.gz 2> /dev/null", 0)   if (! -e "$asm.correctedReads.gz");
        runCommandSilently(".", "gzip -1vc < /dev/null > $asm.trimmedReads.gz   2> /dev/null", 0)   if (! -e "$asm.trimmedReads.gz");

        stashFile("$asm.correctedReads.fasta.gz");
        stashFile("$asm.trimmedReads.fasta.gz");

        generateOutputs($asm);
    }

  finishStage:
    generateReport($asm);
    resetIteration("cor-dumpCorrectedReads");

  allDone:
    stopAfter("correction");
}
