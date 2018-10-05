
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



#  Return the number of correction jobs.
#
sub computeNumberOfCorrectionJobs ($) {
    my $asm     = shift @_;
    my $nJobs   = 0;
    my $nPerJob = 0;

    my $nPart    = getGlobal("corPartitions");
    my $nReads   = getNumberOfReadsInStore($asm, "all");

    $nPerJob     = int($nReads / $nPart + 1);
    $nPerJob     = getGlobal("corPartitionMin")  if ($nPerJob < getGlobal("corPartitionMin"));

    for (my $j=1; $j<=$nReads; $j += $nPerJob) {  #  We could just divide, except for rounding issues....
        $nJobs++;
    }

    return($nJobs, $nPerJob);
}



sub estimateMemoryNeededForCorrectionJobs ($) {
    my $asm     = shift @_;
    my $bin     = getBinDirectory();
    my $cmd;

    my $path    = "correction/2-correction";

    my $readID   = 0;
    my $readLen  = 0;
    my $numOlaps = 0;
    my $alignLen = 0;
    my $memEst   = 0;

    return   if (defined(getGlobal("corMemory")));

    fetchFile("$path/$asm.readsToCorrect.stats");

    if (-e "$path/$asm.readsToCorrect.stats") {
        open(F, "< $path/$asm.readsToCorrect.stats") or caExit("can't open '$path/$asm.readsToCorrect.stats' for reading: $!", undef);
        while (<F>) {
            if (m/Maximum\s+Memory\s+(\d+)/) {
                $memEst = int(4 * $1 / 1073741824.0 + 0.5);
            }
        }
        close(F);
    }

    if ($memEst == 0) {
        $memEst = 12;
    }

    setGlobal("corMemory", $memEst);

    my $err;
    my $all;

    ($err, $all) = getAllowedResources("", "cor", $err, $all, 0);

    print STDERR "--\n";
    print STDERR $all;
    print STDERR "--\n";
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

    estimateMemoryNeededForCorrectionJobs($asm);

    my ($nJobs, $nPerJob)  = computeNumberOfCorrectionJobs($asm);  #  Does math based on number of reads and parameters.

    my $nReads             = getNumberOfReadsInStore($asm, "all");

    open(F, "> $path/correctReads.sh") or caExit("can't open '$path/correctReads.sh' for writing: $!", undef);

    print F "#!" . getGlobal("shell") . "\n";
    print F "\n";
    print F getBinDirectoryShellCode();
    print F "\n";
    print F setWorkDirectoryShellCode($path);
    print F "\n";
    print F getJobIDShellCode();
    print F "\n";
    print F "if [ \$jobid -gt $nJobs ]; then\n";
    print F "  echo Error: Only $nJobs partitions, you asked for \$jobid.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "\n";

    my  $bgnID   = 1;
    my  $endID   = $bgnID + $nPerJob - 1;
    my  $jobID   = 1;

    while ($bgnID < $nReads) {
        $endID  = $bgnID + $nPerJob - 1;
        $endID  = $nReads  if ($endID > $nReads);

        print F "if [ \$jobid -eq $jobID ] ; then\n";
        print F "  bgn=$bgnID\n";
        print F "  end=$endID\n";
        print F "fi\n";

        $bgnID = $endID + 1;
        $jobID++;
    }

    print F "\n";
    print F "jobid=`printf %04d \$jobid`\n";
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
    print F fetchFileShellCode("$base/$asm.corStore", "seqDB.v001.dat", "");
    print F fetchFileShellCode("$base/$asm.corStore", "seqDB.v001.tig", "");
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
    print F "  -r \$bgn-\$end \\\n";
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

    if (defined($stageDir)) {
        print F "rm -rf $stageDir/$asm.seqStore\n";   #  Prevent accidents of 'rm -rf /' if stageDir = "/".
        print F "rmdir  $stageDir\n";
        print F "\n";
    }

    print F stashFileShellCode("$path", "results/\$jobid.cns", "");

    print F "\n";
    print F "exit 0\n";

    close(F);

    makeExecutable("$path/correctReads.sh");
    stashFile("$path/correctReads.sh");


  finishStage:
    generateReport($asm);
    resetIteration("cor-generateCorrectedReadsConfigure");

  allDone:
}



sub generateCorrectedReadsCheck ($) {
    my $asm     = shift @_;
    my $attempt = getGlobal("canuIteration");
    my $bin     = getBinDirectory();

    my $path    = "correction/2-correction";

    #  Compute the size of seqStore for staging

    setGlobal("corStageSpace", getSizeOfSequenceStore($asm));

    #  Compute expected size of jobs, set if not set already.

    estimateMemoryNeededForCorrectionJobs($asm);

    #  Figure out if all the tasks finished correctly.

    fetchFile("$path/correctReads.sh");

    my ($jobs, undef) = computeNumberOfCorrectionJobs($asm);

    my $currentJobID = "0001";
    my @successJobs;
    my @failedJobs;
    my $failureMessage = "";

    for (my $job=1; $job <= $jobs; $job++) {
        if (fileExists("$path/results/$currentJobID.cns")) {
            push @successJobs, "2-correction/results/$currentJobID.cns\n";

        } else {
            $failureMessage .= "--   job 2-correction/results/$currentJobID.cns FAILED.\n";
            push @failedJobs, $job;
        }

        $currentJobID++;
    }

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

    if (runCommand("correction", $cmd)) {
        caExit("failed to load corrected reads into store", "$path/$asm.loadCorrectedReads.err");
    }

    unlink("$path/$asm.loadCorrectedReads.err");

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
    stopAfter("readCorrection");
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
    stopAfter("readCorrection");
}
