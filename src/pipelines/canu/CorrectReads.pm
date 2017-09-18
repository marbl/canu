
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
@EXPORT = qw(setupCorrectionParameters buildCorrectionLayoutsConfigure buildCorrectionLayoutsCheck filterCorrectionLayouts generateCorrectedReadsConfigure generateCorrectedReadsCheck dumpCorrectedReads);

use strict;

use File::Path 2.08 qw(make_path remove_tree);

use canu::Defaults;
use canu::Configure;
use canu::Execution;
use canu::Gatekeeper;
use canu::Report;
use canu::HTML;
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

    my $exp = getExpectedCoverage("correction", $asm);
    my $des = getGlobal("corOutCoverage");

    if (!defined($cov)) {
        $cov = $des;
    } elsif ($cov =~ m/(.*)x/) {
        $cov = int($des * $1);
    }

    return($cov);
}


#  Query gkpStore to find the read types involved.  Return an error rate that is appropriate for
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

    open(L, "< correction/$asm.gkpStore/libraries.txt") or caExit("can't open 'correction/$asm.gkpStore/libraries.txt' for reading: $!", undef);
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
    my $nReads   = getNumberOfReadsInStore("correction", $asm);

    caExit("didn't find any reads in store 'correction/$asm.gkpStore'?", undef)  if ($nReads == 0);

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
                $memEst = int($1 / 1073741824.0 + 0.5) * 2;
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

    ($err, $all) = getAllowedResources("", "cor", $err, $all);

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
        my $cov = getExpectedCoverage("correction", $asm);

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

    goto allDone   if (skipStage($asm, "cor-buildCorrectionLayouts") == 1);
    goto allDone   if (sequenceFileExists("$asm.correctedReads"));              #  Output exists
    goto allDone   if (fileExists("$base/$asm.corStore"));                      #  Jobs all finished
    goto allDone   if (fileExists("$path/generateCorrectionLayouts.sh"));       #  Jobs created

    #  The global filter can be estimated from data saved in ovlStore.  This code will compute it exactly.
    #
    #  IT HAS NOT BEEN UPDATED OR TESTED.

    fetchFile("$path/$asm.globalScores");

    my $computeGlobalScores = 0;

    if ($computeGlobalScores) {
        if (! fileExists("$path/$asm.globalScores")) {
            print STDERR "-- Computing global filter scores '$path/$asm.globalScores'.\n";

            fetchStore("./correction/$asm.ovlStore");

            $cmd  = "$bin/filterCorrectionOverlaps \\\n";
            $cmd .= "  -estimate -nolog \\\n";
            $cmd .= "  -G ../$asm.gkpStore \\\n";
            $cmd .= "  -O ../$asm.ovlStore \\\n";
            $cmd .= "  -S ./$asm.globalScores.WORKING \\\n";
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

    fetchStore("./correction/$asm.gkpStore");
    fetchStore("./correction/$asm.ovlStore");

    print STDERR "-- Computing correction layouts.\n";

    $cmd  = "$bin/generateCorrectionLayouts \\\n";
    $cmd .= "  -G ./$asm.gkpStore \\\n";
    $cmd .= "  -O ./$asm.ovlStore \\\n";
    $cmd .= "  -C ./$asm.corStore.WORKING \\\n";
    $cmd .= "  -S 2-correction/$asm.globalScores \\\n"              if (-e "$path/$asm.globalScores");
    $cmd .= "  -eL " . getGlobal("corMinEvidenceLength") . " \\\n"  if (defined(getGlobal("corMinEvidenceLength")));
    $cmd .= "  -eE " . getGlobal("corMaxEvidenceErate")  . " \\\n"  if (defined(getGlobal("corMaxEvidenceErate")));
    $cmd .= "  -ec " . getGlobal("corMinCoverage") . " \\\n";
    $cmd .= "  -eC " . getCorCov($asm, "Local") . " \\\n";
    $cmd .= "> ./$asm.corStore.err 2>&1\n";

    if (runCommand($base, $cmd)) {
        caExit("failed to generate correction layouts", "$base/$asm.corStore.err");
    }

    rename "$base/$asm.corStore.WORKING", "$base/$asm.corStore";
    unlink "$base/$asm.corStore.err";

    stashStore("./correction/$asm.corStore");

  finishStage:
    emitStage($asm, "cor-buildCorrectionLayoutsConfigure");
    buildHTML($asm, "cor");

  allDone:
}



sub buildCorrectionLayoutsCheck ($) {
    my $asm     = shift @_;
    my $bin     = getBinDirectory();
    my $cmd;

    my $base    = "correction";
    my $path    = "correction/2-correction";

    goto allDone   if (skipStage($asm, "cor-buildCorrectionLayouts") == 1);
    goto allDone   if (sequenceFileExists("$asm.correctedReads"));              #  Output exists
    goto allDone   if (fileExists("$base/$asm.corStore"));                      #  Jobs all finished
    goto allDone   if (fileExists("$path/generateCorrectionLayouts.sh"));       #  Jobs created

    #  Eventually, we'll run generateCorrectionLayouts on the grid.  Then we'll need to load the new
    #  tigs into the corStore here.

  finishStage:
    emitStage($asm, "cor-buildCorrectionLayoutsCheck");
    buildHTML($asm, "cor");

  allDone:
}



sub filterCorrectionLayouts ($) {
    my $asm     = shift @_;
    my $bin     = getBinDirectory();
    my $cmd;

    my $base    = "correction";
    my $path    = "correction/2-correction";

    goto allDone   if (skipStage($asm, "cor-buildCorrectionLayouts") == 1);
    goto allDone   if (sequenceFileExists("$asm.correctedReads"));              #  Output exists
    goto allDone   if (fileExists("$path/$asm.readsToCorrect"));                #  Jobs all finished

    #  Analyze the corStore to decide what reads we want to correct.

    fetchStore("./correction/$asm.gkpStore");
    fetchStore("./correction/$asm.ovlStore");

    my $genomeSize  = getGlobal("genomeSize");
    my $outCoverage = getGlobal("corOutCoverage");

    print STDERR "-- Computing correction layouts.\n";

    $cmd  = "$bin/filterCorrectionLayouts \\\n";
    $cmd .= "  -G ../$asm.gkpStore \\\n";
    $cmd .= "  -C ../$asm.corStore \\\n";
    $cmd .= "  -R ./$asm.readsToCorrect.WORKING \\\n";
    $cmd .= "  -g $genomeSize \\\n";
    $cmd .= "  -c $outCoverage \\\n";
    $cmd .= "> ./$asm.readsToCorrect.err 2>&1\n";

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
    emitStage($asm, "cor-filterCorrectionLayouts");
    buildHTML($asm, "cor");

  allDone:
}



sub generateCorrectedReadsConfigure ($) {
    my $asm     = shift @_;
    my $bin     = getBinDirectory();
    my $cmd;

    my $base    = "correction";
    my $path    = "correction/2-correction";

    make_path("$path/results")  if (! -d "$path/results");

    estimateMemoryNeededForCorrectionJobs($asm);

    my ($nJobs, $nPerJob)  = computeNumberOfCorrectionJobs($asm);  #  Does math based on number of reads and parameters.

    my $nReads             = getNumberOfReadsInStore("correction", $asm);

    open(F, "> $path/correctReads.sh") or caExit("can't open '$path/correctReads.sh' for writing: $!", undef);

    print F "#!" . getGlobal("shell") . "\n";
    print F "\n";
    print F getBinDirectoryShellCode();
    print F "\n";
    print F setWorkDirectoryShellCode($path);
    print F fetchStoreShellCode("$base/$asm.gkpStore", "$base/2-correction", "");
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
    print F "if [ -e \"./results/\$jobid.fasta\" ] ; then\n";
    print F "  echo Job finished successfully.\n";
    print F "  exit 0\n";
    print F "fi\n";
    print F "\n";
    print F "if [ ! -d \"./results\" ] ; then\n";
    print F "  mkdir -p \"./results\"\n";
    print F "fi\n";
    print F "\n";

    print F fetchStoreShellCode("correction/$asm.gkpStore", "correction/3-correction", "");
    print F "\n";
    print F fetchStoreShellCode("correction/$asm.ovlStore", "correction/3-correction", "");
    print F "\n";
    print F fetchFileShellCode("correction/2-correction", "$asm.readsToCorrect", "");
    print F "\n";
    print F fetchFileShellCode("correction/2-correction", "$asm.globalScores", "");
    print F "\n";

    print F "gkpStore=\"../$asm.gkpStore\"\n";
    print F "\n";

    my $stageDir = getGlobal("stageDirectory");

    if (defined($stageDir)) {
        print F "if [ ! -d $stageDir ] ; then\n";
        print F "  mkdir -p $stageDir\n";
        print F "fi\n";
        print F "\n";
        print F "mkdir -p $stageDir/$asm.gkpStore\n";
        print F "\n";
        print F "echo Start copy at `date`\n";
        print F "cp -p \$gkpStore/info      $stageDir/$asm.gkpStore/info\n";
        print F "cp -p \$gkpStore/libraries $stageDir/$asm.gkpStore/libraries\n";
        print F "cp -p \$gkpStore/reads     $stageDir/$asm.gkpStore/reads\n";
        print F "cp -p \$gkpStore/blobs     $stageDir/$asm.gkpStore/blobs\n";
        print F "echo Finished   at `date`\n";
        print F "\n";
        print F "gkpStore=\"$stageDir/$asm.gkpStore\"\n";
        print F "\n";
    }

    print F "\n";
    print F "\$bin/falconsense \\\n";
    print F "  -G \$gkpStore \\\n";
    print F "  -C ../$asm.corStore \\\n";
    print F "  -b \$bgn -e \$end -r ./$asm.readsToCorrect \\\n"     if (  -e "$path/$asm.readsToCorrect");
    print F "  -b \$bgn -e \$end \\\n"                              if (! -e "$path/$asm.readsToCorrect");
    print F "  -t  " . getGlobal("corThreads") . " \\\n";
    print F "  -ci " . getCorIdentity($asm) . "\\\n";
    print F "  -cl " . getGlobal("minReadLength") . "\\\n";
    print F "  -cc " . getGlobal("corMinCoverage") . " \\\n";
    print F "  > ./results/\$jobid.fasta.WORKING \\\n";
    print F " 2> ./results/\$jobid.err \\\n";
    print F "&& \\\n";
    print F "mv ./results/\$jobid.fasta.WORKING ./results/\$jobid.fasta \\\n";
    print F "\n";

    if (defined($stageDir)) {
        print F "rm -rf $stageDir/$asm.gkpStore\n";   #  Prevent accidents of 'rm -rf /' if stageDir = "/".
        print F "rmdir  $stageDir\n";
        print F "\n";
    }

    print F stashFileShellCode("$path", "results/\$jobid.fasta", "");

    print F "\n";
    print F "exit 0\n";

    close(F);

    makeExecutable("$path/correctReads.sh");
    stashFile("$path/correctReads.sh");


  finishStage:
    emitStage($asm, "cor-generateCorrectedReadsConfigure");
    buildHTML($asm, "cor");

  allDone:
}



sub generateCorrectedReadsCheck ($) {
    my $asm     = shift @_;
    my $attempt = getGlobal("canuIteration");
    my $bin     = getBinDirectory();

    my $path    = "correction/2-correction";

    goto allDone   if (skipStage($asm, "cor-generateCorrectedReads", $attempt) == 1);
    goto allDone   if (sequenceFileExists("$asm.correctedReads"));

    #  Compute the size of gkpStore for staging

    {
        my $size = 0;

        $size += -s "correction/$asm.gkpStore/info";
        $size += -s "correction/$asm.gkpStore/libraries";
        $size += -s "correction/$asm.gkpStore/reads";
        $size += -s "correction/$asm.gkpStore/blobs";

        $size = int($size / 1024 / 1024 / 1024 + 1.5);

        setGlobal("corStageSpace", $size);
    }

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
        if (fileExists("$path/results/$currentJobID.fasta")) {
            push @successJobs, "$path/results/$currentJobID.fasta\n";

        } else {
            $failureMessage .= "--   job $path/results/$currentJobID.fasta FAILED.\n";
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

        emitStage($asm, "cor-generateCorrectedReads", $attempt);
        buildHTML($asm, "cor");

        submitOrRunParallelJob($asm, "cor", $path, "correctReads", @failedJobs);
        return;
    }

  finishStage:
    print STDERR "-- Found ", scalar(@successJobs), " read correction output files.\n";

    open(L, "> $path/corjob.files") or caExit("failed to open '$path/corjob.files'", undef);
    print L @successJobs;
    close(L);

    stashFile("$path/corjob.files");

    emitStage($asm, "cor-generateCorrectedReadsCheck");
    buildHTML($asm, "cor");

  allDone:
}



sub dumpCorrectedReads ($) {
    my $asm     = shift @_;
    my $bin     = getBinDirectory();

    my $path    = "correction/2-correction";

    goto allDone   if (skipStage($asm, "cor-dumpCorrectedReads") == 1);
    goto allDone   if (sequenceFileExists("$asm.correctedReads"));

    print STDERR "-- Concatenating correctReads output.\n";

    my $files = 0;
    my $reads = 0;

    stashFile("$path/corjob.files");

    open(F, "< $path/corjob.files")                      or caExit("can't open '$path/corjob.files' for reading: $!", undef);
    open(N, "< correction/$asm.gkpStore/readNames.txt")  or caExit("can't open 'correction/$asm.gkpStore/readNames.txt' for reading: $!", undef);
    open(O, "| gzip -1c > $asm.correctedReads.fasta.gz") or caExit("can't open '$asm.correctedReads.fasta.gz' for writing: $!", undef);
    open(L, "> $asm.correctedReads.length")              or caExit("can't open '$asm.correctedReads.length' for writing: $!", undef);

    while (<F>) {
        chomp;

        fetchFile($_);

        open(R, "< $_") or caExit("can't open correction output '$_' for reading: $!\n", undef);

        my $h;   #  Current header line
        my $s;   #  Current sequence
        my $n;   #  Next header line

        my $nameid;  #  Currently loaded id and name from gkpStore/readNames
        my $name;

        $n = <R>;  chomp $n;  #  Read the first line, the first header.

        while (!eof(R)) {
            $h = $n;              #  Read name.
            $s = undef;           #  No sequence yet.
            $n = <R>;  chomp $n;  #  Sequence, or the next header.

            #  Read sequence until the next header or we stop reading lines.  Perl seems to be
            #  setting EOF when the last line is read, which is early, IMHO.  We loop until
            #  we both are EOF and have an empty line.

            while (($n !~ m/^>/) && ((length($n) > 0) || !eof(R))) {
                $s .= $n;

                $n = <R>;  chomp $n;
            }

            #  Parse the header of the corrected read to find the IID and the split piece number

            my $rid = undef;  #  Read ID
            my $pid = undef;  #  Piece ID

            if ($h =~ m/read(\d+)_(\d+)/) {
                $rid = $1;
                $pid = $2;
            }

            #  Load the next line from the gatekeeper ID map file until we find the
            #  correct one.

            while (!eof(N) && ($rid != $nameid)) {
                my $rn = <N>;

                ($rn =~ m/^(\d+)\s+(.*)$/);

                $nameid = $1;
                $name   = $2;
            }

            #  If a match, replace the header with the actual read id.  If no match, use the bogus
            #  corrected read name as is.

            if ($rid eq $nameid) {
                $h = ">$name id=${rid}_${pid}";
            }

            #  And write the read to the output as FASTA.

            print O "$h", "\n";
            print O "$s", "\n";

            #  Or as FASTQ

            #my $q = $s;
            #$n =~ s/^>/\@/;
            #$q =~ tr/[A-Z][a-z]/*/;

            #print O $h, "\n";
            #print O $s, "\n";
            #print O "+\n";
            #print O $q, "\n";

            print L "$rid\t$pid\t", length($s), "\n";

            $reads++;
        }

        $files++;
        close(R);
    }

    close(O);
    close(F);

    stashFile("$asm.correctedReads.fasta.gz");
    stashFile("$asm.correctedReads.length");

    #  Analyze the results.

    print STDERR "-- Analyzing correctReads output.\n";

    my @origLengthHist;
    my @expcLengthHist;
    my @corrLengthHist;

    #my @origExpcDiffLengthHist;
    #my @origCorrDiffLengthHist;
    #my @expcCorrDiffLengthHist;

    my @corrPiecesHist;

    my $inputReads     = 0;
    my $inputBases     = 0;

    my $failedReads    = 0;

    my $correctedReads = 0;
    my $correctedBases = 0;

    my $maxReadLen     = 0;

    my $minDiff = 0;
    my $maxDiff = 0;

    fetchFile("$path/$asm.readsToCorrect");
    fetchFile("$asm.correctedReads.length");   #  Already here, we just wrote it.

    open(R, "< $path/$asm.readsToCorrect")    or caExit("can't open '$path/$asm.readsToCorrect' for reading: $!", undef);
    open(L, "< $asm.correctedReads.length")   or caExit("can't open '$asm.correctedReads.length' for reading: $!", undef);

    open(S, "> $path/$asm.original-expected-corrected-length.dat") or caExit("", undef);

    my $moreData = 1;

    my $r = <R>;  chomp $r;    #  Read first 'read to correct' (input read)
    my @r = split '\s+', $r;

    if ($r[0] == 0) {
        $r = <R>;  chomp $r;    #  First line was header, read again.
        @r = split '\s+', $r;
    }

    my $l = <L>;  chomp $l;    #  Read first 'corrected read' (output read)
    my @l = split '\s+', $l;

    if ($l[0] == 0) {
        $l = <L>;  chomp $l;    #  First line was header, read again.
        @l = split '\s+', $l;
    }

    while ($moreData) {

        #  The 'read to correct' should be a superset of the 'corrected reads'.
        #  There will be only one 'read to correct', but possibly multiple (or even zero)
        #  'corrected reads'.  Read until the IDs differ.  For this, we only want to
        #  count the number of pieces and the total length.

        my $numPieces = 0;
        my $numBases  = 0;

        while ($l[0] == $r[0]) {
            $numPieces++;
            $numBases += $l[2];

            last   if (eof(L));  #  Ugly, break out if we just processed the last 'corrected read'.

            $l = <L>;  chomp $l;    #  Read next 'corrected read'
            @l = split '\s+', $l;
        }

        #  Add to the length scatter plot (original length, expected length, actual length).

        print S "$r[0]\t$r[1]\t$r[2]\t$numBases\t", $r[1] - $r[2], "\t", $r[1] - $numBases, "\t", $r[2] - $numBases, "\n";

        $minDiff = $r[2] - $numBases  if ($r[2] - $numBases < $minDiff);
        $maxDiff = $r[2] - $numBases  if ($r[2] - $numBases > $minDiff);

        #  Add to the histograms.

        $origLengthHist[$r[1]]++;
        $expcLengthHist[$r[2]]++;

        $corrLengthHist[$numBases]++;

        #$origExpcDiffLengthHist[$r[1] - $r[2]]++;
        #$origCorrDiffLengthHist[$r[1] - $numBases]++;
        #$expcCorrDiffLengthHist[$r[2] - $numBases]++;

        $corrPiecesHist[$numPieces]++;

        $inputReads++;
        $inputBases += $r[1];

        $failedReads++  if ($numBases == 0);

        $correctedReads += $numPieces;
        $correctedBases += $numBases;

        $maxReadLen = $r[1]      if ($maxReadLen < $r[1]);
        $maxReadLen = $r[2]      if ($maxReadLen < $r[2]);
        $maxReadLen = $numBases  if ($maxReadLen < $numBases);

        #  Read next 'read to correct'.  If this reads the last line, eof(R) is true, so can't
        #  loop on that.  Instead, we set a flag after we process that last line.

        if (!eof(R)) {
            $r = <R>;  chomp $r;
            @r = split '\s+', $r;
        } else {
            $moreData = 0;
        }
    }

    close(S);
    close(R);
    close(L);

    stashFile("$path/$asm.original-expected-corrected-length.dat");

    #  Write a summary of the corrections.

    my $report = getFromReport("corrections");

    $report .= "-- Actual corrected reads:\n";
    $report .= "--   $correctedReads reads\n";
    $report .= "--   $correctedBases bp\n";

    for (my $ii=0; $ii<scalar(@corrPiecesHist); $ii++) {
        $corrPiecesHist[$ii] += 0;
        $report .= "--   $corrPiecesHist[$ii] reads with $ii corrected block" . (($ii == 1) ? "" : "s") . "\n";
    }

    addToReport("corrections", $report);


#OBSOLETE
    open(F, "> $path/$asm.correction.summary") or caExit("", undef);
    print F "CORRECTION INPUTS:\n";
    print F "-----------------\n";
    print F "$inputReads (input reads)\n";
    print F "$inputBases (input bases)\n";
    print F "\n";
    print F "CORRECTION OUTPUTS:\n";
    print F "-----------------\n";
    print F "$failedReads (reads that failed to generate any corrected bases)\n";
    print F "$correctedReads (corrected read pieces)\n";
    print F "$correctedBases (corrected bases)\n";
    print F "\n";
    print F "PIECES PER READ:\n";
    print F "---------------\n";

    for (my $ii=0; $ii<scalar(@corrPiecesHist); $ii++) {
        $corrPiecesHist[$ii] += 0;
        print F "$ii pieces: $corrPiecesHist[$ii]\n";
    }
    close(F);

    stashFile("$path/$asm.correction.summary");

    #  Scatterplot of lengths.

    my $gnuplot = getGlobal("gnuplot");
    my $format  = getGlobal("gnuplotImageFormat");

    open(F, "> $path/$asm.originalLength-vs-correctedLength.gp") or caExit("", undef);
    print F "\n";
    print F "set pointsize 0.25\n";
    print F "\n";
    print F "set title 'original read length vs expected corrected read length'\n";
    print F "set xlabel 'original read length'\n";
    print F "set ylabel 'expected corrected read length'\n";
    print F "\n";
    print F "set terminal $format size 1024,1024\n";
    print F "set output './$asm.originalLength-vs-expectedLength.lg.$format'\n";
    print F "plot [0:$maxReadLen] [0:$maxReadLen] './$asm.original-expected-corrected-length.dat' using 2:3 title 'original (x) vs expected (y)'\n";
    print F "set terminal $format size 256,256\n";
    print F "set output './$asm.originalLength-vs-expectedLength.sm.$format'\n";
    print F "replot\n";
    print F "\n";
    print F "set title 'original read length vs sum of corrected read lengths'\n";
    print F "set xlabel 'original read length'\n";
    print F "set ylabel 'sum of corrected read lengths'\n";
    print F "\n";
    print F "set terminal $format size 1024,1024\n";
    print F "set output './$asm.originalLength-vs-correctedLength.lg.$format'\n";
    print F "plot [0:$maxReadLen] [0:$maxReadLen] './$asm.original-expected-corrected-length.dat' using 2:4 title 'original (x) vs corrected (y)'\n";
    print F "set terminal $format size 256,256\n";
    print F "set output './$asm.originalLength-vs-correctedLength.sm.$format'\n";
    print F "replot\n";
    print F "\n";
    print F "set title 'expected read length vs sum of corrected read lengths'\n";
    print F "set xlabel 'expected read length'\n";
    print F "set ylabel 'sum of corrected read lengths'\n";
    print F "\n";
    print F "set terminal $format size 1024,1024\n";
    print F "set output './$asm.expectedLength-vs-correctedLength.lg.$format'\n";
    print F "plot [0:$maxReadLen] [0:$maxReadLen] './$asm.original-expected-corrected-length.dat' using 3:4 title 'expected (x) vs corrected (y)'\n";
    print F "set terminal $format size 256,256\n";
    print F "set output './$asm.expectedLength-vs-correctedLength.sm.$format'\n";
    print F "replot\n";
    close(F);

    if (runCommandSilently($path, "$gnuplot ./$asm.originalLength-vs-correctedLength.gp > /dev/null 2>&1", 0)) {
        print STDERR "--\n";
        print STDERR "-- WARNING: gnuplot failed; no plots will appear in HTML output.\n";
        print STDERR "--\n";
        print STDERR "----------------------------------------\n";
    }

    stashFile("$path/$asm.originalLength-vs-correctedLength.gp");
    stashFile("$asm.originalLength-vs-expectedLength.lg.$format");
    stashFile("$asm.originalLength-vs-expectedLength.sm.$format");
    stashFile("$asm.originalLength-vs-correctedLength.lg.$format");
    stashFile("$asm.originalLength-vs-correctedLength.sm.$format");
    stashFile("$asm.expectedLength-vs-correctedLength.lg.$format");
    stashFile("$asm.expectedLength-vs-correctedLength.sm.$format");

    #  Histograms of lengths, including the difference between expected and actual corrected (the
    #  other two difference plots weren't interesting; original-expected was basically all zero, and
    #  so original-actual was nearly the same as expected-actual.

    open(F, "> $path/$asm.length-histograms.gp") or caExit("", undef);
    print F "set title 'read length'\n";
    print F "set ylabel 'number of reads'\n";
    print F "set xlabel 'read length, bin width = 250'\n";
    print F "\n";
    print F "binwidth=250\n";
    print F "set boxwidth binwidth\n";
    print F "bin(x,width) = width*floor(x/width) + binwidth/2.0\n";
    print F "\n";
    print F "set terminal $format size 1024,1024\n";
    print F "set output './$asm.length-histograms.lg.$format'\n";
    print F "plot [1:$maxReadLen] [0:] \\\n";
    print F "  './$asm.original-expected-corrected-length.dat' using (bin(\$2,binwidth)):(1.0) smooth freq with boxes title 'original', \\\n";
    print F "  './$asm.original-expected-corrected-length.dat' using (bin(\$3,binwidth)):(1.0) smooth freq with boxes title 'expected', \\\n";
    print F "  './$asm.original-expected-corrected-length.dat' using (bin(\$4,binwidth)):(1.0) smooth freq with boxes title 'corrected'\n";
    print F "set terminal $format size 256,256\n";
    print F "set output './$asm.length-histograms.sm.$format'\n";
    print F "replot\n";
    print F "\n";
    print F "set xlabel 'difference between expected and corrected read length, bin width = 250, min=$minDiff, max=$maxDiff'\n";
    print F "\n";
    print F "set terminal $format size 1024,1024\n";
    print F "set output './$asm.length-difference-histograms.lg.$format'\n";
    print F "plot [$minDiff:$maxDiff] [0:] \\\n";
    print F "  './$asm.original-expected-corrected-length.dat' using (bin(\$7,binwidth)):(1.0) smooth freq with boxes title 'expected - corrected'\n";
    print F "set terminal $format size 256,256\n";
    print F "set output './$asm.length-difference-histograms.sm.$format'\n";
    print F "replot\n";
    close(F);

    if (runCommandSilently($path, "$gnuplot ./$asm.length-histograms.gp > /dev/null 2>&1", 0)) {
        print STDERR "--\n";
        print STDERR "-- WARNING: gnuplot failed; no plots will appear in HTML output.\n";
        print STDERR "--\n";
        print STDERR "----------------------------------------\n";
    }

    stashFile("$path/$asm.length-histograms.gp");
    stashFile("$asm.length-histograms.lg.$format");
    stashFile("$asm.length-histograms.sm.$format");
    stashFile("$asm.length-difference-histograms.lg.$format");
    stashFile("$asm.length-difference-histograms.sm.$format");

    #  Now that all outputs are (re)written, cleanup the job outputs.

    print STDERR "--\n";
    print STDERR "-- Purging correctReads output after merging to final output file.\n";
    print STDERR "-- Purging disabled by saveReadCorrections=true.\n"  if (getGlobal("saveReadCorrections") == 1);

    my $Nsuccess = 0;
    my $Nerr     = 0;
    my $Nfasta   = 0;
    my $Nlog     = 0;

    if (getGlobal("saveReadCorrections") != 1) {
        open(F, "< $path/corjob.files") or caExit("can't open '$path/corjob.files' for reading: $!", undef);
        while (<F>) {
            chomp;

            if (m/^(.*)\/results\/0*(\d+).fasta$/) {
                my $ID6 = substr("00000" . $2, -6);
                my $ID4 = substr("000"   . $2, -4);
                my $ID0 = $2;

                if (-e "$1/results/$ID4.dump.success") {
                    $Nsuccess++;
                    unlink "$1/results/$ID4.dump.success";
                }
                if (-e "$1/results/$ID4.err") {
                    $Nerr++;
                    unlink "$1/results/$ID4.err";
                }
                if (-e "$1/results/$ID4.fasta") {
                    $Nfasta++;
                    unlink "$1/results/$ID4.fasta";
                }

                if (-e "$1/correctReads.$ID6.out") {
                    $Nlog++;
                    unlink "$1/correctReads.$ID6.out";
                }
                if (-e "$1/correctReads.$ID0.out") {
                    $Nlog++;
                    unlink "$1/correctReads.$ID0.out";
                }

            } else {
                caExit("unknown correctReads job name '$_'\n", undef);
            }
        }
        close(F);

        print STDERR "-- Purged $Nsuccess .dump.success sentinels.\n"   if ($Nsuccess > 0);
        print STDERR "-- Purged $Nfasta .fasta outputs.\n"              if ($Nfasta > 0);
        print STDERR "-- Purged $Nerr .err outputs.\n"                  if ($Nerr > 0);
        print STDERR "-- Purged $Nlog .out job log outputs.\n"          if ($Nlog > 0);
    }

    remove_tree("correction/$asm.ovlStore")   if (getGlobal("saveOverlaps") eq "0");

  finishStage:
    emitStage($asm, "cor-dumpCorrectedReads");
    buildHTML($asm, "cor");

  allDone:
    print STDERR "--\n";
    print STDERR "-- Corrected reads saved in '", sequenceFileExists("$asm.correctedReads"), "'.\n";

    stopAfter("readCorrection");
}
