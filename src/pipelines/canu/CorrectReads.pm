
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
@EXPORT = qw(buildCorrectionLayouts generateCorrectedReads dumpCorrectedReads);

use strict;

use File::Path qw(make_path remove_tree);

use canu::Defaults;
use canu::Execution;
use canu::Gatekeeper;
use canu::HTML;

#  Returns a coverage:
#    If $cov not defined, default to desired output coverage * 1.0.
#    Otherwise, if defined but ends in an 'x', that's desired output coverage * whatever
#    Otherwise, the coverage is as defined.
#
sub getCorCov ($$$) {
    my $wrk     = shift @_;  #  Local work directory
    my $asm     = shift @_;
    my $typ     = shift @_;
    my $cov     = getGlobal("corMaxEvidenceCoverage$typ");

    my $exp = getExpectedCoverage($wrk, $asm);
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
sub getCorErrorRate ($$) {
    my $wrk     = shift @_;  #  Local work directory
    my $asm     = shift @_;
    my $bin     = getBinDirectory();
    my $erate   = getGlobal("corErrorRate");

    if (defined($erate)) {
        print STDERR "-- Using overlaps no worse than $erate fraction error for correcting reads (from corErrorRate parameter).\n";
        return($erate);
    }

    my $numPacBioRaw         = 0;
    my $numPacBioCorrected   = 0;
    my $numNanoporeRaw       = 0;
    my $numNanoporeCorrected = 0;

    open(L, "< $wrk/$asm.gkpStore/libraries.txt") or caExit("can't open '$wrk/$asm.gkpStore/libraries.txt' for reading: $!", undef);
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

    return($erate);
}



#  Return the number of jobs for 'falcon', 'falconpipe' or 'utgcns'
#
sub computeNumberOfCorrectionJobs ($$) {
    my $wrk     = shift @_;  #  Local work directory
    my $asm     = shift @_;
    my $nJobs   = 0;
    my $nPerJob = 0;

    if (getGlobal("corConsensus") eq "falcon" && -e "$wrk/2-correction/correction_inputs" ) {
        open(F, "ls $wrk/2-correction/correction_inputs/ |") or caExit("can't find list of correction_inputs: $!", undef);
        while (<F>) {
            $nJobs++  if (m/^\d\d\d\d$/);
        }
        close(F);

        return($nJobs, undef);
    }

    if ((getGlobal("corConsensus") eq "utgcns") ||
        (getGlobal("corConsensus") eq "falconpipe") ||
        (getGlobal("corConsensus") eq "falcon")) {
        my $nPart    = getGlobal("corPartitions");
        my $nReads   = getNumberOfReadsInStore($wrk, $asm);

        $nPerJob     = int($nReads / $nPart + 1);
        $nPerJob     = getGlobal("corPartitionMin")  if ($nPerJob < getGlobal("corPartitionMin"));

        for (my $j=1; $j<=$nReads; $j += $nPerJob) {  #  We could just divide, except for rounding issues....
            $nJobs++;
        }
    }

    return($nJobs, $nPerJob);
}


#  Generate a corStore, dump files for falcon to process, generate a script to run falcon.
#
sub buildCorrectionLayouts_direct ($$) {
    my $wrk  = shift @_;  #  Local work directory
    my $asm  = shift @_;
    my $bin  = getBinDirectory();
    my $cmd;

    my $path = "$wrk/2-correction";

    #  Outer level buildCorrectionLayouts() ensures the task is not finished.

    my $maxCov = getCorCov($wrk, $asm, "Local");

    if (! -e "$wrk/$asm.corStore") {
        $cmd  = "$bin/generateCorrectionLayouts \\\n";
        $cmd .= "  -rl $path/$asm.readsToCorrect \\\n"                 if (-e "$path/$asm.readsToCorrect");
        $cmd .= "  -G $wrk/$asm.gkpStore \\\n";
        $cmd .= "  -O $wrk/$asm.ovlStore \\\n";
        $cmd .= "  -S $path/$asm.globalScores \\\n"                    if (-e "$path/$asm.globalScores");
        $cmd .= "  -T $wrk/$asm.corStore.WORKING \\\n";
        $cmd .= "  -L " . getGlobal("corMinEvidenceLength") . " \\\n"  if (defined(getGlobal("corMinEvidenceLength")));
        $cmd .= "  -E " . getGlobal("corMaxEvidenceErate")  . " \\\n"  if (defined(getGlobal("corMaxEvidenceErate")));
        $cmd .= "  -C $maxCov \\\n"                                    if (defined($maxCov));
        $cmd .= "  -legacy \\\n"                                       if (defined(getGlobal("corLegacyFilter")));
        $cmd .= "> $wrk/$asm.corStore.err 2>&1";

        if (runCommand($wrk, $cmd)) {
            caExit("failed to generate layouts for correction", "$wrk/$asm.corStore.err");
        }

        rename "$wrk/$asm.corStore.WORKING", "$wrk/$asm.corStore";
    }

    # first we call this function to compute partioning
    my ($jobs, $nPer) = computeNumberOfCorrectionJobs($wrk, $asm);

    make_path("$path/correction_inputs")  if (! -d "$path/correction_inputs");
    make_path("$path/correction_outputs")  if (! -d "$path/correction_outputs");

    if (getGlobal("corConsensus") eq "falcon") {
        $cmd  = "$bin/createFalconSenseInputs \\\n";
        $cmd .= "  -G $wrk/$asm.gkpStore \\\n";
        $cmd .= "  -T $wrk/$asm.corStore 1 \\\n";
        $cmd .= "  -o $path/correction_inputs/ \\\n";
        $cmd .= "  -p " . $jobs . " \\\n";
        $cmd .= "> $path/correction_inputs.err 2>&1";

        if (runCommand($wrk, $cmd)) {
            caExit("failed to generate falcon inputs", "$path/correction_inputs.err");
        }
    }

    if (getGlobal("corConsensus") eq "falcon") {
       #the second call will confirm we have the proper number of output files and set jobs
       ($jobs, $nPer) = computeNumberOfCorrectionJobs($wrk, $asm);
    }

    #getAllowedResources("", "cor");

    open(F, "> $path/correctReads.sh") or caExit("can't open '$path/correctReads.sh'", undef);

    print F "#!" . getGlobal("shell") . "\n";
    print F "\n";
    print F "jobid=\$" . getGlobal("gridEngineTaskID") . "\n";
    print F "if [ x\$jobid = x -o x\$jobid = xundefined -o x\$jobid = x0 ]; then\n";
    print F "  jobid=\$1\n";
    print F "fi\n";
    print F "if [ x\$jobid = x ]; then\n";
    print F "  echo Error: I need " . getGlobal("gridEngineTaskID") . " set, or a job index on the command line.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "\n";
    print F "if [ \$jobid -gt $jobs ]; then\n";
    print F "  echo Error: Only $jobs partitions, you asked for \$jobid.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "\n";
    if (getGlobal("corConsensus") eq "utgcns") {
       print F "bgn=`expr \\( \$jobid - 1 \\) \\* $nPer`\n";
       print F "end=`expr \\( \$jobid + 0 \\) \\* $nPer`\n";
       print F "\n";
    }
    print F "jobid=`printf %04d \$jobid`\n";
    print F "\n";
    print F "if [ -e \"$path/correction_outputs/\$jobid.fasta\" ] ; then\n";
    print F "  echo Job finished successfully.\n";
    print F "  exit 0\n";
    print F "fi\n";
    print F "\n";
    print F "if [ ! -d \"$path/correction_outputs\" ] ; then\n";
    print F "  mkdir -p \"$path/correction_outputs\"\n";
    print F "fi\n";
    print F "\n";

    print F getBinDirectoryShellCode();

    my $erate  = getCorErrorRate($wrk, $asm);
    my $minidt = 1 - $erate;

    #  UTGCNS for correction is writing FASTQ, but needs to write FASTA.  The names below were changed to fasta preemptively.

    if (getGlobal("corConsensus") eq "utgcns") {
        caExit("UTGCNS for correction is writing FASTQ, but needs to write FASTA", undef);
        print F "\n";
        print F "\$bin/utgcns \\\n";
        print F "  -u \$bgn-\$end \\\n";
        print F "  -e $erate \\\n";
        print F "  -G $wrk/$asm.gkpStore \\\n";
        print F "  -T $wrk/$asm.corStore 1 . \\\n";
        print F "  -O $path/correction_outputs/\$jobid.cns.WORKING \\\n";
        print F "  -L $path/correction_outputs/\$jobid.layout.WORKING \\\n";
        print F "  -F $path/correction_outputs/\$jobid.fasta.WORKING \\\n";
        print F "&& \\\n";
        print F "mv $path/correction_outputs/\$jobid.cns.WORKING $path/correction_outputs/\$jobid.cns \\\n";
        print F "&& \\\n";
        print F "mv $path/correction_outputs/\$jobid.layout.WORKING $path/correction_outputs/\$jobid.layout \\\n";
        print F "&& \\\n";
        print F "mv $path/correction_outputs/\$jobid.fasta.WORKING $path/correction_outputs/\$jobid.fasta \\\n";
        print F "\n";
    }

    if (getGlobal("corConsensus") eq "falcon") {
        print F "\n";
        print F getGlobal("falconSense") . " \\\n"  if ( defined(getGlobal("falconSense")));
        print F "\$bin/falcon_sense \\\n"           if (!defined(getGlobal("falconSense")));
        print F "  --min_idt $minidt \\\n";
        print F "  --min_len " . getGlobal("minReadLength") . "\\\n";
        print F "  --max_read_len " . 2 * getMaxReadInStore($wrk, $asm) . "\\\n";
        print F "  --min_ovl_len " . getGlobal("minOverlapLength") . "\\\n";
        print F "  --min_cov " . getGlobal("corMinCoverage") . " \\\n";
        print F "  --n_core " . getGlobal("corThreads") . " \\\n";
        print F "  < $path/correction_inputs/\$jobid \\\n";
        print F "  > $path/correction_outputs/\$jobid.fasta.WORKING \\\n";
        print F " 2> $path/correction_outputs/\$jobid.err \\\n";
        print F "&& \\\n";
        print F "mv $path/correction_outputs/\$jobid.fasta.WORKING $path/correction_outputs/\$jobid.fasta \\\n";
    }

    close(F);

  finishStage:
    ;
  allDone:
}





#  For falcon_sense, using a pipe and no intermediate files
#
sub buildCorrectionLayouts_piped ($$) {
    my $wrk  = shift @_;  #  Local work directory
    my $asm  = shift @_;
    my $bin  = getBinDirectory();
    my $cmd;
    my $path = "$wrk/2-correction";

    #  Outer level buildCorrectionLayouts() ensures the task is not finished.

    make_path("$path/correction_inputs")  if (! -d "$path/correction_inputs");
    make_path("$path/correction_outputs")  if (! -d "$path/correction_outputs");

    my ($nJobs, $nPerJob)  = computeNumberOfCorrectionJobs($wrk, $asm);  #  Does math based on number of reads and parameters.

    my $nReads             = getNumberOfReadsInStore($wrk, $asm);

    #getAllowedResources("", "cor");

    open(F, "> $path/correctReads.sh") or caExit("can't open '$path/correctReads.sh'", undef);

    print F "#!" . getGlobal("shell") . "\n";
    print F "\n";
    print F "jobid=\$" . getGlobal("gridEngineTaskID") . "\n";
    print F "if [ x\$jobid = x -o x\$jobid = xundefined -o x\$jobid = x0 ]; then\n";
    print F "  jobid=\$1\n";
    print F "fi\n";
    print F "if [ x\$jobid = x ]; then\n";
    print F "  echo Error: I need " . getGlobal("gridEngineTaskID") . " set, or a job index on the command line.\n";
    print F "  exit 1\n";
    print F "fi\n";
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
    print F "if [ -e \"$path/correction_outputs/\$jobid.fasta\" ] ; then\n";
    print F "  echo Job finished successfully.\n";
    print F "  exit 0\n";
    print F "fi\n";
    print F "\n";
    print F "if [ ! -d \"$path/correction_outputs\" ] ; then\n";
    print F "  mkdir -p \"$path/correction_outputs\"\n";
    print F "fi\n";
    print F "\n";
    print F getBinDirectoryShellCode();
    print F "\n";

    my $maxCov   = getCorCov($wrk, $asm, "Local");

    my $erate    = getCorErrorRate($wrk, $asm);
    my $minidt   = 1 - $erate;

    print F "\n";
    print F "if [ \"x\$BASH\" != \"x\" ] ; then\n";   #  Needs doublequotes, else shell doesn't expand $BASH
    print F "  set -o pipefail\n";
    print F "fi\n";
    print F "\n";
    print F "( \\\n";
    print F "\$bin/generateCorrectionLayouts -b \$bgn -e \$end \\\n";
    print F "  -rl $path/$asm.readsToCorrect \\\n"                 if (-e "$path/$asm.readsToCorrect");
    print F "  -G $wrk/$asm.gkpStore \\\n";
    print F "  -O $wrk/$asm.ovlStore \\\n";
    print F "  -S $path/$asm.globalScores \\\n"                    if (-e "$path/$asm.globalScores");
    print F "  -L " . getGlobal("corMinEvidenceLength") . " \\\n"  if (defined(getGlobal("corMinEvidenceLength")));
    print F "  -E " . getGlobal("corMaxEvidenceErate")  . " \\\n"  if (defined(getGlobal("corMaxEvidenceErate")));
    print F "  -C $maxCov \\\n"                                    if (defined($maxCov));
    print F "  -legacy \\\n"                                       if (defined(getGlobal("corLegacyFilter")));
    print F "  -F \\\n";
    print F "&& \\\n";
    print F "  touch $path/correction_outputs/\$jobid.dump.success \\\n";
    print F ") \\\n";
    print F "| \\\n";
    print F getGlobal("falconSense") . " \\\n"  if ( defined(getGlobal("falconSense")));
    print F "\$bin/falcon_sense \\\n"           if (!defined(getGlobal("falconSense")));
    print F "  --min_idt $minidt \\\n";
    print F "  --min_len " . getGlobal("minReadLength") . "\\\n";
    print F "  --max_read_len " . 2 * getMaxReadInStore($wrk, $asm) . "\\\n";
    print F "  --min_ovl_len " . getGlobal("minOverlapLength") . "\\\n";
    print F "  --min_cov " . getGlobal("corMinCoverage") . " \\\n";
    print F "  --n_core " . getGlobal("corThreads") . " \\\n";
    print F "  > $path/correction_outputs/\$jobid.fasta.WORKING \\\n";
    print F " 2> $path/correction_outputs/\$jobid.err \\\n";
    print F "&& \\\n";
    print F "mv $path/correction_outputs/\$jobid.fasta.WORKING $path/correction_outputs/\$jobid.fasta \\\n";
    print F "\n";
    print F "if [ ! -e \"$path/correction_outputs/\$jobid.dump.success\" ] ; then\n";
    print F "  echo Read layout generation failed.\n";
    print F "  mv $path/correction_outputs/\$jobid.fasta $path/correction_outputs/\$jobid.fasta.INCOMPLETE\n";
    print F "fi\n";
    print F "\n";
    print F "exit 0\n";

    close(F);

  finishStage:
    ;
  allDone:
}



sub lengthStats (@) {
    my @v = sort { $b <=> $a } @_;

    my $total  = 0;
    my $mean   = 0;
    my $n50    = 0;

    foreach my $v (@v) {
        $total += $v;
        $n50    = $v  if ($total < getGlobal("genomeSize") / 2);
    }

    $mean = int($total / scalar(@v) + 0.5);

    return($mean, $n50);
}


sub quickFilter ($$$) {
    my $wrk  = shift @_;  #  Local work directory
    my $asm  = shift @_;
    my $minTotal = shift @_;
    my $bin  = getBinDirectory();
    my $path = "$wrk/2-correction";

    my $totCorLengthIn  = 0;
    my $totCorLengthOut = 0;
    my $minCorLength    = 0;

    open(O, "> $path/$asm.readsToCorrect.WORKING") or caExit("can't open '$path/$asm.readsToCorrect.WORKING' for writing: $!\n", undef);
    open(F, "$bin/gatekeeperDumpMetaData -G $wrk/$asm.gkpStore -reads | sort -T . -k3nr | ") or caExit("can't dump gatekeeper for read lengths: $!\n", undef);

    print O "read\toriginalLength\tcorrectedLength\n";

    while (<F>) {
        my @v = split '\s+', $_;

        $totCorLengthIn  += $v[2];
        $totCorLengthOut += $v[2];

        print O "$v[0]\t$v[2]\t0\n";

        if ($minTotal != 0 && $totCorLengthIn >= $minTotal) {
            $minCorLength = $v[2];
            last;
        }
    }
    close(F);
    close(O);

    rename "$path/$asm.readsToCorrect.WORKING", "$path/$asm.readsToCorrect";
}




sub expensiveFilter ($$) {
    my $wrk  = shift @_;  #  Local work directory
    my $asm  = shift @_;
    my $bin  = getBinDirectory();
    my $cmd;
    my $path = "$wrk/2-correction";

    my $maxCov = getCorCov($wrk, $asm, "Local");

    if (! -e "$path/$asm.estimate.log") {
        $cmd  = "$bin/generateCorrectionLayouts \\\n";
        $cmd .= "  -G $wrk/$asm.gkpStore \\\n";
        $cmd .= "  -O $wrk/$asm.ovlStore \\\n";
        $cmd .= "  -S $path/$asm.globalScores \\\n"                    if (-e "$path/$asm.globalScores");
        $cmd .= "  -L " . getGlobal("corMinEvidenceLength") . " \\\n"  if (defined(getGlobal("corMinEvidenceLength")));
        $cmd .= "  -E " . getGlobal("corMaxEvidenceErate")  . " \\\n"  if (defined(getGlobal("corMaxEvidenceErate")));
        $cmd .= "  -C $maxCov \\\n"                                    if (defined($maxCov));
        $cmd .= "  -legacy \\\n"                                       if (defined(getGlobal("corLegacyFilter")));
        $cmd .= "  -p $path/$asm.estimate";

        if (runCommand($wrk, $cmd)) {
            rename "$path/$asm.estimate.log", "$path/$asm.estimate.log.FAILED";
            caExit("failed to generate estimated lengths of corrected reads", "$wrk/$asm.corStore.err");
        }
    }

    if (runCommandSilently($path, "sort -T . -k4nr -k2nr < $path/$asm.estimate.log > $path/$asm.estimate.correctedLength.log", 1)) {
        caExit("failed to sort by corrected read length", undef);
    }

    if (runCommandSilently($path, "sort -T . -k2nr -k4nr < $path/$asm.estimate.log > $path/$asm.estimate.originalLength.log",  1)) {
        caExit("failed to sort by original read length", undef);
    }

    my $totRawLengthIn    = 0;  #  Bases in raw reads we correct
    my $totRawLengthOut   = 0;  #  Expected bases in corrected reads
    my $minRawLength      = 0;

    my $totCorLengthIn    = 0;  #  Bases in raw reads we correct
    my $totCorLengthOut   = 0;  #  Expected bases in corrected reads
    my $minCorLength      = 0;

    my $minTotal          = getGlobal("genomeSize") * getGlobal("corOutCoverage");

    #  Lists of reads to correct if we use the raw length or the corrected length as a filter.

    my %rawReads;
    my %corReads;
    my %corReadLen;

    #  The expected length of the corrected reads, based on the filter

    my @corLengthRawFilter;
    my @corLengthCorFilter;

    #  Filter!

    open(F, "< $path/$asm.estimate.originalLength.log");
    while (<F>) {
        my @v = split '\s+', $_;
        next if ($v[0] eq "read");

        $corReadLen{$v[0]} = $v[3];
    }
    close(F);

    open(F, "< $path/$asm.estimate.originalLength.log");
    while (<F>) {
        my @v = split '\s+', $_;
        next if ($v[0] eq "read");

        $totRawLengthIn  += $v[1];
        $totRawLengthOut += $v[3];

        #print O "$v[0]\t$v[1]\t$v[3]\n";

        $rawReads{$v[0]} = 1;

        push @corLengthRawFilter, $v[3];

        if ($totRawLengthIn >= $minTotal) {  #  Compare against raw bases
            $minRawLength = $v[1];
            last;
        }
    }

    open(F, "< $path/$asm.estimate.correctedLength.log");
    open(O, "| sort -T . -k1n > $path/$asm.readsToCorrect.WORKING") or caExit("can't open sort -k1n > '$path/$asm.readsToCorrect.WORKING' for writing: $!\n", undef);

    print O "read\toriginalLength\tcorrectedLength\n";

    while (<F>) {
        my @v = split '\s+', $_;
        next if ($v[0] eq "read");

        $totCorLengthIn  += $v[1];
        $totCorLengthOut += $v[3];

        print O "$v[0]\t$v[1]\t$v[3]\n";

        $corReads{$v[0]} = 1;

        push @corLengthCorFilter, $v[3];

        if ($totCorLengthOut >= $minTotal) {  #  Compare against corrected bases
            $minCorLength = $v[3];
            last;
        }
    }
    close(O);
    close(F);

    rename "$path/$asm.readsToCorrect.WORKING", "$path/$asm.readsToCorrect";


    #  Generate true/false positive/negative lists.

    open(F, "< $path/$asm.estimate.correctedLength.log") or die;

    open(TN, "> $path/$asm.estimate.tn.log") or die;
    open(FN, "> $path/$asm.estimate.fn.log") or die;
    open(FP, "> $path/$asm.estimate.fp.log") or die;
    open(TP, "> $path/$asm.estimate.tp.log") or die;

    my ($tnReads, $tnBasesR, $tnBasesC, $tnBasesRave, $tnBasesCave) = (0, 0, 0, 0, 0 );
    my ($fnReads, $fnBasesR, $fnBasesC, $fnBasesRave, $fnBasesCave) = (0, 0, 0, 0, 0 );
    my ($fpReads, $tpBasesR, $tpBasesC, $tpBasesRave, $tpBasesCave) = (0, 0, 0, 0, 0 );
    my ($tpReads, $fpBasesR, $fpBasesC, $fpBasesRave, $fpBasesCave) = (0, 0, 0, 0, 0 );

    while (<F>) {
        my @v = split '\s+', $_;
        next if ($v[0] eq "read");

        my $er = exists($rawReads{$v[0]});
        my $ec = exists($corReads{$v[0]});

        if (($er == 0) && ($ec == 0)) {  print TN $_;  $tnReads++;  $tnBasesR += $v[1];  $tnBasesC += $corReadLen{$v[0]};  }  #  True negative, yay!
        if (($er == 0) && ($ec == 1)) {  print FN $_;  $fnReads++;  $fnBasesR += $v[1];  $fnBasesC += $corReadLen{$v[0]};  }  #  False negative.  Bad.
        if (($er == 1) && ($ec == 0)) {  print FP $_;  $fpReads++;  $fpBasesR += $v[1];  $fpBasesC += $corReadLen{$v[0]};  }  #  False positive.  Bad.
        if (($er == 1) && ($ec == 1)) {  print TP $_;  $tpReads++;  $tpBasesR += $v[1];  $tpBasesC += $corReadLen{$v[0]};  }  #  True positive, yay!
    }

    close(TP);
    close(FP);
    close(FN);
    close(TN);

    close(F);

    $tnBasesRave = $tnBasesR / $tnReads  if ($tnReads > 0);
    $tnBasesCave = $tnBasesC / $tnReads  if ($tnReads > 0);

    $fnBasesRave = $fnBasesR / $fnReads  if ($fnReads > 0);
    $fnBasesCave = $fnBasesC / $fnReads  if ($fnReads > 0);

    $fpBasesRave = $fpBasesR / $fpReads  if ($fpReads > 0);
    $fpBasesCave = $fpBasesC / $fpReads  if ($fpReads > 0);

    $tpBasesRave = $tpBasesR / $tpReads  if ($tpReads > 0);
    $tpBasesCave = $tpBasesC / $tpReads  if ($tpReads > 0);

    #  Dump a summary of the filter

    my ($rawFilterMean, $rawFilterN50) = lengthStats(@corLengthRawFilter);
    my ($corFilterMean, $corFilterN50) = lengthStats(@corLengthCorFilter);

    open(F, "> $path/$asm.readsToCorrect.summary") or caExit("can't open '$path/$asm.readsToCorrect.summary' for writing: $!\n", undef);
    print F "Corrected read length filter:\n";
    print F "\n";
    print F "  nReads  ", scalar(keys %corReads), "\n";
    print F "  nBases  $totCorLengthIn (input bases)\n";
    print F "  nBases  $totCorLengthOut (corrected bases)\n";
    print F "  Mean    $corFilterMean\n";
    print F "  N50     $corFilterN50\n";
    print F "\n";
    print F "Raw read length filter:\n";
    print F "\n";
    print F "  nReads  ", scalar(keys %rawReads), "\n";
    print F "  nBases  $totRawLengthIn (input bases)\n";
    print F "  nBases  $totRawLengthOut (corrected bases)\n";
    print F "  Mean    $rawFilterMean\n";
    print F "  N50     $rawFilterN50\n";
    print F "\n";
    printf F "TN %9d reads %13d raw bases (%6d ave) %13d corrected bases (%6d ave)\n", $tnReads, $tnBasesR, $tnBasesRave, $tnBasesC, $tnBasesCave;
    printf F "FN %9d reads %13d raw bases (%6d ave) %13d corrected bases (%6d ave)\n", $fnReads, $fnBasesR, $fnBasesRave, $fnBasesC, $fnBasesCave;
    printf F "FP %9d reads %13d raw bases (%6d ave) %13d corrected bases (%6d ave)\n", $fpReads, $fpBasesR, $fpBasesRave, $fpBasesC, $fpBasesCave;
    printf F "TP %9d reads %13d raw bases (%6d ave) %13d corrected bases (%6d ave)\n", $tpReads, $tpBasesR, $tpBasesRave, $tpBasesC, $tpBasesCave;
    close(F);

    #  Plot a scatter plot of the original vs the expected corrected read lengths.  Early versions
    #  also plotted the sorted length vs the other length, but those were not interesting.

    if (! -e "$path/$asm.estimate.original-x-correctedLength.png") {
        open(F, "> $path/$asm.estimate.original-x-correctedLength.gp");
        print F "set title 'original length (x) vs corrected length (y)'\n";
        print F "set xlabel 'original read length'\n";
        print F "set ylabel 'corrected read length (expected)'\n";
        print F "set pointsize 0.25\n";
        print F "\n";
        print F "set terminal png size 1024,1024\n";
        print F "set output '$path/$asm.estimate.original-x-corrected.lg.png'\n";
        print F "plot '$path/$asm.estimate.tn.log' using 2:4 title 'tn', \\\n";
        print F "     '$path/$asm.estimate.fn.log' using 2:4 title 'fn', \\\n";
        print F "     '$path/$asm.estimate.fp.log' using 2:4 title 'fp', \\\n";
        print F "     '$path/$asm.estimate.tp.log' using 2:4 title 'tp'\n";
        print F "set terminal png size 256,256\n";
        print F "set output '$path/$asm.estimate.original-x-corrected.sm.png'\n";
        print F "replot\n";
        close(F);

        if (runCommandSilently($path, "gnuplot < $path/$asm.estimate.original-x-correctedLength.gp > /dev/null 2>&1", 0)) {
            print STDERR "--\n";
            print STDERR "-- WARNING: gnuplot failed; no plots will appear in HTML output.\n";
            print STDERR "--\n";
            print STDERR "----------------------------------------\n";
        }
    }
}



sub buildCorrectionLayouts ($$) {
    my $WRK     = shift @_;           #  Root work directory (the -d option to canu)
    my $wrk     = "$WRK/correction";  #  Local work directory
    my $asm     = shift @_;
    my $bin     = getBinDirectory();
    my $cmd;

    my $path = "$wrk/2-correction";

    #  All we do here is decide if the job is finished, and delegate
    #  to the correct function if not finished.

    #  But first we analyze the layouts to pick a length cutoff - or more precisely, to pick a set
    #  of reads to correct.
    #
    #  If one doesn't want to compute the (expensive) stats, one can just correct more reads.  How many more to correct
    #  is dependent on the particular reads.  An extra 1x to 5x seems reasonable.

    #  A one-pass algorithm would write layouts to tigStore, computing stats as before, then delete
    #  tigs that shouldn't be corrected.  I suspect this will be slower.

    goto allDone   if (skipStage($WRK, $asm, "cor-buildCorrectionLayouts") == 1);
    goto allDone   if (sequenceFileExists("$WRK/$asm.correctedReads"));    #  Output exists
    goto allDone   if (-e "$path/cnsjob.files");                           #  Jobs all finished
    goto allDone   if (-e "$path/correctReads.sh");                        #  Jobs created

    make_path("$path")  if (! -d "$path");

    #  This will eventually get rolled into overlap store creation.  Generate a list of scores for
    #  'global' overlap filtering.

    if (! -e "$path/$asm.globalScores") {
        my $maxCov = getCorCov($wrk, $asm, "Global");
        my $minLen = (defined(getGlobal("corMinEvidenceLength"))) ? getGlobal("corMinEvidenceLength") : 0;

        $cmd  = "$bin/filterCorrectionOverlaps \\\n";
        $cmd .= "  -G $wrk/$asm.gkpStore \\\n";
        $cmd .= "  -O $wrk/$asm.ovlStore \\\n";
        $cmd .= "  -S $path/$asm.globalScores \\\n";
        $cmd .= "  -c $maxCov \\\n";
        $cmd .= "  -l $minLen \\\n";
        $cmd .= "  -e " . getGlobal("corMaxEvidenceErate")  . " \\\n"  if (defined(getGlobal("corMaxEvidenceErate")));
        $cmd .= "  -legacy \\\n"                                       if (defined(getGlobal("corLegacyFilter")));
        $cmd .= "> $path/$asm.globalScores.err 2>&1";

        if (runCommand($path, $cmd)) {
            caExit("failed to globally filter overlaps for correction", "$path/$asm.globalScores.err");
        }

        unlink "$path/$asm.globalScores.err";
    }

    #  For 'quick' filtering, but more reads to correct, sort the reads by length, and correct the
    #  longest Nx of reads.
    #
    #  For 'expensive' filtering, but fewer reads to correct, first estimate the corrected lengths
    #  (requires a pass through the overlaps), then pull out the longest Nx of corrected reads.
    #
    #  Both are required to create a file $asm.readsToCorrect, containing a list of IDs to correct.

    if      (getGlobal("corFilter") eq "quick") {
        quickFilter($wrk, $asm, (getGlobal("genomeSize") * getGlobal("corOutCoverage")));

    } elsif (getGlobal("corFilter") eq "expensive") {
        expensiveFilter($wrk, $asm);

    } elsif (getGlobal("corFilter") eq "none" ) {
        quickFilter($wrk, $asm, 0);

    } else {
        caFailure("unknown corFilter '" . getGlobal("corFilter") . "'", undef);
    }

    #  Set the minimum coverage for a corrected read based on coverage in input reads.

    if (!defined(getGlobal("corMinCoverage"))) {
        my $cov = getExpectedCoverage($wrk, $asm);

        setGlobal("corMinCoverage", 4);
        setGlobal("corMinCoverage", 4)   if ($cov <  60);
        setGlobal("corMinCoverage", 0)   if ($cov <= 20);

        print STDERR "-- Set corMinCoverage=", getGlobal("corMinCoverage"), " based on read coverage of $cov.\n";
    }

    caExit("failed to create list of reads to correct", undef)  if (! -e "$path/$asm.readsToCorrect");

    buildCorrectionLayouts_direct($wrk, $asm)      if (getGlobal("corConsensus") eq "utgcns");
    buildCorrectionLayouts_direct($wrk, $asm)      if (getGlobal("corConsensus") eq "falcon");
    buildCorrectionLayouts_piped($wrk, $asm)       if (getGlobal("corConsensus") eq "falconpipe");

  finishStage:
    emitStage($WRK, $asm, "cor-buildCorrectionLayouts");
    buildHTML($WRK, $asm, "cor");

  allDone:
}




sub generateCorrectedReads ($$) {
    my $WRK     = shift @_;           #  Root work directory (the -d option to canu)
    my $wrk     = "$WRK/correction";  #  Local work directory
    my $asm     = shift @_;
    my $attempt = getGlobal("canuIteration");
    my $bin     = getBinDirectory();

    my $path    = "$wrk/2-correction";

    goto allDone   if (skipStage($WRK, $asm, "cor-generateCorrectedReads", $attempt) == 1);
    goto allDone   if (sequenceFileExists("$WRK/$asm.correctedReads"));

    #  Figure out if all the tasks finished correctly.

    my ($jobs, undef) = computeNumberOfCorrectionJobs($wrk, $asm);

    my $currentJobID = "0001";
    my @successJobs;
    my @failedJobs;
    my $failureMessage = "";

    for (my $job=1; $job <= $jobs; $job++) {
        if (-e "$path/correction_outputs/$currentJobID.fasta") {
            push @successJobs, "$path/correction_outputs/$currentJobID.fasta\n";

        } else {
            $failureMessage .= "--   job $path/correction_outputs/$currentJobID.fasta FAILED.\n";
            push @failedJobs, $job;
        }

        $currentJobID++;
    }

    #  Failed jobs, retry.

    if (scalar(@failedJobs) > 0) {

        #  If not the first attempt, report the jobs that failed, and that we're recomputing.

        if ($attempt > 1) {
            print STDERR "--\n";
            print STDERR "-- ", scalar(@failedJobs), " read correction jobs failed:\n";
            print STDERR $failureMessage;
        }

        #  If too many attempts, give up.

        if ($attempt > getGlobal("canuIterationMax")) {
            caExit("failed to generate corrected reads.  Made " . ($attempt-1) . " attempts, jobs still failed", undef);
        }

        #  Otherwise, run some jobs.

        print STDERR "-- generate corrected reads attempt $attempt begins with ", scalar(@successJobs), " finished, and ", scalar(@failedJobs), " to compute.\n";

        emitStage($WRK, $asm, "cor-generateCorrectedReads", $attempt);
        buildHTML($WRK, $asm, "cor");

        submitOrRunParallelJob($WRK, $asm, "cor", $path, "correctReads", @failedJobs);
        return;
    }

  finishStage:
    print STDERR "-- Found ", scalar(@successJobs), " read correction output files.\n";

    open(L, "> $path/corjob.files") or caExit("failed to open '$path/corjob.files'", undef);
    print L @successJobs;
    close(L);

    setGlobal("canuIteration", 0);
    emitStage($WRK, $asm, "cor-generateCorrectedReads");
    buildHTML($WRK, $asm, "cor");
    stopAfter("readCorrection");

  allDone:
}



sub dumpCorrectedReads ($$) {
    my $WRK     = shift @_;           #  Root work directory (the -d option to canu)
    my $wrk     = "$WRK/correction";  #  Local work directory
    my $asm     = shift @_;
    my $bin     = getBinDirectory();

    my $path = "$wrk/2-correction";

    goto allDone   if (skipStage($WRK, $asm, "cor-dumpCorrectedReads") == 1);
    goto allDone   if (sequenceFileExists("$WRK/$asm.correctedReads"));

    print STDERR "-- Concatenating correctReads output.\n";

    my $files = 0;
    my $reads = 0;

    open(F, "< $path/corjob.files")                           or caExit("can't open '$path/corjob.files' for reading: $!", undef);
    open(N, "< $wrk/$asm.gkpStore/readNames.txt")             or caExit("can't open '$wrk/$asm.gkpStore/readNames.txt' for reading: $!", undef);
    open(O, "| gzip -1c > $WRK/$asm.correctedReads.fasta.gz") or caExit("can't open '$WRK/$asm.correctedReads.fasta.gz' for writing: $!", undef);
    #open(O, "> $WRK/$asm.correctedReads.fasta")               or caExit("can't open '$WRK/$asm.correctedReads.fasta' for writing: $!", undef);
    open(L, "> $WRK/$asm.correctedReads.length")              or caExit("can't open '$WRK/$asm.correctedReads.length' for writing: $!", undef);

    while (<F>) {
        chomp;
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
                $h = ">$name iid=${rid}_${pid}";
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

    open(R, "< $wrk/2-correction/$asm.readsToCorrect") or caExit("can't open '$wrk/2-correction/$asm.readsToCorrect' for reading: $!", undef);
    open(L, "< $WRK/$asm.correctedReads.length")       or caExit("can't open '$WRK./$asm.correctedReads.length' for reading: $!", undef);

    open(S, "> $wrk/2-correction/$asm.original-expected-corrected-length.dat") or caExit("", undef);

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

    #  Write a summary of the corrections.

    open(F, "> $wrk/2-correction/$asm.correction.summary") or caExit("", undef);
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

    #  Scatterplot of lengths.

    open(F, "> $wrk/2-correction/$asm.originalLength-vs-correctedLength.gp") or caExit("", undef);
    print F "\n";
    print F "set pointsize 0.25\n";
    print F "\n";
    print F "set title 'original read length vs expected corrected read length'\n";
    print F "set xlabel 'original read length'\n";
    print F "set ylabel 'expected corrected read length'\n";
    print F "\n";
    print F "set terminal png size 1024,1024\n";
    print F "set output '$wrk/2-correction/$asm.originalLength-vs-expectedLength.lg.png'\n";
    print F "plot [0:$maxReadLen] [0:$maxReadLen] '$wrk/2-correction/$asm.original-expected-corrected-length.dat' using 2:3 title 'original (x) vs expected (y)'\n";
    print F "set terminal png size 256,256\n";
    print F "set output '$wrk/2-correction/$asm.originalLength-vs-expectedLength.sm.png'\n";
    print F "replot\n";
    print F "\n";
    print F "set title 'original read length vs sum of corrected read lengths'\n";
    print F "set xlabel 'original read length'\n";
    print F "set ylabel 'sum of corrected read lengths'\n";
    print F "\n";
    print F "set terminal png size 1024,1024\n";
    print F "set output '$wrk/2-correction/$asm.originalLength-vs-correctedLength.lg.png'\n";
    print F "plot [0:$maxReadLen] [0:$maxReadLen] '$wrk/2-correction/$asm.original-expected-corrected-length.dat' using 2:4 title 'original (x) vs corrected (y)'\n";
    print F "set terminal png size 256,256\n";
    print F "set output '$wrk/2-correction/$asm.originalLength-vs-correctedLength.sm.png'\n";
    print F "replot\n";
    print F "\n";
    print F "set title 'expected read length vs sum of corrected read lengths'\n";
    print F "set xlabel 'expected read length'\n";
    print F "set ylabel 'sum of corrected read lengths'\n";
    print F "\n";
    print F "set terminal png size 1024,1024\n";
    print F "set output '$wrk/2-correction/$asm.expectedLength-vs-correctedLength.lg.png'\n";
    print F "plot [0:$maxReadLen] [0:$maxReadLen] '$wrk/2-correction/$asm.original-expected-corrected-length.dat' using 3:4 title 'expected (x) vs corrected (y)'\n";
    print F "set terminal png size 256,256\n";
    print F "set output '$wrk/2-correction/$asm.expectedLength-vs-correctedLength.sm.png'\n";
    print F "replot\n";
    close(F);

    if (runCommandSilently("$wrk/2-correction", "gnuplot $wrk/2-correction/$asm.originalLength-vs-correctedLength.gp > /dev/null 2>&1", 0)) {
        print STDERR "--\n";
        print STDERR "-- WARNING: gnuplot failed; no plots will appear in HTML output.\n";
        print STDERR "--\n";
        print STDERR "----------------------------------------\n";
    }

    #  Histograms of lengths, including the difference between expected and actual corrected (the
    #  other two difference plots weren't interesting; original-expected was basically all zero, and
    #  so original-actual was nearly the same as expected-actual.

    open(F, "> $wrk/2-correction/$asm.length-histograms.gp") or caExit("", undef);
    print F "set title 'read length'\n";
    print F "set ylabel 'number of reads'\n";
    print F "set xlabel 'read length, bin width = 250'\n";
    print F "\n";
    print F "binwidth=250\n";
    print F "set boxwidth binwidth\n";
    print F "bin(x,width) = width*floor(x/width) + binwidth/2.0\n";
    print F "\n";
    print F "set terminal png size 1024,1024\n";
    print F "set output '$wrk/2-correction/$asm.length-histograms.lg.png'\n";
    print F "plot [1:$maxReadLen] [0:] \\\n";
    print F "  '$wrk/2-correction/$asm.original-expected-corrected-length.dat' using (bin(\$2,binwidth)):(1.0) smooth freq with boxes title 'original', \\\n";
    print F "  '$wrk/2-correction/$asm.original-expected-corrected-length.dat' using (bin(\$3,binwidth)):(1.0) smooth freq with boxes title 'expected', \\\n";
    print F "  '$wrk/2-correction/$asm.original-expected-corrected-length.dat' using (bin(\$4,binwidth)):(1.0) smooth freq with boxes title 'corrected'\n";
    print F "set terminal png size 256,256\n";
    print F "set output '$wrk/2-correction/$asm.length-histograms.sm.png'\n";
    print F "replot\n";
    print F "\n";
    print F "set xlabel 'difference between expected and corrected read length, bin width = 250, min=$minDiff, max=$maxDiff'\n";
    print F "\n";
    print F "set terminal png size 1024,1024\n";
    print F "set output '$wrk/2-correction/$asm.length-difference-histograms.lg.png'\n";
    print F "plot [$minDiff:$maxDiff] [0:] \\\n";
    print F "  '$wrk/2-correction/$asm.original-expected-corrected-length.dat' using (bin(\$7,binwidth)):(1.0) smooth freq with boxes title 'expected - corrected'\n";
    print F "set terminal png size 256,256\n";
    print F "set output '$wrk/2-correction/$asm.length-difference-histograms.sm.png'\n";
    print F "replot\n";
    close(F);

    if (runCommandSilently("$wrk/2-correction", "gnuplot $wrk/2-correction/$asm.length-histograms.gp > /dev/null 2>&1", 0)) {
        print STDERR "--\n";
        print STDERR "-- WARNING: gnuplot failed; no plots will appear in HTML output.\n";
        print STDERR "--\n";
        print STDERR "----------------------------------------\n";
    }

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

            if (m/^(.*)\/correction_outputs\/0*(\d+).fasta$/) {
                my $ID6 = substr("00000" . $2, -6);
                my $ID4 = substr("000"   . $2, -4);
                my $ID0 = $2;

                if (-e "$1/correction_outputs/$ID4.dump.success") {
                    $Nsuccess++;
                    unlink "$1/correction_outputs/$ID4.dump.success";
                }
                if (-e "$1/correction_outputs/$ID4.err") {
                    $Nerr++;
                    unlink "$1/correction_outputs/$ID4.err";
                }
                if (-e "$1/correction_outputs/$ID4.fasta") {
                    $Nfasta++;
                    unlink "$1/correction_outputs/$ID4.fasta";
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
    }

    print STDERR "-- Purged $Nsuccess .dump.success sentinels.\n"   if ($Nsuccess > 0);
    print STDERR "-- Purged $Nfasta .fasta outputs.\n"              if ($Nfasta > 0);
    print STDERR "-- Purged $Nerr .err outputs.\n"                  if ($Nerr > 0);
    print STDERR "-- Purged $Nlog .out job log outputs.\n"          if ($Nlog > 0);

  finishStage:
    emitStage($WRK, $asm, "cor-dumpCorrectedReads");
    buildHTML($WRK, $asm, "cor");

  allDone:
    print STDERR "--\n";
    print STDERR "-- Corrected reads saved in '", sequenceFileExists("$WRK/$asm.correctedReads"), "'.\n";
}
