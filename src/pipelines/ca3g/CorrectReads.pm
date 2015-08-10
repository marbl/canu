package ca3g::CorrectReads;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(buildCorrectionLayouts generateCorrectedReads dumpCorrectedReads);

use strict;

use File::Path qw(make_path remove_tree);

use ca3g::Defaults;
use ca3g::Execution;
use ca3g::Gatekeeper;

#  Returns a coverage:
#    If $cov not defined, default to desired output coverage * 1.0.
#    Otherwise, if defined but ends in an 'x', that's desired output coverage * whatever
#    Otherwise, the coverage is as defined.
#
sub getCorCov ($$$) {
    my $wrk = shift @_;
    my $asm = shift @_;
    my $typ = shift @_;
    my $cov = getGlobal("corMaxEvidenceCoverage$typ");

    my $exp = getExpectedCoverage($wrk, $asm);
    my $des = getGlobal("corOutCoverage");

    if (!defined($cov)) {
        $cov = $des;
    } elsif ($cov =~ m/(.*)x/) {
        $cov = int($des * $1);
    }

    return($cov);
}



#  Return the number of jobs for 'falcon', 'falconpipe' or 'utgcns'
#
sub computeNumberOfCorrectionJobs ($$) {
    my $wrk     = shift @_;
    my $asm     = shift @_;
    my $nJobs   = 0;
    my $nPerJob = 0;

    if (getGlobal("corConsensus") eq "falcon") {
        open(F, "ls $wrk/correction_inputs/ |") or caExit("can't find list of correction_inputs: $!", undef);
        while (<F>) {
            $nJobs++  if (m/^\d\d\d\d$/);
        }
        close(F);
    }

    if ((getGlobal("corConsensus") eq "utgcns") ||
        (getGlobal("corConsensus") eq "falconpipe")) {
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
    my $wrk  = shift @_;
    my $asm  = shift @_;
    my $bin  = getBinDirectory();
    my $cmd;

    my $path = "$wrk/2-correction";

    return  if (-e "$wrk/$asm.correctedReads.fastq");  #  Output exists
    return  if (-e "$path/cnsjob.files");              #  Jobs all finished

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
        $cmd .= "> $wrk/$asm.corStore.err 2>&1";

        if (runCommand($wrk, $cmd)) {
            caExit("failed to generate layouts for correction", "$wrk/$asm.corStore.err");
        }

        rename "$wrk/$asm.corStore.WORKING", "$wrk/$asm.corStore";
    }

    make_path("$path/correction_inputs")  if (! -d "$path/correction_inputs");
    make_path("$path/correction_outputs")  if (! -d "$path/correction_outputs");

    if (getGlobal("corConsensus") eq "falcon") {
        $cmd  = "$bin/createFalconSenseInputs \\\n";
        $cmd .= "  -G $wrk/$asm.gkpStore \\\n";
        $cmd .= "  -T $wrk/$asm.corStore 1 \\\n";
        $cmd .= "  -o $path/correction_inputs/ \\\n";
        $cmd .= "  -p " . getGlobal("corPartitions") . " \\\n";   #  NEED corPartitionMin!
        $cmd .= "> $path/correction_inputs.err 2>&1";

        if (runCommand($wrk, $cmd)) {
            caExit("failed to generate falcon inputs", "$path/correction_inputs.err");
        }
    }

    my ($jobs, $nPer) = computeNumberOfCorrectionJobs($wrk, $asm);

    my $taskID            = getGlobal("gridEngineTaskID");
    my $submitTaskID      = getGlobal("gridEngineArraySubmitID");

    open(F, "> $path/correctReads.sh") or caExit("can't open '$path/correctReads.sh'", undef);

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
    print F "if [ \$jobid -gt $jobs ]; then\n";
    print F "  echo Error: Only $jobs partitions, you asked for \$jobid.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "\n";
    print F "bgn=`expr \\( \$jobid - 1 \\) \\* $nPer`\n";
    print F "end=`expr \\( \$jobid + 0 \\) \\* $nPer`\n";
    print F "\n";
    print F "jobid=`printf %04d \$jobid`\n";
    print F "\n";
    print F "if [ -e \"$path/correction_outputs/\$jobid.fastq\" ] ; then\n";
    print F "  echo Job finished successfully.\n";
    print F "  exit 0\n";
    print F "fi\n";
    print F "\n";

    print F getBinDirectoryShellCode();

    if (getGlobal("corConsensus") eq "utgcns") {
        print F "\n";
        print F "$bin/utgcns \\\n";
        print F "  -u \$bgn-\$end \\\n";
        print F "  -e 0.30 \\\n";
        print F "  -G $wrk/$asm.gkpStore \\\n";
        print F "  -T $wrk/$asm.corStore 1 . \\\n";
        print F "  -O $path/correction_outputs/\$jobid.cns.WORKING \\\n";
        print F "  -L $path/correction_outputs/\$jobid.layout.WORKING \\\n";
        print F "  -F $path/correction_outputs/\$jobid.fastq.WORKING \\\n";
        print F "&& \\\n";
        print F "mv $path/correction_outputs/\$jobid.cns.WORKING $path/correction_outputs/\$jobid.cns \\\n";
        print F "&& \\\n";
        print F "mv $path/correction_outputs/\$jobid.layout.WORKING $path/correction_outputs/\$jobid.layout \\\n";
        print F "&& \\\n";
        print F "mv $path/correction_outputs/\$jobid.fastq.WORKING $path/correction_outputs/\$jobid.fastq \\\n";
        print F "\n";
    }

    if (getGlobal("corConsensus") eq "falcon") {
        print F "\n";
        print F getGlobal("falconSense") . " \\\n";
        print F "  --max_n_read 200 \\\n";
        print F "  --min_idt 0.70 \\\n";
        print F "  --output_multi \\\n";
        print F "  --local_match_count_threshold 2 \\\n";
        print F "  --min_cov 4 \\\n";
        print F "  --n_core " . getGlobal("corThreads") . " \\\n";
        print F "  < $path/correction_inputs/\$jobid \\\n";
        print F "  > $path/correction_outputs/\$jobid.fasta.WORKING \\\n";
        print F " 2> $path/correction_outputs/\$jobid.err \\\n";
        print F "&& \\\n";
        print F "mv $path/correction_outputs/\$jobid.fasta.WORKING $path/correction_outputs/\$jobid.fasta \\\n";
    }

    close(F);
}





#  For falcon_sense, using a pipe and no intermediate files
#
sub buildCorrectionLayouts_piped ($$) {
    my $wrk  = shift @_;
    my $asm  = shift @_;
    my $bin  = getBinDirectory();
    my $cmd;
    my $path = "$wrk/2-correction";

    return  if (-e "$wrk/$asm.correctedReads.fastq");  #  Output exists
    return  if (-e "$path/cnsjob.files");              #  Jobs all finished

    make_path("$path/correction_inputs")  if (! -d "$path/correction_inputs");
    make_path("$path/correction_outputs")  if (! -d "$path/correction_outputs");

    my ($nJobs, $nPerJob)  = computeNumberOfCorrectionJobs($wrk, $asm);  #  Does math based on number of reads and parameters.

    my $nReads             = getNumberOfReadsInStore($wrk, $asm);

    my $taskID             = getGlobal("gridEngineTaskID");
    my $submitTaskID       = getGlobal("gridEngineArraySubmitID");

    open(F, "> $path/correctReads.sh") or caExit("can't open '$path/correctReads.sh'", undef);

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
    print F "if [ \$jobid -gt $nJobs ]; then\n";
    print F "  echo Error: Only $nJobs partitions, you asked for \$jobid.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "\n";

    my  $bgnID   = 1;
    my  $endID   = $bgnID + $nPerJob;
    my  $jobID   = 1;

    while ($bgnID < $nReads) {
        $endID  = $bgnID + $nPerJob;
        $endID  = $nReads  if ($endID > $nReads);

        print F "if [ \$jobid -eq $jobID ] ; then\n";
        print F "  bgn=$bgnID\n";
        print F "  end=$endID\n";
        print F "fi\n";

        $bgnID = $endID;
        $jobID++;
    }

    print F "\n";
    print F "jobid=`printf %04d \$jobid`\n";
    print F "\n";
    print F "\n";
    print F getBinDirectoryShellCode();
    print F "\n";

    my $maxCov   = getCorCov($wrk, $asm, "Local");

    print F "$bin/generateCorrectionLayouts -b \$bgn -e \$end \\\n";
    print F "  -rl $path/$asm.readsToCorrect \\\n"                 if (-e "$path/$asm.readsToCorrect");
    print F "  -G $wrk/$asm.gkpStore \\\n";
    print F "  -O $wrk/$asm.ovlStore \\\n";
    print F "  -S $path/$asm.globalScores \\\n"                    if (-e "$path/$asm.globalScores");
    print F "  -L " . getGlobal("corMinEvidenceLength") . " \\\n"  if (defined(getGlobal("corMinEvidenceLength")));
    print F "  -E " . getGlobal("corMaxEvidenceErate")  . " \\\n"  if (defined(getGlobal("corMaxEvidenceErate")));
    print F "  -C $maxCov \\\n"                                    if (defined($maxCov));
    print F "  -F \\\n";
    print F "| \\\n";
    print F getGlobal("falconSense") . " \\\n";
    print F "  --max_n_read 200 \\\n";
    print F "  --min_idt 0.70 \\\n";
    print F "  --output_multi \\\n";
    print F "  --local_match_count_threshold 2 \\\n";
    print F "  --min_cov 4 \\\n";
    print F "  --n_core " . getGlobal("corThreads") . " \\\n";
    print F "  > $path/correction_outputs/\$jobid.fasta.WORKING \\\n";
    print F " 2> $path/correction_outputs/\$jobid.err \\\n";
    print F "&& \\\n";
    print F "mv $path/correction_outputs/\$jobid.fasta.WORKING $path/correction_outputs/\$jobid.fasta \\\n";
    print F "\n";
    print F "exit 0\n";

    close(F);
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


sub quickFilter ($$) {
    my $wrk  = shift @_;
    my $asm  = shift @_;
    my $bin  = getBinDirectory();
    my $path = "$wrk/2-correction";

    my $totCorLengthIn  = 0;
    my $totCorLengthOut = 0;
    my $minCorLength    = 0;

    my $minTotal  = getGlobal("genomeSize") * getGlobal("corOutCoverage");

    open(O, "> $path/$asm.readsToCorrect.WORKING") or caExit("can't open '$path/$asm.readsToCorrect.WORKING' for writing: $!\n", undef);
    open(F, "$bin/gatekeeperDumpMetaData -G $wrk/$asm.gkpStore -reads | sort -k3nr | ") or caExit("can't dump gatekeeper for read lengths: $!\n", undef);

    print O "read\toriginalLength\tcorrectedLength\n";

    while (<F>) {
        my @v = split '\s+', $_;

        $totCorLengthIn  += $v[2];
        $totCorLengthOut += $v[2];

        print O "$v[0]\t$v[2]\t0\n";

        if ($totCorLengthIn >= $minTotal) {
            $minCorLength = $v[2];
            last;
        }
    }
    close(F);
    close(O);

    rename "$path/$asm.readsToCorrect.WORKING", "$path/$asm.readsToCorrect";
}




sub expensiveFilter ($$) {
    my $wrk  = shift @_;
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
        $cmd .= "  -p $path/$asm.estimate";

        if (runCommand($wrk, $cmd)) {
            rename "$path/$asm.estimate.log", "$path/$asm.estimate.log.FAILED";
            caExit("failed to generate estimated lengths of corrected reads", "$wrk/$asm.corStore.err");
        }
    }

    runCommandSilently($path, "sort -k4nr -k2nr < $path/$asm.estimate.log > $path/$asm.estimate.correctedLength.log");
    runCommandSilently($path, "sort -k2nr -k4nr < $path/$asm.estimate.log > $path/$asm.estimate.originalLength.log");

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

    open(O, "> $path/$asm.readsToCorrect.WORKING") or caExit("can't open '$path/$asm.readsToCorrect.WORKING' for writing: $!\n", undef);
    open(F, "< $path/$asm.estimate.correctedLength.log");

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
    close(F);
    close(O);

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

    if (! -e "$path/$asm.estimate.original-x-correctedLength.png") {
        open(F, "> $path/$asm.estimate.original-x-correctedLength.gp");
        print F "set title 'original length (x) vs corrected length (y)'\n";
        print F "set xlabel 'original read length'\n";
        print F "set ylabel 'corrected read length (expected)'\n";
        print F "set pointsize 0.25\n";
        print F "set terminal png size 1024,1024\n";
        print F "set output '$path/$asm.estimate.original-x-corrected.png'\n";
        print F "plot '$path/$asm.estimate.tn.log' using 2:4 title 'tn', \\\n";
        print F "     '$path/$asm.estimate.fn.log' using 2:4 title 'fn', \\\n";
        print F "     '$path/$asm.estimate.fp.log' using 2:4 title 'fp', \\\n";
        print F "     '$path/$asm.estimate.tp.log' using 2:4 title 'tp'\n";
        close(F);

        runCommandSilently($path, "gnuplot < $path/$asm.estimate.original-x-correctedLength.gp > /dev/null 2>&1");
    }

    #  These not so interesting

    #if (! -e "$path/$asm.estimate.correctedLength.png") {
    #    open(F, "> $path/$asm.estimate.correctedLength.gp");
    #    print F "set title 'sorted corrected length vs original length'\n";
    #    print F "set pointsize 0.25\n";
    #    print F "set terminal png size 1024,1024\n";
    #    print F "set output '$path/$asm.estimate.correctedLength.png'\n";
    #    print F "plot '$path/$asm.estimate.correctedLength.log' using 2 title 'original length', '$path/$asm.estimate.correctedLength.log' using 4 with lines title 'corrected length'\n";
    #    close(F);
    #
    #    runCommandSilently($path, "gnuplot < $path/$asm.estimate.correctedLength.gp > /dev/null 2>&1");
    #}

    #if (! -e "$path/$asm.estimate.originalLength.png") {
    #    open(F, "> $path/$asm.estimate.originalLength.gp");
    #    print F "set title 'sorted original length vs corrected length'\n";
    #    print F "set pointsize 0.25\n";
    #    print F "set terminal png size 1024,1024\n";
    #    print F "set output '$path/$asm.estimate.originalLength.png'\n";
    #    print F "plot '$path/$asm.estimate.originalLength.log' using 4 title 'corrected length', '$path/$asm.estimate.originalLength.log' using 2 with lines title 'original length'\n";
    #    close(F);
    #
    #    runCommandSilently($path, "gnuplot < $path/$asm.estimate.originalLength.gp > /dev/null 2>&1");
    #}
}





sub buildCorrectionLayouts ($$) {
    my $WRK  = shift @_;
    my $wrk  = "$WRK/correction";
    my $asm  = shift @_;
    my $bin  = getBinDirectory();
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

    return  if (skipStage($WRK, $asm, "cor-buildCorrectionLayouts") == 1);
    return  if (-e "$wrk/$asm.correctedReads.fastq");  #  Output exists
    return  if (-e "$path/cnsjob.files");              #  Jobs all finished
    return  if (-e "$path/correctReads.sh");           #  Jobs created

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
        $cmd .= "  -logfile $path/$asm.globalScores.log \\\n";
        $cmd .= "> $path/$asm.globalScores.err 2>&1";

        if (runCommand($path, $cmd)) {
            caExit("failed to globally filter overlaps for correction", "$path/$asm.globalScores.err");
        }
    }

    #  For 'quick' filtering, but more reads to correct, sort the reads by length, and correct the
    #  longest Nx of reads.
    #
    #  For 'expensive' filtering, but fewer reads to correct, first estimate the corrected lengths
    #  (requires a pass through the overlaps), then pull out the longest Nx of corrected reads.
    #
    #  Both are required to create a file $asm.readsToCorrect, containing a list of IDs to correct.

    if      (getGlobal("corFilter") eq "quick") {
        quickFilter($wrk, $asm);

    } elsif (getGlobal("corFilter") eq "expensive") {
        expensiveFilter($wrk, $asm);

    } else {
        caFailure("unknown corFilter '" . getGlobal("corFilter") . "'", undef);
    }

    caExit("failed to create list of reads to correct", undef)  if (! -e "$path/$asm.readsToCorrect");

    buildCorrectionLayouts_direct($wrk, $asm)      if (getGlobal("corConsensus") eq "utgcns");
    buildCorrectionLayouts_direct($wrk, $asm)      if (getGlobal("corConsensus") eq "falcon");
    buildCorrectionLayouts_piped($wrk, $asm)       if (getGlobal("corConsensus") eq "falconpipe");

    emitStage($WRK, $asm, "cor-buildCorrectionLayouts");
}




sub generateCorrectedReads ($$$) {
    my $WRK     = shift @_;
    my $wrk     = "$WRK/correction";
    my $asm     = shift @_;
    my $attempt = shift @_;
    my $bin     = getBinDirectory();

    my $path    = "$wrk/2-correction";

    return  if (skipStage($WRK, $asm, "cor-generateCorrectedReads", $attempt) == 1);
    return  if (-e "$wrk/$asm.correctedReads.fastq");

    my ($jobs, undef) = computeNumberOfCorrectionJobs($wrk, $asm);

    my $currentJobID = "0001";
    my @successJobs;
    my @failedJobs;
    my $failureMessage = "";

    for (my $job=1; $job <= $jobs; $job++) {
        if (-e "$path/correction_outputs/$currentJobID.fasta") {
            push @successJobs, "$path/correction_outputs/$currentJobID.fasta\n";

        } else {
            $failureMessage .= "  job $path/correction_outputs/$currentJobID.fasta FAILED.\n";
            push @failedJobs, $job;
        }

        $currentJobID++;
    }

    #  No failed jobs?  Success!

    if (scalar(@failedJobs) == 0) {
        open(L, "> $path/corjob.files") or caExit("failed to open '$path/corjob.files'", undef);
        print L @successJobs;
        close(L);
        setGlobal("ca3gIteration", 0);
        emitStage($WRK, $asm, "cor-generateCorrectedReads");
        return;
    }

    #  If not the first attempt, report the jobs that failed, and that we're recomputing.

    if ($attempt > 1) {
        print STDERR "\n";
        print STDERR scalar(@failedJobs), " read correction jobs failed:\n";
        print STDERR $failureMessage;
        print STDERR "\n";
    }

    #  If too many attempts, give up.

    if ($attempt > 2) {
        caExit("failed to generate corrected reads.  Made " . ($attempt-1) . " attempts, jobs still failed", undef);
    }

    #  Otherwise, run some jobs.

    print STDERR "generateCorrectedReads() -- attempt $attempt begins with ", scalar(@successJobs), " finished, and ", scalar(@failedJobs), " to compute.\n";

    emitStage($WRK, $asm, "cor-generateCorrectedReads", $attempt);

    submitOrRunParallelJob($wrk, $asm, "cor", $path, "correctReads", @failedJobs);
}



sub dumpCorrectedReads ($$) {
    my $WRK  = shift @_;
    my $wrk  = "$WRK/correction";
    my $asm  = shift @_;
    my $bin  = getBinDirectory();

    my $path = "$wrk/2-correction";

    return  if (skipStage($WRK, $asm, "cor-dumpCorrectedReads") == 1);
    return  if (-e "$wrk/$asm.correctedReads.fastq");

    my $files = 0;
    my $reads = 0;

    open(F, "< $path/corjob.files")             or caExit("can't open '$path/corjob.files' for reading: $!", undef);
    open(O, "> $wrk/$asm.correctedReads.fastq") or caExit("can't open '$wrk/$asm.correctedReads.fastq' for writing: $!", undef);

    while (<F>) {
        chomp;

        open(R, "< $_") or caExit("can't open correction output '$_' for reading: $!\n", undef);

        while (!eof(R)) {
            my $n = <R>;
            my $s = <R>;
            my $q = $s;

            $n =~ s/^>/\@/;

            $q =~ tr/[A-Z][a-z]/*/;

            print O $n;
            print O $s;
            print O "+\n";
            print O $q;

            $reads++;
        }

        $files++;
        close(R);
    }

    close(O);
    close(F);

    print STDERR "dumpCorrectedReads()-- wrote $reads corrected reads from $files files into '$wrk/$asm.correctedReads.fastq'.\n";

    emitStage($WRK, $asm, "cor-dumpCorrectedReads");
}
