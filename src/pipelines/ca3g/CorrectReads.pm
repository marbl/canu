package ca3g::CorrectReads;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(buildCorrectionLayouts generateCorrectedReads dumpCorrectedReads);

use strict;

use File::Path qw(make_path remove_tree);

use ca3g::Defaults;
use ca3g::Execution;
use ca3g::Gatekeeper;


#  Return the number of jobs for 'falcon', 'falconpipe' or 'utgcns'
#
sub computeNumberOfCorrectionJobs ($$) {
    my $wrk     = shift @_;
    my $asm     = shift @_;
    my $nJobs   = 0;
    my $nPerJob = 0;

    if (getGlobal("correction") eq "falcon") {
        open(F, "ls $wrk/2-correction/falcon_inputs/ |") or die;
        while (<F>) {
            $nJobs++  if (m/^\d\d\d\d$/);
        }
        close(F);
    }

    if (getGlobal("correction") eq "falconpipe") {
        my $nPart    = getGlobal("corPartitions");
        my $nReads   = getNumberOfReadsInStore($wrk, $asm);

        $nPerJob     = int($nReads / $nPart + 1);
        $nPerJob     = getGlobal("corPartitionMin")  if ($nPerJob < getGlobal("corPartitionMin"));

        for (my $j=1; $j<=$nReads; $j += $nPerJob) {  #  We could just divide, except for rounding issues....
            $nJobs++;
        }
    }

    if (getGlobal("correction") eq "utgcns") {
    }

    return($nJobs, $nPerJob);
}


#  Generate a corStore, dump files for falcon to process, generate a script to run falcon.
#
sub buildCorrectionLayouts_falcon ($$) {
    my $wrk          = shift @_;
    my $asm          = shift @_;
    my $bin          = getBinDirectory();
    my $cmd;
    my $path         = "$wrk/2-correction";

    return  if (-e "$wrk/$asm.correctedReads.fastq");  #  Output exists
    return  if (-e "$path/cnsjob.files");              #  Jobs all finished

    my $maxCov = (defined(getGlobal("corMaxCoverage"))) ? getGlobal("corMaxCoverage") : int(getExpectedCoverage($wrk, $asm) * 1.5);

    if (! -e "$wrk/$asm.corStore") {
        $cmd  = "$bin/generateCorrectionLayouts \\\n";
        $cmd .= "  -G $wrk/$asm.gkpStore \\\n";
        $cmd .= "  -O $wrk/$asm.ovlStore \\\n";
        $cmd .= "  -T $wrk/$asm.corStore.WORKING \\\n";
        $cmd .= "  -L " . getGlobal("corMinLength") . " \\\n"  if (defined(getGlobal("corMinLength")));
        $cmd .= "  -E " . getGlobal("corMaxErate")  . " \\\n"   if (defined(getGlobal("corMaxErate")));
        $cmd .= "  -C $maxCov \\\n";
        $cmd .= "> $wrk/$asm.corStore.err 2>&1";

        if (runCommand($wrk, $cmd)) {
            caFailure("failed to generate layouts for correction", "$wrk/$asm.corStore.err");
        }

        rename "$wrk/$asm.corStore.WORKING", "$wrk/$asm.corStore";
    }

    make_path("$path/falcon_inputs")  if (! -d "$path/falcon_inputs");
    make_path("$path/falcon_outputs")  if (! -d "$path/falcon_outputs");

    $cmd  = "$bin/createFalconSenseInputs \\\n";
    $cmd .= "  -G $wrk/$asm.gkpStore \\\n";
    $cmd .= "  -T $wrk/$asm.corStore 1 \\\n";
    $cmd .= "  -o $path/falcon_inputs/ \\\n";
    $cmd .= "  -p " . getGlobal("corPartitions") . " \\\n";   #  NEED corPartitionMin!
    $cmd .= "> $path/falcon_inputs.err 2>&1";

    if (runCommand($wrk, $cmd)) {
        caFailure("failed to generate falcon inputs", "$path/falcon_inputs.err");
    }

    my ($jobs, undef) = computeNumberOfCorrectionJobs($wrk, $asm);  #  Counts the number of files in falcon_inputs

    my $taskID            = getGlobal("gridEngineTaskID");
    my $submitTaskID      = getGlobal("gridEngineArraySubmitID");

    open(F, "> $path/correctReads.sh") or caFailure("can't open '$path/correctReads.sh'", undef);

    print F "#!" . getGlobal("shell") . "\n";
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
    print F "if [ \$jobid -gt $jobs ]; then\n";
    print F "  echo Error: Only $jobs partitions, you asked for \$jobid.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "\n";
    print F "jobid=`printf %04d \$jobid`\n";
    print F "\n";

    print F getBinDirectoryShellCode();

    if (getGlobal("correction") eq "utgcns") {
        caFailure("correction=utgcns not supported for correction", undef);
    }

    if (getGlobal("correction") eq "falcon") {
        print F "\n";
        print F getGlobal("falconSense") . " \\\n";
        print F "  --max_n_read 200 \\\n";
        print F "  --min_idt 0.70 \\\n";
        print F "  --output_multi \\\n";
        print F "  --local_match_count_threshold 2 \\\n";
        print F "  --min_cov 4 \\\n";
        print F "  --n_core " . getGlobal("corThreads") . " \\\n";
        print F "  < $path/falcon_inputs/\$jobid \\\n";
        print F "  > $path/falcon_outputs/\$jobid.fasta.WORKING \\\n";
        print F " 2> $path/falcon_outputs/\$jobid.err \\\n";
        print F "&& \\\n";
        print F "mv $path/falcon_outputs/\$jobid.fasta.WORKING $path/falcon_outputs/\$jobid.fasta \\\n";
    }

    close(F);
}



#  For falcon_sense, using a pipe and no intermediate files
#
sub buildCorrectionLayouts_falconpipe ($$) {
    my $wrk          = shift @_;
    my $asm          = shift @_;
    my $bin          = getBinDirectory();
    my $cmd;
    my $path         = "$wrk/2-correction";

    return  if (-e "$wrk/$asm.correctedReads.fastq");  #  Output exists
    return  if (-e "$path/cnsjob.files");              #  Jobs all finished

    make_path("$path/falcon_inputs")  if (! -d "$path/falcon_inputs");
    make_path("$path/falcon_outputs")  if (! -d "$path/falcon_outputs");

    my ($nJobs, $nPerJob)  = computeNumberOfCorrectionJobs($wrk, $asm);  #  Does math based on number of reads and parameters.

    my $nReads             = getNumberOfReadsInStore($wrk, $asm);

    my $taskID             = getGlobal("gridEngineTaskID");
    my $submitTaskID       = getGlobal("gridEngineArraySubmitID");

    open(F, "> $path/correctReads.sh") or caFailure("can't open '$path/correctReads.sh'", undef);

    print F "#!" . getGlobal("shell") . "\n";
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

    my $maxCov   = (defined(getGlobal("corMaxCoverage"))) ? getGlobal("corMaxCoverage") : int(getExpectedCoverage($wrk, $asm) * 1.5);

    print F "$bin/generateCorrectionLayouts -b \$bgn -e \$end \\\n";
    print F "  -G $wrk/$asm.gkpStore \\\n";
    print F "  -O $wrk/$asm.ovlStore \\\n";
    print F "  -L " . getGlobal("corMinLength") . " \\\n"  if (defined(getGlobal("corMinLength")));
    print F "  -E " . getGlobal("corMaxErate")  . " \\\n"   if (defined(getGlobal("corMaxErate")));
    print F "  -C $maxCov \\\n";
    print F "  -F \\\n";
    print F "| \\\n";
    print F getGlobal("falconSense") . " \\\n";
    print F "  --max_n_read 200 \\\n";
    print F "  --min_idt 0.70 \\\n";
    print F "  --output_multi \\\n";
    print F "  --local_match_count_threshold 2 \\\n";
    print F "  --min_cov 4 \\\n";
    print F "  --n_core " . getGlobal("corThreads") . " \\\n";
    print F "  > $path/falcon_outputs/\$jobid.fasta.WORKING \\\n";
    print F " 2> $path/falcon_outputs/\$jobid.err \\\n";
    print F "&& \\\n";
    print F "mv $path/falcon_outputs/\$jobid.fasta.WORKING $path/falcon_outputs/\$jobid.fasta \\\n";
    print F "\n";
    print F "exit 0\n";

    close(F);
}



sub buildCorrectionLayouts_utgcns ($$) {
    my $wrk          = shift @_;
    my $asm          = shift @_;
    my $bin          = getBinDirectory();
    my $cmd;
    my $path         = "$wrk/2-correction";

    return  if (-e "$wrk/$asm.correctedReads.fastq");  #  Output exists
    return  if (-e "$path/cnsjob.files");              #  Jobs all finished

    #  Extend generateCorrectionLayouts to call consensus directly
}



sub buildCorrectionLayouts ($$) {
    my $wrk          = shift @_;
    my $asm          = shift @_;
    my $bin          = getBinDirectory();
    my $cmd;
    my $path         = "$wrk/2-correction";

    return  if (-e "$wrk/$asm.correctedReads.fastq");  #  Output exists
    return  if (-e "$path/cnsjob.files");              #  Jobs all finished
    return  if (-e "$path/correctReads.sh");           #  Jobs created

    make_path("$path")  if (! -d "$path");

    buildCorrectionLayouts_falcon($wrk, $asm)      if (getGlobal('correction') eq "falcon");
    buildCorrectionLayouts_falconpipe($wrk, $asm)  if (getGlobal('correction') eq "falconpipe");
    buildCorrectionLayouts_utgcns($wrk, $asm)      if (getGlobal('correction') eq "utgcns");
}




sub generateCorrectedReads ($$$) {
    my $wrk          = shift @_;
    my $asm          = shift @_;
    my $attempt      = shift @_;
    my $bin          = getBinDirectory();

    my $path         = "$wrk/2-correction";

    return  if (-e "$wrk/$asm.correctedReads.fastq");

    my ($jobs, undef) = computeNumberOfCorrectionJobs($wrk, $asm);

    my $currentJobID = "0001";
    my @successJobs;
    my @failedJobs;
    my $failureMessage = "";

    for (my $job=1; $job <= $jobs; $job++) {
        if (-e "$path/falcon_outputs/$currentJobID.fasta") {
            push @successJobs, "$path/falcon_outputs/$currentJobID.fasta\n";

        } else {
            $failureMessage .= "  job $path/falcon_outputs/$currentJobID.fasta FAILED.\n";
            push @failedJobs, $job;
        }

        $currentJobID++;
    }

    if (scalar(@failedJobs) == 0) {
        open(L, "> $path/corjob.files") or caFailure("failed to open '$path/corjob.files'", undef);
        print L @successJobs;
        close(L);

        touch("$path/success");

        return;
    }

    if ($attempt > 0) {
        print STDERR "\n";
        print STDERR scalar(@failedJobs), " read correction jobs failed:\n";
        print STDERR $failureMessage;
        print STDERR "\n";
    }

    print STDERR "generateCorrectedReads() -- attempt $attempt begins with ", scalar(@successJobs), " finished, and ", scalar(@failedJobs), " to compute.\n";

    if ($attempt < 1) {
        submitOrRunParallelJob($wrk, $asm, "cor", $path, "correctReads", @failedJobs);
    } else {
        caFailure("failed to generate corrected reads.  Made $attempt attempts, jobs still failed.", undef);
    }

    stopAfter("correctReads");

}


sub dumpCorrectedReads ($$) {
    my $wrk          = shift @_;
    my $asm          = shift @_;
    my $bin          = getBinDirectory();

    my $path         = "$wrk/2-correction";

    my $files = 0;
    my $reads = 0;

    return  if (-e "$wrk/$asm.correctedReads.fastq");

    open(F, "< $path/corjob.files") or die;
    open(O, "> $wrk/$asm.correctedReads.fastq") or die;

    while (<F>) {
        chomp;

        open(R, "< $_") or die "Failed to open read correction output '$_': $_\n";

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
}
