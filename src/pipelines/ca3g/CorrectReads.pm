package ca3g::CorrectReads;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(buildCorrectionLayouts generateCorrectedReads dumpCorrectedReads);

use strict;

use File::Path qw(make_path remove_tree);

use ca3g::Defaults;
use ca3g::Execution;
use ca3g::Gatekeeper;


sub computeNumberOfCorrectionJobs ($$) {
    my $wrk    = shift @_;
    my $asm    = shift @_;
    my $nJobs  = 0;

    if (getGlobal("consensus") eq "falcon") {
        open(F, "ls $wrk/2-correction/falcon_inputs/ |") or die;
        while (<F>) {
            $nJobs++  if (m/^\d\d\d\d$/);
        }
        close(F);
    }

    #my $jobs = 0;
    #open(F, "< $wrk/4-unitigger/$asm.partitioningInfo") or caFailure("can't open '$wrk/4-unitigger/$asm.partitioningInfo'", undef);
    #while (<F>) {
    #    if (m/Partition\s+(\d+)\s+has\s+(\d+)\s+unitigs\sand\s+(\d+)\s+fragments./) {
    #        $jobs = $1;
    #    }
    #}
    #close(F);

    return($nJobs);
}


sub buildCorrectionLayouts ($$) {
    my $wrk          = shift @_;
    my $asm          = shift @_;
    my $bin          = getBinDirectory();
    my $cmd;
    my $path         = "$wrk/2-correction";

    return  if (-e "$wrk/2-correction/consensus.sh");
    return  if (-e "$wrk/correctedReads.fastq");
    return  if (-e "$path/cnsjob.files");

    make_path("$path")  if (! -d "$path");

    #  Assumes we have partial overlaps.  If mhap, quality might improve if mhapReAlign=raw is used.

    if (! -e "$wrk/$asm.corStore") {
        $cmd  = "$bin/generateCorrectionLayouts \\\n";
        $cmd .= "  -G $wrk/$asm.gkpStore \\\n";
        $cmd .= "  -O $wrk/$asm.ovlStore \\\n";
        $cmd .= "  -T $wrk/$asm.corStore.WORKING \\\n";
        $cmd .= "  -L 1000 \\\n";
        $cmd .= "  -E 0.18 \\\n";
        $cmd .= "> $wrk/$asm.corStore.err 2>&1";

        if (runCommand($wrk, $cmd)) {
            caFailure("failed to generate layouts for correction", "$wrk/$asm.corStore.err");
        }

        rename "$wrk/$asm.corStore.WORKING", "$wrk/$asm.corStore";
    }


    #  If we're doing falcon_sense, dump the inputs first.


    if ((getGlobal("consensus") eq "falcon") &&
        (! -e "$path/falcon_inputs/success")) {
        make_path("$path/falcon_inputs")  if (! -d "$path/falcon_inputs");
        make_path("$path/falcon_outputs")  if (! -d "$path/falcon_outputs");

        $cmd  = "$bin/createFalconSenseInputs \\\n";
        $cmd .= "  -G $wrk/$asm.gkpStore \\\n";
        $cmd .= "  -T $wrk/$asm.corStore 1 \\\n";
        $cmd .= "  -o $wrk/2-correction/falcon_inputs/ \\\n";
        $cmd .= "  -p " . getGlobal("cnsPartitions") . " \\\n";
        $cmd .= "> $wrk/2-correction/falcon_inputs.err 2>&1";

        if (runCommand($wrk, $cmd)) {
            caFailure("failed to generate falcon inputs", "$wrk/2-correction/falcon_inputs.err");
        }

        touch("$path/falcon_inputs/success");
    }

    my $jobs = computeNumberOfCorrectionJobs($wrk, $asm);


    #  And a script to run consensus on those layouts.  If we're doing falcon_sense, dump the layouts first.

    if (! -e "$wrk/2-correction/consensus.sh") {
        my $taskID            = getGlobal("gridEngineTaskID");
        my $submitTaskID      = getGlobal("gridEngineArraySubmitID");

        open(F, "> $wrk/2-correction/consensus.sh") or caFailure("can't open '$wrk/2-correction/consensus.sh'", undef);

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
        #print F "if [ -e $wrk/5-consensus/${asm}_\$jobid.success ] ; then\n";
        #print F "  exit 0\n";
        #print F "fi\n";
        print F "\n";

        print F getBinDirectoryShellCode();

        if (getGlobal("consensus") eq "utgcns") {
            print F "\n";
            print F "\$bin/utgcns \\\n";
            print F "  -G $wrk/$asm.gkpStore \\\n";
            print F "  -T $wrk/$asm.tigStore 1 \$jobid \\\n";
            print F "  -O $wrk/5-consensus/\$jobid.cns.WORKING \\\n";
            print F "  -L $wrk/5-consensus/\$jobid.layouts.WORKING \\\n";
            print F "  -maxcoverage " . getGlobal('cnsMaxCoverage') . " \\\n";
            print F "> $wrk/5-consensus/\$jobid.cns.err 2>&1 \\\n";
            print F "&& \\\n";
            print F "mv $wrk/5-consensus/\$jobid.cns.WORKING $wrk/5-consensus/\$jobid.cns \\\n";
            print F "&& \\\n";
            print F "mv $wrk/5-consensus/\$jobid.layouts.WORKING $wrk/5-consensus/\$jobid.layouts \\\n";
        }

        if (getGlobal("consensus") eq "falcon") {
            print F "\n";
            print F getGlobal("falcon") . " \\\n";
            print F "  --max_n_read 200 \\\n";
            #print F "  --trim \\\n";
            print F "  --min_idt 0.70 \\\n";
            print F "  --output_multi \\\n";
            print F "  --local_match_count_threshold 2 \\\n";
            print F "  --min_cov 4 \\\n";
            print F "  --n_core 1 \\\n";
            print F "  < $wrk/2-correction/falcon_inputs/\$jobid \\\n";
            print F "  > $wrk/2-correction/falcon_outputs/\$jobid.fasta.WORKING \\\n";
            print F " 2> $wrk/2-correction/falcon_outputs/\$jobid.err \\\n";
            print F "&& \\\n";
            print F "mv $wrk/2-correction/falcon_outputs/\$jobid.fasta.WORKING $wrk/2-correction/falcon_outputs/\$jobid.fasta \\\n";
        }


        close(F);
    }
}




sub generateCorrectedReads ($$$) {
    my $wrk          = shift @_;
    my $asm          = shift @_;
    my $attempt      = shift @_;
    my $bin          = getBinDirectory();

    my $path         = "$wrk/2-correction";
    my $script       = "consensus";
    my $jobType      = "cns";

    return  if (-e "$wrk/correctedReads.fastq");
    return  if (-e "$path/cnsjob.files");

    #  How many partitions?  There should be a classier way...

    my $jobs = computeNumberOfCorrectionJobs($wrk, $asm);

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
        open(L, "> $path/cnsjob.files") or caFailure("failed to open '$path/cnsjob.files'", undef);
        print L @successJobs;
        close(L);

        touch("$wrk/2-correction/success");

        return;
    }

    if ($attempt > 0) {
        print STDERR "\n";
        print STDERR scalar(@failedJobs), " consensus jobs failed:\n";
        print STDERR $failureMessage;
        print STDERR "\n";
    }

    print STDERR "generateCorrectedReads() -- attempt $attempt begins with ", scalar(@successJobs), " finished, and ", scalar(@failedJobs), " to compute.\n";

    if ($attempt < 1) {
        submitOrRunParallelJob($wrk, $asm, $jobType, $path, $script, getGlobal("cnsConcurrency"), @failedJobs);
    } else {
        caFailure("failed to generate corrected reads.  Made $attempt attempts, jobs still failed.", undef);
    }

    stopAfter("consensus");

}


sub dumpCorrectedReads ($$) {
    my $wrk          = shift @_;
    my $asm          = shift @_;
    my $bin          = getBinDirectory();

    my $path         = "$wrk/2-correction";

    my $files = 0;
    my $reads = 0;

    return  if (-e "$wrk/correctedReads.fastq");

    open(F, "< $path/cnsjob.files") or die;
    open(O, "> $wrk/correctedReads.fastq") or die;

    while (<F>) {
        chomp;

        open(R, "< $_") or die "Failed to open consensus output '$_': $_\n";

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

    print STDERR "dumpCorrectedReads()-- wrote $reads corrected reads from $files files into '$wrk/correctedReads.fastq'.\n";
}
