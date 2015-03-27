package ca3g::OverlapInCore;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(overlapConfigure overlap overlapCheck);

use strict;

use File::Path qw(make_path remove_tree);

use ca3g::Defaults;
use ca3g::Execution;


sub overlapConfigure ($$$) {
    my $wrk          = shift @_;
    my $asm          = shift @_;
    my $bin          = getBinDirectory();
    my $cmd;
    my $type         = shift @_;

    caFailure("invalid type '$type'", undef)  if (($type ne "partial") && ($type ne "normal"));

    return  if  (-d "$wrk/$asm.tigStore");

    return  if ((-d "$wrk/$asm.obtStore") && ($type eq "partial"));
    return  if ((-d "$wrk/$asm.ovlStore") && ($type eq "normal"));

    return  if (-e "$wrk/1-overlapper/ovlopt");

    make_path("$wrk/1-overlapper") if (! -d "$wrk/1-overlapper");
    
    if (! -e "$wrk/1-overlapper/ovlopt") {
        my $checkLibrary       = getGlobal("ovlCheckLibrary");
        my $hashLibrary        = getGlobal("ovlHashLibrary");
        my $refLibrary         = getGlobal("ovlRefLibrary");

        my $ovlHashBlockLength = getGlobal("ovlHashBlockLength");
        my $ovlHashBlockSize   = 0;
        my $ovlRefBlockSize    = getGlobal("ovlRefBlockSize");
        my $ovlRefBlockLength  = getGlobal("ovlRefBlockLength");

        if (($ovlRefBlockSize > 0) && ($ovlRefBlockLength > 0)) {
            caFailure("can't set both ovlRefBlockSize and ovlRefBlockLength", undef);
        }

        $cmd  = "$bin/overlapInCorePartition \\\n";
        $cmd .= " -g  $wrk/$asm.gkpStore \\\n";
        $cmd .= " -bl $ovlHashBlockLength \\\n";
        $cmd .= " -bs $ovlHashBlockSize \\\n";
        $cmd .= " -rs $ovlRefBlockSize \\\n";
        $cmd .= " -rl $ovlRefBlockLength \\\n";
        $cmd .= " -H $hashLibrary \\\n" if ($hashLibrary ne "0");
        $cmd .= " -R $refLibrary \\\n"  if ($refLibrary ne "0");
        $cmd .= " -C \\\n" if (!$checkLibrary);
        $cmd .= " -o  $wrk/1-overlapper \\\n";
        $cmd .= "> $wrk/1-overlapper/overlapInCorePartition.err 2>&1";

        if (runCommand($wrk, $cmd)) {
            caFailure("failed partition for overlapper", undef);
        }
    }

    if (! -e "$wrk/1-overlapper/overlap.sh") {
        my $merSize      = getGlobal("ovlMerSize");
        
        my $hashLibrary  = getGlobal("ovlHashLibrary");
        my $refLibrary   = getGlobal("ovlRefLibrary");
        
        #  Create a script to run overlaps.  We make a giant job array for this -- we need to know
        #  hashBeg, hashEnd, refBeg and refEnd -- from that we compute batchName and jobName.

        my $ovlThreads        = getGlobal("ovlThreads");
        my $ovlHashBits       = getGlobal("ovlHashBits");
        my $ovlHashLoad       = getGlobal("ovlHashLoad");

        my $taskID            = getGlobal("gridEngineTaskID");
        my $submitTaskID      = getGlobal("gridEngineArraySubmitID");

        open(F, "> $wrk/1-overlapper/overlap.sh") or caFailure("can't open '$wrk/1-overlapper/overlap.sh'", undef);
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
        print F "bat=`head -n \$jobid $wrk/1-overlapper/ovlbat | tail -n 1`\n";
        print F "job=`head -n \$jobid $wrk/1-overlapper/ovljob | tail -n 1`\n";
        print F "opt=`head -n \$jobid $wrk/1-overlapper/ovlopt | tail -n 1`\n";
        print F "jid=\$\$\n";
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
        print F "if [ x\$bat = x ]; then\n";
        print F "  echo Error: Job index out of range.\n";
        print F "  exit 1\n";
        print F "fi\n";
        print F "\n";

        print F getBinDirectoryShellCode();

        print F "\$bin/overlapInCore \\\n";
        print F "  -G \\\n"  if ($type eq "partial");
        print F "  -t $ovlThreads \\\n";
        print F "  -k $merSize \\\n";
        print F "  -k $wrk/0-mercounts/$asm.ms$merSize.frequentMers.fasta \\\n";
        print F "  --hashbits $ovlHashBits \\\n";
        print F "  --hashload $ovlHashLoad \\\n";
        print F "  --maxerate  ", getGlobal("ovlErrorRate"), " \\\n";
        print F "  --minlength ", getGlobal("ovlMinLen"), " \\\n";
        print F "  \$opt \\\n";
        print F "  -o $wrk/1-overlapper/\$bat/\$job.ovb.WORKING.gz \\\n";
        print F "  -H $hashLibrary \\\n" if ($hashLibrary ne "0");
        print F "  -R $refLibrary \\\n"  if ($refLibrary  ne "0");
        print F "  $wrk/$asm.gkpStore \\\n";
        print F "&& \\\n";
        print F "mv $wrk/1-overlapper/\$bat/\$job.ovb.WORKING.gz $wrk/1-overlapper/\$bat/\$job.ovb.gz\n";
        print F "\n";
        print F "exit 0\n";
        close(F);

        system("chmod +x $wrk/1-overlapper/overlap.sh");
    }

    caFailure("failed to find overlapInCorePartition output $wrk/1-overlapper/ovlbat", undef)  if (! -e "$wrk/1-overlapper/ovlbat");
    caFailure("failed to find overlapInCorePartition output $wrk/1-overlapper/ovljob", undef)  if (! -e "$wrk/1-overlapper/ovljob");
    caFailure("failed to find overlapInCorePartition output $wrk/1-overlapper/ovlopt", undef)  if (! -e "$wrk/1-overlapper/ovlopt");

    open(F, "< $wrk/1-overlapper/ovlbat") or caFailure("failed partition for overlapper: no ovlbat file found", undef);
    my @bat = <F>;
    close(F);

    open(F, "< $wrk/1-overlapper/ovljob") or caFailure("failed partition for overlapper: no ovljob file found", undef);
    my @job = <F>;
    close(F);

    my $jobs      = scalar(@job);
    my $batchName = $bat[$jobs-1];  chomp $batchName;
    my $jobName   = $job[$jobs-1];  chomp $jobName;

    print STDERR "Created $jobs overlap jobs.  Last batch '$batchName', last job '$jobName'.\n";

    stopAfter("overlap-configure");
    
    #submitOrRunParallelJob($wrk, $asm, "ovl", "$wrk/1-overlapper", "overlap", getGlobal("ovlConcurrency"), "1-$jobs");
}




#  Check that the overlapper jobs properly executed.  If not,
#  complain, but don't help the user fix things.
#
sub overlapCheck ($$$$) {
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

    my $currentJobID   = 1;
    my @successJobs;
    my @failedJobs;
    my $failureMessage = "";

    open(B, "< $path/ovlbat")       or caFailure("failed to open '$path/ovlbat'",       undef);
    open(J, "< $path/ovljob")       or caFailure("failed to open '$path/ovljob'",       undef);

    while (!eof(B) && !eof(J)) {
        my $b = <B>;  chomp $b;
        my $j = <J>;  chomp $j;

        if      (-e "$path/$b/$j.ovb.gz") {
            push @successJobs, "$path/$b/$j.ovb.gz\n";

        } elsif (-e "$path/$b/$j.ovb") {
            push @successJobs, "$path/$b/$j.ovb\n";

        } elsif (-e "$path/$b/$j.ovb.bz2") {
            push @successJobs, "$path/$b/$j.ovb.bz2\n";

        } elsif (-e "$path/$b/$j.ovb.xz") {
            push @successJobs, "$path/$b/$j.ovb.xz\n";

        } else {
            $failureMessage .= "  job $path/$b/$j FAILED.\n";
            push @failedJobs, $currentJobID;
        }

        $currentJobID++;
    }

    print STDERR "Partitioning error; '$path/ovlbat' and '$path/ovljob' have an inconsistent number of lines.\n"  if (!eof(B) || !eof(J));

    close(J);
    close(B);

    if (scalar(@failedJobs) == 0) {
        open(L, "> $path/ovljob.files") or caFailure("failed to open '$path/ovljob.files'", undef);
        print L @successJobs;
        close(L);
        return;
    }

    if ($attempt > 0) {
        print STDERR "\n";
        print STDERR scalar(@failedJobs), " overlapper jobs failed:\n";
        print STDERR $failureMessage;
        print STDERR "\n";
    }

    print STDERR "overlapCheck() -- attempt $attempt begins with ", scalar(@successJobs), " finished, and ", scalar(@failedJobs), " to compute.\n";

    if ($attempt < 1) {
        submitOrRunParallelJob($wrk, $asm, $jobType, $path, $script, getGlobal("ovlConcurrency"), @failedJobs);
    } else {
        caFailure("failed to overlap.  Made $attempt attempts, jobs still failed.\n", undef);
    }
}

