
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
 #    src/pipelines/ca3g/OverlapInCore.pm
 #
 #  Modifications by:
 #
 #    Brian P. Walenz from 2015-FEB-27 to 2015-AUG-25
 #      are Copyright 2015 Battelle National Biodefense Institute, and
 #      are subject to the BSD 3-Clause License
 #
 #    Brian P. Walenz beginning on 2015-OCT-19
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

package canu::OverlapInCore;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(overlapConfigure overlap overlapCheck);

use strict;

use File::Path qw(make_path remove_tree);

use canu::Defaults;
use canu::Execution;
use canu::HTML;


sub overlapConfigure ($$$$) {
    my $WRK     = shift @_;  #  Root work directory (the -d option to canu)
    my $wrk     = $WRK;      #  Local work directory
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $type    = shift @_;

    my $bin  = getBinDirectory();
    my $cmd;

    $wrk = "$wrk/correction"  if ($tag eq "cor");
    $wrk = "$wrk/trimming"    if ($tag eq "obt");
    $wrk = "$wrk/unitigging"  if ($tag eq "utg");

    my $path = "$wrk/1-overlapper";

    caFailure("invalid type '$type'", undef)  if (($type ne "partial") && ($type ne "normal"));

    goto allDone   if (skipStage($WRK, $asm, "$tag-overlapConfigure") == 1);
    goto allDone   if (-e "$path/$asm.partition.ovlopt");
    goto allDone   if (-e "$path/ovljob.files");
    goto allDone   if (-e "$wrk/$asm.ovlStore");

    print STDERR "--\n";
    print STDERR "-- OVERLAPPER (normal) (correction) erate=", getGlobal("${tag}OvlErrorRate"), "\n"  if ($tag eq "cor");
    print STDERR "-- OVERLAPPER (normal) (trimming) erate=", getGlobal("${tag}OvlErrorRate"), "\n"    if ($tag eq "obt");
    print STDERR "-- OVERLAPPER (normal) (assembly) erate=", getGlobal("${tag}OvlErrorRate"), "\n"    if ($tag eq "utg");
    print STDERR "--\n";

    make_path("$path") if (! -d "$path");

    if (! -e "$path/$asm.partition.ovlopt") {

        #  These used to be runCA options, but were removed in canu.  They were used mostly for illumina-pacbio correction,
        #  but were also used (or could have been used) during the Salmon assembly when overlaps were computed differently
        #  depending on the libraries involved (and was run manually).  These are left in for documentation.
        #
        #my $checkLibrary       = getGlobal("${tag}CheckLibrary");
        #my $hashLibrary        = getGlobal("${tag}HashLibrary");
        #my $refLibrary         = getGlobal("${tag}RefLibrary");

        my $hashBlockLength = getGlobal("${tag}OvlHashBlockLength");
        my $hashBlockSize   = 0;
        my $refBlockSize    = getGlobal("${tag}OvlRefBlockSize");
        my $refBlockLength  = getGlobal("${tag}OvlRefBlockLength");
        my $minOlapLength   = getGlobal("minOverlapLength");

        if (($refBlockSize > 0) && ($refBlockLength > 0)) {
            caExit("can't set both ${tag}OvlRefBlockSize and ${tag}OvlRefBlockLength", undef);
        }

        $cmd  = "$bin/overlapInCorePartition \\\n";
        $cmd .= " -g  $wrk/$asm.gkpStore \\\n";
        $cmd .= " -bl $hashBlockLength \\\n";
        $cmd .= " -bs $hashBlockSize \\\n";
        $cmd .= " -rs $refBlockSize \\\n";
        $cmd .= " -rl $refBlockLength \\\n";
        #$cmd .= " -H $hashLibrary \\\n" if ($hashLibrary ne "0");
        #$cmd .= " -R $refLibrary \\\n"  if ($refLibrary ne "0");
        #$cmd .= " -C \\\n" if (!$checkLibrary);
        $cmd .= " -ol $minOlapLength \\\n";
        $cmd .= " -o  $path/$asm.partition \\\n";
        $cmd .= "> $path/$asm.partition.err 2>&1";

        if (runCommand($wrk, $cmd)) {
            caExit("failed partition for overlapper", undef);
        }
    }

    open(BAT, "< $path/$asm.partition.ovlbat") or caExit("can't open '$path/$asm.partition.ovlbat' for reading: $!", undef);
    open(JOB, "< $path/$asm.partition.ovljob") or caExit("can't open '$path/$asm.partition.ovljob' for reading: $!", undef);
    open(OPT, "< $path/$asm.partition.ovlopt") or caExit("can't open '$path/$asm.partition.ovlopt' for reading: $!", undef);

    my @bat = <BAT>;  chomp @bat;
    my @job = <JOB>;  chomp @job;
    my @opt = <OPT>;  chomp @opt;

    close(BAT);
    close(JOB);
    close(OPT);

    #getAllowedResources($tag, "ovl");

    if (! -e "$path/overlap.sh") {
        my $merSize      = getGlobal("${tag}OvlMerSize");

        #my $hashLibrary  = getGlobal("${tag}OvlHashLibrary");
        #my $refLibrary   = getGlobal("${tag}OvlRefLibrary");

        #  Create a script to run overlaps.  We make a giant job array for this -- we need to know
        #  hashBeg, hashEnd, refBeg and refEnd -- from that we compute batchName and jobName.

        my $hashBits       = getGlobal("${tag}OvlHashBits");
        my $hashLoad       = getGlobal("${tag}OvlHashLoad");

        open(F, "> $path/overlap.sh") or caExit("can't open '$path/overlap.sh' for writing: $!", undef);
        print F "#!" . getGlobal("shell") . "\n";
        print F "\n";
        print F "perl='/usr/bin/env perl'\n";
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

        for (my $ii=1; $ii<=scalar(@bat); $ii++) {
            print F "if [ \$jobid -eq $ii ] ; then\n";
            print F "  bat=\"$bat[$ii-1]\"\n";            #  Needed to mkdir.
            print F "  job=\"$bat[$ii-1]/$job[$ii-1]\"\n";  #  Needed to simplify overlapCheck() below.
            print F "  opt=\"$opt[$ii-1]\"\n";
            print F "fi\n";
            print F "\n";
        }

        print F "\n";
        print F "if [ ! -d $path/\$bat ]; then\n";
        print F "  mkdir $path/\$bat\n";
        print F "fi\n";
        print F "\n";
        print F "if [ -e $path/\$job.ovb.gz ]; then\n";
        print F "  echo Job previously completed successfully.\n";
        print F "  exit\n";
        print F "fi\n";
        print F "\n";
        print F getBinDirectoryShellCode();
        print F "\n";
        print F "\$bin/overlapInCore \\\n";
        print F "  -G \\\n"  if ($type eq "partial");
        print F "  -t ", getGlobal("${tag}OvlThreads"), " \\\n";
        print F "  -k $merSize \\\n";
        print F "  -k $wrk/0-mercounts/$asm.ms$merSize.frequentMers.fasta \\\n";
        print F "  --hashbits $hashBits \\\n";
        print F "  --hashload $hashLoad \\\n";
        print F "  --maxerate  ", getGlobal("${tag}OvlErrorRate"), " \\\n";
        print F "  --minlength ", getGlobal("minOverlapLength"), " \\\n";
        print F "  \$opt \\\n";
        print F "  -o $path/\$job.ovb.WORKING.gz \\\n";
        print F "  -s $path/\$job.stats \\\n";
        #print F "  -H $hashLibrary \\\n" if ($hashLibrary ne "0");
        #print F "  -R $refLibrary \\\n"  if ($refLibrary  ne "0");
        print F "  $wrk/$asm.gkpStore \\\n";
        print F "&& \\\n";
        print F "mv $path/\$job.ovb.WORKING.gz $path/\$job.ovb.gz\n";
        print F "\n";
        print F "exit 0\n";
        close(F);

        system("chmod +x $path/overlap.sh");
    }

    my $jobs      = scalar(@job);
    my $batchName = $bat[$jobs-1];  chomp $batchName;
    my $jobName   = $job[$jobs-1];  chomp $jobName;

    if (-e "$path/overlap.sh") {
        my $numJobs = 0;
        open(F, "< $path/overlap.sh") or caExit("can't open '$path/overlap.sh' for reading: $!", undef);
        while (<F>) {
            $numJobs++  if (m/^\s+job=/);
        }
        close(F);
        print STDERR "--\n";
        print STDERR "-- Configured $numJobs overlapInCore jobs.\n";
    }

  finishStage:
    emitStage($WRK, $asm, "$tag-overlapConfigure");
    buildHTML($WRK, $asm, $tag);
    stopAfter("overlapConfigure");

  allDone:
}



sub reportSumMeanStdDev (@) {
    my $sum    = 0;
    my $mean   = 0;
    my $stddev = 0;
    my $n      = scalar(@_);
    my $formatted;

    $sum += $_  foreach (@_);

    $mean = $sum / $n;

    $stddev += ($_ - $mean) * ($_ - $mean)  foreach (@_);
    $stddev  = int(1000 * sqrt($stddev / ($n-1))) / 1000  if ($n > 1);

    $formatted  = substr("         $mean", -10) . " +- $stddev";

    return($sum, $formatted);
}


sub reportOverlapStats ($$@) {
    my $wrk       = shift @_;  #  Local work directory
    my $asm       = shift @_;
    my @statsJobs = @_;

    my @hitsWithoutOlaps;
    my @hitsWithOlaps;
    my @multiOlaps;
    my @totalOlaps;
    my @containedOlaps;
    my @dovetailOlaps;
    my @shortReject;
    my @longReject;

    foreach my $s (@statsJobs) {
        open(F, "< $s") or caExit("can't open '$s' for reading: $!", undef);

        $_ = <F>;  push @hitsWithoutOlaps, $1  if (m/^\s*Kmer\shits\swithout\solaps\s=\s(\d+)$/);
        $_ = <F>;  push @hitsWithOlaps, $1     if (m/^\s*Kmer\shits\swith\solaps\s=\s(\d+)$/);
        $_ = <F>;  push @multiOlaps, $1        if (m/^\s*Multiple\soverlaps\/pair\s=\s(\d+)$/);
        $_ = <F>;  push @totalOlaps, $1        if (m/^\s*Total\soverlaps\sproduced\s=\s(\d+)$/);
        $_ = <F>;  push @containedOlaps, $1    if (m/^\s*Contained\soverlaps\s=\s(\d+)$/);
        $_ = <F>;  push @dovetailOlaps, $1     if (m/^\s*Dovetail\soverlaps\s=\s(\d+)$/);
        $_ = <F>;  push @shortReject, $1       if (m/^\s*Rejected\sby\sshort\swindow\s=\s(\d+)$/);
        $_ = <F>;  push @longReject, $1        if (m/^\s*Rejected\sby\slong\swindow\s=\s(\d+)$/);

        close(F);
    }

    printf STDERR "--\n";
    printf STDERR "-- overlapInCore compute '$wrk/1-overlapper':\n";
    printf STDERR "--   kmer hits\n";
    printf STDERR "--     with no overlap     %12d  %s\n", reportSumMeanStdDev(@hitsWithoutOlaps);
    printf STDERR "--     with an overlap     %12d  %s\n", reportSumMeanStdDev(@hitsWithOlaps);
    printf STDERR "--\n";
    printf STDERR "--   overlaps              %12d  %s\n", reportSumMeanStdDev(@totalOlaps);
    printf STDERR "--     contained           %12d  %s\n", reportSumMeanStdDev(@containedOlaps);
    printf STDERR "--     dovetail            %12d  %s\n", reportSumMeanStdDev(@dovetailOlaps);
    printf STDERR "--\n";
    printf STDERR "--   overlaps rejected\n";
    printf STDERR "--     multiple per pair   %12d  %s\n", reportSumMeanStdDev(@multiOlaps);
    printf STDERR "--     bad short window    %12d  %s\n", reportSumMeanStdDev(@shortReject);
    printf STDERR "--     bad long window     %12d  %s\n", reportSumMeanStdDev(@longReject);
}


#  Check that the overlapper jobs properly executed.  If not,
#  complain, but don't help the user fix things.
#
sub overlapCheck ($$$$) {
    my $WRK     = shift @_;  #  Root work directory (the -d option to canu)
    my $wrk     = $WRK;      #  Local work directory
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $type    = shift @_;
    my $attempt = getGlobal("canuIteration");

    $wrk = "$wrk/correction"  if ($tag eq "cor");
    $wrk = "$wrk/trimming"    if ($tag eq "obt");
    $wrk = "$wrk/unitigging"  if ($tag eq "utg");

    my $path    = "$wrk/1-overlapper";

    goto allDone   if (skipStage($WRK, $asm, "$tag-overlapCheck", $attempt) == 1);
    goto allDone   if (-e "$path/ovljob.files");
    goto allDone   if (-e "$wrk/$asm.ovlStore");

    #  Figure out if all the tasks finished correctly.

    my $currentJobID   = 1;
    my @successJobs;
    my @statsJobs;
    my @miscJobs;
    my @failedJobs;
    my $failureMessage = "";

    open(F, "< $path/overlap.sh") or caExit("can't open '$path/overlap.sh' for reading: $!", undef);

    while (<F>) {
        if (m/^\s+job=\"(\d+\/\d+)\"$/) {
            if      (-e "$path/$1.ovb.gz") {
                push @successJobs, "$path/$1.ovb.gz\n";   #  Dumped to a file, so include \n
                push @statsJobs,   "$path/$1.stats";      #  Used here, don't include \n
                push @miscJobs,    "$path/$1.stats\n";
                push @miscJobs,    "$path/$1.counts\n";

            } elsif (-e "$path/$1.ovb") {
                push @successJobs, "$path/$1.ovb\n";
                push @statsJobs,   "$path/$1.stats";
                push @miscJobs,    "$path/$1.stats\n";
                push @miscJobs,    "$path/$1.counts\n";

            } elsif (-e "$path/$1.ovb.bz2") {
                push @successJobs, "$path/$1.ovb.bz2\n";
                push @statsJobs,   "$path/$1.stats";
                push @miscJobs,    "$path/$1.stats\n";
                push @miscJobs,    "$path/$1.counts\n";

            } elsif (-e "$path/$1.ovb.xz") {
                push @successJobs, "$path/$1.ovb.xz\n";
                push @statsJobs,   "$path/$1.stats";
                push @miscJobs,    "$path/$1.stats\n";
                push @miscJobs,    "$path/$1.counts\n";

            } else {
                $failureMessage .= "--   job $path/$1 FAILED.\n";
                push @failedJobs, $currentJobID;
            }

            $currentJobID++;
        }
    }

    close(F);

    #  Failed jobs, retry.

    if (scalar(@failedJobs) > 0) {

        #  If not the first attempt, report the jobs that failed, and that we're recomputing.

        if ($attempt > 1) {
            print STDERR "--\n";
            print STDERR "-- ", scalar(@failedJobs), " overlapper jobs failed:\n";
            print STDERR $failureMessage;
            print STDERR "--\n";
        }

        #  If too many attempts, give up.

        if ($attempt > getGlobal("canuIterationMax")) {
            caExit("failed to overlap.  Made " . ($attempt-1) . " attempts, jobs still failed", undef);
        }

        #  Otherwise, run some jobs.

        print STDERR "-- overlapInCore attempt $attempt begins with ", scalar(@successJobs), " finished, and ", scalar(@failedJobs), " to compute.\n";

        emitStage($WRK, $asm, "$tag-overlapCheck", $attempt);
        buildHTML($WRK, $asm, $tag);

        submitOrRunParallelJob($WRK, $asm, "${tag}ovl", $path, "overlap", @failedJobs);
        return;
    }

  finishStage:
    print STDERR "-- Found ", scalar(@successJobs), " overlapInCore output files.\n";

    open(L, "> $path/ovljob.files") or caExit("can't open '$path/ovljob.files' for writing: $!", undef);
    print L @successJobs;
    close(L);

    open(L, "> $path/ovljob.more.files") or caExit("can't open '$path/ovljob.more.files' for writing: $!", undef);
    print L @miscJobs;
    close(L);

    reportOverlapStats($wrk, $asm, @statsJobs);

    setGlobal("canuIteration", 0);
    emitStage($WRK, $asm, "$tag-overlapCheck");
    buildHTML($WRK, $asm, $tag);
    stopAfter("overlapper");

  allDone:
}

