
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

package canu::OverlapInCore;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(overlapConfigure overlapCheck);

use strict;
use warnings "all";
no  warnings "uninitialized";

use File::Path 2.08 qw(make_path remove_tree);

use canu::Defaults;
use canu::Execution;

use canu::Report;

use canu::Grid_Cloud;


sub overlapConfigure ($$$) {
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $type    = shift @_;

    my $bin  = getBinDirectory();
    my $cmd;

    my $base;
    my $path;

    $base = "correction"  if ($tag eq "cor");
    $base = "trimming"    if ($tag eq "obt");
    $base = "unitigging"  if ($tag eq "utg");

    $path = "$base/1-overlapper";

    caFailure("invalid type '$type'", undef)  if (($type ne "partial") && ($type ne "normal"));

    goto allDone   if (fileExists("$path/overlap.sh") && fileExists("$path/$asm.partition.ovlbat") && fileExists("$path/$asm.partition.ovljob") && fileExists("$path/$asm.partition.ovlopt"));
    goto allDone   if (fileExists("$path/ovljob.files"));
    goto allDone   if ((-d "$base/$asm.ovlStore") || (fileExists("$base/$asm.ovlStore.tar.gz")));

    print STDERR "--\n";
    print STDERR "-- OVERLAPPER (normal) (correction) erate=", getGlobal("corOvlErrorRate"), "\n"  if ($tag eq "cor");
    print STDERR "-- OVERLAPPER (normal) (trimming) erate=",   getGlobal("obtOvlErrorRate"), "\n"    if ($tag eq "obt");
    print STDERR "-- OVERLAPPER (normal) (assembly) erate=",   getGlobal("utgOvlErrorRate"), "\n"    if ($tag eq "utg");
    print STDERR "--\n";

    make_path("$path") if (! -d "$path");

    #  overlapInCorePartition internally uses 'WORKING' outputs, and renames to the final
    #  version right before it exits.  All we need to do here is check for existence of
    #  the output, and exit if the command fails.

    fetchFile("$path/$asm.partition.ovlbat");
    fetchFile("$path/$asm.partition.ovljob");
    fetchFile("$path/$asm.partition.ovlopt");
    fetchFile("$path/overlap.sh");              #  Fetch early, so we can delete if ovlbat etc don't exist.

    if ((! -e "$path/$asm.partition.ovlbat") ||
        (! -e "$path/$asm.partition.ovljob") ||
        (! -e "$path/$asm.partition.ovlopt")) {

        #  These used to be runCA options, but were removed in canu.  They were used mostly for
        #  illumina-pacbio correction, but were also used (or could have been used) during the
        #  Salmon assembly when overlaps were computed differently depending on the libraries
        #  involved (and was run manually).  These are left in for documentation.
        #
        #my $checkLibrary       = getGlobal("${tag}CheckLibrary");    #  -C
        #my $hashLibrary        = getGlobal("${tag}HashLibrary");     #  -H $hashLibrary
        #my $refLibrary         = getGlobal("${tag}RefLibrary");      #  -R $refLibrary

        $cmd  = "$bin/overlapInCorePartition \\\n";
        $cmd .= " -S  ../../$asm.seqStore \\\n";
        $cmd .= " -hl " . getGlobal("${tag}OvlHashBlockLength") . " \\\n";
        $cmd .= " -rl " . getGlobal("${tag}OvlRefBlockLength")  . " \\\n";
        $cmd .= " -ol " . getGlobal("minOverlapLength") . " \\\n";
        $cmd .= " -o  ./$asm.partition \\\n";
        $cmd .= "> ./$asm.partition.err 2>&1";

        if (runCommand($path, $cmd)) {
            caExit("failed partition for overlapper", undef);
        }

        stashFile("$path/$asm.partition.ovlbat");
        stashFile("$path/$asm.partition.ovljob");
        stashFile("$path/$asm.partition.ovlopt");

        unlink "$path/overlap.sh";
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
        print F getBinDirectoryShellCode();
        print F "\n";
        print F setWorkDirectoryShellCode($path);
        print F fetchSeqStoreShellCode($asm, $path, "");
        print F "\n";
        print F getJobIDShellCode();
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
        print F "if [ ! -d ./\$bat ]; then\n";
        print F "  mkdir ./\$bat\n";
        print F "fi\n";
        print F "\n";
        print F fileExistsShellCode("exists", $path, "\$job.ovb");
        print F "if [ \$exists = true ] ; then\n";
        print F "  echo Job previously completed successfully.\n";
        print F "  exit\n";
        print F "fi\n";
        print F "\n";
        print F "#  Fetch the frequent kmers, if needed.\n";
        print F "if [ ! -e ../0-mercounts/$asm.ms$merSize.dump ] ; then\n";
        print F "  mkdir -p ../0-mercounts\n";
        print F "  cd ../0-mercounts\n";
        print F fetchFileShellCode("$base/0-mercounts", "$asm.ms$merSize.dump", "  ");
        print F "  cd -\n";
        print F "fi\n";
        print F "\n";
        print F "\n";
        print F "\$bin/overlapInCore \\\n";
        print F "  -partial \\\n"  if ($type eq "partial");
        print F "  -t ", getGlobal("${tag}OvlThreads"), " \\\n";
        print F "  -k $merSize \\\n";
        print F "  -k ../0-mercounts/$asm.ms$merSize.dump \\\n";
        print F "  --hashbits $hashBits \\\n";
        print F "  --hashload $hashLoad \\\n";
        print F "  --maxerate  ", getGlobal("corOvlErrorRate"), " \\\n"  if ($tag eq "cor");   #  Explicitly using proper name for grepability.
        print F "  --maxerate  ", getGlobal("obtOvlErrorRate"), " \\\n"  if ($tag eq "obt");
        print F "  --maxerate  ", getGlobal("utgOvlErrorRate"), " \\\n"  if ($tag eq "utg");
        print F "  --minlength ", getGlobal("minOverlapLength"), " \\\n";
        print F "  --minkmers \\\n" if (defined(getGlobal("${tag}OvlFilter")) && getGlobal("${tag}OvlFilter")==1);
        print F "  \$opt \\\n";
        print F "  -o ./\$job.ovb.WORKING \\\n";
        print F "  -s ./\$job.stats \\\n";
        #print F "  -H $hashLibrary \\\n" if ($hashLibrary ne "0");
        #print F "  -R $refLibrary \\\n"  if ($refLibrary  ne "0");
        print F "  ../../$asm.seqStore \\\n";
        print F "&& \\\n";
        print F "mv ./\$job.ovb.WORKING ./\$job.ovb\n";
        print F "\n";
        print F stashFileShellCode("$base/1-overlapper/", "\$job.ovb",   "");
        print F stashFileShellCode("$base/1-overlapper/", "\$job.oc",    "");
        print F stashFileShellCode("$base/1-overlapper/", "\$job.stats", "");
        print F "\n";
        print F "exit 0\n";
        close(F);

        makeExecutable("$path/overlap.sh");
        stashFile("$path/overlap.sh");
    }

    my $jobs      = scalar(@job);
    my $batchName = $bat[$jobs-1];  chomp $batchName;
    my $jobName   = $job[$jobs-1];  chomp $jobName;

    my $numJobs = 0;
    open(F, "< $path/overlap.sh") or caExit("can't open '$path/overlap.sh' for reading: $!", undef);
    while (<F>) {
        $numJobs++  if (m/^\s+job=/);
    }
    close(F);
    print STDERR "--\n";
    print STDERR "-- Configured $numJobs overlapInCore jobs.\n";

  finishStage:
    generateReport($asm);
    resetIteration("$tag-overlapConfigure");

  allDone:
    stopAfter("overlapConfigure");
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
    my $base      = shift @_;
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
        fetchFile("$base/$s");

        open(F, "< $base/$s") or caExit("can't open '$base/$s' for reading: $!", undef);

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
    printf STDERR "-- overlapInCore compute '$base/1-overlapper':\n";
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
sub overlapCheck ($$$) {
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $type    = shift @_;
    my $attempt = getGlobal("canuIteration");

    my $base;
    my $path;

    $base = "correction"  if ($tag eq "cor");
    $base = "trimming"    if ($tag eq "obt");
    $base = "unitigging"  if ($tag eq "utg");

    $path = "$base/1-overlapper";

    goto allDone   if (fileExists("$path/ovljob.files"));
    goto allDone   if ((-d "$base/$asm.ovlStore") || (fileExists("$base/$asm.ovlStore.tar.gz")));

    #  Figure out if all the tasks finished correctly.

    my $currentJobID   = 1;
    my @successJobs;
    my @statsJobs;
    my @miscJobs;
    my @failedJobs;
    my $failureMessage = "";

    fetchFile("$path/overlap.sh");

    open(F, "< $path/overlap.sh") or caExit("can't open '$path/overlap.sh' for reading: $!", undef);

    while (<F>) {
        if (m/^\s+job=\"(\d+\/\d+)\"$/) {
            if      (fileExists("$path/$1.ovb.gz")) {
                push @successJobs, "1-overlapper/$1.ovb.gz\n";   #  Dumped to a file, so include \n
                push @statsJobs,   "1-overlapper/$1.stats";      #  Used here, don't include \n
                push @miscJobs,    "1-overlapper/$1.stats\n";
                push @miscJobs,    "1-overlapper/$1.oc\n";

            } elsif (fileExists("$path/$1.ovb")) {
                push @successJobs, "1-overlapper/$1.ovb\n";
                push @statsJobs,   "1-overlapper/$1.stats";
                push @miscJobs,    "1-overlapper/$1.stats\n";
                push @miscJobs,    "1-overlapper/$1.oc\n";

            } elsif (fileExists("$path/$1.ovb.bz2")) {
                push @successJobs, "1-overlapper/$1.ovb.bz2\n";
                push @statsJobs,   "1-overlapper/$1.stats";
                push @miscJobs,    "1-overlapper/$1.stats\n";
                push @miscJobs,    "1-overlapper/$1.oc\n";

            } elsif (fileExists("$path/$1.ovb.xz")) {
                push @successJobs, "1-overlapper/$1.ovb.xz\n";
                push @statsJobs,   "1-overlapper/$1.stats";
                push @miscJobs,    "1-overlapper/$1.stats\n";
                push @miscJobs,    "1-overlapper/$1.oc\n";

            } else {
                $failureMessage .= "--   job $path/$1.ovb FAILED.\n";
                push @failedJobs, $currentJobID;
            }

            $currentJobID++;
        }
    }

    close(F);

    #  Failed jobs, retry.

    if (scalar(@failedJobs) > 0) {

        #  If too many attempts, give up.

        if ($attempt >= getGlobal("canuIterationMax")) {
            print STDERR "--\n";
            print STDERR "-- Overlap jobs failed, tried $attempt times, giving up.\n";
            print STDERR $failureMessage;
            print STDERR "--\n";
            caExit(undef, undef);
        }

        if ($attempt > 0) {
            print STDERR "--\n";
            print STDERR "-- Overlap jobs failed, retry.\n";
            print STDERR $failureMessage;
            print STDERR "--\n";
        }

        #  Otherwise, run some jobs.

        generateReport($asm);

        submitOrRunParallelJob($asm, "${tag}ovl", $path, "overlap", @failedJobs);
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

    stashFile("$path/ovljob.files");
    stashFile("$path/ovljob.more.files");

    reportOverlapStats($base, $asm, @statsJobs);

    generateReport($asm);
    resetIteration("$tag-overlapCheck");

  allDone:
    stopAfter("overlap");
}

