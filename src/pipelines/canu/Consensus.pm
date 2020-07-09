
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

package canu::Consensus;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(consensusConfigure consensusCheck consensusLoad consensusAnalyze alignGFA);

use strict;
use warnings "all";
no  warnings "uninitialized";

use File::Path 2.08 qw(make_path remove_tree);

use canu::Defaults;
use canu::Configure;

use canu::Execution;
use canu::SequenceStore;
use canu::Unitig;

use canu::Report;

use canu::Grid_Cloud;


sub utgcns ($$$) {
    my $asm     = shift @_;
    my $ctgjobs = shift @_;
    my $utgjobs = shift @_;
    my $jobs    = $ctgjobs + $utgjobs;

    my $path    = "unitigging/5-consensus";

    open(F, "> $path/consensus.sh") or caExit("can't open '$path/consensus.sh' for writing: $!", undef);

    print F "#!" . getGlobal("shell") . "\n";
    print F "\n";
    print F getBinDirectoryShellCode();
    print F "\n";
    print F setWorkDirectoryShellCode($path);
    print F "\n";
    print F getJobIDShellCode();
    print F "\n";
    print F "if [ \$jobid -gt $jobs ]; then\n";
    print F "  echo Error: Only $jobs partitions, you asked for \$jobid.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "\n";
    print F "if [ \$jobid -le $ctgjobs ] ; then\n";
    print F "  tag=\"ctg\"\n";
    print F "else\n";
    print F "  tag=\"utg\"\n";
    print F "  jobid=`expr \$jobid - $ctgjobs`\n";
    print F "fi\n";
    print F "\n";
    print F "jobid=`printf %04d \$jobid`\n";
    print F "\n";
    print F "if [ ! -d ./\${tag}cns ] ; then\n";
    print F "  mkdir -p ./\${tag}cns\n";
    print F "fi\n";
    print F "\n";
    print F "if [ -e ./\${tag}cns/\$jobid.cns ] ; then\n";
    print F "  exit 0\n";
    print F "fi\n";
    print F "\n";
    print F fetchTigStoreShellCode("unitigging/5-consensus", $asm, "\${tag}Store", "001", "");
    print F "\n";
    print F fetchFileShellCode("unitigging/5-consensus", "../$asm.\${tag}Store/partition.\$jobid", "");
    print F "\n";
    print F "\$bin/utgcns \\\n";
    print F "  -R ../$asm.\${tag}Store/partition.\$jobid \\\n";
    print F "  -T ../$asm.\${tag}Store 1 \\\n";
    print F "  -P \$jobid \\\n";
    print F "  -O ./\${tag}cns/\$jobid.cns.WORKING \\\n";
    print F "  -maxcoverage " . getGlobal('cnsMaxCoverage') . " \\\n";
    print F "  -e " . getGlobal("cnsErrorRate") . " \\\n";
    print F "  -quick \\\n"      if (getGlobal("cnsConsensus") eq "quick");
    print F "  -pbdagcon \\\n"   if (getGlobal("cnsConsensus") eq "pbdagcon");
    print F "  -edlib    \\\n"   if (getGlobal("canuIteration") >= 0);
    print F "  -utgcns \\\n"     if (getGlobal("cnsConsensus") eq "utgcns");
    print F "  -threads " . getGlobal("cnsThreads") . " \\\n";
    print F "&& \\\n";
    print F "mv ./\${tag}cns/\$jobid.cns.WORKING ./\${tag}cns/\$jobid.cns \\\n";
    print F "\n";
    print F stashFileShellCode("unitigging/5-consensus", "\${tag}cns/\$jobid.cns", "");
    print F "\n";
    print F "exit 0\n";

    if (getGlobal("canuIteration") < 0) {
        print STDERR "-- Using fast alignment for consensus (iteration '", getGlobal("canuIteration"), "').\n";
    } else {
        print STDERR "-- Using slow alignment for consensus (iteration '", getGlobal("canuIteration"), "').\n";
    }

    close(F);

    makeExecutable("$path/consensus.sh");
    stashFile("$path/consensus.sh");
}



sub partitionTigs ($$) {
    my $asm    = shift @_;
    my $tag    = shift @_;
    my $bin    = getBinDirectory();
    my $cmd;

    return  if (fileExists("unitigging/$asm.${tag}Store/partitioning"));

    fetchTigStore("unitigging", $asm, "${tag}Store", "001");

    # adjust for the fact that compressed contigs will likely expand and thus take more memory/more space
    my $partitionScaling = (defined(getGlobal("homoPolyCompress"))) ? 1.5 : 1.0;

    $cmd  = "$bin/utgcns \\\n";
    $cmd .= "  -S ../$asm.seqStore \\\n";
    $cmd .= "  -T  ./$asm.${tag}Store 1 \\\n";
    #$cmd .= "  -partition " . getGlobal("cnsPartitionSize") . " \\\n"   if (defined(getGlobal("cnsPartitionSize")));
    $cmd .= "  -partition 0.8 $partitionScaling 0.1 \\\n";
    $cmd .= "> ./$asm.${tag}Store/partitioning.log 2>&1";

    if (runCommand("unitigging", $cmd)) {
        caExit("failed to partition the reads", "unitigging/$asm.${tag}Store/partitioning.log");
    }

    stashFile("unitigging/$asm.${tag}Store/partitioning.log");
    stashFile("unitigging/$asm.${tag}Store/partitioning");
}



sub computeNumberOfConsensusJobs ($$) {
    my $asm    = shift @_;
    my $tag    = shift @_;
    my $jobs   = 0;
    my $bin    = getBinDirectory();

    fetchFile("unitigging/$asm.${tag}Store/partitioning");

    open(F, "< unitigging/$asm.${tag}Store/partitioning") or caExit("can't open 'unitigging/$asm.${tag}Store/partitioning' for reading: $!", undef);
    $_ = <F>;
    $_ = <F>;
    while(<F>) {
        s/^\s+//;
        s/\s+$//;

        my @v = split '\s+', $_;

        $jobs = $v[5]   if ($jobs < $v[5]);
    }
    close(F);

    return($jobs);
}



sub estimateMemoryNeededForConsensusJobs ($) {
    my $asm    = shift @_;
    my $minMem = 0;

    fetchFile("unitigging/$asm.ctgStore/partitioning");
    fetchFile("unitigging/$asm.utgStore/partitioning");

    open(F, "< unitigging/$asm.ctgStore/partitioning") or caExit("can't open 'unitigging/$asm.ctgStore/partitioning' for reading: $!", undef);
    $_ = <F>;   #  Header
    $_ = <F>;   #  Header
    while(<F>) {
        s/^\s+//;
        s/\s+$//;

        my @v = split '\s+', $_;

        #  Compute memory needed as 'memory for this tig' + 'memory for read
        #  data' + 0.5 GB, then round to eighths of a GB.  Execution.pm will
        #  massage this into a format acceptable to the scheduler.
        my $thisMem = int(8 * ($v[4] + $v[6] + 0.5)) / 8;

        $minMem = $thisMem   if ($minMem < $thisMem);
    }
    close(F);

    open(F, "< unitigging/$asm.utgStore/partitioning") or caExit("can't open 'unitigging/$asm.utgStore/partitioning' for reading: $!", undef);
    $_ = <F>;   #  Header
    $_ = <F>;   #  Header
    while(<F>) {
        s/^\s+//;
        s/\s+$//;

        #  Compute memory needed as 'memory for this tig' + 'memory for read data' + 0.5 GB, then round to eighths of a GB.
        #  Execution.pm will massage this into a format acceptable to the scheduler.

        my @v = split '\s+', $_;

        my $thisMem = int(8 * ($v[4] + $v[6] + 0.5)) / 8;

        $minMem = $thisMem   if ($minMem < $thisMem);
    }
    close(F);

    #  Warn if the user has overridden our selection.

    my $curMem = getGlobal("cnsMemory");

    if (defined($curMem) && ($curMem > 0)) {
        if ($curMem < $minMem) {
            print STDERR "--\n";
            print STDERR "-- WARNING:\n";
            print STDERR "-- WARNING:  cnsMemory set to $curMem GB, but expected usage is $minMem GB.\n";
            print STDERR "-- WARNING:  Jobs may fail.\n";
            print STDERR "-- WARNING:\n";
        }

    } else {
        setGlobal("cnsMemory", $minMem);

        my $err;
        my $all;

        ($err, $all) = getAllowedResources("", "cns", $err, $all, 0);

        print STDERR "--\n";
        print STDERR $all;
        print STDERR "--\n";
    }

    return($minMem);
}




sub consensusConfigure ($) {
    my $asm    = shift @_;
    my $bin    = getBinDirectory();
    my $cmd;
    my $path   = "unitigging/5-consensus";

    goto allDone   if ((fileExists("unitigging/$asm.ctgStore/seqDB.v002.tig")) &&
                       (fileExists("unitigging/$asm.utgStore/seqDB.v002.tig")));

    make_path($path)  if (! -d $path);

    #  Assign tigs to jobs.

    partitionTigs($asm, "ctg");
    partitionTigs($asm, "utg");

    #  Figure out how many jobs.

    my $ctgjobs = computeNumberOfConsensusJobs($asm, "ctg");
    my $utgjobs = computeNumberOfConsensusJobs($asm, "utg");

    #  This configure is an odd-ball.  Unlike all the other places that write scripts,
    #  we'll rewrite this one every time, so that we can change the alignment algorithm
    #  on the second attempt.

    my $firstTime = (! -e "$path/consensus.sh");

    if ((getGlobal("cnsConsensus") eq "quick") ||
        (getGlobal("cnsConsensus") eq "pbdagcon") ||
        (getGlobal("cnsConsensus") eq "utgcns")) {
        utgcns($asm, $ctgjobs, $utgjobs);

    } else {
        caFailure("unknown consensus style '" . getGlobal("cnsConsensus") . "'", undef);
    }

    print STDERR "-- Configured $ctgjobs contig and $utgjobs unitig consensus jobs.\n";

  finishStage:
    generateReport($asm);
    resetIteration("consensusConfigure")   if ($firstTime);

  allDone:
    stopAfter("consensusConfigure");
}



#  Checks that all consensus jobs are complete, loads them into the store.
#
sub consensusCheck ($) {
    my $asm     = shift @_;
    my $attempt = getGlobal("canuIteration");
    my $path    = "unitigging/5-consensus";

    goto allDone  if ((fileExists("$path/ctgcns.files")) &&
                      (fileExists("$path/utgcns.files")));
    goto allDone  if (fileExists("unitigging/$asm.ctgStore/seqDB.v002.tig"));

    fetchFile("$path/consensus.sh");

    #  Figure out if all the tasks finished correctly.

    my $ctgjobs = computeNumberOfConsensusJobs($asm, "ctg");
    my $utgjobs = computeNumberOfConsensusJobs($asm, "utg");
    my $jobs = $ctgjobs + $utgjobs;

    #  Setup memory and threads and etc.  Complain if not enough memory.

    my $minMem = estimateMemoryNeededForConsensusJobs($asm);

    #  Decide what to run.

    caExit("no consensus jobs found?", undef)   if ($jobs == 0);

    my $currentJobID = "0001";
    my $tag          = "ctgcns";

    my @ctgSuccessJobs;
    my @utgSuccessJobs;
    my @failedJobs;
    my $failureMessage = "";

    for (my $job=1; $job <= $jobs; $job++) {
        if      (fileExists("$path/$tag/$currentJobID.cns")) {
            push @ctgSuccessJobs, "5-consensus/$tag/$currentJobID.cns\n"      if ($tag eq "ctgcns");
            push @utgSuccessJobs, "5-consensus/$tag/$currentJobID.cns\n"      if ($tag eq "utgcns");

        } elsif (fileExists("$path/$tag/$currentJobID.cns.gz")) {
            push @ctgSuccessJobs, "5-consensus/$tag/$currentJobID.cns.gz\n"   if ($tag eq "ctgcns");
            push @utgSuccessJobs, "5-consensus/$tag/$currentJobID.cns.gz\n"   if ($tag eq "utgcns");

        } elsif (fileExists("$path/$tag/$currentJobID.cns.bz2")) {
            push @ctgSuccessJobs, "5-consensus/$tag/$currentJobID.cns.bz2\n"  if ($tag eq "ctgcns");
            push @utgSuccessJobs, "5-consensus/$tag/$currentJobID.cns.bz2\n"  if ($tag eq "utgcns");

        } elsif (fileExists("$path/$tag/$currentJobID.cns.xz")) {
            push @ctgSuccessJobs, "5-consensus/$tag/$currentJobID.cns.xz\n"   if ($tag eq "ctgcns");
            push @utgSuccessJobs, "5-consensus/$tag/$currentJobID.cns.xz\n"   if ($tag eq "utgcns");

        } else {
            $failureMessage .= "--   job $tag/$currentJobID.cns FAILED.\n";
            push @failedJobs, $job;
        }

        $currentJobID++;

        $currentJobID = "0001"    if ($job == $ctgjobs);  #  Reset for first utg job.
        $tag          = "utgcns"  if ($job == $ctgjobs);
    }

    #  Failed jobs, retry.

    if (scalar(@failedJobs) > 0) {

        #  If too many attempts, give up.

        if ($attempt >= getGlobal("canuIterationMax")) {
            print STDERR "--\n";
            print STDERR "-- Consensus jobs failed, tried $attempt times, giving up.\n";
            print STDERR $failureMessage;
            print STDERR "--\n";
            caExit(undef, undef);
        }

        if ($attempt > 0) {
            print STDERR "--\n";
            print STDERR "-- Consensus jobs failed, retry.\n";
            print STDERR $failureMessage;
            print STDERR "--\n";
        }

        #  Otherwise, run some jobs.

        generateReport($asm);

        submitOrRunParallelJob($asm, "cns", $path, "consensus", @failedJobs);
        return;
    }

  finishStage:
    print STDERR "-- All ", scalar(@ctgSuccessJobs) + scalar(@utgSuccessJobs), " consensus jobs finished successfully.\n";

    open(L, "> $path/ctgcns.files") or caExit("can't open '$path/ctgcns.files' for writing: $!", undef);
    print L @ctgSuccessJobs;
    close(L);

    stashFile("$path/ctgcns.files");

    open(L, "> $path/utgcns.files") or caExit("can't open '$path/utgcns.files' for writing: $!", undef);
    print L @utgSuccessJobs;
    close(L);

    stashFile("$path/utgcns.files");

    generateReport($asm);
    resetIteration("consensusCheck");

  allDone:
}



sub purgeFiles ($$$$$$) {
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $Ncns    = shift @_;
    my $Nfastq  = shift @_;
    my $Nlayout = shift @_;
    my $Nlog    = shift @_;

    my $path = "unitigging/5-consensus";

    open(F, "< $path/$tag.files") or caExit("can't open '$path/$tag.files' for reading: $!\n", undef);
    while (<F>) {
        chomp;
        if (m/^(.*)\/0*(\d+).cns$/) {
            my $ID6 = substr("00000" . $2, -6);
            my $ID4 = substr("000"   . $2, -4);
            my $ID0 = $2;

            if (-e "unitigging/$1/$ID4.cns") {
                $Ncns++;
                unlink "unitigging/$1/$ID4.cns";
            }
            if (-e "unitigging/$1/$ID4.fastq") {
                $Nfastq++;
                unlink "unitigging/$1/$ID4.fastq";
            }
            if (-e "unitigging/$1/$ID4.layout") {
                $Nlayout++;
                unlink "unitigging/$1/$ID4.layout";
            }
            if (-e "unitigging/$1/consensus.$ID6.out") {
                $Nlog++;
                unlink "unitigging/$1/consensus.$ID6.out";
            }
            if (-e "unitigging/$1/consensus.$ID0.out") {
                $Nlog++;
                unlink "unitigging/$1/consensus.$ID0.out";
            }

        } else {
            caExit("unknown consensus job name '$_'\n", undef);
        }
    }
    close(F);

    unlink "$path/$tag.files";
    rmdir  "$path/$tag";

    return($Ncns, $Nfastq, $Nlayout, $Nlog);
}



sub consensusLoad ($) {
    my $asm     = shift @_;
    my $bin     = getBinDirectory();
    my $cmd;
    my $path    = "unitigging/5-consensus";

    goto allDone    if ((fileExists("unitigging/$asm.ctgStore/seqDB.v002.tig")) &&
                        (fileExists("unitigging/$asm.utgStore/seqDB.v002.tig")));

    #  Expects to have a list of output files from the consensusCheck() function.

    fetchFile("$path/ctgcns.files");
    fetchFile("$path/utgcns.files");

    caExit("can't find '$path/ctgcns.files' for loading tigs into store: $!", undef)  if (! -e "$path/ctgcns.files");
    caExit("can't find '$path/utgcns.files' for loading tigs into store: $!", undef)  if (! -e "$path/utgcns.files");

    #  Now just load them.

    if (! fileExists("unitigging/$asm.ctgStore/seqDB.v002.tig")) {
        fetchTigStore("unitigging", $asm, "ctgStore", "001");

        open(F, "< $path/ctgcns.files");
        while (<F>) {
            chomp;
            fetchFile("unitigging/$_");
        }
        close(F);

        $cmd  = "$bin/tgStoreLoad \\\n";
        $cmd .= "  -S ../$asm.seqStore \\\n";
        $cmd .= "  -T  ./$asm.ctgStore 2 \\\n";
        $cmd .= "  -L ./5-consensus/ctgcns.files \\\n";
        $cmd .= "> ./5-consensus/ctgcns.files.ctgStoreLoad.err 2>&1";

        if (runCommand("unitigging", $cmd)) {
            caExit("failed to load unitig consensus into ctgStore", "$path/ctgcns.files.ctgStoreLoad.err");
        }
        unlink "$path/ctgcns.files.ctgStoreLoad.err";

        stashFile("unitigging/$asm.ctgStore/seqDB.v002.dat");
        stashFile("unitigging/$asm.ctgStore/seqDB.v002.tig");
    }

    if (! fileExists("unitigging/$asm.utgStore/seqDB.v002.tig")) {
        fetchTigStore("unitigging", $asm, "utgStore", "001");

        open(F, "< $path/utgcns.files");
        while (<F>) {
            chomp;
            fetchFile("unitigging/$_");
        }
        close(F);

        $cmd  = "$bin/tgStoreLoad \\\n";
        $cmd .= "  -S ../$asm.seqStore \\\n";
        $cmd .= "  -T  ./$asm.utgStore 2 \\\n";
        $cmd .= "  -L ./5-consensus/utgcns.files \\\n";
        $cmd .= "> ./5-consensus/utgcns.files.utgStoreLoad.err 2>&1";

        if (runCommand("unitigging", $cmd)) {
            caExit("failed to load unitig consensus into utgStore", "$path/utgcns.files.utgStoreLoad.err");
        }
        unlink "$path/utgcns.files.utgStoreLoad.err";

        stashFile("unitigging/$asm.utgStore/seqDB.v002.dat");
        stashFile("unitigging/$asm.utgStore/seqDB.v002.tig");
    }

    #  Remvoe consensus outputs

    if ((-e "$path/ctgcns.files") ||
        (-e "$path/utgcns.files")) {
        print STDERR "-- Purging consensus output after loading to ctgStore and/or utgStore.\n";

        my $Ncns    = 0;
        my $Nfastq  = 0;
        my $Nlayout = 0;
        my $Nlog    = 0;

        ($Ncns, $Nfastq, $Nlayout, $Nlog) = purgeFiles($asm, "ctgcns", $Ncns, $Nfastq, $Nlayout, $Nlog);
        ($Ncns, $Nfastq, $Nlayout, $Nlog) = purgeFiles($asm, "utgcns", $Ncns, $Nfastq, $Nlayout, $Nlog);

        print STDERR "-- Purged $Ncns .cns outputs.\n"        if ($Ncns > 0);
        print STDERR "-- Purged $Nfastq .fastq outputs.\n"    if ($Nfastq > 0);
        print STDERR "-- Purged $Nlayout .layout outputs.\n"  if ($Nlayout > 0);
        print STDERR "-- Purged $Nlog .err log outputs.\n"    if ($Nlog > 0);
    }

    reportUnitigSizes($asm, 2, "after consensus generation");

  finishStage:
    generateReport($asm);
    resetIteration("consensusLoad");

  allDone:
}




sub consensusAnalyze ($) {
    my $asm     = shift @_;
    my $bin     = getBinDirectory();
    my $cmd;

    #  Left in as a template for any future analysis.

    goto allDone;   #if (fileExists("unitigging/$asm.ctgStore.coverageStat.log"));

    #fetchTigStore("unitigging", $asm, "ctgStore", "001");
    #fetchTigStore("unitigging", $asm, "ctgStore", "002");
    #
    #$cmd  = "$bin/tgStoreCoverageStat \\\n";
    #$cmd .= "  -S ../$asm.seqStore \\\n";
    #$cmd .= "  -T  ./$asm.ctgStore 2 \\\n";
    #$cmd .= "  -s " . getGlobal("genomeSize") . " \\\n";
    #$cmd .= "  -o ./$asm.ctgStore.coverageStat \\\n";
    #$cmd .= "> ./$asm.ctgStore.coverageStat.err 2>&1";
    #
    #if (runCommand("unitigging", $cmd)) {
    #    caExit("failed to compute coverage statistics", "unitigging/$asm.ctgStore.coverageStat.err");
    #}
    #
    #unlink "unitigging/$asm.ctgStore.coverageStat.err";
    #
    #stashFile("unitigging/$asm.ctgStore.coverageStat.stats");
    #stashFile("unitigging/$asm.ctgStore.coverageStat.log");

  finishStage:
    generateReport($asm);
    resetIteration("consensusAnalyze");

  allDone:
    stopAfter("consensus");
}




sub alignGFA ($) {
    my $asm     = shift @_;
    my $attempt = getGlobal("canuIteration");
    my $path    = "unitigging/4-unitigger";

    #  Decide if this is small enough to run right now, or if we should submit to the grid.

    #my $bin     = getBinDirectory();

    #  This is just big enough to not fit comfortably in the canu process itself.

    goto allDone   if (fileExists("unitigging/4-unitigger/$asm.unitigs.aligned.gfa") &&
                       fileExists("unitigging/4-unitigger/$asm.unitigs.aligned.bed"));

    #  If a large genome, run this on the grid, else, run in the canu process itself.
    my $runGrid = (getGlobal("genomeSize") >= 40000000);

    make_path($path);                  #  In cloud mode, 4-unitigger doesn't exist when we get here.
    fetchFile("$path/alignGFA.sh");

    if (! -e "$path/alignGFA.sh") {
        open(F, "> $path/alignGFA.sh") or caExit("can't open '$path/alignGFA.sh' for writing: $!\n", undef);
        print F "#!" . getGlobal("shell") . "\n";
        print F "\n";
        print F getBinDirectoryShellCode();
        print F "\n";
        print F setWorkDirectoryShellCode($path);
        print F "\n";
        print F fetchTigStoreShellCode("unitigging/4-unitigger", $asm, "utgStore", "001", "");
        print F "\n";
        print F fetchTigStoreShellCode("unitigging/4-unitigger", $asm, "utgStore", "002", "");
        print F "\n";
        print F fetchTigStoreShellCode("unitigging/4-unitigger", $asm, "ctgStore", "001", "");
        print F "\n";
        print F fetchTigStoreShellCode("unitigging/4-unitigger", $asm, "ctgStore", "002", "");
        print F "\n";

        print F "if [ ! -e ./$asm.unitigs.aligned.gfa ] ; then\n";
        print F    fetchFileShellCode("unitigging/4-unitigger", "$asm.unitigs.gfa", "  ");
        print F "\n";
        print F "  \$bin/alignGFA \\\n";
        print F "    -T ../$asm.utgStore 2 \\\n";
        print F "    -i ./$asm.unitigs.gfa \\\n";
        print F "    -e " . getGlobal("cnsErrorRate") . " \\\n";
        print F "    -o ./$asm.unitigs.aligned.gfa \\\n";
        print F "    -t " . getGlobal("gfaThreads") . " \\\n";
        print F "  > ./$asm.unitigs.aligned.gfa.err 2>&1";
        print F "\n";
        print F    stashFileShellCode("$path", "$asm.unitigs.aligned.gfa", "  ");
        print F "fi\n";
        print F "\n";
        print F "\n";

        print F "if [ ! -e ./$asm.unitigs.aligned.bed ] ; then\n";
        print F    fetchFileShellCode("unitigging/4-unitigger", "$asm.unitigs.bed", "  ");
        print F "\n";
        print F "  \$bin/alignGFA -bed \\\n";
        print F "    -T ../$asm.utgStore 2 \\\n";
        print F "    -C ../$asm.ctgStore 2 \\\n";
        print F "    -i ./$asm.unitigs.bed \\\n";
        print F "    -o ./$asm.unitigs.aligned.bed \\\n";
        print F "    -t " . getGlobal("gfaThreads") . " \\\n";
        print F "  > ./$asm.unitigs.aligned.bed.err 2>&1";
        print F "\n";
        print F    stashFileShellCode("$path", "$asm.unitigs.aligned.bed", "  ");
        print F "fi\n";
        print F "\n";
        print F "\n";

        print F "if [ -e ./$asm.unitigs.aligned.gfa -a \\\n";
        print F "     -e ./$asm.unitigs.aligned.bed ] ; then\n";
        print F "  echo GFA alignments updated.\n";
        print F "  exit 0\n";
        print F "else\n";
        print F "  echo GFA alignments failed.\n";
        print F "  exit 1\n";
        print F "fi\n";
        close(F);

        makeExecutable("$path/alignGFA.sh");
        stashFile("$path/alignGFA.sh");
    }

    #  Since there is only one job, if we get here, we're not done.  Any other 'check' function
    #  shows how to process multiple jobs.  This only checks for the existence of the final outputs.
    #  (meryl-process and unitig are the same)

    #  If too many attempts, give up.

    if ($attempt >= getGlobal("canuIterationMax")) {
        print STDERR "--\n";
        print STDERR "-- Graph alignment jobs failed, tried $attempt times, giving up.\n";
        print STDERR "--\n";
        caExit(undef, undef);
    }

    if ($attempt > 0) {
        print STDERR "--\n";
        print STDERR "-- Graph alignment jobs failed, retry.\n";
        print STDERR "--\n";
    }

    #  Otherwise, run some jobs.

    generateReport($asm);

    if ($runGrid) {
        submitOrRunParallelJob($asm, "gfa", $path, "alignGFA", (1));
    } else {
        if (runCommand($path, "./alignGFA.sh > alignGFA.err 2>&1")) {
            caExit("failed to create aligned graphs", "./alignGFA.err");
        }
    }

    return;

  finishStage:
    generateReport($asm);
    resetIteration("alignGFA");

  allDone:
}
