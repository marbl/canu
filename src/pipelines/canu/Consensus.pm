
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
@EXPORT = qw(consensusConfigure consensusCheck consensusLoad consensusAnalyze);

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


sub utgcns ($$) {
    my $asm     = shift @_;
    my $jobs    = shift @_;

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
    print F "jobid=`printf %04d \$jobid`\n";
    print F "\n";
    print F "if [ ! -d ./ctgcns ] ; then\n";
    print F "  mkdir -p ./ctgcns\n";
    print F "fi\n";
    print F "\n";
    print F "if [ -e ./ctgcns/\$jobid.cns ] ; then\n";
    print F "  exit 0\n";
    print F "fi\n";
    print F "\n";
    print F fetchTigStoreShellCode("unitigging/5-consensus", $asm, "ctgStore", "001", "");
    print F "\n";
    print F fetchFileShellCode("unitigging/5-consensus", "../$asm.ctgStore/partition.\$jobid", "");
    print F "\n";
    print F "\$bin/utgcns \\\n";
    print F "  -R ../$asm.ctgStore/partition.\$jobid \\\n";
    print F "  -T ../$asm.ctgStore 1 \\\n";
    print F "  -P \$jobid \\\n";
    print F "  -O ./ctgcns/\$jobid.cns.WORKING \\\n";
    print F "  -maxcoverage " . getGlobal('cnsMaxCoverage') . " \\\n";
    print F "  -e " . getGlobal("cnsErrorRate") . " \\\n";
    print F "  -quick \\\n"      if (getGlobal("cnsConsensus") eq "quick");
    print F "  -pbdagcon \\\n"   if (getGlobal("cnsConsensus") eq "pbdagcon");
    print F "  -edlib    \\\n"   if (getGlobal("canuIteration") >= 0);
    print F "  -utgcns \\\n"     if (getGlobal("cnsConsensus") eq "utgcns");
    print F "  -threads " . getGlobal("cnsThreads") . " \\\n";
    print F "&& \\\n";
    print F "mv ./ctgcns/\$jobid.cns.WORKING ./ctgcns/\$jobid.cns \\\n";
    print F "\n";
    print F stashFileShellCode("unitigging/5-consensus", "ctgcns/\$jobid.cns", "");
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



sub partitionTigs ($) {
    my $asm    = shift @_;
    my $bin    = getBinDirectory();
    my $cmd;

    return  if (fileExists("unitigging/$asm.ctgStore/partitioning"));

    fetchTigStore("unitigging", $asm, "ctgStore", "001");

    my $pSize  = 0.8;   #  Make partitions be about 80% the size of the biggest tig.
    my $pScale = 1.0;   #  Don't pre-scale tig lengths (unless homopoly compressed, below).
    my $pReads = 0.1;   #  Allow up to 10% of the reads to be in a single partition.

    $pScale = 1.5       if (defined(getGlobal("homoPolyCompress")));

    #  If cnsPartitions set, switch to using the number of reads as the only
    #  partition sizing.  The +1 is to ensure that we don't end up with a
    #  tiny partition of left over reads.  pReads will be #.####.

    my $np = getGlobal("cnsPartitions");

    if ($np > 0) {
        $pSize   = 0.0;
        $pReads  = int(10000 / $np + 1) / 10000;
    }

    $cmd  = "$bin/utgcns \\\n";
    $cmd .= "  -S ../$asm.seqStore \\\n";
    $cmd .= "  -T  ./$asm.ctgStore 1 \\\n";
    $cmd .= "  -partition $pSize $pScale $pReads \\\n";
    $cmd .= "> ./$asm.ctgStore/partitioning.log 2>&1";

    if (runCommand("unitigging", $cmd)) {
        caExit("failed to partition the reads", "unitigging/$asm.ctgStore/partitioning.log");
    }

    stashFile("unitigging/$asm.ctgStore/partitioning.log");
    stashFile("unitigging/$asm.ctgStore/partitioning");
}



sub computeNumberOfConsensusJobs ($) {
    my $asm    = shift @_;
    my $jobs   = 0;
    my $bin    = getBinDirectory();

    fetchFile("unitigging/$asm.ctgStore/partitioning");

    open(F, "< unitigging/$asm.ctgStore/partitioning") or caExit("can't open 'unitigging/$asm.ctgStore/partitioning' for reading: $!", undef);
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

        my ($err, $all) = getAllowedResources("", "cns", "", "", 0);

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

    goto allDone   if (fileExists("unitigging/$asm.ctgStore/seqDB.v002.tig"));

    make_path($path)  if (! -d $path);

    #  Assign tigs to jobs.

    partitionTigs($asm);

    #  Figure out how many jobs.

    my $ctgjobs = computeNumberOfConsensusJobs($asm);

    #  This configure is an odd-ball.  Unlike all the other places that write scripts,
    #  we'll rewrite this one every time, so that we can change the alignment algorithm
    #  on the second attempt.

    my $firstTime = (! -e "$path/consensus.sh");

    if ((getGlobal("cnsConsensus") eq "quick") ||
        (getGlobal("cnsConsensus") eq "pbdagcon") ||
        (getGlobal("cnsConsensus") eq "utgcns")) {
        utgcns($asm, $ctgjobs);

    } else {
        caFailure("unknown consensus style '" . getGlobal("cnsConsensus") . "'", undef);
    }

    print STDERR "-- Configured $ctgjobs consensus jobs.\n";

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

    goto allDone  if (fileExists("$path/ctgcns.files"));
    goto allDone  if (fileExists("unitigging/$asm.ctgStore/seqDB.v002.tig"));

    fetchFile("$path/consensus.sh");

    #  Figure out if all the tasks finished correctly.

    my $jobs = computeNumberOfConsensusJobs($asm);

    #  Setup memory and threads and etc.  Complain if not enough memory.

    my $minMem = estimateMemoryNeededForConsensusJobs($asm);

    #  Decide what to run.

    caExit("no consensus jobs found?", undef)   if ($jobs == 0);

    my $currentJobID = "0001";

    my @successJobs;
    my @failedJobs;
    my $failureMessage = "";

    for (my $job=1; $job <= $jobs; $job++) {
        if      (fileExists("$path/ctgcns/$currentJobID.cns")) {
            push @successJobs, "5-consensus/ctgcns/$currentJobID.cns\n";

        } elsif (fileExists("$path/ctgcns/$currentJobID.cns.gz")) {
            push @successJobs, "5-consensus/ctgcns/$currentJobID.cns.gz\n";

        } elsif (fileExists("$path/ctgcns/$currentJobID.cns.bz2")) {
            push @successJobs, "5-consensus/ctgcns/$currentJobID.cns.bz2\n";

        } elsif (fileExists("$path/ctgcns/$currentJobID.cns.xz")) {
            push @successJobs, "5-consensus/ctgcns/$currentJobID.cns.xz\n";

        } else {
            $failureMessage .= "--   job ctgcns/$currentJobID.cns FAILED.\n";
            push @failedJobs, $job;
        }

        $currentJobID++;
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
    print STDERR "-- All ", scalar(@successJobs), " consensus jobs finished successfully.\n";

    open(L, "> $path/ctgcns.files") or caExit("can't open '$path/ctgcns.files' for writing: $!", undef);
    print L @successJobs;
    close(L);

    stashFile("$path/ctgcns.files");

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

    goto allDone    if (fileExists("unitigging/$asm.ctgStore/seqDB.v002.tig"));

    #  Expects to have a list of output files from the consensusCheck() function.

    fetchFile("$path/ctgcns.files");

    caExit("can't find '$path/ctgcns.files' for loading tigs into store: $!", undef)  if (! -e "$path/ctgcns.files");

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

    #  Remvoe consensus outputs

    if (-e "$path/ctgcns.files") {
        print STDERR "-- Purging consensus output after loading to ctgStore.\n";

        my $Ncns    = 0;
        my $Nfastq  = 0;
        my $Nlayout = 0;
        my $Nlog    = 0;

        ($Ncns, $Nfastq, $Nlayout, $Nlog) = purgeFiles($asm, "ctgcns", $Ncns, $Nfastq, $Nlayout, $Nlog);

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

