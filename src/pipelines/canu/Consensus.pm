
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
 #    src/pipelines/ca3g/Consensus.pm
 #
 #  Modifications by:
 #
 #    Brian P. Walenz beginning on 2015-MAR-06
 #      are Copyright 2015 Battelle National Biodefense Institute, and
 #      are subject to the BSD 3-Clause License
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

package canu::Consensus;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(consensusConfigure consensusCheck consensusLoad consensusAnalyze consensusFilter);

use strict;

use File::Path qw(make_path remove_tree);

use canu::Defaults;
use canu::Execution;
use canu::Gatekeeper;
use canu::Unitig;



sub computeNumberOfConsensusJobs ($$) {
    my $wrk    = shift @_;
    my $asm    = shift @_;

    my $jobs = 0;
    open(F, "< $wrk/4-unitigger/$asm.partitioningInfo") or caExit("can't open '$wrk/4-unitigger/$asm.partitioningInfo' for reading: $!", undef);
    while (<F>) {
        if (m/Partition\s+(\d+)\s+has\s+(\d+)\s+unitigs\sand\s+(\d+)\s+fragments./) {
            $jobs = $1;
        }
    }
    close(F);

    return($jobs);
}



sub utgcns ($$$) {
    my $wrk    = shift @_;
    my $asm    = shift @_;
    my $jobs   = shift @_;

    my $taskID            = getGlobal("gridEngineTaskID");
    my $submitTaskID      = getGlobal("gridEngineArraySubmitID");

    open(F, "> $wrk/5-consensus/consensus.sh") or caExit("can't open '$wrk/5-consensus/consensus.sh' for writing: $!", undef);

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
    print F "jobid=`printf %04d \$jobid`\n";
    print F "\n";
    print F "if [ -e $wrk/5-consensus/\$jobid.cns ] ; then\n";
    print F "  exit 0\n";
    print F "fi\n";
    print F "\n";
    print F getBinDirectoryShellCode();
    print F "\n";
    print F "\$bin/utgcns \\\n";
    print F "  -G $wrk/$asm.gkpStore \\\n";
    print F "  -T $wrk/$asm.tigStore 1 \$jobid \\\n";
    print F "  -O $wrk/5-consensus/\$jobid.cns.WORKING \\\n";
    print F "  -L $wrk/5-consensus/\$jobid.layout.WORKING \\\n";
    print F "  -F $wrk/5-consensus/\$jobid.fastq.WORKING \\\n";
    print F "  -maxcoverage " . getGlobal('cnsMaxCoverage') . " \\\n";
    print F "&& \\\n";
    print F "mv $wrk/5-consensus/\$jobid.cns.WORKING $wrk/5-consensus/\$jobid.cns \\\n";
    print F "&& \\\n";
    print F "mv $wrk/5-consensus/\$jobid.layout.WORKING $wrk/5-consensus/\$jobid.layout \\\n";
    print F "&& \\\n";
    print F "mv $wrk/5-consensus/\$jobid.fastq.WORKING $wrk/5-consensus/\$jobid.fastq\n";
    print F "\n";
    print F "exit 0\n";

    close(F);
}



sub make_consensus ($$$) {
    my $wrk    = shift @_;
    my $asm    = shift @_;
    my $jobs   = shift @_;

    die "No make_consensus yet.\n";
}



sub falcon_sense ($$$) {
    my $wrk    = shift @_;
    my $asm    = shift @_;
    my $jobs   = shift @_;

    die "No falcon_sense yet.\n";
}



sub pbdagcon ($$$) {
    my $wrk    = shift @_;
    my $asm    = shift @_;
    my $jobs   = shift @_;

    my $blasr = "";

    my $taskID            = getGlobal("gridEngineTaskID");
    my $submitTaskID      = getGlobal("gridEngineArraySubmitID");

    open(F, "> $wrk/5-consensus/consensus.sh") or caExit("can't open '$wrk/5-consensus/consensus.sh' for writing: $!", undef);

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
    print F "jobid=`printf %04d \$jobid`\n";
    print F "\n";
    print F "if [ -e $wrk/5-consensus/\$jobid.cns ] ; then\n";
    print F "  exit 0\n";
    print F "fi\n";
    print F "\n";
    print F getBinDirectoryShellCode();
    print F "\n";
    print F "cat $wrk/4-unitigger/$asm.partitioning | awk -v JOB=\$jobid '{if (\$1 == JOB) print \$NF}' > $wrk/5-consensus/$asm.\$jobid.iid\n";
    print F "\n";
    print F "\$bin/gatekeeper -dumpfasta $wrk/5-consensus/$asm.\$jobid -iid $wrk/5-consensus/$asm.\$jobid.iid $wrk/$asm.gkpStore\n";
    print F "\n";
    print F "rm -f $wrk/5-consensus/$asm.\$jobid.q*\n";
    print F "\n";
    print F "\$bin/tigStore -d layout -U -t $wrk/$asm.tigStore 1 -up \$jobid -g $wrk/$asm.gkpStore > $wrk/5-consensus/$asm.\$jobid.lay\n";
    print F "\n";
    print F "\$bin/convertToPBCNS -path $blasr -consensus pbdagcon -coverage 1 -threads " . getGlobal("cnsConcurrency") . " -prefix $wrk/5-consensus/$asm.\$jobid.tmp -length 500 -sequence $wrk/5-consensus/$asm.\$jobid.fasta -input $wrk/5-consensus/$asm.\$jobid.lay -output $wrk/5-consensus/$asm.\$jobid.fa\n";
    print F "\n";
    print F "\$bin/addCNSToStore -path \$bin -input $wrk/5-consensus/$asm.\$jobid.fa -lay $wrk/5-consensus/$asm.\$jobid.lay -output $wrk/5-consensus/$asm.\$jobid.cns -prefix $wrk/$asm -sequence $wrk/5-consensus/$asm.\$jobid.fasta -partition \$jobid && \$bin/utgcnsfix -g $wrk/$asm.gkpStore  -t $wrk/$asm.tigStore 2 \$jobid -o $wrk/5-consensus/${asm}_\$jobid.fixes > $wrk/5-consensus/${asm}_\$jobid.fix.err 2>&1 && touch $wrk/5-consensus/${asm}_\$jobid.cns\n";
    print F "\n";
    print F "if [ -e $wrk/5-consensus/${asm}_\$jobid.cns ]; then\n";
    print F "   rm -f $wrk/5-consensus/${asm}.\$jobid.fasta*\n";
    print F "   rm -f $wrk/5-consensus/${asm}.\$jobid.lay\n";
    print F "fi\n";

    close(F);

    setGlobal("utgConcurrency", 1);

}





sub consensusConfigure ($$) {
    my $wrk     = shift @_;
    my $asm     = shift @_;
    my $bin     = getBinDirectory();
    my $cmd;
    my $path    = "$wrk/5-consensus";

    goto allDone   if (skipStage($wrk, $asm, "consensusConfigure") == 1);
    goto allDone   if (-e "$wrk/$asm.tigStore/seqDB.v002.tig");

    make_path("$path")  if (! -d "$path");

    #  Make sure the gkpStore is partitioned.

    if (! -e "$wrk/$asm.gkpStore/partitions/map") {
        $cmd  = "$bin/gatekeeperPartition \\\n";
        $cmd .= "  -G $wrk/$asm.gkpStore \\\n";
        $cmd .= "  -P $wrk/4-unitigger/$asm.partitioning \\\n";
        $cmd .= "> $path/$asm.partitioned.err 2>&1";

        stopBefore("consensusConfigure", $cmd);

        if (runCommand("$path", $cmd)) {
            caExit("failed to partition the reads", "$path/$asm.partitioned.err");
        }
    }

    my $jobs = computeNumberOfConsensusJobs($wrk, $asm);

    if      (getGlobal("cnsConsensus") eq "utgcns") {
        utgcns($wrk, $asm, $jobs);

    } elsif (getGlobal("cnsConsensus") eq "make_consensus") {
        make_consensus($wrk, $asm, $jobs);

    } elsif (getGlobal("cnsConsensus") eq "falcon_sense") {
        falcon_sense($wrk, $asm, $jobs);

    } elsif (getGlobal("cnsConsensus") eq "pbdagcon") {
        pbdagcon($wrk, $asm, $jobs);

    } else {
        caFailure("unknown consensus style '" . getGlobal("cnsConsensus") . "'", undef);
    }

  finishStage:
    emitStage($wrk, $asm, "consensusConfigure");
    stopAfter("consensusConfigure");

  allDone:
    print STDERR "--  Configured ", computeNumberOfConsensusJobs($wrk, $asm), " consensus jobs.\n";
}





#  Checks that all consensus jobs are complete, loads them into the store.
#
sub consensusCheck ($$$) {
    my $wrk     = shift @_;
    my $asm     = shift @_;
    my $attempt = shift @_;
    my $path    = "$wrk/5-consensus";

    goto allDone  if (skipStage($wrk, $asm, "consensusCheck", $attempt) == 1);
    goto allDone  if (-e "$path/cnsjob.files");
    goto allDone  if (-e "$wrk/$asm.tigStore/seqDB.v002.tig");

    my $jobs = computeNumberOfConsensusJobs($wrk, $asm);

    my $currentJobID = "0001";
    my @successJobs;
    my @failedJobs;
    my $failureMessage = "";

    for (my $job=1; $job <= $jobs; $job++) {
        if      (-e "$path/$currentJobID.cns") {
            push @successJobs, "$path/$currentJobID.cns\n";

        } elsif (-e "$path/$currentJobID.cns.gz") {
            push @successJobs, "$path/$currentJobID.cns.gz\n";

        } elsif (-e "$path/$currentJobID.cns.bz2") {
            push @successJobs, "$path/$currentJobID.cns.bz2\n";

        } elsif (-e "$path/$currentJobID.cns.xz") {
            push @successJobs, "$path/$currentJobID.cns.xz\n";

        } else {
            $failureMessage .= "--    job $path/$currentJobID.cns FAILED.\n";
            push @failedJobs, $job;
        }

        $currentJobID++;
    }

    #  No failed jobs?  Success!

    if (scalar(@failedJobs) == 0) {
        open(L, "> $path/cnsjob.files") or caExit("can't open '$path/cnsjob.files' for writing: $!", undef);
        print L @successJobs;
        close(L);
        setGlobal("canuIteration", 0);
        emitStage($wrk, $asm, "consensusCheck");
        return;
    }

    #  If not the first attempt, report the jobs that failed, and that we're recomputing.

    if ($attempt > 1) {
        print STDERR "--\n";
        print STDERR "--  ", scalar(@failedJobs), " consensus jobs failed:\n";
        print STDERR $failureMessage;
        print STDERR "--\n";
    }

    #  If too many attempts, give up.

    if ($attempt > 2) {
        caExit("failed to generate consensus.  Made " . ($attempt-1) . " attempts, jobs still failed", undef);
    }

    #  Otherwise, run some jobs.

    print STDERR "--  Consensus attempt $attempt begins with ", scalar(@successJobs), " finished, and ", scalar(@failedJobs), " to compute.\n";

    emitStage($wrk, $asm, "consensusCheck", $attempt);

    submitOrRunParallelJob($wrk, $asm, "cns", $path, "consensus", @failedJobs);

  finishStage:
    emitStage($wrk, $asm, "consensusCheck");
    stopAfter("consensusCheck");
  allDone:
    print STDERR "--  All ", computeNumberOfConsensusJobs($wrk, $asm), " consensus jobs finished successfully.\n";
}




sub consensusLoad ($$) {
    my $wrk     = shift @_;
    my $asm     = shift @_;
    my $bin     = getBinDirectory();
    my $cmd;
    my $path    = "$wrk/5-consensus";

    goto allDone    if (skipStage($wrk, $asm, "consensusLoad") == 1);
    goto allDone    if (-e "$wrk/$asm.tigStore/seqDB.v002.tig");

    #  Expects to have a cnsjob.files list of output files from the consensusCheck() function.

    caExit("can't find '$path/cnsjob.files' for loading tigs into store: $!", undef)  if (! -e "$path/cnsjob.files");

    #  Now just load them.

    $cmd  = "$bin/tgStoreLoad \\\n";
    $cmd .= "  -G $wrk/$asm.gkpStore \\\n";
    $cmd .= "  -T $wrk/$asm.tigStore 2 \\\n";
    $cmd .= "  -L $path/cnsjob.files \\\n";
    $cmd .= "> $path/cnsjobs.files.tigStoreLoad.err 2>&1";

    if (runCommand($path, $cmd)) {
        caExit("failed to load unitig consensus into tigStore", "$path/cnsjobs.files.tigStoreLoad.err");
    }

    #  Remvoe consensus outputs

    if (-e "$path/cnsjob.files") {
        print STDERR "--  Purging consensus output after loading to tigStore.\n";

        my $Ncns    = 0;
        my $Nfastq  = 0;
        my $Nlayout = 0;
        my $Nlog    = 0;

        open(F, "< $path/cnsjob.files") or caExit("can't open '$path/cnsjob.files' for reading: $!\n", undef);
        while (<F>) {
            chomp;
            if (m/^(.*)\/0*(\d+).cns$/) {
                my $ID6 = substr("00000" . $2, -6);
                my $ID4 = substr("000"   . $2, -4);
                my $ID0 = $2;

                if (-e "$1/$ID4.cns") {
                    $Ncns++;
                    unlink "$1/$ID4.cns";
                }
                if (-e "$1/$ID4.fastq") {
                    $Nfastq++;
                    unlink "$1/$ID4.fastq";
                }
                if (-e "$1/$ID4.layout") {
                    $Nlayout++;
                    unlink "$1/$ID4.layout";
                }
                if (-e "$1/consensus.$ID6.out") {
                    $Nlog++;
                    unlink "$1/consensus.$ID6.out";
                }
                if (-e "$1/consensus.$ID0.out") {
                    $Nlog++;
                    unlink "$1/consensus.$ID0.out";
                }

            } else {
                caExit("unknown consensus job name '$_'\n", undef);
            }
        }
        close(F);

        print STDERR "--  Purged $Ncns .cns outputs.\n"        if ($Ncns > 0);
        print STDERR "--  Purged $Nfastq .fastq outputs.\n"    if ($Nfastq > 0);
        print STDERR "--  Purged $Nlayout .layout outputs.\n"  if ($Nlayout > 0);
        print STDERR "--  Purged $Nlog .err log outputs.\n"    if ($Nlog > 0);
    }

  finishStage:
    emitStage($wrk, $asm, "consensusLoad");
    stopAfter("consensusLoad");
  allDone:
    reportUnitigSizes($wrk, $asm, 2, "after consenss generation");
}




sub consensusAnalyze ($$) {
    my $wrk     = shift @_;
    my $asm     = shift @_;
    my $bin     = getBinDirectory();
    my $cmd;
    my $path    = "$wrk/5-consensus";

    goto allDone   if (skipStage($wrk, $asm, "consensusAnalyze") == 1);
    goto allDone   if (-e "$wrk/$asm.tigStore/status.coverageStat");

    $cmd  = "$bin/tgStoreCoverageStat \\\n";
    $cmd .= "  -G       $wrk/$asm.gkpStore \\\n";
    $cmd .= "  -T       $wrk/$asm.tigStore 2 \\\n";
    $cmd .= "  -s       " . getGlobal("genomeSize") . " \\\n";
    $cmd .= "  -o       $wrk/$asm.tigStore.coverageStat \\\n";
    $cmd .= "> $wrk/$asm.tigStore.coverageStat.err 2>&1";

    if (runCommand($path, $cmd)) {
        caExit("failed to compute coverage statistics", "$wrk/$asm.tigStore.coverageStat.err");
    }

    unlink "$wrk/$asm.tigStore.coverageStat.err";

  finishStage:
    emitStage($wrk, $asm, "consensusAnalyze");
    touch("$wrk/$asm.tigStore/status.coverageStat");
    stopAfter("consensusAnalyze");
  allDone:
}


sub consensusFilter ($$) {
    my $wrk     = shift @_;
    my $asm     = shift @_;
    my $bin     = getBinDirectory();
    my $cmd;
    my $path    = "$wrk/5-consensus";

    goto allDone   if (skipStage($wrk, $asm, "consensusFilter") == 1);
    goto allDone   if (-e "$wrk/$asm.tigStore/status.filter");

    my $msrs = getGlobal("maxSingleReadSpan");
    my $lca  = getGlobal("lowCoverageAllowed");  #  checkParameters() ensures that this and
    my $lcd  = getGlobal("lowCoverageDepth");    #  this are both set or both unset
    my $mru  = getGlobal("minReadsUnique");
    my $mul  = getGlobal("minUniqueLength");
    my $mrl  = getGlobal("maxRepeatLength");

    $cmd  = "$bin/tgStoreFilter \\\n";
    $cmd .= "  -G       $wrk/$asm.gkpStore \\\n";
    $cmd .= "  -T       $wrk/$asm.tigStore 2 \\\n";
    $cmd .= "  -o       $wrk/$asm.tigStore.filter \\\n";
    $cmd .= "  -span    $msrs \\\n"       if (defined($msrs));;
    $cmd .= "  -lowcov  $lca $lcd \\\n"   if (defined($lca) && defined($lcd));;
    $cmd .= "  -reads   $mru \\\n"        if (defined($mru));
    $cmd .= "  -long    $mul \\\n"        if (defined($mul));
    $cmd .= "  -short   $mrl \\\n"        if (defined($mrl));
    $cmd .= "> $wrk/$asm.tigStore.filter.err 2>&1";

    if (runCommand($path, $cmd)) {
        caExit("failed to filter unitigs", "$wrk/$asm.tigStore.filter.err");
    }

    unlink "$wrk/$asm.tigStore.filter.err";

  finishStage:
    emitStage($wrk, $asm, "consensusFilter");
    touch("$wrk/$asm.tigStore/status.filter");
    stopAfter("consensusFilter");
  allDone:
}
