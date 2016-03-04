
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
 #    Brian P. Walenz from 2015-MAR-06 to 2015-AUG-25
 #      are Copyright 2015 Battelle National Biodefense Institute, and
 #      are subject to the BSD 3-Clause License
 #
 #    Brian P. Walenz beginning on 2015-NOV-03
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #    Sergey Koren beginning on 2015-DEC-16
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

package canu::Consensus;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(consensusConfigure consensusCheck consensusLoad consensusAnalyze);

use strict;

use File::Path qw(make_path remove_tree);

use canu::Defaults;
use canu::Execution;
use canu::Gatekeeper;
use canu::Unitig;
use canu::HTML;


sub computeNumberOfConsensusJobs ($$) {
    my $wrk    = shift @_;  #  Local work directory
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
    my $wrk    = shift @_;  #  Local work directory
    my $asm    = shift @_;
    my $jobs   = shift @_;

    #getAllowedResources("", "cns");

    open(F, "> $wrk/5-consensus/consensus.sh") or caExit("can't open '$wrk/5-consensus/consensus.sh' for writing: $!", undef);

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
    #print F "  -L $wrk/5-consensus/\$jobid.layout.WORKING \\\n";
    #print F "  -Q $wrk/5-consensus/\$jobid.fastq.WORKING \\\n";
    print F "  -maxcoverage " . getGlobal('cnsMaxCoverage') . " \\\n";
    print F "  -e " . getGlobal("cnsErrorRate") . " \\\n";
    print F "  -quick \\\n"      if (getGlobal("cnsConsensus") eq "quick");
    print F "  -pbdagcon \\\n"   if (getGlobal("cnsConsensus") eq "pbdagcon");
    print F "  -utgcns \\\n"     if (getGlobal("cnsConsensus") eq "utgcns");
    print F "  -threads " . getGlobal("cnsThreads") . " \\\n";
    print F "&& \\\n";
    print F "mv $wrk/5-consensus/\$jobid.cns.WORKING $wrk/5-consensus/\$jobid.cns \\\n";
    #print F "&& \\\n";
    #print F "mv $wrk/5-consensus/\$jobid.layout.WORKING $wrk/5-consensus/\$jobid.layout \\\n";
    #print F "&& \\\n";
    #print F "mv $wrk/5-consensus/\$jobid.fastq.WORKING $wrk/5-consensus/\$jobid.fastq\n";
    print F "\n";
    print F "exit 0\n";

    close(F);
}



sub consensusConfigure ($$) {
    my $WRK    = shift @_;           #  Root work directory
    my $wrk    = "$WRK/unitigging";  #  Local work directory
    my $asm    = shift @_;
    my $bin    = getBinDirectory();
    my $cmd;
    my $path   = "$wrk/5-consensus";

    goto allDone   if (skipStage($WRK, $asm, "consensusConfigure") == 1);
    goto allDone   if (-e "$wrk/$asm.tigStore/seqDB.v002.tig");

    make_path("$path")  if (! -d "$path");

    #  If the gkpStore partitions are older than the tigStore unitig output, assume the unitigs have
    #  changed and remove the gkpStore partition.  -M is (annoyingly) 'file age', so we need to
    #  rebuild if gkp is older (larger) than tig.

    if (-e "$wrk/$asm.gkpStore/partitions/map") {
        my $gkpTime = -M "$wrk/$asm.gkpStore/partitions/map";
        my $tigTime = -M "$wrk/$asm.tigStore/seqDB.v001.tig";

        if ($gkpTime > $tigTime) {
            print STDERR "-- Partitioned gkpStore is older than tigs, rebuild partitioning (gkpStore $gkpTime days old; tigStore $tigTime days old).\n";

            if (runCommandSilently($wrk, "rm -rf $wrk/$asm.gkpStore/partitions", 1)) {
                caExit("failed to remove old partitions ($wrk/$asm.gkpStore/partitions), can't continue until these are removed", undef);
            }
        }
    }

    #  Partition gkpStore if needed.

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

    #  Set up the consensus compute.  It's in a useless if chain because there used to be
    #  different executables; now they're all rolled into utgcns itself.

    my $jobs = computeNumberOfConsensusJobs($wrk, $asm);

    if ((getGlobal("cnsConsensus") eq "quick") ||
        (getGlobal("cnsConsensus") eq "pbdagcon") ||
        (getGlobal("cnsConsensus") eq "utgcns")) {
        utgcns($wrk, $asm, $jobs);

    } else {
        caFailure("unknown consensus style '" . getGlobal("cnsConsensus") . "'", undef);
    }

  finishStage:
    emitStage($WRK, $asm, "consensusConfigure");
    buildHTML($WRK, $asm, "utg");
    stopAfter("consensusConfigure");

  allDone:
    print STDERR "-- Configured ", computeNumberOfConsensusJobs($wrk, $asm), " consensus jobs.\n";
}





#  Checks that all consensus jobs are complete, loads them into the store.
#
sub consensusCheck ($$) {
    my $WRK     = shift @_;           #  Root work directory
    my $wrk     = "$WRK/unitigging";  #  Local work directory
    my $asm     = shift @_;
    my $attempt = getGlobal("canuIteration");
    my $path    = "$wrk/5-consensus";

    goto allDone  if (skipStage($WRK, $asm, "consensusCheck", $attempt) == 1);
    goto allDone  if (-e "$path/cnsjob.files");
    goto allDone  if (-e "$wrk/$asm.tigStore/seqDB.v002.tig");

    #  Figure out if all the tasks finished correctly.

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
            $failureMessage .= "--   job $path/$currentJobID.cns FAILED.\n";
            push @failedJobs, $job;
        }

        $currentJobID++;
    }

    #  Failed jobs, retry.

    if (scalar(@failedJobs) > 0) {

        #  If not the first attempt, report the jobs that failed, and that we're recomputing.

        if ($attempt > 1) {
            print STDERR "--\n";
            print STDERR "-- ", scalar(@failedJobs), " consensus jobs failed:\n";
            print STDERR $failureMessage;
            print STDERR "--\n";
        }

        #  If too many attempts, give up.

        if ($attempt > getGlobal("canuIterationMax")) {
            caExit("failed to generate consensus.  Made " . ($attempt-1) . " attempts, jobs still failed", undef);
        }

        #  Otherwise, run some jobs.

        print STDERR "-- Consensus attempt $attempt begins with ", scalar(@successJobs), " finished, and ", scalar(@failedJobs), " to compute.\n";

        emitStage($WRK, $asm, "consensusCheck", $attempt);
        buildHTML($WRK, $asm, "utg");

        submitOrRunParallelJob($WRK, $asm, "cns", $path, "consensus", @failedJobs);
        return;
    }

  finishStage:
    print STDERR "-- All ", scalar(@successJobs), " consensus jobs finished successfully.\n";

    open(L, "> $path/cnsjob.files") or caExit("can't open '$path/cnsjob.files' for writing: $!", undef);
    print L @successJobs;
    close(L);

    setGlobal("canuIteration", 0);
    emitStage($WRK, $asm, "consensusCheck");
    buildHTML($WRK, $asm, "utg");
    stopAfter("consensusCheck");

  allDone:
}




sub consensusLoad ($$) {
    my $WRK     = shift @_;           #  Root work directory
    my $wrk     = "$WRK/unitigging";  #  Local work directory
    my $asm     = shift @_;
    my $bin     = getBinDirectory();
    my $cmd;
    my $path    = "$wrk/5-consensus";

    goto allDone    if (skipStage($WRK, $asm, "consensusLoad") == 1);
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
        print STDERR "-- Purging consensus output after loading to tigStore.\n";

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

        print STDERR "-- Purged $Ncns .cns outputs.\n"        if ($Ncns > 0);
        print STDERR "-- Purged $Nfastq .fastq outputs.\n"    if ($Nfastq > 0);
        print STDERR "-- Purged $Nlayout .layout outputs.\n"  if ($Nlayout > 0);
        print STDERR "-- Purged $Nlog .err log outputs.\n"    if ($Nlog > 0);
    }

  finishStage:
    emitStage($WRK, $asm, "consensusLoad");
    buildHTML($WRK, $asm, "utg");
    stopAfter("consensusLoad");
  allDone:
    reportUnitigSizes($wrk, $asm, 2, "after consenss generation");
}




sub consensusAnalyze ($$) {
    my $WRK     = shift @_;           #  Root work directory
    my $wrk     = "$WRK/unitigging";  #  Local work directory
    my $asm     = shift @_;
    my $bin     = getBinDirectory();
    my $cmd;
    my $path    = "$wrk/5-consensus";

    goto allDone   if (skipStage($WRK, $asm, "consensusAnalyze") == 1);
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
    emitStage($WRK, $asm, "consensusAnalyze");
    buildHTML($WRK, $asm, "utg");
    touch("$wrk/$asm.tigStore/status.coverageStat");
    stopAfter("consensusAnalyze");
  allDone:
}
