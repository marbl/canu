
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

package canu::Grid_Slurm;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(detectSlurm configureSlurm);

use strict;
use warnings "all";
no  warnings "uninitialized";

use canu::Defaults;
use canu::Execution;

use canu::Grid "formatAllowedResources";



sub detectSlurm () {

    return   if ( defined(getGlobal("gridEngine")));

    my $sinfo = findExecutable("sinfo");

    return   if (!defined($sinfo));

    if (getGlobal("useGrid") eq "0") {
        print STDERR "--\n";
        print STDERR "-- Detected Slurm with 'sinfo' binary in $sinfo.\n";
        print STDERR "--          Slurm disabled by useGrid=false\n";
    }
    else {
        print STDERR "--\n";
        print STDERR "-- Detected Slurm with 'sinfo' binary in $sinfo.\n";

        setGlobal("gridEngine", "SLURM");
    }
}


sub configureSlurm () {

    return   if (uc(getGlobal("gridEngine")) ne "SLURM");

    my $maxArraySize  = 65535;
    my $maxArrayTasks = 65535;

    #  Slurm has one or two configuration variables we need to check.
    #    MaxArraySize    is (one more than) the maximum index value of any array task.
    #    max_array_tasks is the maximum number of array tasks per array job.
    #
    #  Usually they're the same value.  We just need to pick the smaller of the two.

    open(F, "scontrol show config |") or caExit("can't run 'scontrol' to get SLURM config", undef);
    while (<F>) {
        if (m/MaxArraySize\s*=\s*(\d+)/) {
            $maxArraySize  = $1 - 1;
            print STDERR "-- Detected Slurm with task IDs up to $maxArraySize allowed.\n";
        }
        if (m/max_array_tasks\s*=\s*(\d+)/) {
            $maxArrayTasks = $1 - 1;
            print STDERR "-- Detected Slurm with up to $maxArrayTasks tasks per job.\n";
        }
    }
    close(F);

    if ($maxArraySize > $maxArrayTasks) {   #  Now just pick the smaller.
        $maxArraySize = $maxArrayTasks;
    }

    setGlobalIfUndef("gridEngineSubmitCommand",              "sbatch");
    setGlobalIfUndef("gridEngineNameOption",                 "-D `pwd` -J");
    setGlobalIfUndef("gridEngineArrayOption",                "-a ARRAY_JOBS");
    setGlobalIfUndef("gridEngineArrayName",                  "ARRAY_NAME");
    setGlobalIfUndef("gridEngineArrayMaxJobs",               $maxArraySize);
    setGlobalIfUndef("gridEngineOutputOption",               "-o");                                        ## NB: SLURM default joins STDERR & STDOUT if no -e specified
    setGlobalIfUndef("gridEngineResourceOption",             "--cpus-per-task=THREADS --mem-per-cpu=MEMORY");
    #setGlobalIfUndef("gridEngineMemoryPerJob",              "0");   #  Do NOT set; default set below
    setGlobalIfUndef("gridEngineNameToJobIDCommand",         "squeue -h -o\%F -n \"WAIT_TAG\" | uniq");    ## TODO: manually verify this in all cases
    setGlobalIfUndef("gridEngineNameToJobIDCommandNoArray",  "squeue -h -o\%i -n \"WAIT_TAG\"");     ## TODO: manually verify this in all cases
    setGlobalIfUndef("gridEngineTaskID",                     "SLURM_ARRAY_TASK_ID");
    setGlobalIfUndef("gridEngineArraySubmitID",              "%A_%a");
    setGlobalIfUndef("gridEngineJobID",                      "SLURM_JOB_ID");

    #  Set memory-per-job if configured and the user hasn't supplied a value.
    #  Otherwise, leave it in the default off state.
    if (!defined(getGlobal("gridEngineMemoryPerJob"))) {
        if (getGlobal("gridEngineResourceOption") !~ m/mem-per-cpu/i) {
            print STDERR "-- Enable memory-per-job mode.\n";
            setGlobalfUnset("gridEngineMemoryPerJob", "1");
        }
    }


    #  Build a list of the resources available in the grid.  This will contain a list with keys
    #  of "#CPUs-#GBs" and values of the number of nodes With such a config.  Later on, we'll use this
    #  to figure out what specific settings to use for each algorithm.
    #
    #  The list is saved in global{"availableHosts"}
    #
    my %hosts;

    #NODELIST NODES CPUS MEMORY
    open(F, "sinfo --exact -o '%N %D %c %m' | grep -v drained | grep -v interactive |");
    my $h = <F>;  #  header

    my @h = split '\s+', $h;

    my $nodeIdx = 1;
    my $cpuIdx  = 4;
    my $memIdx  = 6;

    for (my $ii=0; ($ii < scalar(@h)); $ii++) {
        $nodeIdx = $ii  if ($h[$ii] eq "NODES");
        $cpuIdx  = $ii  if ($h[$ii] eq "CPUS");
        $memIdx  = $ii  if ($h[$ii] eq "MEMORY");
    }

    while (<F>) {
        my @v = split '\s+', $_;

        my $cpus  = int($v[$cpuIdx]);
        my $mem   = int($v[$memIdx] / 1024 - 0.5);   #  Round down.
        my $nodes = int($v[$nodeIdx]);

        $hosts{"$cpus-$mem"}+=int($nodes)    if ($cpus gt 0);
    }
    close(F);
    setGlobal("availableHosts", formatAllowedResources(%hosts, "Slurm"));
}
