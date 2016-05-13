
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
 #  Modifications by:
 #
 #    Brian P. Walenz beginning on 2015-NOV-27
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #    Sergey Koren beginning on 2015-NOV-30
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

package canu::Grid_Slurm;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(detectSlurm configureSlurm);

use strict;

use canu::Defaults;
use canu::Execution;
use canu::Grid;

sub detectSlurm () {

    return   if ( defined(getGlobal("gridEngine")));

    my $sinfo = findExecutable("sinfo");

    return   if (!defined($sinfo));

    print STDERR "-- Detected Slurm with 'sinfo' binary in $sinfo.\n";
    setGlobal("gridEngine", "SLURM");
}


sub configureSlurm () {

    return   if (uc(getGlobal("gridEngine")) ne "SLURM");

    my $maxArraySize = 65535;

    open(F, "scontrol show config |") or caExit("can't run 'scontrol' to get SLURM config", undef);
    while (<F>) {
        if (m/MaxArraySize\s+=\s+(\d+)/) {
            $maxArraySize = $1;
        }
    }
    close(F);

    setGlobalIfUndef("gridEngineSubmitCommand",              "sbatch");
    setGlobalIfUndef("gridEngineHoldOption",                 "--depend=afterany:WAIT_TAG");
    setGlobalIfUndef("gridEngineHoldOptionNoArray",          "--depend=afterany:WAIT_TAG");
    setGlobalIfUndef("gridEngineSyncOption",                 "");                                          ## TODO: SLURM may not support w/out wrapper; See LSF bsub manpage to compare
    setGlobalIfUndef("gridEngineNameOption",                 "-D `pwd` -J");
    setGlobalIfUndef("gridEngineArrayOption",                "-a ARRAY_JOBS");
    setGlobalIfUndef("gridEngineArrayName",                  "ARRAY_NAME");
    setGlobalIfUndef("gridEngineArrayMaxJobs",               $maxArraySize);
    setGlobalIfUndef("gridEngineOutputOption",               "-o");                                        ## NB: SLURM default joins STDERR & STDOUT if no -e specified
    setGlobalIfUndef("gridEngineThreadsOption",              "--cpus-per-task=THREADS");
    setGlobalIfUndef("gridEngineMemoryOption",               "--mem-per-cpu=MEMORY");
    setGlobalIfUndef("gridEnginePropagateCommand",           "scontrol update job=\"WAIT_TAG\"");          ## TODO: manually verify this in all cases
    setGlobalIfUndef("gridEngineNameToJobIDCommand",         "squeue -h -o\%F -n \"WAIT_TAG\" | uniq");    ## TODO: manually verify this in all cases
    setGlobalIfUndef("gridEngineNameToJobIDCommandNoArray",  "squeue -h -o\%i -n \"WAIT_TAG\"");     ## TODO: manually verify this in all cases
    setGlobalIfUndef("gridEngineTaskID",                     "SLURM_ARRAY_TASK_ID");
    setGlobalIfUndef("gridEngineArraySubmitID",              "%A_%a");
    setGlobalIfUndef("gridEngineJobID",                      "SLURM_JOB_ID");


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

        my $cpus  = $v[$cpuIdx];
        my $mem   = $v[$memIdx] / 1024;
        my $nodes = $v[$nodeIdx];

        $hosts{"$cpus-$mem"}+=int($nodes)    if ($cpus gt 0);
    }
    close(F);
    setGlobal("availableHosts", formatAllowedResources(%hosts, "Slurm"));
}
