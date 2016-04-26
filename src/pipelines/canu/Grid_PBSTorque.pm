
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

package canu::Grid_PBSTorque;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(detectPBSTorque configurePBSTorque);

use strict;

use canu::Defaults;
use canu::Execution;
use canu::Grid;

sub detectPBSTorque () {

    return   if ( defined(getGlobal("gridEngine")));

    my $pbsnodes = findExecutable("pbsnodes");

    return   if (!defined($pbsnodes));

    print STDERR "-- Detected PBS/Torque with 'pbsnodes' binary in $pbsnodes.\n";
    setGlobal("gridEngine", "PBS");
}


sub configurePBSTorque () {

    return   if (uc(getGlobal("gridEngine")) ne "PBS");

    setGlobalIfUndef("gridEngineSubmitCommand",              "qsub");
    setGlobalIfUndef("gridEngineHoldOption",                 "-W depend=afteranyarray:WAIT_TAG");
    setGlobalIfUndef("gridEngineHoldOptionNoArray",          "-W depend=afterany:WAIT_TAG");
    setGlobalIfUndef("gridEngineSyncOption",                 "");
    setGlobalIfUndef("gridEngineNameOption",                 "-d `pwd` -N");
    setGlobalIfUndef("gridEngineArrayOption",                "-t ARRAY_JOBS");
    setGlobalIfUndef("gridEngineArrayName",                  "ARRAY_NAME");
    setGlobalIfUndef("gridEngineArrayMaxJobs",               65535);
    setGlobalIfUndef("gridEngineOutputOption",               "-j oe -o");
    setGlobalIfUndef("gridEngineThreadsOption",              "-l nodes=1:ppn=THREADS");
    setGlobalIfUndef("gridEngineMemoryOption",               "-l mem=MEMORY");
    setGlobalIfUndef("gridEnginePropagateCommand",           "qalter -W depend=afterany:\"WAIT_TAG\"");
    setGlobalIfUndef("gridEngineNameToJobIDCommand",         "qstat -f |grep -F -B 1 WAIT_TAG | grep Id: | grep -F [] |awk '{print \$NF}'");
    setGlobalIfUndef("gridEngineNameToJobIDCommandNoArray",  "qstat -f |grep -F -B 1 WAIT_TAG | grep Id: |awk '{print \$NF}'");
    setGlobalIfUndef("gridEngineTaskID",                     "PBS_ARRAYID");
    setGlobalIfUndef("gridEngineArraySubmitID",              "\\\$PBS_ARRAYID");
    setGlobalIfUndef("gridEngineJobID",                      "PBS_JOBID");

    #  Build a list of the resources available in the grid.  This will contain a list with keys
    #  of "#CPUs-#GBs" and values of the number of nodes With such a config.  Later on, we'll use this
    #  to figure out what specific settings to use for each algorithm.
    #
    #  The list is saved in global{"availableHosts"}
    #
    my %hosts;

    open(F, "pbsnodes |");

    while (<F>) {
        my $cpus = 0;
        my $mem = 0;
        if ($_ =~ m/status/) {
            my @stats = split ',', $_;
            for my $stat (@stats) {
                if ($stat =~ m/physmem/) {
                    $mem = ( split '=', $stat )[-1];
                } elsif ($stat =~ m/ncpus/) {
                    $cpus = int(( split '=', $stat )[-1]);
                }
            }
            $mem  = $1 * 1024         if ($mem =~ m/(\d+.*\d+)[tT]/);
            $mem  = $1 * 1            if ($mem =~ m/(\d+.*\d+)[gG]/);
            $mem  = $1 / 1024         if ($mem =~ m/(\d+.*\d+)[mM]/);
            $mem  = $1 / 1024 / 1024  if ($mem =~ m/(\d+.*\d+)[kK]/);
            $mem  = int($mem);
            $hosts{"$cpus-$mem"}++    if ($cpus gt 0);
        }
    }
    close(F);

    setGlobal("availableHosts", formatAllowedResources(%hosts, "PBS/Torque"));
}
