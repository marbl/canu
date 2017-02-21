
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


sub detectPBSVersion () {
    my $isPro   = 0;
    my $version = "";

    open(F, "pbsnodes --version 2>&1 |");
    while (<F>) {
        if (m/pbs_version\s+=\s+(.*)/) {
            $isPro   =  1;
            $version = $1;
        }
        if (m/Version:\s+(.*)/) {
            $version = $1;
        }
    }
    close(F);

    return($version, $isPro);
}


sub detectPBSTorque () {

    return   if ( defined(getGlobal("gridEngine")));

    my $pbsnodes = findExecutable("pbsnodes");

    return   if (!defined($pbsnodes));

    my ($version, $isPro) = detectPBSVersion();

    if ($isPro == 0) {
        print STDERR "-- Detected PBS/Torque '$version' with 'pbsnodes' binary in $pbsnodes.\n";
        setGlobal("gridEngine", "PBS");
    } else {
        print STDERR "-- Detected PBSPro '$version' with 'pbsnodes' binary in $pbsnodes.\n";
        setGlobal("gridEngine", "PBSPRO");
    }
}



sub configurePBSTorqueNodes () {
    my %hosts;

    print STDERR "-- Detecting PBS/Torque resources.\n";

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



sub configurePBSProNodes () {
    my %hosts;
    my $mem  = 0;
    my $cpus = 0;

    print STDERR "-- Detecting PBSPro resources.\n";

    open(F, "pbsnodes -av |");
    while (<F>) {
        if (m/resources_available.mem\s*=\s*(\d+)kb/) {
            $mem = int($1 / 1024 / 1024);
        }
        if (m/resources_available.mem\s*=\s*(\d+)mb/) {
            $mem = int($1 / 1024);
        }
        if (m/resources_available.mem\s*=\s*(\d+)gb/) {
            $mem = int($1);
        }

        if (m/resources_available.ncpus\s*=\s*(\d+)/) {
            $cpus = $1;
        }

        if (($cpus > 0) && ($mem > 0)) {
            $hosts{"$cpus-$mem"}++;
            $cpus = 0;
            $mem  = 0;
        }
    }
    close(F);

    setGlobal("availableHosts", formatAllowedResources(%hosts, "PBSPro"));
}



sub configurePBSTorque () {

    return   if ((uc(getGlobal("gridEngine")) ne "PBS") &&
                 (uc(getGlobal("gridEngine")) ne "PBSPRO"));

    my $isPro = (uc(getGlobal("gridEngine")) eq "PBSPRO");

    #  PBSPro, again, throws a curve ball at us.  There is no way to set the output of array jobs
    #  to someting reasonable like name.TASK_ID.err, even though that is basically the default.
    #  So, we unset gridEngineArraySubmitID to get the default name, but then need to move the '-j oe'
    #  somewhere else - and put it in the submit command.

    setGlobalIfUndef("gridEngineSubmitCommand",              "qsub -j oe -d `pwd`")   if ($isPro == 0);
    setGlobalIfUndef("gridEngineSubmitCommand",              "qsub -j oe")            if ($isPro == 1);
    setGlobalIfUndef("gridEngineNameOption",                 "-N");
    setGlobalIfUndef("gridEngineArrayOption",                "-t ARRAY_JOBS")  if ($isPro == 0);
    setGlobalIfUndef("gridEngineArrayOption",                "-J ARRAY_JOBS")  if ($isPro == 1);
    setGlobalIfUndef("gridEngineArrayName",                  "ARRAY_NAME");
    setGlobalIfUndef("gridEngineArrayMaxJobs",               268435456);  #  Effectively unlimited.
    setGlobalIfUndef("gridEngineOutputOption",               "-o");
    setGlobalIfUndef("gridEngineThreadsOption",              "-l nodes=1:ppn=THREADS");
    setGlobalIfUndef("gridEngineMemoryOption",               "-l mem=MEMORY");
    setGlobalIfUndef("gridEnginePropagateCommand",           "qalter -W depend=afterany:\"WAIT_TAG\"");
    setGlobalIfUndef("gridEngineNameToJobIDCommand",         "qstat -f |grep -F -B 1 WAIT_TAG | grep Id: | grep -F [] |awk '{print \$NF}'");
    setGlobalIfUndef("gridEngineNameToJobIDCommandNoArray",  "qstat -f |grep -F -B 1 WAIT_TAG | grep Id: |awk '{print \$NF}'");
    setGlobalIfUndef("gridEngineTaskID",                     "PBS_ARRAYID")           if ($isPro == 0);
    setGlobalIfUndef("gridEngineTaskID",                     "PBS_ARRAY_INDEX")       if ($isPro == 1);
    setGlobalIfUndef("gridEngineArraySubmitID",              "\\\$PBS_ARRAYID")       if ($isPro == 0);
    setGlobalIfUndef("gridEngineArraySubmitID",              undef)                   if ($isPro == 1);   #  Was "\\\$PBS_ARRAY_INDEX"
    setGlobalIfUndef("gridEngineJobID",                      "PBS_JOBID");

    #  Build a list of the resources available in the grid.  This will contain a list with keys
    #  of "#CPUs-#GBs" and values of the number of nodes With such a config.  Later on, we'll use this
    #  to figure out what specific settings to use for each algorithm.
    #
    #  The list is saved in global{"availableHosts"}

    configurePBSTorqueNodes()   if ($isPro == 0);
    configurePBSProNodes()      if ($isPro == 1);
}
