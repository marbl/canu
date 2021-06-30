
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

package canu::Grid_Local;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(detectLocal
             configureLocal);

use strict;
use warnings "all";
no  warnings "uninitialized";

use List::Util qw(min max);

use canu::Defaults;
use canu::Execution;

use canu::Grid "formatAllowedResources";



#  The following two get() functions should only be used in this file.
#  Everyone else should use getGlobal(localMemory) and
#  getGlobal(localThreads).

sub getNumberOfCPUs () {
    my $os   = $^O;
    my $ncpu = 1;

    #  See http://stackoverflow.com/questions/6481005/obtain-the-number-of-cpus-cores-in-linux

    if ($os eq "freebsd") {
        $ncpu = int(`/sbin/sysctl -n hw.ncpu`);
    }

    if ($os eq "darwin") {
        $ncpu = int(`/usr/bin/getconf _NPROCESSORS_ONLN`);
    }

    if ($os eq "linux" || $os eq "cygwin") {
        $ncpu = int(`getconf _NPROCESSORS_ONLN`);
    }

    #  See src/utility/src/utility/system.C

    if (exists($ENV{"OMP_NUM_THREADS"})) {
        $ncpu = $ENV{"OMP_NUM_THREADS"};
    }

    if (exists($ENV{"SLURM_JOB_CPUS_PER_NODE"})) {
        $ncpu = $ENV{"SLURM_JOB_CPUS_PER_NODE"};
    }

    if (exists($ENV{"PBS_NCPUS"})) {
        $ncpu = $ENV{"PBS_NCPUS"};
    }

    if (exists($ENV{"PBS_NUM_PPN"})) {
        $ncpu = $ENV{"PBS_NUM_PPN"};
    }

    if (exists($ENV{"NSLOTS"})) {
        $ncpu = $ENV{"NSLOTS"};
    }

    return($ncpu);
}



sub getPhysicalMemorySize () {
    my $os     = $^O;
    my $memory = 1;

    if ($os eq "freebsd") {
        $memory = `/sbin/sysctl -n hw.physmem` / 1024 / 1024 / 1024;
    }

    if ($os eq "darwin") {
        $memory = `/usr/sbin/sysctl -n hw.memsize` / 1024 / 1024 / 1024;
    }

    if ($os eq "linux" || $os eq "cygwin") {
        open(F, "< /proc/meminfo");        #  Way to go, Linux!  Make it easy on us!
        while (<F>) {
            if (m/MemTotal:\s+(\d+)/) {
                $memory = $1 / 1024 / 1024;
            }
        }
        close(F);
    }

    #  See src/utility/src/utility/system.C

    if (exists($ENV{"SLURM_MEM_PER_CPU"}) &&
        exists($ENV{"SLURM_JOB_CPUS_PER_NODE"})) {
        $memory = $ENV{"SLURM_MEM_PER_CPU"} * $ENV{"SLURM_JOB_CPUS_PER_NODE"};
    }

    if (exists($ENV{"SLURM_MEM_PER_NODE"})) {
        $memory = $ENV{"SLURM_MEM_PER_NODE"};
    }

    if (exists($ENV{"PBS_RESC_MEM"})) {
        $memory = $ENV{"PBS_RESC_MEM"} / 1024 / 1024 / 1024;
    }

    return(int($memory + 0.5));  #  Poor man's rounding
}



sub detectLocal () {
    my $ncpu   = getNumberOfCPUs();
    my $maxm   = getPhysicalMemorySize();

    setGlobal("localThreads", $ncpu);       #  Set available resources on local machine.
    setGlobal("localMemory",  $maxm);

    print STDERR "--\n";
    print STDERR "-- Detected $ncpu CPUs and $maxm gigabytes of memory on the local machine.\n";

    my $maxcpu = getGlobal("maxThreads");   #  Get user-supplied limit.
    my $maxmem = getGlobal("maxMemory");

    if (($ncpu < $maxcpu) ||                #  Warn if user limit is more than available resource
        ($maxm < $maxmem)) {                #  but don't adjust.
        print STDERR "--\n";
        print STDERR "-- WARNING: maxThreads=$maxcpu has no effect when only $ncpu CPUs present.\n"              if ($ncpu < $maxcpu);
        print STDERR "-- WARNING: maxMemory=$maxmem has no effect when only $maxm gigabytes memory present.\n"   if ($maxm < $maxmem);
    }
}



sub configureLocal () {

    #  If a grid found, don't setup for local.

    return   if (defined(getGlobal("gridEngine")));

    #  Otherwise, setup for local mode...only there's nothing we need to do
    #  at the moment.  This used to disable stageDirectory (before 29 June
    #  2021) but that was because there was no way to prevent jobs from using
    #  the same stage directory on top of each other.  See
    #  CorrectReads.pm/OverlapMhap.pm for how to do it 'properly'.

    print STDERR "--\n";
    print STDERR "-- Local machine mode enabled; grid support not detected or not allowed.\n";
}
