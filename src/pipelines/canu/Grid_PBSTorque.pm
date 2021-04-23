
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

package canu::Grid_PBSTorque;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(detectPBSTorque configurePBSTorque);

use strict;
use warnings "all";
no  warnings "uninitialized";

use canu::Defaults;
use canu::Execution;

use canu::Grid "formatAllowedResources";



sub detectPBSVersion () {
    my $name    = "";
    my $version = "";

    open(F, "pbsnodes --version 2>&1 |");
    while (<F>) {
        if (m/pbs_version\s+=\s+(.*)/) {
            $name    = "PBSPro";
            $version = $1;
        }
        if (m/Version:\s+(.*)/) {
            $name    = "PBS/Torque";
            $version = $1;
        }
    }
    close(F);

    return($version, $name);
}


sub detectPBSTorque () {

    return   if ( defined(getGlobal("gridEngine")));

    my $pbsnodes = findExecutable("pbsnodes");

    return   if (!defined($pbsnodes));

    my ($version, $name) = detectPBSVersion();

    if (getGlobal("useGrid") eq "0") {
        print STDERR "--\n";
        print STDERR "-- Detected $name '$version' with 'pbsnodes' binary in $pbsnodes.\n";
        print STDERR "--          $name disabled by useGrid=false\n";
    }
    else {
        print STDERR "--\n";
        print STDERR "-- Detected $name '$version' with 'pbsnodes' binary in $pbsnodes.\n";

        setGlobal("gridEngine", ($name eq "PBSPro") ? "PBSPRO" : "PBS");
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

    #  For Torque, see if there is a max array size.
    #  For Pro, set to 1000.

    my $maxArraySize = getGlobal("gridEngineArrayMaxJobs");

    if (!defined($maxArraySize)) {
        $maxArraySize = 1000;

        open(F, "qmgr -c 'p s' |");
        while (<F>) {
            if (m/max_job_array_size\s+=\s+(\d+)/) {   #  Torque
                $maxArraySize = $1;
            }
            if (m/max_array_size\s+=\s+(\d+)/) {       #  PBSPro
                $maxArraySize = $1;
            }
        }
        close(F);
    }

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
    setGlobalIfUndef("gridEngineArrayMaxJobs",               $maxArraySize);
    setGlobalIfUndef("gridEngineOutputOption",               "-o");
    setGlobalIfUndef("gridEngineResourceOption",             "-l nodes=1:ppn=THREADS,mem=MEMORY")      if ($isPro == 0);
    setGlobalIfUndef("gridEngineResourceOption",             "-l select=1:ncpus=THREADS:mem=MEMORY")   if ($isPro == 1);
    setGlobalIfUndef("gridEngineMemoryPerJob",               "1");
    setGlobalIfUndef("gridEngineNameToJobIDCommand",         "qstat -f |grep -F -B 1 WAIT_TAG | grep Id: | grep -F [] |awk '{print \$NF}'");
    setGlobalIfUndef("gridEngineNameToJobIDCommandNoArray",  "qstat -f |grep -F -B 1 WAIT_TAG | grep Id: |awk '{print \$NF}'");
    setGlobalIfUndef("gridEngineTaskID",                     "PBS_ARRAYID")           if ($isPro == 0);
    setGlobalIfUndef("gridEngineTaskID",                     "PBS_ARRAY_INDEX")       if ($isPro == 1);
    setGlobalIfUndef("gridEngineArraySubmitID",              undef);
    setGlobalIfUndef("gridEngineJobID",                      "PBS_JOBID");

    #  Build a list of the resources available in the grid.  This will contain a list with keys
    #  of "#CPUs-#GBs" and values of the number of nodes With such a config.  Later on, we'll use this
    #  to figure out what specific settings to use for each algorithm.
    #
    #  The list is saved in global{"availableHosts"}

    configurePBSTorqueNodes()   if ($isPro == 0);
    configurePBSProNodes()      if ($isPro == 1);
}
