
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

package canu::Grid_LSF;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(detectLSF configureLSF);

use strict;

use canu::Defaults;
use canu::Execution;
use canu::Grid;

sub detectLSF () {

    return   if ( defined(getGlobal("gridEngine")));

    my $bsub = findExecutable("bsub");

    return   if (!defined($bsub));

    print STDERR "-- Detected LSF with 'bsub' binary in $bsub.\n";
    setGlobal("gridEngine", "LSF");
}


sub configureLSF () {

    return   if (uc(getGlobal("gridEngine")) ne "LSF");

    setGlobalIfUndef("gridEngineSubmitCommand",              "bsub");
    setGlobalIfUndef("gridEngineHoldOption",                 "-w \"ended\(\"WAIT_TAG\"\)\"");
    setGlobalIfUndef("gridEngineHoldOptionNoArray",          "-w \"done\(\"WAIT_TAG\"\)\"");
    setGlobalIfUndef("gridEngineSyncOption",                 "-K");
    setGlobalIfUndef("gridEngineNameOption",                 "-J");
    setGlobalIfUndef("gridEngineArrayOption",                "");
    setGlobalIfUndef("gridEngineArrayName",                  "ARRAY_NAME\[ARRAY_JOBS\]");
    setGlobalIfUndef("gridEngineArrayMaxJobs",               65535);
    setGlobalIfUndef("gridEngineOutputOption",               "-o");
    setGlobalIfUndef("gridEngineThreadsOption",              "-n THREADS");
    setGlobalIfUndef("gridEngineMemoryOption",               "-M MEMORY");
    setGlobalIfUndef("gridEnginePropagateCommand",           "bmodify -w \"done\(\"WAIT_TAG\"\)\"");
    setGlobalIfUndef("gridEngineNameToJobIDCommand",         "bjobs -A -J \"WAIT_TAG\" | grep -v JOBID");
    setGlobalIfUndef("gridEngineNameToJobIDCommandNoArray",  "bjobs -J \"WAIT_TAG\" | grep -v JOBID");
    setGlobalIfUndef("gridEngineTaskID",                     "LSB_JOBINDEX");
    setGlobalIfUndef("gridEngineArraySubmitID",              "%I");
    setGlobalIfUndef("gridEngineJobID",                      "LSB_JOBID");

    #
    #  LSF has variation in the units used to request memory
    #  They are defined by the LSF_UNIT_FOR_LIMITS variable in lsf.conf
    #  Poll and see if we can find it
    #

    my $memUnits = undef;

    open(F, "lsadmin showconf lim |");

    my $s = <F>;  #  cluster name
    my $d = <F>;  #  dat/time

    while (<F>) {
        my @v = split '=', $_;
        if ($v[0] =~ m/LSF_UNIT_FOR_LIMITS/) {
            $memUnits = "t" if ($v[1] =~ m/[tT]/);
            $memUnits = "g" if ($v[1] =~ m/[gG]/);
            $memUnits = "m" if ($v[1] =~ m/[mM]/);
            $memUnits = "k" if ($v[1] =~ m/[kK]/);
        }
    }

    close(F);

    if (!defined($memUnits)) {
        print STDERR "-- Warning: unknown memory units for grid engine LSF assuming KB\n";
        $memUnits = "k";
    }

    #  Build a list of the resources available in the grid.  This will contain a list with keys
    #  of "#CPUs-#GBs" and values of the number of nodes With such a config.  Later on, we'll use this
    #  to figure out what specific settings to use for each algorithm.
    #
    #  The list is saved in global{"availableHosts"}
    #
    #  !!! UNTESTED !!
    #
    my %hosts;

    open(F, "lshosts |");

    my $h = <F>;  #  header

    my @h = split '\s+', $h;

    my $cpuIdx  = 4;
    my $memIdx  = 5;

    for (my $ii=0; ($ii < scalar(@h)); $ii++) {
        $cpuIdx  = $ii  if ($h[$ii] eq "ncpus");
        $memIdx  = $ii  if ($h[$ii] eq "maxmem");
    }

    while (<F>) {
        my @v = split '\s+', $_;

        my $cpus  = $v[$cpuIdx];
        my $mem   = $v[$memIdx];

        # if we failed to find the units from the configuration, inherit it from the lshosts output
        if (!defined($memUnits)) {
            $memUnits = "t" if ($mem =~ m/(\d+.*\d+)[tT]/);
            $memUnits = "g" if ($mem =~ m/(\d+.*\d+)[gG]/);
            $memUnits = "m" if ($mem =~ m/(\d+.*\d+)[mM]/);
            $memUnits = "k" if ($mem =~ m/(\d+.*\d+)[kK]/);
        }

        $mem  = $1 * 1024         if ($mem =~ m/(\d+.*\d+)[tT]/);
        $mem  = $1 * 1            if ($mem =~ m/(\d+.*\d+)[gG]/);
        $mem  = $1 / 1024         if ($mem =~ m/(\d+.*\d+)[mM]/);
        $mem  = $1 / 1024 / 1024  if ($mem =~ m/(\d+.*\d+)[kK]/);
        $mem  = int($mem);

        $hosts{"$cpus-$mem"}++    if ($cpus gt 0);
    }
    close(F);
    setGlobal("availableHosts", formatAllowedResources(%hosts, "LSF"));
    setGlobal("gridEngineMemoryUnits", $memUnits);
    print STDERR "-- \n";
    print STDERR "-- On LSF detected memory is requested in " . uc(${memUnits}) . "B\n";
    print STDERR "-- \n";
}

