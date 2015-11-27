
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
 #    Brian P. Walenz beginning on 2015-MOV-27
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
    setGlobalIfUndef("gridEngineHoldOption",                 "-w \"numended\(\"WAIT_TAG\", \*\)\"");
    setGlobalIfUndef("gridEngineHoldOptionNoArray",          "-w \"done\(\"WAIT_TAG\"\)\"");
    setGlobalIfUndef("gridEngineSyncOption",                 "-K");
    setGlobalIfUndef("gridEngineNameOption",                 "-J");
    setGlobalIfUndef("gridEngineArrayOption",                "");
    setGlobalIfUndef("gridEngineArrayName",                  "ARRAY_NAME\[ARRAY_JOBS\]");
    setGlobalIfUndef("gridEngineOutputOption",               "-o");
    setGlobalIfUndef("gridEngineThreadsOption",              undef);
    setGlobalIfUndef("gridEngineMemoryOption",               undef);
    setGlobalIfUndef("gridEnginePropagateCommand",           "bmodify -w \"done\(\"WAIT_TAG\"\)\"");
    setGlobalIfUndef("gridEngineNameToJobIDCommand",         "bjobs -A -J \"WAIT_TAG\" | grep -v JOBID");
    setGlobalIfUndef("gridEngineNameToJobIDCommandNoArray",  "bjobs -J \"WAIT_TAG\" | grep -v JOBID");
    setGlobalIfUndef("gridEngineTaskID",                     "\$LSB_JOBINDEX");
    setGlobalIfUndef("gridEngineArraySubmitID",              "%I");
    setGlobalIfUndef("gridEngineJobID",                      "LSB_JOBID");

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

        $mem  = $1 * 1024  if ($mem =~ m/(\d+.*\d+)[tT]/);
        $mem  = $1 * 1     if ($mem =~ m/(\d+.*\d+)[gG]/);
        $mem  = $1 / 1024  if ($mem =~ m/(\d+.*\d+)[mM]/);
        $mem  = int($mem);

        $hosts{"$cpus-$mem"}++;
    }
    close(F);
    setGlobal("availableHosts", formatAllowedResources(%hosts, "LSF"));
}

