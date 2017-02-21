
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
 #    Brian P. Walenz beginning on 2015-FEB-11
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

package canu::Grid_DNANexus;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(detectDNANexus configureDNANexus);

use strict;

use canu::Defaults;
use canu::Grid;
use canu::Execution;

sub detectDNANexus () {

    return   if ( defined(getGlobal("gridEngine")));   #  Grid not requested.
    return   if (!defined($ENV{'DNA_NEXUS'}));         #  Not a DNA Nexus grid

    print STDERR "-- Detected DNA Nexus '...some-version...'.\n";
    setGlobal("gridEngine", "DNANEXUS");

    #  DNANexus mode doesn't support (easily) the gatekeeper check on short reads.
    #  The issue is that we'd need to save the store, ask the user to accept it (and rename),
    #  then continue.  Nothing super tricky, just not done.

    setGlobal("stopOnReadQuality", 0);
}


sub configureDNANexus () {

    return   if (uc(getGlobal("gridEngine")) ne "DNANEXUS");

    my $maxArraySize = 65535;

    #  Probe for the maximum array job size

    setGlobalIfUndef("gridEngineSubmitCommand",              "");
    setGlobalIfUndef("gridEngineNameOption",                 "");
    setGlobalIfUndef("gridEngineArrayOption",                "");
    setGlobalIfUndef("gridEngineArrayName",                  "");
    setGlobalIfUndef("gridEngineArrayMaxJobs",               $maxArraySize);
    setGlobalIfUndef("gridEngineOutputOption",               "");
    setGlobalIfUndef("gridEnginePropagateCommand",           "");
    setGlobalIfUndef("gridEngineThreadsOption",              undef);
    setGlobalIfUndef("gridEngineMemoryOption",               undef);
    setGlobalIfUndef("gridEngineNameToJobIDCommand",         undef);
    setGlobalIfUndef("gridEngineNameToJobIDCommandNoArray",  undef);
    setGlobalIfUndef("gridEngineTaskID",                     "");
    setGlobalIfUndef("gridEngineArraySubmitID",              "");
    setGlobalIfUndef("gridEngineJobID",                      "");


    my %hosts;

    #  Probe for how to request multiple CPUs on each node, set 

    #  .
    #  .
    #  .

    #  Probe for how to reserve memory on each node

    #  .
    #  .
    #  .

    #  Build a list of the resources available in the grid.  This will contain a list with keys of
    #  "#CPUs-#GBs" and values of the number of nodes With such a config.  Later on, we'll use this
    #  to figure out what specific settings to use for each algorithm.
    #
    #  The list is saved in global{"availableHosts"}

    #  .
    #  $hosts{"4-32"} = 15;   #  15 machines with 4 CPUs and 32gb memory
    #  .

    setGlobal("availableHosts", formatAllowedResources(%hosts, "DNA Nexus"));
}
