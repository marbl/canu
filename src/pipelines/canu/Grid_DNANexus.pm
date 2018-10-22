
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
 #    Brian P. Walenz beginning on 2017-FEB-11
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #    Sergey Koren beginning on 2018-MAY-09
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

package canu::Grid_DNANexus;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(detectDNANexus configureDNANexus getDNANexusInstance);

use strict;
use warnings "all";
no  warnings "uninitialized";

use canu::Defaults;
use canu::Execution;

use canu::Grid "formatAllowedResources";



sub buildAvailableNodeList() {
    my %hosts;

    open(F, "dx-get-instance-info.py -p -s 128 | iconv -c -f UTF-8 -t ASCII//TRANSLIT | ");

    my $cpuIdx = 0;
    my $memIdx = 0;
    my $nameIdx = 0;

    while (<F>) {
       if (m/Name\s+|\s+Memory_GB/) {
          my @h = split '\t';

          for (my $ii=0; ($ii < scalar(@h)); $ii++) {
             $cpuIdx  = $ii  if ($h[$ii] eq "CPU_Cores\n");
             $memIdx  = $ii  if ($h[$ii] eq "Memory_GB");
             $nameIdx = $ii  if ($h[$ii] =~ m/Name/);
          }
          last;
       }
    }

    while (<F>) {
        my @v = split '\t', $_;

        my $cpus = $v[$cpuIdx];
        my $mem  = $v[$memIdx] * 0.85; # the machines don't have swap so save some space for OS overhead
        my $name = $v[$nameIdx];

        $mem  = int($mem);
        $hosts{"$cpus-$mem-$name"}++ if ($mem > 0 && $cpus > 0);
     }
     close(F);

     return %hosts;
}

sub getDNANexusInstance($$) {
    my $requestMem       = shift @_;
    my $requestCPU       = shift @_;
    my $instance         = "";
    my $currCPU          = 99999;
    my $currMem          = 99999;

    my %hosts = buildAvailableNodeList();

    foreach my $host (keys %hosts) {
        my @v = split '-', $host;

        my $cpus = $v[0];
        my $mem  = $v[1];

        if ($mem >= $requestMem && $cpus >= $requestCPU && $mem <= $currMem) {
            $instance = $v[2];
            $currCPU = $cpus;
            $currMem = $mem;
         }
    }

    die "-- ERROR:    cannot fine a node to satisfy request for CPU=$requestCPU and MEM=${requestMem}G\n" if (! defined($instance));
    print STDERR "-- Selected $instance for request of ${requestMem}G and $requestCPU.\n";
    return $instance;
}


sub detectDNANexus () {

    return   if ( defined(getGlobal("gridEngine")));   #  Grid not requested.
    return   if (!defined($ENV{'DNANEXUS_HOME'}));     #  Not a DNA Nexus grid

    my $dnanodes = findExecutable("dx-jobutil-new-job");

    return   if (!defined($dnanodes));

    my ($version) = `dx --version`;
    chomp $version;

    print STDERR "-- Detected DNA Nexus '$version' in '$ENV{'DNANEXUS_HOME'}'.\n";
    setGlobal("gridEngine", "DNANEXUS");

    #  DNANexus mode doesn't support (easily) the sequence store check on short reads.
    #  The issue is that we'd need to save the store, ask the user to accept it (and rename),
    #  then continue.  Nothing super tricky, just not done.

    setGlobal("stopOnReadQuality", 0);
}


sub configureDNANexus () {

    return   if (uc(getGlobal("gridEngine")) ne "DNANEXUS");

    my $maxArraySize = 1;    # no array jobs

    #  Probe for the maximum array job size

    setGlobalIfUndef("gridEngineSubmitCommand",              "dx-jobutil-new-job");
    setGlobalIfUndef("gridEngineNameOption",                 "--name");
    setGlobalIfUndef("gridEngineArrayOption",                "-iDX_ARRAY_ID:int=ARRAY_JOBS");
    setGlobalIfUndef("gridEngineArrayName",                  "ARRAY_NAME");
    setGlobalIfUndef("gridEngineArrayMaxJobs",               $maxArraySize);
    setGlobalIfUndef("gridEngineOutputOption",               undef);
    setGlobalIfUndef("gridEnginePropagateCommand",           undef);
    setGlobalIfUndef("gridEngineThreadsOption",              undef);
    setGlobalIfUndef("gridEngineMemoryOption",               "--instance-type MEMORY");
    setGlobalIfUndef("gridEngineNameToJobIDCommand",         undef);
    setGlobalIfUndef("gridEngineNameToJobIDCommandNoArray",  undef);
    setGlobalIfUndef("gridEngineTaskID",                     "DX_ARRAY_ID");
    setGlobalIfUndef("gridEngineArraySubmitID",              undef);
    setGlobalIfUndef("gridEngineJobID",                      "DX_JOB_ID");


    #  Build a list of the resources available in the grid.  This will contain a list with keys of
    #  "#CPUs-#GBs" and values of the number of nodes With such a config.  Later on, we'll use this
    #  to figure out what specific settings to use for each algorithm.
    #
    #  The list is saved in global{"availableHosts"}


    my %hosts = buildAvailableNodeList();


    foreach my $host (keys %hosts) {
        my @v = split '-', $host;

        my $cpus = $v[0];
        my $mem  = $v[1];

        $hosts{"$cpus-$mem"}++    if ($cpus gt 0);
    }

    if (scalar(keys(%hosts)) == 0) {
        my $mm = getGlobal("maxMemory");
        my $mt = getGlobal("maxThreads");

        print STDERR "--\n";
        print STDERR "-- WARNING:  No hosts found in 'instance' report.\n";
        print STDERR "-- WARNING:  Will use maxMemory=$mm and maxThreads=$mt instead.\n";
        print STDERR "-- ERROR:    maxMemory not defined!\n"   if (!defined($mm));
        print STDERR "-- ERROR:    maxThreads not defined!\n"  if (!defined($mt));

        caExit("maxMemory or maxThreads not defined", undef)  if (!defined($mm) || !defined($mt));

        $hosts{"$mt-$mm"}++;
    }

    setGlobal("availableHosts", formatAllowedResources(%hosts, "DNA Nexus"));
}
