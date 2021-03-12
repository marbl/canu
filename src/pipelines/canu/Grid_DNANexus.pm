
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
    my $nameIdx = 0;
    my $memIdx  = 0;
    my $stoIdx  = 0;
    my $cpuIdx  = 0;

    #  Get all instances valid for this project (-p) that
    #  have at least 128 GB storage.

    open(F, "dx-get-instance-info.py -p -s 128 | iconv -c -f UTF-8 -t ASCII//TRANSLIT | ");
    while (<F>) {
        s/^\s+//;
        s/\s+$//;

        my @v = split '\s+', $_;

        #  Parse the header to figure out which columns have which data.
        #    'Name	Memory_GB	Storage_GB	CPU_Cores'
        if (m/Memory_GB/) {
            for (my $ii=0; ($ii < scalar(@v)); $ii++) {
                $nameIdx = $ii  if ($v[$ii] eq "Name");
                $memIdx  = $ii  if ($v[$ii] eq "Memory_GB");
                $stoIdx  = $ii  if ($v[$ii] eq "Storage_GB");
                $cpuIdx  = $ii  if ($v[$ii] eq "CPU_Cores");
            }
        }

        #  Parse a line for data, then save it if valid.
        else {
            my $name =     $v[$nameIdx];
            my $mem  = int($v[$memIdx] * 0.85);   #  No swap; leave space for the OS
            my $sto  =     $v[$stoIdx];
            my $cpus =     $v[$cpuIdx];

            if (($mem > 0) && ($cpus > 0)) {
                $hosts{"$cpus-$mem-$sto-$name"}++;
            }
        }
    }
    close(F);

    return %hosts;
}



sub getDNANexusInstance($$) {
    my $reqMem = shift @_;
    my $reqCPU = shift @_;
    my $sel   = undef;
    my $selCPU;
    my $selMem;
    my $selSto;

    my %hosts = buildAvailableNodeList();

    foreach my $host (keys %hosts) {
        my ($cpu, $mem, $sto, $name) = split '-', $host;

        if (($cpu < $reqCPU) ||     #  Instance doesn't have enough
            ($mem < $reqMem)) {     #  CPUs or memory.
            #print STDERR "-- Instance '$name' with $mem GB and $cpu CPUs rejected; too small.\n";
            next;
        }

        #  If nothing defined yet, accept the first valid instance.

        if (!defined($sel)) {
            #print STDERR "-- Instance '$name' with $mem GB and $cpu CPUs accepted.\n";
            ($sel, $selCPU, $selMem, $selSto) = ($name, $cpu, $mem, $sto);
            next;
        }

        #  If the instance is bigger than the selected, skip it.

        if (($cpu > $selCPU) ||
            ($mem > $selMem)) {
            #print STDERR "-- Instance '$name' with $mem GB and $cpu CPUs rejected; extra CPU or memory.\n";
            next;
        }

        if (($cpu == $selCPU) &&
            ($mem == $selMem) &&
            ($sto  > $selSto)) {
            print STDERR "-- Instance '$name' with $mem GB and $cpu CPUs rejected; extra storage.\n";
            next;
        }

        #  Otherwise, use this instance.

        print STDERR "-- Instance '$name' with $mem GB and $cpu CPUs accepted.\n";
        ($sel, $selCPU, $selMem, $selSto) = ($name, $cpu, $mem, $sto);
    }

    if (!defined($sel)) {
        print STDERR "--\n";
        print STDERR "-- ERROR: Cannot find an instance to satisfy request for $reqMem GB and $reqCPU CPUs.\n";
        print STDERR "--\n";

        caExit("no appropriate instance found", undef);
    }

    print STDERR "-- Selected $sel ($selMem GB and $selCPU CPUs) for request of $reqMem GB and $reqCPU CPUs.\n";

    return($sel);
}



sub detectDNANexus () {

    return   if ( defined(getGlobal("gridEngine")));   #  Grid not requested.
    return   if (!defined($ENV{'DNANEXUS_HOME'}));     #  Not a DNA Nexus grid

    my $dnanodes = findExecutable("dx-jobutil-new-job");

    return   if (!defined($dnanodes));

    my ($version) = `dx --version`;
    chomp $version;

    print STDERR "--\n";
    print STDERR "-- Detected DNA Nexus '$version' in '$ENV{'DNANEXUS_HOME'}'.\n";

    setGlobal("gridEngine", "DNANEXUS");
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
    setGlobalIfUndef("gridEngineResourceOption",             "--instance-type MEMORY");
    setGlobalIfUndef("gridEngineNameToJobIDCommand",         undef);
    setGlobalIfUndef("gridEngineNameToJobIDCommandNoArray",  undef);
    setGlobalIfUndef("gridEngineTaskID",                     "DX_ARRAY_ID");
    setGlobalIfUndef("gridEngineArraySubmitID",              undef);
    setGlobalIfUndef("gridEngineJobID",                      "DX_JOB_ID");

    #  Force this to be enabled, regardless of what the user thinks it should
    #  be.  We always need to have memory per job.
    setGlobal("gridEngineMemoryPerJob", 1);

    #  Build a list of the resources available in the grid.  This will contain a list with keys of
    #  "#CPUs-#GBs" and values of the number of nodes With such a config.  Later on, we'll use this
    #  to figure out what specific settings to use for each algorithm.
    #
    #  The list is saved in global{"availableHosts"}

    my %hosts = buildAvailableNodeList();

    if (scalar(keys(%hosts)) == 0) {
        print STDERR "--\n";
        print STDERR "-- ERROR:  No instances found in output of 'dx-get-instance-info.py -p -s 128'.\n";
        print STDERR "--\n";

        caExit("no instances found", undef);
    }

    setGlobal("availableHosts", formatAllowedResources(%hosts, "DNA Nexus"));
}
