
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

package canu::Grid_SGE;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(detectSGE configureSGE);

use strict;

use canu::Defaults;
use canu::Grid;

sub detectSGE () {

    return   if ( defined(getGlobal("gridEngine")));
    return   if (!defined($ENV{'SGE_ROOT'}));

    print STDERR "-- Detected Sun Grid Engine in '$ENV{'SGE_ROOT'}/$ENV{'SGE_CELL'}'.\n";
    setGlobal("gridEngine", "SGE");
}


sub configureSGE () {

    return   if (uc(getGlobal("gridEngine")) ne "SGE");

    setGlobalIfUndef("gridEngineSubmitCommand",              "qsub");
    setGlobalIfUndef("gridEngineHoldOption",                 "-hold_jid \"WAIT_TAG\"");
    setGlobalIfUndef("gridEngineHoldOptionNoArray",          undef);
    setGlobalIfUndef("gridEngineSyncOption",                 "-sync y");
    setGlobalIfUndef("gridEngineNameOption",                 "-cwd -N");
    setGlobalIfUndef("gridEngineArrayOption",                "-t ARRAY_JOBS");
    setGlobalIfUndef("gridEngineArrayName",                  "ARRAY_NAME");
    setGlobalIfUndef("gridEngineOutputOption",               "-j y -o");
    setGlobalIfUndef("gridEnginePropagateCommand",           "qalter -hold_jid \"WAIT_TAG\"");
    setGlobalIfUndef("gridEngineThreadsOption",              undef);  #"-pe threads THREADS");
    setGlobalIfUndef("gridEngineMemoryOption",               undef);  #"-l mem_free=MEMORY");
    setGlobalIfUndef("gridEngineNameToJobIDCommand",         undef);
    setGlobalIfUndef("gridEngineNameToJobIDCommandNoArray",  undef);
    setGlobalIfUndef("gridEngineTaskID",                     "SGE_TASK_ID");
    setGlobalIfUndef("gridEngineArraySubmitID",              "\\\$TASK_ID");
    setGlobalIfUndef("gridEngineJobID",                      "JOB_ID");

    #  Try to figure out the name of the threaded job execution environment.
    #  It's the one with allocation_rule of $pe_slots.

    if (!defined(getGlobal("gridEngineThreadsOption"))) {
        my @env = `qconf -spl`;  chomp @env;
        my $thr = undef;
        my $nth = 0;

        foreach my $env (@env) {
            my $ar = 0;
            my $cs = 0;
            my $jf = 0;

            open(F, "qconf -sp $env |");
            while (<F>) {
                $ar = 1  if (m/allocation_rule.*pe_slots/);
                $cs = 1  if (m/control_slaves.*FALSE/);
                $jf = 1  if (m/job_is_first_task.*TRUE/);
            }
            close(F);

            if (($ar == 1) && ($cs == 1) && ($jf == 1)) {
                $thr = $env;
                $nth++;
            }
        }

        if ($nth == 1) {
            print STDERR "-- Detected Grid Engine environment '$thr'.\n";

            setGlobal("gridEngineThreadsOption", "-pe $thr THREADS");

        } else {
            print STDERR "-- WARNING:  Couldn't determine the SGE parallel environment to run multi-threaded codes.\n";
            print STDERR "--          Set 'gridEngineThreadsOption' manually; example: '-pe threaded THREADS'.\n";
        }
    } else {
        if (getGlobal("gridEngineThreadsOption") =~ m/^-pe\s+(.*)$/) {
            print STDERR "-- User supplied Grid Engine environment '$1'.\n";
        } else {
            caFailure("Couldn't parse gridEngineThreadsOption='" . getGlobal("gridEngineThreadsOption") . "'", undef);
        }
    }

    #  Try to figure out the name of the memory resource.

    if (!defined(getGlobal("gridEngineMemoryOption"))) {
        my $mem = undef;
        my $nmm = 0;

        open(F, "qconf -sc |");
        while (<F>) {
            my @vals = split '\s+', $_;

            next  if ($vals[5] ne "YES");        #  Not a consumable resource.
            next  if ($vals[2] ne "MEMORY");     #  Not a memory resource.
            next  if ($vals[0] =~ m/swap/);      #  Don't care about swap.
            next  if ($vals[0] =~ m/virtual/);   #  Don't care about vm space.

            $mem .= " " if (defined($mem));
            $mem .= $vals[0];
            $nmm++;
        }
        close(F);

        if      ($nmm == 1) {
            print STDERR "-- Detected Grid Engine consumable '$mem'.\n";

            setGlobal("gridEngineMemoryOption", "-l $mem=MEMORY");

        } elsif ($nmm > 1) {
            print STDERR "-- WARNING:  Couldn't determine the SGE resource to request memory.\n";
            print STDERR "--          Found $nmm choices: $mem\n";
            print STDERR "--          Set 'gridEngineMemoryOption' manually; example: '-l mem_free=MEMORY'.\n";

        } else {
            print STDERR "-- WARNING:  Couldn't determine the SGE resource to request memory.\n";
            print STDERR "--          Set 'gridEngineMemoryOption' manually; example: '-l mem_free=MEMORY'.\n";
        }
    } else {
        if (getGlobal("gridEngineMemoryOption") =~ m/^-l\s+(.*)=MEMORY$/) {
            print STDERR "-- User supplied Grid Engine consumable '$1'.\n";
        } else {
            caFailure("Couldn't parse gridEngineMemoryOption='" . getGlobal("gridEngineMemoryOption") . "'", undef);
        }
    }

    #  Build a list of the resources available in the grid.  This will contain a list with keys
    #  of "#CPUs-#GBs" and values of the number of nodes With such a config.  Later on, we'll use this
    #  to figure out what specific settings to use for each algorithm.
    #
    #  The list is saved in global{"availableHosts"}

    my %hosts;
    my $hosts = "";

    open(F, "qhost |");

    my $h = <F>;  #  Header
    my $b = <F>;  #  Table bar

    my @h = split '\s+', $h;

    my $cpuIdx = 2;
    my $memIdx = 4;

    for (my $ii=0; ($ii < scalar(@h)); $ii++) {
        $cpuIdx = $ii  if ($h[$ii] eq "NCPU");
        $memIdx = $ii  if ($h[$ii] eq "MEMTOT");
    }

    while (<F>) {
        my @v = split '\s+', $_;

        next if ($v[3] eq "-");  #  Node disabled or otherwise not available

        my $cpus = $v[$cpuIdx];
        my $mem  = $v[$memIdx];

        $mem  = $1 * 1024  if ($mem =~ m/(\d+.*\d+)[tT]/);
        $mem  = $1 * 1     if ($mem =~ m/(\d+.*\d+)[gG]/);
        $mem  = $1 / 1024  if ($mem =~ m/(\d+.*\d+)[mM]/);
        $mem  = int($mem);

        $hosts{"$cpus-$mem"}++    if ($cpus gt 0);
    }
    close(F);
    setGlobal("availableHosts", formatAllowedResources(%hosts, "Sun Grid Engine"));
}
