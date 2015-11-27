
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
 #  This file is derived from:
 #
 #    src/pipelines/ca3g/Defaults.pm
 #
 #  Modifications by:
 #
 #    Brian P. Walenz from 2015-FEB-27 to 2015-SEP-21
 #      are Copyright 2015 Battelle National Biodefense Institute, and
 #      are subject to the BSD 3-Clause License
 #
 #    Brian P. Walenz beginning on 2015-OCT-21
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

package canu::Configure;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(configureAssembler);

use strict;
use Carp qw(cluck);
use Sys::Hostname;

use canu::Defaults;


sub expandRange ($$) {
    my $var = shift @_;
    my $val = shift @_;

    my @v = split ',', $val;
    my @r;

    foreach my $v (@v) {
        if      ($v =~ m/^(\d+)$/) {
            push @r, $1;

        } elsif ($v =~ m/^(\d+)([kKmMgGtT]{0,1})-(\d+)([kKmMgGtT]{0,1})$/) {
            my $b = adjustMemoryValue("$1$2");
            my $e = adjustMemoryValue("$3$4");
            my $s = adjustMemoryValue("1$4");     #  Yup, not a bug: 1-4g == 1g-4g == 1g,2g,3g,4g

            for (my $ii=$b; $ii<=$e; $ii += $s) {
                push @r, $ii;
            }

        } elsif ($v =~ m/^(\d+[kKmMgGtT]{0,1})-(\d+[kKmMgGtT]{0,1}):(\d+[kKmMgGtT]{0,1})$/) {
            my $b = adjustMemoryValue($1);  #  Unlike above, all three need to have multipliers
            my $e = adjustMemoryValue($2);
            my $s = adjustMemoryValue($3);

            for (my $ii=$b; $ii<=$e; $ii += $s) {
                push @r, $ii;
            }

        } else {
            caExit("can't parse '$var' entry '$v'", undef);
        }
    }

    #print "$var = ";
    #foreach my $r (@r) {
    #    print "$r ";
    #}
    #print "\n";

    return(@r);
}


#  Side effect!  This will RESET the $global{} parameters to the computed value.  This lets
#  the rest of canu - in particular, the part that runs the jobs - use the correct value.  Without
#  resetting, I'd be making code changes all over the place to support the values returned.

sub getAllowedResources ($$$$) {
    my $tag  = shift @_;  #  Variant, e.g., "cor", "utg"
    my $alg  = shift @_;  #  Algorithm, e.g., "mhap", "ovl"
    my $err  = shift @_;  #  Report of things we can't run.
    my $all  = shift @_;  #  Report of things we can run.

    #  If no grid, or grid not enabled, everything falls under 'lcoal'.

    my $class = ((getGlobal("useGrid") == 0) || (getGlobal("gridEngine") eq undef)) ? "local" : "grid";

    #  Figure out limits.

    my $maxMemory    = getGlobal("maxMemory");
    my $maxThreads   = getGlobal("maxThreads");

    my $taskMemory   = getGlobal("${tag}${alg}Memory");   #  Algorithm limit, "utgovlMemory", etc.
    my $taskThreads  = getGlobal("${tag}${alg}Threads");  #

    #  The task limits MUST be defined.

    caExit("${tag}${alg}Memory is not defined", undef)   if (!defined($taskMemory));
    caExit("${tag}${alg}Threads is not defined", undef)  if (!defined($taskThreads));

    #  If the maximum limits aren't set, default to 'unlimited' (for the grid; we'll effectively filter
    #  by the number of jobs we can fit on the hosts) or to the current hardware limits.

    $maxMemory  = (($class eq "grid") ? 1024 * 1024 : getPhysicalMemorySize())  if (!defined($maxMemory));    #  1 PB memory!
    $maxThreads = (($class eq "grid") ? 1024        : getNumberOfCPUs())        if (!defined($maxThreads));   #  1 k  cores!

    #  Build a list of the available hardware configurations we can run on.  If grid, we get this
    #  from the list of previously discovered hosts.  If local, it's just this machine.

    my @gridCor;  #  Number of cores
    my @gridMem;  #  GB's of memory
    my @gridNum;  #  Number of nodes

    if ($class eq "grid") {
        my @grid = split '\0', getGlobal("availableHosts");

        foreach my $g (@grid) {
            my ($cpu, $mem, $num) = split '-', $g;

            push @gridCor, $cpu;
            push @gridMem, $mem;
            push @gridNum, $num;
        }
    } else {
        push @gridCor, $maxThreads;
        push @gridMem, $maxMemory;
        push @gridNum, 1;
    }

    #  The task usually has multiple choices, and we have a little optimization problem to solve.  For each
    #  pair of memory/threads, compute three things:
    #    a) how many processes we can get running
    #    b) how many cores we can get running
    #    c) how much memory we can consume
    #  We then (typically) want to maximize the number of cores we can get running.
    #  Other options would be number of cores * amount of memory.

    my @taskMemory  = expandRange("${tag}${alg}Memory",  $taskMemory);
    my @taskThreads = expandRange("${tag}${alg}Threads", $taskThreads);

    #  Filter out task settings that can't be run based on the gridMemory/gridThreads or masterMemory/masterThreads setting.
    #  (actually, this just reports those that would be filtered; the actual filtering is inline in the algorithm)

    my $ignoreM;
    my $ignoreT;

    foreach my $m (@taskMemory) {
        $m = adjustMemoryValue($m);
    }

    foreach my $m (@taskMemory) {
        next  if ($m <= $maxMemory);
        $ignoreM .= ","  if (defined($ignoreM));
        $ignoreM .= "${m}g";
    }
    foreach my $t (@taskThreads) {
        next  if ($t <= $maxThreads);
        $ignoreT .= ","  if (defined($ignoreT));
        $ignoreT .= "$t";
    }

    #  Too verbose with long value lists
    #
    #if      (defined($ignoreM) && defined($ignoreT)) {
    #    $err .= "-- Can't use ${tag}${alg}Memory=$ignoreM and ${tag}${alg}Threads=$ignoreT because of maxMemory=${maxMemory}g and maxThreads=$maxThreads limits.\n";
    #
    #} elsif (defined($ignoreM)) {
    #    $err .= "-- Can't use ${tag}${alg}Memory=$ignoreM because of maxMemory=${maxMemory}g limit.\n";
    #
    #} elsif (defined($ignoreT)) {
    #    $err .= "-- Can't use ${tag}${alg}Threads=$ignoreT because of maxThreads=$maxThreads limit.\n";
    #}

    #  Find task memory/thread settings that will maximize the number of cores running.  This used
    #  to also compute best as 'cores * memory' but that is better handled by ordering the task
    #  settings parameters.  The example below will pick the largest (last) configuration that
    #  maximizes core utilization:
    #
    #    taskThreads = 4,8,32,64
    #    taskMemory  = 16g,32g,64g

    my ($bestCores,  $bestCoresM,  $bestCoresT)  = (0, undef, undef);

    foreach my $m (@taskMemory) {
        foreach my $t (@taskThreads) {

            next  if ($m > $maxMemory);   #  Bail if either of the suggest settings are
            next  if ($t > $maxThreads);  #  larger than the maximum allowed.

            my $processes = 0;
            my $cores     = 0;
            my $memory    = 0;

            #  For a job using $m GB memory and $t threads, we can compute how many processes will
            #  fit on each node in our set of available machines.  The smaller of the two is then
            #  the number of processes we can run on this node.

            for (my $ii=0; $ii<scalar(@gridCor); $ii++) {
                my $np_cpu = $gridNum[$ii] * int($gridCor[$ii] / $t);  #  Each process uses $t cores, node has $gridCor[$ii] cores available.
                my $np_mem = $gridNum[$ii] * int($gridMem[$ii] / $m);  #  Same idea.

                my $np = ($np_cpu < $np_mem) ? $np_cpu : $np_mem;

                $processes += $np;
                $cores     += $np * $t;
                $memory    += $np * $m;
            }

            #  Save the best one seen so far.

            if ($bestCores <= $cores) {
                $bestCores  = $cores;
                $bestCoresM = $m;
                $bestCoresT = $t;
            }
        }
    }

    if (!defined($bestCoresM)) {
        print STDERR "--\n";
        print STDERR "-- Task $tag$alg can't run on any available machines.\n";
        print STDERR "-- It is requesting ", getGlobal("${tag}${alg}Memory"), " GB memory and ", getGlobal("${tag}${alg}Threads"), " threads.\n";
        print STDERR "-- See above for hardware limits.\n";
        print STDERR "--\n";

        caExit("task $tag$alg failed to find a configuration to run on", undef);
    }

    $taskMemory  = $bestCoresM;
    $taskThreads = $bestCoresT;

    #  Check for stupidity.

    caExit("invalid taskThread=$taskMemory; maxMemory=$maxMemory", undef)     if ($taskMemory > $maxMemory);
    caExit("invalid taskThread=$taskThreads; maxThreads=$maxThreads", undef)  if ($taskThreads > $maxThreads);

    #  Reset the global values for later use.

    setGlobal("${tag}${alg}Memory",  $taskMemory);
    setGlobal("${tag}${alg}Threads", $taskThreads);

    #  Finally, reset the concurrency (if we're running locally) so we don't swamp our poor workstation.

    my $concurrent = undef;

    if ($class eq "local") {
        my $nc = int($maxThreads / $taskThreads);

        if (($taskThreads * getGlobal("${tag}${alg}Concurrency") > $maxThreads)) {
            $err .= "-- Reset concurrency from ", getGlobal("${tag}${alg}Concurrency"), " to $nc.\n";
            setGlobal("${tag}${alg}Concurrency", $nc);
        }

        if (!defined(getGlobal("${tag}${alg}Concurrency"))) {
            setGlobal("${tag}${alg}Concurrency", $nc);
        }

        $concurrent = getGlobal("${tag}${alg}Concurrency");
    }

    #  And report.

    my $nam;

    if    ($alg eq "bat")      {  $nam = "bogart (unitigger)"; }
    elsif ($alg eq "cns")      {  $nam = "utgcns (consensus"; }
    elsif ($alg eq "cor")      {  $nam = "falcon_sense (read correction)"; }
    elsif ($alg eq "meryl")    {  $nam = "meryl (k-mer counting)"; }
    elsif ($alg eq "oea")      {  $nam = "overlap error adjustment"; }
    elsif ($alg eq "ovb")      {  $nam = "overlap store parallel bucketizer"; }
    elsif ($alg eq "ovlStore") {  $nam = "overlap store sequential building"; }
    elsif ($alg eq "ovs")      {  $nam = "overlap store parallel sorting"; }
    elsif ($alg eq "red")      {  $nam = "read error detection (overlap error adjustment)"; }
    elsif ($alg eq "mhap")     {  $nam = "mhap (overlapper)"; }
    elsif ($alg eq "ovl")      {  $nam = "overlapper"; }
    else {
        caFailure("unknown task '$alg' in getAllowedResources().", undef);
    }

    $all .= "-- Allowed to";
    $all .= " run " . substr("   $concurrent", -3) . " job" . (($concurrent == 1) ? " " : "s") . " concurrently,"  if (defined($concurrent));
    $all .= " run under grid control,"                                                                             if (!defined($concurrent));
    $all .= " and use up to " . substr("   $taskThreads", -3) . " compute thread" . (($taskThreads == 1) ? " " : "s");
    $all .= " and " . substr("   $taskMemory", -4) . " GB memory for stage '$nam'.\n";

    return($err, $all);
}




#  Converts number with units to gigabytes.  If no units, gigabytes is assumed.
sub adjustMemoryValue ($) {
    my $val = shift @_;

    return(undef)                     if (!defined($val));

    return($1)                        if ($val =~ m/^(\d+\.{0,1}\d*)$/);
    return($1 / 1024 / 1024)          if ($val =~ m/^(\d+\.{0,1}\d*)[kK]$/);
    return($1 / 1024)                 if ($val =~ m/^(\d+\.{0,1}\d*)[mM]$/);
    return($1)                        if ($val =~ m/^(\d+\.{0,1}\d*)[gG]$/);
    return($1 * 1024)                 if ($val =~ m/^(\d+\.{0,1}\d*)[tT]$/);
    return($1 * 1024 * 1024)          if ($val =~ m/^(\d+\.{0,1}\d*)[pP]$/);

    die "Invalid memory value '$val'\n";
}


#  Converts gigabytes to number with units.
sub displayMemoryValue ($) {
    my $val = shift @_;

    return(($val * 1024 * 1024)        . "k")   if ($val < adjustMemoryValue("1m"));
    return(($val * 1024)               . "m")   if ($val < adjustMemoryValue("1g"));
    return(($val)                      . "g")   if ($val < adjustMemoryValue("1t"));
    return(($val / 1024)               . "t");
}


#  Converts number with units to bases.
sub adjustGenomeSize ($) {
    my $val = shift @_;

    return(undef)               if (!defined($val));

    return($1)                  if ($val =~ m/^(\d+\.{0,1}\d*)$/i);
    return($1 * 1000)           if ($val =~ m/^(\d+\.{0,1}\d*)[kK]$/i);
    return($1 * 1000000)        if ($val =~ m/^(\d+\.{0,1}\d*)[mM]$/i);
    return($1 * 1000000000)     if ($val =~ m/^(\d+\.{0,1}\d*)[gG]$/i);
    return($1 * 1000000000000)  if ($val =~ m/^(\d+\.{0,1}\d*)[tT]$/i);

    die "Invalid genome size '$val'\n";
}


#  Converts bases to number with units.
sub displayGenomeSize ($) {
    my $val = shift @_;

    return(($val))                        if ($val < adjustGenomeSize("1k"));
    return(($val / 1000)          . "k")  if ($val < adjustGenomeSize("1m"));
    return(($val / 1000000)       . "m")  if ($val < adjustGenomeSize("1g"));
    return(($val / 1000000000)    . "g")  if ($val < adjustGenomeSize("1t"));
    return(($val / 1000000000000) . "t");
}





#
#  If minMemory or minThreads isn't defined, pick a reasonable pair based on genome size.
#

sub configureAssembler () {

    #  First, parse units on things the user possibly set.

    setGlobal("genomeSize", adjustGenomeSize(getGlobal("genomeSize")));

    setGlobal("minMemory", adjustMemoryValue(getGlobal("minMemory")));
    setGlobal("maxMemory", adjustMemoryValue(getGlobal("maxMemory")));

    #

    if (!defined(getGlobal("minMemory")) || !defined(getGlobal("minThreads"))) {
        my  $minMemoryD = "";
        my  $minMemory  = 0;
        my  $minThreads = 0;

        #  Based on genome size, pick some arbitrary minimums to meet.

        if      (getGlobal("genomeSize") < adjustGenomeSize("40m")) {
            $minMemoryD = "8g";
            $minThreads =  4;

        } elsif (getGlobal("genomeSize") < adjustGenomeSize("500m")) {
            $minMemoryD = "32g";
            $minThreads =  4;

        } elsif (getGlobal("genomeSize") < adjustGenomeSize("2g")) {
            $minMemoryD = "128g";
            $minThreads =  8;

        } elsif (getGlobal("genomeSize") < adjustGenomeSize("5g")) {
            $minMemoryD = "256g";
            $minThreads =  16;

        } else {
            $minMemoryD = "512g";
            $minThreads =  32;
        }

        $minMemory = adjustMemoryValue($minMemoryD);

        print STDERR "--\n";
        print STDERR "-- For genomeSize ", displayGenomeSize(getGlobal("genomeSize")), ", $minMemoryD memory and $minThreads threads seems reasonable.\n";

        #  Set the minMemory/minThreads to the minimum if not defined.

        if (!defined(getGlobal("minMemory"))) {
            print STDERR "--   minMemory  not defined, set to $minMemoryD\n";
            setGlobal("minMemory", $minMemory);
        } else {
            print STDERR "--   minMemory  defined, current value ", displayMemoryValue(getGlobal("minMemory")) . "\n";
        }

        if (!defined(getGlobal("minThreads"))) {
            print STDERR "--   minThreads not defined, set to $minThreads\n";
            setGlobal("minThreads", $minThreads);
        } else {
            print STDERR "--   minThreads defined, current value ", getGlobal("minThreads") . "\n";
        }

        my $reqMemory  = getGlobal("minMemory");
        my $reqThreads = getGlobal("minThreads");

        if     (($reqMemory < $minMemory) && ($reqThreads < $minThreads)) {
            print STDERR "--\n";
            print STDERR "--  WARNING:  supplied minMemory and minThreads are both smaller than suggested!\n";
        } elsif ($reqMemory < $minMemory) {
            print STDERR "--\n";
            print STDERR "--  WARNING:  supplied minMemory is smaller than suggested!\n";
        } elsif ($reqThreads < $minThreads) {
            print STDERR "--\n";
            print STDERR "--  WARNING:  supplied minMemory is smaller than suggested!\n";
        }

        #
        #  Now, tricky.  We need to decide if we can run a job with reqMemory/reqThreads.  We've
        #  already decided that these are sufficient to run the job, we just need to decide if we
        #  can actually run the job.
        #

        #  If we're not running jobs on the grid, make sure the minimums are below what the machine has.

        if (getGlobal("useGrid") == 0) {
            if ((getPhysicalMemorySize() < $reqMemory) ||
                (getNumberOfCPUs()       < $reqThreads)) {
                print STDERR "--\n";
                print STDERR "-- WARNING: For genome size ", displayGenomeSize(getGlobal("genomeSize")), " I suggest minimum values:\n";
                print STDERR "-- WARNING:   minMemory=$minMemoryD\n";
                print STDERR "-- WARNING:   minThreads=$minThreads\n";
                print STDERR "-- WARNING:\n";
                print STDERR "-- WARNING: This machine has only:\n";
                print STDERR "-- WARNING:   maxMemory=", displayMemoryValue(getPhysicalMemorySize()), "\n";
                print STDERR "-- WARNING:   maxThreads=", getNumberOfCPUs(), "\n";
                print STDERR "-- WARNING:\n";
                print STDERR "-- WARNING: If the values are close, the assembly might be possible, but you need to set minMemory,\n";
                print STDERR "-- WARNING: minThreads, maxMemory and/or maxThreads manually to have:\n";
                print STDERR "-- WARNING:   minMemory  <= maxMemory\n";
                print STDERR "-- WARNING:   minThreads <= maxThreads\n";
                print STDERR "-- WARNING: (setting minMemory/minThreads is suggested; changing the max could over-commit resources)\n";
                print STDERR "--\n";
                caExit("machine limits exceeded", undef);
            }

            print STDERR "--\n";
            print STDERR "-- Local host has ", displayMemoryValue(getPhysicalMemorySize()), " memory and ", getNumberOfCPUs(), " CPUs, and can run jobs with ", displayMemoryValue($minMemory), " memory and $minThreads threads.\n";
        }

        #  Otherwise, we are running jobs on the grid, either under full grid control (useGrid=1) or manually (useGrid=remote).
        #  Make sure there is at least one host on the grid that can run the job.

        else {
            my @grid   = split '\0', getGlobal("availableHosts");
            my $nHosts = 0;
            my $tHosts = 0;

            foreach my $g (@grid) {
                my ($cpu, $mem, $num) = split '-', $g;

                $tHosts += $num;
                $nHosts += $num   if (($reqMemory <= $mem) && ($reqThreads <= $cpu));
            }

            if ($nHosts == 0) {
                print STDERR "--\n";
                print STDERR "-- WARNING: For genome size ", getGlobal("genomeSize"), " I suggest minimum values:\n";
                print STDERR "-- WARNING:   minMemory=$minMemoryD\n";
                print STDERR "-- WARNING:   minThreads=$minThreads\n";
                print STDERR "-- WARNING:\n";
                print STDERR "-- WARNING: There are no suitable hosts in the grid (listed above).\n";
                print STDERR "-- WARNING:\n";
                print STDERR "-- WARNING: If the values are close, the assembly might be possible, but you need to\n";
                print STDERR "-- WARNING: set minMemory and minThreads manually.\n";
                print STDERR "--\n";
                caExit("machine limits exceeded", undef);
            }

            print STDERR "--\n";
            print STDERR "-- Found $nHosts ", (($nHosts == 1) ? "host" : "hosts"), " that we can run jobs with ", displayMemoryValue($minMemory), " memory and $minThreads threads.\n";
        }
    }

    my $err;
    my $all;

    ($err, $all) = getAllowedResources("",    "bat",      $err, $all);
    ($err, $all) = getAllowedResources("cor", "mhap",     $err, $all);
    ($err, $all) = getAllowedResources("obt", "mhap",     $err, $all);
    ($err, $all) = getAllowedResources("utg", "mhap",     $err, $all);
    ($err, $all) = getAllowedResources("",    "red",      $err, $all);
    ($err, $all) = getAllowedResources("",    "oea",      $err, $all);
    ($err, $all) = getAllowedResources("",    "cns",      $err, $all);
    ($err, $all) = getAllowedResources("",    "ovlStore", $err, $all);
    ($err, $all) = getAllowedResources("",    "ovb",      $err, $all);
    ($err, $all) = getAllowedResources("",    "ovs",      $err, $all);
    ($err, $all) = getAllowedResources("cor", "ovl",      $err, $all);
    ($err, $all) = getAllowedResources("obt", "ovl",      $err, $all);
    ($err, $all) = getAllowedResources("utg", "ovl",      $err, $all);
    ($err, $all) = getAllowedResources("",    "meryl",    $err, $all);
    ($err, $all) = getAllowedResources("",    "cor",      $err, $all);

    print STDERR "--\n" if (defined($err));
    print STDERR $err   if (defined($err));
    print STDERR "--\n";
    print STDERR $all;
}


