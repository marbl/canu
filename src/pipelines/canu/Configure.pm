
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
 #    Sergey Koren beginning on 2015-DEC-02
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

package canu::Configure;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(displayMemoryValue displayGenomeSize configureAssembler);

use strict;
use Carp qw(cluck);
use Sys::Hostname;

use canu::Defaults;
use canu::Execution;

#  This is called to expand parameter ranges for memory and thread parameters.
#  Examples of valid ranges:
#
#  no units  - 1-4:2                      - assumes 'g' in the adjust (if memory)
#  one unit  - 1g-4:2 1-4g:2 1-4:2g       - all the others are set to 'g'
#  two units - 1g-4g:2 1g-4:2g 1-4g:2g    - bgn/end are the same, stp uses end
#  all three - 1g-4g:2g                   - use as is
#
#  Quirks:  1g-2000m will increment every 1m.
#           1g-2000m:1g only adds 1g.
#           1g-2048m:1g adds 1 and 2g.

sub expandRange ($$$$) {
    my $var = shift @_;
    my $val = shift @_;
    my $min = shift @_;  #  limit the minimum to be above this
    my $max = shift @_;  #  limit the maximum to be below this

    my @v = split ',', $val;
    my @r;

    foreach my $v (@v) {
        my $bgn;  my $bgnu;
        my $end;  my $endu;
        my $stp;  my $stpu;

        #  Decode the range.

        if      ($v =~ m/^(\d+\.{0,1}\d*)([kKmMgGtT]{0,1})$/) {
            $bgn = $1;  $bgnu = $2;
            $end = $1;  $endu = $2;
            $stp =  1;  $stpu = $2;
        } elsif ($v =~ m/^(\d+\.{0,1}\d*)([kKmMgGtT]{0,1})-(\d+\.{0,1}\d*)([kKmMgGtT]{0,1})$/) {
            $bgn = $1;  $bgnu = $2;
            $end = $3;  $endu = $4;
            $stp =  1;  $stpu = $4;
        } elsif ($v =~ m/^(\d+\.{0,1}\d*)([kKmMgGtT]{0,1})-(\d+\.{0,1}\d*)([kKmMgGtT]{0,1}):(\d+\.{0,1}\d*)([kKmMgGtT]{0,1})$/) {
            $bgn = $1;  $bgnu = $2;
            $end = $3;  $endu = $4;
            $stp = $5;  $stpu = $6;
        } else {
            caExit("can't parse '$var' entry '$v'", undef);
        }

        #  Undef things that are null.  The code that follows this was written assuming undef.

        $bgnu = undef   if ($bgnu eq "");
        $endu = undef   if ($endu eq "");
        $stpu = undef   if ($stpu eq "");

        #  Process the range

        my $def = defined($bgnu) + defined($endu) + defined($stpu);

        #  If no units, this could be a memory or a thread setting.  Don't use units.
        if      ($def == 0) {
        }

        #  If only one unit specified, set the others to the same.
        elsif ($def == 1) {
            if    (defined($bgnu))  { $endu = $stpu = $bgnu;  }
            elsif (defined($endu))  { $bgnu = $stpu = $endu;  }
            elsif (defined($stpu))  { $bgnu = $endu = $stpu;  }
        }

        #  If two units specified, set the unset as:
        #    bgn or end unset - set based on the other range
        #    stp unset        - set on end if stp<end otherwise bgn
        elsif ($def == 2) {

            if ((!defined($bgnu) && ($endu ne $stpu)) ||
                (!defined($endu) && ($bgnu ne $stpu)) ||
                (!defined($stpu) && ($bgnu ne $endu))) {
                print STDERR "--\n";
                print STDERR "-- WARNING: incomplete and inconsistent units on '$var=$val'.\n";
            }

            $bgnu = $endu  if (!defined($bgnu));
            $endu = $bgnu  if (!defined($endu));
            $stpu = $endu  if (!defined($stpu) && ($stp <= $end));
            $stpu = $bgnu  if (!defined($stpu) && ($stp >  $end));
        }

        #  Nothing to do if all three are set!
        elsif ($def == 3) {
        }

        #  Convert the value and unit to gigabytes.

        my $b = adjustMemoryValue("$bgn$bgnu");
        my $e = adjustMemoryValue("$end$endu");
        my $s = adjustMemoryValue("$stp$stpu");

        #  Enforce the user supplied minimum and maximum.  We cannot 'decrease min to user supplied
        #  maximum' because this effectively ignores the task setting.  For, e.g., batMemory=64-128
        #  and maxMemory=32, we want it to fail.

        $b = $min   if ((defined($min)) && ($b < $min));    #  Increase min to user supplied minimum.
        $e = $min   if ((defined($min)) && ($e < $min));    #  Increase max to user supplied minimum.

        #$b = $max   if ((defined($max)) && ($b > $max));    #  Decrease min to use supplied maximum.
        $e = $max   if ((defined($max)) && ($e > $max));    #  Decrease max to use supplied maximum.

        #  Iterate over the range, push values to test onto the array.

        for (my $ii=$b; $ii<=$e; $ii += $s) {
            push @r, $ii;
        }
    }

    #print "$var = ";
    #foreach my $r (@r) {
    #    print "$r ";
    #}
    #print "\n";

    return(@r);
}


sub findGridMaxMemoryAndThreads () {
    my @grid   = split '\0', getGlobal("availableHosts");
    my $maxmem = 0;
    my $maxcpu = 0;

    foreach my $g (@grid) {
        my ($cpu, $mem, $num) = split '-', $g;

        $maxmem = ($maxmem < $mem) ? $mem : $maxmem;
        $maxcpu = ($maxcpu < $cpu) ? $cpu : $maxcpu;
    }

    return($maxmem, $maxcpu);
}


#  Side effect!  This will RESET the $global{} parameters to the computed value.  This lets
#  the rest of canu - in particular, the part that runs the jobs - use the correct value.  Without
#  resetting, I'd be making code changes all over the place to support the values returned.

sub getAllowedResources ($$$$@) {
    my $tag  = shift @_;  #  Variant, e.g., "cor", "utg"
    my $alg  = shift @_;  #  Algorithm, e.g., "mhap", "ovl"
    my $err  = shift @_;  #  Report of things we can't run.
    my $all  = shift @_;  #  Report of things we can run.
    my $dbg  = shift @_;  #  Optional, report debugging stuff

    #  If no grid, or grid not enabled, everything falls under 'lcoal'.

    my $class = ((getGlobal("useGrid") ne "0") && (defined(getGlobal("gridEngine")))) ? "grid" : "local";

    #  If grid, but no hosts, fail.

    if (($class eq "grid") && (!defined(getGlobal("availableHosts")))) {
        caExit("invalid useGrid (" . getGlobal("useGrid") . ") and gridEngine (" . getGlobal("gridEngine") . "); found no execution hosts - is grid available from this host?", undef);
    }

    #  Figure out limits.

    my $minMemory    = getGlobal("minMemory");
    my $minThreads   = getGlobal("minThreads");

    my $maxMemory    = getGlobal("maxMemory");
    my $maxThreads   = getGlobal("maxThreads");

    my $taskMemory   = getGlobal("${tag}${alg}Memory");   #  Algorithm limit, "utgovlMemory", etc.
    my $taskThreads  = getGlobal("${tag}${alg}Threads");  #

    #  The task limits MUST be defined.

    caExit("${tag}${alg}Memory is not defined", undef)   if (!defined($taskMemory));
    caExit("${tag}${alg}Threads is not defined", undef)  if (!defined($taskThreads));

    #  If the maximum limits aren't set, default to 'unlimited' (for the grid; we'll effectively filter
    #  by the number of jobs we can fit on the hosts) or to the current hardware limits.

    if ($dbg) {
        print STDERR "--\n";
        print STDERR "-- ERROR\n";
        print STDERR "-- ERROR  Limited to at least $minMemory GB memory via minMemory option\n"   if (defined($minMemory));
        print STDERR "-- ERROR  Limited to at least $minThreads threads via minThreads option\n"   if (defined($minThreads));
        print STDERR "-- ERROR  Limited to at most $maxMemory GB memory via maxMemory option\n"    if (defined($maxMemory));
        print STDERR "-- ERROR  Limited to at most $maxThreads threads via maxThreads option\n"    if (defined($maxThreads));
    }

    #  Figure out the largest memory and threads that could ever be supported.  This lets us short-circuit
    #  the loop below.

    my ($gridMaxMem, $gridMaxThr) = findGridMaxMemoryAndThreads();

    $maxMemory  = (($class eq "grid") ? $gridMaxMem : getPhysicalMemorySize())  if (!defined($maxMemory));
    $maxThreads = (($class eq "grid") ? $gridMaxThr : getNumberOfCPUs())        if (!defined($maxThreads));

    #  Build a list of the available hardware configurations we can run on.  If grid, we get this
    #  from the list of previously discovered hosts.  If local, it's just this machine.

    my @gridCor;  #  Number of cores
    my @gridMem;  #  GB's of memory
    my @gridNum;  #  Number of nodes

    if ($class eq "grid") {
        my @grid = split '\0', getGlobal("availableHosts");

        foreach my $g (@grid) {
            my ($cpu, $mem, $num) = split '-', $g;

            if (($cpu > 0) && ($mem > 0) && ($num > 0)) {
                push @gridCor, $cpu;
                push @gridMem, $mem;
                push @gridNum, $num;
            }
        }
    } else {
        push @gridCor, $maxThreads;
        push @gridMem, $maxMemory;
        push @gridNum, 1;
    }

    if ($dbg) {
        print STDERR "-- ERROR\n";
        print STDERR "-- ERROR  Found ", scalar(@gridCor), " machine ", ((scalar(@gridCor) == 1) ? "configuration:\n" : "configurations:\n");
        for (my $ii=0; $ii<scalar(@gridCor); $ii++) {
            print STDERR "-- ERROR    class$ii - $gridNum[$ii] machines with $gridCor[$ii] cores with $gridMem[$ii] GB memory each.\n";
        }
    }

    #  The task usually has multiple choices, and we have a little optimization problem to solve.  For each
    #  pair of memory/threads, compute three things:
    #    a) how many processes we can get running
    #    b) how many cores we can get running
    #    c) how much memory we can consume
    #  We then (typically) want to maximize the number of cores we can get running.
    #  Other options would be number of cores * amount of memory.

    my @taskMemory  = expandRange("${tag}${alg}Memory",  $taskMemory,  $minMemory,  $maxMemory);
    my @taskThreads = expandRange("${tag}${alg}Threads", $taskThreads, $minThreads, $maxThreads);

    #  Find task memory/thread settings that will maximize the number of cores running.  This used
    #  to also compute best as 'cores * memory' but that is better handled by ordering the task
    #  settings parameters.  The example below will pick the largest (last) configuration that
    #  maximizes core utilization:
    #
    #    taskThreads = 4,8,32,64
    #    taskMemory  = 16g,32g,64g

    my ($bestCores, $bestCoresM, $bestCoresT, $availMemoryMin, $availMemoryMax)  = (0, undef, undef, undef, undef);

    foreach my $m (@taskMemory) {
        foreach my $t (@taskThreads) {
            #if ($dbg && (($m > $maxMemory) || ($t > $maxThreads))) {
            #    print STDERR "-- ERROR Tested $tag$alg requesting $t cores and ${m}GB memory - rejected: limited to ${maxMemory}GB and $maxThreads cores.\n";
            #}
            next  if ($m > $maxMemory);   #  Bail if either of the suggest settings are
            next  if ($t > $maxThreads);  #  larger than the maximum allowed.

            #  Save this memory size.  ovsMemory uses a list of possible memory sizes to
            #  pick the smallest one that results in an acceptable number of files.

            $availMemoryMin = $m    if (!defined($availMemoryMin) || ($m < $availMemoryMin));
            $availMemoryMax = $m    if (!defined($availMemoryMax) || ($availMemoryMax < $m));

            #  For a job using $m GB memory and $t threads, we can compute how many processes will
            #  fit on each node in our set of available machines.  The smaller of the two is then
            #  the number of processes we can run on this node.

            my $processes = 0;
            my $cores     = 0;
            my $memory    = 0;

            for (my $ii=0; $ii<scalar(@gridCor); $ii++) {
                my $np_cpu = $gridNum[$ii] * int($gridCor[$ii] / $t);  #  Each process uses $t cores, node has $gridCor[$ii] cores available.
                my $np_mem = $gridNum[$ii] * int($gridMem[$ii] / $m);  #  Same idea.

                my $np = ($np_cpu < $np_mem) ? $np_cpu : $np_mem;

                if ($dbg) {
                    print STDERR "-- ERROR  for $t threads and $m memory - class$ii can support $np_cpu jobs(cores) and $np_mem jobs(memory), so $np jobs.\n";
                }

                $processes += $np;
                $cores     += $np * $t;
                $memory    += $np * $m;
            }

            if ($dbg) {
                print STDERR "-- ERROR  Tested $tag$alg requesting $t cores and ${m}GB memory and found $cores could be used.\n";
            }

            #  If no cores, then all machines were too small.

            next if ($cores == 0);

            #  Save the best one seen so far.

            if ($bestCores <= $cores) {
                $bestCores  = $cores;
                $bestCoresM = $m;
                $bestCoresT = $t;
            }
        }
    }

    if (!defined($bestCoresM)) {
        getAllowedResources($tag, $alg, $err, $all, 1)  if (!defined($dbg));

        print STDERR "-- ERROR\n";
        print STDERR "-- ERROR  Task $tag$alg can't run on any available machines.\n";
        print STDERR "-- ERROR  It is requesting:\n";
        print STDERR "-- ERROR    ${tag}${alg}Memory=", getGlobal("${tag}${alg}Memory"), " memory (gigabytes)\n";
        print STDERR "-- ERROR    ${tag}${alg}Threads=", getGlobal("${tag}${alg}Threads"), " threads\n";
        print STDERR "-- ERROR\n";
        print STDERR "-- ERROR  No available machine configuration can run this task.\n";
        print STDERR "-- ERROR\n";
        print STDERR "-- ERROR  Possible solutions:\n";
        print STDERR "-- ERROR    Increase maxMemory\n"  if (defined(getGlobal("maxMemory")));
        print STDERR "-- ERROR    Change ${tag}${alg}Memory and/or ${tag}${alg}Threads\n";
        print STDERR "-- ERROR\n";

        caExit("task $tag$alg failed to find a configuration to run on", undef);
    }

    #  Reset the global values for later use.  SPECIAL CASE!  For ovsMemory, we just want the list
    #  of valid memory sizes.

    if ("$alg" ne "ovs") {
        $taskMemory  = $bestCoresM;
        $taskThreads = $bestCoresT;

        setGlobal("${tag}${alg}Memory",  $taskMemory);
        setGlobal("${tag}${alg}Threads", $taskThreads);

    } else {
        $taskMemory  = $availMemoryMax;
        $taskThreads = $bestCoresT;

        setGlobal("${tag}${alg}Memory",  "$availMemoryMin-$availMemoryMax");
        setGlobal("${tag}${alg}Threads",  $taskThreads);
    }

    #  Check for stupidity.

    caExit("invalid taskMemory=$taskMemory; maxMemory=$maxMemory", undef)     if ($taskMemory > $maxMemory);
    caExit("invalid taskThread=$taskThreads; maxThreads=$maxThreads", undef)  if ($taskThreads > $maxThreads);

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

    if    ($alg eq "meryl")    {  $nam = "(k-mer counting)"; }
    elsif ($alg eq "mhap")     {  $nam = "(overlap detection with mhap)"; }
    elsif ($alg eq "mmap")     {  $nam = "(overlap detection with minimap)"; }
    elsif ($alg eq "ovl")      {  $nam = "(overlap detection)"; }
    elsif ($alg eq "cor")      {  $nam = "(read correction)"; }
    elsif ($alg eq "ovb")      {  $nam = "(overlap store bucketizer)"; }
    elsif ($alg eq "ovs")      {  $nam = "(overlap store sorting)"; }
    elsif ($alg eq "red")      {  $nam = "(read error detection)"; }
    elsif ($alg eq "oea")      {  $nam = "(overlap error adjustment)"; }
    elsif ($alg eq "bat")      {  $nam = "(contig construction)"; }
    elsif ($alg eq "cns")      {  $nam = "(consensus)"; }
    elsif ($alg eq "gfa")      {  $nam = "(GFA alignment and processing)"; }
    else {
        caFailure("unknown task '$alg' in getAllowedResources().", undef);
    }

    my $job = substr("    $concurrent",  -3) . " job" . (($concurrent == 1) ? " " : "s");
    my $thr = substr("    $taskThreads", -3) . " CPU" . (($taskThreads == 1) ? " " : "s");
    my $mem = substr("    $taskMemory",  -4) . " GB";

    my $t = substr("$tag$alg     ", 0, 7);

    if (!defined($all)) {
        #$all .= "-- Memory, Threads and Concurrency configuration:\n"  if ( defined($concurrent));
        #$all .= "-- Memory and Threads configuration:\n"               if (!defined($concurrent));

        if (defined($concurrent)) {
            $all .= "--                            (tag)Concurrency\n";
            $all .= "--                     (tag)Threads          |\n";
            $all .= "--            (tag)Memory         |          |\n";
            $all .= "--        (tag)         |         |          |  algorithm\n";
            $all .= "--        -------  ------  --------   --------  -----------------------------\n";
        } else {
            $all .= "--                     (tag)Threads\n";
            $all .= "--            (tag)Memory         |\n";
            $all .= "--        (tag)         |         |  algorithm\n";
            $all .= "--        -------  ------  --------  -----------------------------\n";
        }
    }
    $all .= "-- Local: $t $mem  $thr x $job  $nam\n"      if ( defined($concurrent));
    $all .= "-- Grid:  $t $mem  $thr  $nam\n"             if (!defined($concurrent));

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

    #  Parse units on things the user possibly set.

    setGlobal("genomeSize", adjustGenomeSize(getGlobal("genomeSize")));

    setGlobal("minMemory",  adjustMemoryValue(getGlobal("minMemory")));
    setGlobal("maxMemory",  adjustMemoryValue(getGlobal("maxMemory")));


    #  For overlapper and mhap, allow larger maximums for larger genomes.  More memory won't help
    #  smaller genomes, and the smaller minimums won't hurt larger genomes (which are probably being
    #  run on larger machines anyway, so the minimums won't be used).

    #  For uncorrected overlapper, both memory and thread count is reduced.  Memory because it is
    #  very CPU bound, and thread count because it can be quite unbalanced.

    if      (getGlobal("genomeSize") < adjustGenomeSize("40m")) {
        setGlobalIfUndef("corOvlMemory", "2-6");     setGlobalIfUndef("corOvlThreads", "1");
        setGlobalIfUndef("obtOvlMemory", "4-8");     setGlobalIfUndef("obtOvlThreads", "1-8");
        setGlobalIfUndef("utgOvlMemory", "4-8");     setGlobalIfUndef("utgOvlThreads", "1-8");

        setGlobalIfUndef("corMhapMemory", "4-6");   setGlobalIfUndef("corMhapThreads", "1-16");
        setGlobalIfUndef("obtMhapMemory", "4-6");   setGlobalIfUndef("obtMhapThreads", "1-16");
        setGlobalIfUndef("utgMhapMemory", "4-6");   setGlobalIfUndef("utgMhapThreads", "1-16");

        setGlobalIfUndef("corMMapMemory", "4-6");   setGlobalIfUndef("corMMapThreads", "1-16");
        setGlobalIfUndef("obtMMapMemory", "4-6");   setGlobalIfUndef("obtMMapThreads", "1-16");
        setGlobalIfUndef("utgMMapMemory", "4-6");   setGlobalIfUndef("utgMMapThreads", "1-16");


    } elsif (getGlobal("genomeSize") < adjustGenomeSize("500m")) {
        setGlobalIfUndef("corOvlMemory", "2-6");     setGlobalIfUndef("corOvlThreads", "1");
        setGlobalIfUndef("obtOvlMemory", "4-8");     setGlobalIfUndef("obtOvlThreads", "1-8");
        setGlobalIfUndef("utgOvlMemory", "4-8");     setGlobalIfUndef("utgOvlThreads", "1-8");

        setGlobalIfUndef("corMhapMemory", "8-13");   setGlobalIfUndef("corMhapThreads", "1-16");
        setGlobalIfUndef("obtMhapMemory", "8-13");   setGlobalIfUndef("obtMhapThreads", "1-16");
        setGlobalIfUndef("utgMhapMemory", "8-13");   setGlobalIfUndef("utgMhapThreads", "1-16");

        setGlobalIfUndef("corMMapMemory", "8-13");   setGlobalIfUndef("corMMapThreads", "1-16");
        setGlobalIfUndef("obtMMapMemory", "8-13");   setGlobalIfUndef("obtMMapThreads", "1-16");
        setGlobalIfUndef("utgMMapMemory", "8-13");   setGlobalIfUndef("utgMMapThreads", "1-16");

    } elsif (getGlobal("genomeSize") < adjustGenomeSize("2g")) {
        setGlobalIfUndef("corOvlMemory", "2-8");     setGlobalIfUndef("corOvlThreads", "1");
        setGlobalIfUndef("obtOvlMemory", "4-12");    setGlobalIfUndef("obtOvlThreads", "1-8");
        setGlobalIfUndef("utgOvlMemory", "4-12");    setGlobalIfUndef("utgOvlThreads", "1-8");

        setGlobalIfUndef("corMhapMemory", "16-32");   setGlobalIfUndef("corMhapThreads", "4-16");
        setGlobalIfUndef("obtMhapMemory", "16-32");   setGlobalIfUndef("obtMhapThreads", "4-16");
        setGlobalIfUndef("utgMhapMemory", "16-32");   setGlobalIfUndef("utgMhapThreads", "4-16");

        setGlobalIfUndef("corMMapMemory", "16-32");   setGlobalIfUndef("corMMapThreads", "1-16");
        setGlobalIfUndef("obtMMapMemory", "16-32");   setGlobalIfUndef("obtMMapThreads", "1-16");
        setGlobalIfUndef("utgMMapMemory", "16-32");   setGlobalIfUndef("utgMMapThreads", "1-16");

    } elsif (getGlobal("genomeSize") < adjustGenomeSize("5g")) {
        setGlobalIfUndef("corOvlMemory", "2-8");     setGlobalIfUndef("corOvlThreads", "1");
        setGlobalIfUndef("obtOvlMemory", "4-16");    setGlobalIfUndef("obtOvlThreads", "1-8");
        setGlobalIfUndef("utgOvlMemory", "4-16");    setGlobalIfUndef("utgOvlThreads", "1-8");

        setGlobalIfUndef("corMhapMemory", "16-48");  setGlobalIfUndef("corMhapThreads", "4-16");
        setGlobalIfUndef("obtMhapMemory", "16-48");  setGlobalIfUndef("obtMhapThreads", "4-16");
        setGlobalIfUndef("utgMhapMemory", "16-48");  setGlobalIfUndef("utgMhapThreads", "4-16");

        setGlobalIfUndef("corMMapMemory", "16-48");  setGlobalIfUndef("corMMapThreads", "1-16");
        setGlobalIfUndef("obtMMapMemory", "16-48");  setGlobalIfUndef("obtMMapThreads", "1-16");
        setGlobalIfUndef("utgMMapMemory", "16-48");  setGlobalIfUndef("utgMMapThreads", "1-16");

    } else {
        setGlobalIfUndef("corOvlMemory", "2-8");     setGlobalIfUndef("corOvlThreads", "1");
        setGlobalIfUndef("obtOvlMemory", "4-16");    setGlobalIfUndef("obtOvlThreads", "1-8");
        setGlobalIfUndef("utgOvlMemory", "4-16");    setGlobalIfUndef("utgOvlThreads", "1-8");

        setGlobalIfUndef("corMhapMemory", "32-64");  setGlobalIfUndef("corMhapThreads", "4-16");
        setGlobalIfUndef("obtMhapMemory", "32-64");  setGlobalIfUndef("obtMhapThreads", "4-16");
        setGlobalIfUndef("utgMhapMemory", "32-64");  setGlobalIfUndef("utgMhapThreads", "4-16");

        setGlobalIfUndef("corMMapMemory", "32-64");  setGlobalIfUndef("corMMapThreads", "1-16");
        setGlobalIfUndef("obtMMapMemory", "32-64");  setGlobalIfUndef("obtMMapThreads", "1-16");
        setGlobalIfUndef("utgMMapMemory", "32-64");  setGlobalIfUndef("utgMMapThreads", "1-16");
    }

    #  Overlapper block sizes probably don't need to be modified based on genome size.

    setGlobalIfUndef("corOvlHashBlockLength",   2500000);   setGlobalIfUndef("corOvlRefBlockSize",   20000);   setGlobalIfUndef("corOvlRefBlockLength", 0);
    setGlobalIfUndef("obtOvlHashBlockLength", 100000000);   setGlobalIfUndef("obtOvlRefBlockSize", 2000000);   setGlobalIfUndef("obtOvlRefBlockLength", 0);
    setGlobalIfUndef("utgOvlHashBlockLength", 100000000);   setGlobalIfUndef("utgOvlRefBlockSize", 2000000);   setGlobalIfUndef("utgOvlRefBlockLength", 0);

    #  Overlap store construction should be based on the number of overlaps, but we obviously don't
    #  know that until much later.  If we set memory too large, we risk (in the parallel version for sure)
    #  inefficiency; too small and we run out of file handles.

    if      (getGlobal("genomeSize") < adjustGenomeSize("300m")) {
        setGlobalIfUndef("ovsMethod", "sequential");
        setGlobalIfUndef("ovbMemory",   "2-4");     setGlobalIfUndef("ovbThreads",   "1");
        setGlobalIfUndef("ovsMemory",   "2-8");     setGlobalIfUndef("ovsThreads",   "1");

    } elsif (getGlobal("genomeSize") < adjustGenomeSize("1g")) {
        setGlobalIfUndef("ovsMethod", "parallel");
        setGlobalIfUndef("ovbMemory",   "2-4");     setGlobalIfUndef("ovbThreads",   "1");
        setGlobalIfUndef("ovsMemory",   "4-16");    setGlobalIfUndef("ovsThreads",   "1");

    } else {
        setGlobalIfUndef("ovsMethod", "parallel");
        setGlobalIfUndef("ovbMemory",   "2-4");     setGlobalIfUndef("ovbThreads",   "1");
        setGlobalIfUndef("ovsMemory",   "4-32");    setGlobalIfUndef("ovsThreads",   "1");
    }

    #  Correction and consensus are somewhat invariant.  Correction memory is set based on read length
    #  in CorrectReads.pm.

    if      (getGlobal("genomeSize") < adjustGenomeSize("40m")) {
        setGlobalIfUndef("cnsMemory",     "8-32");     setGlobalIfUndef("cnsThreads",      "1-4");
        setGlobalIfUndef("corMemory",     undef);      setGlobalIfUndef("corThreads",      "4");
        setGlobalIfUndef("cnsPartitions", "8");        setGlobalIfUndef("cnsPartitionMin", "15000");
        setGlobalIfUndef("corPartitions", "256");      setGlobalIfUndef("corPartitionMin", "5000");

    } elsif (getGlobal("genomeSize") < adjustGenomeSize("1g")) {
        setGlobalIfUndef("cnsMemory",     "16-48");    setGlobalIfUndef("cnsThreads",      "2-8");
        setGlobalIfUndef("corMemory",     undef);      setGlobalIfUndef("corThreads",      "4");
        setGlobalIfUndef("cnsPartitions", "64");       setGlobalIfUndef("cnsPartitionMin", "20000");
        setGlobalIfUndef("corPartitions", "512");      setGlobalIfUndef("corPartitionMin", "10000");

    } else {
        setGlobalIfUndef("cnsMemory",     "64-128");   setGlobalIfUndef("cnsThreads",      "2-8");
        setGlobalIfUndef("corMemory",     undef);      setGlobalIfUndef("corThreads",      "4");
        setGlobalIfUndef("cnsPartitions", "256");      setGlobalIfUndef("cnsPartitionMin", "25000");
        setGlobalIfUndef("corPartitions", "1024");     setGlobalIfUndef("corPartitionMin", "15000");
    }

    #  Meryl too, basically just small or big.  This should really be using the number of bases
    #  reported from gatekeeper.

    if      (getGlobal("genomeSize") < adjustGenomeSize("100m")) {
        setGlobalIfUndef("merylMemory", "4-8");     setGlobalIfUndef("merylThreads", "1-4");

    } elsif (getGlobal("genomeSize") < adjustGenomeSize("1g")) {
        setGlobalIfUndef("merylMemory", "16-64");    setGlobalIfUndef("merylThreads", "1-16");

    } else {
        setGlobalIfUndef("merylMemory", "64-256");   setGlobalIfUndef("merylThreads", "1-32");
    }

    #  Overlap error adjustment
    #
    #  Configuration is primarily done though memory size.  This blows up when there are many
    #  short(er) reads and large memory machines are available.
    #
    #  The limit is arbitrary.
    #   On medicago,   with 740,000 reads (median len  ~1,500bp), this will result in about 150 jobs.
    #   The memory-only limit generated only 7 jobs.
    #
    #   On drosophila, with 270,000 reads (median len ~17,000bp), this will result in about  50 jobs.
    #   The memory-only limit generated 36 jobs.
    #
    setGlobalIfUndef("redBatchSize",   "5000");
    setGlobalIfUndef("redBatchLength", "");

    setGlobalIfUndef("oeaBatchSize",   "25000");
    setGlobalIfUndef("oeaBatchLength", "");

    if      (getGlobal("genomeSize") < adjustGenomeSize("40m")) {
        setGlobalIfUndef("redMemory",   "1-2");    setGlobalIfUndef("redThreads",   "1-4");
        setGlobalIfUndef("oeaMemory",   "1");      setGlobalIfUndef("oeaThreads",   "1");

    } elsif (getGlobal("genomeSize") < adjustGenomeSize("500m")) {
        setGlobalIfUndef("redMemory",   "2-6");    setGlobalIfUndef("redThreads",   "1-6");
        setGlobalIfUndef("oeaMemory",   "2");       setGlobalIfUndef("oeaThreads",   "1");

    } elsif (getGlobal("genomeSize") < adjustGenomeSize("2g")) {
        setGlobalIfUndef("redMemory",   "2-8");    setGlobalIfUndef("redThreads",   "1-8");
        setGlobalIfUndef("oeaMemory",   "2");       setGlobalIfUndef("oeaThreads",   "1");

    } elsif (getGlobal("genomeSize") < adjustGenomeSize("5g")) {
        setGlobalIfUndef("redMemory",   "2-16");    setGlobalIfUndef("redThreads",   "1-8");
        setGlobalIfUndef("oeaMemory",   "4");       setGlobalIfUndef("oeaThreads",   "1");

    } else {
        setGlobalIfUndef("redMemory",   "2-16");    setGlobalIfUndef("redThreads",   "1-8");
        setGlobalIfUndef("oeaMemory",   "4");       setGlobalIfUndef("oeaThreads",   "1");
    }

    #  And bogart and GFA alignment/processing.
    #
    #  GFA for genomes less than 40m is run in the canu process itself.

    if      (getGlobal("genomeSize") < adjustGenomeSize("40m")) {
        setGlobalIfUndef("batMemory",   "2-16");        setGlobalIfUndef("batThreads",   "1-4");
        setGlobalIfUndef("gfaMemory",   "2-8");         setGlobalIfUndef("gfaThreads",   "1-4");

    } elsif (getGlobal("genomeSize") < adjustGenomeSize("500m")) {
        setGlobalIfUndef("batMemory",   "16-64");       setGlobalIfUndef("batThreads",   "2-8");
        setGlobalIfUndef("gfaMemory",   "4-8");         setGlobalIfUndef("gfaThreads",   "2-8");

    } elsif (getGlobal("genomeSize") < adjustGenomeSize("2g")) {
        setGlobalIfUndef("batMemory",   "32-256");      setGlobalIfUndef("batThreads",   "4-16");
        setGlobalIfUndef("gfaMemory",   "8-16");         setGlobalIfUndef("gfaThreads",  "4-16");

    } elsif (getGlobal("genomeSize") < adjustGenomeSize("5g")) {
        setGlobalIfUndef("batMemory",   "128-512");     setGlobalIfUndef("batThreads",   "8-32");
        setGlobalIfUndef("gfaMemory",   "16-32");        setGlobalIfUndef("gfaThreads",  "8-32");

    } else {
        setGlobalIfUndef("batMemory",   "256-1024");    setGlobalIfUndef("batThreads",   "16-64");
        setGlobalIfUndef("gfaMemory",   "32-64");       setGlobalIfUndef("gfaThreads",   "16-64");
    }




    #  Finally, use all that setup to pick actual values for each component.
    #
    #  ovsMemory needs to be configured here iff the sequential build method is used.  This runs in
    #  the canu process, and needs to have a single memory size.  The parallel method will pick a
    #  memory size based on the number of overlaps and submit jobs using that size.

    my $err;
    my $all;

    ($err, $all) = getAllowedResources("",    "meryl",    $err, $all);

    ($err, $all) = getAllowedResources("cor", "mhap",     $err, $all)   if (getGlobal("corOverlapper") eq "mhap");
    ($err, $all) = getAllowedResources("cor", "mmap",     $err, $all)   if (getGlobal("corOverlapper") eq "minimap");
    ($err, $all) = getAllowedResources("cor", "ovl",      $err, $all)   if (getGlobal("corOverlapper") eq "ovl");

    ($err, $all) = getAllowedResources("obt", "mhap",     $err, $all)   if (getGlobal("obtOverlapper") eq "mhap");
    ($err, $all) = getAllowedResources("obt", "mmap",     $err, $all)   if (getGlobal("obtOverlapper") eq "minimap");
    ($err, $all) = getAllowedResources("obt", "ovl",      $err, $all)   if (getGlobal("obtOverlapper") eq "ovl");

    ($err, $all) = getAllowedResources("utg", "mhap",     $err, $all)   if (getGlobal("utgOverlapper") eq "mhap");
    ($err, $all) = getAllowedResources("utg", "mmap",     $err, $all)   if (getGlobal("utgOverlapper") eq "minimap");
    ($err, $all) = getAllowedResources("utg", "ovl",      $err, $all)   if (getGlobal("utgOverlapper") eq "ovl");

    #  Usually set based on read length in CorrectReads.pm.  If defined by the user, run through configuration.
    ($err, $all) = getAllowedResources("",    "cor",      $err, $all)   if (getGlobal("corMemory") ne undef);

    ($err, $all) = getAllowedResources("",    "ovb",      $err, $all);
    ($err, $all) = getAllowedResources("",    "ovs",      $err, $all);

    ($err, $all) = getAllowedResources("",    "red",      $err, $all);
    ($err, $all) = getAllowedResources("",    "oea",      $err, $all);

    ($err, $all) = getAllowedResources("",    "bat",      $err, $all);

    ($err, $all) = getAllowedResources("",    "cns",      $err, $all);

    ($err, $all) = getAllowedResources("",    "gfa",      $err, $all);

    #  Check some minimums.

    if ((getGlobal("ovsMemory") =~ m/^([0123456789.]+)-*[0123456789.]*$/) &&
        ($1 < 0.25)) {
        caExit("ovsMemory must be at least 0.25g or 256m", undef);
    }

    #  2017-02-21 -- not sure why $err is being reported here if it doesn't stop.  What's in it?

    print STDERR "--\n" if (defined($err));
    print STDERR $err   if (defined($err));
    print STDERR "--\n";
    print STDERR $all;
}


