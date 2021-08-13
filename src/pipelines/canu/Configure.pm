
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

package canu::Configure;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(getAllowedResources
             displayMemoryValue
             displayGenomeSize
             configureAssembler);

use strict;
use warnings "all";
no  warnings "uninitialized";
use Carp qw(cluck);
use Sys::Hostname;

use List::Util qw(min max);

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


sub findAbsoluteMax ($$) {
    my $type   = shift @_;
    my $class  = shift @_;

    #  Find the local resources.

    my $localmem  = getGlobal("localMemory");
    my $localcpu  = getGlobal("localThreads");

    #  Find the max memory available on any grid host.

    my $gridmem = 0;
    my $gridcpu = 0;
    my @grid    = split '\0', getGlobal("availableHosts");

    foreach my $g (@grid) {
        my ($cpu, $mem, $num) = split '-', $g;

        $gridmem = ($gridmem < $mem) ? $mem : $gridmem;
        $gridcpu = ($gridcpu < $cpu) ? $cpu : $gridcpu;
    }

    #  Threshold to any max limit.

    my $maxmem  = getGlobal("maxMemory");
    my $maxcpu  = getGlobal("maxThreads");

    if (defined($maxmem)) {
        $localmem = min($localmem, $maxmem);
        $gridmem = min($gridmem, $maxmem);
    }
    if (defined($maxcpu)) {
        $localcpu = min($localcpu, $maxcpu);
        $gridcpu = min($gridcpu, $maxcpu);
    }

    #  Return one or the other or fail.

    return(($class eq "grid") ? $gridmem : $localmem)   if ($type eq "memory");
    return(($class eq "grid") ? $gridcpu : $localcpu)   if ($type eq "threads");

    caExit("Invalid type '$type' supplied to findAbsoluteMax().", undef)
}


#  Side effect!  This will RESET the $global{} parameters to the computed value.  This lets
#  the rest of canu - in particular, the part that runs the jobs - use the correct value.  Without
#  resetting, I'd be making code changes all over the place to support the values returned.

sub getAllowedResources ($$$$$@);  #  Recursive call to getAllowedResources() wants the prototype.

sub getAllowedResources ($$$$$@) {
    my $tag  = shift @_;  #  Variant, e.g., "cor", "utg"
    my $alg  = shift @_;  #  Algorithm, e.g., "mhap", "ovl"
    my $err  = shift @_;  #  Report of things we can't run.
    my $all  = shift @_;  #  Report of things we can run.
    my $uni  = shift @_;  #  There's only one task to run (bogart)
    my $dbg  = shift @_;  #  Optional, report debugging stuff

    #  If no grid, or grid not enabled, everything falls under 'lcoal'.

    my $class = ((getGlobal("useGrid") ne "0") && (defined(getGlobal("gridEngine")))) ? "grid" : "local";

    #  If grid, but no hosts, fail.

    if (($class eq "grid") && (!defined(getGlobal("availableHosts")))) {
        caExit("invalid useGrid (" . getGlobal("useGrid") . ") and gridEngine (" . getGlobal("gridEngine") . "); found no execution hosts - is grid available from this host?", undef);
    }

    #  Figure out absolute limits of job sizes.  These are just limits on
    #  what values we'll iterate over when finding a job size.

    my $minMemory    = getGlobal("minMemory");
    my $minThreads   = getGlobal("minThreads");

    my $maxMemory    = findAbsoluteMax("memory",  $class);  #  This is the min of any user-defined max and the
    my $maxThreads   = findAbsoluteMax("threads", $class);  #  local/grid host maximum.  Not param maxMemory, etc.

    caExit("maxMemory is not defined", undef)   if (!defined($maxMemory));
    caExit("maxThreads is not defined", undef)  if (!defined($maxThreads));

    my $taskMemory   = getGlobal("${tag}${alg}Memory");   #  Algorithm limit, "utgovlMemory", etc.
    my $taskThreads  = getGlobal("${tag}${alg}Threads");  #

    caExit("${tag}${alg}Memory is not defined", undef)   if (!defined($taskMemory));
    caExit("${tag}${alg}Threads is not defined", undef)  if (!defined($taskThreads));

    if ($dbg) {
        my $minm  = getGlobal("minMemory");
        my $mint  = getGlobal("minThreads");

        my $maxm  = getGlobal("maxMemory");
        my $maxt  = getGlobal("maxThreads");

        undef $minm   if ($minm == 0);   #  The default is "0"; undef for easier testing below.
        undef $mint   if ($mint == 0);

        print STDERR "--\n"                              if (defined($minm) || defined($mint) || defined($maxm) || defined($maxt));
        print STDERR "-- ERROR  minMemory=${minm}g\n"    if (defined($minm));
        print STDERR "-- ERROR  maxMemory=${maxm}g\n"    if (defined($maxm));
        print STDERR "-- ERROR  minThreads=${mint}\n"    if (defined($mint));
        print STDERR "-- ERROR  maxThreads=${maxt}\n"    if (defined($maxt));
    }

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

    my $bestCores      = 0;
    my $bestMemory     = 16 * 1024 * 1024;   #  16 petabytes.
    my $bestCoresM     = undef;
    my $bestCoresT     = undef;
    my $availMemoryMin = undef;
    my $availMemoryMax = undef;

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

            for (my $ii=0; $ii<scalar(@gridCor); $ii++) {                                 #  Each process uses:
                my $np_cpu = $gridNum[$ii] * (($t == 0) ? 1 : int($gridCor[$ii] / $t));   #    $t cores, node has $gridCor[$ii] available.
                my $np_mem = $gridNum[$ii] * (($m == 0) ? 1 : int($gridMem[$ii] / $m));   #    $m GBs,   node ame $gridMem[$ii] available.

                my $np = ($np_cpu < $np_mem) ? $np_cpu : $np_mem;      #  Number of processes we can fit on this machine.

                if ($dbg) {
                    print STDERR "-- ERROR  for $t threads and $m memory - class$ii can support $np_cpu jobs(cores) and $np_mem jobs(memory), so $np jobs.\n";
                }

                #  If we only need to run one task, just remember if we can run the task on any machine here.

                if ($uni) {
                    if ($np > 0) {
                        $processes  = 1;
                        $cores      = $t;
                        $memory     = $m;
                    }
                }

                #  Otherwise, sum the number of processes we can run on the entire grid.

                else {
                    $processes += $np;        #  Total number of processes running
                    $cores     += $np * $t;   #  Total cores in use
                    $memory    += $np * $m;   #  Total memory in use
                }
            }

            if ($dbg) {
                print STDERR "-- ERROR  Tested $tag$alg requesting $t cores and ${m}GB memory and found $cores could be used.\n";
            }

            #  If no cores, then all machines were too small.

            next if ($cores == 0);

            #  Save the best one seen so far.  Break ties by selecting the one with the most memory.

            if (($bestCores <  $cores) ||
                ($bestCores <= $cores) && ($bestMemory > $memory)) {
                $bestCores  = $cores;
                $bestCoresT = $t;
                $bestCoresM = $m;
            }
        }
    }

    if (!defined($bestCoresM)) {
        getAllowedResources($tag, $alg, $err, $all, $uni, 1)  if (!defined($dbg));

        print STDERR "-- ERROR\n";
        print STDERR "-- ERROR  Task $tag$alg can't run on any available machines.\n";
        print STDERR "-- ERROR  It is requesting:\n";
        print STDERR "-- ERROR    ${tag}${alg}Memory=", getGlobal("${tag}${alg}Memory"), " gigabytes\n";
        print STDERR "-- ERROR    ${tag}${alg}Threads=", getGlobal("${tag}${alg}Threads"), " threads\n";
        print STDERR "-- ERROR\n";
        print STDERR "-- ERROR  No available machine configuration can run this task.\n";
        print STDERR "-- ERROR\n";
        print STDERR "-- ERROR  Possible solutions:\n";
        print STDERR "-- ERROR    Change maxMemory or maxThreads\n";
        print STDERR "-- ERROR    Change ${tag}${alg}Memory and/or ${tag}${alg}Threads\n";
        print STDERR "-- ERROR\n";

        caExit("task $tag$alg failed to find a configuration to run on", undef);
    }

    #  Reset the global values for later use.  SPECIAL CASE!  For ovsMemory, we just want the list
    #  of valid memory sizes.

    if ("$alg" eq "ovs") {
        $taskMemory  = $availMemoryMax;
        $taskThreads = $bestCoresT;

        setGlobal("${tag}${alg}Memory",  "$availMemoryMin-$availMemoryMax");
        setGlobal("${tag}${alg}Threads",  $taskThreads);

    } else {
        $taskMemory  = $bestCoresM;
        $taskThreads = $bestCoresT;

        setGlobal("${tag}${alg}Memory",  $taskMemory);
        setGlobal("${tag}${alg}Threads", $taskThreads);
    }

    #  Check for stupidity.

    caExit("invalid taskMemory=$taskMemory; maxMemory=$maxMemory", undef)     if ($taskMemory  > $maxMemory);
    caExit("invalid taskThread=$taskThreads; maxThreads=$maxThreads", undef)  if ($taskThreads > $maxThreads);

    #  Finally, reset the concurrency (if we're running locally) so we don't swamp our poor workstation.

    my $concurrent = undef;  #  Undef if in grid mode.

    if ($class eq "local") {
        my $nct = ($taskThreads == 0) ? 1 : int($maxThreads / $taskThreads);
        my $ncm = ($taskMemory  == 0) ? 1 : int($maxMemory  / $taskMemory);

        my $nc  = ($nct < $ncm) ? $nct : $ncm;

        $nc = 1  if (($uni) && ($nc > 0));

        #  If already set (on the command line), reset if too big.

        if (getGlobal("${tag}${alg}Concurrency") > $nc) {
            $err .= "-- Reset ${tag}${alg}Concurrency from " . getGlobal("${tag}${alg}Concurrency") . " to $nc.\n";
            setGlobal("${tag}${alg}Concurrency", $nc);
        }

        #  If not set, or set but zero, set it.

        if (!defined(getGlobal("${tag}${alg}Concurrency")) ||
            (getGlobal("${tag}${alg}Concurrency") == 0)) {
            setGlobal("${tag}${alg}Concurrency", $nc);
        }

        #  But if no memory set, set concurrency to zero.
        #  This is mostly to get the report correct; if set to undef, we think this is a grid job.

        if ($taskMemory == 0) {
            setGlobal("${tag}${alg}Concurrency", 0);
        }

        #  Update the local variable for the report.

        $concurrent = getGlobal("${tag}${alg}Concurrency");
    }

    #  And report.

    my $nam;

    if    ($alg eq "meryl")     {  $nam = "(k-mer counting)"; }
    elsif ($alg eq "hap")       {  $nam = "(read-to-haplotype assignment)"; }
    elsif ($alg eq "mhap")      {  $nam = "(overlap detection with mhap)"; }
    elsif ($alg eq "mmap")      {  $nam = "(overlap detection with minimap)"; }
    elsif ($alg eq "ovl")       {  $nam = "(overlap detection)"; }
    elsif ($alg eq "cor")       {  $nam = "(read correction)"; }
    elsif ($alg eq "ovb")       {  $nam = "(overlap store bucketizer)"; }
    elsif ($alg eq "ovs")       {  $nam = "(overlap store sorting)"; }
    elsif ($alg eq "red")       {  $nam = "(read error detection)"; }
    elsif ($alg eq "oea")       {  $nam = "(overlap error adjustment)"; }
    elsif ($alg eq "bat")       {  $nam = "(contig construction with bogart)"; }
    elsif ($alg eq "cns")       {  $nam = "(consensus)"; }
    else {
        caFailure("unknown task '$alg' in getAllowedResources().", undef);
    }

    my $mem  = sprintf("%7.3f GB",  $taskMemory);
    my $thr  = sprintf("%3d CPU%s", $taskThreads, (($taskThreads == 1) ? " " : "s"));
    my $job  = sprintf("%3d job%s", $concurrent,  (($concurrent  == 1) ? " " : "s"));

    my $memt = sprintf("%8.3f GB",  $taskMemory  * $concurrent);
    my $thrt = sprintf("%3d CPU%s", $taskThreads * $concurrent, (($taskThreads * $concurrent == 1) ? " " : "s"));

    $mem  = "  -.--- GB"    if ($taskMemory  == 0);
    $thr  = "  - CPUs"      if ($taskThreads == 0);
    $job  = "  - jobs"      if ($concurrent  == 0);

    $memt = "   -.--- GB"   if ($taskMemory  * $concurrent == 0);
    $thrt = "  - CPUs"      if ($taskThreads * $concurrent == 0);

    my $t = substr("$tag$alg     ", 0, 7);

    if (!defined($all)) {
        #$all .= "-- Memory, Threads and Concurrency configuration:\n"  if ( defined($concurrent));
        #$all .= "-- Memory and Threads configuration:\n"               if (!defined($concurrent));

        if (defined($concurrent)) {
            $all .= "--                                (tag)Concurrency\n";
            $all .= "--                         (tag)Threads          |\n";
            $all .= "--                (tag)Memory         |          |\n";
            $all .= "--        (tag)             |         |          |       total usage      algorithm\n";
            $all .= "--        -------  ----------  --------   --------  --------------------  -----------------------------\n";
        } else {
            $all .= "--                         (tag)Threads\n";
            $all .= "--                (tag)Memory         |\n";
            $all .= "--        (tag)             |         |  algorithm\n";
            $all .= "--        -------  ----------  --------  -----------------------------\n";
        }
    }
    $all .= "-- Local: $t  $mem  $thr x $job  $memt $thrt  $nam\n"      if ( defined($concurrent));
    $all .= "-- Grid:  $t  $mem  $thr  $nam\n"                          if (!defined($concurrent));

    return($err, $all);
}



#
#  If minMemory or minThreads isn't defined, pick a reasonable pair based on genome size.
#

sub configureAssembler () {

    #  For overlapper and mhap, allow larger maximums for larger genomes.  More memory won't help
    #  smaller genomes, and the smaller minimums won't hurt larger genomes (which are probably being
    #  run on larger machines anyway, so the minimums won't be used).

    #  For uncorrected overlapper, both memory and thread count is reduced.  Memory because it is
    #  very CPU bound, and thread count because it can be quite unbalanced.

    #  22 bits ->   864 MB table structure,   64 million kmers
    #  23 bits ->  1728 MB table structure,  128 million kmers
    #  24 bits ->  3456 MB table structure,  256 million kmers
    #  25 bits ->  6912 MB table structure,  512 million kmers
    #  26 bits -> 13824 MB table structure, 1024 million kmers
    #
    #    sequence generate -min 5000 -max 25000 -bases 10000000000                                   > random.fasta
    #    sequence generate -min 5000 -max 25000 -bases 10000000000 -a 0.9 -c 0.033 -g 0.033 -t 0.033 > repeat.fasta
    #
    #               TABLE    W/DATA
    #    bits 20   216 MB -  2500 MB random -   16 million kmers at 75% load
    #    bits 21   432 MB -  2750 MB random -   32 million kmers
    #
    #    bits 22   864 MB -  3000 MB random -   64 million kmers
    #
    #    bits 23  1728 MB -  3500 MB random -  128 million kmers
    #
    #    bits 24  3456 MB -  5750 MB random -  200 million kmers at 56% load
    #             3456 MB -  6500 MB random -  256 million kmers at 75% load
    #             3456 MB -  8200 MB repeat -  256 million kmers at  5% load
    #
    #    bits 25  6912 MB - 10000 MB random -  300 million kmers at 42% load
    #             6912 MB - 12750 MB random -  512 million kmers at 75% load
    #             6912 MB - 21000 MB repeat -  600 million kmers at  5% load
    #
    #    bits 26 13824 MB - 25000 MB random - 1024 million kmers at 75% load
    #
    #
    #  An expansion factor for the bases to load into the hash table.  Each table size will hold a
    #  little more than the above number of different kmers.  The additional third in the expansion
    #  is to adjust for repeated kmers in the reads.

    my $hx = 1.25 * 1000000;

    if      (getGlobal("genomeSize") < adjustGenomeSize("40m")) {
        setGlobalIfUndef("corOvlHashBlockLength",     2500000);    setGlobalIfUndef("obtOvlHashBlockLength",    64 * $hx);    setGlobalIfUndef("utgOvlHashBlockLength",    64 * $hx);
        setGlobalIfUndef("corOvlRefBlockLength",      2000000);    setGlobalIfUndef("obtOvlRefBlockLength",   1000000000);    setGlobalIfUndef("utgOvlRefBlockLength",   1000000000);   #    1 Gbp

        setGlobalIfUndef("corOvlMemory", "2");       setGlobalIfUndef("corOvlThreads", "1");      setGlobalIfUndef("corOvlHashBits", 22);
        setGlobalIfUndef("obtOvlMemory", "4");       setGlobalIfUndef("obtOvlThreads", "2-8");    setGlobalIfUndef("obtOvlHashBits", 22);
        setGlobalIfUndef("utgOvlMemory", "4");       setGlobalIfUndef("utgOvlThreads", "2-8");    setGlobalIfUndef("utgOvlHashBits", 22);

        setGlobalIfUndef("corMhapMemory", "4-6");    setGlobalIfUndef("corMhapThreads", "1-16");
        setGlobalIfUndef("obtMhapMemory", "4-6");    setGlobalIfUndef("obtMhapThreads", "1-16");
        setGlobalIfUndef("utgMhapMemory", "4-6");    setGlobalIfUndef("utgMhapThreads", "1-16");

        setGlobalIfUndef("corMMapMemory", "4-6");    setGlobalIfUndef("corMMapThreads", "1-16");
        setGlobalIfUndef("obtMMapMemory", "4-6");    setGlobalIfUndef("obtMMapThreads", "1-16");
        setGlobalIfUndef("utgMMapMemory", "4-6");    setGlobalIfUndef("utgMMapThreads", "1-16");

    } elsif (getGlobal("genomeSize") < adjustGenomeSize("500m")) {
        setGlobalIfUndef("corOvlHashBlockLength",     2500000);    setGlobalIfUndef("obtOvlHashBlockLength",   128 * $hx);    setGlobalIfUndef("utgOvlHashBlockLength",   128 * $hx);
        setGlobalIfUndef("corOvlRefBlockLength",      2000000);    setGlobalIfUndef("obtOvlRefBlockLength",   5000000000);    setGlobalIfUndef("utgOvlRefBlockLength",   5000000000);   #    5 Gbp

        setGlobalIfUndef("corOvlMemory", "2");       setGlobalIfUndef("corOvlThreads", "1");      setGlobalIfUndef("corOvlHashBits", 23);
        setGlobalIfUndef("obtOvlMemory", "8");       setGlobalIfUndef("obtOvlThreads", "2-8");    setGlobalIfUndef("obtOvlHashBits", 23);
        setGlobalIfUndef("utgOvlMemory", "8");       setGlobalIfUndef("utgOvlThreads", "2-8");    setGlobalIfUndef("utgOvlHashBits", 23);

        setGlobalIfUndef("corMhapMemory", "8-13");   setGlobalIfUndef("corMhapThreads", "1-16");
        setGlobalIfUndef("obtMhapMemory", "8-13");   setGlobalIfUndef("obtMhapThreads", "1-16");
        setGlobalIfUndef("utgMhapMemory", "8-13");   setGlobalIfUndef("utgMhapThreads", "1-16");

        setGlobalIfUndef("corMMapMemory", "8-13");   setGlobalIfUndef("corMMapThreads", "1-16");
        setGlobalIfUndef("obtMMapMemory", "8-13");   setGlobalIfUndef("obtMMapThreads", "1-16");
        setGlobalIfUndef("utgMMapMemory", "8-13");   setGlobalIfUndef("utgMMapThreads", "1-16");

    } elsif (getGlobal("genomeSize") < adjustGenomeSize("2g")) {
        setGlobalIfUndef("corOvlHashBlockLength",     2500000);    setGlobalIfUndef("obtOvlHashBlockLength",   256 * $hx);    setGlobalIfUndef("utgOvlHashBlockLength",   256 * $hx);
        setGlobalIfUndef("corOvlRefBlockLength",      2000000);    setGlobalIfUndef("obtOvlRefBlockLength",  15000000000);    setGlobalIfUndef("utgOvlRefBlockLength",  15000000000);   #   15 Gbp

        setGlobalIfUndef("corOvlMemory", "8");       setGlobalIfUndef("corOvlThreads", "1");      setGlobalIfUndef("corOvlHashBits", 24);
        setGlobalIfUndef("obtOvlMemory", "16");      setGlobalIfUndef("obtOvlThreads", "4-16");   setGlobalIfUndef("obtOvlHashBits", 24);
        setGlobalIfUndef("utgOvlMemory", "16");      setGlobalIfUndef("utgOvlThreads", "4-16");   setGlobalIfUndef("utgOvlHashBits", 24);

        setGlobalIfUndef("corMhapMemory", "16-32");  setGlobalIfUndef("corMhapThreads", "4-16");
        setGlobalIfUndef("obtMhapMemory", "16-32");  setGlobalIfUndef("obtMhapThreads", "4-16");
        setGlobalIfUndef("utgMhapMemory", "16-32");  setGlobalIfUndef("utgMhapThreads", "4-16");

        setGlobalIfUndef("corMMapMemory", "16-32");  setGlobalIfUndef("corMMapThreads", "1-16");
        setGlobalIfUndef("obtMMapMemory", "16-32");  setGlobalIfUndef("obtMMapThreads", "1-16");
        setGlobalIfUndef("utgMMapMemory", "16-32");  setGlobalIfUndef("utgMMapThreads", "1-16");

    } elsif (getGlobal("genomeSize") < adjustGenomeSize("5g")) {
        setGlobalIfUndef("corOvlHashBlockLength",     2500000);    setGlobalIfUndef("obtOvlHashBlockLength",   512 * $hx);    setGlobalIfUndef("utgOvlHashBlockLength",   512 * $hx);
        setGlobalIfUndef("corOvlRefBlockLength",      2000000);    setGlobalIfUndef("obtOvlRefBlockLength",  20000000000);    setGlobalIfUndef("utgOvlRefBlockLength",  20000000000);   #   20 Gbp

        setGlobalIfUndef("corOvlMemory", "8");       setGlobalIfUndef("corOvlThreads", "1");      setGlobalIfUndef("corOvlHashBits", 25);
        setGlobalIfUndef("obtOvlMemory", "24");      setGlobalIfUndef("obtOvlThreads", "4-16");   setGlobalIfUndef("obtOvlHashBits", 25);
        setGlobalIfUndef("utgOvlMemory", "24");      setGlobalIfUndef("utgOvlThreads", "4-16");   setGlobalIfUndef("utgOvlHashBits", 25);

        setGlobalIfUndef("corMhapMemory", "16-48");  setGlobalIfUndef("corMhapThreads", "4-16");
        setGlobalIfUndef("obtMhapMemory", "16-48");  setGlobalIfUndef("obtMhapThreads", "4-16");
        setGlobalIfUndef("utgMhapMemory", "16-48");  setGlobalIfUndef("utgMhapThreads", "4-16");

        setGlobalIfUndef("corMMapMemory", "16-48");  setGlobalIfUndef("corMMapThreads", "1-16");
        setGlobalIfUndef("obtMMapMemory", "16-48");  setGlobalIfUndef("obtMMapThreads", "1-16");
        setGlobalIfUndef("utgMMapMemory", "16-48");  setGlobalIfUndef("utgMMapThreads", "1-16");

    } else {
        setGlobalIfUndef("corOvlHashBlockLength",     2500000);    setGlobalIfUndef("obtOvlHashBlockLength",   512 * $hx);    setGlobalIfUndef("utgOvlHashBlockLength",   512 * $hx);
        setGlobalIfUndef("corOvlRefBlockLength",      2000000);    setGlobalIfUndef("obtOvlRefBlockLength",  30000000000);    setGlobalIfUndef("utgOvlRefBlockLength",  30000000000);   #   30 Gbp

        setGlobalIfUndef("corOvlMemory", "8");       setGlobalIfUndef("corOvlThreads", "1");      setGlobalIfUndef("corOvlHashBits", 25);
        setGlobalIfUndef("obtOvlMemory", "24");      setGlobalIfUndef("obtOvlThreads", "4-16");   setGlobalIfUndef("obtOvlHashBits", 25);
        setGlobalIfUndef("utgOvlMemory", "24");      setGlobalIfUndef("utgOvlThreads", "4-16");   setGlobalIfUndef("utgOvlHashBits", 25);

        setGlobalIfUndef("corMhapMemory", "32-64");  setGlobalIfUndef("corMhapThreads", "4-16");
        setGlobalIfUndef("obtMhapMemory", "32-64");  setGlobalIfUndef("obtMhapThreads", "4-16");
        setGlobalIfUndef("utgMhapMemory", "32-64");  setGlobalIfUndef("utgMhapThreads", "4-16");

        setGlobalIfUndef("corMMapMemory", "32-64");  setGlobalIfUndef("corMMapThreads", "1-16");
        setGlobalIfUndef("obtMMapMemory", "32-64");  setGlobalIfUndef("obtMMapThreads", "1-16");
        setGlobalIfUndef("utgMMapMemory", "32-64");  setGlobalIfUndef("utgMMapThreads", "1-16");
    }

    #  Overlap store construction should be based on the number of overlaps, but we obviously don't
    #  know that until much later.  If we set memory too large, we risk (in the parallel version for sure)
    #  inefficiency; too small and we run out of file handles.

    if      (getGlobal("genomeSize") < adjustGenomeSize("300m")) {
        setGlobalIfUndef("ovbMemory",   "4");       setGlobalIfUndef("ovbThreads",   "1");
        setGlobalIfUndef("ovsMemory",   "4-8");     setGlobalIfUndef("ovsThreads",   "1");

    } elsif (getGlobal("genomeSize") < adjustGenomeSize("1g")) {
        setGlobalIfUndef("ovbMemory",   "4");       setGlobalIfUndef("ovbThreads",   "1");
        setGlobalIfUndef("ovsMemory",   "8-16");    setGlobalIfUndef("ovsThreads",   "1");

    } else {
        setGlobalIfUndef("ovbMemory",   "4");       setGlobalIfUndef("ovbThreads",   "1");
        setGlobalIfUndef("ovsMemory",   "16-32");   setGlobalIfUndef("ovsThreads",   "1");
    }

    #  Correction and consensus are somewhat invariant.
    #    Correction memory is set based on read length in CorrectReads.pm.
    #    Consensus memory is set based on tig size in Consensus.pm.
    #  Both are set to zero here, a special case that will configure only the thread component.

    setGlobalIfUndef("cnsMemory",  "0");
    setGlobalIfUndef("corMemory",  "0");

    if      (getGlobal("genomeSize") < adjustGenomeSize("40m")) {
        setGlobalIfUndef("corThreads", "4");
        setGlobalIfUndef("cnsThreads", "1-4");
    } else {
        setGlobalIfUndef("corThreads", "4");
        setGlobalIfUndef("cnsThreads", "2-8");
    }

    if      (getGlobal("genomeSize") < adjustGenomeSize("40m")) {
        setGlobalIfUndef("corPartitions",   "64");
        setGlobalIfUndef("corPartitionMin", "10000");

    } elsif (getGlobal("genomeSize") < adjustGenomeSize("1g")) {
        setGlobalIfUndef("corPartitions",   "128");
        setGlobalIfUndef("corPartitionMin", "20000");

    } else {
        setGlobalIfUndef("corPartitions",   "256");
        setGlobalIfUndef("corPartitionMin", "40000");
    }

    #  Meryl too, basically just small or big.  This should really be using the number of bases
    #  reported from sqStore.

    if      (getGlobal("genomeSize") < adjustGenomeSize("100m")) {
        setGlobalIfUndef("merylMemory", "4-12");        setGlobalIfUndef("merylThreads", "1-4");

    } elsif (getGlobal("genomeSize") < adjustGenomeSize("1g")) {
        setGlobalIfUndef("merylMemory", "12-24");       setGlobalIfUndef("merylThreads", "1-8");

    } else {
        setGlobalIfUndef("merylMemory", "24-64");       setGlobalIfUndef("merylThreads", "1-8");
    }

    #  Total guesses on read-to-haplotype assignment.  A well-behaved nanopore human
    #  was running in 2GB memory.

    if      (getGlobal("genomeSize") < adjustGenomeSize("100m")) {
        setGlobalIfUndef("hapMemory", "4-8");     setGlobalIfUndef("hapThreads", "1-4");

    } elsif (getGlobal("genomeSize") < adjustGenomeSize("1g")) {
        setGlobalIfUndef("hapMemory", "6-12");    setGlobalIfUndef("hapThreads", "8-24");

    } else {
        setGlobalIfUndef("hapMemory", "8-16");    setGlobalIfUndef("hapThreads", "16-64");
    }


    if (getGlobal("genomeSize") < 10000000) {
        setGlobal("corOvlMerDistinct",  "0.9999")   if (!defined(getGlobal("corOvlMerThreshold")) && !defined(getGlobal("corOvlMerDistinct")));
        setGlobal("obtOvlMerDistinct",  "0.9999")   if (!defined(getGlobal("obtOvlMerThreshold")) && !defined(getGlobal("obtOvlMerDistinct")));
        setGlobal("utgOvlMerDistinct",  "0.9999")   if (!defined(getGlobal("utgOvlMerThreshold")) && !defined(getGlobal("utgOvlMerDistinct")));
    } else {
        setGlobal("corOvlMerDistinct",  "0.9990")   if (!defined(getGlobal("corOvlMerThreshold")) && !defined(getGlobal("corOvlMerDistinct")));
        setGlobal("obtOvlMerDistinct",  "0.9990")   if (!defined(getGlobal("obtOvlMerThreshold")) && !defined(getGlobal("obtOvlMerDistinct")));
        setGlobal("utgOvlMerDistinct",  "0.9990")   if (!defined(getGlobal("utgOvlMerThreshold")) && !defined(getGlobal("utgOvlMerDistinct")));
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
    setGlobalIfUndef("redBatchSize",   undef);
    setGlobalIfUndef("redBatchLength", "500000000");

    setGlobalIfUndef("oeaBatchSize",   undef);
    setGlobalIfUndef("oeaBatchLength", "300000000");

    if      (getGlobal("genomeSize") < adjustGenomeSize("40m")) {
        setGlobalIfUndef("redMemory", "8-16");        setGlobalIfUndef("redThreads", "2-4");
        setGlobalIfUndef("oeaMemory", "8");           setGlobalIfUndef("oeaThreads", "1");

    } elsif (getGlobal("genomeSize") < adjustGenomeSize("500m")) {
        setGlobalIfUndef("redMemory", "8-16");        setGlobalIfUndef("redThreads", "4-6");
        setGlobalIfUndef("oeaMemory", "8");           setGlobalIfUndef("oeaThreads", "1");

    } elsif (getGlobal("genomeSize") < adjustGenomeSize("2g")) {
        setGlobalIfUndef("redMemory", "16-32");       setGlobalIfUndef("redThreads", "4-8");
        setGlobalIfUndef("oeaMemory", "8");           setGlobalIfUndef("oeaThreads", "1");

    } elsif (getGlobal("genomeSize") < adjustGenomeSize("5g")) {
        setGlobalIfUndef("redMemory", "32-48");       setGlobalIfUndef("redThreads", "4-8");
        setGlobalIfUndef("oeaMemory", "8");           setGlobalIfUndef("oeaThreads", "1");

    } else {
        setGlobalIfUndef("redMemory", "32-64");       setGlobalIfUndef("redThreads", "6-10");
        setGlobalIfUndef("oeaMemory", "8");           setGlobalIfUndef("oeaThreads", "1");
    }

    #  And bogart.

    if      (getGlobal("genomeSize") < adjustGenomeSize("40m")) {
        setGlobalIfUndef("batMemory", "4-16");        setGlobalIfUndef("batThreads", "2-4");

    } elsif (getGlobal("genomeSize") < adjustGenomeSize("500m")) {
        setGlobalIfUndef("batMemory", "16-64");       setGlobalIfUndef("batThreads", "2-8");

    } elsif (getGlobal("genomeSize") < adjustGenomeSize("2g")) {
        setGlobalIfUndef("batMemory", "32-256");      setGlobalIfUndef("batThreads", "4-16");

    } elsif (getGlobal("genomeSize") < adjustGenomeSize("5g")) {
        setGlobalIfUndef("batMemory", "128-512");     setGlobalIfUndef("batThreads", "8-32");

    } else {
        setGlobalIfUndef("batMemory", "256-1024");    setGlobalIfUndef("batThreads", "16-64");
    }

    #  Log maxMemory setting.

    {
        my $maxcpu = getGlobal("maxThreads");            #  Get user-supplied limit.
        my $maxmem = getGlobal("maxMemory");

        if (($maxcpu > 0) ||
            ($maxmem > 0)) {
            printf STDERR "--\n";
            printf STDERR "-- Job limits:\n";
            printf(STDERR "--   %4d gigabytes memory  (maxMemory option).\n", $maxmem)    if (defined($maxmem));
            printf(STDERR "--   %4d CPUs              (maxThreads option).\n", $maxcpu)   if (defined($maxcpu));
        }
    }

    #  Finally, use all that setup to pick actual values for each component.
    #
    #  ovsMemory needs to be configured here iff the sequential build method is used.  This runs in
    #  the canu process, and needs to have a single memory size.  The parallel method will pick a
    #  memory size based on the number of overlaps and submit jobs using that size.

    my $err;
    my $all;

    ($err, $all) = getAllowedResources("",    "meryl",     $err, $all, 0);

    ($err, $all) = getAllowedResources("",    "hap",       $err, $all, 0);

    ($err, $all) = getAllowedResources("cor", "mhap",      $err, $all, 0)   if (getGlobal("corOverlapper") eq "mhap");
    ($err, $all) = getAllowedResources("cor", "mmap",      $err, $all, 0)   if (getGlobal("corOverlapper") eq "minimap");
    ($err, $all) = getAllowedResources("cor", "ovl",       $err, $all, 0)   if (getGlobal("corOverlapper") eq "ovl");

    ($err, $all) = getAllowedResources("obt", "mhap",      $err, $all, 0)   if (getGlobal("obtOverlapper") eq "mhap");
    ($err, $all) = getAllowedResources("obt", "mmap",      $err, $all, 0)   if (getGlobal("obtOverlapper") eq "minimap");
    ($err, $all) = getAllowedResources("obt", "ovl",       $err, $all, 0)   if (getGlobal("obtOverlapper") eq "ovl");

    ($err, $all) = getAllowedResources("utg", "mhap",      $err, $all, 0)   if (getGlobal("utgOverlapper") eq "mhap");
    ($err, $all) = getAllowedResources("utg", "mmap",      $err, $all, 0)   if (getGlobal("utgOverlapper") eq "minimap");
    ($err, $all) = getAllowedResources("utg", "ovl",       $err, $all, 0)   if (getGlobal("utgOverlapper") eq "ovl");

    ($err, $all) = getAllowedResources("",    "cor",       $err, $all, 0);

    ($err, $all) = getAllowedResources("",    "ovb",       $err, $all, 0);
    ($err, $all) = getAllowedResources("",    "ovs",       $err, $all, 0);

    ($err, $all) = getAllowedResources("",    "red",       $err, $all, 0);
    ($err, $all) = getAllowedResources("",    "oea",       $err, $all, 0);

    ($err, $all) = getAllowedResources("",    "bat",       $err, $all, 1)   if (uc(getGlobal("unitigger")) eq "BOGART");

    ($err, $all) = getAllowedResources("",    "cns",       $err, $all, 0);

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
