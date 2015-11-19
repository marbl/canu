
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

package canu::Defaults;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(getCommandLineOptions addCommandLineOption writeLog caExit caFailure getNumberOfCPUs getPhysicalMemorySize getAllowedResources diskSpace printHelp setParametersFromFile setParametersFromCommandLine checkParameters getGlobal setGlobal showErrorRates setErrorRate setDefaults);

use strict;
use Carp qw(cluck);
use Sys::Hostname;

use Filesys::Df;  #  for diskSpace()

my %global;
my %synops;
my %synnam;

my $cLineOpts = "";
my $specLog   = "";



#  Return the second argument, unless the first argument is found in
#  %global, in which case return that.
#
sub getGlobal ($) {
    my $var = shift @_;

    $var =~ tr/A-Z/a-z/;

    caFailure("parameter '$var' is not known", undef) if (!exists($global{$var}));

    return($global{$var});
}


sub globalExists ($) {
    my $var = shift @_;

    $var =~ tr/A-Z/a-z/;

    return(exists($global{$var}));
}


sub setGlobalSpecialization ($@) {
    my $val = shift @_;

    foreach my $var (@_) {
        #print STDERR "set specialization $var = $val\n";
        $global{$var} = $val;
    }

    return(1);
}

sub setGlobal ($$) {
    my $var = shift @_;
    my $val = shift @_;
    my $set = 0;

    $var =~ tr/A-Z/a-z/;
    $val = undef  if ($val eq "");  #  Set to undefined, the default for many of the options.

    #  Map 'true'/'false' et al. to 0/1.

    $val = 0  if (($val =~ m/^false$/i) || ($val =~ m/^f$/i));
    $val = 1  if (($val =~ m/^true$/i)  || ($val =~ m/^t$/i));

    #  Translate from generic to specialized var

    foreach my $alg ("ovl", "mhap") {
        foreach my $opt ("usegrid", "gridoptions") {
            $set += setGlobalSpecialization($val, ("${opt}cor${alg}", "${opt}obt${alg}", "${opt}utg${alg}"))  if ($var eq "${opt}${alg}");
        }

        foreach my $opt ("memory", "threads", "concurrency") {
            $set += setGlobalSpecialization($val, ("cor${alg}${opt}", "obt${alg}${opt}", "utg${alg}${opt}"))  if ($var eq "${alg}${opt}");
        }
    }

    foreach my $opt ("overlapper") {
        $set += setGlobalSpecialization($val, ("cor${opt}", "obt${opt}", "utg${opt}"))  if ($var eq "${opt}");
    }

    #  e.g., corOvlHashBlockLength
    foreach my $opt ("ovlerrorrate", "ovlhashblocklength", "ovlrefblocksize", "ovlrefblocklength", "ovlhashbits", "ovlhashload", "ovlmersize", "ovlmerthreshold", "ovlmerdistinct", "ovlmertotal", "ovlfrequentmers") {
        $set += setGlobalSpecialization($val, ("cor${opt}", "obt${opt}", "utg${opt}"))  if ($var eq "${opt}");
    }

    #  e.g., corMhapBlockSize
    foreach my $opt ("mhapblocksize", "mhapmersize", "mhaprealign", "mhapsensitivity") {
        $set += setGlobalSpecialization($val, ("cor${opt}", "obt${opt}", "utg${opt}"))  if ($var eq "${opt}");
    }

    return  if ($set > 0);

    if ($var eq "errorrate") {
        setErrorRate($val);
        return;
    }

    caFailure("paramter '$var' is not known", undef) if (!exists($global{$var}));

    $global{$var} = $val;
}



sub setGlobalIfUndef ($$) {
    my $var = shift @_;
    my $val = shift @_;

    $var =~ tr/A-Z/a-z/;
    $val = undef  if ($val eq "");  #  Set to undefined, the default for many of the options.

    return  if (defined($global{$var}));

    $global{$var} = $val;
}



sub getCommandLineOptions () {
    return($cLineOpts);
}



sub addCommandLineOption ($) {
    if ($cLineOpts =~ m/\s$/) {
        $cLineOpts .= "$_[0]";
    } else {
        $cLineOpts .= " $_[0]";
    }
}



sub writeLog ($) {
    my $wrk = shift @_;

    my $time = time();
    my $host = hostname();
    my $pid  = $$;

    open(F, "> $wrk/canu-logs/${time}_${host}_${pid}_canu");
    print F $specLog;
    close(F);
}



#  Use caExit() for transient errors, like not opening files, processes that die, etc.
sub caExit ($$) {
    my  $msg   = shift @_;
    my  $log   = shift @_;

    print STDERR "================================================================================\n";
    print STDERR "Don't panic, but a mostly harmless error occurred and canu failed.\n";
    print STDERR "\n";

    #  Really should pass in $wrk
    if (defined($log)) {
        my  $df = diskSpace($log);

        print STDERR "Disk space available:  $df GB\n";
        print STDERR "\n";
    }

    if (-e $log) {
        print STDERR "Last 50 lines of the relevant log file ($log):\n";
        print STDERR "\n";
        system("tail -n 50 $log");
        print STDERR "\n";
    }

    print STDERR "canu failed with '$msg'.\n";
    print STDERR "\n";

    exit(1);
}


#  Use caFailure() for errors that definitely will require code changes to fix.
sub caFailure ($$) {
    my  $msg   = shift @_;
    my  $log   = shift @_;

    print STDERR "================================================================================\n";
    print STDERR "Please panic.  canu failed, and it shouldn't have.\n";
    print STDERR "\n";
    print STDERR "Stack trace:\n";
    print STDERR "\n";
    cluck;
    print STDERR "\n";

    if (-e $log) {
        print STDERR "Last few lines of the relevant log file ($log):\n";
        print STDERR "\n";
        system("tail -n 50 $log");
    }

    print STDERR "\n";
    print STDERR "canu failed with '$msg'.\n";

    exit(1);
}


#
#  Host management - these really belong in 'Execution.pm' (or 'Utilities.pm') but can't go there
#  (Execution.pm) and be used here too.
#

sub getNumberOfCPUs () {
    my $os   = $^O;
    my $ncpu = 1;

    #  See http://stackoverflow.com/questions/6481005/obtain-the-number-of-cpus-cores-in-linux

    if ($os eq "freebsd") {
        $ncpu = int(`/sbin/sysctl -n hw.ncpu`);
    }

    if ($os eq "darwin") {
        $ncpu = int(`/usr/bin/getconf _NPROCESSORS_ONLN`);
    }

    if ($os eq "linux") {
        $ncpu = int(`getconf _NPROCESSORS_ONLN`);
    }

    return($ncpu);
}


sub getPhysicalMemorySize () {
    my $os     = $^O;
    my $memory = 1;

    if ($os eq "freebsd") {
        $memory = `/sbin/sysctl -n hw.physmem` / 1024 / 1024 / 1024;
    }

    if ($os eq "darwin") {
        $memory = `/usr/sbin/sysctl -n hw.memsize` / 1024 / 1024 / 1024;
    }

    if ($os eq "linux") {
        open(F, "< /proc/meminfo");        #  Way to go, Linux!  Make it easy on us!
        while (<F>) {
            if (m/MemTotal:\s+(\d+)/) {
                $memory = $1 / 1024 / 1024;
            }
        }
        close(F);
    }

    return(int($memory + 0.5));  #  Poor man's rounding
}



#  Side effect!  This will RESET the $global{} parameters to the computed value.  This lets
#  the rest of canu - in particular, the part that runs the jobs - use the correct value.  Without
#  resetting, I'd be making code changes all over the place to support the values returned.

sub getAllowedResources ($$$$) {
    my $tag  = shift @_;  #  Variant, e.g., "cor", "utg"
    my $alg  = shift @_;  #  Algorithm, e.g., "mhap", "ovl"
    my $cls  = undef;     #  Class, "grid" or "master"
    my $nam  = undef;     #  Human readable name
    my $err  = shift @_;  #  Report of things we can't run.
    my $all  = shift @_;  #  Report of things we can run.

    #  Decide if the task is 'grid' or 'master'.

    if    ($alg eq "master")   {  $cls = "master";  $nam = "sequential stages"; }
    elsif ($alg eq "bat")      {  $cls = "master";  $nam = "bogart (unitigger)"; }
    elsif ($alg eq "cns")      {  $cls = "grid";    $nam = "utgcns (consensus"; }
    elsif ($alg eq "cor")      {  $cls = "grid";    $nam = "falcon_sense (read correction)"; }
    elsif ($alg eq "meryl")    {  $cls = "master";  $nam = "meryl (k-mer counting)"; }
    elsif ($alg eq "oea")      {  $cls = "grid";    $nam = "overlap error adjustment"; }
    elsif ($alg eq "ovb")      {  $cls = "grid";    $nam = "overlap store parallel bucketizer"; }
    elsif ($alg eq "ovlStore") {  $cls = "master";  $nam = "overlap store sequential building"; }
    elsif ($alg eq "ovs")      {  $cls = "grid";    $nam = "overlap store parallel sorting"; }
    elsif ($alg eq "red")      {  $cls = "grid";    $nam = "read error detection (overlap error adjustment)"; }
    elsif ($alg eq "mhap")     {  $cls = "grid";    $nam = "mhap (overlapper)"; }
    elsif ($alg eq "ovl")      {  $cls = "grid";    $nam = "overlapper"; }
    else {
        caFailure("unknown task '$alg' in getAllowedResources().", undef);
    }

    #  If no grid, or grid not enabled, everything falls under 'master'.

    $cls = "master"  if ((getGlobal("useGrid") == 0) || (getGlobal("gridEngine") eq undef));

    #  Figure out limits.

    my $hostMemory   = getGlobal("${cls}Memory");       #  Host limit, "gridMemory", "masterThreads", etc.
    my $hostThreads  = getGlobal("${cls}Threads");      #

    my $taskMemory   = getGlobal("${tag}${alg}Memory");   #  Algorithm limit, "utgovlMemory", etc.
    my $taskThreads  = getGlobal("${tag}${alg}Threads");  #

    #  If the host limits aren't set, default to 'unlimited' (for the grid; we'll effectively filter
    #  by the number of jobs we can fit on the hosts) or to the current hardware limits.

    $hostMemory  = (($cls eq "grid") ? 1024 * 1024 : getPhysicalMemorySize())  if (!defined($hostMemory));    #  1 PB memory!
    $hostThreads = (($cls eq "grid") ? 1024        : getNumberOfCPUs())        if (!defined($hostThreads));   #  1 k  cores!

    #  If the two task parameters is undefined, they are left up to the algorithm.  This causes
    #  problems here.  They should only be undefined for 'master' jobs (bogart and meryl).  We reset
    #  to the limits.  To be safe, ALL types are reset.

    $taskMemory  = $hostMemory    if (!defined($taskMemory));
    $taskThreads = $hostThreads   if (!defined($taskThreads));

    #  Build a list of the available hardware configurations we can run on.

    my @gridCor;  #  Number of cores
    my @gridMem;  #  GB's of memory
    my @gridNum;  #  Number of nodes

    if ($cls eq "grid") {
        my @grid = split '\0', getGlobal("availableHosts");

        foreach my $g (@grid) {
            my ($cpu, $mem, $num) = split '-', $g;

            push @gridCor, $cpu;
            push @gridMem, $mem;
            push @gridNum, $num;
        }
    } else {
        push @gridCor, $hostThreads;
        push @gridMem, $hostMemory;
        push @gridNum, 1;
    }

    #  If the task has multiple choices, we have a little optimization problem to solve.  Foreach
    #  pair of memory/threads, compute three things:
    #    a) how many processes we can get running
    #    b) how many cores we can get running
    #    c) how much memory we can consume
    #  We then (typically) want to maximize the number of cores we can get running.
    #  Other options would be number of cores * amount of memory.

    #  Filter out task settings that can't be run based on the gridMemory/gridThreads or masterMemory/masterThreads setting.
    #  (actually, this just reports those that would be filtered; the actual filtering is inline in the algorithm)

    my @taskMemory  = split ',', $taskMemory;
    my @taskThreads = split ',', $taskThreads;

    my $ignoreM;
    my $ignoreT;

    foreach my $m (@taskMemory) {
        $m = adjustMemoryValue($m);
    }

    foreach my $m (@taskMemory) {
        next  if ($m <= $hostMemory);
        $ignoreM .= ","  if (defined($ignoreM));
        $ignoreM .= "${m}g";
    }
    foreach my $t (@taskThreads) {
        next  if ($t <= $hostThreads);
        $ignoreT .= ","  if (defined($ignoreT));
        $ignoreT .= "$t";
    }

    if      (defined($ignoreM) && defined($ignoreT)) {
        $err .= "-- Can't use ${tag}${alg}Memory=$ignoreM and ${tag}${alg}Threads=$ignoreT because of ${cls}Memory=" . getGlobal("${cls}Memory") . "g and ${cls}Threads=" . getGlobal("${cls}Threads") . " limits.\n";

    } elsif (defined($ignoreM)) {
        $err .= "-- Can't use ${tag}${alg}Memory=$ignoreM because of ${cls}Memory=" . getGlobal("${cls}Memory") . "g limit.\n";

    } elsif (defined($ignoreT)) {
        $err .= "-- Can't use ${tag}${alg}Threads=$ignoreT because of ${cls}Threads=" . getGlobal("${cls}Threads") . " limit.\n";
    }

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

            next  if ($m > $hostMemory);   #  Bail if either of the suggest settings are
            next  if ($t > $hostThreads);  #  larger than the maximum allowed.

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

    $taskThreads = $hostThreads   if ($taskThreads > $hostThreads);
    $taskMemory  = $hostMemory    if ($taskMemory  > $hostMemory);

    #  Reset the global values for later use.

    setGlobal("${tag}${alg}Memory",  $taskMemory);
    setGlobal("${tag}${alg}Threads", $taskThreads);

    #  Finally, reset the concurrency (if we're master) so we don't swamp our poor workstation.

    my $concurrent = undef;

    if (($cls eq "master") &&
        (globalExists("${tag}${alg}Concurrency"))) {
        my $nc = int($hostThreads / $taskThreads);

        if (($taskThreads * getGlobal("${tag}${alg}Concurrency") > $hostThreads)) {
            $err .= "-- Reset concurrency from ", getGlobal("${tag}${alg}Concurrency"), " to $nc.\n";
            setGlobal("${tag}${alg}Concurrency", $nc);
        }

        if (!defined(getGlobal("${tag}${alg}Concurrency"))) {
            setGlobal("${tag}${alg}Concurrency", $nc);
        }

        $concurrent = getGlobal("${tag}${alg}Concurrency");
    }

    #  And report.

    $all .= "-- Allowed to";
    $all .= " run " . substr("   $concurrent", -3) . " job" . (($concurrent == 1) ? " " : "s") . " concurrently,"  if (defined($concurrent));
    $all .= " run    under grid control,"                                                                          if (!defined($concurrent));
    $all .= " and use up to " . substr("   $taskThreads", -3) . " compute thread" . (($taskThreads == 1) ? " " : "s");
    $all .= " and " . substr("   $taskMemory", -4) . " GB memory for stage '$nam'.\n";

    return($err, $all);
}




sub dirname ($) {
    my $d = shift @_;

    return($d)  if (-d $d);

    my @d = split '/', $d;
    pop @d;

    $d = join('/', @d);

    return($d);
}


sub diskSpace ($) {
    my $wrk   = dirname($_[0]);
    my $df    = df($wrk, 1024);

    my $total = int(10 * $df->{blocks} / 1048576) / 10;
    my $used  = int(10 * $df->{used}   / 1048576) / 10;
    my $free  = int(10 * $df->{bfree}  / 1048576) / 10;
    my $avail = int(10 * $df->{bavail} / 1048576) / 10;

    #print STDERR "Disk space: total $total GB, used $used GB, free $free GB, available $avail GB\n";

    return (wantarray) ? ($total, $used, $free, $avail) : $avail;
}




sub printHelp ($) {
    my $bin = shift @_;  #  Can't include canu::Execution without a loop.

    if (getGlobal("version")) {
        system("$bin/gatekeeperCreate --version");
        system("$bin/overlapInCore    --version");
        system("$bin/bogart           --version");
        system("$bin/utgcns           --version");
        exit(0);
    }

    if (getGlobal("options")) {
        foreach my $k (sort values %synnam) {
            my $o = substr("$k                                    ", 0, 35);
            my $d = substr($global{$k}   . "                      ", 0, 20);
            my $u = $synops{$k};

            if (!defined($global{$k})) {
                $d = substr("<unset>                    ", 0, 20);
            }

            print "$o$d($u)\n";
        }
        exit(0);
    }

    if (getGlobal("help") ne "") {
        print "\n";
        print "usage: canu [-correct | -trim | -assemble] \\\n";
        print "            [-s <assembly-specifications-file>] \\\n";
        print "             -p <assembly-prefix> \\\n";
        print "             -d <assembly-directory> \\\n";
        print "             genomeSize=<number>[g|m|k] \\\n";
        print "             errorRate=0.X \\\n";
        print "            [other-options] \\\n";
        print "            [-pacbio-raw | -pacbio-corrected | -nanopore-raw | -nanopore-corrected] *fastq\n";
        print "\n";
        print "  By default, all three stages (correct, trim, assemble) are computed.\n";
        print "  To compute only a single stage, use:\n";
        print "    -correct  - generate corrected reads\n";
        print "    -trim     - generate trimmed reads\n";
        print "    -assemble - generate an assembly\n";
        print "\n";
        print "  The assembly is computed in the (created) -d <assembly-directory>, with most\n";
        print "  files named using the -p <assembly-prefix>.\n";
        print "\n";
        print "  The genome size is your best guess of the genome size of what is being assembled.\n";
        print "  It is used mostly to compute coverage in reads.  Fractional values are allowed: '4.7m'\n";
        print "  is the same as '4700k' and '4700000'\n";
        print "\n";
        print "  The errorRate is not used correctly (we're working on it).  Set it to 0.06 and\n";
        print "  use the various utg*ErrorRate options.\n";
        print "\n";
        print "  A full list of options can be printed with '-options'.  All options\n";
        print "  can be supplied in an optional sepc file.\n";
        print "\n";
        print "  Reads can be either FASTA or FASTQ format, uncompressed, or compressed\n";
        print "  with gz, bz2 or xz.  Reads are specified by the technology they were\n";
        print "  generated with:\n";
        print "    -pacbio-raw         <files>\n";
        print "    -pacbio-corrected   <files>\n";
        print "    -nanopore-raw       <files>\n";
        print "    -nanopore-corrected <files>\n";
        print "\n";
        print "Complete documentation at http://canu.github.io/\n";
        print "\n";
        print $global{"help"};
        exit(0);
    }

    undef $global{"version"};
    undef $global{"options"};
    undef $global{"help"};
}



sub makeAbsolute ($) {
    my $var = shift @_;
    my $val = getGlobal($var);
    if (defined($val) && ($val !~ m!^/!)) {
        $val = "$ENV{'PWD'}/$val";
        setGlobal($var, $val);
        $val =~ s/\\\"/\"/g;
        $val =~ s/\"/\\\"/g;
        $val =~ s/\\\$/\$/g;
        $val =~ s/\$/\\\$/g;

        addCommandLineOption("\"$var=$val\"");
    }
}



sub fixCase ($) {
    my $var = shift @_;
    my $val = getGlobal($var);

    if (defined($val)) {
        $val =~ tr/A-Z/a-z/;
        setGlobal($var, $val);
    }
}



sub setParametersFromFile ($@) {
    my $specFile  = shift @_;
    my @fragFiles = @_;

    #  Client should be ensuring that the file exists before calling this function.
    die "specFile '$specFile' not found.\n"  if (! -e "$specFile");

    $specLog .= "\n";
    $specLog .= "###\n";
    $specLog .= "###  Reading options from '$specFile'\n";
    $specLog .= "###\n";
    $specLog .= "\n";

    open(F, "< $specFile") or caExit("can't open '$specFile' for reading: $!", undef);

    while (<F>) {
        $specLog .= $_;

        s/^\s+//;
        s/\s+$//;

        next if (m/^\s*\#/);
        next if (m/^\s*$/);

        if (-e $_) {
            my $xx = $_;
            $xx = "$ENV{'PWD'}/$xx" if ($xx !~ m!^/!);
            if (-e $xx) {
                push @fragFiles, $xx;
            } else {
                setGlobal("help", getGlobal("help") . "File not found '$_' after appending absolute path.\n");
            }
        } elsif (m/\s*(\w*)\s*=([^#]*)#*.*$/) {
            my ($var, $val) = ($1, $2);
            $var =~ s/^\s+//; $var =~ s/\s+$//;
            $val =~ s/^\s+//; $val =~ s/\s+$//;
            undef $val if ($val eq "undef");
            setGlobal($var, $val);
        } else {
            setGlobal("help", getGlobal("help") . "File not found or unknown specFile option line '$_'.\n");
        }
    }
    close(F);

    return(@fragFiles);
}



sub setParametersFromCommandLine(@) {
    my @specOpts = @_;

    if (scalar(@specOpts) > 0) {
        $specLog .= "\n";
        $specLog .= "###\n";
        $specLog .= "###  Reading options from the command line.\n";
        $specLog .= "###\n";
        $specLog .= "\n";
    }

    foreach my $s (@specOpts) {
        $specLog .= "$s\n";

        if ($s =~ m/\s*(\w*)\s*=(.*)/) {
            my ($var, $val) = ($1, $2);
            $var =~ s/^\s+//; $var =~ s/\s+$//;
            $val =~ s/^\s+//; $val =~ s/\s+$//;
            setGlobal($var, $val);
        } else {
            setGlobal("help", getGlobal("help") . "Misformed command line option '$s'.\n");
        }
    }
}



sub adjustMemoryValue ($) {
    my $val = shift @_;

    $val = $1                if ($val =~ m/(\d+.*\d*)g/);
    $val = $1 / 1024         if ($val =~ m/(\d+.*\d*)m/);
    $val = $1 / 1024 / 1024  if ($val =~ m/(\d+.*\d*)k/);

    return($val);
}


sub checkParameters ($) {
    my $bin = shift @_;  #  Can't include canu::Execution without a loop.

    #
    #  Pick a nice looking set of binaries, and check them.
    #

    caExit("can't find 'gatekeeperCreate' program in $bin.  Possibly incomplete installation", undef) if (! -x "$bin/gatekeeperCreate");
    caExit("can't find 'meryl' program in $bin.  Possibly incomplete installation", undef)            if (! -x "$bin/meryl");
    caExit("can't find 'overlapInCore' program in $bin.  Possibly incomplete installation", undef)    if (! -x "$bin/overlapInCore");
    caExit("can't find 'bogart' program in $bin.  Possibly incomplete installation", undef)           if (! -x "$bin/bogart");
    caExit("can't find 'utgcns' program in $bin.  Possibly incomplete installation", undef)           if (! -x "$bin/utgcns");

    #
    #  Fiddle with filenames to make them absolute paths.
    #

    makeAbsolute("pathMap");

    makeAbsolute("corOvlFrequentMers");
    makeAbsolute("obtOvlFrequentMers");
    makeAbsolute("utgOvlFrequentMers");

    #
    #  Adjust case on some of them
    #

    fixCase("corOverlapper");
    fixCase("obtOverlapper");
    fixCase("utgOverlapper");

    fixCase("corConsensus");
    fixCase("cnsConsensus");

    fixCase("corFilter");

    fixCase("unitigger");
    fixCase("stopBefore");
    fixCase("stopAfter");

    #
    #  Check for inconsistent parameters
    #

    if (getGlobal("minReadLength") < getGlobal("minOverlapLength")) {
        my $mr = getGlobal("minReadLength");
        my $mo = getGlobal("minOverlapLength");

        caExit("minReadLength=$mr must be at least minOverlapLength=$mo", undef);

        print STDERR "-- WARNING: minReadLength reset from $mr to $mo (limited by minOverlapLength)\n";

        setGlobal("minOverlapLength", $mo);
    }

    #
    #  Adjust memory to be gigabytes.
    #

    foreach my $key (keys %global) {
        next  if ($key =~ m/gridEngineMemoryOption/i);
        next  if ($key !~ m/Memory/i);

        setGlobal($key, adjustMemoryValue(getGlobal($key)));
    }

    #
    #  Check for invalid usage
    #

    foreach my $tag ("cor", "obt", "utg") {
        if ((getGlobal("${tag}Overlapper") ne "mhap") &&
            (getGlobal("${tag}Overlapper") ne "ovl")) {
            caExit("invalid '${tag}Overlapper' specified (" . getGlobal("${tag}Overlapper") . "); must be 'mhap' or 'ovl'", undef);
        }
    }

    if ((getGlobal("unitigger") ne "unitigger") &&
        (getGlobal("unitigger") ne "bogart")) {
        caExit("invalid 'unitigger' specified (" . getGlobal("unitigger") . "); must be 'unitigger' or 'bogart'", undef);
    }

    foreach my $tag ("cor", "cns") {
        if ((getGlobal("${tag}consensus") ne "utgcns") &&
            (getGlobal("${tag}consensus") ne "falcon") &&
            (getGlobal("${tag}consensus") ne "falconpipe") &&
            (getGlobal("${tag}consensus") ne "pbdagcon") &&
            (getGlobal("${tag}consensus") ne "pbutgcns")) {
            caExit("invalid 'consensus' specified (" . getGlobal("${tag}consensus") . "); must be 'utgcns' or 'falcon' or 'falconpipe' or 'pbdagcon' or 'pbutgcns'", undef);
        }
    }


    if ((!defined("lowCoverageAllowed") &&  defined("lowCoverageDepth")) ||
        ( defined("lowCoverageAllowed") && !defined("lowCoverageDepth"))) {
        caExit("invalid 'lowCoverageAllowed' and 'lowCoverageDepth' specified; both must be set", undef);
    }


    #if ((getGlobal("cleanup") ne "none") &&
    #    (getGlobal("cleanup") ne "light") &&
    #    (getGlobal("cleanup") ne "heavy") &&
    #    (getGlobal("cleanup") ne "aggressive")) {
    #    caExit("invalid cleaup specified (" . getGlobal("cleanup") . "); must be 'none', 'light', 'heavy' or 'aggressive'", undef);
    #}

    if ((getGlobal("corFilter") ne "quick") &&
        (getGlobal("corFilter") ne "expensive")) {
        caExit("invalid 'corFilter' specified (" . getGlobal("corFilter") . "); must be 'quick' or 'expensive'", undef);
    }

    if (defined(getGlobal("stopBefore"))) {
        my $ok = 0;
        my $st = getGlobal("stopBefore");
        $st =~ tr/A-Z/a-z/;

        my $failureString = "Invalid stopBefore specified (" . getGlobal("stopBefore") . "); must be one of:\n";

        my @stopBefore = ("gatekeeper",
                          "meryl",
                          "trimReads",
                          "splitReads",
                          "unitig",
                          "consensusConfigure");

        foreach my $sb (@stopBefore) {
            $failureString .= "    '$sb'\n";
            $sb =~ tr/A-Z/a-z/;
            if ($st eq $sb) {
                $ok++;
                setGlobal('stopBefore', $st);
            }
        }

        caExit($failureString, undef) if ($ok == 0);
    }

    if (defined(getGlobal("stopAfter"))) {
        my $ok = 0;
        my $st = getGlobal("stopAfter");
        $st =~ tr/A-Z/a-z/;

        my $failureString = "Invalid stopAfter specified (" . getGlobal("stopAfter") . "); must be one of:\n";

        my @stopAfter = ("gatekeeper",
                         "meryl",
                         "mhapConfigure",
                         "overlapStoreConfigure",
                         "overlapStore",
                         "unitig",
                         "consensusConfigure",
                         "consensusCheck",
                         "consensusLoad",
                         "consensusFilter");

        foreach my $sa (@stopAfter) {
            $failureString .= "    '$sa'\n";
            $sa =~ tr/A-Z/a-z/;
            if ($st eq $sa) {
                $ok++;
                setGlobal('stopAfter', $st);
            }
        }

        caExit($failureString, undef) if ($ok == 0);
    }

    caExit("required parameter 'errorRate' is not set", undef)   if (! defined(getGlobal("errorRate")));
    caExit("required parameter 'genomeSize' is not set", undef)  if (! defined(getGlobal("genomeSize")));

    setGlobal("genomeSize", $1 * 1000)        if (getGlobal("genomeSize") =~ m/(\d+.*\d*)k/i);
    setGlobal("genomeSize", $1 * 1000000)     if (getGlobal("genomeSize") =~ m/(\d+.*\d*)m/i);
    setGlobal("genomeSize", $1 * 1000000000)  if (getGlobal("genomeSize") =~ m/(\d+.*\d*)g/i);

    #
    #  Java?  Need JRE 1.8.
    #

    if ((getGlobal("corOverlapper") eq "mhap") ||
        (getGlobal("obtOverlapper") eq "mhap") ||
        (getGlobal("utgOverlapper") eq "mhap")) {
        my $java       = getGlobal("java");
        my $versionStr = "unknown";
        my $version    = 0;

        open(F, "$java -showversion 2>&1 |");
        while (<F>) {
            #  First word is either "java" or "openjdk" or ...
            if (m/^.*\s+version\s+\"(\d+.\d+)(.*)\"$/) {
                $versionStr = "$1$2";
                $version    =  $1;
            }
        }
        close(F);

        print STDERR "-- Detected Java(TM) Runtime Environment '$versionStr' (from '$java').\n";

        caExit("mhap overlapper requires java version at least 1.8.0; you have $versionStr", undef)  if ($version < 1.8);
    }

    #
    #  Falcon?  Need to find it.
    #

    if ((getGlobal("corConsensus") eq "falcon") ||
        (getGlobal("corConsensus") eq "falconpipe")) {
        my $falcon = getGlobal("falconSense");

        caExit("didn't find falcon program with option falconSense='$falcon'", undef)  if ((defined($falcon)) && (! -e $falcon));
    }


    #
    #  Finish grid configuration.  If any of these are set, they were set by the user.
    #

    #  Handle special cases.

    if (uc(getGlobal("gridEngine")) eq "SGE") {
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
        setGlobalIfUndef("gridEngineMemoryOption",               undef);  #"-l mem=MEMORY");
        setGlobalIfUndef("gridEngineNameToJobIDCommand",         undef);
        setGlobalIfUndef("gridEngineNameToJobIDCommandNoArray",  undef);
        setGlobalIfUndef("gridEngineTaskID",                     "\$SGE_TASK_ID");
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
                print STDERR "--          Set 'gridEngineMemoryOption' manually; example: '-l mem=MEMORY'.\n";

            } else {
                print STDERR "-- WARNING:  Couldn't determine the SGE resource to request memory.\n";
                print STDERR "--          Set 'gridEngineMemoryOption' manually; example: '-l mem=MEMORY'.\n";
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

        open(F, "qhost -ncb |");
        $_ = <F>;  #  Header
        $_ = <F>;  #  Table bar
        while (<F>) {
            my @v = split '\s+', $_;

            next if ($v[3] eq "-");  #  Node disabled or otherwise not available

            my $cpus = $v[2];
            my $mem  = $v[4];

            $mem  = $1 * 1024  if ($mem =~ m/(\d+.*\d+)[tT]/);
            $mem  = $1 * 1     if ($mem =~ m/(\d+.*\d+)[gG]/);
            $mem  = $1 / 1024  if ($mem =~ m/(\d+.*\d+)[mM]/);
            $mem  = int($mem);

            $hosts{"$cpus-$mem"}++;
        }
        close(F);

        print STDERR "--\n";

        foreach my $c (keys %hosts) {
            my ($cpus, $mem) = split '-', $c;
            my  $nodes       = $hosts{$c};

            printf(STDERR "-- Found %3d host%s with %3d core%s and %4d GB memory under Sun Grid Engine control.\n",
                   $nodes, ($nodes == 1) ? " " : "s",
                   $cpus,  ($cpus  == 1) ? " " : "s",
                   $mem);

            $hosts .= "\0"                  if (defined($hosts));
            $hosts .= "$cpus-$mem-$nodes";
        }

        setGlobal("availableHosts", $hosts);
    }

    if (uc(getGlobal("gridEngine")) eq "PBS") {
        setGlobalIfUndef("gridEngineSubmitCommand",              "qsub");
        setGlobalIfUndef("gridEngineHoldOption",                 "-W depend=afterany:\"WAIT_TAG\"");
        setGlobalIfUndef("gridEngineHoldOptionNoArray",          undef);
        setGlobalIfUndef("gridEngineSyncOption",                 "");
        setGlobalIfUndef("gridEngineNameOption",                 "-d `pwd` -N");
        setGlobalIfUndef("gridEngineArrayOption",                "-t ARRAY_JOBS");
        setGlobalIfUndef("gridEngineArrayName",                  "ARRAY_NAME\[ARRAY_JOBS\]");
        setGlobalIfUndef("gridEngineOutputOption",               "-j oe -o");
        setGlobalIfUndef("gridEnginePropagateCommand",           "qalter -W depend=afterany:\"WAIT_TAG\"");
        setGlobalIfUndef("gridEngineThreadsOption",              undef);
        setGlobalIfUndef("gridEngineMemoryOption",               undef);
        setGlobalIfUndef("gridEngineNameToJobIDCommand",         undef);
        setGlobalIfUndef("gridEngineNameToJobIDCommandNoArray",  undef);
        setGlobalIfUndef("gridEngineTaskID",                     "\$PBS_TASKNUM");
        setGlobalIfUndef("gridEngineArraySubmitID",              "\\\$PBS_TASKNUM");
        setGlobalIfUndef("gridEngineJobID",                      "PBS_JOBID");
    }

    if (uc(getGlobal("gridEngine")) eq "LSF") {
        setGlobalIfUndef("gridEngineSubmitCommand",              "bsub");
        setGlobalIfUndef("gridEngineHoldOption",                 "-w \"numended\(\"WAIT_TAG\", \*\)\"");
        setGlobalIfUndef("gridEngineHoldOptionNoArray",          "-w \"done\(\"WAIT_TAG\"\)\"");
        setGlobalIfUndef("gridEngineSyncOption",                 "-K");
        setGlobalIfUndef("gridEngineNameOption",                 "-J");
        setGlobalIfUndef("gridEngineArrayOption",                "");
        setGlobalIfUndef("gridEngineArrayName",                  "ARRAY_NAME\[ARRAY_JOBS\]");
        setGlobalIfUndef("gridEngineOutputOption",               "-o");
        setGlobalIfUndef("gridEnginePropagateCommand",           "bmodify -w \"done\(\"WAIT_TAG\"\)\"");
        setGlobalIfUndef("gridEngineThreadsOption",              undef);
        setGlobalIfUndef("gridEngineMemoryOption",               undef);
        setGlobalIfUndef("gridEngineNameToJobIDCommand",         "bjobs -A -J \"WAIT_TAG\" | grep -v JOBID");
        setGlobalIfUndef("gridEngineNameToJobIDCommandNoArray",  "bjobs -J \"WAIT_TAG\" | grep -v JOBID");
        setGlobalIfUndef("gridEngineTaskID",                     "\$LSB_JOBINDEX");
        setGlobalIfUndef("gridEngineArraySubmitID",              "%I");
        setGlobalIfUndef("gridEngineJobID",                      "LSB_JOBID");
    }

    #
    #  Set default error rates based on the per-read error rate.
    #

    setGlobalIfUndef("corOvlErrorRate",      3.0 * getGlobal("errorRate"));
    setGlobalIfUndef("obtOvlErrorRate",      3.0 * getGlobal("errorRate"));
    setGlobalIfUndef("utgOvlErrorRate",      3.0 * getGlobal("errorRate"));

    setGlobalIfUndef("ovlErrorRate",	     2.5 * getGlobal("errorRate"));
    setGlobalIfUndef("utgGraphErrorRate",    2.0 * getGlobal("errorRate"));
    setGlobalIfUndef("utgBubbleErrorRate",   2.0 * getGlobal("errorRate") + 0.5 * getGlobal("errorRate"));
    setGlobalIfUndef("utgMergeErrorRate",    2.0 * getGlobal("errorRate") - 0.5 * getGlobal("errorRate"));
    setGlobalIfUndef("utgRepeatErrorRate",   1.0 * getGlobal("errorRate"));

    setGlobalIfUndef("corsErrorRate",        8.0 * getGlobal("errorRate"));
    setGlobalIfUndef("cnsErrorRate",         3.0 * getGlobal("errorRate"));

    #
    #  Finally, process all the memory/thread settings to make them compatible with our
    #  hardware.  These wanted to be closer to their use (like right before the various
    #  scripts are written), but it needs to be done everytime canu starts, otherwise,
    #  the canu invocation that writes scripts get the fixes, and the canu invocation
    #  that runs the scripts does not.
    #

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


sub setExecDefaults ($$$$) {
    my $tag         = shift @_;
    my $name        = shift @_;
    my $memory      = shift @_;
    my $threads     = shift @_;

    $global{"useGrid${tag}"}       = 1;
    $synops{"useGrid${tag}"}       = "Use grid engine for $name computes";

    $global{"gridOptions${tag}"}   = undef;
    $synops{"gridOptions${tag}"}   = "Grid engine options applied to $name jobs";

    $global{"${tag}Memory"}        = $memory;
    $synops{"${tag}Memory"}        = "Amount of memory, in gigabytes, to use for $name jobs";

    $global{"${tag}Threads"}       = $threads;
    $synops{"${tag}Threads"}       = "Number of threads to use for $name jobs";

    $global{"${tag}Concurrency"}   = undef;
    $synops{"${tag}Concurrency"}   = "If grid not enabled, number of $name jobs to run at the same time; default is n_proc / n_threads";
}



sub showErrorRates ($) {
    my $prefix = shift @_;

    print STDERR "${prefix}\n";
    print STDERR "${prefix}genomeSize          -- ", getGlobal("genomeSize"), "\n";
    print STDERR "${prefix}errorRate           -- ", getGlobal("errorRate"), "\n";
    print STDERR "${prefix}\n";
    print STDERR "${prefix}corOvlErrorRate     -- ", getGlobal("corOvlErrorRate"), "\n";
    print STDERR "${prefix}obtOvlErrorRate     -- ", getGlobal("obtOvlErrorRate"), "\n";
    print STDERR "${prefix}utgOvlErrorRate     -- ", getGlobal("utgOvlErrorRate"), "\n";
    print STDERR "${prefix}\n";
    print STDERR "${prefix}obtErrorRate        -- ", getGlobal("obtErrorRate"), "\n";
    print STDERR "${prefix}\n";
    print STDERR "${prefix}utgGraphErrorRate   -- ", getGlobal("utgGraphErrorRate"), "\n";
    print STDERR "${prefix}utgBubbleErrorRate  -- ", getGlobal("utgBubbleErrorRate"), "\n";
    print STDERR "${prefix}utgMergeErrorRate   -- ", getGlobal("utgMergeErrorRate"), "\n";
    print STDERR "${prefix}utgRepeatErrorRate  -- ", getGlobal("utgRepeatErrorRate"), "\n";
    print STDERR "${prefix}\n";
    print STDERR "${prefix}corErrorRate        -- ", getGlobal("corErrorRate"), "\n";
    print STDERR "${prefix}cnsErrorRate        -- ", getGlobal("cnsErrorRate"), "\n";
}



#  Defaults are set for yeast:
#    trimming   errorRate = 0.009  obtOvlErrorRate = 0.06  obtErrorRate = 0.035
#    assembly   errorRate = 0.009  utgOvlErrorRate = 0.06  bogart 0.035
#
sub setErrorRate ($@) {
    my $er      = shift @_;
    my $verbose = shift @_;

    print STDERR "-- Set errorRate to $er (verbose='$verbose')\n"  if (defined($verbose));

    #  Can NOT call setGlobal() for this, because it calls setErrorRate()!.
    $global{"errorrate"} = $er;
    setGlobal("corOvlErrorRate",    $er * 3);  #  Not used, except for realigning
    setGlobal("obtOvlErrorRate",    $er * 3);  #  Generally must be smaller than utgGraphErrorRate
    setGlobal("utgOvlErrorRate",    $er * 3);

    setGlobal("obtErrorRate",       $er * 2.5);

    setGlobal("utgGraphErrorRate",  $er * 2);
    setGlobal("utgBubbleErrorRate", $er * 2 + 0.5 * $er);  #  Not tested!
    setGlobal("utgMergeErrorRate",  $er * 2 - 0.5 * $er);
    setGlobal("utgRepeatErrorRate", $er * 2);

    setGlobal("corErrorRate",       $er * 10);  #  Erorr rate used for raw sequence alignment/consensus
    setGlobal("cnsErrorRate",       $er);

    showErrorRates("--  ")  if (defined($verbose));
}



sub setOverlapDefaults ($$$) {
    my $tag     = shift @_;  #  If 'cor', some parameters are loosened for raw pacbio reads
    my $name    = shift @_;
    my $default = shift @_;  #  Sets ${tag}Overlapper

    #  Which overlapper to use.

    $global{"${tag}Overlapper"}                  = $default;
    $synops{"${tag}Overlapper"}                  = "Which overlap algorithm to use for $name";

    #  OverlapInCore parameters.

    $global{"${tag}OvlHashBlockLength"}          = ($tag eq "cor") ? 2500000 : 100000000;
    $synops{"${tag}OvlHashBlockLength"}          = "Amount of sequence (bp) to load into the overlap hash table";

    $global{"${tag}OvlRefBlockSize"}             = ($tag eq "cor") ? 20000 : 2000000;
    $synops{"${tag}OvlRefBlockSize"}             = "Number of reads to search against the hash table per batch";

    $global{"${tag}OvlRefBlockLength"}           = 0;
    $synops{"${tag}OvlRefBlockLength"}           = "Amount of sequence (bp) to search against the hash table per batch";

    $global{"${tag}OvlHashBits"}                 = ($tag eq "cor") ? 18 : 22;
    $synops{"${tag}OvlHashBits"}                 = "Width of the kmer hash.  Width 22=1gb, 23=2gb, 24=4gb, 25=8gb.  Plus 10b per ${tag}OvlHashBlockLength";

    $global{"${tag}OvlHashLoad"}                 = 0.75;
    $synops{"${tag}OvlHashLoad"}                 = "Maximum hash table load.  If set too high, table lookups are inefficent; if too low, search overhead dominates run time";

    $global{"${tag}OvlMerSize"}                  = ($tag eq "cor") ? 19 : 22;
    $synops{"${tag}OvlMerSize"}                  = "K-mer size for seeds in overlaps";

    $global{"${tag}OvlMerThreshold"}             = "auto";
    $synops{"${tag}OvlMerThreshold"}             = "K-mer frequency threshold; mers more frequent than this count are ignored";

    $global{"${tag}OvlMerDistinct"}              = undef;
    $synops{"${tag}OvlMerDistinct"}              = "K-mer frequency threshold; the least frequent fraction of distinct mers can seed overlaps";

    $global{"${tag}OvlMerTotal"}                 = undef;
    $synops{"${tag}OvlMerTotal"}                 = "K-mer frequency threshold; the least frequent fraction of all mers can seed overlaps";

    $global{"${tag}OvlFrequentMers"}             = undef;
    $synops{"${tag}OvlFrequentMers"}             = "Do not seed overlaps with these kmers (fasta format)";

    #  Mhap parameters.

    $global{"${tag}MhapBlockSize"}            = 20000;
    $synops{"${tag}MhapBlockSize"}            = "Number of reads per block; one block is loaded into memory per job";

    $global{"${tag}MhapMerSize"}              = ($tag eq "cor") ? 16 : 22;
    $synops{"${tag}MhapMerSize"}              = "K-mer size for seeds in mhap";

    $global{"${tag}MhapReAlign"}              = undef;
    $synops{"${tag}MhapReAlign"}              = "Compute actual alignments from mhap overlaps; 'raw' from mhap output, 'final' from overlap store; uses either obtErrorRate or ovlErrorRate, depending on which overlaps are computed";

    $global{"${tag}MhapSensitivity"}          = "normal";
    $synops{"${tag}MhapSensitivity"}          = "Coarse sensitivity level: 'normal' or 'high'";
}



sub setDefaults () {

    #####  General Configuration Options (aka miscellany)

    $global{"canuIteration"}               = 0;  #  See documentation in Execution.pm
    $global{"canuIterationMax"}            = 2;

    $global{"showNext"}                    = undef;
    $synops{"showNext"}                    = "Don't run any commands, just report what would run";

    $global{"pathMap"}                     = undef;
    $synops{"pathMap"}                     = "File with a hostname to binary directory map";

    $global{"shell"}                       = "/bin/sh";
    $synops{"shell"}                       = "Command interpreter to use; sh-compatible (e.g., bash), NOT C-shell (csh or tcsh)";

    $global{"java"}                        = "java";
    #$global{"java"}                        = "/usr/bin/java"                               if (-e "/usr/bin/java");
    #$global{"java"}                        = "/usr/local/bin/java"                         if (-e "/usr/local/bin/java");
    $global{"java"}                        = "/nbacc/local/packages/jdk1.8.0_25/bin/java"  if (-e "/nbacc/local/packages/jdk1.8.0_25/bin/java");
    $synops{"java"}                        = "Java interpreter to use; at least version 1.8";

    #####  Cleanup options

    $global{"saveOverlaps"}                = 0;
    $synops{"saveOverlaps"}                = "Save intermediate overlap files, almost never a good idea";

    $global{"saveMerCounts"}               = 0;
    $synops{"saveMerCounts"}               = "Save full mer counting results, sometimes useful";

    #####  Error Rates

    $global{"errorRate"}                   = undef;
    $synops{"errorRate"}                   = "The expected error rate in the input reads";

    $global{"corOvlErrorRate"}             = undef;
    $synops{"corOvlErrorRate"}             = "Overlaps above this error rate are not computed";

    $global{"obtOvlErrorRate"}             = undef;
    $synops{"obtOvlErrorRate"}             = "Overlaps at or below this error rate are used to trim reads";

    $global{"utgOvlErrorRate"}             = undef;
    $synops{"utgOvlErrorRate"}             = "Overlaps at or below this error rate are used to trim reads";

    #$global{"utgErrorRate"}                = undef;
    #$synops{"utgErrorRate"}                = "Overlaps at or below this error rate are used to construct unitigs (BOG and UTG)";

    $global{"utgGraphErrorRate"}           = undef;
    $synops{"utgGraphErrorRate"}           = "Overlaps at or below this error rate are used to construct unitigs (BOGART)";

    $global{"utgBubbleErrorRate"}          = undef;
    $synops{"utgBubbleErrorRate"}          = "Overlaps at or below this error rate are used to construct unitigs (BOGART)";

    $global{"utgMergeErrorRate"}           = undef;
    $synops{"utgMergeErrorRate"}           = "Overlaps at or below this error rate are used to construct unitigs (BOGART)";

    $global{"utgRepeatErrorRate"}          = undef;
    $synops{"utgRepeatErrorRate"}          = "Overlaps at or below this error rate are used to construct unitigs (BOGART)";

    $global{"corErrorRate"}                = undef;
    $synops{"corErrorRate"}                = "Only use raw alignments below this error rate to construct corrected reads";

    $global{"cnsErrorRate"}                = undef;
    $synops{"cnsErrorRate"}                = "Consensus expects alignments at about this error rate";

    #####  Minimums

    $global{"minReadLength"}               = 1000;
    $synops{"minReadLength"}               = "Reads shorter than this length are not loaded into the assembler";

    $global{"minOverlapLength"}            = 500;
    $synops{"minOverlapLength"}            = "Overlaps shorter than this length are not computed";

    #####  Stopping conditions

    $global{"stopBefore"}                  = undef;
    $synops{"stopBefore"}                  = "Tell canu when to halt execution";

    $global{"stopAfter"}                   = undef;
    $synops{"stopAfter"}                   = "Tell canu when to halt execution";

    #####  Grid Engine configuration, internal parameters

    $global{"availableHosts"}                       = undef;  #  Internal list of cpus-memory-nodes describing the grid

    $global{"gridEngine"}                           = undef;
    $global{"gridEngineSubmitCommand"}              = undef;
    $global{"gridEngineHoldOption"}                 = undef;
    $global{"gridEngineHoldOptionNoArray"}          = undef;
    $global{"gridEngineSyncOption"}                 = undef;
    $global{"gridEngineNameOption"}                 = undef;
    $global{"gridEngineArrayOption"}                = undef;
    $global{"gridEngineArrayName"}                  = undef;
    $global{"gridEngineOutputOption"}               = undef;
    $global{"gridEnginePropagateCommand"}           = undef;
    $global{"gridEngineThreadsOption"}              = undef;
    $global{"gridEngineMemoryOption"}               = undef;
    $global{"gridEngineNameToJobIDCommand"}         = undef;
    $global{"gridEngineNameToJobIDCommandNoArray"}  = undef;
    $global{"gridEngineTaskID"}                     = undef;
    $global{"gridEngineArraySubmitID"}              = undef;
    $global{"gridEngineJobID"}                      = undef;

    #  Try to decide which grid engine we have.  If this isn't set, the execution methods
    #  submitScript() and submitOrRunParallelJob() will return without submitting, or run locally
    #  (respectively).  This means that we can just trivially change the defaults for useGrid and
    #  useGridMaster to 'enabled' and it'll do the right thing when SGE isn't present.

    if (defined($ENV{'SGE_ROOT'})) {
        print STDERR "-- Detected Sun Grid Engine in '$ENV{'SGE_ROOT'}/$ENV{'SGE_CELL'}'.\n";
        $global{"gridEngine"} = "SGE";
    } else {
        print STDERR "-- No grid engine detected, grid disabled.\n";
    }

    #####  Grid Engine Pipeline

    $global{"useGrid"}                     = 1;
    $synops{"useGrid"}                     = "Enable SGE; if unset, no grid will be used";

    #####  Grid Engine configuration, for each step of the pipeline

    $global{"gridOptions"}                 = undef;
    $synops{"gridOptions"}                 = "Grid engine options applied to all jobs";

    $global{"gridOptionsMaster"}           = undef;
    $synops{"gridOptionsMaster"}           = "Grid engine options applied to sequential (usually high memory) jobs";

    $global{"gridOptionsJobName"}          = undef;
    $synops{"gridOptionsJobName"}          = "Grid jobs job-name suffix";

    #####  Grid Engine configuration and parameters, for each step of the pipeline (memory, threads)

    setExecDefaults("cns",    "unitig consensus",                       "4,6,8,10,12,14,16",  "1");
    setExecDefaults("cor",    "read correction",                        "4,6,8",              "2,4,6,8");
    setExecDefaults("red",    "read error detection",                   "4,6,8,10,12,14,16",  "2,4,6,8");
    setExecDefaults("oea",    "overlap error adjustment",               "4,6,8,10,12,14,16",  "1");

    setExecDefaults("corovl",  "overlaps for correction",               "4,6,8,12", "2,4,6");
    setExecDefaults("obtovl",  "overlaps for trimming",                 "4,6,8,12", "2,4,6");
    setExecDefaults("utgovl",  "overlaps for unitig construction",      "4,6,8,12", "2,4,6");

    setExecDefaults("cormhap", "mhap overlaps for correction",          "8,12,14,16,20,24,28,32,36,40,44,48,52,56,60,64", "4,8,12,16");
    setExecDefaults("obtmhap", "mhap overlaps for trimming",            "8,12,14,16,20,24,28,32,36,40,44,48,52,56,60,64", "4,8,12,16");
    setExecDefaults("utgmhap", "mhap overlaps for unitig construction", "8,12,14,16,20,24,28,32,36,40,44,48,52,56,60,65", "4,8,12,16");

    setExecDefaults("ovb",    "overlap store bucketizing",              "2,4",   "1");
    setExecDefaults("ovs",    "overlap store sorting",                  "4,6,8", "1");

    setExecDefaults("grid",   "parallel jobs",                          undef, undef);  #  gridThreads and gridMemory, other settingss are unused
    setExecDefaults("master", "master script",                          undef, undef);  #  masterThreads and masterMemory and masterConcurrency

    $global{"gridMemory"}         = undef;                     #  Unset.  User can limit if wanted.
    $global{"gridThreads"}        = undef;

    $global{"masterMemory"}       = getPhysicalMemorySize();   #  Set to this machine.
    $global{"masterThreads"}      = getNumberOfCPUs();
    #$global{"masterConcurrency"}  = undef;

    $global{"useGrid"}            = 1;
    $synops{"useGrid"}            = "Enable use of a grid, if one exists";

    $global{"useGridMaster"}      = 1;
    $synops{"useGridMaster"}      = "Run the canu pipeline entirely under grid control";

    #####  Overlapper

    setOverlapDefaults("cor", "correction",             "mhap");  #  Overlaps computed for correction
    setOverlapDefaults("obt", "overlap based trimming", "ovl");   #  Overlaps computed for trimming
    setOverlapDefaults("utg", "unitig construction",    "ovl");   #  Overlaps computed for unitigging

    ##### Overlap Store

    $global{"ovlStoreMemory"}              = 4;
    $synops{"ovlStoreMemory"}              = "How much memory, in gigabytes, to use when constructing overlap stores";

    $global{"ovlStoreThreads"}             = 1;
    $synops{"ovlStoreThreads"}             = "Unused, only one thread supported";

    $global{"ovlStoreConcurrency"}         = 1;
    $synops{"ovlStoreConcurrency"}         = "Unused, only one process supported";

    $global{"ovlStoreMethod"}              = "sequential";
    $synops{"ovlStoreMethod"}              = "Use the 'sequential' or 'parallel' algorithm for constructing an overlap store";

    $global{"ovlStoreSlices"}              = 128;
    $synops{"ovlStoreSlices"}              = "How many pieces to split the sorting into, for the parallel store build";

    #####  Mers

    $global{"merylMemory"}                 = undef;
    $synops{"merylMemory"}                 = "Amount of memory, in gigabytes, to use for mer counting";

    $global{"merylThreads"}                = undef;
    $synops{"merylThreads"}                = "Number of threads to use for mer counting";

    $global{"merylConcurrency"}            = "1";
    $synops{"merylConcurrency"}            = "Unused, there is only one process";

    #####  Overlap Based Trimming

    $global{"obtErrorRate"}                = undef;
    $synops{"obtErrorRate"}                = "Stringency of overlaps to use for trimming";

    $global{"trimReadsOverlap"}            = 1;
    $synops{"trimReadsOverlap"}            = "Minimum overlap between evidence to make contiguous trim";

    $global{"trimReadsCoverage"}           = 1;
    $synops{"trimReadsCoverage"}           = "Minimum depth of evidence to retain bases";

    #$global{"splitReads..."}               = 1;
    #$synops{"splitReads..."}               = "";

    #####  Fragment/Overlap Error Correction

    $global{"enableOEA"}                   = 1;
    $synops{"enableOEA"}                   = "Do overlap error adjustment - comprises two steps: read error detection (RED) and overlap error adjustment (OEA)";

    $global{"redBatchSize"}                = 0;
    $synops{"redBatchSize"}                = "Number of reads per fragment error detection batch";

    $global{"redBatchLength"}              = 0;
    $synops{"redBatchLength"}              = "Number of bases per fragment error detection batch";

    $global{"oeaBatchSize"}                = 0;
    $synops{"oeaBatchSize"}                = "Number of reads per overlap error correction batch";

    $global{"oeaBatchLength"}              = 0;
    $synops{"oeaBatchLength"}              = "Number of bases per overlap error correction batch";

    #####  Unitigger & BOG & bogart Options

    $global{"unitigger"}                   = "bogart";
    $synops{"unitigger"}                   = "Which unitig algorithm to use; utg or bogart (defalut)";

    $global{"genomeSize"}                  = undef;
    $synops{"genomeSize"}                  = "An estimate of the size of the genome";

    $global{"batOptions"}                  = undef;
    $synops{"batOptions"}                  = "Advanced options to bogart";

    $global{"batMemory"}                   = undef;
    $synops{"batMemory"}                   = "Approximate maximum memory usage for loading overlaps, in gigabytes, default is unlimited";

    $global{"batThreads"}                  = undef;
    $synops{"batThreads"}                  = "Number of threads to use in the Merge/Split/Join phase; default is whatever OpenMP wants";

    $global{"batConcurrency"}              = 1;
    $synops{"batConcurrency"}              = "Unused, only one process supported";

    #####  Unitig Filtering Options



    #####  Unitig Repeat/Unique Options (formerly in scaffolder)

    $global{"maxSingleReadSpan"}           = undef;  #  1.0 (default as in the binary)
    $synops{"maxSingleReadSpan"}           = "Unitigs with a single read spanning more than this fraction of the unitig are never labeled unique";

    $global{"lowCoverageAllowed"}          = undef;  #  1.0
    $synops{"lowCoverageAllowed"}          = "Unitigs with more than fraction lowCoverageAllowed bases at depth at most lowCoverageDepth bases are never labeled unique";

    $global{"lowCoverageDepth"}            = undef;  #  2
    $synops{"lowCoverageDepth"}            = "Unitigs with more than fraction lowCoverageAllowed bases at depth at most lowCoverageDepth bases are never labeled unique";

    $global{"minReadsUnique"}              = undef;  #  2
    $synops{"minReadsUnique"}              = "Unitigs with fewer reads that this are never labeled unique";

    $global{"minUniqueLength"}             = undef;  #  1000
    $synops{"minUniqueLength"}             = "Unitigs shorter than this are always labeled non-unique";

    $global{"maxRepeatLength"}             = undef;  #  max_int
    $synops{"maxRepeatLength"}             = "Unitigs longer than this are always labeled unique";

    #####  Consensus Options

    $global{"cnsPartitions"}               = 128;
    $synops{"cnsPartitions"}               = "Partition consensus into N jobs";

    $global{"cnsPartitionMin"}             = 75000;
    $synops{"cnsPartitionMin"}             = "Don't make a consensus partition with fewer than N reads";

    $global{"cnsMaxCoverage"}              = 0;
    $synops{"cnsMaxCoverage"}              = "Limit unitig consensus to at most this coverage";

    $global{"cnsConsensus"}                = "utgcns";
    $synops{"cnsConsensus"}                = "Which consensus algorithm to use; only 'utgcns' is supported";

    #####  Correction Options

    $global{"corPartitions"}               = 128;
    $synops{"corPartitions"}               = "Partition read correction into N jobs";

    $global{"corPartitionMin"}             = 25000;
    $synops{"corPartitionMin"}             = "Don't make a read correction partition with fewer than N reads";

    $global{"corMinEvidenceLength"}        = undef;
    $synops{"corMinEvidenceLength"}        = "Limit read correction to only overlaps longer than this; default: unlimited";

    $global{"corMaxEvidenceErate"}         = undef;
    $synops{"corMaxEvidenceErate"}         = "Limit read correction to only overlaps at or below this fraction error; default: unlimited";

    $global{"corMaxEvidenceCoverageGlobal"}= "1.0x";
    $synops{"corMaxEvidenceCoverageGlobal"}= "Limit reads used for correction to supporting at most this coverage; default: 1.0 * estimated coverage";

    $global{"corMaxEvidenceCoverageLocal"} = "2.0x";
    $synops{"corMaxEvidenceCoverageLocal"} = "Limit reads being corrected to at most this much evidence coverage; default: 10 * estimated coverage";

    $global{"corOutCoverage"}              = 40;
    $synops{"corOutCoverage"}              = "Only correct the longest reads up to this coverage; default 40";

    $global{"corMinCoverage"}              = 4;
    $synops{"corMinCoverage"}              = "Minimum number of bases supporting each corrected base, if less than this sequences are split";

    $global{"corFilter"}                   = "expensive";
    $synops{"corFilter"}                   = "Method to filter short reads from correction; 'quick' or 'expensive'";

    $global{"corConsensus"}                = "falconpipe";
    $synops{"corConsensus"}                = "Which consensus algorithm to use; only 'falcon' and 'falconpipe' are supported";

    $global{"falconSense"}                 = undef;
    $synops{"falconSense"}                 = "Path to fc_consensus.py or falcon_sense.bin";



    #####  Ugly, command line options passed to printHelp()

    $global{"help"}                        = "";
    $synops{"help"}                        = undef;

    $global{"version"}                     = 0;
    $synops{"version"}                     = undef;

    $global{"options"}                     = 0;
    $synops{"options"}                     = undef;

    #  Convert all the keys to lowercase, and remember the case-sensitive version

    foreach my $k (keys %global) {
        (my $l = $k) =~ tr/A-Z/a-z/;

        next  if ($k eq "version");
        next  if ($k eq "options");
        next  if ($k eq "help");

        if (! exists($synnam{$l})) {
            $synnam{$l} = $k;

            if (!exists($global{$l})) {
                $global{$l} = $global{$k};
                delete $global{$k};
            }

            #print "$k -> $l\n";
        }
    }

    #  If this is set, it breaks the consensus.sh and overlap.sh scripts.  Good grief!  Why
    #  are you running this in a task array!?

    if (exists($ENV{getGlobal("gridEngineTaskID")})) {
        undef $ENV{getGlobal("gridEngineTaskID")};
        print STDERR "ENV: ", getGlobal("gridEngineTaskID"), " needs to be unset, done.\n";
    }

    #  Finally, set the global default error rate

    setErrorRate(0.01);
}

1;
