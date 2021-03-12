
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

package canu::Grid_SGE;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(detectSGE configureSGE);

use strict;
use warnings "all";
no  warnings "uninitialized";

use canu::Defaults;
use canu::Execution;

use canu::Grid "formatAllowedResources";



sub detectSGE () {

    return   if ( defined(getGlobal("gridEngine")));
    return   if (!defined($ENV{'SGE_ROOT'}));

    if (getGlobal("useGrid") eq "0") {
        print STDERR "--\n";
        print STDERR "-- Detected Sun Grid Engine in '$ENV{'SGE_ROOT'}/$ENV{'SGE_CELL'}'.\n";
        print STDERR "--          Sun Grid Engine disabled by useGrid=false\n";
    }
    else {
        print STDERR "--\n";
        print STDERR "-- Detected Sun Grid Engine in '$ENV{'SGE_ROOT'}/$ENV{'SGE_CELL'}'.\n";

        setGlobal("gridEngine", "SGE");
    }
}


sub discoverThreadsOption () {
    my @env = `qconf -spl 2> /dev/null`;  chomp @env;
    my @thr;
    my $bestThr   = undef;
    my $bestSlots = 0;

    my @peErrorList;

    foreach my $env (@env) {
        my $ns = 0;
        my $ar = "unknown";

        open(F, "qconf -sp $env 2> /dev/null |");
        while (<F>) {
            $ns = $1   if (m/slots\s+(\d+)/);               #  How many slots can we use?
            $ar = $1   if (m/allocation_rule\s+(.*)/);
        }
        close(F);

        if ($ar ne "\$pe_slots") {
            push @peErrorList, "-- WARNING: Parallel Environment '$env' ignored; uses allocation_rule '$ar'\n";
            next;
        }

        push @thr, $env;

        if ($ns > $bestSlots) {
            $bestThr   = $env;
            $bestSlots = $ns;
        }
    }

    if (scalar(@thr) == 0) {
        print STDERR "--\n";

        foreach my $err (@peErrorList) {
            print STDERR $err;
        }

        print STDERR "--\n";
        print STDERR "-- ERROR:  No Sun Grid Engine Parallel Environment is using allocation rule '\$pe_threads'.\n";
        print STDERR "-- ERROR:  Multi-thread execution not supported on this grid.\n";
        print STDERR "--\n";

        caExit("can't configure for SGE", undef);
    }

    elsif (scalar(@thr) == 1) {
        print STDERR "-- Detected Sun Grid Engine parallel environment '$thr[0]'.\n";
        return($thr[0]);
    }

    else {
        print STDERR "--\n";
        print STDERR "-- WARNING:  Multiple Sun Grid Engine Parallel Environments can run multi-threaded codes.\n";
        print STDERR "-- WARNING:  Valid choices are:\n";
        print STDERR "-- WARNING:    $_\n"   foreach (@thr);
        print STDERR "-- WARNING:\n";
        print STDERR "-- WARNING:  Using Parallel Environment '$bestThr'.\n";
        print STDERR "--\n";

        return($bestThr);
    }
}


sub discoverMemoryOption () {
    my @mem;

    open(F, "qconf -sc 2> /dev/null |");
    while (<F>) {
        my @vals = split '\s+', $_;

        next  if (($vals[5] ne "YES") &&      #  Not a consumable resource.
                  ($vals[5] ne "JOB"));       #    ('JOB' is processed later)
        next  if  ($vals[2] ne "MEMORY");     #  Not a memory resource.
        next  if  ($vals[0] =~ m/swap/);      #  Don't care about swap.
        next  if  ($vals[0] =~ m/virtual/);   #  Don't care about vm space.

        push @mem, $vals[0];
    }
    close(F);

    if (scalar(@mem) == 1) {
        print STDERR "-- Detected Sun Grid Engine memory resource '$mem[0]'.\n";
        return($mem[0]);
    }

    if (scalar(@mem) > 1) {
        print STDERR "--\n";
        print STDERR "-- WARNING:  Couldn't determine the SGE resource to request memory.\n";
        print STDERR "-- WARNING:  Valid choices are:\n";
        print STDERR "-- WARNING:    $_\n"   foreach (@mem);
        print STDERR "--\n";

        return(undef);
    }

    return(undef);
}


sub configureSGE () {

    return   if (uc(getGlobal("gridEngine")) ne "SGE");

    my $maxArraySize = getGlobal("gridEngineArrayMaxJobs");

    if (!defined($maxArraySize)) {
        $maxArraySize = 65535;

        open(F, "qconf -sconf 2> /dev/null |") or caExit("can't run 'qconf' to get SGE config", undef);
        while (<F>) {
            if (m/max_aj_tasks\s+(\d+)/) {
                $maxArraySize = $1;
            }
        }
        close(F);
    }

    setGlobalIfUndef("gridEngineSubmitCommand",              "qsub");
    setGlobalIfUndef("gridEngineNameOption",                 "-cwd -N");
    setGlobalIfUndef("gridEngineArrayOption",                "-t ARRAY_JOBS");
    setGlobalIfUndef("gridEngineArrayName",                  "ARRAY_NAME");
    setGlobalIfUndef("gridEngineArrayMaxJobs",               $maxArraySize);
    setGlobalIfUndef("gridEngineOutputOption",               "-j y -o");
    setGlobalIfUndef("gridEngineResourceOption",             "");     #"-pe threads THREADS -l mem=MEMORY");
    #setGlobalIfUndef("gridEngineMemoryPerJob",              "0");    #  Do NOT set; default set below
    setGlobalIfUndef("gridEngineNameToJobIDCommand",         undef);
    setGlobalIfUndef("gridEngineNameToJobIDCommandNoArray",  undef);
    setGlobalIfUndef("gridEngineTaskID",                     "SGE_TASK_ID");
    setGlobalIfUndef("gridEngineArraySubmitID",              "\\\$TASK_ID");
    setGlobalIfUndef("gridEngineJobID",                      "JOB_ID");

    #  Try to figure out the name of the threaded job execution environment
    #  (it's the one with allocation_rule of $pe_slots),
    #  Try to figure out the name of the memory resource.

    if (!defined(getGlobal("gridEngineResourceOption"))) {
        my $thr = discoverThreadsOption();
        my $mem = discoverMemoryOption();

        if (!defined($thr) ||
            !defined($mem)) {
            print STDERR "--\n";
            print STDERR "-- ERROR:  Couldn't determine how to request resources.\n";
            print STDERR "-- ERROR:\n";
            print STDERR "-- ERROR:  Find an appropriate Parallel Environment (qconf -spl) and\n";
            print STDERR "-- ERROR:       an appropriate Memory Resource (qconf -sc) and set:\n";
            print STDERR "-- ERROR:\n";
            print STDERR "-- ERROR:    gridEngineResourceOption=\"-pe <pe_name> THREADS -l <mem_name>=MEMORY\"\n";
            print STDERR "-- ERROR:\n";
            print STDERR "-- ERROR:  Replace <pe_name> and <mem_name> with your choices.  Canu will\n";
            print STDERR "-- ERROR:  substitute the correct values for THREADS and MEMORY.\n";
            print STDERR "--\n";

            caExit("can't configure for SGE", undef);
        }

        setGlobal("gridEngineResourceOption", "-pe $thr THREADS -l $mem=MEMORY");
    }

    #  Otherwise, just log what the user requested.

    else {
        my $opt = getGlobal("gridEngineResourceOption");
        my $thr = $1  if ($opt =~ m/-pe\s+(.*)\s+THREADS/);
        my $mem = $1  if ($opt =~ m/-l\s+(.*)=MEMORY/);

        if (defined($thr)) {
            print STDERR "-- User supplied Parallel Environment '$thr'.\n";
        } else {
            print STDERR "-- No Sun Grid Engine parallel environment detected in gridEngineResourceOption.\n";
        }

        if (defined($mem)) {
            print STDERR "-- User supplied Memory Resource      '$mem'.\n";
        } else {
            print STDERR "-- No Sun Grid Engine memory resource detected in gridEngineResourceOption.\n";
        }
    }

    #  Check that the threads and memory options look sane.

    {
        my $opt = getGlobal("gridEngineResourceOption");
        my $thr = $1  if ($opt =~ m/-pe\s+(.*)\s+THREADS/);
        my $mem = $1  if ($opt =~ m/-l\s+(.*)=MEMORY/);

        if (defined($thr)) {
            my $ar = "\$pe_slots";   #  Default to $pe_slots state for Mick's busted grid config.
            my @lines;

            open(F, "qconf -sp $thr 2> /dev/null |");
            while (<F>) {
                push @lines, $_;

                $ar = $1  if (m/^allocation_rule\s+(.*)$/);
            }
            close(F);

            if ($ar ne "\$pe_slots") {
                print STDERR "\n";
                print STDERR "ERROR:  Sun Grid Engine parallel environment '$thr' is using 'allocation_rule' other than '\$pe_slots'.\n"  if ($ar == 0);
                print STDERR "ERROR:\n";

                if (scalar(@lines) == 0) {
                    print STDERR "ERROR:    (no output from qconf -sp $thr)\n";
                } else {
                    foreach my $l (@lines) {
                        print STDERR "ERROR:    $l";
                    }
                }

                caExit("can't configure for SGE", undef);
            }
        }

        #  Set memory-per-job if configured and the user hasn't supplied a value.
        #  Otherwise, leave it in the default off state.
        if (defined($mem) && !defined(getGlobal("gridEngineMemoryPerJob"))) {
            open(F, "qconf -sc 2> /dev/null |");
            while (<F>) {
                my @vals = split '\s+', $_;

                if (($vals[0] eq $mem) &&
                    ($vals[5] eq "JOB")) {
                    print STDERR "-- Enable memory-per-job mode.\n";
                    setGlobalIfUndef("gridEngineMemoryPerJob", "1");
                }
            }
            close(F);
        }
    }

    #  Check that SGE is setup to use the #! line instead of the (stupid) defaults.

    my %start_mode;
    my %start_shell;

    if (getGlobal('gridOptions') !~ m/-S/) {
        open(Q, "qconf -sql 2> /dev/null |");
        while (<Q>) {
            chomp;
            my $q = $_;

            $start_mode{$q}  = "na";
            $start_shell{$q} = "na";

            open(F, "qconf -sq $q 2> /dev/null |");
            while (<F>) {
                $start_mode{$q}  = $1   if (m/shell_start_mode\s+(\S+)/);
                $start_shell{$q} = $1   if (m/shell\s+(\S+)/);
            }
            close(F);
        }

        my $startBad = undef;

        foreach my $q (keys %start_mode) {
            if (($start_mode{$q}  ne "unix_behavior") &&
                ($start_shell{$q} =~ m/csh$/)) {
                $startBad .= "-- WARNING:  Queue '$q' has start mode set to 'posix_behavior' and shell set to '$start_shell{$q}'.\n";
            }
        }

        if (defined($startBad)) {
            my $bash = findCommand("bash");
            my $sh   = findCommand("sh");
            my $shell;

            $shell = $bash   if ($bash ne "");
            $shell = $sh     if ($sh   ne "");

            print STDERR "--\n";
            print STDERR "-- WARNING:\n";
            print STDERR "$startBad";
            print STDERR "-- WARNING:\n";
            print STDERR "-- WARNING:  Some queues in your configuration will fail to start jobs correctly.\n";
            print STDERR "-- WARNING:  Jobs will be submitted with option:\n";
            print STDERR "-- WARNING:    gridOptions=-S $shell\n";
            print STDERR "-- WARNING:\n";
            print STDERR "-- WARNING:  If jobs fail to start, modify the above option to use a valid shell\n";
            print STDERR "-- WARNING:  and supply it directly to canu.\n";
            print STDERR "-- WARNING:\n";

            if (!defined(getGlobal('gridOptions'))) {
                setGlobal('gridOptions', "-S $shell");
            } else {
                setGlobal('gridOptions', getGlobal('gridOptions') . " -S $shell");
            }
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

    if (scalar(keys(%hosts)) == 0) {
        my $mm = getGlobal("maxMemory");
        my $mt = getGlobal("maxThreads");

        print STDERR "--\n";
        print STDERR "-- WARNING:  No hosts found in 'qhost' report.\n";
        print STDERR "-- WARNING:  Will use maxMemory=$mm and maxThreads=$mt instead.\n";
        print STDERR "-- ERROR:    maxMemory not defined!\n"   if (!defined($mm));
        print STDERR "-- ERROR:    maxThreads not defined!\n"  if (!defined($mt));

        caExit("maxMemory or maxThreads not defined", undef)  if (!defined($mm) || !defined($mt));

        $hosts{"$mt-$mm"}++;
    }

    setGlobal("availableHosts", formatAllowedResources(%hosts, "Sun Grid Engine"));
}
