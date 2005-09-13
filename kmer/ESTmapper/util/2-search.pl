use strict;

sub checkFinished {
    my ($path, $cdnaInInput, @segList) = @_;
    my @searchesToRun;

    open(S, "> $path/1-search/run.sh");
    foreach my $s (@segList) {
        my $c = 0;

        #  If the hits file is NOT found, remove the count file.
        #
        unlink "$path/1-search/$s.count" if (! -e "$path/1-search/$s.hits");

        #  If there is a count file, count the number of lines in it.
        #
        $c = int(`wc -l < $path/1-search/$s.count`) if (-e "$path/1-search/$s.count");

        #  There should be one for each cDNA.
        #
        if ($c != $cdnaInInput) {
            print STDERR "ESTmapper/search-- Found signal for $c cDNA out of $cdnaInInput for segment $s.\n";
            print S "$path/1-search/$s.sh\n";
            push @searchesToRun, "$path/1-search/$s.sh";
        }
    }
    close(S);

    return(@searchesToRun);
}



sub search {
    my $startTime = time();
    my @ARGS         = @_;
    my $path         = "";
    my $opts         = "";
    my $mersize      = 20;
    my $maskFile     = undef;
    my $verbose      = "";
    my $stats        = 1;

    my $farmname     = undef;
    my $farmcode     = undef;
    my $farmqueue    = undef;
    my $filtqueue    = undef;

    my $sgename      = undef;
    my $sgeaccount   = undef;
    my $sgepriority  = undef;

    my $runtype      = "local";

    my $hitMemory    = 600;    #  Don't change the value without 3-filter

    my $numthread    = 2;
    my $numproc      = 4;
    my $species      = "";

    print STDERR "ESTmapper: Performing a search.\n";

    while (scalar @ARGS > 0) {
        my $arg = shift @ARGS;

        $verbose = "-verbose" if ($arg eq "-verbose");

        if ($arg eq "-searchest") {
            $path = shift @ARGS;
            $opts = "-maxintron 2000000 -singlelength 20 -multiplelength 30 -smallsequence 100";
        }
        if ($arg eq "-searchmrna") {
            $path = shift @ARGS;
            $opts = "-maxintron 2000000 -singlelength 30 -multiplelength 50 -smallsequence 0";
        }
        if ($arg eq "-searchsnp") {
            $path = shift @ARGS;
            $opts = "-maxintron 2000000 -singlecoverage 0.3 -multiplecoverage 0.3 -smallsequence 10000000 -extendminimum 100 -extendweight 2";
        }
        if ($arg eq "-searchopts") {
            $opts .= " " . shift @ARGS;
        }
        if ($arg eq "-mersize") {
            $mersize = shift @ARGS;
        }

        $stats = 1 if ($arg eq "-stats");
        $stats = 0 if ($arg eq "-nostats");

        $farmname   = shift @ARGS if ($arg eq "-lsfjobname");
        $farmcode   = shift @ARGS if ($arg eq "-lsfproject");
        $farmqueue  = shift @ARGS if ($arg eq "-lsfsearchqueue");
        $filtqueue  = shift @ARGS if ($arg eq "-lsffilterqueue");

        $sgename     = shift @ARGS         if ($arg eq "-sge");
        $sgeaccount  = "-A " . shift @ARGS if ($arg eq "-sgeaccount");
        $sgepriority = "-p " . shift @ARGS if ($arg eq "-sgepriority");

        $runtype    = "later"     if ($arg eq "-runlater");
        $runtype    = "lsf"       if (defined($farmname));
        $runtype    = "sge"       if (defined($sgename));

        $numproc   = shift @ARGS if ($arg eq "-localsearches");
        $numthread = shift @ARGS if ($arg eq "-searchthreads");

        $hitMemory = shift @ARGS if ($arg eq "-hitsortmemory");

        if ($arg eq "-species") {
            $species  = shift @ARGS;
            $maskFile = undef;
        }

        #  Use carefully or not at all.
        $maskFile = shift @ARGS if ($arg eq "-maskmers");
    }

    if (!defined($maskFile)) {
        if ($species eq "none") {
            $maskFile = "";
        } elsif (($species eq "human") ||
                 ($species eq "hum")) {
            $maskFile = "$exechome/data/frequentMers-hum-$mersize.fasta";
        } elsif (($species eq "mouse") ||
                 ($species eq "mus")) {
            $maskFile = "$exechome/data/frequentMers-mus-$mersize.fasta";
        } elsif (($species eq "rattus") ||
                 ($species eq "rat")) {
            $maskFile = "$exechome/data/frequentMers-rat-$mersize.fasta";
        } else {
            die "ESTmapper/search-- Unknown species '$species'.\n";
        }

        if (($maskFile ne "") && (! -e "$maskFile")) {
            print STDERR "ESTmapper/search-- Can't find mer mask file '$maskFile'.\n";
            print STDERR "ESTmapper/search-- Try a different mersize.\n";
            exit(1);
        }
    }


    #  If we're all done, skip the searches.
    #
    if (-e "$path/1-search/allDone") {
        print STDERR "ESTmapper/search-- this stage all done, and results checked previously.\n";
        return;
    }

    ($path eq "")                              and die "ERROR: ESTmapper/search-- no directory given.\n";
    (! -d "$path")                             and die "ERROR: ESTmapper/search-- no directory '$path' found!\n";

    mkdir "$path/1-search" if (! -d "$path/1-search");

    my $cdnaInInput  = int(`$leaff -F $path/0-input/cDNA.fasta -d`);

    #  Read the list of segments the searches used
    #
    open(F, "< $path/0-input/genome/segments") or die "Can't open genome segments list!\n";
    my @segList = <F>;
    close(F);
    chomp(@segList);



    #  Create a bunch of scripts to process
    #
    #  Rewrite the command everytime.  This fixes the problem where
    #  we would, say, change the stats or number of threads, and is
    #  needed for switching from -buildtables to -usetables.
    #
    foreach my $s (@segList) {

        open(F, "> $path/1-search/$s.sh");
        print F "$searchGENOME \\\n";
        print F " $verbose -binary -mersize $mersize $opts -numthreads $numthread \\\n";
        print F " -cdna $path/0-input/cDNA.fasta \\\n";
        print F " -genomic $path/0-input/genome/genome.fasta \\\n";
        print F " -positions $path/0-input/genome/seg$s.posDB \\\n";
        print F " -output $path/1-search/$s.hits \\\n";
        print F " -count $path/1-search/$s.count \\\n";
        print F " -stats $path/1-search/$s.stats \\\n" if ($stats == 1);
        print F " -mask $maskFile" if ($maskFile ne "");
        close(F);

        chmod 0755, "$path/1-search/$s.sh";
    }


    #  We build a list of jobs to run, and put it into
    #  $path/1-search/run.sh.
    #
    my @searchesToRun = checkFinished($path, $cdnaInInput, @segList);



    #  Run searches.  If the search terminated properly, the
    #  hit-counts file should exist.  Run (maybe re-run) the search if
    #  it isn't there.
    #
    if ($runtype eq "later") {
        print STDERR "ESTmapper/search-- Please run the jobs in\n";
        print STDERR "ESTmapper/search--   $path/1-search/run.sh\n";
        exit(0);
    } elsif ($runtype eq "local") {
        print STDERR "ESTmapper/search-- Local mode requested; ", scalar @segList, " processes to compute,\n";
        print STDERR "ESTmapper/search-- Local mode requested; $numproc concurrent processes,\n";
        print STDERR "ESTmapper/search-- Local mode requested; each with $numthread threads.\n";

        #  Run the searches.  We use the scheduler, then check
        #  everything at the end.  This is a little less friendly
        #  to the user, but much easier for the implementor.
        #
        if (scalar(@searchesToRun) > 0) {
            &scheduler::schedulerSetNumberOfProcesses($numproc);
            &scheduler::schedulerSetShowCommands(0);
            &scheduler::schedulerSetShowStatus(1);
            foreach my $s (@searchesToRun) {
                &scheduler::schedulerSubmit("sh $s");
            }
            &scheduler::schedulerFinish();

            #  See if anything failed.
            #
            print STDERR "ESTmapper/search-- checking search output.  All should have $cdnaInInput cDNA.\n";
            @searchesToRun = checkFinished($path, $cdnaInInput, @segList);
        }

        if (scalar(@searchesToRun) > 0) {
            print STDERR "ESTmapper/search-- Searches failed.\n";
            foreach my $s (@searchesToRun) {
                if ($s =~ m/(\d+).sh/) {
                    print STDERR "ESTmapper/search-- Search $1 failed.  Output saved as *.CRASH\n";
                    rename "$path/1-search/$1.count", "$path/1-search/$1.count.CRASH";
                    rename "$path/1-search/$1.hits", "$path/1-search/$1.hits.CRASH";
                } else {
                    print STDERR "ESTmapper/search-- Search $s failed.\n";
                }
            }
        }
    } elsif ($runtype eq "lsf") {
        open(F, "< $path/0-input/memoryLimit");
        my $farmMemory = <F>;
        close(F);
        chomp $farmMemory;

        my $cmd;
        my $jobid;

        if (scalar(@searchesToRun) > 0) {
            print STDERR "ESTmapper/search-- LSF mode requested; ", scalar @searchesToRun, " processes to compute,\n";
            print STDERR "ESTmapper/search-- LSF mode requested; each with $numthread threads,\n";
            print STDERR "ESTmapper/search-- LSF mode requested; $farmMemory MB per process.\n";

            $cmd  = "bsub -q $farmqueue -P $farmcode -n $numthread -R \"span[hosts=1]select[mem>$farmMemory]rusage[mem=$farmMemory]\" -o $path/1-search/%J-%I.lsfout ";
            $cmd .= " -J \"s$farmname\[1-" . scalar(@searchesToRun) . "\]\" ";
            $cmd .= "$exechome/util/jobarray.pl $path search";

            $jobid = runLSF($cmd);

            #  Submit the filter, and make it wait for the searches, if they were submitted.
            #
            $cmd  = "bsub -q $filtqueue -P $farmcode -R \"select[mem>$hitMemory]rusage[mem=$hitMemory]\" -o $path/stage2.lsfout ";
            $cmd .= " -w \"ended($jobid)\"" if (defined($jobid));
            $cmd .= " -J f$farmname ";
            $cmd .= " $ESTmapper -restart $path";

            if (runCommand($cmd)) {
                die "Failed to submit job to SGE.\n";
            }

            print STDERR "ESTmapper/search-- Searches submitted.   Rest of run is on the farm.\n";

            exit(0);
        }
    } elsif ($runtype eq "sge") {
        open(F, "< $path/0-input/memoryLimit");
        my $farmMemory = <F>;
        close(F);
        chomp $farmMemory;

        my $cmd;
        my $jobid;

        if (scalar(@searchesToRun) > 0) {
            print STDERR "ESTmapper/search-- SGE mode requested; ", scalar @searchesToRun, " processes to compute,\n";
            print STDERR "ESTmapper/search-- SGE mode requested; each with $numthread threads,\n";
            print STDERR "ESTmapper/search-- SGE mode requested; $farmMemory MB per process.\n";

            #  Still need to set number of threads, use parallel execution environment
            #  -p priority

            $cmd  = "qsub -cwd $sgeaccount $sgepriority ";
            $cmd .= " -pe thread $numthread ";
            $cmd .= " -j y -o $path/1-search/seagen-\\\$TASK_ID.sgeout ";
            $cmd .= " -N \"s$sgename\" ";
            $cmd .= " -t 1-" . scalar(@searchesToRun) . " ";
            $cmd .= "$exechome/util/jobarray.sh $exechome/util/jobarray.pl $path search";

            if (runCommand($cmd)) {
                die "Failed to submit job to SGE.\n";
            }

            #  Submit the filter, and make it wait for the searches, if they were submitted.
            #
            $cmd  = "qsub -cwd $sgeaccount $sgepriority ";
            $cmd .= " -j y -o $path/stage2.sgeout ";
            $cmd .= " -hold_jid \"s$sgename\" ";
            $cmd .= " -N \"f$sgename\" ";
            $cmd .= " $ESTmappersh $ESTmapper -restart $path";

            if (runCommand($cmd)) {
                die "Failed to submit job to SGE.\n";
            }

            print STDERR "ESTmapper/search-- Searches submitted.   Rest of run is on the farm.\n";

            exit(0);
        }
    } else {
        print STDERR "ESTmapper/search-- Unknown execution method $runtype.\n";
        exit(1);
    }



    #  See if anything is still broken.  If it did, holler!
    #
    my $died = 0;
    foreach my $s (@segList) {
        if (! -e "$path/1-search/$s.count") {
            print STDERR "ERROR: ESTmapper/search-- SEARCH $s DIED!\n";
            $died++;
        }
    }
    die "ERROR: ESTmapper/search-- Searches have a problem ($died).\n" if ($died > 0);



    #  Summarize compute time used.
    #
    my $sysTimeBuild  = 0;
    my $usrTimeBuild  = 0;
    my $sysTimeSearch = 0;
    my $usrTimeSearch = 0;
    foreach my $s (@segList) {
        if (-e "$path/1-search/$s.stats") {
            open(F, "< $path/1-search/$s.stats");
            while (<F>) {
                $sysTimeSearch += $1 if (m/^systemTime:\s+(\d+\.\d+)$/);
                $usrTimeSearch += $1 if (m/^userTime:\s+(\d+\.\d+)$/);
            }
            close(F);
        }
    }


    #  Rather lazy way to inform the next step (and future calls to this step) that we're all done.
    #
    open(F, "> $path/1-search/allDone");
    close(F);

    print STDERR "ESTmapper/search-- Used $sysTimeBuild seconds system time to precompute tables.\n";
    print STDERR "ESTmapper/search-- Used $usrTimeBuild seconds user time to precompute tables.\n";
    print STDERR "ESTmapper/search-- Used $sysTimeSearch seconds system time.\n";
    print STDERR "ESTmapper/search-- Used $usrTimeSearch seconds user time.\n";
    print STDERR "ESTmapper/search-- Script finished in ", time() - $startTime, " wall-clock seconds.\n" if (time() > $startTime + 5);
}


1;
