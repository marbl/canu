use strict;

sub checkFinished {
    my ($path, $cdnaInInput, @scafList) = @_;
    my @searchesToRun;

    open(S, "> $path/1-search/run.sh");
    foreach my $s (@scafList) {
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

            print S "$path/1-search/$s.cmd";
            push @searchesToRun, "$path/1-search/$s.cmd";
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
    my $maskFile     = "";
    my $builddir     = undef;
    my $buildprefix  = undef;
    my $verbose      = "";
    my $stats        = 1;

    my $farmname;
    my $farmcode;
    my $farmqueue;
    my $filtqueue;

    my $hitMemory    = 600;    #  Don't change the value without 3-filter

    my $numthread    = 2;
    my $numproc      = 4;
    my $runnow       = 1;
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

        if ($arg eq "-prebuild") {
            $builddir    = shift @ARGS;
            $buildprefix = shift @ARGS;
        }

        $stats = 1 if ($arg eq "-stats");
        $stats = 0 if ($arg eq "-nostats");

        $farmname  = shift @ARGS if ($arg eq "-lsfjobname");
        $farmcode  = shift @ARGS if ($arg eq "-lsfproject");
        $farmqueue = shift @ARGS if ($arg eq "-lsfsearchqueue");
        $filtqueue = shift @ARGS if ($arg eq "-lsffilterqueue");

        $numproc   = shift @ARGS if ($arg eq "-localsearches");
        $numthread = shift @ARGS if ($arg eq "-searchthreads");

        $hitMemory = shift @ARGS if ($arg eq "-hitsortmemory");

        $runnow    = 0           if ($arg eq "-runlater");

        if ($arg eq "-species") {
            $species  = shift @ARGS;
            $maskFile = "";
            if ($species eq "none") {
                $maskFile = "";
            } elsif (($species eq "human") ||
                     ($species eq "hum")) {
                $maskFile = "$exechome/data/frequentMers-hum-20.fasta";
            } elsif (($species eq "mouse") ||
                     ($species eq "mus")) {
                $maskFile = "$exechome/data/frequentMers-mus-20.fasta";
            } elsif (($species eq "rattus") ||
                     ($species eq "rat")) {
                $maskFile = "$exechome/data/frequentMers-rat-20.fasta";
            } else {
                print STDERR "ESTmapper/search-- Unknown species '$species'.\n";
                $maskFile = "";
            }
        }

        #  Use carefully or not at all.
        $maskFile = shift @ARGS if ($arg eq "-maskmers");
    }

    #  If we're all done, skip the searches.
    #
    if (-e "$path/1-search/allDone") {
        print STDERR "ESTmapper/search-- this stage all done, and results checked previously.\n";
        return;
    }

    ($path eq "")                              and die "ERROR: ESTmapper/search-- no directory given.\n";
    (! -d "$path")                             and die "ERROR: ESTmapper/search-- no directory '$path' found!\n";
    (! -f "$path/0-input/scaffolds-list")      and die "ERROR: ESTmapper/search-- no scaffolds-list?\n";
    (($maskFile ne "") && (! -e "$maskFile"))  and die "ERROR: ESTmapper/search-- Can't find mer mask file '$maskFile'!\n";

    system("mkdir $path/1-search")                  if (! -d "$path/1-search");

    my $cdnaInInput  = int(`$leaff -F $path/0-input/cDNA.fasta -d`);

    #  Read the list of segments the searches used
    #
    open(F, "< $path/0-input/scaffolds-list");
    my @scafList = <F>;
    close(F);
    chomp @scafList;

    #  Create a bunch of scripts to process
    #
    #  Rewrite the command everytime.  This fixes the problem where
    #  we would, say, change the stats or number of threads.
    #
    foreach my $s (@scafList) {
        open(F, "> $path/1-search/$s.cmd");
        print F "$searchGENOME $verbose -binary -mersize $mersize $opts -numthreads $numthread";
        print F " -buildonly $builddir/$buildprefix.$mersize.$s" if (defined($buildprefix));
        print F " -cdna $path/0-input/cDNA.fasta";
        print F " -genomic $path/0-input/genomic.fasta";
        print F " -use $path/0-input/scaffolds-$s";
        print F " -output $path/1-search/$s.hits";
        print F " -count $path/1-search/$s.count";
        print F " -stats $path/1-search/$s.stats" if ($stats == 1);
        print F " -mask $maskFile" if ($maskFile ne "");
        close(F);
        system("chmod 755 $path/1-search/$s.cmd");
    }

    #  We build a list of jobs to run, and put it into
    #  $path/1-search/run.sh.
    #
    my @searchesToRun = checkFinished($path, $cdnaInInput, @scafList);

    #  Run searches.  If the search terminated properly, the
    #  hit-counts file should exist.  Run (maybe re-run) the search if
    #  it isn't there.
    #
    if ($runnow) {
        if (!defined($farmqueue)) {
            print STDERR "ESTmapper/search-- Local mode requested; ", scalar @scafList, " processes to compute,\n";
            print STDERR "ESTmapper/search-- Local mode requested; $numproc concurrent processes,\n";
            print STDERR "ESTmapper/search-- Local mode requested; each with $numthread threads.\n";

            my $numTries = 0;
          again:
            $numTries++;

            #  Run the searches.  We use the scheduler, then check
            #  everything at the end.  This is a little less friendly
            #  to the user, but much easier for the implementor.
            #
            if (scalar(@searchesToRun) > 0) {
                &libBri::schedulerSetNumberOfProcesses($numproc);
                &libBri::schedulerSetShowCommands(0);
                &libBri::schedulerSetShowStatus(1);
                foreach my $s (@searchesToRun) {
                    &libBri::schedulerSubmit("sh $s");
                }
                &libBri::schedulerFinish();

                #  See if anything failed.
                #
                print STDERR "ESTmapper/search-- checking search output.  All should have $cdnaInInput cDNA.\n";
                @searchesToRun = checkFinished($path, $cdnaInInput, @scafList);
            }

            if (($numTries < 3) && (scalar(@searchesToRun) > 0)) {
                print STDERR "ESTmapper/search-- ", scalar(@searchesToRun), " searches failed.  Retrying.\n";
                goto again;
            }

            if (scalar(@searchesToRun) > 0) {
                print STDERR "ESTmapper/search-- Searches failed.\n";
                foreach my $s (@searchesToRun) {
                    if ($s =~ m/(\d+).cmd/) {
                        print STDERR "ESTmapper/search-- Search $1 failed.  Output saved as *.CRASH\n";
                        rename "$path/1-search/$1.count", "$path/1-search/$1.count.CRASH";
                        rename "$path/1-search/$1.hits", "$path/1-search/$1.hits.CRASH";
                    } else {
                        print STDERR "ESTmapper/search-- Search $s failed.\n";
                    }
                }
            }
        } else {
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

                system($cmd);

                print STDERR "ESTmapper/search-- Searches submitted.   Rest of run is on the farm.\n";

                exit(0);
            }
        }
    } else {
        print STDERR "ESTmapper/search-- Please run the jobs in\n";
        print STDERR "ESTmapper/search--   $path/1-search/run.sh\n";
        exit(0);
    }


    #  See if anything is still broken.  If it did, holler!
    #  We also sum the time used here.
    #
    my $died    = 0;
    my $sysTime = 0;
    my $usrTime = 0;
    foreach my $s (@scafList) {
        if (! -e "$path/1-search/$s.count") {
            print STDERR "ERROR: ESTmapper/search-- SEARCH $s DIED!\n";
            $died++;
        }
        if (-e "$path/1-search/$s.stats") {
            open(F, "< $path/1-search/$s.stats");
            while (<F>) {
                if (m/^systemTime:\s+(\d+\.\d+)$/) {
                    $sysTime += $1;
                }
                if (m/^userTime:\s+(\d+\.\d+)$/) {
                    $usrTime += $1;
                }
            }
            close(F);
        }
    }
    die "ERROR: ESTmapper/search-- Searches have a problem ($died).\n" if ($died > 0);

    #  Rather lazy way to inform the next step (and future calls to this step) that we're all done.
    #
    system("touch $path/1-search/allDone");

    print STDERR "ESTmapper: searchGENOME required $sysTime seconds system time.\n";
    print STDERR "ESTmapper: searchGENOME required $usrTime seconds user time.\n";
    print STDERR "ESTmapper: Search script finished in ", time() - $startTime, " wall-clock seconds.\n" if (time() > $startTime + 5);
}


1;
