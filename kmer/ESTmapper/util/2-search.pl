sub countTheCountFile {
    my ($path, $s) = @_;
    my $c = 0;

    #  If the hits file is NOT found, remove the count file.
    #
    if (! -e "$path/1-search/$s.hits") {
        unlink "$path/1-search/$s.count";
    }

    #  If there is a count file, count the number of lines in it.
    #  There should be one for each cDNA.
    #
    if (-e "$path/1-search/$s.count") {
        $c = int(`wc -l < $path/1-search/$s.count`);
    }

    print STDERR "ESTmapper/search-- found $c cDNA for search $s\n";

    return($c);
}



sub search {
    my $startTime = time();
    my $errHdr       = "ERROR: ESTmapper/search--";
    my @ARGS         = @_;
    my $path         = "";
    my $opts         = "";
    my $cdna         = "";
    my $mersize      = 20;
    my $maskFile     = "";
    my $builddir     = undef;
    my $buildprefix  = undef;
    my $verbose      = "";
    my $stats        = 1;
    my $farm         = 0;
    my $farmqueue    = "";
    my $farmcode     = "";
    my $local        = 1;
    my $numthread    = 2;
    my $numproc      = 4;
    my $runnow       = 1;
    my $maxintron    = "-maxintron 2000000";
    my $species      = "";

    print STDERR "ESTmapper: Performing a search.\n";

    while (scalar @ARGS > 0) {
        my $arg = shift @ARGS;

        if ($arg eq "-searchest") {
            $path = shift @ARGS;
            $opts = "-singlelength 20 -multiplelength 30 -smallsequence 100";
        }
        if ($arg eq "-searchmrna") {
            $path = shift @ARGS;
            $opts = "-singlelength 30 -multiplelength 50 -smallsequence 0";
        }
        if ($arg eq "-searchsnp") {
            $path = shift @ARGS;
            $opts = "-singlecoverage 0.3 -multiplecoverage 0.3 -smallsequence 10000000 -extendminimum 100 -extendweight 2";
        }
        if ($arg eq "-searchopts") {
            $opts .= " " . shift @ARGS;
        }
        if ($arg eq "-cdna") {
            $cdna = shift @ARGS;
        }
        if ($arg eq "-mersize") {
            $mersize = shift @ARGS;
        }
        if ($arg eq "-maskmers") {
            $species  = "";
            $maskFile = shift @ARGS;
        }
        if ($arg eq "-prebuild") {
            $builddir    = shift @ARGS;
            $buildprefix = shift @ARGS;
        }
        if ($arg eq "-nomaskmers") {
            $maskFile = "";
        }
        if ($arg eq "-verbose") {
            $verbose = "-verbose";
        }
        if ($arg eq "-stats") {
            $stats = 1;
        }
        if ($arg eq "-nostats") {
            $stats = 0;
        }
        if ($arg eq "-farmsearches") {
            $farm     = 1;
            $local    = 0;
            $farmqueue = shift @ARGS;
            $farmcode  = shift @ARGS
        }
        if ($arg eq "-localsearches") {
            $farm  = 0;
            $local = 1;
            $numproc = shift @ARGS;
        }
        if ($arg eq "-searchthreads") {
            $numthread = shift @ARGS;
        }
        if ($arg eq "-runlater") {
            $runnow = 0;
        }
        if ($arg eq "-maxintron") {
            $maxintron = "-maxintron " . shift @ARGS;
        }
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
    }

    ($path eq "") and die "$errHdr no directory given.\n";
    (! -d "$path") and die "$errHdr no directory '$path' found!\n";
    (! -f "$path/0-input/scaffolds-list") and die "$errHdr no scaffolds-list?\n";
    (($maskFile ne "") && (! -e "$maskFile")) and die "$errHdr Can't find -maskmers '$maskFile'!\n";
    ($cdna eq "") and die "$errHdr Please supply some cDNA to search for!\n";
    (! -f $cdna) and die "$errHdr cDNA sequences '$cdna' not found!\n";


    if ($maskFile eq "") {
        print STDERR "ESTmapper/search-- No file of maskmers - mer masking disabled.\n";
    } else {
        print STDERR "ESTmapper/search-- Using '$maskFile' for mer masking.\n";
    }

    system("mkdir $path/1-search") if (! -d "$path/1-search");

    if (-e "$path/1-search/allDone") {
        print STDERR "ESTmapper/search-- this stage all done, and results checked previously.\n";
        return;
    }

    #  Save a pointer to the cDNA used
    #
    if (! -e "$path/0-input/cDNA.fasta") {
        system("ln -s $cdna      $path/0-input/cDNA.fasta");
    }

    if (! -e "$path/0-input/cDNA.fastaidx") {
        if (-f "${cdna}idx") {
            system("ln -s ${cdna}idx $path/0-input/cDNA.fastaidx");
        } else {
            print STDERR "Building index for cDNA sequences.\n";
            system("$leaff -F $path/0-input/cDNA.fasta");
        }
    }

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
        print F "$searchGENOME $verbose -binary -mersize $mersize $maxintron $opts -numthreads $numthread";
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

    my $cdnaInInput  = int(`$leaff -F $path/0-input/cDNA.fasta -d`);

    #  Run searches.  If the search terminated properly, the
    #  hit-counts file should exist.  Run (maybe re-run) the search if
    #  it isn't there.
    #
    if ($runnow) {
        if ($local) {
            print STDERR "ESTmapper/search-- Local mode requested; ", scalar @scafList, " processes to compute,\n";
            print STDERR "ESTmapper/search-- Local mode requested; $numproc concurrent processes,\n";
            print STDERR "ESTmapper/search-- Local mode requested; each with $numthread threads.\n";

            #  Build a list of things to run.
            #
            my @searchesToRun;

            foreach my $s (@scafList) {
                if (countTheCountFile($path, $s) != $cdnaInInput) {
                    push @searchesToRun, $s;
                }
            }

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
                    &libBri::schedulerSubmit("sh $path/1-search/$s.cmd");
                }
                &libBri::schedulerFinish();

                #  See if anything failed.
                #
                undef @searchesToRun;
                print STDERR "ESTmapper/search-- checking search output.  All should have $cdnaInInput cDNA.\n";
                foreach my $s (@scafList) {
                    if (countTheCountFile($path, $s) != $cdnaInInput) {
                        push @searchesToRun, $s;
                    }
                }
            }

            if (($numTries < 3) && (scalar(@searchesToRun) > 0)) {
                print STDERR "ESTmapper/search-- ", scalar(@searchesToRun), " searches failed.  Retrying.\n";
                goto again;
            }

            if (scalar(@searchesToRun) > 0) {
                print STDERR "ESTmapper/search-- Searches failed.\n";
                foreach my $s (@searchesToRun) {
                    print STDERR "ESTmapper/search-- Search $s failed.  Output saved as *.CRASH\n";
                    rename "$path/1-search/$s.count", "$path/1-search/$s.count.CRASH";
                    rename "$path/1-search/$s.hits", "$path/1-search/$s.hits.CRASH";
                }
            }
        }

        if ($farm) {
            open(F, "< $path/0-input/memoryLimit");
            my $farmMemory = <F>;
            close(F);
            chomp $farmMemory;

            print STDERR "ESTmapper/search-- Farm mode requested; ", scalar @scafList, " processes to compute,\n";
            print STDERR "ESTmapper/search-- Farm mode requested; each with $numthread threads,\n";
            print STDERR "ESTmapper/search-- Farm mode requested; $farmMemory MB per process.\n";

            for (my $x=10; $x>0; $x--) {
                print STDERR "ESTmapper/search-- Abort in $x seconds to change. \r";
                sleep(1);
            }
            print STDERR "ESTmapper/search-- STARTING FARM RUN!           \n";

            my $jobsToRun = 0;

            foreach my $s (@scafList) {
                if (countTheCountFile($path, $s) != $cdnaInInput) {
                    print STDERR "ESTmapper/search-- submitting search $s\n";
                    if (`uname` =~ m/aix/i) {
                        system("bsub -q $farmqueue -o $path/1-search/$s.stdout -R \"select[mem>$farmMemory]rusage[mem=$farmMemory]\" -P $farmcode $path/1-search/$s.cmd");
                    } else {
                        system("bsub -q $farmqueue -o $path/1-search/$s.stdout -R \"select[physmem>$farmMemory]rusage[physmem=$farmMemory]\" -P $farmcode $path/1-search/$s.cmd");
                    }
                    $jobsToRun++;
                } else {
                    print STDERR "ESTmapper/search-- search $s finished successfully!\n";
                }
            }

            #  Hang around waiting for them to finish.
            if ($jobsToRun > 0) {
                die "ESTmapper/search-- I don't know how to monitor LSF jobs!\nESTmapper/search-- Please restart when they're all done.\n";
            }
        }
    }


    #  See if anything is still broken.  If it did, holler!
    #  We also sum the time used here.
    #
    my $died    = 0;
    my $sysTime = 0;
    my $usrTime = 0;
    foreach my $s (@scafList) {
        if (! -e "$path/1-search/$s.count") {
            print STDERR "$errHdr SEARCH $s DIED!\n";
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
    die "Searches have a problem ($died).\n" if ($died > 0);


    #  Rather lazy way to inform the next step (and future calls to this step) that we're all done.
    #
    system("touch $path/1-search/allDone");

    print STDERR "ESTmapper: searchGENOME required $sysTime seconds system time.\n";
    print STDERR "ESTmapper: searchGENOME required $usrTime seconds user time.\n";
    print STDERR "ESTmapper: Search script finished in ", time() - $startTime, " wall-clock seconds.\n" if (time() > $startTime + 5);
}


1;
