use strict;

sub filter {
    my $startTime = time();
    my $verbose   = "";
    my @ARGS      = @_;
    my $path      = "";
    my $type      = "";
    my $farmname;
    my $farmcode;
    my $filtqueue;
    my $extrafilter = "cat";

    #  Don't change the value without 2-search
    my $hitMemory = "600";


    print STDERR "ESTmapper: Performing a filter.\n";

    while (scalar @ARGS > 0) {
        my $arg = shift @ARGS;

        if ($arg eq "-filterest") {
            $path = shift @ARGS;
            $type  = "est";
        }
        if ($arg eq "-filtersnp") {
            $path = shift @ARGS;
            $type  = "snp";
        }
        if ($arg eq "-filtermrna") {
            $path = shift @ARGS;
            $type  = "mrna";
        }
        if ($arg eq "-filternone") {
            $type  = "none";
        }
        if ($arg eq "-hitsortmemory") {
            $hitMemory = shift @ARGS;
        }
        if ($arg eq "-verbose") {
            $verbose = "-verbose";
        }

        #  If present, this will be called after the real filtering
        #  gets done.  If -filternone, then this is the only filtering
        #  done.  Thanks to Alex Levitsky for suggesting this.
        #
        if ($arg eq "-extrahitfilter") {
            $extrafilter = shift @ARGS;
        }

        $farmname  = shift @ARGS if ($arg eq "-lsfjobname");
        $farmcode  = shift @ARGS if ($arg eq "-lsfproject");
        $filtqueue = shift @ARGS if ($arg eq "-lsffilterqueue");
    }

    ($path eq "")                         and die "FATAL ERROR: ESTmapper/filter-- No directory given.\n";
    (! -d "$path")                        and die "FATAL ERROR: ESTmapper/filter-- No directory '$path' found!\n";
    (! -f "$path/0-input/scaffolds-list") and die "FATAL ERROR: ESTmapper/filter-- No scaffolds-list?\n";
    (! -f "$path/1-search/allDone")       and die "FATAL ERROR: ESTmapper/filter-- The searches failed to complete successfully.\n";

    mkdir "$path/2-filter" if (! -d "$path/2-filter");


    #  If we're supposed to be running on LSF, but we aren't, restart.
    #  This can occur if the searches have finished, but the filter
    #  didn't, and we restart.  (also in 5-assemble.pl)
    #
    if (! -e "$path/2-filter/filteredHits" &&
        defined($filtqueue) &&
        !defined($ENV{'LSB_JOBID'})) {
        my $cmd;
        $cmd  = "bsub -q $filtqueue -P $farmcode -R \"select[mem>$hitMemory]rusage[mem=$hitMemory]\" -o $path/stage2.lsfout ";
        $cmd .= " -J f$farmname ";
        $cmd .= " $ESTmapper -restart $path";

        if (runCommand($cmd)) {
            die "ESTmapper/filter-- Failed to restart LSF execution.\n";
        }

        print STDERR "ESTmapper/filter-- Restarted LSF execution.\n";

        exit;
    }


    #  Merge all the hit counts into one list
    #
    if (! -e "$path/2-filter/hitCounts") {
        print STDERR "ESTmapper/filter-- Merging counts.\n";

        my $cmd = "$mergeCounts";

        open(F, "< $path/0-input/scaffolds-list");
        while (<F>) {
            chomp;
            $cmd .= " $path/1-search/$_.count";
        }
        close(F);

        $cmd .= "> $path/2-filter/hitCounts";

        if (runCommand($cmd)) {
            die "Failed.\n";
        }
    }


    #  Setup the filtering and sorting
    #
    if (! -e "$path/2-filter/filteredHits") {
        my $fcmd;
        my $scmd;

        #  The sort command used to depend on the type, but no more.
        #
        $scmd = "$sortHits $verbose -m $hitMemory -t $path/2-filter $path/2-filter/filtHits > $path/2-filter/filteredHits";

        if      ($type eq "est") {
            #  Original settings, but $filterEST was rewritten to fix
            #  a bug, and now these settings stink.
            #
            #$fcmd = "$filterEST -u 200 -r 200 -q 0.2 -log $path/2-filter/filterLog $path/1-search/*hits | $extrafilter > $path/2-filter/filtHits";

            #  bpw, 20051005, this isn't the perfect filter, but it
            #  does nearly as good as the best filter I've seen, and
            #  produces significantly fewer false positives.
            #
            $fcmd = "$filterEST -u 200000000000 -r 0 -log $path/2-filter/filterLog $path/1-search/*hits | $extrafilter > $path/2-filter/filtHits";
        } elsif ($type eq "snp") {
            $fcmd = "$filterMRNA $verbose -c $path/2-filter/hitCounts $path/1-search/*hits | $extrafilter > $path/2-filter/filtHits";
        } elsif ($type eq "mrna") {
            $fcmd = "$filterMRNA $verbose -c $path/2-filter/hitCounts $path/1-search/*hits | $extrafilter > $path/2-filter/filtHits";
        } elsif ($type eq "none") {
            $fcmd = "cat $path/1-search/*hits | $extrafilter > $path/2-filter/filtHits";
        }

        print STDERR "ESTmapper/filter-- Filtering.\n";
        if (runCommand($fcmd)) {
            unlink "$path/2-filter/filtHits";
            die "Failed.\n";
        }

        print STDERR "ESTmapper/filter-- Sorting.\n";
        if (runCommand($scmd)) {
            unlink "$path/2-filter/filteredHits";
            die "Failed.\n";
        }
    }

    die "ESTmapper/filter-- FATAL: filter and sort produced no hits?\n" if (-z "$path/2-filter/filteredHits");

    print STDERR "ESTmapper: Filter script finished in ", time() - $startTime, " wall-clock seconds.\n" if (time() > $startTime + 5);
}


1;
