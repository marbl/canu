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

        $farmname  = shift @ARGS if ($arg eq "-lsfjobname");
        $farmcode  = shift @ARGS if ($arg eq "-lsfproject");
        $filtqueue = shift @ARGS if ($arg eq "-lsffilterqueue");
    }

    ($path eq "")                         and die "FATAL ERROR: ESTmapper/filter-- No directory given.\n";
    (! -d "$path")                        and die "FATAL ERROR: ESTmapper/filter-- No directory '$path' found!\n";
    (! -f "$path/0-input/scaffolds-list") and die "FATAL ERROR: ESTmapper/filter-- No scaffolds-list?\n";
    (! -f "$path/1-search/allDone")       and die "FATAL ERROR: ESTmapper/filter-- The searches failed to complete successfully.\n";

    system("mkdir $path/2-filter") if (! -d "$path/2-filter");



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

        print STDERR "ESTmapper/filter-- Restarted LSF execution.\n";

        system($cmd);

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

        system($cmd);
    }


    #  Setup the filtering and sorting
    #
    if (! -e "$path/2-filter/filteredHits") {
        my $fcmd;
        my $scmd;

        if      ($type eq "est") {
            $fcmd = "$filterEST -u 200 -r 200 -q 0.2 -log $path/2-filter/filterLog $path/1-search/*hits > $path/2-filter/filtHits";
            $scmd = "$sortHits $verbose -m $hitMemory -t $path/2-filter $path/2-filter/filtHits > $path/2-filter/filteredHits";
        } elsif ($type eq "snp") {
            $fcmd = "$filterMRNA $verbose -c $path/2-filter/hitCounts $path/1-search/*hits > $path/2-filter/filtHits";
            $scmd = "$sortHits $verbose -m $hitMemory -t $path/2-filter $path/2-filter/filtHits > $path/2-filter/filteredHits";
        } elsif ($type eq "mrna") {
            $fcmd = "$filterMRNA $verbose -c $path/2-filter/hitCounts $path/1-search/*hits > $path/2-filter/filtHits";
            $scmd = "$sortHits $verbose -m $hitMemory -t $path/2-filter $path/2-filter/filtHits > $path/2-filter/filteredHits";
        } elsif ($type eq "none") {
            $fcmd = "echo ESTmapper/filter-- No filter requested.";
            $scmd = "$sortHits $verbose -m $hitMemory -t $path/2-filter $path/1-search/*hits > $path/2-filter/filteredHits";
        }

        print STDERR "ESTmapper/filter-- Filtering.\n";
        system($fcmd) == 0 or die "FATAL ERROR: ESTmapper/filter-- Failed to filter!\n! = $!\n? = $?\n";

        print STDERR "ESTmapper/filter-- Sorting.\n";
        system($scmd) == 0 or die "FATAL ERROR: ESTmapper/filter-- Failed to sort!\n! = $!\n? = $?\n";
    }

    die "ESTmapper/filter-- FATAL: filter and sort produced no hits?\n" if (-z "$path/2-filter/filteredHits");

    print STDERR "ESTmapper: Filter script finished in ", time() - $startTime, " wall-clock seconds.\n" if (time() > $startTime + 5);
}


1;
