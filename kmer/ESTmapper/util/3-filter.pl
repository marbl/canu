use strict;

sub filter {
    my $startTime = time();
    my $errHdr    = "";
    my @ARGS      = @_;
    my $path      = "";
    my $type      = "";
    my $hitMemory = "600";
    my $verbose   = "";

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
        if ($arg eq "-skiphitfilter") {
            $type  = "none";
        }
        if ($arg eq "-hitsortmemory") {
            $hitMemory = shift @ARGS;
        }
        if ($arg eq "-verbose") {
            $verbose = "-verbose";
        }
    }

    ($path eq "")                         and die "FATAL ERROR: ESTmapper/filter-- No directory given.\n";
    (! -d "$path")                        and die "FATAL ERROR: ESTmapper/filter-- No directory '$path' found!\n";
    (! -f "$path/0-input/scaffolds-list") and die "FATAL ERROR: ESTmapper/filter-- No scaffolds-list?\n";
    (! -f "$path/1-search/allDone")       and die "FATAL ERROR: ESTmapper/filter-- The searches failed to complete successfully.\n";

    system("mkdir $path/2-filter") if (! -d "$path/2-filter");


    #  Read the list of segments the searches used
    #
    open(F, "< $path/0-input/scaffolds-list");
    my @scafList = <F>;
    close(F);
    chomp @scafList;


    #  Merge all the hit counts into one list
    #
    if (! -e "$path/2-filter/hitCounts") {
        print STDERR "ESTmapper/search-- Merging counts.\n";
        my $cmd = "$mergeCounts";
        foreach my $s (@scafList) {
            $cmd .= " $path/1-search/$s.count";
        }
        $cmd .= "> $path/2-filter/hitCounts";
        system($cmd);
    }

    #  Setup the filtering and sorting
    #
    if (! -e "$path/2-filter/filteredHits") {
        my $fcmd;
        my $scmd;

        if ($type eq "est") {
            $fcmd  = "$filterEST -u 200 -r 200 -q 0.2 -log $path/2-filter/filterLog $path/1-search/*hits > $path/2-filter/filtHits";
            $scmd = "$sortHits $verbose -m $hitMemory -t $path/2-filter $path/2-filter/filtHits > $path/2-filter/filteredHits";
        }

        if ($type eq "snp") {
            $fcmd  = "$filterMRNA $verbose -c $path/2-filter/hitCounts $path/1-search/*hits > $path/2-filter/filtHits";
            $scmd = "$sortHits $verbose -m $hitMemory -t $path/2-filter $path/2-filter/filtHits > $path/2-filter/filteredHits";
        }

        if ($type eq "mrna") {
            $fcmd  = "$filterMRNA $verbose -c $path/2-filter/hitCounts $path/1-search/*hits > $path/2-filter/filtHits";
            $scmd = "$sortHits $verbose -m $hitMemory -t $path/2-filter $path/2-filter/filtHits > $path/2-filter/filteredHits";
        }

        if ($type eq "none") {
            $fcmd  = "echo";
            $scmd = "$sortHits $verbose -m $hitMemory -t $path/2-filter $path/1-search/*hits > $path/2-filter/filteredHits";
        }

        print STDERR "ESTmapper/filter-- Filtering.\n";
        system($fcmd) == 0 or die "FATAL ERROR: ESTmapper/filter-- Failed to filter!\n! = $!\n? = $?\n";

        print STDERR "ESTmapper/filter-- Sorting.\n";
        system($scmd) == 0 or die "FATAL ERROR: ESTmapper/filter-- Failed to sort!\n! = $!\n? = $?\n";
    }

    die "FATAL ERROR: ESTmapper/filter-- filter and sort produced no hits?\n" if (-z "$path/2-filter/filteredHits");

    print STDERR "ESTmapper: Filter script finished in ", time() - $startTime, " wall-clock seconds.\n" if (time() > $startTime + 5);
}


1;
