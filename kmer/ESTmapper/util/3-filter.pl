use strict;

sub filter {
    my $startTime = time();

    #  If we're all done, just get outta here.
    return if (-e "$args{'path'}/2-filter/filteredHits");

    #  If we're supposed to be running on the grid, but we aren't, restart.
    #  This can occur if the searches have finished, but the filter
    #  didn't, and we restart.  (also in 5-assemble.pl)
    #
    if (defined($args{'sgename'}) && !defined($ENV{'SGE_TASK_ID'})) {
        submitFilter();
        print STDERR "ESTmapper/filter-- Restarted LSF execution.\n";
        exit;
    }

    my $path         = $args{'path'};
    my $verbose      = ($args{'verbose'}) ? "-verbose" : "";

    my $hitMemory    = ($args{'hitsortmemory'} or 600);    #  Don't change the value without 2-search

    print STDERR "ESTmapper: Performing a filter.\n";


    #  Merge all the hit counts into one list -- this is needed for output filtering!
    #
    if (! -e "$path/2-filter/hitCounts") {
        print STDERR "ESTmapper/filter-- Merging counts.\n";
        if (runCommand("$prog{'mergeCounts'} $path/1-search/[0-9]*[0-9].count > $path/2-filter/hitCounts")) {
            unlink "$path/2-filter/hitCounts";
            die "Failed.\n";
        }
    }

    #
    #  Setup the filtering and sorting
    #

    #  No verbose for filterNULL!
    #
    my $fcmd;

    #  bpw, 20051005, this isn't the perfect EST filter, but it does
    #  nearly as good as the best filter I've seen, and produces
    #  significantly fewer false positives.

    if      ($args{'nofilter'} eq 1) {
        $fcmd = "$prog{'filterNULL'} $path/1-search/*hits > $path/2-filter/filtHits";
    } elsif ($args{'runstyle'} eq "est") {
        $fcmd = "$prog{'filterEST'} -u 200000000000 -r 0 -log $path/2-filter/filterLog $path/1-search/*hits > $path/2-filter/filtHits";
    } elsif ($args{'runstyle'} eq "snp") {
        $fcmd = "$prog{'filterMRNA'} $verbose $path/1-search/*hits > $path/2-filter/filtHits";
    } elsif ($args{'runstyle'} eq "mrna") {
        $fcmd = "$prog{'filterMRNA'} $verbose $path/1-search/*hits > $path/2-filter/filtHits";
    } else {
        print STDERR "ESTmapper/filter-- nofilter = $args{'nofilter'}\n";
        print STDERR "ESTmapper/filter-- runstyle = $args{'runstyle'}\n";
        die "ESTmapper/filter-- Don't know how to filter!\n";
    }

    print STDERR "ESTmapper/filter-- Filtering.\n";
    if (runCommand($fcmd)) {
        unlink "$path/2-filter/filtHits";
        die "Failed.\n";
    }

    my $scmd = "$prog{'sortHits'} $verbose -m $hitMemory -t $path/2-filter $path/2-filter/filtHits > $path/2-filter/filteredHits";
    
    print STDERR "ESTmapper/filter-- Sorting.\n";
    if (runCommand($scmd)) {
        unlink "$path/2-filter/filteredHits";
        die "Failed.\n";
    }

    die "ESTmapper/filter-- FATAL: filter and sort produced no hits?\n" if (-z "$path/2-filter/filteredHits");

    print STDERR "ESTmapper: Filter script finished in ", time() - $startTime, " wall-clock seconds.\n" if (time() > $startTime + 5);
}


1;
