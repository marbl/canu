use strict;

#  A wild over-estimation of the cost per base for the search.  Last
#  thing we want to do is thrash.
#
my $scaleFactor = 12;


sub configure_pack {
    my $path      = shift @_;
    my $segments  = shift @_;
    my $maxmemory = 0;
    my $segmentID = "000";

    open(L, "> $path/0-input/scaffolds-list");
    open(F, "$leaff -F $path/0-input/genomic.fasta --partition $segments |");
    $segments = <F>;
    while(<F>) {
        my $segments = "";
        my @pieces   = split '\s+', $_;
        my $memory   = shift @pieces;

        if ($memory =~ m/^\d+\]\((\d+)\)$/) {
            $memory = $1;
        } else {
            die "Error parsing memory for segment: $_\n";
        }

        foreach my $piece (@pieces) {
            if ($piece =~ m/(\d+)\(\d+\)/) {
                $segments .= "$1\n";
            } else {
                die "Error parsing segment: $piece\n";
            }
        }

        printf STDERR "ESTmapper/configure-- Created group $segmentID with $memory bases (%8.3fMB of memory).\n", $memory * $scaleFactor / 1024.0 / 1024.0;

        open(S, "> $path/0-input/scaffolds-$segmentID");
        print S $segments;
        close(S);

        print L "$segmentID\n";

        $maxmemory = $memory if ($memory > $maxmemory);

        $segmentID++;
    }
    close(F);
    close(L);

    return $scaleFactor * $maxmemory;
}


sub configure {
    my $startTime = time();
    my $errHdr    = "ERROR: ESTmapper/configure--";
    my @ARGS      = @_;
    my $path      = "";
    my $genomic   = "";
    my $memory    = 800;
    my $segments  = 0;

    print STDERR "ESTmapper: Performing a configure.\n";

    #
    #  Parse the args to find the path, then read any saved
    #  configuration, then reparse the args.
    #

    while (scalar @ARGS > 0) {
        my $arg = shift @ARGS;

        if ($arg eq "-configure") {
            $path = shift @ARGS;
        }
        if ($arg eq "-genomic") {
            $genomic = shift @ARGS;
        }
        if ($arg eq "-memory") {
            $memory   = shift @ARGS;
            $segments = 0;
        }
        if ($arg eq "-segments") {
            $memory   = 0;
            $segments = shift @ARGS;
        }
    }

    ($path eq "") and die "$errHdr no directory given.\n";
    ($genomic eq "") and die "$errHdr no genomic sequence given.\n";
    (! -f $genomic) and die "$errHdr can't find the genomic sequence '$genomic'\n";

    #  Make a place for us to work
    #
    system("mkdir $path")          if (! -d "$path");
    system("mkdir $path/0-input")  if (! -d "$path/0-input");
    system("mkdir $path/1-search") if (! -d "$path/1-search");
    system("mkdir $path/2-filter") if (! -d "$path/2-filter");
    system("mkdir $path/3-polish") if (! -d "$path/3-polish");


    #  Remember the genomic file for later
    #
    system("ln -s ${genomic}    $path/0-input/genomic.fasta")    if ((-e "${genomic}")    && (! -e "$path/0-input/genomic.fasta"));
    system("ln -s ${genomic}idx $path/0-input/genomic.fastaidx") if ((-e "${genomic}idx") && (! -e "$path/0-input/genomic.fastaidx"));

    if (! -f "$path/0-input/genomic.fasta") {
        die "$errHdr can't find the genomic sequence '$path/0-input/genomic.fasta'\n";
    }

    if (! -f "$path/0-input/genomic.fastaidx") {
        print STDERR "ESTmapper/configure-- Generating the index for '$path/0-input/genomic.fasta'\n";
        print STDERR "ESTmapper/configure-- WARNING:  This is done in the work directory!\n";
        system("$leaff -F $path/0-input/genomic.fasta");
    }


    #
    #  Partition the genome into itty-bitty pieces
    #
    if (! -e "$path/0-input/scaffolds-list") {
        if ($memory > 0) {
            print STDERR "ESTmapper/configure-- packing to preserve ${memory}MB memory limit\n";
            $memory /= $scaleFactor;
            $memory .= "mbp";
            $memory = configure_pack($path, $memory);
        }

        if ($segments > 0) {
            print STDERR "ESTmapper/configure-- packing to preserve $segments processor limit\n";
            $memory = configure_pack($path, $segments);
        }

        #  Argh!  Looks like int() is really trunc()
        $memory = int($memory / 1024.0 / 1024.0 + 1);

        open(F, "> $path/0-input/memoryLimit") or die "Can't write $path/0-input/memoryLimit\n";
        print F "$memory\n";
        close(F);

        print STDERR "ESTmapper/configure-- Created groups with maximum memory requirement of ${memory}MB.\n";
    }


    print STDERR "ESTmapper: Configure script finished in ", time() - $startTime, " wall-clock seconds.\n" if (time() > $startTime + 5);

    exit;
}


1;
