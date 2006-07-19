use strict;


sub setExecutables {
    my $exechome      = "$FindBin::Bin";

    $prog{'ESTmapper'}           = "$exechome/ESTmapper.pl";
    $prog{'seagen'}              = "$exechome/seagen";
    $prog{'mergeCounts'}         = "$exechome/mergeCounts";
    $prog{'filterEST'}           = "$exechome/filterEST";
    $prog{'filterMRNA'}          = "$exechome/filterMRNA";
    $prog{'filterNULL'}          = "$exechome/filterNULL";
    $prog{'sim4db'}              = "$exechome/sim4th";
    $prog{'leaff'}               = "$exechome/leaff";
    $prog{'cleanPolishes'}       = "$exechome/cleanPolishes";
    $prog{'toFILTER'}            = "$exechome/filterPolishes";
    $prog{'sortHits'}            = "$exechome/sortHits";
    $prog{'sortPolishes'}        = "$exechome/sortPolishes";
    $prog{'parseSNPs'}           = "$exechome/parseSNP";
    $prog{'pickBest'}            = "$exechome/pickBestPolish";
    $prog{'positionDB'}          = "$exechome/positionDB";
    $prog{'mersInMerStreamFile'} = "$exechome/mersInMerStreamFile";
    $prog{'terminate'}           = "$exechome/terminate";

    foreach my $e (keys %prog) {
        die "Can't find/execute $e ('$prog{$e}')\n" if (! -e $prog{$e});
    }
}


sub parseArgs (@) {
    my @ARGS = @_;

    $args{'scriptVersion'} = "10";
    $args{'startTime'}     = time();

    while (scalar(@ARGS) > 0) {
        my $arg = shift @ARGS;

        if      (($arg =~ m/^-dir/) || ($arg =~ m/^-path/)) {
            $args{'path'} = shift @ARGS;
        } elsif ($arg =~ m/^-genome/) {
            $args{'genome'} = shift @ARGS;

        } elsif (($arg =~ m/^-map(est)/) ||
                 ($arg =~ m/^-map(mrna)/) ||
                 ($arg =~ m/^-map(snp)/)) {
            $args{'runstyle'} = $1;
            $args{'queries'}  = shift @ARGS;
        } elsif ($arg =~ m/^-restart/) {
            $args{'runstyle'} = "restart";
            $args{'path'}  = shift @ARGS;
        } elsif ($arg =~ m/^-help/) {
            $args{'runstyle'} = "help";
        } elsif ($arg =~ m/^-time/) {
            $args{'runstyle'} = "time";

        } elsif ($arg =~ m/^-verbose/) {
            $args{'verbose'} = 1;
        } elsif ($arg =~ m/^-stats/) {
            $args{'stats'}   = 1;
        }


        #
        #  RUN options
        #
        elsif   ($arg =~ m/^-runlater/) {
            $args{'runlater'} = 1;
        }

        #
        #  LSF options
        #

        #
        #  SGE options
        #
        elsif ($arg =~ m/^-sge$/) {
            $args{'sgename'} = shift @ARGS;
        } elsif (($arg =~ m/^-(sgeoptions)/) ||
                 ($arg =~ m/^-(sgesearch)/)  ||
                 ($arg =~ m/^-(sgefilter)/)  ||
                 ($arg =~ m/^-(sgepolish)/)  ||
                 ($arg =~ m/^-(sgefinish)/)) {
            $args{$1} = shift @ARGS;
        }


        #
        #  search options
        #
        elsif  (($arg =~ m/^-(mersize)/)       ||
                ($arg =~ m/^-(merskip)/)       ||
                ($arg =~ m/^-(searchopts)/)    ||
                ($arg =~ m/^-(localsearches)/) ||
                ($arg =~ m/^-(searchthreads)/) ||
                ($arg =~ m/^-(hitsortmemory)/) ||
                ($arg =~ m/^-(maskmers)/)) {
            $args{$1}       = shift @ARGS;
        }

        #
        #  filter options
        #

        #
        #  polish options
        #
        elsif  (($arg =~ m/^-(mincoverage)/) ||
                ($arg =~ m/^-(minidentity)/) ||
                ($arg =~ m/^-(minlength)/) ||
                ($arg =~ m/^-(minsim4coverage)/) ||
                ($arg =~ m/^-(minsim4identity)/) ||
                ($arg =~ m/^-(minsim4length)/) ||
                ($arg =~ m/^-(relink)/) ||
                ($arg =~ m/^-(alwaysprint)/) ||
                ($arg =~ m/^-(batchsize)/) ||
                ($arg =~ m/^-(numbatches)/) ||
                ($arg =~ m/^-(localpolishes)/)) {
            $args{$1} = shift @ARGS;
        } elsif ($arg =~ m/^-interspecies/) {
            $args{'interspecies'} = 1;
        } elsif ($arg =~ m/^-aligns/) {
            $args{'aligns'} = 1;
        } elsif ($arg =~ m/^-noaligns/) {
            delete $args{'aligns'};
        } elsif ($arg =~ m/^-abort/) {
            $args{'abort'} = 1;
        } elsif ($arg =~ m/^-nofilter/) {
            $args{'nofilter'} = 1;
        } elsif ($arg =~ m/^-yn/) {
            $args{'nofilter'} = 1;
            $args{'sim4-yn'} = 1;
        }

        #
        #  finish options
        #
        elsif   ($arg =~ m/^-cleanup/) {
            $args{'cleanup'} = shift @ARGS;
        } elsif ($arg =~ m/^-nocleanup/) {
            delete $args{'cleanup'};
        } elsif ($arg =~ m/^-savetemporary/) {
            $args{'savetemporary'} = 1;
        }

        else {
            die "ESTmapper/configure-- unknown option '$arg'\n";
        }
    }

    #  Check we have a path!
    #
    ($args{'path'} eq "") and die "ERROR: ESTmapper/configure-- no directory given.\n";


    #  If told to restart, suck in the original configration, but
    #  don't overwrite things already defined.
    #
    if ($args{'runstyle'} eq "restart") {
        if (! -e "$args{'path'}/.runOptions") {
            print STDERR "ESTmapper/restart-- Nothing to restart!\n";
            exit;
        }

        delete $args{'runstyle'};

        open(F, "< $args{'path'}/.runOptions");
        while (<F>) {
            chomp;

            if (m/\s*(\S+)\s*=\s*(.*)\s*$/) {
                $args{$1} = $2 if (!defined($args{$1}));
            } else {
                die "Invalid runOption line '$_'\n";
            }
        }
        close(F);
    }


    foreach my $k (keys %args) {
        print STDERR "DEBUG $k=$args{$k}\n";
    }


    #  Make some organization
    #
    mkdir "$args{'path'}"          if (! -d "$args{'path'}");
    mkdir "$args{'path'}/0-input"  if (! -d "$args{'path'}/0-input");
    mkdir "$args{'path'}/1-search" if (! -d "$args{'path'}/1-search");
    mkdir "$args{'path'}/2-filter" if (! -d "$args{'path'}/2-filter");
    mkdir "$args{'path'}/3-polish" if (! -d "$args{'path'}/3-polish");


    #  Write the current set of args to the runOptions file
    #
    open(F, "> $args{'path'}/.runOptions") or die "Failed to open '$args{'path'}/.runOptions' to save options.\n";
    foreach my $k (keys %args) {
        print F "$k=$args{$k}\n";
    }
    close(F);
}




sub configure {
    my $path      = $args{'path'};

    print STDERR "ESTmapper: Performing a configure.\n";

    ($args{'genome'}  eq "") and die "ERROR: ESTmapper/configure-- no genomic sequences given.\n";
    ($args{'queries'} eq "") and die "ERROR: ESTmapper/configure-- no cDNA sequences given.\n";

    (! -f $args{'queries'})     and die "ERROR: ESTmapper/configure-- can't find the cdna sequence '$args{'queries'}'\n";

    #  XXX:  We should check that the genome dir is valid and complete.
    #
    symlink "$args{'genome'}", "$path/0-input/genome"  if (! -d "$path/0-input/genome");

    #  Check the input files exist, create symlinks to them, and find/build index files
    #
    symlink "$args{'queries'}",       "$path/0-input/cDNA.fasta"       if ((! -f "$path/0-input/cDNA.fasta"));
    symlink "$args{'queries'}idx",    "$path/0-input/cDNA.fastaidx"    if ((! -f "$path/0-input/cDNA.fastaidx") && (-f "$args{'queries'}idx"));

    if (! -f "$path/0-input/cDNA.fastaidx") {
        print STDERR "ESTmapper/configure-- Generating the index for '$path/0-input/cDNA.fasta'\n";
        runCommand("$prog{'leaff'} -F $path/0-input/cDNA.fasta") and die "Failed.\n";
    }

    #  Create a .runInformaiton file, containing supposedly useful information
    #  about this run.
    #
    my $time = time();
    $args{'runInfoFile'} = "$args{'path'}/.runInformation.$time";

    #  Write some information and the args to a run info file
    #
    open(F, "> $args{'runInfoFile'}");
    print F "startTime:  $time (", scalar(localtime($time)), ")\n";
    print F "operator:   $ENV{'USER'}\n";
    print F "host:       " . `uname -a`;
    print F "version:    $args{'scriptVersion'}\n";
    print F "parameters:";
    foreach my $k (keys %args) {
        print F "$k=$args{$k}\n";
    }
    close(F);

    unlink "$args{'path'}/.runInformation";
    symlink "$args{'path'}/.runInformation.$time", "$args{'path'}/.runInformation";

    print STDERR "ESTmapper: configured.\n";
}

1;
