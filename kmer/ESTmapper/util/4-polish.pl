use strict;

sub polish {
    my $startTime = time();
    my @ARGS         = @_;
    my $path         = "";
    my $mini         = "95";
    my $minc         = "50";
    my $minl         = "0";
    my $minsim4i     = "90";
    my $minsim4c     = "45";
    my $minsim4l     = "0";
    my $always       = "";
    my $relink       = "";
    my $batchsize    = 0;
    my $numbatches   = 256;
    my $farm         = 0;
    my $farmqueue    = "";
    my $farmcode     = "";
    my $local        = 1;
    my $numcpus      = 4;
    my $runnow       = 1;
    my $aligns       = "-align";
    my $stats        = 1;
    my $abort        = "";
    my $interspecies = "";

    print STDERR "ESTmapper: Performing a polish.\n";

    while (scalar @ARGS > 0) {
        my $arg = shift @ARGS;

        $path = shift @ARGS                      if ($arg eq "-polish");
        $minc = shift @ARGS                      if ($arg eq "-mincoverage");
        $mini = shift @ARGS                      if ($arg eq "-minidentity");
        $minl = shift @ARGS                      if ($arg eq "-minlength");
        $minsim4c = shift @ARGS                  if ($arg eq "-minsim4coverage");
        $minsim4i = shift @ARGS                  if ($arg eq "-minsim4identity");
        $minsim4l = shift @ARGS                  if ($arg eq "-minsim4length");
        $always = "-alwaysprint " . shift @ARGS  if ($arg eq "-alwaysprint");
        $relink = "-H " . shift @ARGS            if ($arg eq "-relink");
        $batchsize = int(shift @ARGS)            if ($arg eq "-batchsize");
        $numbatches = int(shift @ARGS)           if ($arg eq "-numbatches");
        $runnow = 0                              if ($arg eq "-runlater");
        $aligns = "-align"                       if ($arg eq "-aligns");
        $aligns = ""                             if ($arg eq "-noaligns");
        $stats = 1                               if ($arg eq "-stats");
        $stats = 0                               if ($arg eq "-nostats");
        $abort = "-Mp 0.25 -Ma 10000"            if ($arg eq "-abort");
        $interspecies = "-interspecies"          if ($arg eq "-interspecies");

        if ($arg eq "-farmpolishes") {
            $farm     = 1;
            $local    = 0;
            $farmqueue = shift @ARGS;
            $farmcode  = shift @ARGS
        }

        if ($arg eq "-localpolishes") {
            $farm  = 0;
            $local = 1;
            $numcpus = shift @ARGS;
        }
    }

    ($path eq "")   and die "ERROR: ESTmapper/polish-- no directory given.\n";
    (! -d "$path")  and die "ERROR: ESTmapper/polish-- no directory '$path' found!\n";

    system("mkdir $path/3-polish") if (! -d "$path/3-polish");


    #  Save the parameters, these are used on later invocations of
    #  polish, and in filter to make sure the user isn't an idiot.
    #
    if (-e "$path/3-polish/parameters") {
        print STDERR "ESTmapper/polish-- Using original parameters.\n";

        open(F, "< $path/3-polish/parameters");
        $numbatches   = int(<F>);
        $batchsize    = int(<F>);
        $mini         = <F>;      chomp $mini;
        $minc         = <F>;      chomp $minc;
        $minl         = <F>;      chomp $minl;
        $minsim4i     = <F>;      chomp $minsim4i;
        $minsim4c     = <F>;      chomp $minsim4c;
        $minsim4l     = <F>;      chomp $minsim4l;
        $relink       = <F>;      chomp $relink;
        $always       = <F>;      chomp $always;
        $aligns       = <F>;      chomp $aligns;
        $abort        = <F>;      chomp $abort;
        $interspecies = <F>;      chomp $interspecies;
        close(F);

        print STDERR "ESTmapper/polish-- Polish quality suitable for $minsim4i percent identity and\n";
        print STDERR "ESTmapper/polish--                             $minsim4c percent coverage\n";
        print STDERR "ESTmapper/polish-- To rerun polishes at a different quality level,\n";
        print STDERR "ESTmapper/polish-- remove the 3-polish directory.\n";
    } else {

        #  Do a little error checking; if both $batchsize and
        #  $numbatches are zero, set $batchsize to make 256 batches.
        #
        if (($batchsize == 0) && ($numbatches == 0)) {
            $numbatches = 256;
        }

        #  If $batchsize is not specified, compute it.
        #
        if ($batchsize == 0) {
            $batchsize = int(`wc -l < $path/2-filter/filteredHits` / $numbatches) + 1;
            $batchsize = 500 if ($batchsize < 500);
        }

        #  Adjust the sim4 qualities based on the final quality desired
        #
        $mini = 0 if ($mini < 0);
        $minc = 0 if ($minc < 0);
        $minl = 0 if ($minl < 0);

        $minsim4i = $mini - 5 if ($mini - 5 < $minsim4i);
        $minsim4c = $minc - 5 if ($minc - 5 < $minsim4c);
        $minsim4l = $minl     if ($minl     < $minsim4l);

        $minsim4i = 0 if ($minsim4i < 0);
        $minsim4c = 0 if ($minsim4c < 0);
        $minsim4l = 0 if ($minsim4l < 0);

        #  Save the parameters
        #
        open(F, "> $path/3-polish/parameters");
        print F "$numbatches\n$batchsize\n";
        print F "$mini\n$minc\n$minl\n";
        print F "$minsim4i\n$minsim4c\n$minsim4l\n";
        print F "$relink\n$always\n$aligns\n$abort\n$interspecies\n";
        close(F);
    }


    #  Display what parameters we are using
    #
    print STDERR "ESTmapper/polish--   minidentity   = $mini ($minsim4i)\n";
    print STDERR "ESTmapper/polish--   mincoverage   = $minc ($minsim4c)\n";
    print STDERR "ESTmapper/polish--   minlength     = $minl ($minsim4l)\n";
    print STDERR "ESTmapper/polish--   relink        = $relink\n";
    print STDERR "ESTmapper/polish--   always        = $always\n";
    print STDERR "ESTmapper/polish--   aligns        = $aligns\n";
    print STDERR "ESTmapper/polish--   abort         = $abort\n";
    print STDERR "ESTmapper/polish--   interspecies  = $interspecies\n";


    #  Splits the filteredHits into several pieces, and outputs a script
    #  that runs sim4db on those pieces.
    #
    if (! -e "$path/3-polish/run-script") {
        print STDERR "ESTmapper/polish-- Creating scripts with $batchsize lines in each.\n";

        my $idx = "0000";

        open(H, "< $path/2-filter/filteredHits");
        open(S, "> $path/3-polish/run-script");
        while (!eof(H)) {
            my $c = 0;

            open(F, "> $path/3-polish/$idx.scr");
            while (($c < $batchsize) && (!eof(H))) {
                $_ = <H>;
                print F $_;
                $c++;
            }
            close(F);

            #  The run-script is composed of three lines per command:
            #    the first line is the index of the run
            #    the second is the command
            #
            print S "$idx\n";
            if (`uname` =~ m/aix/i) {
                print S "bsub -q $farmqueue -o $path/3-polish/$idx.stdout -R \"select[mem>600]rusage[mem=500]\" -P $farmcode " if ($farm);
            } else {
                print S "bsub -q $farmqueue -o $path/3-polish/$idx.stdout -R \"select[physmem>600]rusage[physmem=500]\" -P $farmcode " if ($farm);
            }
            print S "$sim4db -cdna $path/0-input/cDNA.fasta -genomic $path/0-input/genomic.fasta ";
            print S "$aligns $always $relink $abort $interspecies -cut 0.6 ";
            print S "-mincoverage $minsim4c ";
            print S "-minidentity $minsim4i ";
            print S "-minlength $minsim4l ";
            print S "-script $path/3-polish/$idx.scr ";
            print S "-output $path/3-polish/$idx.polished ";
            print S "-stats  $path/3-polish/$idx.stats " if ($stats == 1);
            print S "-touch  $path/3-polish/$idx.touch\n";

            $idx++;
        }
        close(S);
        close(H);

        print STDERR "ESTmapper/polish-- Created $idx scripts.\n";
    }

    #  Builds a list of things to run by looking at the run-script, and
    #  seeing if the output for a given command exists.
    #
    my $polishesToPerform = 0;

    open(F, "< $path/3-polish/run-script");
    open(S, "> $path/3-polish/run.sh");
    while (!eof(F)) {
        my $idx = <F>;  chomp $idx;
        my $cmd = <F>;

        if (! -e "$path/3-polish/$idx.touch") {
            $polishesToPerform = 1;
            print S $cmd;
        }
    }
    close(S);
    close(F);

    #  Wipe any summaries, cDNA-* and polished files if we need to polish more stuff.
    #
    if ($polishesToPerform == 1) {
        print STDERR "ESTmapper/polish-- more polishes to compute - removing old output.\n";
        unlink "$path/cDNA-good.fasta";
        unlink "$path/cDNA-goodshort.fasta";
        unlink "$path/cDNA-lowquality.fasta";
        unlink "$path/cDNA-missing.fasta";
        unlink "$path/cDNA-repeat.fasta";
        unlink "$path/cDNA-zero.fasta";
        unlink "$path/polishes-aborted";
        unlink "$path/polishes-good";
        unlink "$path/polishes-goodshort";
        unlink "$path/polishes-lowquality";
        unlink "$path/summary";
    }


    #  Run things, or tell the user to do it for us.
    #
    if ($runnow) {
        if ($local && $polishesToPerform) {
            print STDERR "ESTmapper/polish-- Running locally, $numcpus at a time.\n";

            &libBri::schedulerSetNumberOfProcesses($numcpus);
            &libBri::schedulerSetShowCommands(0);
            &libBri::schedulerSetShowStatus(1);
            open(F, "< $path/3-polish/run.sh");
            while (<F>) {
                chomp;
                &libBri::schedulerSubmit($_);
            }
            close(F);
            &libBri::schedulerFinish();
        }

        if ($farm && $polishesToPerform) {
            print STDERR "ESTmapper/polish-- Submitting to the farm.\n";

            open(F, "< $path/3-polish/run.sh");
            while (<F>) {
                chomp;
                system($_);
            }
            close(F);
            
            #  Hang around waiting for them to finish.
            die "ESTmapper/polish-- I don't know how to monitor LSF jobs!\nESTmapper/polish-- Please restart when they're all done.\n";
        }

        unlink "$path/3-polish/run.sh";
    } else {
        print STDERR "ESTmapper/polish-- Please run the jobs in\n";
        print STDERR "ESTmapper/polish--   $path/3-polish/run.sh\n";
        exit(0);
    }

    #
    #  Summarize run-time performance of the jobs
    #
    my $clkTime = 0;
    my $sysTime = 0;
    my $usrTime = 0;
    open(F, "< $path/3-polish/run-script");
    while (!eof(F)) {
        my $idx = <F>;  chomp $idx;
        my $cmd = <F>;

        if (-e "$path/3-polish/$idx.stats") {
            open(X, "< $path/3-polish/$idx.stats");
            while (<X>) {
                $clkTime += $1 if (m/^clockTime:\s+(\d+\.\d+)$/);
                $sysTime += $1 if (m/^systemTime:\s+(\d+\.\d+)$/);
                $usrTime += $1 if (m/^userTime:\s+(\d+\.\d+)$/);
            }
            close(X);
        }
    }
    close(F);

    print STDERR "ESTmapper: sim4db required $clkTime seconds wall-clock time.\n";
    print STDERR "ESTmapper: sim4db required $sysTime seconds system time.\n";
    print STDERR "ESTmapper: sim4db required $usrTime seconds user time.\n";
    print STDERR "ESTmapper: Polish script finished in ", time() - $startTime, " wall-clock seconds.\n" if (time() > $startTime + 5);
}

1;
