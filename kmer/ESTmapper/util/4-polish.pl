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

    my $farmname     = undef;
    my $farmcode     = undef;
    my $farmqueue    = undef;
    my $finiqueue    = undef;

    my $sgename      = undef;
    my $sgeaccount   = undef;
    my $sgepriority  = undef;

    my $runtype      = "local";

    my $numproc      = 4;

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
        $aligns = "-align"                       if ($arg eq "-aligns");
        $aligns = ""                             if ($arg eq "-noaligns");

        $abort = "-Mp 0.25 -Ma 10000"            if ($arg eq "-abort");
        $interspecies = "-interspecies"          if ($arg eq "-interspecies");

        $stats = 1                               if ($arg eq "-stats");
        $stats = 0                               if ($arg eq "-nostats");

        $farmname  = shift @ARGS if ($arg eq "-lsfjobname");
        $farmcode  = shift @ARGS if ($arg eq "-lsfproject");
        $farmqueue = shift @ARGS if ($arg eq "-lsfpolishqueue");
        $finiqueue = shift @ARGS if ($arg eq "-lsffinishqueue");

        $sgename     = shift @ARGS         if ($arg eq "-sge");
        $sgeaccount  = "-A " . shift @ARGS if ($arg eq "-sgeaccount");
        $sgepriority = "-p " . shift @ARGS if ($arg eq "-sgepriority");

        $runtype    = "later"     if ($arg eq "-runlater");
        $runtype    = "lsf"       if (defined($farmname));
        $runtype    = "sge"       if (defined($sgename));

        $numproc   = shift @ARGS if ($arg eq "-localpolishes");
    }

    ($path eq "")   and die "ERROR: ESTmapper/polish-- no directory given.\n";
    (! -d "$path")  and die "ERROR: ESTmapper/polish-- no directory '$path' found!\n";

    mkdir "$path/3-polish" if (! -d "$path/3-polish");


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
            $batchsize = 10000 if ($batchsize < 10000);
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


    #  Splits the filteredHits into several pieces, and outputs a script
    #  that runs sim4db on those pieces.
    #
    if (! -e "$path/3-polish/splitDone") {
        print STDERR "ESTmapper/polish-- Creating scripts with $batchsize lines in each.\n";

        my @idxs;
        my $idx  = "0000";

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

            push @idxs, "$idx\n";
            $idx++;
        }
        close(S);
        close(H);

        print STDERR "ESTmapper/polish-- Created $idx scripts.\n";

        open(S, "> $path/3-polish/splitDone");
        print S @idxs;
        close(S);
    }

    #  Build a list of things to run.
    #
    my @jobsToRun;
    
    open(F, "< $path/3-polish/splitDone");
    open(S, "> $path/3-polish/run.sh");
    while (<F>) {
        my $idx = $_;
        chomp $idx;

        if (! -e "$path/3-polish/$idx.touch") {
            my $cmd;
            $cmd  = "$sim4db -cdna $path/0-input/cDNA.fasta -genomic $path/0-input/genome/genome.fasta ";
            $cmd .= "$aligns $always $relink $abort $interspecies -cut 0.6 ";
            $cmd .= "-mincoverage $minsim4c ";
            $cmd .= "-minidentity $minsim4i ";
            $cmd .= "-minlength $minsim4l ";
            $cmd .= "-script $path/3-polish/$idx.scr ";
            $cmd .= "-output $path/3-polish/$idx.polished ";
            $cmd .= "-stats  $path/3-polish/$idx.stats " if ($stats == 1);
            $cmd .= "-touch  $path/3-polish/$idx.touch";

            print S "$cmd\n";

            push @jobsToRun, $cmd;
        }
    }
    close(S);
    close(F);


    #  Wipe any summaries, cDNA-* and polished files if we need to polish more stuff.
    #
    if (scalar(@jobsToRun) > 0) {
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

        #  Display what parameters we are using
        #
        print STDERR "ESTmapper/polish-- more polishes to compute.\n";
        print STDERR "ESTmapper/polish--   minidentity   = $mini ($minsim4i)\n";
        print STDERR "ESTmapper/polish--   mincoverage   = $minc ($minsim4c)\n";
        print STDERR "ESTmapper/polish--   minlength     = $minl ($minsim4l)\n";
        print STDERR "ESTmapper/polish--   relink        = $relink\n";
        print STDERR "ESTmapper/polish--   always        = $always\n";
        print STDERR "ESTmapper/polish--   aligns        = $aligns\n";
        print STDERR "ESTmapper/polish--   abort         = $abort\n";
        print STDERR "ESTmapper/polish--   interspecies  = $interspecies\n";


        #  Run things, or tell the user to do it for us.
        #
        if      ($runtype eq "later") {
            print STDERR "ESTmapper/polish-- Please run the jobs in\n";
            print STDERR "ESTmapper/polish--   $path/3-polish/run.sh\n";
            exit(0);
        } elsif ($runtype eq "local") {
                print STDERR "ESTmapper/polish-- Running locally, $numproc at a time.\n";

                &scheduler::schedulerSetNumberOfProcesses($numproc);
                &scheduler::schedulerSetShowCommands(0);
                &scheduler::schedulerSetShowStatus(1);

                foreach my $cmd (@jobsToRun) {
                    &scheduler::schedulerSubmit($cmd);
                }

                &scheduler::schedulerFinish();

                unlink "$path/3-polish/run.sh";
        } elsif ($runtype eq "lsf") {
                print STDERR "ESTmapper/polish-- Submitting to LSF.\n";

                my $cmd;
                $cmd  = "bsub -q $farmqueue -P $farmcode -R \"select[mem>300]rusage[mem=600]\" -o $path/3-polish/%J-%I.lsfout ";
                $cmd .= " -J \"p$farmname\[1-" . scalar(@jobsToRun) . "\]\" ";
                $cmd .= "$exechome/util/jobarray.pl $path polish";

                my $jobid = runLSF($cmd);

                $cmd  = "bsub -q $finiqueue -P $farmcode -R \"select[mem>200]rusage[mem=200]\" -o $path/stage3.lsfout ";
                $cmd .= " -w \"ended($jobid)\"";
                $cmd .= " -J \"o$farmname\" ";
                $cmd .= " $ESTmapper -restart $path";

                if (runCommand($cmd)) {
                    die "Failed.\n";
                }

                print STDERR "ESTmapper/polish-- Finish submitted.   See ya later!\n";

                exit(0);
        } elsif ($runtype eq "sge") {
                print STDERR "ESTmapper/polish-- Submitting to SGE.\n";

                my $cmd;
                $cmd  = "qsub -cwd $sgeaccount $sgepriority ";
                $cmd .= " -pe thread 2 ";
                $cmd .= " -j y -o $path/3-polish/sim4db-\\\$TASK_ID.sgeout ";
                $cmd .= " -N \"p$sgename\" ";
                $cmd .= " -t 1-" . scalar(@jobsToRun) . " ";
                $cmd .= "$exechome/util/jobarray.sh $exechome/util/jobarray.pl $path polish";

                if (runCommand($cmd)) {
                    die "Failed.\n";
                }

                $cmd  = "qsub -cwd $sgeaccount $sgepriority ";
                $cmd .= " -j y -o $path/stage3.sgeout ";
                $cmd .= " -hold_jid \"p$sgename\" ";
                $cmd .= " -N \"o$sgename\" ";
                $cmd .= " $ESTmappersh $ESTmapper -restart $path";

                if (runCommand($cmd)) {
                    die "Failed.\n";
                }

                print STDERR "ESTmapper/polish-- Finish submitted.   See ya later!\n";

                exit(0);
        } else {
        }
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
