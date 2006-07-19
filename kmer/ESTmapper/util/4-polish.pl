use strict;

sub polish {
    my $startTime = time();

    #  If we're all done, just get outta here.
    return if (-e "$args{'path'}/3-polish/allDone");

    my $path         = $args{'path'};

    my $mini         = ($args{'minidentity'} or 95);
    my $minc         = ($args{'mincoverage'} or 50);
    my $minl         = ($args{'minlength'}   or 0);

    my $minsim4i     = ($args{'minsim4identity'} or 90);
    my $minsim4c     = ($args{'minsim4coverage'} or 45);
    my $minsim4l     = ($args{'minsim4length'}   or 0);

    my $relink       = "-H $args{'relink'}" if ($args{'relink'});
    my $always       = "-alwaysprint $args{'alwaysprint'}" if ($args{'alwaysprint'});

    my $batchsize    = ($args{'batchsize'} or 0);
    my $numbatches   = ($args{'numbatches'} or 256);

    my $numproc      = ($args{'localpolishes'} or 4);

    my $aligns       = "-aligns" if ($args{'aligns'});
    my $stats        = 1;
    my $abort        = "-Mp 0.25 -Ma 10000" if ($args{'abort'});
    my $interspecies = "-interspecies"      if ($args{'interspecies'});



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


    #  Build the sim4 command
    #
    open(F, "> $path/3-polish/polish.sh");
    print F "#!/bin/sh\n";
    print F "\n";
    print F "jid=\$SGE_TASK_ID\n";
    print F "if [ x\$jid = x -o x\$jid = xundefined ] ; then\n";
    print F "  if [ x\$1 = x ] ; then\n";
    print F "    echo \"ERROR: I need a job-id on the command line or in \$SGE_TASK_ID\"\n";
    print F "    exit 1\n";
    print F "  fi\n";
    print F "  jid=`expr \$1 + 1`\n";;
    print F "fi\n";
    print F "\n";
    print F "jid=`head -\$jid $path/3-polish/partitions | tail -1`\n";
    print F "\n";
    print F "if [ -e \"$path/3-polish/\$jid.success\" ] ; then\n";
    print F "  exit\n";
    print F "fi\n";
    print F "\n";
    print F "$prog{'sim4db'} \\\n";
    print F "  -cdna $path/0-input/cDNA.fasta \\\n";
    print F "  -genomic $path/0-input/genome/genome.fasta \\\n";
    print F "  $aligns \\\n"         if ($aligns ne "");
    print F "  $always \\\n"         if ($always ne "");
    print F "  $relink \\\n"         if ($relink ne "");
    print F "  $abort \\\n"          if ($abort  ne "");
    print F "  $interspecies \\\n"   if ($interspecies ne "");
    print F "  -cut 0.6 \\\n";
    print F "  -mincoverage $minsim4c \\\n";
    print F "  -minidentity $minsim4i \\\n";
    print F "  -minlength   $minsim4l \\\n";
    print F "  -script      $path/3-polish/\$jid.sim4script \\\n";
    print F "  -output      $path/3-polish/\$jid.sim4db \\\n";
    print F "  -stats       $path/3-polish/\$jid.stats \\\n" if ($stats == 1);
    print F "  -YN          $path/3-polish/\$jid.yn \\\n" if ($args{'sim4-yn'} == 1);
    print F "&& \\\n";
    print F "touch $path/3-polish/\$jid.success\n";
    close(F);


    #  Splits the filteredHits into several pieces, and outputs a script
    #  that runs sim4db on those pieces.
    #
    if (! -e "$path/3-polish/partitions") {
        print STDERR "ESTmapper/polish-- Creating scripts with $batchsize lines in each.\n";

        my @idxs;
        my $idx  = "0000";

        open(H, "< $path/2-filter/filteredHits");
        while (!eof(H)) {
            my $c = 0;

            open(F, "> $path/3-polish/$idx.sim4script");
            while (($c < $batchsize) && (!eof(H))) {
                $_ = <H>;
                print F $_;
                $c++;
            }
            close(F);

            push @idxs, "$idx\n";
            $idx++;
        }
        close(H);

        print STDERR "ESTmapper/polish-- Created $idx scripts.\n";

        open(S, "> $path/3-polish/partitions");
        print S @idxs;
        close(S);
    }


    #  Build a list of things to run.
    #
    my @jobsToRun;
    
    open(F, "< $path/3-polish/partitions");
    while (<F>) {
        chomp;
        push @jobsToRun, $_ if (! -e "$path/3-polish/$_.success");
    }
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
        if      (defined($args{'runlater'})) {
            print STDERR "ESTmapper/polish-- Please run the jobs in\n";
            print STDERR "ESTmapper/polish--   $path/3-polish/run.sh\n";
            exit(0);
        } elsif (defined($args{'sgename'})) {
            print STDERR "ESTmapper/polish-- Submitting to SGE.\n";

            #  Don't resubmit jobs that are already done, and do
            #  submit the smallest number of jobs to finish.
            #  Bugs here should be fixed in 2-search.pl as well.

            my @watchJobs;

            my $fJob = shift @jobsToRun;
            my $lJob = $fJob;

            while (defined($lJob)) {
                my $nJob = shift @jobsToRun;

                if (($lJob + 1 != $nJob) || (!defined($nJob))) {

                    #  SGE expects jobs to start at 1, but we start at 0.
                    $fJob++;
                    $lJob++;

                    print STDERR "Sumbit $fJob - $lJob (njob=$nJob)\n";

                    my $cmd;
                    $cmd  = "qsub -cwd -j y -o $path/3-polish/sgeout-\\\$TASK_ID ";
                    $cmd .= " $args{'sgeoptions'} " if (defined($args{'sgeoptions'}));;
                    $cmd .= " $args{'sgepolish'} "  if (defined($args{'sgepolish'}));
                    $cmd .= " -N \"p$args{'sgename'}$fJob\" ";
                    $cmd .= " -t $fJob-$lJob ";
                    $cmd .= "$path/3-polish/polish.sh";

                    push @watchJobs, "p$args{'sgename'}$fJob";

                    die "Failed to submit job to SGE.\n" if (runCommand($cmd));

                    $fJob = $nJob;
                }
                $lJob = $nJob;
            }

            submitFinish(@watchJobs);

            print STDERR "ESTmapper/polish-- Finish submitted.   See ya later!\n";

            exit(0);
        } else {
            print STDERR "ESTmapper/polish-- Running locally, $numproc at a time.\n";

            &scheduler::schedulerSetNumberOfProcesses($numproc);
            &scheduler::schedulerSetShowCommands(0);
            &scheduler::schedulerSetShowStatus(1);

            foreach my $cmd (@jobsToRun) {
                &scheduler::schedulerSubmit("/bin/sh $path/3-polish/polish.sh $cmd");
            }

            &scheduler::schedulerFinish();

            #unlink "$path/3-polish/run.sh";
        }
    }


    #  Make sure that all the polishes are finished and OK.
    #  If not, print dire warnings and exit.
    #
    my $fail = 0;

    open(F, "< $path/3-polish/partitions") or die "Failed to open '$path/3-polish/partitions'\n";;
    while (<F>) {
        chomp;
        if (! -e "$path/3-polish/$_.success") {
            $fail++;
            print STDERR "ESTmapper/polish-- segment $_ failed.\n";
        }
    }
    close(F);

    die "Dang." if ($fail);

    #  Hooray!  Now we're all done!

    open(F, "> $args{'path'}/3-polish/allDone");
    close(F);

    print STDERR "ESTmapper: Polish script finished in ", time() - $startTime, " wall-clock seconds.\n" if (time() > $startTime + 5);
}

1;
