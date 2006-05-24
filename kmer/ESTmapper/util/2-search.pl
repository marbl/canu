use strict;


sub submitFilter (@) {
    my $watch = join ",", @_;
    my $path  = $args{'path'};

    print STDERR "submitFilter\n";

    open(F, "> $path/1-search/filter-restart.sh");
    print F "#!/bin/sh\n";
    print F "/usr/bin/perl $prog{'ESTmapper'} -restart $path\n";
    close(F);

    my $cmd;
    $cmd  = "qsub -cwd -j y -o $path/stage2.sgeout ";
    $cmd .= " $args{'sgeoptions'} " if (defined($args{'sgeoptions'}));;
    $cmd .= " $args{'sgefilter'} "  if (defined($args{'sgefilter'}));
    $cmd .= " -N \"f$args{'sgename'}\" ";
    $cmd .= " -hold_jid $watch " if ($watch ne "");
    $cmd .= " $path/1-search/filter-restart.sh";

    die "Failed to submit job to SGE.\n" if (runCommand($cmd));
}


sub submitFinish (@) {
    my $watch = join ",", @_;
    my $path  = $args{'path'};

    print STDERR "submitFinish\n";

    open(F, "> $path/3-polish/finish-restart.sh");
    print F "#!/bin/sh\n";
    print F "/usr/bin/perl $prog{'ESTmapper'} -restart $path\n";
    close(F);

    my $cmd;
    $cmd  = "qsub -cwd -j y -o $path/stage3.sgeout ";
    $cmd .= " $args{'sgeoptions'} " if (defined($args{'sgeoptions'}));;
    $cmd .= " $args{'sgefinish'} "  if (defined($args{'sgefinish'}));
    $cmd .= " -N \"o$args{'sgename'}\" ";
    $cmd .= " -hold_jid $watch " if ($watch ne "");
    $cmd .= " $path/3-polish/finish-restart.sh";

    die "Failed to submit job to SGE.\n" if (runCommand($cmd));
}


sub search {
    my $startTime = time();
    my $path         = $args{'path'};

    #  If we're all done, just get outta here.
    return if (-e "$path/1-search/allDone");

    my $mersize      = ($args{'mersize'} or 20);
    my $merskip      = ($args{'merskip'} or 0);
    my $searchopts   = "";

    $searchopts = "-maxintron 2000000 -singlelength 20 -multiplelength 30 -smallsequence 100" if ($args{'runstyle'} eq "est");
    $searchopts = "-maxintron 2000000 -singlelength 30 -multiplelength 50 -smallsequence 0"  if ($args{'runstyle'} eq "mrna");
    $searchopts = "-maxintron 2000000 -singlecoverage 0.3 -multiplecoverage 0.3 -smallsequence 10000000 -extendminimum 100 -extendweight 2"  if ($args{'runstyle'} eq "snp");

    $searchopts .= $args{'searchopts'};

    my $numproc      = ($args{'localsearches'} or 4);
    my $numthread    = ($args{'searchthreads'} or 2);

    my $hitMemory    = ($args{'hitsortmemory'} or 600);    #  Don't change the value without 3-filter

    my $cdnaInInput  = int(`$prog{'leaff'} -F $path/0-input/cDNA.fasta -d`);


    #  Look for a mer masking file, or use the one supplied.
    #
    if (!defined($args{'maskmers'})) {
        $args{'maskmers'} = "$args{'genome'}/frequentMers-$mersize.fasta";
    }
    if (($args{'maskmers'} ne "none") && (! -e $args{'maskmers'})) {
        print STDERR "ESTmapper/search-- Can't find mer mask file '$args{'maskmers'}'.\n";
        print STDERR "ESTmapper/search-- Perhaps your genome isn't installed correctly?\n";
        print STDERR "ESTmapper/search-- Try a different mersize?\n";
        exit(1);
    }


    open(F, "< $path/0-input/memoryLimit");
    my $farmMemory = <F>;
    close(F);
    chomp $farmMemory;


    #  Create a bunch of scripts to process
    #
    #  Rewrite the command everytime.  This fixes the problem where
    #  we would, say, change the stats or number of threads...
    #
    open(F, "> $path/1-search/search.sh");
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
    print F "jid=`head -\$jid $path/0-input/genome/segments | tail -1`\n";
    print F "\n";
    print F "if [ -e \"$path/1-search/\$jid.success\" ] ; then\n";
    print F "  exit\n";
    print F "fi\n";
    print F "\n";
    print F "$prog{'seagen'} \\\n";
    print F "  -verbose \\\n"  if ($args{'verbose'});
    print F "  -binary \\\n";
    print F "  -mersize $mersize \\\n";
    print F "  -numthreads $numthread \\\n";
    print F "  $searchopts \\\n";
    print F "  -cdna      $path/0-input/cDNA.fasta \\\n";
    print F "  -genomic   $path/0-input/genome/genome.fasta \\\n";
    print F "  -positions $path/0-input/genome/seg\$jid.posDB \\\n";
    print F "  -mask      $args{'maskmers'} \\\n" if ($args{'maskmers'} ne "none");
    print F "  -output    $path/1-search/\$jid.hits \\\n";
    print F "  -count     $path/1-search/\$jid.count \\\n";
    print F "  -stats     $path/1-search/\$jid.stats \\\n" if ($args{'stats'} == 1);
    print F "&& \\\n";
    print F "touch $path/1-search/\$jid.success\n";
    close(F);

    chmod 0755, "$path/1-search/search.sh";


    #  Read the list of segments to figure out which segments we need to run.
    #
    my @searchesToRun;

    open(F, "< $path/0-input/genome/segments") or die "Can't open genome segments list!\n";
    while (<F>) {
        chomp;
        if (! -e "$path/1-search/$_.success") {
            print STDERR "ESTmapper/search-- search segment $_ not computed.\n";
            push @searchesToRun, $_;
        }
    }
    close(F);



    #  Run searches.  If the search terminated properly, the
    #  hit-counts file should exist.  Run (maybe re-run) the search if
    #  it isn't there.
    #
    if (defined($args{'runlater'})) {
        my $jobs = join " ", @searchesToRun;
        print STDERR "ESTmapper/search-- Please run the jobs:\n";
        print STDERR "ESTmapper/search--   $jobs\n";
        print STDERR "ESTmapper/search--  using:\n";
        print STDERR "ESTmapper/search--   $path/1-search/search.sh\n";
        exit(0);
    } elsif (defined($args{'sgename'})) {

        if (scalar(@searchesToRun) > 0) {
            print STDERR "ESTmapper/search-- SGE mode requested; ", scalar @searchesToRun, " processes to compute,\n";
            print STDERR "ESTmapper/search-- SGE mode requested; each with $numthread threads,\n";
            print STDERR "ESTmapper/search-- SGE mode requested; $farmMemory MB per process.\n";

            #  Don't resubmit jobs that are already done, and do
            #  submit the smallest number of jobs to finish.
            #  Bugs here should be fixed in 2-search.pl as well.

            my @watchJobs;

            my $fJob = shift @searchesToRun;
            my $lJob = $fJob;

            while (defined($lJob)) {
                my $nJob = shift @searchesToRun;

                if (($lJob + 1 != $nJob) || (!defined($nJob))) {

                    #  SGE expects jobs to start at 1, but we start at 0.
                    $fJob++;
                    $lJob++;

                    print STDERR "Sumbit $fJob - $lJob (njob=$nJob)\n";

                    my $cmd;
                    $cmd  = "qsub -cwd -j y -o $path/1-search/sgeout-\\\$TASK_ID ";
                    $cmd .= " $args{'sgeoptions'} " if (defined($args{'sgeoptions'}));
                    $cmd .= " $args{'sgesearch'} "  if (defined($args{'sgesearch'}));
                    $cmd .= " -N \"s$args{'sgename'}$fJob\" ";
                    $cmd .= " -t $fJob-$lJob ";
                    $cmd .= "$path/1-search/search.sh";

                    push @watchJobs, "s$args{'sgename'}$fJob";

                    die "Failed to submit job to SGE.\n" if (runCommand($cmd));

                    $fJob = $nJob;
                }
                $lJob = $nJob;
            }

            #  Submit the filter, and make it wait for the searches, if they were submitted.
            #
            submitFilter(@watchJobs);

            print STDERR "ESTmapper/search-- Searches submitted.   Rest of run is on the farm.\n";

            exit(0);
        }
    } else {
        print STDERR "ESTmapper/search-- Local mode requested; ", scalar @searchesToRun, " processes to compute,\n";
        print STDERR "ESTmapper/search-- Local mode requested; $numproc concurrent processes,\n";
        print STDERR "ESTmapper/search-- Local mode requested; each with $numthread threads.\n";

        #  Run the searches.  We use the scheduler, then check
        #  everything at the end.  This is a little less friendly
        #  to the user, but much easier for the implementor.
        #
        if (scalar(@searchesToRun) > 0) {
            &scheduler::schedulerSetNumberOfProcesses($numproc);
            &scheduler::schedulerSetShowCommands(0);
            &scheduler::schedulerSetShowStatus(1);
            foreach my $s (@searchesToRun) {
                print STDERR "sh $path/1-search/search.sh $s\n";
                &scheduler::schedulerSubmit("sh $path/1-search/search.sh $s");
            }
            &scheduler::schedulerFinish();
        }
    }


    #  See if anything failed.
    #
    print STDERR "ESTmapper/search-- checking search output.  All should have $cdnaInInput cDNA.\n";

    my $fail         = 0;

    open(F, "< $path/0-input/genome/segments") or die "Can't open genome segments list!\n";
    while (<F>) {
        chomp;

        #  If the hits file is NOT found, remove the count file.  Then
        #  figure out how many ESTs we have signals for, and fail if
        #  it's not what we expect.

        unlink "$path/1-search/$_.count" if (! -e "$path/1-search/$_.hits");

        my $c = int(`wc -l < $path/1-search/$_.count`) if (-e "$path/1-search/$_.count");

        if ($c != $cdnaInInput) {
            print STDERR "ESTmapper/search-- Search $_ failed, only $c signals.  Output saved as *.CRASH\n";
            rename "$path/1-search/$_.count", "$path/1-search/$_.count.CRASH";
            rename "$path/1-search/$_.hits",  "$path/1-search/$_.hits.CRASH";
            $fail++;
        }
    }
    close(F);

    die "Dang." if ($fail);

    #  Hooray!  Now we're all done!

    open(F, "> $path/1-search/allDone");
    close(F);

    print STDERR "ESTmapper/search-- Script finished in ", time() - $startTime, " wall-clock seconds.\n" if (time() > $startTime + 5);
}


1;
