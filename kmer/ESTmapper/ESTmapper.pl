#!/usr/bin/perl

#  Copyright (c) 2002 PE Corporation (NY) through the Celera Genomics Group
#  Copyright (c) 2003, 2004 Applied Biosystems
#  Copyright (c) 2004, 2005, 2006 Brian Walenz

$| = 1;

# Perl version 5.005_03 is too old, it requires two args to mkdir.

use strict;
use FindBin;
use Config;  #  for @signame
use lib "$FindBin::Bin/util";

use scheduler;

my %prog;
my %args;


################################################################################
#
#  Utility to run a command and check the exit status (sadly, duplicated
#  in configureESTmapper.pl).
#
################################################################################


sub runCommand {
    my $cmd = shift @_;

    print STDERR "$cmd\n";

    my $rc = 0xffff & system($cmd);

    #  Pretty much copied from Programming Perl page 230

    return(0) if ($rc == 0);

    #  Bunch of busy work to get the names of signals.  Is it really worth it?!
    #
    my @signame;
    if (defined($Config{sig_name})) {
        my $i = 0;
        foreach my $n (split('\s+', $Config{sig_name})) {
            $signame[$i] = $n;
            $i++;
        }
    }

    my $error = "ERROR: $cmd\n        failed with ";

    if ($rc == 0xff00) {
        $error .= "$!\n";
    } elsif ($rc > 0x80) {
        $rc >>= 8;
        $error .= "exit status $rc\n";
    } else {
        if ($rc & 0x80) {
            $rc &= ~0x80;
            $error .= "coredump from ";
        }
        if (defined($signame[$rc])) {
            $error .= "signal $signame[$rc]\n";
        } else {
            $error .= "signal $rc\n";
        }
    }

    print STDERR $error;

    return(1);
}


################################################################################
#
#  Command line parsing and configuration
#
################################################################################


sub setExecutables {
    my $exechome      = "$FindBin::Bin";

    $prog{'ESTmapper'}           = "$exechome/ESTmapper.pl";
    $prog{'seagen'}              = "$exechome/seagen";
    $prog{'mergeCounts'}         = "$exechome/mergeCounts";
    $prog{'filterEST'}           = "$exechome/filterEST";
    $prog{'filterMRNA'}          = "$exechome/filterMRNA";
    $prog{'filterNULL'}          = "$exechome/filterNULL";
    $prog{'sim4db'}              = "$exechome/sim4db";
    $prog{'leaff'}               = "$exechome/leaff";
    $prog{'meryl'}               = "$exechome/meryl";
    $prog{'cleanPolishes'}       = "$exechome/cleanPolishes";
    $prog{'toFILTER'}            = "$exechome/filterPolishes";
    $prog{'sortHits'}            = "$exechome/sortHits";
    $prog{'sortPolishes'}        = "$exechome/sortPolishes";
    $prog{'parseSNPs'}           = "$exechome/parseSNP";
    $prog{'pickBest'}            = "$exechome/pickBestPolish";
    $prog{'positionDB'}          = "$exechome/positionDB";
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

        if      (($arg =~ m/^-dir/) ||   #  depricated
                 ($arg =~ m/^-path/) ||  #  depricated
                 ($arg =~ m/^-outputdir/) ||
                 ($arg =~ m/^-mapdir/)) {
            $args{'path'} = shift @ARGS;
        } elsif (($arg =~ m/^-genomedir/) ||
                 ($arg =~ m/-genome/)) {  #  depricated
            $args{'genomedir'} = shift @ARGS;

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
        elsif  (($arg =~ m/^-(searchopts)/)    ||
                ($arg =~ m/^-(localsearches)/) ||
                ($arg =~ m/^-(searchthreads)/) ||
                ($arg =~ m/^-(hitsortmemory)/) ||
                ($arg =~ m/^-(mermaskfile)/) ||
                ($arg =~ m/^-(merignore)/)) {
            $args{$1}       = shift @ARGS;
        }

        #
        #  filter options
        #
        elsif  (($arg =~ m/^-(hitsortmemory)/)) {
            $args{$1}       = shift @ARGS;
        } elsif ($arg =~ m/^-nofilter/) {
            $args{'nofilter'} = 1;
        }

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

        #
        #  Are we installed correctly?
        #
        elsif   ($arg =~ m/-justtestingifitworks/) {
            exit(0);
        }

        else {
            die "ESTmapper/configure-- unknown option '$arg'\n";
        }
    }

    #  Check we have a path!
    #
    ($args{'path'} eq "") and die "ERROR: ESTmapper/configure-- no directory given.\n";

    #print STDERR "CONF $args{'genomedir'}\n";
    #print STDERR "CONF $args{'queries'}\n";
    #print STDERR "CONF $args{'path'}\n";


    #  Be tolerant of relative paths, but don't use them!
    #
    $args{'genomedir'} = "$ENV{'PWD'}/$args{'genomedir'}" if (defined($args{'genomedir'}) && ($args{'genomedir'} !~ m!^/!));
    $args{'queries'}   = "$ENV{'PWD'}/$args{'queries'}"   if (defined($args{'queries'})   && ($args{'queries'}   !~ m!^/!));
    $args{'path'}      = "$ENV{'PWD'}/$args{'path'}"      if (defined($args{'path'})      && ($args{'path'}      !~ m!^/!));


    #  Make some organization
    #
    mkdir "$args{'path'}"          if (! -d "$args{'path'}");
    mkdir "$args{'path'}/0-input"  if (! -d "$args{'path'}/0-input");
    mkdir "$args{'path'}/1-search" if (! -d "$args{'path'}/1-search");
    mkdir "$args{'path'}/2-filter" if (! -d "$args{'path'}/2-filter");
    mkdir "$args{'path'}/3-polish" if (! -d "$args{'path'}/3-polish");


    #  If told to restart, suck in the original configration, but
    #  don't overwrite things already defined.
    #
    if ($args{'runstyle'} eq "restart") {
        if (! -e "$args{'path'}/.runOptions") {
            print STDERR "ESTmapper/restart-- Nothing to restart!\n";
            exit;
        }

        delete $args{'runstyle'};

        open(F, "< $args{'path'}/.runOptions") or die "Failed to open '$args{'path'}/.runOptions' to read options.\n";
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

    #  Write the current set of args to the runOptions file
    #
    open(F, "> $args{'path'}/.runOptions") or die "Failed to open '$args{'path'}/.runOptions' to save options.\n";
    foreach my $k (keys %args) {
        #print STDERR "DEBUG $k=$args{$k}\n";
        print F "$k=$args{$k}\n";
    }
    close(F);
}


sub configure {
    my $path      = $args{'path'};

    print STDERR "ESTmapper: Performing a configure.\n";

    ($args{'genomedir'} eq "") and die "ERROR: ESTmapper/configure-- no genomic sequences given.\n";
    ($args{'queries'}   eq "") and die "ERROR: ESTmapper/configure-- no cDNA sequences given.\n";

    (! -f $args{'queries'})     and die "ERROR: ESTmapper/configure-- can't find the cdna sequence '$args{'queries'}'\n";

    #  XXX:  We should check that the genome dir is valid and complete.
    #
    symlink "$args{'genomedir'}", "$path/0-input/genome"  if (! -d "$path/0-input/genome");

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


################################################################################
#
#  Signal Finding
#
################################################################################


sub submitFilter (@) {
    my $watch = join ",", @_;
    my $path  = $args{'path'};

    open(F, "> $path/1-search/filter-restart.sh");
    print F "#!/bin/sh\n";
    print F "#\n";
    print F "#  Attempt to (re)configure SGE.  For reasons Bri doesn't know,\n";
    print F "#  jobs submitted to SGE, and running under SGE, fail to read his\n";
    print F "#  .tcshrc (or .bashrc, limited testing), and so they don't setup\n";
    print F "#  SGE (or ANY other paths, etc) properly.  For the record,\n";
    print F "#  interactive SGE logins (qlogin, etc) DO set the environment.\n";
    print F "#\n";
    print F ". \$SGE_ROOT/\$SGE_CELL/common/settings.sh\n";
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

    open(F, "> $path/3-polish/finish-restart.sh");
    print F "#!/bin/sh\n";
    print F "#\n";
    print F "#  Attempt to (re)configure SGE.  For reasons Bri doesn't know,\n";
    print F "#  jobs submitted to SGE, and running under SGE, fail to read his\n";
    print F "#  .tcshrc (or .bashrc, limited testing), and so they don't setup\n";
    print F "#  SGE (or ANY other paths, etc) properly.  For the record,\n";
    print F "#  interactive SGE logins (qlogin, etc) DO set the environment.\n";
    print F "#\n";
    print F ". \$SGE_ROOT/\$SGE_CELL/common/settings.sh\n";
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
    if (!defined($args{'mermaskfile'})) {
        $args{'merignore'} = 1000  if (!defined($args{'merignore'}));
        $args{'merignore'} = substr("000000$args{'merignore'}", -4);
        $args{'mermaskfile'} = "$args{'genomedir'}/frequentMers-ge$args{'merignore'}.fasta";
    }
    if (($args{'mermaskfile'} ne "none") && (! -e $args{'mermaskfile'})) {
        print STDERR "ESTmapper/search-- Didn't find mer mask file '$args{'mermaskfile'}', attempting\n";
        print STDERR "ESTmapper/search-- create it.\n";
        my $cmd;
        $cmd  = "$prog{'meryl'}";
        $cmd .= " -Dt -n $args{'merignore'} ";
        $cmd .= " -s \"$args{'genomedir'}//genome\"";
        $cmd .= " > \"$args{'genomedir'}/frequentMers-ge$args{'merignore'}.fasta\"";
        if (runCommand($cmd)) {
            die "ESTmapper/search-- Failed to create mask file.\n";
        }
    }
    if (($args{'mermaskfile'} ne "none") && (! -e $args{'mermaskfile'})) {
        print STDERR "ESTmapper/search-- Can't find mer mask file '$args{'mermaskfile'}'.\n";
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
    #  we would, say, change the number of threads...
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
    print F "  -genomic   $path/0-input/genome/genome.seqStore \\\n";
    print F "  -positions $path/0-input/genome/seg\$jid.posDB \\\n";
    print F "  -mask      $args{'mermaskfile'} \\\n" if ($args{'mermaskfile'} ne "none");
    print F "  -output    $path/1-search/\$jid.hits \\\n";
    print F "  -count     $path/1-search/\$jid.count \\\n";
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
                    $cmd .= " -N \"s$args{'sgename'}.$fJob\" ";
                    $cmd .= " -t $fJob-$lJob ";
                    $cmd .= "$path/1-search/search.sh";

                    push @watchJobs, "s$args{'sgename'}.$fJob";

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


################################################################################
#
#  Signal Filtering
#
################################################################################


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


################################################################################
#
#  Signal Polishing
#
################################################################################


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
    print F "  -genomic $path/0-input/genome/genome.seqStore \\\n";
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
                    $cmd .= " -N \"p$args{'sgename'}.$fJob\" ";
                    $cmd .= " -t $fJob-$lJob ";
                    $cmd .= "$path/3-polish/polish.sh";

                    push @watchJobs, "p$args{'sgename'}.$fJob";

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


################################################################################
#
#  Output
#
################################################################################


#  This is way too complicated.
#
#  1) Collect output from 4-polish, put into polishes-good
#  2) Filter -> polishes-best
#
#  Given as input a single polishes file and a cdna file,
#  we need an executable that: 
#    Generate stats on mapping, good and best, missing, zero
#    Filter cDNA to good, missing, zero


sub assembleOutput {
    my $startTime   = time();

    my $path        = $args{'path'};
    my $mini        = ($args{'minidentity'} or 95);
    my $minc        = ($args{'mincoverage'} or 50);
    my $minl        = ($args{'minlength'}   or 0);

    my $intronLimit = $args{'cleanup'} or 100000;

    print STDERR "ESTmapper: Performing an assembleOutput.\n";

    (($mini < 0) || ($mini > 100))  and die "ERROR: ESTmapper/assembleOutput-- supply a value 0 <= x <= 100 for minidentity!\n";
    (($minc < 0) || ($minc > 100))  and die "ERROR: ESTmapper/assembleOutput-- supply a value 0 <= x <= 100 for mincoverage!\n";
    ($minl < 0)                     and die "ERROR: ESTmapper/assembleOutput-- supply a value x >= 0 for minlength!\n";



    #  Check that the filtering is compatable with the polishing.
    #
    if (-e "$path/3-polish/parameters") {
        open(F, "< $path/3-polish/parameters");
        $_ = <F>;
        $_ = <F>;

        my $miniL = int(<F>);  #  Quality values used for last filtering
        my $mincL = int(<F>);
        my $minlL = int(<F>);

        my $miniP = int(<F>);  #  Quality values used for polishing
        my $mincP = int(<F>);
        my $minlP = int(<F>);
        close(F);

        if ($mini < $miniP) {
            printf STDERR "ESTmapper/assembleOutput-- WARNING: Percent identity quality level too low for existing polishing!\n";
            printf STDERR "ESTmapper/assembleOutput-- WARNING: Polished at percent align-sequence identity = %3d, requested filtration at %3d.\n", $miniP, $mini;
        }
        if ($minc < $mincP) {
            printf STDERR "ESTmapper/assembleOutput-- WARNING: Coverage quality level too low for existing polishing!\n";
            printf STDERR "ESTmapper/assembleOutput-- WARNING: Polished at percent query-sequence identity = %3d, requested filtration at %3d.\n", $mincP, $minc;
        }
        if ($minl < $minlP) {
            printf STDERR "ESTmapper/assembleOutput-- WARNING: Length quality level too low for existing polishing!\n";
            printf STDERR "ESTmapper/assembleOutput-- WARNING: Polished at length = %3d, requested filtration at %3d.\n", $minlP, $minl;
        }

        #  If the filter quality has changed, we need to refilter.  Nuke
        #  the filterLevel file, print a message.
        #
        if (($mini != $miniL) ||
            ($minc != $mincL) ||
            ($minl != $minlL)) {
            print STDERR "ESTmapper/assembleOutput-- filtering criteria changed; refiltering.\n";

            printf STDERR "ESTmapper/assembleOutput-- identity:  percent align-sequence identity: old=%3d new=%3d\n", $miniL, $mini;
            printf STDERR "ESTmapper/assembleOutput-- coverage:  percent query-sequence identity: old=%3d new=%3d\n", $mincL, $minc;
            printf STDERR "ESTmapper/assembleOutput-- length:    length in bp of match:           old=%3d new=%3d\n", $minlL, $minl;

            unlink "$path/polishes-good";
            unlink "$path/polishes-best";
            unlink "$path/polishes-lowquality";
            unlink "$path/summary";
        }
    } else {
        die "ESTmapper/assemblyOutput-- ERROR: Couldn't find polishing parameters.  Script error.\n";
    }



    #  If we're supposed to be running on LSF, but we aren't, restart.
    #  This can occur if the searches have finished, but the filter
    #  didn't, and we restart.  (also in 3-filter.pl)
    #
    if (defined($args{'sgename'}) && !defined($ENV{'SGE_TASK_ID'})) {
        submitFinish();
        print STDERR "ESTmapper/filter-- Restarted LSF execution.\n";
        exit;
    }



    if (! -e "$path/polishes-good") {
        print STDERR "ESTmapper/assembleOutput-- filtering polishes by quality.\n";

        print STDERR "ESTmapper/assembleOutput-- identity:  percent align-sequence identity: $mini\n";
        print STDERR "ESTmapper/assembleOutput-- coverage:  percent query-sequence identity: $minc\n";
        print STDERR "ESTmapper/assembleOutput-- length:    length in bp of match:           $minl\n";

        #  Find all the polishes, run them through the cleaner, and filter by quality.
        #
        my $cmd;
        $cmd  = "find $path/3-polish/ -name '*.sim4db' -print | sort | xargs -n 100 cat | ";
        $cmd .= "$prog{'cleanPolishes'} -threshold $intronLimit -savejunk | " if (defined($args{'cleanup'}));
        $cmd .= "$prog{'toFILTER'} -c $minc -i $mini -l $minl -o $path/polishes-good -j $path/polishes-aborted > /dev/null";

        if (runCommand($cmd)) {
            unlink "$path/polishes-good";
            unlink "$path/polishes-aborted";
            die "Failed.\n";
        }

        unlink "$path/polishes-best";
        unlink "$path/cDNA-good.fasta";
        unlink "$path/cDNA-missing.fasta";
        unlink "$path/cDNA-repeat.fasta";
        unlink "$path/cDNA-zero.fasta";
        unlink "$path/summary";
    }


    if (! -e "$path/polishes-best") {
        if      ($args{'runstyle'} eq "mrna") {
            print STDERR "ESTmapper/assembleOutput--  Picking the best mRNA polish.\n";
            if (runCommand("$prog{'sortPolishes'} -m 400 -c < $path/polishes-good | $prog{'pickBest'} -mrna > $path/polishes-best")) {
                unlink "$path/polishes-best";
                die "Failed.";
            }
        } elsif ($args{'runstyle'} eq "est") {
            print STDERR "ESTmapper/assembleOutput--  Picking the best EST polish.\n";
            if (runCommand("$prog{'sortPolishes'} -m 400 -c < $path/polishes-good | $prog{'pickBest'} -est > $path/polishes-best")) {
                unlink "$path/polishes-best";
                die "Failed.";
            }
        } else {
            print STDERR "ESTmapper/assembleOutput--  Not mRNA and not EST, so not picking the best polish.\n";
        }
    }

    #
    #  Segregate the sequences
    #

    #  XXXX  if the filter prints a list of repeats, we should add those here!

    if (! -e "$path/cDNA-good.fasta") {
        my $iid = 0;
        open(F, "< $path/2-filter/hitCounts");
        open(G, "> $path/zero-hit-iid");
        while (<F>) {
            if ($_ == 0) {
                print G "$iid\n";
            }
            $iid++;
        }
        close(G);
        close(F);

        my $cmd;
        $cmd  = "$prog{'terminate'}";
        $cmd .= " -P $path/polishes-best $path/cDNA-best.fasta";
        $cmd .= " -P $path/polishes-good $path/cDNA-good.fasta";
        $cmd .= " -I $path/zero-hit-iid  $path/cDNA-zero.fasta";
        $cmd .= " -O $path/cDNA-missing.fasta";
        $cmd .= " -i $path/0-input/cDNA.fasta";
        print $cmd;
        if (runCommand($cmd)) {
            rename "$path/cDNA-good.fasta", "$path/cDNA-good.fasta.FAILED";
            rename "$path/cDNA-missing.fasta", "$path/cDNA-missing.fasta.FAILED";
            rename "$path/cDNA-zero.fasta", "$path/cDNA-zero.fasta.FAILED";
            die "Failed.\n";
        }

        unlink "zero-hit-iid";
    }

    #
    #  Summarize
    #

    if ((! -e "$path/summary") || (-z "$path/summary")) {
        my ($mat, $est, $scf);

        open(F, "> $path/summary");

        print STDERR "ESTmapper/assembleOutput-- counting 'good' matches.\n";
        ($mat, $est, $scf) = summarizePolishes("$path/polishes-good");
        print F "GOOD: >= $mini% identity, >= $minc% composite, >= $minl bp\n"; 
        if ($mat > 0) {
            print F "cDNA-genomic matches  $mat matches ($est different cDNA and $scf genomic)\n";
            print F "Matches per cDNA      ", int(10000 * $mat / $est) / 10000.0, " matches/cDNA\n";
            print F "Matches per genomic   ", int(10000 * $mat / $scf) / 10000.0, " matches/genomic\n";
        } else {
            print F "cDNA-genomic matches  None.\n";
        }
        print F "\n";

        print STDERR "ESTmapper/assembleOutput-- counting cDNA.\n";
        print F "cDNA COUNTS:\n";
        my $cnttotl = int(`grep -c '^>' $path/0-input/cDNA.fasta`);
        my $cntgood = int(`grep -c '^>' $path/cDNA-good.fasta`);
        my $cntmiss = int(`grep -c '^>' $path/cDNA-missing.fasta`);
        my $cntrept = int(`grep -c '^>' $path/cDNA-repeat.fasta`) if (-e "$path/cDNA-repeat.fasta");
        my $cntzero = int(`grep -c '^>' $path/cDNA-zero.fasta`);

        printf F "cDNA:            %8d\n", $cnttotl, "\n";
        printf F "cDNA-good:       %8d (%8.4f%%)\n", $cntgood, 100 * $cntgood / $cnttotl;
        printf F "cDNA-missing:    %8d (%8.4f%%)\n", $cntmiss, 100 * $cntmiss / $cnttotl;
        printf F "cDNA-repeat:     %8d (%8.4f%%)\n", $cntrept, 100 * $cntrept / $cnttotl  if (-e "$path/cDNA-repeat.fasta");
        printf F "cDNA-zero:       %8d (%8.4f%%)\n", $cntzero, 100 * $cntzero / $cnttotl;
    }


    #
    #  All done!
    #
    if ($args{'savetemporary'} != 1) {
        if (runCommand("rm -rf $path/1-search $path/2-filter $path/3-polish")) {
            print STDERR "ESTmapper/assembleOutput-- WARNING: Failed to remove temporary directories.\n";
        }
    }


    print STDERR "ESTmapper: assembleOutput script finished in ", time() - $startTime, " wall-clock seconds.\n" if (time() > $startTime + 5);
}



######################################################################
#
#  Generates a report on a set of polishes.
#
#  number of cDNA-scaffold matches
#  number of different cDNA sequences in the set
#  number of different scaffolds in the set
#
sub summarizePolishes {
    my (@files) = @_;

    my %est;
    my %scf;
    my $mat = 0;
    my $ests = 0;
    my $scfs = 0;

    foreach my $infile (@files) {
        open(INPUT, "< $infile");

        while (<INPUT>) {
            if (m/^sim4begin$/) {
                $mat++;
            } elsif (m/^edef=/) {
                $ests++;
                $est{$_} = 1;
            } elsif (m/^ddef=/) {
                $scfs++;
                $scf{$_} = 1;
            }
        }

        close(INPUT);
    }

    if (($ests != $mat) || ($scfs != $mat)) {
        print STDERR "WARNING: summarizePolishes counted\n";
        print STDERR "           $mat matches\n";
        print STDERR "           $ests cDNA deflines\n";
        print STDERR "           $scfs scaffold deflines\n";
        print STDERR "         The number of deflines and the number of matches should be the same!\n";
    }

    return($mat, (scalar (keys %est)), (scalar (keys %scf)));
}


################################################################################
#
#  Utilities for Main
#
################################################################################


sub parseSNP {
    #  Parse the SNPs out
    #
    if (! -e "$args{'path'}/snps-parsed") {
        print STDERR "ESTmapper--  Parsing the SNPs\n";

        #  Sort, if needed.
        #
        if (! -e "$args{'path'}/polishes-good.sorted") {
            print STDERR "ESTmapper--  Sorting polishes by sequence ID; using 2GB memory maximum.\n";
            if (runCommand("$prog{'sortPolishes'} -m 2000 -c < $args{'path'}/polishes-good > $args{'path'}/polishes-good.sorted")) {
                unlink "$args{'path'}/polishes-good.sorted";
                die "Failed to sort the polishes.\n";
            }
        }

        #  Parse the options, looking for SNP specific ones
        #
        my @ARGS = @ARGV;
        my $snpdelimiter = "";
        my $snpsizetag   = "";
        my $snppostag    = "";
        my $snpoffset    = "";
        my $snpoutformat = "";

        while (scalar @ARGS > 0) {
            my $arg = shift @ARGS;

            if ($arg eq "-snpdelimiter") {
                $arg = shift @ARGS;
                $snpdelimiter = "-d \"$arg\"";
            } elsif ($arg eq "-snpsizetag") {
                $arg = shift @ARGS;
                $snpsizetag = "-s \"$arg\"";
            } elsif ($arg eq "-snppostag") {
                $arg = shift @ARGS;
                $snppostag = "-p \"$arg\"";
            } elsif ($arg eq "-snpoffset") {
                $arg = shift @ARGS;
                $snpoffset = "-o $arg";
            } elsif ($arg eq "-snpoutformat") {
                $arg = shift @ARGS;
                $snpoutformat = "-format $arg";
            }
        }

        #  PARSE!
        #
        if (runCommand("$prog{'parseSNPs'} $snpdelimiter $snpsizetag $snppostag $snpoffset $snpoutformat -F $args{'path'}/snps-failed -O $args{'path'}/snps-parsed < $args{'path'}/polishes-good.sorted > $args{'path'}/summary-snps")) {
            unlink "$args{'path'}/snps-failed";
            unlink "$args{'path'}/snps-parsed";
            unlink "$args{'path'}/summary-snps";
            die "Failed to parse SNP locations from polishes.\n";
        }
    }
}


sub sTOhms ($) {
    my ($s, $m, $h) = @_;
    $h = $s / 3600;
    $m = int(($h - int($h)) * 60);
    $h = int($h);
    $s = int($s);
    return($h,$m,$s);
}


################################################################################
#
#  Main
#
################################################################################


setExecutables();
parseArgs(@ARGV);

if      ($args{'runstyle'} eq "est") {
    configure();
    search();
    filter();
    polish();
    assembleOutput();
} elsif ($args{'runstyle'} eq "mrna") {
    $args{'relink'} = 1000;
    $args{'abort'}  = 1;

    configure();
    search();
    filter();
    polish();
    assembleOutput();
} elsif ($args{'runstyle'} eq "snp") {
    $args{'minidentity'} = 95;
    $args{'mincoverage'} = 80;

    configure();
    search();
    filter();
    polish();
    assembleOutput();
    parseSNP();
} else {
    print STDERR "Basic help N/A.\n";
}

print STDERR "ESTmapper: script finished everything in ", time() - $args{'startTime'}, " wall-clock seconds.\n" if (time() != $args{'startTime'});


if (-e $args{'runInforFile'}) {
    my $time = time();

    open(F, ">> $args{'runInforFile'}");
    print F "endTime:  $time (", scalar(localtime($time)), ")\n";
    close(F);
}

exit(0);
