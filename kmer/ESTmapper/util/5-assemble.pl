use strict;


#  Make sure that all the polishes are finished and OK.
#  Returns 0 if all done, 1 if not done.
#
sub polishesNotDone {
    my ($path) = @_;

    my $failed = 0;

    open(F, "< $path/3-polish/run-script");
    while (!eof(F)) {
        my $idx = <F>;  chomp $idx;
        my $cmd  = <F>;
        if (! -e "$path/3-polish/$idx.touch") {
            print STDERR "Polish $idx failed to finish.\n";
            $failed = 1;
        }
    }
    close(F);

    return $failed;
}


sub assembleOutput {
    my $startTime   = time();
    my $errHdr      = "";
    my @ARGS        = @_;
    my $path        = "";
    my $minc        = 50;
    my $mini        = 95;
    my $minl        = 0;
    my $intronLimit = 100000;
    my $deletetemp  = 1;

    my $farmname;
    my $farmcode;
    my $finiqueue;

    print STDERR "ESTmapper: Performing an assembleOutput.\n";

    while (scalar @ARGS > 0) {
        my $arg = shift @ARGS;

        $path        = shift @ARGS            if ($arg eq "-assembleoutput");
        $minc        = shift @ARGS            if ($arg eq "-mincoverage");
        $mini        = shift @ARGS            if ($arg eq "-minidentity");
        $minl        = shift @ARGS            if ($arg eq "-minlength");
        $intronLimit = int(shift @ARGS)       if ($arg eq "-cleanup");
        $intronLimit = 0                      if ($arg eq "-nocleanup");
        $deletetemp  = 0                      if ($arg eq "-savetemporary");

        $farmname  = shift @ARGS if ($arg eq "-lsfjobname");
        $farmcode  = shift @ARGS if ($arg eq "-lsfproject");
        $finiqueue = shift @ARGS if ($arg eq "-lsffinishqueue");
    }

    ($path eq "")                   and die "ERROR: ESTmapper/assembleOutput-- no directory given.\n";
    (! -d "$path")                  and die "ERROR: ESTmapper/assembleOutput-- no directory '$path' found!\n";
    (($mini < 0) || ($mini > 100))  and die "ERROR: ESTmapper/assembleOutput-- supply a value 0 <= x <= 100 for minidentity!\n";
    (($minc < 0) || ($minc > 100))  and die "ERROR: ESTmapper/assembleOutput-- supply a value 0 <= x <= 100 for mincoverage!\n";
    ($minl < 0)                     and die "ERROR: ESTmapper/assembleOutput-- supply a value x >= 0 for minlength!\n";

    (polishesNotDone($path) > 0) and die "There are unfinished polishing jobs.\n";


    #  If we're supposed to be running on LSF, but we aren't, restart.
    #  This can occur if the searches have finished, but the filter
    #  didn't, and we restart.  (also in 3-filter.pl)
    #
    if (((! -e "$path/summary") || (-z "$path/summary")) &&
        defined($finiqueue) &&
        !defined($ENV{'LSB_JOBID'})) {
        my $cmd;
        $cmd  = "bsub -q $finiqueue -P $farmcode -R \"select[mem>200]rusage[mem=200]\" -o $path/stage3.lsfout ";
        $cmd .= " -J \"o$farmname\" ";
        $cmd .= " $ESTmapper -restart $path";

        if (runCommand($cmd)) {
            die "ESTmapper/assembleOutput-- Failed to restart LSF execution.";
        }

        print STDERR "ESTmapper/assembleOutput-- Restarted LSF execution.\n";

        exit;
    }



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



    if (! -e "$path/polishes-good") {
        print STDERR "ESTmapper/assembleOutput-- filtering polishes by quality.\n";

        print STDERR "ESTmapper/assembleOutput-- identity:  percent align-sequence identity: $mini\n";
        print STDERR "ESTmapper/assembleOutput-- coverage:  percent query-sequence identity: $minc\n";
        print STDERR "ESTmapper/assembleOutput-- length:    length in bp of match:           $minl\n";

        #  Find all the polishes, run them through the cleaner, and filter by quality.
        #
        my $cmd;
        $cmd  = "find $path/3-polish/ -name '*.polished' -print | sort | xargs -n 100 cat | ";
        $cmd .= "$cleanPolishes -threshold $intronLimit -savejunk | " if ($intronLimit);
        $cmd .= "$toFILTER -c $minc -i $mini -l $minl -o $path/polishes-good -j $path/polishes-aborted > /dev/null";

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
        if      ($personality eq "-mapmrna") {
            print STDERR "ESTmapper/assembleOutput--  Picking the best mRNA polish.\n";
            if (runCommand("$sortPolishes -m 400 -c < $path/polishes-good | $pickBest -mrna > $path/polishes-best")) {
                unlink "$path/polishes-best";
                die "Failed.";
            }
        } elsif ($personality eq "-mapest") {
            print STDERR "ESTmapper/assembleOutput--  Picking the best EST polish.\n";
            if (runCommand("$sortPolishes -m 400 -c < $path/polishes-good | $pickBest -est > $path/polishes-best")) {
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

    if (! -e "$path/cDNA-missing.fasta") {
        unlink "$path/summary";

        #  For each polished file, figure out the ESTs that belong with those polishes.
        #  cDNA-good.fasta, cDNA-lowquality.fasta, cDNA-zero.fasta, cDNA-missing.fasta
        #
        print STDERR "ESTmapper/assembleOutput-- finding 'good' cDNA.\n";
        splitFastABasedOnPolishes("$path/0-input/cDNA.fasta",
                                  "$path/polishes-good",
                                  "$path/cDNA-good.fasta",
                                  "$path/cDNA-lost.fasta");

        if (-e "$path/2-filter/repeats") {
            print STDERR "ESTmapper/assembleOutput-- finding 'repeat' cDNA.\n";
            if (runCommand("$leaff -F $path/0-input/cDNA.fasta -q $path/2-filter/repeats > $path/cDNA-repeat.fasta")) {
                unlink "$path/cDNA-repeat.fasta";
                die "Failed.";
            }
        }

        print STDERR "ESTmapper/assembleOutput-- finding 'zero hit' cDNA.\n";
        copyZeroFastA("$path/0-input/cDNA.fasta",
                      "$path/2-filter/hitCounts",
                      "$path/cDNA-zero.fasta");

        #  subtractFastAfromFastA checks for input file existance, so we need
        #  not worry that cDNA-repeat.fasta might not exist
        #
        print STDERR "ESTmapper/assembleOutput-- finding 'missing' cDNA.\n";
        subtractFastAfromFastA("$path/cDNA-lost.fasta",
                               "$path/cDNA-repeat.fasta",
                               "$path/cDNA-zero.fasta",
                               "$path/cDNA-missing.fasta");

        #  Remove the temporary files
        #
        unlink "$path/cDNA-lost.fasta";
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
        my $cnttotl = int(`grep -c '>' $path/0-input/cDNA.fasta`);
        my $cntgood = int(`grep -c '>' $path/cDNA-good.fasta`);
        my $cntmiss = int(`grep -c '>' $path/cDNA-missing.fasta`);
        my $cntrept = int(`grep -c '>' $path/cDNA-repeat.fasta`) if (-e "$path/cDNA-repeat.fasta");
        my $cntzero = int(`grep -c '>' $path/cDNA-zero.fasta`);

        printf F "cDNA:            %8d\n", $cnttotl, "\n";
        printf F "cDNA-good:       %8d (%8.4f%%)\n", $cntgood, 100 * $cntgood / $cnttotl;
        printf F "cDNA-missing:    %8d (%8.4f%%)\n", $cntmiss, 100 * $cntmiss / $cnttotl;
        printf F "cDNA-repeat:     %8d (%8.4f%%)\n", $cntrept, 100 * $cntrept / $cnttotl  if (-e "$path/cDNA-repeat.fasta");
        printf F "cDNA-zero:       %8d (%8.4f%%)\n", $cntzero, 100 * $cntzero / $cnttotl;
    }


    #
    #  All done!
    #
    if ($deletetemp) {
        if (runCommand("rm -rf $path/1-search $path/2-filter $path/3-polish")) {
            print STDERR "ESTmapper/assembleOutput-- WARNING: Failed to remove temporary directories.\n";
        }
    }


    print STDERR "ESTmapper: assembleOutput script finished in ", time() - $startTime, " wall-clock seconds.\n" if (time() > $startTime + 5);
}


1;
