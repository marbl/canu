use strict;

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



1;
