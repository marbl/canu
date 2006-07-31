use strict;

sub overlapTrim {

    return if (getGlobal("doOverlapTrimming") == 0);

    goto alldone if (-e "$wrk/0-overlaptrim/overlaptrim.success");

    system("mkdir $wrk/0-overlaptrim")         if (! -d "$wrk/0-overlaptrim");
    system("mkdir $wrk/0-overlaptrim-overlap") if (! -d "$wrk/0-overlaptrim-overlap");

    #  We use this in a couple of places, so just make it more global.
    #
    my $im = getGlobal("immutableFrags");
    die "immutableFrags '$im' supplied, but not found!\n" if (defined($im) && (! -e $im));

    #  Do an initial overly-permissive quality trimming, intersected
    #  with any known vector trimming.
    #
    if ((! -e "$wrk/0-overlaptrim/$asm.initialTrimLog") &&
        (! -e "$wrk/0-overlaptrim/$asm.initialTrimLog.bz2")) {

        my $vi = getGlobal("vectorIntersect");
        my $im = getGlobal("immutableFrags");
        die "vectorIntersect '$vi' supplied, but not found!\n" if (defined($vi) && (! -e $vi));
        die "immutableFrags '$im' supplied, but not found!\n" if (defined($im) && (! -e $im));

        backupFragStore("beforeInitialTrim");

        my $cmd;
        $cmd  = "$bin/initialTrim -update -q 12 ";
        $cmd .= " -vector $vi "    if (defined($vi));
        $cmd .= " -immutable $im " if (defined($im));
        $cmd .= " -log $wrk/0-overlaptrim/$asm.initialTrimLog ";
        $cmd .= " -frg $wrk/$asm.frgStore ";
        $cmd .= " > $wrk/0-overlaptrim/initialTrim.err 2>&1";

        if (runCommand($cmd)) {
            rename "$wrk/0-overlaptrim/$asm.initialTrimLog", "$wrk/0-overlaptrim/$asm.initialTrimLog.failed";
            die "Failed.\n";
        }
    }


    #  Compute overlaps, if we don't have them already

    if ((! -e "$wrk/0-overlaptrim/$asm.ovl.sorted") &&
        (! -e "$wrk/0-overlaptrim/$asm.ovl.sorted.bz2")) {

        #  Run meryl on the now quality-trimmed frags.  This hopefully
        #  will get around some sequencing centers' habit of having N's in
        #  the low quality region, and is generally just a good idea.
        #
        meryl();

        #  Filter the standard set of nmers, throw out things below 100.
        #  If you change 100, you should also change meryl.pl.

        if (! -e "$wrk/0-overlaptrim-overlap/$asm.nmers.fasta") {
            open(F, "< $wrk/0-preoverlap/$asm.nmers.fasta")  or die "Failed to open $wrk/0-preoverlap/$asm.nmers.fasta for reading.\n";
            open(G, "> $wrk/0-overlaptrim-overlap/$asm.nmers.fasta") or die "Failed to open $wrk/0-overlaptrim-overlap/$asm.nmers.fasta for writing.\n";
            while (!eof(F)) {
                my $def = <F>;
                my $mer = <F>;
                if ($def =~ m/^>(\d+)$/) {
                    print G "$def$mer" if ($1 > 100);
                } else {
                    chomp $def;
                    print STDERR "ERROR:  Got '$def' for a defline!\n";
                }
            }
            close(G);
            close(F);
        }

        createOverlapJobs("trim");
        checkOverlap("trim");

        #  Sort the overlaps -- this also duplicates each overlap so that
        #  all overlaps for a fragment A are localized.

        if (runCommand("find $wrk/0-overlaptrim-overlap -follow -name \\*ovb -print > $wrk/0-overlaptrim/all-overlaps-trim.ovllist")) {
            print STDERR "Failed to generate a list of all the overlap files.\n";
            exit(1);
        }

        my $cmd;
        $cmd  = "cd $wrk && ";
        $cmd .= "$bin/sort-overlaps";
        $cmd .= " -memory " . getGlobal('ovlSortMemory') . " ";
        $cmd .= " -maxiid $numFrags ";
        $cmd .= " -L $wrk/0-overlaptrim/all-overlaps-trim.ovllist";
        $cmd .= " > $wrk/0-overlaptrim/$asm.ovl.sorted";
        $cmd .= " 2> $wrk/0-overlaptrim/$asm.ovl.sorted.err";

        if (runCommand($cmd)) {
            unlink "$wrk/0-overlaptrim/$asm.ovl.sorted";
            die "Failed to sort.\n";
        }
    }

    #  Consolidate the overlaps, listing all overlaps for a single
    #  fragment on a single line.  These are still iid's.

    if ((! -e "$wrk/0-overlaptrim/$asm.ovl.consolidated") &&
        (! -e "$wrk/0-overlaptrim/$asm.ovl.consolidated.bz2")) {

        if (runCommand("$bin/consolidate < $wrk/0-overlaptrim/$asm.ovl.sorted > $wrk/0-overlaptrim/$asm.ovl.consolidated")) {
          unlink "$wrk/0-overlaptrim/$asm.ovl.consolidated";
          die "Failed to consolidate.\n";
        }
    }


    #  We need to have all the overlaps squashed already, in particular so
    #  that we can get the mode of the 5'mode.  We could do this all in
    #  core, but that would take lots of space.

    if ((! -e "$wrk/0-overlaptrim/$asm.mergeLog") &&
        (! -e "$wrk/0-overlaptrim/$asm.mergeLog.bz2")) {

        backupFragStore("beforeTrimMerge");

        my $cmd;
        $cmd  = "$bin/merge-trimming ";
        $cmd .= "-immutable $im " if (defined($im));
        $cmd .= "-log $wrk/0-overlaptrim/$asm.mergeLog ";
        $cmd .= "-frg $wrk/$asm.frgStore ";
        $cmd .= "-ovl $wrk/0-overlaptrim/$asm.ovl.consolidated";
        $cmd .= "> $wrk/0-overlaptrim/$asm.ovl.consolidated.err 2>&1";

        if (runCommand($cmd)) {
            unlink "$wrk/0-overlaptrim/$asm.mergeLog";
            unlink "$wrk/0-overlaptrim/$asm.mergeLog.stats";
            die "Failed to merge trimming.\n";
        }
    }


    #  Be nice, and generate a report of our trimming done.
    #
    if (0) {

        if ((! -e "$wrk/0-overlaptrim/$asm.report") &&
            (! -e "$wrk/0-overlaptrim/$asm.report.bz2")) {


            #  Clean up stuff
            #   - add missing fragments to $wrk/0-overlaptrim/$asm.ovl.consolidated
            #
            open(F, "< $wrk/0-overlaptrim/$asm.ovl.consolidated");
            open(G, "> $wrk/0-overlaptrim/$asm.ovl.consolidated.full");
            my $inId = 0;
            my $otId = 0;
            while (<F>) {
                ($inId) = split '\s+', $_;
                $otId++;
                while ($otId < $inId) {
                    #print STDERR "$otId has no overlaps (but $inId does).\n";
                    print G "$otId  0 0 0 0 0  0 0 0 0 0  0\n";
                    $otId++;
                }
                print G $_;
                $otId = $inId;
            }
            close(F);

            $otId++;
            while ($otId <= $numFrags) {
                print G "$otId  0 0 0 0 0  0 0 0 0 0  0\n";
                $otId++;
            }
            close(G);


            open(A, "< $wrk/0-overlaptrim/$asm.qualityLog") or die "Failed to open $wrk/0-overlaptrim/$asm.qualityLog\n";
            open(B, "< $wrk/0-overlaptrim/$asm.mergeLog") or die "Failed to open $wrk/0-overlaptrim/$asm.mergeLog\n";
            open(C, "< $wrk/0-overlaptrim/$asm.ovl.consolidated.full") or die "Failed to open $wrk/0-overlaptrim/$asm.ovl.consolidated.full\n";
            open(F, "> $wrk/0-overlaptrim/$asm.report") or die "Failed to open $wrk/0-overlaptrim/$asm.report\n";

            while (!eof(A) || !eof(B) || !eof(C)) {
                my $a = <A>; chomp $a;
                my $b = <B>; chomp $b;
                my $c = <C>; chomp $c;

                my @av = split '\s+', $a;
                my @bv = split '\s+', $b;
                my @cv = split '\s+', $c;

                if (($av[0] != $bv[0]) || ($bv[0] != $cv[0]) || ($av[0] != $cv[0])) {
                    print STDERR "ERROR: ID MISMATCH!\n";
                    print STDERR "A: $a\nB: $b\nC: $c\n";
                    die;
                }

                printf(F "%6d : TI: %4d %4d Q1: %4d %4d Q2: %4d %4d TF: %4d %4d : %s\n",
                       $av[0],
                       $av[1], $av[2],  #  TI
                       $av[4], $av[5],  #  Q1
                       $bv[1], $bv[2],  #  Q2
                       $bv[3], $bv[4],  #  TF
                       $c);
            }

            close(C);
            close(B);
            close(A);
            close(F);

            unlink "$wrk/0-overlaptrim/$asm.ovl.consolidated.full";
        }
    }


    #  Add "-delete" to remove, instead of fix, chimera and spurs.
    #
    if ((! -e "$wrk/0-overlaptrim/$asm.chimera.report") &&
        (! -e "$wrk/0-overlaptrim/$asm.chimera.report.bz2")) {

        backupFragStore("beforeChimera");

        my $cmd;
        $cmd  = "$bin/chimera ";
        $cmd .= " -frg $wrk/$asm.frgStore ";
        $cmd .= " -immutable $im " if (defined($im));
        $cmd .= " -summary $wrk/0-overlaptrim/$asm.chimera.summary ";
        $cmd .= " -report  $wrk/0-overlaptrim/$asm.chimera.report ";
        $cmd .= " < $wrk/0-overlaptrim/$asm.ovl.sorted ";
        $cmd .= " 2> $wrk/0-overlaptrim/$asm.chimera.err ";
        if (runCommand($cmd)) {
            rename "$wrk/0-overlaptrim/$asm.chimera.report", "$wrk/0-overlaptrim/$asm.chimera.report.FAILED";
            die "Failed.\n";
        }
    }


    #  Finally, fix up gatekeeper, delete any mate links for fragments that we've deleted.
    #
    if (! -e "$wrk/0-overlaptrim/$asm.deletelinks.out") {
        my $cmd;
        $cmd  = "$bin/deleteLinks ";
        $cmd .= " -f $wrk/$asm.frgStore ";
        $cmd .= " -g $wrk/$asm.gkpStore ";
        $cmd .= " > $wrk/0-overlaptrim/$asm.deletelinks.out 2>&1";
        if (runCommand($cmd)) {
            rename "$wrk/0-overlaptrim/$asm.deletelinks.out", "$wrk/0-overlaptrim/$asm.deletelinks.out.FAILED";
            die "Failed.\n";
        }
    }

    touch("$wrk/0-overlaptrim/overlaptrim.success");

  alldone:
    stopAfter("overlapBasedTrimming");
    stopAfter("OBT");
}

1;
