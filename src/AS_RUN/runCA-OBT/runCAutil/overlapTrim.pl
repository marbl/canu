use strict;

########################################
#
#  Do overlap based trimming
#
#  Do a leniant quality filter.  Run overlapper with the Granger
#  option (-G).  We used to fiddle with the sequences to convert
#  any N into a random base with low quality.

sub overlapTrim {

    system("mkdir $wrk/0-overlaptrim")         if (! -d "$wrk/0-overlaptrim");
    system("mkdir $wrk/0-overlaptrim-overlap") if (! -d "$wrk/0-overlaptrim-overlap");

    if ((! -e "$wrk/0-overlaptrim/$asm.trim.qualityLog") &&
        (! -e "$wrk/0-overlaptrim/$asm.trim.qualityLog.bz2")) {
        #die "Failed test - quality\n";
        print STDERR "Starting -- overlap trimming - quality trimming\n";

        if (($doBackups) && (! -e "$wrk/$asm.frgStore/db.frg.beforeQualityTrim")) {
            print STDERR "Backing up the frgStore.\n";
            if (runCommand("cp -p $wrk/$asm.frgStore/db.frg $wrk/$asm.frgStore/db.frg.beforeQualityTrim")) {
                unlink "$wrk/$asm.frgStore/db.frg.beforeQualityTrim";
                die "Failed to backup frgStore.\n";
            }
        }

        if (runCommand("$bin/qualityTrim -update -log $wrk/0-overlaptrim/$asm.trim.qualityLog -q 12 -frg $wrk/$asm.frgStore")) {
            rename "$wrk/0-overlaptrim/$asm.trim.quailtyLog", "$wrk/0-overlaptrim/$asm.trim.qualityLog.failed";
            die "Failed.\n";
        }
    }

    #  Do the _optional_ vector intersection

    if ((-e $vectorIntersect) &&
        (! -e "$wrk/0-overlaptrim/$asm.trim.vectorIntersectionLog") &&
        (! -e "$wrk/0-overlaptrim/$asm.trim.vectorIntersectionLog.bz2")) {
        #die "Failed test - intersect\n";
        print STDERR "Starting -- overlap trimming - vector intersection\n";

        if (($doBackups) && (! -e "$wrk/$asm.frgStore/db.frg.beforeOverlapIntersection")) {
            print STDERR "Backing up the frgStore.\n";
            if (runCommand("cp -p $wrk/$asm.frgStore/db.frg $wrk/$asm.frgStore/db.frg.beforeOverlapIntersection")) {
                unlink "$wrk/$asm.frgStore/db.frg.beforeOverlapIntersection";
                die "Failed to backup frgStore.\n";
            }
        }

        if (runCommand("$bin/intersectTrim -update -intersect $vectorIntersect -log $wrk/0-overlaptrim/$asm.trim.vectorIntersectionLog -frg $wrk/$asm.frgStore")) {
            rename "$wrk/0-overlaptrim/$asm.trim.vectorIntersectionLog", "$wrk/0-overlaptrim/$asm.trim.vectorIntersectionLog.failed";
            die "Failed.\n";
        }
    }

    #  NEW!  Run meryl on the now quality-trimmed frags.  This
    #  hopefully will get around BCM's habit of having N's in the low
    #  quality region, and is generally just a good idea.
    #
    meryl();



    #  Filter the standard set of nmers, throw out things below 100.
    #  If you change 100, you should also change meryl.pl.

    if (! -e "$wrk/0-overlaptrim-overlap/$asm.nmers.fasta") {
        print STDERR "Starting -- overlap trimming - meryl\n";

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

    if ((! -e "$wrk/0-overlaptrim/$asm.trim.ovl.sorted") &&
        (! -e "$wrk/0-overlaptrim/$asm.trim.ovl.sorted.bz2")) {
        #die "Failed test - sort\n";
        print STDERR "Starting -- overlap trimming - sorting\n";

        if (runCommand("find $wrk/0-overlaptrim-overlap -follow -name \\*ovb -print > $wrk/0-overlaptrim/all-overlaps-trim.ovllist")) {
            print STDERR "Failed to generate a list of all the overlap files.\n";
            exit(1);
        }

        if (runCommand("$bin/sort-overlaps -memory 16000 -maxiid $numFrags -L $wrk/0-overlaptrim/all-overlaps-trim.ovllist > $wrk/0-overlaptrim/$asm.trim.ovl.sorted")) {
            unlink "$wrk/0-overlaptrim/$asm.trim.ovl.sorted";
            die "Failed to sort.\n";
        }
    }

    #  Consolidate the overlaps, listing all overlaps for a single
    #  fragment on a single line.  These are still iid's.

    if ((! -e "$wrk/0-overlaptrim/$asm.trim.ovl.consolidated") &&
        (! -e "$wrk/0-overlaptrim/$asm.trim.ovl.consolidated.bz2")) {
        #die "Failed test - consolidate\n";
        print STDERR "Starting -- overlap trimming - consolidation\n";

        if (runCommand("$bin/consolidate < $wrk/0-overlaptrim/$asm.trim.ovl.sorted > $wrk/0-overlaptrim/$asm.trim.ovl.consolidated.1")) {
          unlink "$wrk/0-overlaptrim/$asm.ovl.trim.consolidated.1";
          die "Failed to sort.\n";
        }

        print STDERR "Starting -- overlap trimming - consolidation - cleaning\n";

        #  Clean up stuff
        #   - add missing fragments to $wrk/0-overlaptrim/$asm.trim.ovl.consolidated
        #
        open(G, "> $wrk/0-overlaptrim/$asm.trim.ovl.consolidated");
        open(F, "< $wrk/0-overlaptrim/$asm.trim.ovl.consolidated.1");
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

        print STDERR "$otId $numFrags\n";

        $otId++;
        while ($otId <= $numFrags) {
            print G "$otId  0 0 0 0 0  0 0 0 0 0  0\n";
            $otId++;
        }

        close(G);

        unlink "$wrk/0-overlaptrim/$asm.trim.ovl.consolidated.1";
    }


    #  We need to have all the overlaps squashed already, in particular so
    #  that we can get the mode of the 5'mode.  We could do this all in
    #  core, but that would take lots of space.

    #  This is for restarting -- I always seem to remove the *.ofg, and
    #  forget to rename the original back.
    if ((! -e "$wrk/0-preoverlap/$asm.ofg") && (-e "$wrk/0-preoverlap/$asm.ofg.orig")) {
        rename "$wrk/0-preoverlap/$asm.ofg.orig", "$wrk/0-preoverlap/$asm.ofg";
    }

    if ((! -e "$wrk/0-preoverlap/$asm.ofg.orig") &&
        (! -e "$wrk/0-preoverlap/$asm.ofg.orig.bz2") &&
        (! -e "$wrk/0-overlaptrim/$asm.trim.mergeLog") &&
        (! -e "$wrk/0-overlaptrim/$asm.trim.mergeLog.bz2")) {
        #die "Failed test - ofg - merge\n";
        print STDERR "Starting -- overlap trimming - merging\n";

        if (($doBackups) && (! -e "$wrk/$asm.frgStore/db.frg.beforeTrimMerge")) {
            print STDERR "Backing up the frgStore.\n";
            if (runCommand("cp -p $wrk/$asm.frgStore/db.frg $wrk/$asm.frgStore/db.frg.beforeTrimMerge")) {
                unlink "$wrk/$asm.frgStore/db.frg.beforeTrimMerge";
                die "Failed to backup frgStore.\n";
            }
        }

        runCommand("$bin/merge-trimming -log $wrk/0-overlaptrim/$asm.trim.mergeLog -frg $wrk/$asm.frgStore -ovl $wrk/0-overlaptrim/$asm.trim.ovl.consolidated") and die;
        runCommand("$bin/make_OFG_from_FragStore $wrk/$asm.frgStore > $wrk/0-preoverlap/$asm.2.ofg") and die;
        rename "$wrk/0-preoverlap/$asm.ofg", "$wrk/0-preoverlap/$asm.ofg.orig";
        rename "$wrk/0-preoverlap/$asm.2.ofg", "$wrk/0-preoverlap/$asm.ofg";
    }


    #  Be nice, and generate a report for Granger
    #
    if ((! -e "$wrk/0-overlaptrim/$asm.trim.report") &&
        (! -e "$wrk/0-overlaptrim/$asm.trim.report.bz2")) {
        #die "Failed test - report\n";
        print STDERR "Starting -- overlap trimming - reporting\n";

        open(A, "< $wrk/0-overlaptrim/$asm.trim.qualityLog") or die "Failed to open $wrk/0-overlaptrim/$asm.trim.qualityLog\n";
        open(B, "< $wrk/0-overlaptrim/$asm.trim.mergeLog") or die "Failed to open $wrk/0-overlaptrim/$asm.trim.mergeLog\n";
        open(C, "< $wrk/0-overlaptrim/$asm.trim.ovl.consolidated") or die "Failed to open $wrk/0-overlaptrim/$asm.trim.ovl.consolidated\n";
        open(F, "> $wrk/0-overlaptrim/$asm.trim.report") or die "Failed to open $wrk/0-overlaptrim/$asm.trim.report\n";

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
    }


    if ((! -e "$wrk/0-overlaptrim/$asm.trim.chimera.report") &&
        (! -e "$wrk/0-overlaptrim/$asm.trim.chimera.report.bz2")) {
        #die "Failed test - chimera\n";
        print STDERR "Starting -- overlap trimming - chimera\n";

        #  Add "-delete" to remove, instead of fix, chimera and spurs.

        runCommand("$bin/chimera -frg $wrk/$asm.frgStore < $wrk/0-overlaptrim/$asm.trim.ovl.sorted > $wrk/0-overlaptrim/$asm.trim.chimera.report") and die;
    }
}



1;
