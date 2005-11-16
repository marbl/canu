use strict;

########################################
#
#  Do overlap based trimming
#
#  Do a leniant quality filter.  Run overlapper with the Granger
#  option (-G).  We used to fiddle with the sequences to convert
#  any N into a random base with low quality.

sub backupFragStore ($) {
    my $backupName = shift @_;
    my $doBackups  = getGlobal("doBackupFragStore", 1);

    return if ($doBackups == 0);

    if (-e "$wrk/$asm.frgStore/db.frg.$backupName") {
        print STDERR "Found a backup for $backupName!  Restoring!\n";
        unlink "$wrk/$asm.frgStore/db.frg";
        if (runCommand("cp -p $wrk/$asm.frgStore/db.frg.$backupName $wrk/$asm.frgStore/db.frg")) {
            unlink "$wrk/$asm.frgStore/db.frg";
            die "Failed to restore frgStore from backup.\n";
        }
    }
    if (! -e "$wrk/$asm.frgStore/db.frg.$backupName") {
        print STDERR "Backing up the frgStore to $backupName.\n";
        if (runCommand("cp -p $wrk/$asm.frgStore/db.frg $wrk/$asm.frgStore/db.frg.$backupName")) {
            unlink "$wrk/$asm.frgStore/db.frg.$backupName";
            die "Failed to backup frgStore.\n";
        }
    }
}


sub overlapTrim {

    return if (getGlobal("doOverlapTrimming", 1) == 0);

    system("mkdir $wrk/0-overlaptrim")         if (! -d "$wrk/0-overlaptrim");
    system("mkdir $wrk/0-overlaptrim-overlap") if (! -d "$wrk/0-overlaptrim-overlap");

    #  Do an initial overly-permissive quality trimming, intersected
    #  with any known vector trimming.
    #
    if ((! -e "$wrk/0-overlaptrim/$asm.initialTrimLog") &&
        (! -e "$wrk/0-overlaptrim/$asm.initialTrimLog.bz2")) {
        print STDERR "Starting -- overlap trimming - initial trimming\n";

        my $vi = getGlobal("vectorIntersect", undef);
        my $im = getGlobal("immutableFrags", undef);
        die "vectorIntersect '$vi' supplied, but not found!\n" if (defined($vi) && (! -e $vi));
        die "immutableFrags '$im' supplied, but not found!\n" if (defined($im) && (! -e $im));

        backupFragStore("beforeInitialTrim");

        my $cmd;
        $cmd  = "$bin/initialTrim -update -q 12 ";
        $cmd .= " -vector $vi " if (defined($vi));
        $cmd .= " -immutable $im " if (defined($im));
        $cmd .= " -log $wrk/0-overlaptrim/$asm.initialTrimLog ";
        $cmd .= " -frg $wrk/$asm.frgStore ";
        $cmd .= " 2>&1 > $wrk/0-overlaptrim/initialTrim.err";

        if (runCommand($cmd)) {
            rename "$wrk/0-overlaptrim/$asm.initialTrimLog", "$wrk/0-overlaptrim/$asm.initialTrimLog.failed";
            die "Failed.\n";
        }
    }


    #  Run meryl on the now quality-trimmed frags.  This hopefully
    #  will get around some sequencing centers' habit of having N's in
    #  the low quality region, and is generally just a good idea.
    #
    meryl();

    #  Filter the standard set of nmers, throw out things below 100.
    #  If you change 100, you should also change meryl.pl.

    if (! -e "$wrk/0-overlaptrim-overlap/$asm.nmers.fasta") {
        print STDERR "Starting -- overlap trimming - meryl (filter)\n";

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

    if ((! -e "$wrk/0-overlaptrim/$asm.ovl.sorted") &&
        (! -e "$wrk/0-overlaptrim/$asm.ovl.sorted.bz2")) {
        print STDERR "Starting -- overlap trimming - sorting\n";

        if (runCommand("find $wrk/0-overlaptrim-overlap -follow -name \\*ovb -print > $wrk/0-overlaptrim/all-overlaps-trim.ovllist")) {
            print STDERR "Failed to generate a list of all the overlap files.\n";
            exit(1);
        }

        my $cmd;
        $cmd  = "$bin/sort-overlaps";
        $cmd .= " -memory " . getGlobal('ovlSortMemory', 1000) . " ";
        $cmd .= " -maxiid $numFrags ";
        $cmd .= " -L $wrk/0-overlaptrim/all-overlaps-trim.ovllist";
        $cmd .= " > $wrk/0-overlaptrim/$asm.ovl.sorted";

        if (runCommand($cmd)) {
            unlink "$wrk/0-overlaptrim/$asm.ovl.sorted";
            die "Failed to sort.\n";
        }
    }

    #  Consolidate the overlaps, listing all overlaps for a single
    #  fragment on a single line.  These are still iid's.

    if ((! -e "$wrk/0-overlaptrim/$asm.ovl.consolidated") &&
        (! -e "$wrk/0-overlaptrim/$asm.ovl.consolidated.bz2")) {
        print STDERR "Starting -- overlap trimming - consolidation\n";

        if (runCommand("$bin/consolidate < $wrk/0-overlaptrim/$asm.ovl.sorted > $wrk/0-overlaptrim/$asm.ovl.consolidated.1")) {
          unlink "$wrk/0-overlaptrim/$asm.ovl.consolidated.1";
          die "Failed to sort.\n";
        }

        print STDERR "Starting -- overlap trimming - consolidation - cleaning\n";

        #  Clean up stuff
        #   - add missing fragments to $wrk/0-overlaptrim/$asm.ovl.consolidated
        #
        open(G, "> $wrk/0-overlaptrim/$asm.ovl.consolidated");
        open(F, "< $wrk/0-overlaptrim/$asm.ovl.consolidated.1");
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

        unlink "$wrk/0-overlaptrim/$asm.ovl.consolidated.1";
    }


    #  We need to have all the overlaps squashed already, in particular so
    #  that we can get the mode of the 5'mode.  We could do this all in
    #  core, but that would take lots of space.

    if ((! -e "$wrk/0-overlaptrim/$asm.mergeLog") &&
        (! -e "$wrk/0-overlaptrim/$asm.mergeLog.bz2")) {
        print STDERR "Starting -- overlap trimming - merging\n";

        backupFragStore("beforeTrimMerge");

        my $im = getGlobal("immutableFrags", undef);
        die "immutableFrags '$im' supplied, but not found!\n" if (defined($im) && (! -e $im));

        my $cmd;
        $cmd  = "$bin/merge-trimming ";
        $cmd .= "-immutable $im " if (defined($im));
        $cmd .= "-log $wrk/0-overlaptrim/$asm.mergeLog ";
        $cmd .= "-frg $wrk/$asm.frgStore ";
        $cmd .= "-ovl $wrk/0-overlaptrim/$asm.ovl.consolidated";

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
            print STDERR "Starting -- overlap trimming - reporting\n";

            open(A, "< $wrk/0-overlaptrim/$asm.qualityLog") or die "Failed to open $wrk/0-overlaptrim/$asm.qualityLog\n";
            open(B, "< $wrk/0-overlaptrim/$asm.mergeLog") or die "Failed to open $wrk/0-overlaptrim/$asm.mergeLog\n";
            open(C, "< $wrk/0-overlaptrim/$asm.ovl.consolidated") or die "Failed to open $wrk/0-overlaptrim/$asm.ovl.consolidated\n";
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
        }
    }


    #  Add "-delete" to remove, instead of fix, chimera and spurs.
    #
    if ((! -e "$wrk/0-overlaptrim/$asm.chimera.report") &&
        (! -e "$wrk/0-overlaptrim/$asm.chimera.report.bz2")) {
        print STDERR "Starting -- overlap trimming - chimera\n";

        backupFragStore("beforeChimera");

        my $cmd;
        $cmd  = "$bin/chimera ";
        $cmd .= " -frg $wrk/$asm.frgStore ";
        $cmd .= " -summary $wrk/0-overlaptrim/$asm.chimera.summary ";
        $cmd .= " -report  $wrk/0-overlaptrim/$asm.chimera.report ";
        $cmd .= " < $wrk/0-overlaptrim/$asm.ovl.sorted ";
        $cmd .= " 2> $wrk/0-overlaptrim/$asm.chimera.err ";
        if (runCommand($cmd)) {
            rename "$wrk/0-overlaptrim/$asm.chimera.report", "$wrk/0-overlaptrim/$asm.chimera.report.FAILED";
            die "Failed.\n";
        }
    }


    #  Finally, for any UID's in the "do not trim" file, restore their
    #  original trimming, and possibly undelete them.


    #  Unitigger needs an input .ofg file, which is pretty much just a dump of the fragstore.
    #
    if ((! -e "$wrk/0-preoverlap/$asm.ofg.orig") &&
        (! -e "$wrk/0-preoverlap/$asm.ofg.orig.bz2")) {
        if (runCommand("$bin/make_OFG_from_FragStore $wrk/$asm.frgStore > $wrk/0-preoverlap/$asm.2.ofg")) {
            unlink "$wrk/0-preoverlap/$asm.2.ofg";
            die "Failed.\n";
        }
        rename "$wrk/0-preoverlap/$asm.ofg", "$wrk/0-preoverlap/$asm.ofg.orig";
        rename "$wrk/0-preoverlap/$asm.2.ofg", "$wrk/0-preoverlap/$asm.ofg";
    }
}



1;
