use strict;

sub overlapTrim {

    return if (getGlobal("doOverlapTrimming") == 0);
    return if (getGlobal("ovlOverlapper") eq "umd");

    goto alldone if (-e "$wrk/0-overlaptrim/overlaptrim.success");

    system("mkdir $wrk/0-overlaptrim")         if (! -d "$wrk/0-overlaptrim");
    system("mkdir $wrk/0-overlaptrim-overlap") if (! -d "$wrk/0-overlaptrim-overlap");

    #  Do an initial overly-permissive quality trimming, intersected
    #  with any known vector trimming.
    #
    if ((! -e "$wrk/0-overlaptrim/$asm.initialTrimLog") &&
        (! -e "$wrk/0-overlaptrim/$asm.initialTrimLog.bz2")) {
        my $bin = getBinDirectory();
        my $cmd;
        $cmd  = "$bin/initialTrim ";
        $cmd .= " -log $wrk/0-overlaptrim/$asm.initialTrimLog ";
        $cmd .= " -frg $wrk/$asm.gkpStore ";
        $cmd .= " >  $wrk/0-overlaptrim/$asm.initialTrim.report ";
        $cmd .= " 2> $wrk/0-overlaptrim/$asm.initialTrim.err ";

        if (runCommand("$wrk/0-overlaptrim", $cmd)) {
            rename "$wrk/0-overlaptrim/$asm.initialTrimLog", "$wrk/0-overlaptrim/$asm.initialTrimLog.FAILED";
            caFailure("initial trimming failed", "$wrk/0-overlaptrim/$asm.initialTrim.err");
        }

        unlink "0-overlaptrim/$asm.initialTrim.err";
    }

    #  Compute overlaps, if we don't have them already

    if (! -e "$wrk/0-overlaptrim/$asm.obtStore") {

        createOverlapJobs("trim");
        checkOverlap("trim");

        #  Sort the overlaps -- this also duplicates each overlap so that
        #  all overlaps for a fragment A are localized.

        if (runCommand("$wrk/0-overlaptrim",
                       "find $wrk/0-overlaptrim-overlap -follow -name \\*ovb.gz -print > $wrk/0-overlaptrim/$asm.obtStore.list")) {
            caFailure("failed to generate a list of all the overlap files", undef);
        }

        my $bin = getBinDirectory();
        my $cmd;
        $cmd  = "$bin/overlapStore ";
        $cmd .= " -O ";
        $cmd .= " -c $wrk/0-overlaptrim/$asm.obtStore.BUILDING ";
        $cmd .= " -g $wrk/$asm.gkpStore ";
        $cmd .= " -M " . getGlobal('ovlStoreMemory');
        $cmd .= " -L $wrk/0-overlaptrim/$asm.obtStore.list";
        $cmd .= " > $wrk/0-overlaptrim/$asm.obtStore.err 2>&1";

        if (runCommand("$wrk/0-overlaptrim", $cmd)) {
            caFailure("failed to build the obt store", "$wrk/0-overlaptrim/$asm.obtStore.err");
        }

        rename "$wrk/0-overlaptrim/$asm.obtStore.BUILDING", "$wrk/0-overlaptrim/$asm.obtStore";

        rmrf("$asm.obtStore.list");
        rmrf("$asm.obtStore.err");
    }

    if (getGlobal("doDeDuplication") != 0) {
    if (! -e "$wrk/0-overlaptrim/$asm.deduplicate.summary") {
        my $bin = getBinDirectory();
        my $cmd;

        if (! -e "$wrk/0-overlaptrim/$asm.dupStore") {
            if (runCommand("$wrk/0-overlaptrim",
                           "find $wrk/0-overlaptrim-overlap -follow -name \\*ovb.gz -print > $wrk/0-overlaptrim/$asm.dupStore.list")) {
                caFailure("failed to generate a list of all the overlap files", undef);
            }

            $cmd  = "$bin/overlapStore ";
            $cmd .= " -O -O ";
            $cmd .= " -c $wrk/0-overlaptrim/$asm.dupStore.BUILDING ";
            $cmd .= " -g $wrk/$asm.gkpStore ";
            $cmd .= " -M " . getGlobal('ovlStoreMemory');
            $cmd .= " -L $wrk/0-overlaptrim/$asm.dupStore.list";
            $cmd .= " > $wrk/0-overlaptrim/$asm.dupStore.err 2>&1";

            if (runCommand("$wrk/0-overlaptrim", $cmd)) {
                caFailure("failed to build the dup store", "$wrk/0-overlaptrim/$asm.dupStore.err");
            }

            rename "$wrk/0-overlaptrim/$asm.dupStore.BUILDING", "$wrk/0-overlaptrim/$asm.dupStore";

            rmrf("$asm.dupStore.list");
            rmrf("$asm.dupStore.err");
        }

        $cmd  = "$bin/deduplicate ";
        $cmd .= "-gkp     $wrk/$asm.gkpStore ";
        $cmd .= "-ovs     $wrk/0-overlaptrim/$asm.obtStore ";
        $cmd .= "-ovs     $wrk/0-overlaptrim/$asm.dupStore ";
        $cmd .= "-report  $wrk/0-overlaptrim/$asm.deduplicate.report ";
        $cmd .= "-summary $wrk/0-overlaptrim/$asm.deduplicate.summary ";
        $cmd .= "> $wrk/0-overlaptrim/$asm.deduplicate.err 2>&1";

        if (runCommand("$wrk/0-overlaptrim", $cmd)) {
            unlink "$wrk/0-overlaptrim/$asm.deduplicate.summary";
            caFailure("failed to deduplicate the reads", "$wrk/0-overlaptrim/$asm.deduplicate.err");
        }
    }
    }

    #  Consolidate the overlaps, listing all overlaps for a single
    #  fragment on a single line.  These are still iid's.

    if ((! -e "$wrk/0-overlaptrim/$asm.ovl.consolidated") &&
        (! -e "$wrk/0-overlaptrim/$asm.ovl.consolidated.bz2")) {

        my $bin = getBinDirectory();
        my $cmd;
        $cmd  = "$bin/consolidate ";
        $cmd .= " -ovs $wrk/0-overlaptrim/$asm.obtStore ";
        $cmd .= " > $wrk/0-overlaptrim/$asm.ovl.consolidated ";
        $cmd .= "2> $wrk/0-overlaptrim/$asm.ovl.consolidated.err";

        if (runCommand("$wrk/0-overlaptrim", $cmd)) {
          unlink "$wrk/0-overlaptrim/$asm.ovl.consolidated";
          caFailure("failed to consolidate overlaps", "$wrk/0-overlaptrim/$asm.ovl.consolidated.err");
        }
        unlink "$wrk/0-overlaptrim/$asm.ovl.consolidated.err";
    }


    #  We need to have all the overlaps squashed already, in particular so
    #  that we can get the mode of the 5'mode.  We could do this all in
    #  core, but that would take lots of space.

    if ((! -e "$wrk/0-overlaptrim/$asm.mergeLog") &&
        (! -e "$wrk/0-overlaptrim/$asm.mergeLog.bz2")) {
        my $bin = getBinDirectory();
        my $cmd;
        $cmd  = "$bin/merge-trimming ";
        $cmd .= "-log $wrk/0-overlaptrim/$asm.mergeLog ";
        $cmd .= "-frg $wrk/$asm.gkpStore ";
        $cmd .= "-ovl $wrk/0-overlaptrim/$asm.ovl.consolidated ";
        $cmd .= "> $wrk/0-overlaptrim/$asm.merge.err 2>&1";

        if (runCommand("$wrk/0-overlaptrim", $cmd)) {
            unlink "$wrk/0-overlaptrim/$asm.mergeLog";
            unlink "$wrk/0-overlaptrim/$asm.mergeLog.stats";
            caFailure("failed to merge trimming", "$wrk/0-overlaptrim/$asm.merge.err");
        }
    }

    if (getGlobal("doChimeraDetection") != 0) {
        if ((! -e "$wrk/0-overlaptrim/$asm.chimera.report") &&
            (! -e "$wrk/0-overlaptrim/$asm.chimera.report.bz2")) {
            my $bin = getBinDirectory();
            my $cmd;
            $cmd  = "$bin/chimera ";
            $cmd .= " -gkp $wrk/$asm.gkpStore ";
            $cmd .= " -ovs $wrk/0-overlaptrim/$asm.obtStore ";
            $cmd .= " -summary $wrk/0-overlaptrim/$asm.chimera.summary ";
            $cmd .= " -report  $wrk/0-overlaptrim/$asm.chimera.report ";
            $cmd .= " > $wrk/0-overlaptrim/$asm.chimera.err 2>&1";
            if (runCommand("$wrk/0-overlaptrim", $cmd)) {
                rename "$wrk/0-overlaptrim/$asm.chimera.report", "$wrk/0-overlaptrim/$asm.chimera.report.FAILED";
                caFailure("chimera cleaning failed", "$wrk/0-overlaptrim/$asm.chimera.err");
            }
        }
    }

    rmrf("$asm.obtStore");

    touch("$wrk/0-overlaptrim/overlaptrim.success");

  alldone:
    stopAfter("overlapBasedTrimming");
    stopAfter("OBT");
}

1;
