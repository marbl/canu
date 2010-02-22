use strict;

sub overlapTrim {

    return if (getGlobal("doOverlapBasedTrimming") == 0);
    return if (getGlobal("doMerBasedTrimming") == 1);
    return if (getGlobal("ovlOverlapper") eq "umd");

    #  Skip overlap based trimming if it is done, or if the ovlStore already exists.
    #
    goto alldone if (-e "$wrk/0-overlaptrim/overlaptrim.success");
    goto alldone if (-d "$wrk/$asm.ovlStore");

    system("mkdir $wrk/0-overlaptrim")         if (! -d "$wrk/0-overlaptrim");
    system("mkdir $wrk/0-overlaptrim-overlap") if (! -d "$wrk/0-overlaptrim-overlap");


    #  Disable dedup, unless reads request it.  This avoids an expensive ovlStore build.
    #
    if (getGlobal("doDeDuplication") != 0) {
        my $bin = getBinDirectory();

        setGlobal("doDeDuplication", 0);

        if (system("$bin/gatekeeper -isfeatureset 0 doRemoveDuplicateReads $wrk/$asm.gkpStore") == 0) {
           setGlobal("doDeDuplication", 1);   
        }
    }

    #  Do an initial overly-permissive quality trimming, intersected
    #  with any known vector trimming.  This is skipped if mer based trimming is enabled.
    #
#
#  TESTING:  DO NOT DO THIS IF WE HAVE MER TRIMMED
#
if (getGlobal("doMerBasedTrimming") == 0) {
    if ((! -e "$wrk/0-overlaptrim/$asm.initialTrimLog") &&
        (! -e "$wrk/0-overlaptrim/$asm.initialTrimLog.bz2")) {
        my $bin = getBinDirectory();
        my $cmd;
        $cmd  = "$bin/initialTrim \\\n";
        $cmd .= " -log $wrk/0-overlaptrim/$asm.initialTrimLog \\\n";
        $cmd .= " -frg $wrk/$asm.gkpStore \\\n";
        $cmd .= " >  $wrk/0-overlaptrim/$asm.initialTrim.report \\\n";
        $cmd .= " 2> $wrk/0-overlaptrim/$asm.initialTrim.err ";

        stopBefore("initialTrim", $cmd);

        if (runCommand("$wrk/0-overlaptrim", $cmd)) {
            rename "$wrk/0-overlaptrim/$asm.initialTrimLog", "$wrk/0-overlaptrim/$asm.initialTrimLog.FAILED";
            caFailure("initial trimming failed", "$wrk/0-overlaptrim/$asm.initialTrim.err");
        }

        unlink "0-overlaptrim/$asm.initialTrim.err";
    }
}

    #  Compute overlaps, if we don't have them already
    #
    if (! -e "$wrk/0-overlaptrim/$asm.obtStore") {
        createOverlapJobs("trim");
        checkOverlap("trim");

        #  Sort the overlaps -- this also duplicates each overlap so that
        #  all overlaps for a fragment A are localized.

        if (runCommand("$wrk/0-overlaptrim",
                       "find $wrk/0-overlaptrim-overlap -follow \\( -name \\*ovb.gz -or -name \\*ovb \\) -print > $wrk/0-overlaptrim/$asm.obtStore.list")) {
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

        #  Delete overlaps unless we're told to save them, or we need to dedup.
        if ((getGlobal("saveOverlaps") == 0) && (getGlobal("doDeDuplication") == 0)) {
            open(F, "< $wrk/0-overlaptrim/$asm.obtStore.list");
            while (<F>) {
                chomp;
                unlink $_;
            }
            close(F);
        }

        rmrf("$wrk/0-overlaptrim/$asm.obtStore.list");
        rmrf("$wrk/0-overlaptrim/$asm.obtStore.err");
    }


    #  Deduplicate?
    #
#
#  TESTING:  DO NOT DO THIS IF WE HAVE MER TRIMMED
#
if (getGlobal("doMerBasedTrimming") == 0) {
    if ((getGlobal("doDeDuplication") != 0) &&
        (! -e "$wrk/0-overlaptrim/$asm.deduplicate.summary")) {
        my $bin = getBinDirectory();
        my $cmd;

        if (! -e "$wrk/0-overlaptrim/$asm.dupStore") {
            if (runCommand("$wrk/0-overlaptrim",
                           "find $wrk/0-overlaptrim-overlap -follow \\( -name \\*ovb.gz -or -name \\*ovb \\) -print > $wrk/0-overlaptrim/$asm.dupStore.list")) {
                caFailure("failed to generate a list of all the overlap files", undef);
            }

            $cmd  = "$bin/overlapStore \\\n";
            $cmd .= " -O -O \\\n";
            $cmd .= " -c $wrk/0-overlaptrim/$asm.dupStore.BUILDING \\\n";
            $cmd .= " -g $wrk/$asm.gkpStore \\\n";
            $cmd .= " -M \\\n" . getGlobal('ovlStoreMemory');
            $cmd .= " -L $wrk/0-overlaptrim/$asm.dupStore.list \\\n";
            $cmd .= " > $wrk/0-overlaptrim/$asm.dupStore.err 2>&1";

            if (runCommand("$wrk/0-overlaptrim", $cmd)) {
                caFailure("failed to build the dup store", "$wrk/0-overlaptrim/$asm.dupStore.err");
            }

            rename "$wrk/0-overlaptrim/$asm.dupStore.BUILDING", "$wrk/0-overlaptrim/$asm.dupStore";

            #  Delete overlaps unless we're told to save them
            if (getGlobal("saveOverlaps") == 0) {
                open(F, "< $wrk/0-overlaptrim/$asm.dupStore.list");
                while (<F>) {
                    chomp;
                    unlink $_;
                }
                close(F);
            }

            rmrf("$asm.dupStore.list");
            rmrf("$asm.dupStore.err");
        }

        $cmd  = "$bin/deduplicate \\\n";
        $cmd .= "-gkp     $wrk/$asm.gkpStore \\\n";
        $cmd .= "-ovs     $wrk/0-overlaptrim/$asm.obtStore \\\n";
        $cmd .= "-ovs     $wrk/0-overlaptrim/$asm.dupStore \\\n";
        $cmd .= "-report  $wrk/0-overlaptrim/$asm.deduplicate.report \\\n";
        $cmd .= "-summary $wrk/0-overlaptrim/$asm.deduplicate.summary \\\n";
        $cmd .= "> $wrk/0-overlaptrim/$asm.deduplicate.err 2>&1";

        stopBefore("deDuplication", $cmd);

        if (runCommand("$wrk/0-overlaptrim", $cmd)) {
            unlink "$wrk/0-overlaptrim/$asm.deduplicate.summary";
            caFailure("failed to deduplicate the reads", "$wrk/0-overlaptrim/$asm.deduplicate.err");
        }
    }

    #  Consolidate the overlaps, listing all overlaps for a single
    #  fragment on a single line.  These are still iid's.

    if ((! -e "$wrk/0-overlaptrim/$asm.ovl.consolidated") &&
        (! -e "$wrk/0-overlaptrim/$asm.ovl.consolidated.bz2")) {

        my $bin = getBinDirectory();
        my $cmd;
        $cmd  = "$bin/consolidate \\\n";
        $cmd .= " -ovs $wrk/0-overlaptrim/$asm.obtStore \\\n";
        $cmd .= " > $wrk/0-overlaptrim/$asm.ovl.consolidated \\\n";
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
        $cmd  = "$bin/merge-trimming \\\n";
        $cmd .= "-log $wrk/0-overlaptrim/$asm.mergeLog \\\n";
        $cmd .= "-frg $wrk/$asm.gkpStore \\\n";
        $cmd .= "-ovl $wrk/0-overlaptrim/$asm.ovl.consolidated \\\n";
        $cmd .= "> $wrk/0-overlaptrim/$asm.merge.err 2>&1";

        stopBefore("mergeTrimming", $cmd);

        if (runCommand("$wrk/0-overlaptrim", $cmd)) {
            unlink "$wrk/0-overlaptrim/$asm.mergeLog";
            unlink "$wrk/0-overlaptrim/$asm.mergeLog.stats";
            caFailure("failed to merge trimming", "$wrk/0-overlaptrim/$asm.merge.err");
        }
    }
}

    if (getGlobal("doChimeraDetection") ne 'off') {
        if ((! -e "$wrk/0-overlaptrim/$asm.chimera.report") &&
            (! -e "$wrk/0-overlaptrim/$asm.chimera.report.bz2")) {
            my $bin = getBinDirectory();
            my $cmd;
            $cmd  = "$bin/chimera \\\n";
            $cmd .= " -gkp $wrk/$asm.gkpStore \\\n";
            $cmd .= " -ovs $wrk/0-overlaptrim/$asm.obtStore \\\n";
            $cmd .= " -summary $wrk/0-overlaptrim/$asm.chimera.summary \\\n";
            $cmd .= " -report  $wrk/0-overlaptrim/$asm.chimera.report \\\n";
            $cmd .= " -mininniepair 0 -minoverhanging 0 \\\n" if (getGlobal("doChimeraDetection") eq "aggressive");
            $cmd .= " > $wrk/0-overlaptrim/$asm.chimera.err 2>&1";

            stopBefore("chimeraDetection", $cmd);

            if (runCommand("$wrk/0-overlaptrim", $cmd)) {
                rename "$wrk/0-overlaptrim/$asm.chimera.report", "$wrk/0-overlaptrim/$asm.chimera.report.FAILED";
                caFailure("chimera cleaning failed", "$wrk/0-overlaptrim/$asm.chimera.err");
            }
        }
    }

    #rmrf("$asm.obtStore");

    touch("$wrk/0-overlaptrim/overlaptrim.success");

  alldone:
    stopAfter("overlapBasedTrimming");
    stopAfter("OBT");
}

1;
