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

        #  OBT needs to backup the frag store because it doesn't have
        #  enough entries to be non-destructive.  In particular, the
        #  merge (might) and chimera (definitely) read/write the same
        #  entry (OBT).
        #
        backupFragStore("beforeTrimming");

        my $bin = getBinDirectory();
        my $cmd;
        $cmd  = "$bin/initialTrim ";
        $cmd .= " -log $wrk/0-overlaptrim/$asm.initialTrimLog ";
        $cmd .= " -frg $wrk/$asm.gkpStore ";
        $cmd .= " >  $wrk/0-overlaptrim/$asm.initialTrim.report ";
        $cmd .= " 2> $wrk/0-overlaptrim/$asm.initialTrim.err ";

        if (runCommand("$wrk/0-overlaptrim", $cmd)) {
            restoreFragStoreBackup("beforeTrimming");
            rename "$wrk/0-overlaptrim/$asm.initialTrimLog", "$wrk/0-overlaptrim/$asm.initialTrimLog.FAILED";
            caFailure("initial trimming failed.\n");
        }

        rmrf("0-overlaptrim/$asm.initialTrim.err");
    }

    #  Look for any exact prefix reads -- only does something if
    #  you've got 454 reads.
    #
    if ((! -e "$wrk/0-overlaptrim/$asm.prefixDeleteLog") &&
        (! -e "$wrk/0-overlaptrim/$asm.prefixDeleteLog.bz2")) {

        #  See comment on first backupFragStore() call.
        backupFragStore("beforePrefixDelete");

        my $bin = getBinDirectory();
        my $cmd;
        $cmd  = "$bin/prefixDelete ";
        $cmd .= " -log $wrk/0-overlaptrim/$asm.prefixDeleteLog ";
        $cmd .= " -frg $wrk/$asm.gkpStore ";
        $cmd .= " >  $wrk/0-overlaptrim/$asm.prefixDelete.report ";
        $cmd .= " 2> $wrk/0-overlaptrim/$asm.prefixDelete.err ";

        if (runCommand("$wrk/0-overlaptrim", $cmd)) {
            restoreFragStoreBackup("beforePrefixDelete");
            rename "$wrk/0-overlaptrim/$asm.initialTrimLog", "$wrk/0-overlaptrim/$asm.initialTrimLog.FAILED";
            caFailure("prefix delete failed.\n");
        }

        rmrf("0-overlaptrim/$asm.prefixDelete.err");
    }

    #  Compute overlaps, if we don't have them already

    if (! -e "$wrk/$asm.obtStore") {

        createOverlapJobs("trim");
        checkOverlap("trim");

        #  Sort the overlaps -- this also duplicates each overlap so that
        #  all overlaps for a fragment A are localized.

        if (runCommand("$wrk/0-overlaptrim",
                       "find $wrk/0-overlaptrim-overlap -follow -name \\*ovb.gz -print > $wrk/$asm.obtStore.list")) {
            caFailure("Failed to generate a list of all the overlap files.\n");
        }

        my $bin = getBinDirectory();
        my $cmd;
        $cmd  = "$bin/overlapStore ";
        $cmd .= " -c $wrk/$asm.obtStore.BUILDING ";
        $cmd .= " -g $wrk/$asm.gkpStore ";
        $cmd .= " -M " . getGlobal('ovlStoreMemory');
        $cmd .= " -L $wrk/$asm.obtStore.list";
        $cmd .= " > $wrk/$asm.obtStore.err 2>&1";

        if (runCommand("$wrk/0-overlaptrim", $cmd)) {
            caFailure("Failed to build the obt store.\n");
        }

        rename "$wrk/$asm.obtStore.BUILDING", "$wrk/$asm.obtStore";

        rmrf("$asm.obtStore.list");
        rmrf("$asm.obtStore.err");
    }

    #  Consolidate the overlaps, listing all overlaps for a single
    #  fragment on a single line.  These are still iid's.

    if ((! -e "$wrk/0-overlaptrim/$asm.ovl.consolidated") &&
        (! -e "$wrk/0-overlaptrim/$asm.ovl.consolidated.bz2")) {

        my $bin = getBinDirectory();
        my $cmd;
        $cmd  = "$bin/consolidate ";
        $cmd .= " -ovs $wrk/$asm.obtStore";
        $cmd .= " > $wrk/0-overlaptrim/$asm.ovl.consolidated";

        if (runCommand("$wrk/0-overlaptrim", $cmd)) {
          unlink "$wrk/0-overlaptrim/$asm.ovl.consolidated";
          caFailure("Failed to consolidate.\n");
        }
    }


    #  We need to have all the overlaps squashed already, in particular so
    #  that we can get the mode of the 5'mode.  We could do this all in
    #  core, but that would take lots of space.

    if ((! -e "$wrk/0-overlaptrim/$asm.mergeLog") &&
        (! -e "$wrk/0-overlaptrim/$asm.mergeLog.bz2")) {

        #  See comment on first backupFragStore() call.
        backupFragStore("beforeTrimMerge");

        my $bin = getBinDirectory();
        my $cmd;
        $cmd  = "$bin/merge-trimming ";
        $cmd .= "-log $wrk/0-overlaptrim/$asm.mergeLog ";
        $cmd .= "-frg $wrk/$asm.gkpStore ";
        $cmd .= "-ovl $wrk/0-overlaptrim/$asm.ovl.consolidated ";
        $cmd .= "> $wrk/0-overlaptrim/$asm.merge.err 2>&1";

        if (runCommand("$wrk/0-overlaptrim", $cmd)) {
            restoreFragStoreBackup("beforeTrimMerge");
            unlink "$wrk/0-overlaptrim/$asm.mergeLog";
            unlink "$wrk/0-overlaptrim/$asm.mergeLog.stats";
            caFailure("Failed to merge trimming.\n");
        }
    }

    if ((! -e "$wrk/0-overlaptrim/$asm.chimera.report") &&
        (! -e "$wrk/0-overlaptrim/$asm.chimera.report.bz2")) {

        #  See comment on first backupFragStore() call.
        backupFragStore("beforeChimera");

        my $bin = getBinDirectory();
        my $cmd;
        $cmd  = "$bin/chimera ";
        $cmd .= " -gkp $wrk/$asm.gkpStore ";
        $cmd .= " -ovs $wrk/$asm.obtStore ";
        $cmd .= " -summary $wrk/0-overlaptrim/$asm.chimera.summary ";
        $cmd .= " -report  $wrk/0-overlaptrim/$asm.chimera.report ";
        $cmd .= " > $wrk/0-overlaptrim/$asm.chimera.err 2>&1";
        if (runCommand("$wrk/0-overlaptrim", $cmd)) {
            restoreFragStoreBackup("beforeChimera");
            rename "$wrk/0-overlaptrim/$asm.chimera.report", "$wrk/0-overlaptrim/$asm.chimera.report.FAILED";
            caFailure("chimera cleaning failed.\n");
        }
    }

    removeFragStoreBackup("beforePrefixDelete");
    removeFragStoreBackup("beforeTrimMerge");
    removeFragStoreBackup("beforeChimera");

    #  Well, except we never get here if UMD.  Needs some tweaking.

    removeFragStoreBackup("beforeVectorTrim");
    removeFragStoreBackup("beforeUMDOverlapper");

    backupFragStore("afterTrimming");

    rmrf("$asm.obtStore");

    touch("$wrk/0-overlaptrim/overlaptrim.success");

  alldone:
    stopAfter("overlapBasedTrimming");
    stopAfter("OBT");
}

1;
