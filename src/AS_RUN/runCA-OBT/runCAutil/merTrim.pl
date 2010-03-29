use strict;

sub merTrim {

    return if (getGlobal("doMerBasedTrimming") == 0);
    return if (getGlobal("doOverlapBasedTrimming") == 1);
    return if (getGlobal("ovlOverlapper") eq "umd");

    #  Skip mer based trimming if it is done, or if the ovlStore already exists.
    #
    goto alldone if (-e "$wrk/0-overlaptrim/mertrim.success");
    goto alldone if (-d "$wrk/$asm.ovlStore");

    system("mkdir $wrk/0-overlaptrim")         if (! -d "$wrk/0-overlaptrim");
    system("mkdir $wrk/0-overlaptrim-overlap") if (! -d "$wrk/0-overlaptrim-overlap");


    if (! -e "$wrk/0-overlaptrim/$asm.merTrimLog") {
        meryl();

        my $merSize     = getGlobal("obtMerSize");
        my $merComp     = 0;  # getGlobal("merCompression");

        my $bin = getBinDirectory();
        my $cmd;
        $cmd  = "$bin/merTrim \\\n";
        $cmd .= " -g  $wrk/$asm.gkpStore \\\n";
        $cmd .= " -m  $merSize \\\n";
        $cmd .= " -c  $merComp \\\n";
        $cmd .= " -mc $wrk/0-mercounts/$asm-C-ms$merSize-cm$merComp \\\n";
        $cmd .= " -l  $wrk/0-overlaptrim/$asm.merTrimLog \\\n";
        $cmd .= " >  $wrk/0-overlaptrim/$asm.merTrim.err 2>&1\n";

        stopBefore("initialTrim", $cmd);

        if (runCommand("$wrk/0-overlaptrim", $cmd)) {
            caFailure("mer trimming failed", "$wrk/0-overlaptrim/$asm.merTrim.err");
        }

        touch("$wrk/0-overlaptrim/$asm.merTrimLog");
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




    if (! -e "$wrk/0-overlaptrim/$asm.overlapMask.err") {
        my $bin = getBinDirectory();
        my $cmd;
        $cmd  = "$bin/overlapMask \\\n";
        $cmd .= " -gkp $wrk/$asm.gkpStore \\\n";
        $cmd .= " -ovs $wrk/0-overlaptrim/$asm.obtStore \\\n";
        $cmd .= " > $wrk/0-overlaptrim/$asm.overlapMask.err 2>&1";

        if (runCommand("$wrk/0-overlaptrim", $cmd)) {
            caFailure("mer trim masking failed", "$wrk/0-overlaptrim/$asm.overlapMask.err");
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



    touch("$wrk/0-overlaptrim/mertrim.success");

  alldone:
    stopAfter("overlapBasedTrimming");
    stopAfter("OBT");
}

1;
