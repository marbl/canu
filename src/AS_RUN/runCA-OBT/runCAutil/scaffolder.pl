use strict;

#  Don't do interleaved merging unless we are throwing stones.

sub CGW ($$$$$$) {
    my $thisDir     = shift @_;
    my $lastDir     = shift @_;
    my $cgiFile     = shift @_;
    my $stoneLevel  = shift @_;
    my $logickp     = shift @_;
    my $finalRun    = shift @_;

    return($thisDir) if (-e "$wrk/$thisDir/cgw.success");

    my $lastckp = 0;
    my $ckp     = "";

    $lastckp = findLastCheckpoint($lastDir)  if (defined($lastDir));
    $ckp     = "-y -R $lastckp -N $logickp"  if (defined($lastckp) && defined($logickp));

    #  If there is a timing file here, assume we are restarting.  Not
    #  all restarts are possible, but we try hard to make it so.
    #
    if (-e "$wrk/$thisDir/$asm.timing") {
        $ckp = "";

        open(F, "< $wrk/$thisDir/$asm.timing");
        while (<F>) {
            if (m/Done\swith/) {
                print $_;
            }
            if (m/Done\swith\scheckpoint\s(\d+)\s\(logical\s(\d+)\)/) {
                $ckp = "-y -R $1 -N $2";
            }
        }
        close(F);

        die "ERROR:  Found a timing file, but didn't find the checkpoint information!\n" if (!defined($ckp));
        print STDERR "Found a timing file, restarting: $ckp\n";
    }

    system("mkdir $wrk/$thisDir")               if (! -d "$wrk/$thisDir");
    system("mkdir $wrk/$asm.SeqStore")          if (! -d "$wrk/$asm.SeqStore");

    $cgiFile = "../5-consensus/$asm.cgi" if (!defined($cgiFile));

    system("ln -s $cgiFile          $wrk/$thisDir/$asm.cgi")          if (! -e "$wrk/$thisDir/$asm.cgi");
    system("ln -s ../$asm.SeqStore  $wrk/$thisDir/$asm.SeqStore")     if (! -e "$wrk/$thisDir/$asm.SeqStore");

    system("ln -s ../$lastDir/$asm.ckp.$lastckp $wrk/$thisDir/$asm.ckp.$lastckp") if (defined($lastDir));

    my $sampleSize = getGlobal("cgwDistanceSampleSize");
    
    my $cmd;
    $cmd  = "$bin/cgw $ckp -c -j 1 -k 5 -r 5 -s $stoneLevel ";
    $cmd .= " -S 0 " if (($finalRun == 0)   || (getGlobal("doResolveSurrogates") == 0));
    $cmd .= " -G "   if (($finalRun == 0)   && (getGlobal("cgwOutputIntermediate") == 0));
    $cmd .= " -M "   if (($stoneLevel == 0) && (getGlobal("delayInterleavedMerging") == 1));
    $cmd .= " -z "   if (getGlobal("cgwDemoteRBP") == 1);
    $cmd .= " -m $sampleSize";
    $cmd .= " -g $wrk/$asm.gkpStore ";
    $cmd .= " -o $wrk/$thisDir/$asm ";
    $cmd .= " $wrk/$thisDir/$asm.cgi ";
    $cmd .= " > $wrk/$thisDir/cgw.out 2>&1";
    if (runCommand("$wrk/$thisDir", $cmd)) {
        print STDERR "Failed.\n";
        exit(1);
    }


    open(F, "ls -1 $wrk/$thisDir |");
    while (<F>) {
        chomp;

        if (m/\.log$/) {
            system("mkdir $wrk/$thisDir/log")        if (! -d "$wrk/$thisDir/log");
            rename "$wrk/$thisDir/$_", "$wrk/$thisDir/log/$_";
        }

        if (m/\.analysis$/) {
            system("mkdir $wrk/$thisDir/analysis")   if (! -d "$wrk/$thisDir/analysis");
            rename "$wrk/$thisDir/$_", "$wrk/$thisDir/analysis/$_";
        }
    }
    close(F);


    if (getGlobal("cgwPurgeCheckpoints") != 0) {
        my $f = findFirstCheckpoint($thisDir);
        my $l = findLastCheckpoint($thisDir);

        while ($f < $l) {
            #print STDERR "Purging $wrk/$thisDir/$asm.ckp.$f\n";
            unlink "$wrk/$thisDir/$asm.ckp.$f";
            $f++;
        }
    }

    touch("$wrk/$thisDir/cgw.success");

    return $thisDir;
}


sub eCR ($$$) {
    my $thisDir = shift @_;
    my $lastDir = shift @_;
    my $iter    = shift @_;

    return $thisDir if (-e "$wrk/$thisDir/extendClearRanges.success");

    my $lastckp = findLastCheckpoint($lastDir);

    system("mkdir $wrk/$thisDir") if (! -d "$wrk/$thisDir");

    system("ln -s ../$lastDir/$asm.ckp.$lastckp $wrk/$thisDir/$asm.ckp.$lastckp")  if (! -e "$wrk/$thisDir/$asm.ckp.$lastckp");
    system("ln -s ../$asm.SeqStore              $wrk/$thisDir/$asm.SeqStore")      if (! -e "$wrk/$thisDir/$asm.SeqStore");

    #  Run eCR in smaller batches, hopefully making restarting from a failure both
    #  faster and easier.

    my $curScaffold  = 0;
    my $endScaffold  = 0;
    my $numScaffolds = findNumScaffoldsInCheckpoint($thisDir, $lastckp);
    my $stepSize     = getGlobal("extendClearRangesStepSize");

    my $substrlen = length("$numScaffolds");

    while ($curScaffold < $numScaffolds) {
        $endScaffold = $curScaffold + $stepSize;
        $endScaffold = $numScaffolds if ($endScaffold > $numScaffolds);

        $curScaffold = substr("000000000$curScaffold", -$substrlen);

        if (! -e "$wrk/$thisDir/extendClearRanges-scaffold.$curScaffold.success") {
            backupFragStore("before-$thisDir-scaffold.$curScaffold");

            $lastckp = findLastCheckpoint($thisDir);

            my $cmd;
            $cmd  = "$bin/extendClearRanges ";
            $cmd .= " -g $wrk/$asm.gkpStore ";
            $cmd .= " -n $lastckp ";
            $cmd .= " -c $asm ";
            $cmd .= " -b $curScaffold -e $endScaffold ";
            $cmd .= " -i $iter ";
            $cmd .= " > $wrk/$thisDir/extendClearRanges-scaffold.$curScaffold.err 2>&1";

            open(F, "> $wrk/$thisDir/extendClearRanges-scaffold.$curScaffold.sh");
            print F "#!/bin/sh\n\n";
            print F "$cmd\n";
            close(F);

            if (runCommand("$wrk/$thisDir", $cmd)) {
                print STDERR "Failed.\n";
                print STDERR "fragStore restored:    frg -> frg.during.$thisDir-scaffold.$curScaffold.FAILED\n";
                print STDERR "                       frg.before-$thisDir-scaffold.$curScaffold -> frg\n";
                rename "$wrk/$asm.gkpStore/frg", "$wrk/$asm.gkpStore/frg.during.$thisDir-scaffold.$curScaffold.FAILED";
                rename "$wrk/$asm.gkpStore/frg.before-$thisDir-scaffold.$curScaffold", "$wrk/$asm.gkpStore/frg";
                exit(1);
            }
            touch("$wrk/$thisDir/extendClearRanges-scaffold.$curScaffold.success");
        }

        $curScaffold = $endScaffold;
    }

    touch("$wrk/$thisDir/extendClearRanges.success");

    return $thisDir;
}


sub updateDistanceRecords ($) {
    my $thisDir = shift @_;

    return if (getGlobal("doUpdateDistanceRecords") == 0);
    return if (-e "$wrk/$thisDir/cgw.distupdate.success");

    #  Older versions needed to actually compute the updated
    #  distances.  Now, cgw outputs it!  Yay!

    my $cmd;
    $cmd  = "$bin/gatekeeper ";
    $cmd .= " -a -o $wrk/$asm.gkpStore ";
    $cmd .= " $wrk/$thisDir/stat/scaffold_final.distupdate.dst ";
    $cmd .= " $wrk/$thisDir/stat/contig_final.distupdate.dst ";
    $cmd .= " > $wrk/$thisDir/cgw.distupdate.err 2>&1";
    if (runCommand("$wrk/$thisDir", $cmd)) {
        die "Gatekeeper Failed.\n";
    }

    touch("$wrk/$thisDir/cgw.distupdate.success");
}


sub scaffolder ($) {
    my $cgiFile    = shift @_;
    my $lastDir    = undef;
    my $thisDir    = 0;
    my $stoneLevel = getGlobal("stoneLevel");

    goto alldone if (-e "$wrk/7-CGW/cgw.success");

    #  Do an initial CGW to update distances, then update the
    #  gatekeeper.  This initial run shouldn't be used for later
    #  CGW'ing.  We need to check explicitly for
    #  doUpdateDistanceRecords, otherwise that CGW() is run, _then_ we
    #  check if we should update.
    #
    if (getGlobal("doUpdateDistanceRecords")) {
        updateDistanceRecords(CGW("6-clonesize", undef, $cgiFile, $stoneLevel, undef, 0));
    }


    #  If we're not doing eCR, we just do a single scaffolder run, and
    #  get the heck outta here!  OK, we'll do resolveSurrogates(), maybe.
    #
    if (getGlobal("doExtendClearRanges") == 0) {
        $lastDir = CGW("7-$thisDir-CGW", $lastDir, $cgiFile, $stoneLevel, undef, 1);
        $thisDir++;
    } else {

        #  Do the initial CGW, making sure to not throw stones.
        #
        $lastDir = CGW("7-$thisDir-CGW", $lastDir, $cgiFile, 0, undef, 0);
        $thisDir++;

        #  Followed by at least one eCR
        #
        $lastDir = eCR("7-$thisDir-ECR", $lastDir, 1);
        $thisDir++;

        #  Iterate eCR: do another scaffolder still without stones,
        #  then another eCR.  Again, and again, until we get dizzy and
        #  fall over.
        #
        my $iterationMax = getGlobal("doExtendClearRanges") + 1;
        for (my $iteration = 2; $iteration < $iterationMax; $iteration++) {
            $lastDir = CGW("7-$thisDir-CGW", $lastDir, $cgiFile, 0, 3, 0);
            $thisDir++;

            $lastDir = eCR("7-$thisDir-ECR", $lastDir, $iteration);
            $thisDir++;
        }

        #  Then another scaffolder, chucking stones into the big holes,
        #  filling in surrogates, and writing output.
        #
        $lastDir = CGW("7-$thisDir-CGW", $lastDir, $cgiFile, $stoneLevel, 3, 1);
        $thisDir++;
    }


    #  And, finally, hold on, we're All Done!  Point to the correct output directory.
    #
    system("ln -s $lastDir $wrk/7-CGW") if (! -d "$wrk/7-CGW");

  alldone:
    stopAfter("scaffolder");
}


1;
