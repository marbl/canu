use strict;

#  Don't do interleaved merging unless we are throwing stones.

sub CGW ($$$$$) {
    my $thisDir    = shift @_;
    my $lastDir    = shift @_;
    my $cgiFile    = shift @_;
    my $stoneLevel = shift @_;
    my $logickp    = shift @_;

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

    my $cmd;
    $cmd  = "$bin/cgw $ckp -c -j 1 -k 5 -r 5 -s $stoneLevel -w 0 -T -P ";
    $cmd .= " -M " if (($stoneLevel == 0) && (getGlobal("delayInterleavedMerging") == 1));
    $cmd .= " -f $wrk/$asm.frgStore ";
    $cmd .= " -g $wrk/$asm.gkpStore ";
    $cmd .= " -o $wrk/$thisDir/$asm ";
    $cmd .= " $wrk/$thisDir/$asm.cgi ";
    $cmd .= " > $wrk/$thisDir/cgw.out 2>&1";
    if (runCommand("$wrk/$thisDir", $cmd)) {
        print STDERR "Failed.\n";
        exit(1);
    }

    system("mkdir $wrk/$thisDir/log")        if (! -d "$wrk/$thisDir/log");
    system("mkdir $wrk/$thisDir/analysis")   if (! -d "$wrk/$thisDir/analysis");
    system("mv $wrk/$thisDir/*.log      $wrk/$thisDir/log");
    system("mv $wrk/$thisDir/*.analysis $wrk/$thisDir/analysis");

    touch("$wrk/$thisDir/cgw.success");

    return $thisDir;
}


sub eCR ($$) {
    my $thisDir = shift @_;
    my $lastDir = shift @_;

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
    my $stepSize     = 5000;

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
            $cmd .= " -B ";
            $cmd .= " -b $curScaffold -e $endScaffold ";
            $cmd .= " -f $wrk/$asm.frgStore ";
            $cmd .= " -g $wrk/$asm.gkpStore ";
            $cmd .= " -c $asm ";
            $cmd .= " -n $lastckp ";
            $cmd .= " > $wrk/$thisDir/extendClearRanges-scaffold.$curScaffold.err 2>&1";

            open(F, "> $wrk/$thisDir/extendClearRanges-scaffold.$curScaffold.sh");
            print F "#!/bin/sh\n\n";
            print F "$cmd\n";
            close(F);

            if (runCommand("$wrk/$thisDir", $cmd)) {
                print STDERR "Failed.\n";
                print STDERR "fragStore restored:    db.frg -> db.frg.during.$thisDir-scaffold.$curScaffold.FAILED\n";
                print STDERR "                       db.frg.before-$thisDir-scaffold.$curScaffold -> db.frg\n";
                rename "$wrk/$asm.frgStore/db.frg", "$wrk/$asm.frgStore/db.frg.during.$thisDir-scaffold.$curScaffold.FAILED";
                rename "$wrk/$asm.frgStore/db.frg.before-$thisDir-scaffold.$curScaffold", "$wrk/$asm.frgStore/db.frg";
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
    return if (-e "$wrk/$thisDir/distupdate/distupdate.success");

    my $lastckp    = findLastCheckpoint($thisDir);
    my $distupdate = "$wrk/$thisDir/distupdate";

    system("mkdir $distupdate") if (! -d "$distupdate");

    system("ln -s ../$asm.ckp.$lastckp   $distupdate/$asm.ckp.$lastckp")  if (! -e "$distupdate/$asm.ckp.$lastckp");
    system("ln -s ../$asm.SeqStore       $distupdate/$asm.SeqStore")      if (! -e "$distupdate/$asm.SeqStore");

    #  Instead of suffering through extreme pain making nice fake
    #  UIDs, we just use the checkpoint number to offset from some
    #  arbitrary UID.
    #
    my $terminateFakeUID = getGlobal("fakeUIDs");
    my $uidServer        = getGlobal("uidServer");

    #  Some funny baloney to get around perl wanting to print large
    #  numbers in scientific notation.  Makes you just want to hate and
    #  love perl, don't it?
    #
    $terminateFakeUID = "987654312198765" . (4320 + $lastckp + $terminateFakeUID) if ($terminateFakeUID > 0);

    my $cmd;
    $cmd  = "$bin/dumpDistanceEstimates ";
    $cmd .= "    -u "                    if ($terminateFakeUID == 0);
    $cmd .= "    -s $terminateFakeUID "  if ($terminateFakeUID  > 0);
    $cmd .= "    $uidServer "            if (defined($uidServer));
    $cmd .= " -f $wrk/$asm.frgStore ";
    $cmd .= " -g $wrk/$asm.gkpStore ";
    $cmd .= " -p $asm ";
    $cmd .= " -n $lastckp ";
    $cmd .= " > $distupdate/update.dst ";
    $cmd .= " 2> $distupdate/dumpDistanceEstimates.err ";
    if (runCommand("$distupdate", $cmd)) {
        rename "$distupdate/update.dst", "$distupdate/update.dst.FAILED";
        die "dumpDistanceEstimates Failed.\n";
    }

    $cmd  = "$bin/gatekeeper ";
    $cmd .= " -X -Q -C -P -a $wrk/$asm.gkpStore $distupdate/update.dst ";
    $cmd .= " > $distupdate/gatekeeper.err 2>&1";
    if (runCommand("$distupdate", $cmd)) {
        die "Gatekeeper Failed.\n";
    }

    touch("$distupdate/distupdate.success");
}


sub resolveSurrogates ($$) {
    my $thisDir = shift @_;
    my $lastDir = shift @_;

    return $thisDir if (-e "$wrk/$thisDir/resolveSurrogates.success");

    my $lastckp = findLastCheckpoint($lastDir);

    system("mkdir $wrk/$thisDir")  if (! -d "$wrk/$thisDir");

    system("ln -s ../$lastDir/$asm.ckp.$lastckp $wrk/$thisDir/$asm.ckp.$lastckp")  if (! -e "$wrk/$thisDir/$asm.ckp.$lastckp");
    system("ln -s ../$asm.SeqStore              $wrk/$thisDir/$asm.SeqStore")      if (! -e "$wrk/$thisDir/$asm.SeqStore");

    my $cmd;
    $cmd  = "$bin/resolveSurrogates ";
    $cmd .= " -f $wrk/$asm.frgStore ";
    $cmd .= " -g $wrk/$asm.gkpStore ";
    $cmd .= " -c $asm ";
    $cmd .= " -n $lastckp ";
    $cmd .= " -S 0.666 ";
    $cmd .= " > $wrk/$thisDir/resolveSurrogates.err 2>&1";
    if (runCommand("$wrk/$thisDir", $cmd)) {
        die "Failed.\n";
    }

    touch("$wrk/$thisDir/resolveSurrogates.success");

    return $thisDir;
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
    if ((getGlobal("updateDistanceType") eq "pre") && (getGlobal("doUpdateDistanceRecords"))) {
        updateDistanceRecords(CGW("6-clonesize", undef, $cgiFile, $stoneLevel, undef));
    }


    #  If we're not doing eCR, we just do a single scaffolder run, and
    #  get the heck outta here!  OK, we'll do resolveSurrogates(), maybe.
    #
    if (getGlobal("doExtendClearRanges") == 0) {
        $lastDir = CGW("7-$thisDir-CGW", $lastDir, $cgiFile, $stoneLevel, undef);
        $thisDir++;
    } else {

        #  Do the initial CGW, making sure to not throw stones.
        #
        $lastDir = CGW("7-$thisDir-CGW", $lastDir, $cgiFile, 0, undef);
        $thisDir++;

        #  Followed by at least one eCR, and a distance update.
        #
        $lastDir = eCR("7-$thisDir-ECR", $lastDir);
        $thisDir++;

        if (getGlobal("updateDistanceType") ne "pre") {
            updateDistanceRecords($lastDir);
        }

        #  Iterate eCR: do another scaffolder still without stones,
        #  then another eCR.  Again, and again, until we get dizzy and
        #  fall over.
        #
        for (my $rounds = getGlobal("doExtendClearRanges") - 1; $rounds > 0; $rounds--) {
            $lastDir = CGW("7-$thisDir-CGW", $lastDir, $cgiFile, 0, 3);
            $thisDir++;

            $lastDir = eCR("7-$thisDir-ECR", $lastDir);
            $thisDir++;

            if (getGlobal("updateDistanceType") ne "pre") {
                updateDistanceRecords($lastDir);
            }
        }

        #  Then another scaffolder, chucking stones into the big holes.
        #
        $lastDir = CGW("7-$thisDir-CGW", $lastDir, $cgiFile, $stoneLevel, 3);
        $thisDir++;
    }

    #  Then resolve surrogates.  And yet another scaffolder, this time
    #  to just do output
    #
    if (getGlobal("doResolveSurrogates")) {
        $lastDir = resolveSurrogates("7-$thisDir-resolveSurrogates", $lastDir);
        $thisDir++;
        
        $lastDir = CGW("7-$thisDir-CGW", $lastDir, $cgiFile, $stoneLevel, 14);
        $thisDir++;
    }


    #  And, finally, hold on, we're All Done!  Point to the correct output directory.
    #
    system("ln -s $lastDir $wrk/7-CGW") if (! -d "$wrk/7-CGW");

  alldone:
    stopAfter("scaffolder");
}


1;
