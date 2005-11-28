use strict;

#  should we require stoneLevel be passed into scaffolder?

sub CGW ($$$$) {
    my $thisDir    = shift @_;
    my $lastDir    = shift @_;
    my $stoneLevel = shift @_;
    my $logickp    = shift @_;

    return $thisDir if (-e "$wrk/$thisDir/cgw.success");

    my $lastckp = findLastCheckpoint($lastDir)  if (defined($lastDir));
    my $ckp     = "-y -R $lastckp -N $logickp"  if (defined($lastckp) && defined($logickp));

    #  If there is a timing file here, assume we are restarting.  Not
    #  all restarts are possible, but we try hard to make it so.
    #
    if (-e "$wrk/$thisDir/$asm.timing") {
        undef $ckp;

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

    system("ln -s $wrk/5-consensus/$asm.cgi $wrk/$thisDir/$asm.cgi")          if (! -e "$wrk/$thisDir/$asm.cgi");
    system("ln -s $wrk/$asm.SeqStore        $wrk/$thisDir/$asm.SeqStore")     if (! -e "$wrk/$thisDir/$asm.SeqStore");

    system("ln -s $wrk/$lastDir/$asm.ckp.$lastckp $wrk/$thisDir/$asm.ckp.$lastckp") if (defined($lastDir));

    my $cmd;
    $cmd  = "cd $wrk/$thisDir && ";
    $cmd .= "$bin/cgw $ckp -c -j 1 -k 5 -r 4 -s $stoneLevel -w 0 -T -P ";
    $cmd .= " -f $wrk/$asm.frgStore ";
    $cmd .= " -g $wrk/$asm.gkpStore ";
    $cmd .= " -o $wrk/$thisDir/$asm ";
    $cmd .= " $wrk/$thisDir/$asm.cgi ";
    $cmd .= " > $wrk/$thisDir/cgw.out 2>&1";
    if (runCommand($cmd)) {
        print STDERR "Failed.\n";
        exit(1);
    }

    touch("$wrk/$thisDir/cgw.success");

    return $thisDir;
}


sub eCR ($$) {
    my $thisDir = shift @_;
    my $lastDir = shift @_;

    return $thisDir if (-e "$wrk/$thisDir/extendClearRanges.success");

    my $lastckp = findLastCheckpoint($lastDir);

    system("mkdir $wrk/$thisDir") if (! -d "$wrk/$thisDir");

    system("ln -s $wrk/$lastDir/$asm.ckp.$lastckp $wrk/$thisDir/$asm.ckp.$lastckp")  if (! -e "$wrk/$thisDir/$asm.ckp.$lastckp");
    system("ln -s $wrk/$asm.SeqStore              $wrk/$thisDir/$asm.SeqStore")      if (! -e "$wrk/$thisDir/$asm.SeqStore");

    backupFragStore("before-$thisDir");

    my $cmd;
    $cmd  = "cd $wrk/$thisDir && ";
    $cmd .= "$bin/extendClearRanges ";
    $cmd .= " -B ";
    $cmd .= " -f $wrk/$asm.frgStore ";
    $cmd .= " -g $wrk/$asm.gkpStore ";
    $cmd .= " -c $asm ";
    $cmd .= " -n $lastckp ";
    $cmd .= " -s -1 ";
    $cmd .= " > $wrk/$thisDir/extendClearRanges.err 2>&1";
    if (runCommand($cmd)) {
        print STDERR "Failed.\n";
        exit(1);
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

    system("ln -s $wrk/$thisDir/$asm.ckp.$lastckp $distupdate/$asm.ckp.$lastckp")  if (! -e "$distupdate/$asm.ckp.$lastckp");
    system("ln -s $wrk/$asm.SeqStore              $distupdate/$asm.SeqStore")      if (! -e "$distupdate/$asm.SeqStore");

    #  Instead of suffering through extreme pain making nice fake
    #  UIDs, we just use the checkpoint number to offset from some
    #  arbitrary UID.
    #
    my $terminateFakeUID = getGlobal("fakeUIDs");
    my $uidServer        = getGlobal("uidServer");

    #  Some funny baloney to get around perl wanting to print large
    #  numbers in scientific notation.  Make you just want to hate and
    #  love perl, don't it?
    #
    $terminateFakeUID = "987654312198765" . (4321 + $lastckp) if ($terminateFakeUID > 0);

    my $cmd;
    $cmd  = "cd $distupdate && $bin/dumpDistanceEstimates ";
    $cmd .= "    -u "                    if ($terminateFakeUID == 0);
    $cmd .= "    -s $terminateFakeUID "  if ($terminateFakeUID  > 0);
    $cmd .= "    $uidServer "            if (defined($uidServer));
    $cmd .= " -f $wrk/$asm.frgStore ";
    $cmd .= " -g $wrk/$asm.gkpStore ";
    $cmd .= " -p $asm ";
    $cmd .= " -n $lastckp ";
    $cmd .= " > $distupdate/update.dst ";
    $cmd .= " 2> $distupdate/dumpDistanceEstimates.err ";
    if (runCommand($cmd)) {
        rename "$distupdate/update.dst", "$distupdate/update.dst.FAILED";
        print STDERR "dumpDistanceEstimates Failed.\n";
        exit(1);
    }

    $cmd  = "cd $distupdate && $bin/gatekeeper ";
    $cmd .= " -X -Q -C -P -a $wrk/$asm.gkpStore $distupdate/update.dst ";
    $cmd .= " > gatekeeper.err 2>&1";
    if (runCommand($cmd)) {
        print STDERR "Gatekeeper Failed.\n";
        exit(1);
    }

    touch("$distupdate/distupdate.success");
}


sub resolveSurrogates ($$) {
    my $thisDir = shift @_;
    my $lastDir = shift @_;

    return $thisDir if (-e "$wrk/$thisDir/resolveSurrogates.success");

    my $lastckp = findLastCheckpoint($lastDir);

    system("mkdir $wrk/$thisDir")  if (! -d "$wrk/$thisDir");

    system("ln -s $wrk/$lastDir/$asm.ckp.$lastckp $wrk/$thisDir/$asm.ckp.$lastckp")  if (! -e "$wrk/$thisDir/$asm.ckp.$lastckp");
    system("ln -s $wrk/$asm.SeqStore              $wrk/$thisDir/$asm.SeqStore")      if (! -e "$wrk/$thisDir/$asm.SeqStore");

    my $cmd;
    $cmd  = "cd $wrk/$thisDir && ";
    $cmd .= "$bin/resolveSurrogates ";
    $cmd .= " -f $wrk/$asm.frgStore ";
    $cmd .= " -g $wrk/$asm.gkpStore ";
    $cmd .= " -c $asm ";
    $cmd .= " -n $lastckp ";
    $cmd .= " -1 ";
    $cmd .= " > $wrk/$thisDir/resolveSurrogates.err 2>&1";
    if (runCommand($cmd)) {
        print STDERR "Failed.\n";
        exit(1);
    }

    touch("$wrk/$thisDir/resolveSurrogates.success");

    return $thisDir;
}


sub scaffolder {
    my $lastDir    = undef;
    my $thisDir    = 0;
    my $stoneLevel = getGlobal("throwStones");

    #  If we're not doing eCR, we just do a single scaffolder run, and
    #  get the heck outta here!  OK, we'll do resolveSurrogates() too.
    #
    if (getGlobal("doExtendClearRanges") == 0) {
        $lastDir = CGW("7-CGW", $lastDir, $stoneLevel, undef);
        $thisDir++;
    } else {

        #  Do the initial CGW, making sure to not throw stones.
        #
        $lastDir = CGW("7-$thisDir-CGW", $lastDir, 0, undef);
        $thisDir++;

        #  Followed by at least one eCR, and a distance update.
        #
        $lastDir = eCR("7-$thisDir-ECR", $lastDir);
        $thisDir++;

        updateDistanceRecords($lastDir);

        #  Iterate eCR: do another scaffolder still without stones,
        #  then another eCR.  Again, and again, until we get dizzy and
        #  fall over.
        #
        for (my $rounds = getGlobal("doExtendClearRanges") - 1; $rounds > 0; $rounds--) {
            $lastDir = CGW("7-$thisDir-CGW", $lastDir, 0, 3);
            $thisDir++;

            $lastDir = eCR("7-$thisDir-ECR", $lastDir);
            $thisDir++;

            updateDistanceRecords($lastDir);
        }
    }

    #  Then another scaffolder, chucking stones into the big holes.
    #
    $lastDir = CGW("7-$thisDir-CGW", $lastDir, $stoneLevel, 3);
    $thisDir++;

    #  Then resolve surrogates.
    #
    $lastDir = resolveSurrogates("7-$thisDir-resolveSurrogates", $lastDir);
    $thisDir++;
        
    #  And yet another scaffolder, this time to just do output
    #
    $lastDir = CGW("7-$thisDir-CGW", $lastDir, $stoneLevel, 14);
    $thisDir++;

    #  And, finally, hold on, we're All Done!  Point to the correct output directory.
    #
    system("ln -s $wrk/$lastDir $wrk/7-CGW") if (! -d "$wrk/7-CGW");
}


1;
