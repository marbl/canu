use strict;

#  should we require stoneLevel be passed into scaffolder?

sub CGW ($$$$) {
    my $thisDir    = shift @_;
    my $lastDir    = shift @_;
    my $stoneLevel = shift @_;
    my $logickp    = shift @_;

    return $thisDir if (-e "$wrk/$thisDir/cgw.success");

    my $lastckp = findLastCheckpoint($lastDir)  if (defined($lastDir));
    my $ckp     = "-R $lastckp -N $logickp"     if (defined($lastckp) && defined($logickp));

    system("mkdir $wrk/$thisDir")               if (! -d "$wrk/$thisDir");
    system("mkdir $wrk/$asm.SeqStore")          if (! -d "$wrk/$asm.SeqStore");

    system("ln -s $wrk/5-consensus/$asm.cgi $wrk/$thisDir/$asm.cgi")          if (! -e "$wrk/$thisDir/$asm.cgi");
    system("ln -s $wrk/$asm.SeqStore        $wrk/$thisDir/$asm.SeqStore")     if (! -e "$wrk/$thisDir/$asm.SeqStore");

    system("ln -s $wrk/$lastDir/$asm.ckp.$lastckp $wrk/$thisDir/$asm.ckp.$lastckp") if (defined($lastDir));

    my $cmd;
    $cmd  = "cd $wrk/$thisDir && ";
    $cmd .= "$bin/cgw $ckp -c -j 1 -k 5 -r 4 -s $stoneLevel -w 0 -T -P ";
    $cmd .= "-f $wrk/$asm.frgStore ";
    $cmd .= "-g $wrk/$asm.gkpStore ";
    $cmd .= "-o $wrk/$thisDir/$asm ";
    $cmd .= "$wrk/$thisDir/$asm.cgi ";
    $cmd .= "> $wrk/$thisDir/cgw.out ";
    $cmd .= "2> $wrk/$thisDir/cgw.err";
        
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

    my $cmd;
    $cmd  = "cd $wrk/$thisDir && ";
    $cmd .= "$bin/extendClearRanges ";
    $cmd .= "-f $wrk/$asm.frgStore ";
    $cmd .= "-g $wrk/$asm.gkpStore ";
    $cmd .= "-c $asm ";
    $cmd .= "-n $lastckp ";
    $cmd .= "-s -1 ";
    $cmd .= "> $wrk/$thisDir/extendClearRanges.out ";
    $cmd .= "2> $wrk/$thisDir/extendClearRanges.err ";

    if (runCommand($cmd)) {
        print STDERR "Failed.\n";
        exit(1);
    }

    touch("$wrk/$thisDir/extendClearRanges.success");

    return $thisDir;
}


sub updateDistanceRecords ($) {
    my $thisDir = shift @_;

    return if (getGlobal("doUpdateDistanceRecords", 0) == 0);
    return if (-e "$wrk/$thisDir/distupdate/distupdate.success");

    my $lastckp    = findLastCheckpoint($thisDir);
    my $distupdate = "$wrk/$thisDir/distupdate";

    system("mkdir $distupdate") if (! -d "$distupdate");

    system("ln -s $wrk/$thisDir/$asm.ckp.$lastckp $distupdate/$asm.ckp.$lastckp")  if (! -e "$distupdate/$asm.ckp.$lastckp");
    system("ln -s $wrk/$asm.SeqStore              $distupdate/$asm.SeqStore")      if (! -e "$distupdate/$asm.SeqStore");

    my $cmd;
    $cmd  = "cd $distupdate && $bin/dumpDistanceEstimates ";
    $cmd .= " -f $wrk/$asm.frgStore ";
    $cmd .= " -g $wrk/$asm.gkpStore ";
    $cmd .= " -p $asm ";
    $cmd .= " -n $lastckp ";
    $cmd .= " > $distupdate/update.dst ";
    $cmd .= "2> $distupdate/update.err ";
    if (runCommand($cmd)) {
        rename "$distupdate/update.dst", "$distupdate/update.dst.FAILED";
        print STDERR "dumpDistanceEstimates Failed.\n";
        exit(1);
    }

    $cmd  = "cd $distupdate && $bin/gatekeeper ";
    $cmd .= " -X -Q -C -P -a $wrk/$asm.gkpStore $distupdate/update.dst";
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
    $cmd .= "-f $wrk/$asm.frgStore ";
    $cmd .= "-g $wrk/$asm.gkpStore ";
    $cmd .= "-c $asm ";
    $cmd .= "-n $lastckp ";
    $cmd .= "-1 ";
    $cmd .= "> $wrk/$thisDir/resolveSurrogates.out ";
    $cmd .= "2> $wrk/$thisDir/resolveSurrogates.err ";

    if (runCommand($cmd)) {
        print STDERR "Failed.\n";
        exit(1);
    }

    touch("$wrk/$thisDir/resolveSurrogates.success");

    return $thisDir;
}


sub scaffolder {
    my $lastDir = undef;
    my $thisDir = 0;


    #  If we're not doing eCR, we just do a single scaffolder run, and
    #  get the heck outta here!
    #
    if ($doExtendClearRanges == 0) {
        CGW("7-CGW", $lastDir, $stoneLevel, undef);
        return;
    }


    #  Otherwise, we want to do eCR (and surrogates and a bunch of
    #  scaffolder runs), so run the first scaffolder WITHOUT stones.
    #
    $lastDir = CGW("7-$thisDir-CGW", $lastDir, 0, undef);
    $thisDir++;


    #  Then do a round of eCR
    #
    $lastDir = eCR("7-$thisDir-ECR", $lastDir);
    $thisDir++;

    updateDistanceRecords($lastDir);


    #  If we're iterating over eCR, do another scaffolder WITHOUT
    #  stones, then another eCR.  Again, and again.
    #
    for (my $eCRRounds = getGlobal("eCRRounds", 1)-1; $eCRRounds > 0; $eCRRounds--) {
        $lastDir = CGW("7-$thisDir-CGW", $lastDir, 0, 3);
        $thisDir++;

        $lastDir = eCR("7-$thisDir-ECR", $lastDir);
        $thisDir++;

        updateDistanceRecords($lastDir);
    }


    #  Finally, finish it off with a scaffolder from checkpoint level
    #  3, this time throwing stones!
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


    #  And, finally, we're All Done!  Point to the correct output directory.
    #
    system("ln -s $wrk/$lastDir $wrk/7-CGW") if (! -d "$wrk/7-CGW");
}


1;
