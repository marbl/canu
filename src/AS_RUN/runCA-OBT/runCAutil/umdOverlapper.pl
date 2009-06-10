use strict;

sub getUMDOverlapperClearRange ($) {
    my $dir     = shift @_;
    my $fileName = "$asm.obtClrRange";

    open(F, "ls -1 -d $wrk/$dir/*overlapperRunDir* |");
    open(G, ">$wrk/$dir/$fileName") or caFailure("failed to write '$wrk/$dir/$fileName'", undef);
    while (<F>) {
        chomp;

        open(T, "< $_/revisedOrigTrimsForReads.txt") or caFailure("failed to open '$_/revisedOrigTrimsForReads.txt'", undef);
        while (<T>) {
           my @trimData = split(/\s+/,$_);
           my $uid = $trimData[0];
           my $bgn = $trimData[1];
           my $end = $trimData[2];

           if ($bgn < $end) {
             print G "frg uid $uid obt all $bgn $end\n";
           } else {
             print G "frg uid $uid obt all $end $bgn\n";
           }
        }
        close(T);
    }
    close(F);
    close(G);

    return $fileName;
}

sub UMDoverlapper () {
    goto alldone if (-d "$wrk/$asm.ovlStore");
    goto alldone if (getGlobal("ovlOverlapper") ne "umd");

    my $outDir  = "1-overlapper";
    system("mkdir $wrk/$outDir") if (! -d "$wrk/$outDir");

    my $jobID = "0000001";
    system("mkdir $wrk/$outDir/$jobID") if (! -d "$wrk/$outDir/$jobID");

    my $vi = getGlobal("vectorIntersect");

    my $bin = getBinDirectory();

    #dump the frag file from gkp if it does not exist already
    # should check if vector clear then dump vec range else dump this range
    if (defined($vi)) {
       if (runCommand($wrk, "$bin/gatekeeper -clear VEC -dumpfrg $wrk/$asm.gkpStore 2> $wrk/gatekeeper.err | grep -v 'No source' > $wrk/$asm.vec.frg")) {
          caFailure("failed to dump gatekeeper store for UMD overlapper", "$wrk/gatekeeper.err");
       }
    }
    elsif ( ! -s "$wrk/$asm.frg" ) {
       if (runCommand($wrk, "$bin/gatekeeper -dumpfrg $wrk/$asm.gkpStore 2> $wrk/gatekeeper.err | grep -v 'No source' > $wrk/$asm.frg")) {
          caFailure("failed to dump gatekeeper store for UMD overlapper", "$wrk/gatekeeper.err");
       }
    }

    # create a job list (we have only one job for right now)
    open(SUB, "> $wrk/$outDir/ovljobs.dat") or caFailure("failed to open '$wrk/$outDir/ovljobs.dat'", undef);
    print SUB "$jobID ";   print SUB "\n";
    print SUB "$jobID ";   print SUB "\n";
    close(SUB);

    # run frg file command
    #
    my $cmd  = "$bin/runUMDOverlapper ";
    $cmd .= getGlobal("umdOverlapperFlags") . " ";

    # when we have vector clear, pass it to the overlapper, otherwise tell the overlapper to figure it out
    if (defined($vi)) {
       $cmd .= "-vector-trim-file $wrk/$asm.vec.frg $wrk/$asm.vec.frg "
    } else {
       $cmd .= "-calculate-trims $wrk/$asm.frg ";
    }

    $cmd .= "$wrk/$outDir/$jobID/$asm.umd.frg ";
    $cmd .= " > $wrk/$outDir/$jobID/overlapper.out 2>$wrk/$outDir/$jobID/overlapper.err";

    if (runCommand("$wrk/$outDir", $cmd)) {
      caFailure("failed to run UMD overlapper", "$wrk/$outDir/$jobID/overlapper.err");
    }

    my $trimFile = getUMDOverlapperClearRange($outDir);
    $cmd = "";
    $cmd .= "$bin/gatekeeper --edit ";
    $cmd .= "$wrk/$outDir/$trimFile $wrk/$asm.gkpStore";
    if (runCommand("$wrk/$outDir", $cmd)) {
      caFailure("failed to update OBT trims", "undef");
    }

    # now create the binary overlaps
    $cmd = "";
    $cmd .= "cat $wrk/$outDir/$jobID/$asm.umd.reliable.overlaps | ";
    $cmd .= "awk '{print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6\"\\t\"\$7}' | ";
    $cmd .= "$bin/convertOverlap ";
    $cmd .= "-b -ovldump ";
    $cmd .= " > $wrk/$outDir/$jobID/$jobID.ovb";
    if (runCommand("$wrk/$outDir", $cmd)) {
      caFailure("failed to create overlaps", undef);
    }

    #cleanup
    rmrf("$asm.vec.frg");

    touch("$wrk/$outDir/$jobID/$jobID.success");
    stopAfter("overlapper");

  alldone:
}

1;
