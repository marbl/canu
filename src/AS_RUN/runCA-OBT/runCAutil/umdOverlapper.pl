use strict;

sub findUMDOverlapper() {
   # double check we can find the executable in path/bin/exec
   caFailure("Unable to find the UMD Overlapper executable.\n") if (! -e "$bin/runUMDOverlapper");
   
   return $bin;
}

sub getUMDClearRange ($) {
    my $dir     = shift @_;
    my $fileName = "$asm.obtClrRange";

    open(F, "ls -1 -d $wrk/$dir/*overlapperRunDir* |");    
    open(G, ">$wrk/$dir/$fileName") or caFailure("Failed to write '$wrk/$dir/$fileName'\n");    
    while (<F>) {
        chomp;

        open(T, "< $_/revisedOrigTrimsForReads.txt") or caFailure("Failed to open '$_/revisedOrigTrimsForReads.txt'\n");
        while (<T>) {
           my @trimData = split(/\s+/,$_);
           my $uid = $trimData[0];
           my $bgn = $trimData[1];
           my $end = $trimData[2];
           
           if ($bgn < $end) {
             print G "frg uid $uid obt $bgn $end\n";
           } else {
             print G "frg uid $uid obt $end $bgn\n";
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
    
    my $umdOvlBin = findUMDOverlapper();

    my $jobID = "0000001";
    system("mkdir $wrk/$outDir/$jobID") if (! -d "$wrk/$outDir/$jobID");

    #dump the frag file from gkp if it does not exist already
    if ( ! -s "$wrk/$asm.frg" ) {
       if (runCommand($wrk, "$bin/gatekeeper -dumpfrg $wrk/$asm.gkpStore 2> $wrk/gatekeeper.err | grep -v 'No source' > $wrk/$asm.frg")) {
          caFailure("Failed to dump gatekeeper store for UMD overlapper");
       }    
    }
    
    # create a job list (we have only one job for right now)
    open(SUB, "> $wrk/$outDir/ovljobs.dat") or caFailure("Failed to open '$wrk/$outDir/ovljobs.dat'\n");
    print SUB "$jobID ";   print SUB "\n";
    print SUB "$jobID ";   print SUB "\n";
    close(SUB);
    
    # run frg file command
    #
    my $cmd  = "$umdOvlBin/runUMDOverlapper ";
    $cmd .= getGlobal("umdOverlapperFlags") . " ";
    $cmd .= "$wrk/$asm.frg $wrk/$outDir/$jobID/$asm.umd.frg ";
    $cmd .= " > $wrk/$outDir/$jobID/overlapper.out 2>$wrk/$outDir/$jobID/overlapper.err";
   
    if (runCommand("$wrk/$outDir", $cmd)) {
      caFailure("Failed to run UMD overlapper.\n");
    }   

    # update the gkpStore with newly computed clear ranges
    my $trimFile = getUMDClearRange($outDir);
    $cmd = "";
    $cmd .= "$bin/gatekeeper --edit ";
    $cmd .= "$wrk/$outDir/$trimFile $wrk/$asm.gkpStore";
    if (runCommand("$wrk/$outDir", $cmd)) {
      caFailure("Failed to update OBT trims.\n");
    }   
    
    # now create the binary overlaps
    $cmd = "";
    $cmd .= "cat $wrk/$outDir/$jobID/$asm.umd.reliable.overlaps | ";
    $cmd .= "awk '{print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6\"\\t\"\$7}' | ";
    $cmd .= "$bin/convertOverlap ";
    $cmd .= "-b -ovldump ";
    $cmd .= " > $wrk/$outDir/$jobID/$jobID.ovb";
    if (runCommand("$wrk/$outDir", $cmd)) {
      caFailure("Failed to create overlaps .\n");
    }   
        
    touch("$wrk/$outDir/$jobID/$jobID.success");
    stopAfter("overlapper");
   
  alldone:
}

1;
