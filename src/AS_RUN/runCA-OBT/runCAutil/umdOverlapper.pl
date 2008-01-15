use strict;

sub findUMDOverlapper() {
   # double check we can find the executable in path/bin/exec
   caFailure("Unable to find the UMD Overlapper executable.\n") if (! -e "$bin/runUMDOverlapper");
   
   return $bin;
}
sub UMDOverlapper (@) {
    my @fragFiles  = @_;

    caFailure("Either doUseUMDOverlapper or doUseOverlapper must be true.\n") if (getGlobal("doUseUMDOverlapper") == 0 && getGlobal("doUseOverlapper") == 0);
    
    goto alldone if (getGlobal("doUseUMDOverlapper") == 0);
    goto alldone if (-d "$wrk/$asm.ovlStore");

    # when using UMD overlapper, skip OBT (umd overlap has OBT built-in)
    setGlobal("doUseOverlapper", 0);
    setGlobal("doOverlapTrimming", 0);
    setGlobal("doMeryl", 0);

    my $outDir  = "1-overlapper";
    system("mkdir $wrk/$outDir") if (! -d "$wrk/$outDir");
    
    my $umdOvlBin = findUMDOverlapper();

    my $jobID = "0000001";
    system("mkdir $wrk/$outDir/$jobID") if (! -d "$wrk/$outDir/$jobID");

    #cat the frg files together
    my $frgFile = "";
    my $failedFiles = 0;
    
    my $cmd = "";
    if (@fragFiles > 1) {
       $frgFile = "$wrk/$outDir/$asm.frg";
       $cmd = "cat ";
       foreach my $frg (@fragFiles) {
         if (! -e $frg) {
            print STDERR "MISSING: $frg\n";
            next;
            $failedFiles++;
         }
         $cmd .= "$frg ";
       }
       caFailure("Files supplied on command line not found.\n") if ($failedFiles);    
       $cmd .= " > $frgFile";
       runCommand("$wrk/$outDir", $cmd);
    } else {
       $frgFile = $fragFiles[0];
    }

    # create a job list (we have only one job for right now)
    open(SUB, "> $wrk/$outDir/ovljobs.dat") or caFailure("Failed to open '$wrk/$outDir/ovljobs.dat'\n");
    print SUB "$jobID ";   print SUB "\n";
    print SUB "$jobID ";   print SUB "\n";
    close(SUB);
    
    # run frg file command
    #
    $cmd = "";
   
    $cmd  = "$umdOvlBin/runUMDOverlapper ";
    $cmd .= getGlobal("umdOverlapperFlags") . " ";
    $cmd .= "$frgFile $wrk/$outDir/$jobID/$asm.umd.frg ";
    $cmd .= " > $wrk/$outDir/$jobID/overlapper.out 2>$wrk/$outDir/$jobID/overlapper.err";
   
    if (runCommand("$wrk/$outDir", $cmd)) {
      caFailure("Failed to run UMD overlapper.\n");
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
 
    @fragFiles = ( "$wrk/$outDir/$jobID/$asm.umd.frg" );
    print STDERR "The frag files are now @fragFiles\n";
   
    touch("$wrk/$outDir/$jobID/$jobID.success");
    stopAfter("overlapper");
   
  alldone:    
    return @fragFiles;
}

1;
