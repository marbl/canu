use strict;

#  Assembly all done, toggle the unitigs and re-run CGW and subsequent steps of the assembly.

sub toggler () {
   my $toggledDir = "10-toggledAsm";
   my $ecrEdits = "frg.ECREdits.txt";
   
   return if (-d "$wrk/$toggledDir/$asm.asm");
   return if (getGlobal("doToggle") == 0);

   my $minLength = getGlobal("toggleUnitigLength");
   my $numInstances = getGlobal("toggleNumInstances");
   my $maxDistance = getGlobal("toggleMaxDistance"); 
   my $bin = getBinDirectory();
   my $cmd = "";
   my $scaffoldDir;

   system("mkdir $wrk/$toggledDir") if (! -d "$wrk/$toggledDir");

   # link the stores for space savings
   if (! -e "$wrk/$toggledDir/$asm.ovlStore") {
      system("ln -s $wrk/$asm.ovlStore $wrk/$toggledDir/$asm.ovlStore") if (! -e "$wrk/$toggledDir/$asm.ovlStore");
   }

   if (! -e "$wrk/$toggledDir/$asm.gkpStore") {
      system("mkdir $wrk/$toggledDir/$asm.gkpStore") if (! -d "$wrk/$toggledDir/$asm.gkpStore");
      system("ln -s $wrk/$asm.gkpStore/* $wrk/$toggledDir/$asm.gkpStore") if (! -e "$wrk/$toggledDir/$asm.gkpStore/frg");

      # but the frg store is rewritten by cgw, so reset the ECR clear-ranges
      system("rm -rf $wrk/$toggledDir/$asm.gkpStore/frg");
      system("cp $wrk/$asm.gkpStore/frg $wrk/$toggledDir/$asm.gkpStore/frg");
      
      # back out the ECR changes from the gkp store   
      $cmd  = "$bin/gatekeeper ";
      $cmd .= " -dumpfragments -tabular";
      $cmd .= " -allreads -clear OBT ";
      $cmd .= " $wrk/$asm.gkpStore ";
      $cmd .= " | grep -v \"UID\" ";
      $cmd .= " | awk '{print \"frg uid \"\$1\" ECR1 ALL \"\$12\" \"\$13}' ";
      $cmd .= " > $wrk/$toggledDir/$asm.gkpStore/$ecrEdits 2> $wrk/$toggledDir/$asm.gkpStore/$ecrEdits.err";   
      if (runCommand("$wrk/$toggledDir", $cmd)) {
         caFailure("failed to get pre-ECR clear-ranges for toggling", "$wrk/$toggledDir/$asm.gkpStore/$ecrEdits.err");
      }
      
      $cmd  = "$bin/gatekeeper ";
      $cmd .= " --edit $wrk/$toggledDir/$asm.gkpStore/$ecrEdits";
      $cmd .= " $wrk/$toggledDir/$asm.gkpStore";
      $cmd .= " > $wrk/$toggledDir/$asm.gkpStore/gkpEdit.err 2>&1";
      if (runCommand("$wrk/$toggledDir", $cmd)) {
         caFailure("failed to edit gatekeeper to set ECR clear-ranges for toggling", "$wrk/$toggledDir/$asm.gkpStore/gkpEdit.err");
      }
   }

   system("ln -s $wrk/4-unitigger $wrk/$toggledDir") if (! -e "$wrk/$toggledDir/4-unitigger");
   system("mkdir $wrk/$toggledDir/5-consensus") if (! -d "$wrk/$toggledDir/5-consensus");
   
   my $cgiFile;
   open(F, "ls $wrk/5-consensus |");
   while (<F>) {
      chomp;
      if (m/cgi$/) {
         $cgiFile .= " $wrk/5-consensus/$_";
      }
   }
   close(F);
   
   # create the toggled cgi file
   if (! -e "$wrk/$toggledDir/toggled.success") {
      $cmd  = "$bin/markUniqueUnique ";
      $cmd .= " -a $wrk/9-terminator/$asm.asm ";
      $cmd .= " -l $minLength ";
      $cmd .= " -n $numInstances ";
      $cmd .= " -d $maxDistance ";
      $cmd .= " $cgiFile";
      $cmd .= " > $wrk/$toggledDir/5-consensus/$asm.cgi 2> $wrk/$toggledDir/toggle.err";
      if (runCommand("$wrk/$toggledDir", $cmd)) {
         caFailure("failed to toggle unitigs ", "$wrk/$toggledDir/toggle.err");
      }
      
      touch("$wrk/$toggledDir/toggled.success");
   }

   my $numToggles = `tail -n 1 $wrk/$toggledDir/toggle.err | awk '{print \$2}'`;
   if ($numToggles == 0) {
       print "No toggling occured. Finished.\n";
   }
   else {
      $wrk = "$wrk/$toggledDir";
      $cgiFile = "$wrk/5-consensus/$asm.cgi";

      scaffolder($cgiFile);
      postScaffolderConsensus($scaffoldDir);
      terminate($scaffoldDir);
      cleaner();
   }
}
