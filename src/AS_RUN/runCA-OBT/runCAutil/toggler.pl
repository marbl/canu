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

   system("mkdir $wrk/$toggledDir") if (! -d "$wrk/$toggledDir");

   # link the stores for space savings
   if (! -e "$wrk/$toggledDir/$asm.ovlStore") {
      system("ln -s $wrk/$asm.ovlStore $wrk/$toggledDir/$asm.ovlStore") if (! -e "$wrk/$toggledDir/$asm.ovlStore");
   }

   if (! -e "$wrk/$toggledDir/$asm.gkpStore") {
      system("mkdir $wrk/$toggledDir/$asm.gkpStore") if (! -d "$wrk/$toggledDir/$asm.gkpStore");
      system("ln -s $wrk/$asm.gkpStore/* $wrk/$toggledDir/$asm.gkpStore") if (! -e "$wrk/$toggledDir/$asm.gkpStore/frg");

      # but the frg store is rewritten by cgw, so reset the ECR clear-ranges
      system("rm -rf $wrk/$toggledDir/$asm.gkpStore/*00*");
      system("rm -rf $wrk/$toggledDir/$asm.gkpStore/fnm");
      system("cp $wrk/$asm.gkpStore/fnm $wrk/$toggledDir/$asm.gkpStore/fnm");
      
      # back out the ECR changes from the gkp store   
      $cmd  = "$bin/gatekeeper ";
      $cmd .= " --revertclear OBTCHIMERA $wrk/$toggledDir/$asm.gkpStore";
      $cmd .= " > $wrk/$toggledDir/$asm.gkpStore/$ecrEdits.err 2> $wrk/$toggledDir/$asm.gkpStore/$ecrEdits.err";   
      if (runCommand("$wrk/$toggledDir", $cmd)) {
         caFailure("failed to get pre-ECR clear-ranges for toggling", "$wrk/$toggledDir/$asm.gkpStore/$ecrEdits.err");
      }
   }

   system("ln -s $wrk/4-unitigger $wrk/$toggledDir") if (! -e "$wrk/$toggledDir/4-unitigger");
   system("ln -s $wrk/5-consensus $wrk/$toggledDir") if (! -e "$wrk/$toggledDir/5-consensus");
   
   # copy the tigStore
   if (! -e "$wrk/$toggledDir/$asm.tigStore") {
      system("mkdir $wrk/$toggledDir/$asm.tigStore") ;
      system("cp -rf $wrk/$asm.tigStore/*v001* $wrk/$toggledDir/$asm.tigStore");
      system("cp -rf $wrk/$asm.tigStore/*v002* $wrk/$toggledDir/$asm.tigStore");
   }
   
   # create the toggled cgi file
   if (! -e "$wrk/$toggledDir/toggled.success") {
      $cmd  = "$bin/markUniqueUnique ";
      $cmd .= " -a $wrk/9-terminator/$asm.asm ";
      $cmd .= " -l $minLength ";
      $cmd .= " -n $numInstances ";
      $cmd .= " -d $maxDistance ";
      $cmd .= " $wrk/$toggledDir/$asm.tigStore";
      $cmd .= " > $wrk/$toggledDir/toggle.err 2> $wrk/$toggledDir/toggle.err";
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

      scaffolder();
      postScaffolderConsensus();
      terminate();
      cleaner();
   }
}
