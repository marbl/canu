use strict;

#  Assembly all done, toggle the unitigs and re-run CGW and subsequent steps of the assembly.

sub toggler () {
    my $toggledDir = "10-toggledAsm";
    my $ecrEdits = "";

    return if (-d "$wrk/$toggledDir/$asm.asm");
    return if (getGlobal("doToggle") == 0);

    my $minLength = getGlobal("toggleUnitigLength");
    my $numInstances = getGlobal("toggleNumInstances");
    my $maxDistance = getGlobal("toggleMaxDistance");

    my $bin = getBinDirectory();

    system("mkdir $wrk/$toggledDir") if (! -d "$wrk/$toggledDir");

    #  A simple link to the ovlStore suffices.
    #
    if (! -e "$wrk/$toggledDir/$asm.ovlStore") {
        system("ln -s $wrk/$asm.ovlStore $wrk/$toggledDir/$asm.ovlStore");
    }

    #  The gatekeeper store must be paritally copied, so that we can first undo
    #  clear range changes made by any previous ECR, and allow clear range changes
    #  to be made by future ECR.
    #
    if (! -e "$wrk/$toggledDir/$asm.gkpStore") {
        system("mkdir $wrk/$toggledDir/$asm.gkpStore");

        system("ln -s $wrk/$asm.gkpStore/[qsu]?? $wrk/$toggledDir/$asm.gkpStore/");
        system("cp $wrk/$asm.gkpStore/[filp]??   $wrk/$toggledDir/$asm.gkpStore/");
        system("cp $wrk/$asm.gkpStore/clr-*      $wrk/$toggledDir/$asm.gkpStore/");

        my $cmd;
        $cmd  = "$bin/gatekeeper";
        $cmd .= " --revertclear OBTCHIMERA $wrk/$toggledDir/$asm.gkpStore";
        $cmd .= " > $wrk/$toggledDir/$asm.gkpStore.resetClearRange.err 2>&1";

        if (runCommand("$wrk/$toggledDir", $cmd)) {
            caFailure("failed to get pre-ECR clear-ranges for toggling", "$wrk/$toggledDir/$asm.gkpStore.resetClearRange.err");
        }
    }

    #  The tigStore needs only a partial copy, and links suffice.
    #
    if (! -e "$wrk/$toggledDir/$asm.tigStore") {
        system("mkdir $wrk/$toggledDir/$asm.tigStore") ;

        system("ln -s $wrk/$asm.tigStore/*v00[12]* $wrk/$toggledDir/$asm.tigStore/");
    }

    system("ln -s $wrk/4-unitigger $wrk/$toggledDir") if (! -e "$wrk/$toggledDir/4-unitigger");
    system("ln -s $wrk/5-consensus $wrk/$toggledDir") if (! -e "$wrk/$toggledDir/5-consensus");

    #  Update the tigStore, flipping repeat untigs to unique unitigs.
    #
    if (! -e "$wrk/$toggledDir/toggled.success") {
        my $cmd;
        $cmd  = "$bin/markUniqueUnique ";
        $cmd .= " -a $wrk/9-terminator/$asm.asm ";
        $cmd .= " -t $wrk/$toggledDir/$asm.tigStore 2 ";
        $cmd .= " -l $minLength ";
        $cmd .= " -n $numInstances ";
        $cmd .= " -d $maxDistance ";
        $cmd .= " > $wrk/$toggledDir/toggle.err 2>&1";

        if (runCommand("$wrk/$toggledDir", $cmd)) {
            caFailure("failed to toggle unitigs ", "$wrk/$toggledDir/toggle.err");
        }

        touch("$wrk/$toggledDir/toggled.success");
    }

    my $numToggles = `tail -n 1 $wrk/$toggledDir/toggle.err | awk '{print \$2}'`;

    if ($numToggles == 0) {
        print "No toggling occured. Finished.\n";
        return;
    }

    my $oldwrk = $wrk;

    $wrk = "$wrk/$toggledDir";

    scaffolder();
    postScaffolderConsensus();
    terminate();
    cleaner();

    $wrk = $oldwrk;
}
