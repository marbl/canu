use strict;

#  aka, CGW

sub scaffolder {
    my $dir = shift @_;
    my $rkp = shift @_;
    my $lkp = shift @_;
    my $ckp = undef;

    if (!defined($dir)) {
        $dir = "7-CGW";
    }
    if (defined($rkp) && defined($lkp)) {
        $ckp = "-R $rkp -N $lkp";
    }

    system("mkdir $wrk/$dir") if (! -d "$wrk/$dir");

    if (! -e "$wrk/$dir/$asm.cgi") {
        system("ln -s $wrk/5-consensus/$asm.cgi $wrk/$dir/$asm.cgi");
    }

    if (! -e "$wrk/$dir/cgw.success") {
        print STDERR "Starting f -- CGW\n";

        my $cmd;
        $cmd  = "cd $wrk/$dir && ";
        $cmd .= "$bin/cgw $ckp -c -j 1 -k 5 -r 4 -s $stoneLevel -w 0 -T -P ";
        $cmd .= "-f $wrk/$asm.frgStore ";
        $cmd .= "-g $wrk/$asm.gkpStore ";
        $cmd .= "-o $wrk/$dir/$asm ";
        $cmd .= "$wrk/$dir/$asm.cgi ";
        $cmd .= "> $wrk/$dir/cgw.out ";
        $cmd .= "2> $wrk/$dir/cgw.err";
        
        if (runCommand($cmd)) {
            print STDERR "Failed.\n";
            exit(1);
        }

        touch("$wrk/$dir/cgw.success");
    }
}


sub extendClearRanges {

    ########################################
    #
    #  Run the first scaffolder, WITHOUT STONES ALWAYS!
    #
    if (! -e "$wrk/7-0-CGW/cgw.success") {
        my $stoneLevelSave = $stoneLevel;
        $stoneLevel = 0;
        scaffolder("7-0-CGW");
        $stoneLevel = $stoneLevelSave;
    }

    ########################################

    if (! -e "$wrk/7-1-ECR/extendClearRanges.success") {
        my $lastckp = findLastCheckpoint("7-0-CGW");

        system("mkdir $wrk/7-1-ECR") if (! -d "$wrk/7-1-ECR");
        system("ln -s $wrk/7-0-CGW/$asm.ckp.$lastckp $wrk/7-1-ECR/$asm.ckp.$lastckp");
        system("ln -s $wrk/7-0-CGW/$asm.SeqStore     $wrk/7-1-ECR/$asm.SeqStore");

        my $cmd;
        $cmd  = "cd $wrk/7-1-ECR && ";
        $cmd .= "$bin/extendClearRanges ";
        $cmd .= "-f $wrk/$asm.frgStore ";
        $cmd .= "-g $wrk/$asm.gkpStore ";
        $cmd .= "-c $asm ";
        $cmd .= "-n $lastckp ";
        $cmd .= "-s -1 ";
        $cmd .= "> $wrk/7-1-ECR/extendClearRanges.out ";
        $cmd .= "2> $wrk/7-1-ECR/extendClearRanges.err ";

        if (runCommand($cmd)) {
            print STDERR "Failed.\n";
            exit(1);
        }

        touch("$wrk/7-1-ECR/extendClearRanges.success");
    }


    ########################################

    if (! -e "$wrk/7-2-CGW/cgw.success") {
        my $lastckp = findLastCheckpoint("7-1-ECR");

        system("mkdir $wrk/7-2-CGW") if (! -d "$wrk/7-2-CGW");
        system("ln -s $wrk/7-1-ECR/$asm.ckp.$lastckp $wrk/7-2-CGW/$asm.ckp.$lastckp");
        system("ln -s $wrk/7-0-CGW/$asm.SeqStore     $wrk/7-2-CGW/$asm.SeqStore");
   
        scaffolder("7-2-CGW", findLastCheckpoint("7-1-ECR"), 3);
    }

    ########################################
    #
    #
    #
    if (! -e "$wrk/7-3-resolveSurrogates/resolveSurrogates.success") {
        my $lastckp = findLastCheckpoint("7-2-CGW");

        system("mkdir $wrk/7-3-resolveSurrogates") if (! -d "$wrk/7-3-resolveSurrogates");
        system("ln -s $wrk/7-2-CGW/$asm.ckp.$lastckp $wrk/7-3-resolveSurrogates/$asm.ckp.$lastckp");
        system("ln -s $wrk/7-0-CGW/$asm.SeqStore     $wrk/7-3-resolveSurrogates/$asm.SeqStore");

        my $cmd;
        $cmd  = "cd $wrk/7-3-resolveSurrogates && ";
        $cmd .= "$bin/resolveSurrogates ";
        $cmd .= "-f $wrk/$asm.frgStore ";
        $cmd .= "-g $wrk/$asm.gkpStore ";
        $cmd .= "-c $asm ";
        $cmd .= "-n $lastckp ";
        $cmd .= "-1 ";
        $cmd .= "> $wrk/7-3-resolveSurrogates/resolveSurrogates.out ";
        $cmd .= "2> $wrk/7-3-resolveSurrogates/resolveSurrogates.err ";

        if (runCommand($cmd)) {
            print STDERR "Failed.\n";
            exit(1);
        }

        touch("$wrk/7-3-resolveSurrogates/resolveSurrogates.success");
    }



    ########################################
    #
    #  Another scaffolder, this time to just do output
    #
    if (! -e "$wrk/7-4-CGW/cgw.success") {
        my $lastckp = findLastCheckpoint("7-3-resolveSurrogates");

        system("mkdir $wrk/7-4-CGW") if (! -d "$wrk/7-4-CGW");
        system("ln -s $wrk/7-3-resolveSurrogates/$asm.ckp.$lastckp $wrk/7-4-CGW/$asm.ckp.$lastckp");
        system("ln -s $wrk/7-0-CGW/$asm.SeqStore                   $wrk/7-4-CGW/$asm.SeqStore");

        scaffolder("7-4-CGW", $lastckp, 14);
    }


    ########################################
    #
    #  All done, point to the correct output
    #
    if (! -d "$wrk/7-CGW") {
        system("ln -s $wrk/7-4-CGW $wrk/7-CGW");
    }
}




1;
