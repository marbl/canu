use strict;

sub getFigaroClearRange ($) {
    my $outDir     = shift @_;
    my $fileName = "$asm.clv";

    # the figaro output is UID,IID CLR_BGN
    # first reformat is as UID CLR_BGN
    runCommand("$wrk/$outDir", "awk '{print substr(\$1, 1, index(\$1, \",\")-1)\" \"\$2}' $wrk/$outDir/$asm.vectorcuts > $wrk/$outDir/$asm.clrBgn");

    # sort by UID and join it together with the read end to form the full vector clear range
    runCommand("$wrk/$outDir", "sort -nk 1 -T $wrk/$outDir $wrk/$outDir/$asm.clrBgn > $wrk/$outDir/$asm.clrBgn.sorted");
    runCommand("$wrk/$outDir", "join $wrk/$outDir/$asm.clrBgn.sorted $wrk/$asm.untrimmed -o 1.1,1.2,2.3 > $wrk/$outDir/$fileName");

    # clean up
    rmrf("$outDir/$asm.clrBgn");
    rmrf("$outDir/$asm.clrBgn.sorted");

    return $fileName;
}

sub generateFigaroTrim($) {
    my $outDir = shift @_;
    my $bin = getBinDirectory();

    return if (-e "$wrk/$outDir/trim.success");

    # run command
    #
    my $cmd  = "$bin/figaro ";
    $cmd .= getGlobal("figaroFlags") . " ";
    $cmd .= "-F $wrk/$asm.fasta -P $asm ";
    $cmd .= " > $wrk/$outDir/figaro.out 2>$wrk/$outDir/figaro.err";

    if (runCommand("$wrk/$outDir", $cmd)) {
      caFailure("figaro died", "$wrk/$outDir/figaro.err");
    }

    # update the gkpStore with newly computed clear ranges
    return getFigaroClearRange($outDir);
}

sub getUMDTrimClearRange($) {
   my $outDir = shift @_;
   my $fileName = "$asm.clv";

   # the umd output is CLR_BGN (in the same order as the input)
   # to join it with the UID we first number both the list of UIDs in the fasta file and the CLR_BGN
   runCommand("$wrk/$outDir", "cat $wrk/$asm.fasta | grep \">\" | awk '{print NR\" \"substr(\$1, 2, index(\$1, \",\")-2)}' > $wrk/$outDir/$asm.numberedUids");
   runCommand("$wrk/$outDir", "awk '{print NR\" \"\$0}' $asm.vectorcuts > $asm.numberedCuts");

   # now we join them together
   runCommand("$wrk/$outDir", "join $wrk/$outDir/$asm.numberedUids $wrk/$outDir/$asm.numberedCuts -o 1.2,2.2 > $wrk/$outDir/$asm.clrBgn");

   # now we can join together the UID CLR_BGN with the read-end information for the full clear range
   runCommand("$wrk/$outDir", "sort -nk 1 -T $wrk/$outDir $wrk/$outDir/$asm.clrBgn > $wrk/$outDir/$asm.clrBgn.sorted");
   runCommand("$wrk/$outDir", "join $wrk/$outDir/$asm.clrBgn.sorted $wrk/$asm.untrimmed -o 1.1,1.2,2.3 > $wrk/$outDir/$fileName");

   # clean up
   rmrf("$outDir/$asm.numberedUids");
   rmrf("$outDir/$asm.numberedCuts");
   rmrf("$outDir/$asm.clrBgn");
   rmrf("$outDir/$asm.clrBgn.sorted");
   rmrf("$outDir/vectorTrimIntermediateFile001.*");

   return $fileName;
}

sub generateUMDTrim($) {
    my $outDir = shift @_;
    my $bin = getBinDirectory();

    return if (-e "$wrk/$outDir/trim.success");

    # run command
    #
    my $cmd  = "$bin/dataWorkReduced/findVectorTrimPoints.perl ";
    $cmd .= "$wrk/$asm.fasta $wrk/$outDir/$asm.vectorcuts ";
    $cmd .= " > $wrk/$outDir/umd.out 2>$wrk/$outDir/umd.err";

    if (runCommand("$wrk/$outDir", $cmd)) {
      caFailure("UMD overlapper dataWorkReduced/findVectorTrimPoints.perl died",
                "$wrk/$outDir/umd.err");
    }

    return getUMDTrimClearRange($outDir);
}

sub generateVectorTrim ($) {
    my $vi = getGlobal("vectorIntersect");
    my $trimmer = getGlobal("vectorTrimmer");
    my $outDir  = "0-preoverlap";
    my $bin = getBinDirectory();
    my $trimFile = undef;

    # when vector insersect is specified or no external trimming is requested, do nothing
    return if (defined($vi));
    return if ($trimmer eq "ca");
    return if (-e "$wrk/$outDir/trim.success");

    #dump the fasta file from gkp
    if ( ! -e "$wrk/$asm.fasta" ) {
       if (runCommand($wrk, "$bin/gatekeeper -dumpfastaseq -clear UNTRIM $wrk/$asm.gkpStore 2> $wrk/$outDir/gatekeeper.err > $wrk/$asm.fasta")) {
           caFailure("failed to dump gatekeeper store for figaro trimmer",
                     "$wrk/$outDir/gatekeeper.err");
       }
    }
    #dump the clr range
    if ( ! -e "$wrk/$asm.untrimmed" ) {
       if (runCommand($wrk, "$bin/gatekeeper -dumpfragments -tabular -clear UNTRIM $wrk/$asm.gkpStore 2> $wrk/$outDir/gatekeeper.err | grep -v 'UID' |awk '{print \$1\" \"\$12\" \"\$13}' | sort -nk 1 -T $wrk/ > $wrk/$asm.untrimmed")) {
           caFailure("failed to dump gatekeeper quality trim points for figaro trimmer",
                     "$wrk/$outDir/gatekeeper.err");
       }
    }

    if ($trimmer eq "figaro") {
       $trimFile = generateFigaroTrim($outDir);
    } elsif($trimmer eq "umd") {
       $trimFile = generateUMDTrim($outDir);
    } else {
       caFailure("unknown vector trimmer $trimmer", undef);
    }

    # set the global vector trim file so that the subsequent code will update the gkp for us
    setGlobal("vectorIntersect", "$wrk/$outDir/$trimFile");

    #cleanup
    rmrf("$asm.fasta");
    rmrf("$asm.untrimmed");

    touch("$wrk/$outDir/trim.success");

    return;
}

1;
