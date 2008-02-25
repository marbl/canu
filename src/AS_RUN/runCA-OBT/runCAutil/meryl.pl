use strict;

sub runMeryl ($$$$$$) {
    my $merSize      = shift @_;
    my $merComp      = shift @_;
    my $merCanonical = shift @_;
    my $merThresh    = shift @_;
    my $merType      = shift @_;
    my $merDump      = shift @_;
    my $bin          = getBinDirectory();
    my $cmd;

    #  The fasta file we should be creating.
    my $ffile = "$wrk/0-mercounts/$asm.nmers.$merType.fasta";

    if ($merThresh == 0) {
        touch $ffile;
        return;
    }

    if (-e $ffile) {
        return;
    }

    if (merylVersion() eq "Mighty") {

        #  Use the better meryl!  This is straightforward.  We count,
        #  then we dump.

        #  Intermediate file
        my $ofile = "$wrk/0-mercounts/$asm$merCanonical-ms$merSize-cm$merComp";

        if (! -e "$ofile.mcdat") {
            my $merylMemory = getGlobal("merylMemory");
	    my $merylThreads = getGlobal("merylThreads");

            $cmd .= "$bin/meryl ";
            $cmd .= " -B $merCanonical -v -m $merSize -memory $merylMemory -threads $merylThreads -c $merComp ";
            $cmd .= " -s $wrk/$asm.gkpStore:obt ";
            $cmd .= " -o $ofile ";
            $cmd .= "> $wrk/0-mercounts/meryl.out 2>&1";

            if (runCommand("$wrk/0-mercounts", $cmd)) {
                caFailure("Failed.\n");
            }
        }

        #  We only need the ascii dump if we're doing overlapper, mer
        #  overlapper reads meryl directly.
        #
        if ($merDump) {
            if (! -e $ffile) {
                $cmd  = "$bin/meryl ";
                $cmd .= "-Dt -n $merThresh ";
                $cmd .= "-s $ofile ";
                $cmd .= "> $ffile ";

                if (runCommand("$wrk/0-mercounts", $cmd)) {
                    unlink $ffile;
                    caFailure("Failed.\n");
                }
            }
        }
    } elsif (merylVersion() eq "CA") {

        #  Sigh.  The old meryl.  Not as easy.  If we assume the
        #  process, in particular that the Ovl threshold is less than
        #  the Obt threshold, and that we have already computed the
        #  Ovl mers, we could filter the Ovl mers to get the Obt mers.
        #  But that's tough, especially if we allow mer compression.

        my $merSkip = 10;

        #  Intermediate file
        my $ofile = "$wrk/0-mercounts/$asm-ms$merSize-mt$merThresh-mk$merSkip.$merType.fasta";

        if ($merComp > 0) {
            print STDERR "ERROR!  merCompression not supported without installing kmer\n";
            print STDERR "        (http://sourceforge.net/projects/kmer/).\n";
            print STDERR "If you have installed kmer, then your build is broken, as I\n";
            print STDERR "did not find the correct 'meryl' (meryl -V should have said Mighty).\n";
            die;
        }

        if ($merCanonical ne "-C") {
            print STDERR "ERROR!  mer overlapper not supported without installing kmer\n";
            print STDERR "        (http://sourceforge.net/projects/kmer/).\n";
            print STDERR "If you have installed kmer, then your build is broken, as I\n";
            print STDERR "did not find the correct 'meryl' (meryl -V should have said Mighty).\n";
            die;
        }

        if (! -e $ofile) {
            $merThresh /= $merSkip;

            $cmd  = "$bin/meryl ";
            $cmd .= "-s $wrk/$asm.gkpStore -m $merSize -n $merThresh -K $merSkip ";
            $cmd .= " -o $ofile";
            $cmd .= "> $wrk/0-mercounts/meryl.out 2>&1";

            if (runCommand("$wrk/0-mercounts", $cmd)) {
                unlink $ofile;
                caFailure("Failed.\n");
            }
        }

        symlink($ofile, $ffile) if (! -e $ffile);
    } else {
        caFailure("Unknown meryl version.\n");
    }
}

sub meryl {
    system("mkdir $wrk/0-mercounts") if (! -d "$wrk/0-mercounts");

    if (getGlobal("ovlOverlapper") eq "umd") {
        caFailure("meryl() attempted to compute mer counts for the umd overlapper?\n");
    }

    my $ovlc = 0;  #  No compression, unless we're the mer overlapper
    my $obtc = 0;

    my $ovlC = "-C";  #  Canonical, unless we're the mer overlapper
    my $obtC = "-C";  #  (except the mer overlapper now wants canonical)

    my $ovlD = 1;  #  Dump, unless we're the mer overlapper
    my $obtD = 1;

    if (getGlobal("ovlOverlapper") eq "mer") {
        $ovlc = getGlobal("merCompression");
        $ovlC = "-C";
        $ovlD = 0;
    }
    if (getGlobal("obtOverlapper") eq "mer") {
        $obtc = getGlobal("merCompression");
        $obtC = "-C";
        $obtD = 0;
    }

    runMeryl(getGlobal('ovlMerSize'), $ovlc, $ovlC, getGlobal("ovlMerThreshold"), "ovl", $ovlD);
    runMeryl(getGlobal('obtMerSize'), $obtc, $obtC, getGlobal("obtMerThreshold"), "obt", $obtD);
}

1;
