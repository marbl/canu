use strict;

sub runMeryl ($$$$$$) {
    my $merSize      = shift @_;
    my $merComp      = shift @_;
    my $merCanonical = shift @_;
    my $merThresh    = shift @_;
    my $merScale     = 1.0;
    my $merType      = shift @_;
    my $merDump      = shift @_;

    my $bin          = getBinDirectory();
    my $cmd;

    #  The fasta file we should be creating.
    my $ffile = "$wrk/0-mercounts/$asm.nmers.$merType.fasta";

    if ($merThresh =~ m/auto\s*\*\s*(\S+)/) {
        $merThresh = "auto";
        $merScale  = $1;
    }

    if ($merThresh =~ m/auto\s*\/\s*(\S+)/) {
        $merThresh = "auto";
        $merScale  = 1.0 / $1;
    }

    if (($merThresh ne "auto") && ($merThresh == 0)) {
        touch $ffile;
        return;
    }

    #if (-e $ffile) {
    #    print STDERR "runMeryl() would have returned.\n";
    #}
    #print STDERR "$ffile\n";

    if (merylVersion() eq "Mighty") {

        #  Use the better meryl!  This is straightforward.  We count,
        #  then we dump.

        #  Intermediate file
        my $ofile = "$wrk/0-mercounts/$asm$merCanonical-ms$merSize-cm$merComp";

        if (! -e "$ofile.mcdat") {
            my $merylMemory  = getGlobal("merylMemory");
	    my $merylThreads = getGlobal("merylThreads");

            if ($merylMemory !~ m/^-/) {
                $merylMemory = "-memory $merylMemory";
            }

            #  A small optimization we could do if (a) not mer
            #  overlapper, (b) not auto threshold: only save mer
            #  counts above the smaller (of obt & ovl thresholds).
            #  It's complicated, and potentially screws up restarts
            #  (if the threshold is changed after meryl is finished,
            #  for example).  It's only useful on large assemblies,
            #  which we usually assume you know what's going on
            #  anyway.
            #
            #  N.B. the mer overlapper NEEDS all mer counts 2 and
            #  higher.

            $cmd  = "$bin/meryl ";
            $cmd .= " -B $merCanonical -v -m $merSize $merylMemory -threads $merylThreads -c $merComp ";
            $cmd .= " -L 2 ";
            $cmd .= " -s $wrk/$asm.gkpStore:chain ";  #  Or 'latest' to get the original version
            $cmd .= " -o $ofile ";
            $cmd .= "> $wrk/0-mercounts/meryl.err 2>&1";

            stopBefore("meryl", $cmd);

            if (runCommand("$wrk/0-mercounts", $cmd)) {
                caFailure("meryl failed", "$wrk/0-mercounts/meryl.err");
            }
            unlink "$wrk/0-mercounts/meryl.err";
        }

        if ($merThresh eq "auto") {
            if (! -e "$ofile.estMerThresh.out") {
                $cmd  = "$bin/estimate-mer-threshold ";
                $cmd .= " -g $wrk/$asm.gkpStore:chain ";  #  Or 'latest' to get the original version
                $cmd .= " -m $ofile ";
                $cmd .= " > $ofile.estMerThresh.out ";
                $cmd .= "2> $ofile.estMerThresh.err";

                stopBefore("meryl", $cmd);

                if (runCommand("$wrk/0-mercounts", $cmd)) {
                    rename "$ofile.estMerThresh.out", "$ofile.estMerThresh.out.FAILED";
                    caFailure("estimate-mer-threshold failed", "$ofile.estMerThresh.err");
                }
            }

            open(F, "< $ofile.estMerThresh.out") or caFailure("failed to read estimated mer threshold from '$ofile.estMerThresh.out'", undef);
            $merThresh = <F>;
            $merThresh = int($merThresh * $merScale);
            close(F);

            if ($merThresh == 0) {
                caFailure("failed to estimate a mer threshold", "$ofile.estMerThresh.err");
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
                $cmd .= "2> $ffile.err ";

                if (runCommand("$wrk/0-mercounts", $cmd)) {
                    unlink $ffile;
                    caFailure("meryl failed to dump frequent mers", "$ffile.err");
                }
                unlink "$ffile.err";
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

        if ($merThresh eq "auto") {
            print STDERR "WARNING!  auto picking a mer threshold not supported without installing kmer\n";
            print STDERR "          (http://sourceforge.net/projects/kmer/).\n";
            print STDERR "Using historical defaults.\n";

            if ($merType eq "obt") {
                $merThresh = 1000;
            } else {
                $merThresh = 500;
            }
        }

        if (! -e $ofile) {
            my $mt = $merThresh / $merSkip;

            $cmd  = "$bin/meryl ";
            $cmd .= "-s $wrk/$asm.gkpStore -m $merSize -n $mt -K $merSkip ";
            $cmd .= " -o $ofile";
            $cmd .= "> $wrk/0-mercounts/meryl.err 2>&1";

            stopBefore("meryl", $cmd);

            if (runCommand("$wrk/0-mercounts", $cmd)) {
                unlink $ofile;
                caFailure("meryl failed to dump frequent mers", "$wrk/0-mercounts/meryl.err");
            }
            unlink "$wrk/0-mercounts/meryl.err";
        }

        symlink($ofile, $ffile) if (! -e $ffile);
    } else {
        caFailure("unknown meryl version '" . merylVersion() . "'", "");
    }

    return($merThresh);
}

sub meryl {
    system("mkdir $wrk/0-mercounts") if (! -d "$wrk/0-mercounts");

    if (getGlobal("ovlOverlapper") eq "umd") {
        caFailure("meryl attempted to compute mer counts for the umd overlapper", undef);
    }

    my $ovlc = 0;  #  No compression, unless we're the mer overlapper
    my $obtc = 0;

    my $ovlC = "-C";  #  Canonical, unless we're the mer overlapper
    my $obtC = "-C";  #  (except the mer overlapper now wants canonical)

    my $ovlD = 1;  #  Dump, unless we're the mer overlapper
    my $obtD = 1;

    my $obtT = 0;  #  New threshold
    my $ovlT = 0;

    #  If the mer overlapper, we don't care about single-copy mers,
    #  only mers that occur in two or more frags (kind of important
    #  for large assemblies).

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

    $ovlT = runMeryl(getGlobal('ovlMerSize'), $ovlc, $ovlC, getGlobal("ovlMerThreshold"), "ovl", $ovlD);
    $obtT = runMeryl(getGlobal('obtMerSize'), $obtc, $obtC, getGlobal("obtMerThreshold"), "obt", $obtD) if (getGlobal("doOverlapBasedTrimming") || getGlobal("doMerBasedTrimming"));

    if ((getGlobal("obtMerThreshold") ne $obtT) && (getGlobal("doOverlapBasedTrimming"))) {
        print STDERR "Reset OBT mer threshold from ", getGlobal("obtMerThreshold"), " to $obtT.\n";
        setGlobal("obtMerThreshold", $obtT);
    }
    
    if (getGlobal("ovlMerThreshold") ne $ovlT) {
        print STDERR "Reset OVL mer threshold from ", getGlobal("ovlMerThreshold"), " to $ovlT.\n";
        setGlobal("ovlMerThreshold", $ovlT);
    }
}

1;
