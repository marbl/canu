use strict;

sub runMeryl ($$$) {
    my $merSize   = shift @_;
    my $merThresh = shift @_;
    my $merType   = shift @_;
    my $cmd;

    if (getGlobal("doMeryl") == 0) {
        touch "$wrk/0-mercounts/$asm.nmers.$merType.fasta";
        return;
    }

    if (merylVersion() eq "Mighty") {

        #  Use the better meryl!  This is straightforward.  We count,
        #  then we dump.

        if (! -e "$wrk/0-mercounts/$asm-ms$merSize.mcdat") {
            my $merylMemory = getGlobal("merylMemory");

            $cmd .= "$bin/meryl ";
            $cmd .= " -B -C -v -m $merSize -memory $merylMemory ";
            $cmd .= " -s $wrk/$asm.gkpStore ";
            $cmd .= " -o $wrk/0-mercounts/$asm-ms$merSize ";
            $cmd .= "> $wrk/0-mercounts/meryl.out 2>&1";

            if (runCommand("$wrk/0-mercounts", $cmd)) {
                (print "Failed.\n" && return -1);
            }
        }

        if (! -e "$wrk/0-mercounts/$asm.nmers.$merType.fasta") {
            $cmd  = "$bin/meryl ";
            $cmd .= "-Dt -n $merThresh ";
            $cmd .= "-s $wrk/0-mercounts/$asm-ms$merSize ";
            $cmd .= "> $wrk/0-mercounts/$asm.nmers.$merType.fasta ";

            if (runCommand("$wrk/0-mercounts", $cmd)) {
                unlink "$wrk/0-mercounts/$asm.nmers.$merType.fasta";
                (print "Failed.\n" && return -1);
            }
        }
    } elsif (merylVersion() eq "CA") {

        #  Sigh.  The old meryl.  Not as easy.  If we assume the
        #  process, in particular that the Ovl threshold is less than
        #  the Obt threshold, and that we have already computed the
        #  Ovl mers, we could filter the Ovl mers to get the Obt mers.

        if (! -e "$wrk/0-mercounts/$asm.nmers.$merType.fasta") {
            my $merSkip = 10;

            $merThresh /= $merSkip;

            $cmd  = "$bin/meryl ";
            $cmd .= "-s $wrk/$asm.gkpStore -m $merSize -n $merThresh -K $merSkip ";
            $cmd .= " -o $wrk/0-mercounts/$asm.nmers.$merType.fasta";
            $cmd .= "> $wrk/0-mercounts/meryl.out 2>&1";

            if (runCommand("$wrk/0-mercounts", $cmd)) {
                unlink "$wrk/0-mercounts/$asm.nmers.$merType.fasta";
                (print "Failed.\n" && return -1);
            }
        }
    } else {
        (print "Unknown meryl version.\n" && return -1);
    }
}

sub meryl {
    system("mkdir $wrk/0-mercounts") if (! -d "$wrk/0-mercounts");

    runMeryl(getGlobal('merSizeOvl'), getGlobal("merylOvlThreshold"), "ovl");
    runMeryl(getGlobal('merSizeObt'), getGlobal("merylObtThreshold"), "obt");
}

1;
