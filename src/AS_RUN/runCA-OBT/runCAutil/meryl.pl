use strict;

#  This needs the fragstore, but we could use the parallel meryl on
#  *.inp, the output from gatekeeper.

sub meryl {
    my $merSize = getGlobal('merSize');


    #  If we start from prepackaged frg and gkp stores, we don't have
    #  this directory.
    #
    system("mkdir $wrk/0-preoverlap") if (! -d "$wrk/0-preoverlap");

    #  What threshold to use?
    #
    my $merylObtThreshold = getGlobal("merylObtThreshold");
    my $merylOvlThreshold = getGlobal("merylOvlThreshold");

    #  Use the better meryl, .....
    #
    if (merylVersion() eq "Mighty") {
        if (! -e "$wrk/0-preoverlap/$asm.nmers.fasta") {
            my $merylMemory = getGlobal("merylMemory");

            my $cmd;
            $cmd .= "$bin/meryl ";
            $cmd .= " -B -C -v -m $merSize -memory $merylMemory ";
            $cmd .= " -s $wrk/$asm.gkpStore ";
            $cmd .= " -o $wrk/0-preoverlap/$asm-ms$merSize ";
            $cmd .= "> $wrk/0-preoverlap/meryl.out 2>&1";

            if (runCommand("$wrk/0-preoverlap", $cmd)) {
                die "Failed.\n";
            }

            $cmd = "$bin/meryl -Dt -n $merylOvlThreshold -s $wrk/0-preoverlap/$asm-ms$merSize > $wrk/0-preoverlap/$asm.nmers.fasta ";
            if (runCommand("$wrk/0-preoverlap", $cmd)) {
                rename "$wrk/$asm.nmers.fasta", "$wrk/$asm.nmers.fasta.FAILED";
                die "Failed.\n";
            }
        }

        if ((! -e "$wrk/0-overlaptrim-overlap/$asm.nmers.fasta") && (getGlobal("doOverlapTrimming") != 0)) {
            my $cmd;
            $cmd = "$bin/meryl -Dt -n $merylObtThreshold -s $wrk/0-preoverlap/$asm-ms$merSize > $wrk/0-overlaptrim-overlap/$asm.nmers.fasta ";
            
            if (runCommand("$wrk/0-overlaptrim-overlap", $cmd)) {
                rename "$wrk/$asm.nmers.fasta", "$wrk/$asm.nmers.fasta.FAILED";
                die "Failed.\n";
            }
        }
    } elsif (merylVersion() eq "CA") {
        if (! -e "$wrk/0-preoverlap/$asm.nmers.fasta") {

            #  Meryl is run at the most resticted setting (the one used
            #  for assembly); Overlap trimming is responsible for taking
            #  the $asm.nmers.fasta and extracting mers that it cares
            #  about (typically 2x the threshold here).

            my $merylSkip = 10;

            $merylObtThreshold /= $merylSkip;
            $merylOvlThreshold /= $merylSkip;

            my $cmd;
            $cmd  = "$bin/meryl ";
            $cmd .= " -m $merSize -s $wrk/$asm.gkpStore -n $merylOvlThreshold -K $merylSkip ";
            $cmd .= " -o $wrk/0-preoverlap/$asm.nmers.fasta";
            $cmd .= "> $wrk/0-preoverlap/meryl.out 2>&1";

            if (runCommand("$wrk/0-preoverlap", $cmd)) {
                rename "$wrk/$asm.nmers.fasta", "$wrk/$asm.nmers.fasta.FAILED";
                die "Failed.\n";
            }
        }

        if ((! -e "$wrk/0-overlaptrim-overlap/$asm.nmers.fasta") && (getGlobal("doOverlapTrimming") != 0)) {
            open(F, "< $wrk/0-preoverlap/$asm.nmers.fasta")  or die "Failed to open $wrk/0-preoverlap/$asm.nmers.fasta for reading.\n";
            open(G, "> $wrk/0-overlaptrim-overlap/$asm.nmers.fasta") or die "Failed to open $wrk/0-overlaptrim-overlap/$asm.nmers.fasta for writing.\n";
            while (!eof(F)) {
                my $def = <F>;
                my $mer = <F>;
                if ($def =~ m/^>(\d+)$/) {
                    print G "$def$mer" if ($1 > $merylObtThreshold);
                } else {
                    chomp $def;
                    print STDERR "ERROR:  Got '$def' for a defline!\n";
                }
            }
            close(G);
            close(F);
        }
    } else {
        die "Unknown meryl version.\n";
    }
}


1;
