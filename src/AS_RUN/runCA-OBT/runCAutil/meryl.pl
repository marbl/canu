use strict;

#  This needs the fragstore, but we could use the parallel meryl on
#  *.inp, the output from gatekeeper.

sub meryl {
    my $merSize = 22;

    #  If we want to use the better meryl, .....
    #
    if (0) {
        if (! -e "$wrk/0-preoverlap/$asm.nmers.fasta") {
            my $cmd;
            $cmd  = "cd $wrk/0-preoverlap && ";
            $cmd .= "$bin/dumpFragStoreAsFasta -frg $wrk/$asm.frgStore | ";
            $cmd .= "/home/work/src/genomics/meryl/meryl -B -C -v -m $merSize -s - -o $wrk/0-preoverlap/$asm -L 10 -memory 800 ";
            $cmd .= "> $wrk/0-preoverlap/meryl.out ";
            $cmd .= "2> $wrk/0-preoverlap/meryl.err";

            if (runCommand($cmd)) {
                print STDERR "Failed.\n";
                rename "$wrk/$asm.nmers.fasta", "$wrk/$asm.nmers.fasta.FAILED";
                exit(1);
            }

            $cmd  = "cd $wrk/0-preoverlap && ";
            $cmd .= "/home/work/src/genomics/meryl/meryl -Dt -n 500 -s $wrk/0-preoverlap/$asm > $wrk/0-preoverlap/$asm.nmers.fasta ";

            if (runCommand($cmd)) {
                print STDERR "Failed.\n";
                rename "$wrk/$asm.nmers.fasta", "$wrk/$asm.nmers.fasta.FAILED";
                exit(1);
            }


            #  We need a special case for making the overlap-trim nmers, since the thresholds changed with the new meryl

            if ((-d "$wrk/0-overlaptrim-overlap") && (! -e "$wrk/0-overlaptrim-overlap/$asm.nmers.fasta")) {
                my $cmd;
                $cmd  = "cd $wrk/0-overlaptrim-overlap && ";
                $cmd .= "/home/work/src/genomics/meryl/meryl -Dt -n 1000 -s $wrk/0-preoverlap/$asm > $wrk/0-overlaptrim-overlap/$asm.nmers.fasta ";

                if (runCommand($cmd)) {
                    print STDERR "Failed.\n";
                    rename "$wrk/$asm.nmers.fasta", "$wrk/$asm.nmers.fasta.FAILED";
                    exit(1);
                }
            }
        }
    }



    if (! -e "$wrk/0-preoverlap/$asm.nmers.fasta") {

        #  Meryl is run at the most resticted setting (the one used
        #  for assembly); Overlap trimming is responsible for taking
        #  the $asm.nmers.fasta and extracting mers that it cares
        #  about (typically 2x the threshold here).

        #  If you change these thresholds, you should also change
        #  overlapTrim.pl.

        my $cmd;
        $cmd  = "cd $wrk/0-preoverlap && ";
        $cmd .= "$bin/meryl -m $merSize -s $wrk/$asm.frgStore -n 50 -K 10 ";
        $cmd .= "-o $wrk/0-preoverlap/$asm.nmers.fasta";
        $cmd .= "> $wrk/0-preoverlap/meryl.out ";
        $cmd .= "2> $wrk/0-preoverlap/meryl.err";

        if (runCommand($cmd)) {
            print STDERR "Failed.\n";
            rename "$wrk/$asm.nmers.fasta", "$wrk/$asm.nmers.fasta.FAILED";
            exit(1);
        }
    }
}

1;
