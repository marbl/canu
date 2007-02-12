use strict;


sub preoverlap {
    my @fragFiles = @_;

    $numFrags = getNumberOfFragsInStore($bin, $wrk, $asm);

    #  Return if there are fragments in the store, and die if there
    #  are no fragments and no source files.
    #
    goto alldone if ($numFrags > 0);

    die "ERROR: No fragment files specified, and stores not already created.\n" if (scalar(@fragFiles) == 0);

    system("mkdir $wrk/0-preoverlap") if (! -d "$wrk/0-preoverlap");

    if (! -e "$wrk/0-preoverlap/$asm.frg") {
        my $failedFiles = 0;
        foreach my $frg (@fragFiles) {
            if (! -e $frg) {
                print STDERR "MISSING: $frg\n";
                $failedFiles++;
            }
        }
        die if ($failedFiles);

        #  Rather than deal with grabbing the first BAT from the frist
        #  file (which is arbitrary anyway) we just write our own BAT.
        #
        #  The accession is chosen to be one before fakeUIDs used in
        #  updating distance records.
        #
        open(G, "> $wrk/0-preoverlap/$asm.frg") or die;
        print G "{BAT\n";
        print G "bna:$asm frags\n";
        print G "crt:0\n";
        print G "acc:987654321987654320\n";
        print G "com:";
        print G "All fragments contained in:\n";
        foreach my $frg (@fragFiles) {
            print G "  $frg\n";
        }
        print G ".\n";
        print G "}\n";
        close(G);

        #  Remove BAT and ADT messages from the input
        #    o BAT should only be at the start of the file
        #    o ADT (audit) should be OK to have multiple copies, but test
        #    o _major_ performance boost if one frag file, as many things
        #      will copy the store before adding

        foreach my $frg (@fragFiles) {
            print STDERR "$frg\n";

            my $cmd = "$bin/extractmessages -x ADT BAT < $frg >> $wrk/0-preoverlap/$asm.frg";

            if      ($frg =~ m/\.gz$/) {
                $cmd = "gzip -dc $frg | $bin/extractmessages -x ADT BAT >> $wrk/0-preoverlap/$asm.frg";
            } elsif ($frg =~ m/\.bz2$/) {
                $cmd = "bzip2 -dc $frg | $bin/extractmessages -x ADT BAT >> $wrk/0-preoverlap/$asm.frg";
            }

            if (runCommand("$wrk/0-preoverlap", "$cmd 2>> $wrk/0-preoverlap/extract.err")) {
                rename "$wrk/0-preoverlap/$asm.frg", "$wrk/0-preoverlap/$asm.frg.FAILED";
                die "Failed.\n";
            }
        }
    }

    if ((! -d "$wrk/$asm.gkpStore") || (! -e "$wrk/$asm.gkpStore/frg")) {
        my $cmd;
        $cmd  = "$bin/gatekeeper -e 10000000 ";
        $cmd .= "-Q -T " if (getGlobal("doOverlapTrimming"));
        $cmd .= "-o $wrk/$asm.gkpStore ";
        $cmd .= "$wrk/0-preoverlap/$asm.frg ";
        $cmd .= "> $wrk/0-preoverlap/gatekeeper.out ";
        $cmd .= "2> $wrk/0-preoverlap/gatekeeper.err";

        if (runCommand("$wrk/0-preoverlap", $cmd)) {
            print STDERR "Failed.\n";
            rename "$wrk/0-preoverlap/$asm.inp", "$wrk/0-preoverlap/$asm.inp.FAILED";
            rename "$wrk/$asm.gkpStore", "$wrk/$asm.gkpStore.FAILED";
            exit(1);
        }
    }

    $numFrags = getNumberOfFragsInStore($bin, $wrk, $asm);

alldone:
    stopAfter("initialStoreBuilding");
}

1;
