use strict;

sub unitigger {

    system("mkdir $wrk/4-unitigger") if (! -e "$wrk/4-unitigger");

    if (! -e "$wrk/4-unitigger/unitigger.success") {

        #  Unitigger needs an input .ofg file, which is pretty much just a
        #  dump of the fragstore.  We used to save the original, for no
        #  real reason.  Now, just rebuild the .ofg if the fragstore is
        #  newer than the ofg we have.
        #
        #  -M says "start time - mod time" so if that is small, it's newer
        #
        if ((! -e "$wrk/4-unitigger/$asm.ofg") ||
            (-M "$wrk/$asm.frgStore/db.frg" < -M "$wrk/4-unitigger/$asm.ofg")) {
            my $cmd;
            #$cmd  = "[ $wrk/$asm.frgStore/db.frg -nt $wrk/0-preoverlap/$asm.ofg ] && ";
            #$cmd .= "$bin/make_OFG_from_FragStore $wrk/$asm.frgStore > $wrk/0-preoverlap/$asm.ofg || ";
            #$cmd .= "true";
            $cmd .= "$bin/make_OFG_from_FragStore $wrk/$asm.frgStore > $wrk/4-unitigger/$asm.ofg";
            if (runCommand($cmd)) {
                rename "$wrk/4-unitigger/$asm.ofg", "$wrk/4-unitigger/$asm.ofg.FAILED";
                die "Failed.\n";
            }
        }

        if (! -e "$wrk/4-unitigger/$asm.ofgList") {
            if (runCommand("ls -1 $wrk/4-unitigger/*.ofg > $wrk/4-unitigger/$asm.ofgList")) {
                print STDERR "Failed to find the ofg's.\n";
                rename "$wrk/4-unitigger/$asm.ofgList", "$wrk/4-unitigger/$asm.ofgList.FAILED";
                exit(1);
            }
        }

        my $cmd;
        $cmd  = "cd $wrk/4-unitigger && ";
        $cmd .= "$bin/unitigger ";
        $cmd .= " -B 250000 ";

        my $l = getGlobal("utgGenomeSize");
        my $m = getGlobal("utgEdges");
        my $e = getGlobal("utgErrorRate");
        my $n = getGlobal("utgFragments");

        $cmd .= " -l $l " if defined($l);
        $cmd .= " -m $m " if defined($m);
        $cmd .= " -n $n " if defined($n);

        $cmd .= " -c -P -A 1 -d 1 -x 1 -z 10 -j 5 -U 1 -e $e ";
        $cmd .= " -F $wrk/$asm.frgStore ";
        $cmd .= " -f ";
        $cmd .= " -o $wrk/4-unitigger/$asm.fgbStore ";
        $cmd .= " -L $wrk/4-unitigger/$asm.ofgList ";
        $cmd .= " -I $wrk/$asm.ovlStore ";
        $cmd .= " > $wrk/4-unitigger/unitigger.out ";
        $cmd .= " 2> $wrk/4-unitigger/unitigger.err ";

        if (runCommand($cmd)) {
            print STDERR "Failed to unitig.\n";
            exit(1);
        }

        touch("$wrk/4-unitigger/unitigger.success");
    }
}

1;
