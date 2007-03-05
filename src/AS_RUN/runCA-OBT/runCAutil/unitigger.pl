use strict;

sub unitigger (@) {
    my @cgbFiles = @_;

    goto alldone if (scalar(@cgbFiles) > 0);

    if ($global{'useBogUnitig'}) {
        bogUnitigger();
        goto alldone;
    }

    if (! -e "$wrk/4-unitigger/unitigger.success") {
        system("mkdir $wrk/4-unitigger") if (! -e "$wrk/4-unitigger");

        #  Unitigger needs an input .ofg file, which is pretty much just a
        #  dump of the fragstore.  We used to save the original, for no
        #  real reason.  Now, just rebuild the .ofg if the fragstore is
        #  newer than the ofg we have.
        #
        #  -M says "start time - mod time" so if that is small, it's newer
        #
        if ((! -e "$wrk/4-unitigger/$asm.ofg") ||
            (-M "$wrk/$asm.gkpStore/db.frg" < -M "$wrk/4-unitigger/$asm.ofg")) {
            my $cmd;
            $cmd  = "$bin/gatekeeper -ofg $wrk/$asm.gkpStore ";
            $cmd .= " > $wrk/4-unitigger/$asm.ofg ";
            $cmd .= " 2> $wrk/4-unitigger/$asm.ofg.err";

            if (runCommand("$wrk/4-unitigger", $cmd)) {
                rename "$wrk/4-unitigger/$asm.ofg", "$wrk/4-unitigger/$asm.ofg.FAILED";
                die "Failed.\n";
            }
        }

        if (! -e "$wrk/4-unitigger/$asm.ofgList") {
            open(G, "> $wrk/4-unitigger/$asm.ofgList") or die;
            print G "$wrk/4-unitigger/$asm.ofg\n";
            close(G);
        }

        my $l = getGlobal("utgGenomeSize");
        my $m = getGlobal("utgEdges");
        my $e = getGlobal("utgErrorRate");
        my $n = getGlobal("utgFragments");
        my $u = getGlobal("utgBubblePopping");
        my $B = int($numFrags / getGlobal("cnsPartitions"));
        $B = getGlobal("cnsMinFrags") if ($B < getGlobal("cnsMinFrags"));

        my $cmd;
        $cmd  = "$bin/unitigger ";
        $cmd .= " -B $B ";
        $cmd .= " -l $l " if defined($l);
        $cmd .= " -m $m " if defined($m);
        $cmd .= " -n $n " if defined($n);
        $cmd .= " -c -A 1 -d 1 -x 1 -z 10 -j 5 -U $u -e $e ";
        $cmd .= " -F $wrk/$asm.gkpStore ";
        $cmd .= " -f ";
        $cmd .= " -o $wrk/4-unitigger/$asm.fgbStore ";
        $cmd .= " -L $wrk/4-unitigger/$asm.ofg ";
        $cmd .= " -I $wrk/$asm.ovlStore ";
        $cmd .= " > $wrk/4-unitigger/unitigger.err 2>&1";

        if (runCommand("$wrk/4-unitigger", $cmd)) {
            print STDERR "Failed to unitig.\n";
            exit(1);
        }

        touch("$wrk/4-unitigger/unitigger.success");
    }

  alldone:
    #  Other steps (consensus) need the list of cgb files, so we just do it here.
    #
    open(F, "ls $wrk/4-unitigger/*.cgb |") or die;
    @cgbFiles = <F>;
    close(F);
    chomp @cgbFiles;

    stopAfter("unitigger");
    return(@cgbFiles);
}

1;
