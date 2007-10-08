use strict;

sub unitigger (@) {
    my @cgbFiles  = @_;

    goto alldone if (scalar(@cgbFiles) > 0);

    if ($global{'useBogUnitig'}) {
        bogUnitigger();
        goto alldone;
    }

    if (! -e "$wrk/4-unitigger/unitigger.success") {
        system("mkdir $wrk/4-unitigger") if (! -e "$wrk/4-unitigger");

        my $l = getGlobal("utgGenomeSize");
        my $m = getGlobal("utgEdges");
        my $e = getGlobal("utgErrorRate");
        my $n = getGlobal("utgFragments");
        my $u = getGlobal("utgBubblePopping");
        my $B = int($numFrags / getGlobal("cnsPartitions"));
        $B = getGlobal("cnsMinFrags") if ($B < getGlobal("cnsMinFrags"));

        my $cmd;
        $cmd  = "$bin/unitigger ";
	$cmd .= " -k " if (getGlobal("utgRecalibrateGAR") == 1);
        $cmd .= " -B $B ";
        $cmd .= " -l $l " if defined($l);
        $cmd .= " -m $m " if defined($m);
        $cmd .= " -n $n " if defined($n);
        $cmd .= " -d 1 -x 1 -z 10 -j 5 -U $u -e $e ";
        $cmd .= " -F $wrk/$asm.gkpStore ";
        $cmd .= " -o $wrk/4-unitigger/$asm ";
        $cmd .= " -I $wrk/$asm.ovlStore ";
        $cmd .= " > $wrk/4-unitigger/unitigger.err 2>&1";

        if (runCommand("$wrk/4-unitigger", $cmd)) {
            caFailure("Failed to unitig.\n");
        }

        touch("$wrk/4-unitigger/unitigger.success");
    }

  alldone:
    #  Other steps (consensus) need the list of cgb files, so we just do it here.
    #
    open(F, "ls $wrk/4-unitigger/*.cgb |") or caFailure("Failed to ls '$wrk/4-unitigger/*.cgb'\n");
    @cgbFiles = <F>;
    close(F);
    chomp @cgbFiles;

    stopAfter("unitigger");
    return @cgbFiles;
}

1;
