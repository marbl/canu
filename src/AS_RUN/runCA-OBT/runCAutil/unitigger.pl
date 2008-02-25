use strict;

sub unitigger (@) {
    my @cgbFiles  = @_;

    goto alldone if (scalar(@cgbFiles) > 0);

    my $bin = getBinDirectory();

    if (! -e "$wrk/4-unitigger/unitigger.success") {
        system("mkdir $wrk/4-unitigger") if (! -e "$wrk/4-unitigger");

        my $l = getGlobal("utgGenomeSize");
        my $e = getGlobal("utgErrorRate");

        my $B = int($numFrags / getGlobal("cnsPartitions"));
        $B = getGlobal("cnsMinFrags") if ($B < getGlobal("cnsMinFrags"));

        my $unitigger = getGlobal("unitigger");

        my $cmd;

        if ($unitigger eq "bog") {
            my $bmd = getGlobal("bogBadMateDepth");

            $cmd  = "$bin/buildUnitigs ";
            $cmd .= " -O $wrk/$asm.ovlStore ";
            $cmd .= " -G $wrk/$asm.gkpStore ";
            $cmd .= " -B $B ";
            $cmd .= " -e $e ";
            $cmd .= " -s $l "   if (defined($l));
            $cmd .= " -b "      if (getGlobal("bogPromiscuous") == 0);
            $cmd .= " -k "      if (getGlobal("bogEjectUnhappyContain") == 1);
            $cmd .= " -m $bmd " if (defined($bmd));
            $cmd .= " -o $wrk/4-unitigger/$asm ";
            $cmd .= " > $wrk/4-unitigger/unitigger.out 2>$wrk/4-unitigger/unitigger.err";
        } elsif ($unitigger eq "utg") {
            my $m = getGlobal("utgEdges");
            my $n = getGlobal("utgFragments");
            my $u = getGlobal("utgBubblePopping");

            $cmd  = "$bin/unitigger ";
            $cmd .= " -k " if (getGlobal("utgRecalibrateGAR") == 1);
            $cmd .= " -B $B ";
            $cmd .= " -l $l " if defined($l);
            $cmd .= " -m $m " if defined($m);
            $cmd .= " -n $n " if defined($n);
            $cmd .= " -d 1 -x 1 -z 10 -j 5 -U $u ";
            $cmd .= " -e $e ";
            $cmd .= " -F $wrk/$asm.gkpStore ";
            $cmd .= " -o $wrk/4-unitigger/$asm ";
            $cmd .= " -I $wrk/$asm.ovlStore ";
            $cmd .= " > $wrk/4-unitigger/unitigger.err 2>&1";
        } else {
            caFailure("Unknown unitigger $unitigger.\n");
        }

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
