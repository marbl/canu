use strict;

sub partitionCGB($$) {
    my ($cgbFile, $targetFrags) = @_;
    local $/ = "\n{IUM\n";
    my $numFrags;
    my $fileNum=1;
    my $firstIUM = 1;
    my $name = sprintf "%s_%03d.cgb",$asm,$fileNum;
    open(IN,"<$cgbFile") || die "Can't read $cgbFile";
    open(OUT,">$name") || die "Can't write $name";
    while(<IN>) {
        chomp;
        my $lastIMP = 0;
        while ( -1 != ($lastIMP = index $_,"\n{IMP\n",$lastIMP)) {
            $numFrags++;
            $lastIMP++;
        }
        if ($numFrags >= $targetFrags) {
            $numFrags = 0;
            $fileNum++;
            my $name = sprintf "%s_%03d.cgb",$asm,$fileNum;
            open(OUT,">$name") || die "Can't write $name";
            print OUT "{IUM\n";
            $firstIUM = 1;
        }
        print OUT $/ unless $firstIUM;
        print OUT $_;
        $firstIUM = 0;
    }
    close IN;
    close OUT;
}

sub bogUnitigger {

    system("mkdir $wrk/4-unitigger") if (! -e "$wrk/4-unitigger");

    if (! -e "$wrk/4-unitigger/unitigger.success") {

        my $cmd;
        $cmd  = "$bin/buildUnitigs ";
        $cmd .= " $wrk/$asm.ovlStore ";
        $cmd .= " $wrk/$asm.frgStore ";

        my $l = getGlobal("utgGenomeSize");
        $l ||= 0;
        $cmd .= " $l ";

        $cmd .= " > unitigger.out ";
        $cmd .= " 2> unitigger.err ";

        if (runCommand("$wrk/4-unitigger", $cmd)) {
            print STDERR "Failed to unitig.\n";
            exit(1);
        }

        link 'len15.ium',"$asm.cgb";
        partitionCGB( "$asm.cgb", 250000 );
        unlink "$asm.cgb";

        touch('unitigger.success');
    }
}

1;
