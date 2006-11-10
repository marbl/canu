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

    my $workDir = "$wrk/4-unitigger";
    mkdir $workDir unless -e $workDir;

    if (! -e "$workDir/unitigger.success") {

        my $cmd;
        $cmd  = "$bin/buildUnitigs ";
        $cmd .= " $wrk/$asm.ovlStore ";
        $cmd .= " $wrk/$asm.frgStore ";

        my $l = getGlobal("utgGenomeSize");
        $l ||= 0;
        $cmd .= " $l ";

        $cmd .= " > unitigger.out ";
        $cmd .= " 2> unitigger.err ";

        if (runCommand($workDir, $cmd)) {
            print STDERR "Failed to unitig.\n";
            exit(1);
        }

        my $prevPwd = $ENV{PWD};
        chdir $workDir || die "chdir $workDir failed.";
        link('len15.ium',"$asm.cgb") || die "link to $asm.cgb failed in $ENV{PWD}";
        partitionCGB( "$asm.cgb", 250000 );
        unlink "$asm.cgb";
        chdir $prevPwd;

        touch('unitigger.success');
    }
}

1;
