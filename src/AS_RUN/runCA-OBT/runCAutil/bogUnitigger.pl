use strict;

sub partitionCGB($$) {
    my ($cgbFile, $targetFrags) = @_;
    local $/ = "\n{IUM\n";
    my $numFrags;
    my $fileNum=1;
    my $firstIUM = 1;
    my $name = sprintf "%s_%03d.cgb",$asm,$fileNum;
    open(IN,"<$cgbFile") || (print "Can't read $cgbFile" && return -1);
    open(OUT,">$name") || (print "Can't write $name" && return -1);
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
            open(OUT,">$name") || (print "Can't write $name" && return -1);
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
        $cmd .= " -O $wrk/$asm.ovlStore ";
        $cmd .= " -G $wrk/$asm.gkpStore ";

        my $l = getGlobal("utgGenomeSize");
        $l ||= 0;
        $cmd .= " -s $l ";

        my $l = getGlobal("bogPromiscuous");
        $cmd .= " -b " if !$l;

        $cmd .= " > unitigger.out ";
        $cmd .= " 2> unitigger.err ";

        if (runCommand($workDir, $cmd)) {
            print STDERR "Failed to unitig.\n";
            caFailure();
        }

        my $prevPwd = $ENV{PWD};
        chdir $workDir || (print "chdir $workDir failed." && return -1);
        link('len150.ium',"$asm.cgb") || (print "link to $asm.cgb failed in $ENV{PWD}" && return -1);
        partitionCGB( "$asm.cgb", 250000 );
        unlink "$asm.cgb";
        chdir $prevPwd;

        touch("$workDir/unitigger.success");
    }
}

1;
