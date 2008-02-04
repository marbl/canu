use File::Copy;
use FindBin qw($Bin);

my ($username) = getpwuid($>); 

sub localDefaults () {

    setGlobal("sge", ' -P 08010 -b n -l msc');
    setGlobal("sgeOverlap", ' -pe threaded 2');
    setGlobal("useGrid", 1);
    setGlobal("scriptOnGrid", 1);
    setGlobal("cnsPartitions", 1);  # since it's microbes, no need to partition
}



sub localOption($@) {
    my $arg = shift;
    return (0, @_);
}



sub localSetup($) {
    my $numSteps = shift @_;

    #catmap
    if ( -e "$asm.catmap" and !-e "$wrk/$asm.catmap" ) {
      copy("$asm.catmap", "$wrk/$asm.catmap") or die "Could not copy: $asm.catmap\n";
    }

    #seq.features
    if ( -e "$asm.seq.features" and !-e "$wrk/$asm.seq.features" ) {
      copy("$asm.seq.features", "$wrk/$asm.seq.features") or die "Could not copy: $asm.seq.features\n";
    }

    my $vi = getGlobal("vectorIntersect");
    my $clvFile = "in.clv";
    if( !defined($vi) and !-e "$wrk/$clvFile" ) {
        if ( !-e "$asm.seq.features" ) {
            die("ERROR: Unable to create vector intersect file (.clv). Please provide the $asm.seq.features file.\n")
        }
        #create .clv file
        my $clv_cmd = "awk '{print \$1,\$5,\$6}' $asm.seq.features > $wrk/$clvFile";
        system($clv_cmd);
        setGlobal("vectorIntersect", $clvFile);
    }
}



sub localStart ($) {
    my $cmd_name = shift;
}



sub localFinish ($) {
    my $cmd_name = shift;
}



sub localFailure ($) {
    my $msg        = shift @_;

    open(MAIL,"| mail -s 'Assembly Failed' $username");
	print MAIL "$msg\n";
    close MAIL;
}


sub localPostTerminator($) {
    my $cmd_name = shift;
}

sub localFinalize() {

    # now send out a notice e-mail
    system("echo $wrk | mail -s 'Assembly completed' $username");
}


1;
