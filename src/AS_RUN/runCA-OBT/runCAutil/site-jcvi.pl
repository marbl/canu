
use File::Copy;
use Mail::Mailer;
my ($username) = getpwuid($>) . '@jcvi.org'; 

sub localDefaults () {

    #  This sets some options JCVI likes to use for microbes.
    #
    #  These two used to be fixed, but the defaults work:
    #
    #  cnsPartitions     -- Partitions are no smaller than 75,000 fragments.  If you
    #                       have more than that, you'll want to pay the cost of a
    #                       partitioning so you can get more CPU.
    #  computeInsertSize -- Enabled if you have fewer than 1,000,000 fragments.

    setGlobal("sge", ' -P 08010 -b n -l msc');
    setGlobal("sgeOverlap", ' -pe threaded 2');
    setGlobal("useGrid", 1);
    setGlobal("scriptOnGrid", 1);
}



sub localOption($@) {
    my $arg = shift;

    #  If the user gave us a vectorIntersect option, DISABLE the
    #  requirement of seq.features below.  A complete abuse of
    #  quasi-private data.
    #
    if ($arg =~ m/^vectorIntersect/) {
        $global{vectorIntersectSupplied} = 1;
    }
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

    my $clvFile = "in.clv";

    #  If there is no vectorIntersect, and vectorIntersect= was not
    #  set (the noVec.specfile sets this to disable this block) and
    #  the clv file isn't already there, make one.

    if ( !defined(getGlobal("vectorIntersect")) and !defined($global{"vectorIntersectSupplied"}) and !-e "$wrk/$clvFile" ) {
        if ( !-e "$asm.seq.features" ) {
            die("ERROR: Unable to create vector intersect file (.clv). Please provide the $asm.seq.features file.\n")
        }
        #create .clv file
        my $clv_cmd = "awk '{print \$1,\$5,\$6}' $asm.seq.features > $wrk/$clvFile";
        system($clv_cmd);
        setGlobal("vectorIntersect", "$wrk/$clvFile");
        $commandLineOptions .= qq( "vectorIntersect=$wrk/$clvFile" );
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

    my $mailer = new Mail::Mailer;
    $mailer->open({ To => $username, Subject => 'Assembly Failed'});
	print $mailer "$msg\n";
    $mailer->close;
}


sub localPostTerminator($) {
    my $cmd_name = shift;
}

sub localFinalize() {

    # now send out a notice e-mail
    my $mailer = new Mail::Mailer;
    $mailer->open({ To => $username, Subject => 'Assembly completed'});
	print $mailer "$wrk\n";
    $mailer->close;
}


1;
