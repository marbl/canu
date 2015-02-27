package ca3g::Gatekeeper;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(getNumberOfReadsInStore gatekeeper);

use strict;

use ca3g::Defaults;
use ca3g::Execution;




sub getNumberOfReadsInStore ($$) {
    my $wrk = shift @_;
    my $asm = shift @_;
    my $nr  = 0;

    if (-e "$wrk/$asm.gkpStore/info.txt") {
        open(F, "< $wrk/$asm.gkpStore/info.txt") or caFailure("failed to open '$wrk/$asm.gkpStore/info.txt'", undef);
        while (<F>) {
            if (m/numReads\s+=\s+(\d+)/) {
                $nr = $1;
            }
        }
        close(F);
    }

    return($nr);
}



sub gatekeeper ($$@) {
    my $wrk    = shift @_;
    my $asm    = shift @_;
    my $bin    = getBinDirectory();
    my @inputs = @_;

    my $numReads = getNumberOfReadsInStore($wrk, $asm);

    #  Return if there are fragments in the store, and die if there
    #  are no fragments and no source files.

    if ($numReads > 0) {
        stopAfter("gatekeeper");
        return;
    }

    #  Check for inputs, or an incomplete or empty store.

    caFailure("no input files specified, and store not already created", undef)
        if (scalar(@inputs) == 0);

    caFailure("incomplete gatekeeper store $wrk/$asm.gkpStore; fix or remove\n", undef)
        if ((-d "$wrk/$asm.gkpStore") && (! -e "$wrk/$asm.gkpStore/info.txt"));

    caFailure("empty gatekeeper store $wrk/$asm.gkpStore?\n", undef)
        if (($numReads == 0) && (-e "$wrk/$asm.gkpStore/info.txt"));

    #  Make sure all the inputs are here.

    my $failedFiles = undef;
    my $gkpInput = "";
    foreach my $frg (@inputs) {
        if (! -e $frg) {
            if (defined($failedFiles)) {
                $failedFiles .= "; '$frg' not found";
            } else {
                $failedFiles = "'$frg' not found";
            }
        }

        $gkpInput .= " $frg";
    }
    caFailure($failedFiles, undef) if defined($failedFiles);

    #  Load the store.

    my $cmd = "$bin/gatekeeperCreate -o $wrk/$asm.gkpStore.BUILDING $gkpInput > $wrk/$asm.gkpStore.err 2>&1";

    if (runCommand($wrk, $cmd)) {
        caFailure("gatekeeper failed", "$wrk/$asm.gkpStore.err");
    }

    rename "$wrk/$asm.gkpStore.BUILDING",             "$wrk/$asm.gkpStore";
    rename "$wrk/$asm.gkpStore.BUILDING.errorLog",    "$wrk/$asm.gkpStore.errorLog";

 stopafter:
    stopAfter("gatekeeper");
}

1;

