package ca3g::Gatekeeper;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(getNumberOfReadsInStore getNumberOfBasesInStore getExpectedCoverage gatekeeper);

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



sub getNumberOfBasesInStore ($$) {
    my $wrk = shift @_;
    my $asm = shift @_;
    my $nb  = 0;

    if (-e "$wrk/$asm.gkpStore/info.txt") {
        open(F, "< $wrk/$asm.gkpStore/info.txt") or caFailure("failed to open '$wrk/$asm.gkpStore/info.txt'", undef);
        while (<F>) {
            if (m/numBases\s+=\s+(\d+)/) {
                $nb = $1;
            }
        }
        close(F);
    }

    if ($nb == 0) {
        my $bin    = getBinDirectory();
        open(F, "$bin/gatekeeperDumpMetaData -G $wrk/$asm.gkpStore -stats |");
        while (<F>) {
            if (m/^library\s+0\s+reads\s+\d+\s+bases:\s+total\s+(\d+)\s+/) {
                $nb = $1;
            }
        }
        close(F);
    }

    return($nb);
}



sub getExpectedCoverage ($$) {
    my $wrk = shift @_;
    my $asm = shift @_;

    return(getNumberOfBasesInStore($wrk, $asm) / getGlobal("genomeSize"));
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

    caFailure("incomplete gatekeeper store $wrk/$asm.gkpStore; fix or remove", undef)
        if ((-d "$wrk/$asm.gkpStore") && (! -e "$wrk/$asm.gkpStore/info.txt"));

    caFailure("empty gatekeeper store $wrk/$asm.gkpStore", undef)
        if (($numReads == 0) && (-e "$wrk/$asm.gkpStore/info.txt"));

    #  Make sure all the inputs are here.

    my $failedFiles = undef;

    foreach my $iii (@inputs) {
        my $file = $iii;  #  This stupid foreach works by reference!

        $file = $2  if ($file =~ m/^(.*):(.*)/);   #  Handle the raw sequence inputs.

        if (! -e $file) {
            if (defined($failedFiles)) {
                $failedFiles .= "; '$file' not found";
            } else {
                $failedFiles = "'$file' not found";
            }
        }
    }
    caFailure($failedFiles, undef) if defined($failedFiles);

    #  Build a gkp file for all the raw sequence inputs.  For simplicity, we just copy in any gkp
    #  files as is.  This documents what gatekeeper was built with, etc.

    open(F, "> $wrk/$asm.gkpStore.gkp") or caFailure("failed to open '$wrk/$asm.gkpStore.gkp': $0", undef);

    foreach my $iii (@inputs) {
        if ($iii =~ m/^-(.*):(.*)$/) {
            my $tech = $1;
            my $file = $2;
            my @name = split '/', $2;
            my $name = $name[scalar(@name)-1];

            print F "########################################\n";
            print F "#  $tech: $file\n";
            print F "#\n";
            print F "name   $name\n";
            print F "preset $tech\n";
            print F "$file\n";
            print F "\n";

        } elsif (-e $iii) {
            print F "########################################\n";
            print F "#  $iii\n";
            print F "#\n";
            open(I, "< $iii") or caFailure("failed to open gatekeeper input '$iii': $0", undef);
            while (<I>) {
                print F $_;
            }
            close(I);
            print F "\n";

        } else {
            caFailure("don't recognize gatekeeper input '$iii'", undef);
        }
    }

    close(F);

    #  Load the store.

    my $cmd = "$bin/gatekeeperCreate -o $wrk/$asm.gkpStore.BUILDING $wrk/$asm.gkpStore.gkp > $wrk/$asm.gkpStore.err 2>&1";

    if (runCommand($wrk, $cmd)) {
        caFailure("gatekeeper failed", "$wrk/$asm.gkpStore.err");
    }

    rename "$wrk/$asm.gkpStore.BUILDING",             "$wrk/$asm.gkpStore";
    rename "$wrk/$asm.gkpStore.BUILDING.errorLog",    "$wrk/$asm.gkpStore.errorLog";

 stopafter:
    stopAfter("gatekeeper");
}

1;

