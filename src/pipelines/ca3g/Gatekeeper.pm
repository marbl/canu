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
        open(F, "< $wrk/$asm.gkpStore/info.txt") or caExit("can't open '$wrk/$asm.gkpStore/info.txt' for reading: $!", undef);
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
        open(F, "< $wrk/$asm.gkpStore/info.txt") or caExit("can't open '$wrk/$asm.gkpStore/info.txt' for reading: $!", undef);
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

    return(int(getNumberOfBasesInStore($wrk, $asm) / getGlobal("genomeSize")));
}



sub gatekeeper ($$$@) {
    my $WRK    = shift @_;  #  Root work directory
    my $wrk    = $WRK;
    my $asm    = shift @_;
    my $tag    = shift @_;
    my $bin    = getBinDirectory();
    my @inputs = @_;

    $wrk = "$wrk/correction"  if ($tag eq "cor");
    $wrk = "$wrk/trimming"    if ($tag eq "obt");

    my $numReads = getNumberOfReadsInStore($wrk, $asm);

    goto stopAfter  if (skipStage($WRK, $asm, "$tag-gatekeeper") == 1);    #  Finished.
    goto allDone    if ($numReads > 0);                                    #  Store exists, with reads.

    #  Check for inputs, or an incomplete or empty store.

    caExit("empty gatekeeper store $wrk/$asm.gkpStore", undef)
        if (($numReads == 0) && (-e "$wrk/$asm.gkpStore/info.txt"));

    caExit("incomplete gatekeeper store $wrk/$asm.gkpStore; fix or remove", undef)
        if ((-d "$wrk/$asm.gkpStore") && (! -e "$wrk/$asm.gkpStore/info.txt"));

    caExit("no input files specified, and store not already created", undef)
        if (scalar(@inputs) == 0);

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
    caExit($failedFiles, undef) if defined($failedFiles);

    #  Build a gkp file for all the raw sequence inputs.  For simplicity, we just copy in any gkp
    #  files as is.  This documents what gatekeeper was built with, etc.

    open(F, "> $wrk/$asm.gkpStore.gkp") or caExit("cant' open '$wrk/$asm.gkpStore.gkp' for writing: $0", undef);

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
            open(I, "< $iii") or caExit("can't open gatekeeper input '$iii' for reading: $0", undef);
            while (<I>) {
                print F $_;
            }
            close(I);
            print F "\n";

        } else {
            caExit("unrecognized gatekeeper input file '$iii'", undef);
        }
    }

    close(F);

    #  Load the store.

    my $cmd;
    $cmd .= "$bin/gatekeeperCreate \\\n";
    $cmd .= "  -minlength " . getGlobal("minReadLength") . " \\\n";
    $cmd .= "  -o $wrk/$asm.gkpStore.BUILDING \\\n";
    $cmd .= "  $wrk/$asm.gkpStore.gkp \\\n";
    $cmd .= "> $wrk/$asm.gkpStore.err 2>&1";

    if (runCommand($wrk, $cmd)) {
        caExit("gatekeeper failed", "$wrk/$asm.gkpStore.err");
    }

    rename "$wrk/$asm.gkpStore.BUILDING",             "$wrk/$asm.gkpStore";
    rename "$wrk/$asm.gkpStore.BUILDING.errorLog",    "$wrk/$asm.gkpStore.errorLog";

  allDone:
    emitStage($WRK, $asm, "$tag-gatekeeper");
 stopAfter:
    stopAfter("gatekeeper");
}
