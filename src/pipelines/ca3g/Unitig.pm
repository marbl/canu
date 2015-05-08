package ca3g::Unitig;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(unitig);

use strict;

use File::Path qw(make_path remove_tree);

use ca3g::Defaults;
use ca3g::Execution;
use ca3g::Gatekeeper;


sub unitigger ($$$) {
    my $wrk  = shift @_;
    my $asm  = shift @_;
    my $per  = shift @_;
    my $bin  = getBinDirectory();
    my $cmd;

    my $e = " -e" .   getGlobal("utgErrorRate");
    my $u = " -U" .   getGlobal("utgBubblePopping");
    my $k = " -k" if (getGlobal("utgRecalibrateGAR") == 1);

    $cmd  = "$bin/unitigger \\\n";
    $cmd .= " -I $wrk/$asm.ovlStore \\\n";
    $cmd .= " -F $wrk/$asm.gkpStore \\\n";
    $cmd .= " -T $wrk/$asm.tigStore \\\n";
    $cmd .= " -o $wrk/4-unitigger/$asm \\\n";
    $cmd .= " -B $per -d 1 -x 1 -z 10 -j 5 $e$u$k \\\n";  #  $k comes with a space if it is needed
    $cmd .= " > $wrk/4-unitigger/unitigger.err 2>&1";

    return($cmd);
}


sub bogart ($$$) {
    my $wrk  = shift @_;
    my $asm  = shift @_;
    my $per  = shift @_;
    my $bin  = getBinDirectory();
    my $cmd;

    my $eg = getGlobal("utgGraphErrorRate");
    my $em = getGlobal("utgMergeErrorRate");

    my $t  = " -threads " . getGlobal("batThreads")  if defined(getGlobal("batThreads"));
    my $m  = " -M "       . getGlobal("batMemory")   if defined(getGlobal("batMemory"));
    my $o  = " "          . getGlobal("batOptions")  if defined(getGlobal("batOptions"));

    $cmd  = "$bin/bogart \\\n";
    $cmd .= " -O $wrk/$asm.ovlStore \\\n";
    $cmd .= " -G $wrk/$asm.gkpStore \\\n";
    $cmd .= " -T $wrk/$asm.tigStore \\\n";
    $cmd .= " -o $wrk/4-unitigger/$asm \\\n";
    $cmd .= " -B $per -eg $eg -em $em$t$m$o \\\n";
    $cmd .= " > $wrk/4-unitigger/unitigger.err 2>&1";

    return($cmd);
}


sub unitig ($$) {
    my $wrk    = shift @_;
    my $asm    = shift @_;

    goto stopAfter  if (skipStage($wrk, $asm, "unitig") == 1);
    goto allDone    if (-d "$wrk/$asm.tigStore");

    make_path("$wrk/4-unitigger")  if (! -d "$wrk/4-unitigger");

    #  How many reads per partition?  This will change - it'll move to be after unitigs are constructed.

    my $perPart = int(getNumberOfReadsInStore($wrk, $asm) / getGlobal("cnsPartitions"));
    my $minPart = getGlobal("cnsPartitionMin");

    $perPart = ($perPart < $minPart) ? ($perPart) : ($minPart);

    my $cmd;

    if      (getGlobal("unitigger") eq "bogart") {
        $cmd = bogart($wrk, $asm, $perPart);

    } elsif (getGlobal("unitigger") eq "unitigger") {
        $cmd = unitigger($wrk, $asm, $perPart);

    } else {
        caFailure("unknown unitigger '" . getGlobal("unitigger") . "'", undef);
    }

    stopBefore("unitig", $cmd);

    if (runCommand("$wrk/4-unitigger", $cmd)) {
        caExit("failed to unitig", "$wrk/4-unitigger/unitigger.err");
    }

  allDone:
    emitStage($wrk, $asm, "unitig");
  stopAfter:
    stopAfter("unitig");
}
