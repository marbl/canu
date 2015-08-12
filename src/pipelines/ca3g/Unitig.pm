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

    $cmd  = "$bin/bogart \\\n";
    $cmd .= " -G $wrk/$asm.gkpStore \\\n";
    $cmd .= " -O $wrk/$asm.ovlStore \\\n";
    $cmd .= " -T $wrk/$asm.tigStore \\\n";
    $cmd .= " -o $wrk/4-unitigger/$asm \\\n";
    $cmd .= " -B $per \\\n";
    $cmd .= " -eg "      . getGlobal("utgGraphErrorRate")  . " \\\n";
    $cmd .= " -eb "      . getGlobal("utgBubbleErrorRate") . " \\\n";
    $cmd .= " -em "      . getGlobal("utgMergeErrorRate")  . " \\\n";
    $cmd .= " -er "      . getGlobal("utgRepeatErrorRate") . " \\\n";
    $cmd .= " -threads " . getGlobal("batThreads")         . " \\\n"   if defined(getGlobal("batThreads"));
    $cmd .= " -M "       . getGlobal("batMemory")          . " \\\n"   if defined(getGlobal("batMemory"));
    $cmd .= " "          . getGlobal("batOptions")         . " \\\n"   if defined(getGlobal("batOptions"));
    $cmd .= " > $wrk/4-unitigger/unitigger.err 2>&1";

    return($cmd);
}



sub reportSizes ($$$$) {
    my $wrk       = shift @_;
    my $asm       = shift @_;
    my $version   = shift @_;
    my $label     = shift @_;

    my $bin       = getBinDirectory();
    my $cmd       = "";

    my $asmnum    = 0;
    my $asmbases  = 0;
    my $asmsizes  = "";

    my $singnum   = 0;
    my $singbases = 0;

    if (! -e "$wrk/$asm.tigStore.sizes.1.beforeConsensus") {
        $cmd  = "$bin/tgStoreDump \\\n";
        $cmd .= "  -G $wrk/$asm.gkpStore \\\n";
        $cmd .= "  -T $wrk/$asm.tigStore $version \\\n";
        $cmd .= "  -U \\\n";
        $cmd .= "  -d sizes \\\n";
        $cmd .= "  -s " . getGlobal("genomeSize") . " \\\n";
        $cmd .= "> $wrk/$asm.tigStore.sizes.1.beforeConsensus\n";

        if (runCommand($wrk, $cmd)) {
            caExit("failed to generate unitig sizes", undef);
        }
    }

    open(F, "< $wrk/$asm.tigStore.sizes.1.beforeConsensus") or caExit("failed to open '$wrk/$asm.tigStore.sizes.1.beforeConsensus' for reading: $!\n", undef);
    while (<F>) {
        $singbases = $1  if (m/lenSingleton\s+sum\s+(\d+)/);
        $singnum   = $1  if (m/lenSingleton\s+num\s+(\d+)/);
        $asmbases  = $1  if (m/lenAssembled\s+sum\s+(\d+)/);
        $asmnum    = $1  if (m/lenAssembled\s+num\s+(\d+)/);

        $asmsizes .= "--    $_"  if (m/lenAssembled\s+(n\d+)\s+siz/);
    }
    close(F);

    print STDERR "--  Found, in version $version, $label:\n";
    print STDERR "--    unitigs:     $asmnum sequences, total length $asmbases bp.\n";
    print STDERR "--    singletons:  $singnum sequences, total length $singbases bp.\n";
    print STDERR "--\n";
    print STDERR "$asmsizes";
    print STDERR "--\n";
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
    reportSizes($wrk, $asm, 1, "after unitig construction");
    emitStage($wrk, $asm, "unitig");
  stopAfter:
    stopAfter("unitig");
}
