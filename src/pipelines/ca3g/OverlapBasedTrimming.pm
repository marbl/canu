package ca3g::OverlapBasedTrimming;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(qualTrimReads dedupeReads trimReads splitReads dumpReads);

use strict;

use File::Path qw(make_path remove_tree);

use ca3g::Defaults;
use ca3g::Execution;

#  Assumes overlaps exist, in $asm.ovlStore.
#
#  Prior to calling overlapBasedTrimming, we must have computed overlaps (with -G enabled) and built an
#  ovlStore.


#  Step 1:  Remove vector/adapter/other obvious junk.  This MUST chop the reads before loading into
#  gatekeeper.  Gatekeeper, and overlapper, do not know about clear ranges.  The notion of a clear
#  range is only valid during the OBT compute itself; the output of OBT is a new set of chopped
#  reads.  This also GREATLY simplifies the other modules.
#
#  Step 2:  Reads are deduplicated.  This doesn't change the clear range at all, just marks some
#  reads as deleted.  The program does take an input clear range, but it's not ever used in the
#  current pipeline.
#
#  Step 3:  Trimming.
#
#  Step 4:  Splitting.  Chimera and subreads are detected and cleaned.



sub dedupeReads ($$) {
    my $wrk    = shift @_;
    my $asm    = shift @_;
    my $bin    = getBinDirectory();
    my $cmd;
    my $path   = "$wrk/3-overlapbasedtrimming";

    return  if (getGlobal("doDeDuplication") == 0);
    return  if (-e "$path/deduplicated");

    make_path($path)  if (! -d $path);

    $cmd  = "$bin/deduplicate \\\n";
    $cmd .= "  -G  $wrk/$asm.gkpStore \\\n";
    $cmd .= "  -O  $wrk/$asm.ovlStore \\\n";
    $cmd .= "  -e  " . getGlobal("obtErrorRate") . " \\\n";
    $cmd .= "  -Ci $path/$asm.latest.clear \\\n"  if (-e "$path/$asm.latest.clear");
    $cmd .= "  -Co $path/$asm.dedupe.clear \\\n";
    $cmd .= "  -o  $path/$asm.dedupe \\\n";
    $cmd .= "> $path/$asm.dedupe.err 2>&1";

    stopBefore("deDuplication", $cmd);

    if (runCommand($path, $cmd)) {
        caFailure("failed to deduplicate the reads", "$path/$asm.dedupe.err");
    }

    caFailure("dedupe finished, but no '$asm.dedupe.clear' output found", undef)  if (! -e "$path/$asm.dedupe.clear");

    unlink("$path/$asm.latest.clear");
    symlink("$path/$asm.dedupe.clear", "$path/$asm.latest.clear");

    touch("$path/deduplicated");
}



sub trimReads ($$) {
    my $wrk    = shift @_;
    my $asm    = shift @_;
    my $bin    = getBinDirectory();
    my $cmd;
    my $path   = "$wrk/3-overlapbasedtrimming";

    return  if (-e "$path/trimmed");

    make_path($path)  if (! -d $path);

    #  Previously, we'd pick the error rate used by unitigger.  Now, we don't know unitigger here,
    #  and require an obt specific error rate.

    $cmd  = "$bin/finalTrim \\\n";
    $cmd .= "  -G  $wrk/$asm.gkpStore \\\n";
    $cmd .= "  -O  $wrk/$asm.ovlStore \\\n";
    $cmd .= "  -e  " . getGlobal("obtErrorRate") . " \\\n";
    $cmd .= "  -Ci $path/$asm.latest.clear \\\n"       if (-e "$path/$asm.latest.clear");
    $cmd .= "  -Cm $path/$asm.max.clear \\\n"          if (-e "$path/$asm.max.clear");      #  Probably not useful anymore
    $cmd .= "  -Co $path/$asm.finalTrim.clear \\\n";
    $cmd .= "  -o  $path/$asm.finalTrim \\\n";
    $cmd .= "> $path/$asm.finalTrim.err 2>&1";

    stopBefore("finalTrimming", $cmd);

    if (runCommand($path, $cmd)) {
        caFailure("failed to compute final trimming", "$path/$asm.finalTrim.err");
    }

    caFailure("finalTrim finished, but no '$asm.finalTrim.clear' output found", undef)  if (! -e "$path/$asm.finalTrim.clear");

    unlink("$path/$asm.latest.clear");
    symlink("$path/$asm.finalTrim.clear", "$path/$asm.latest.clear");

    touch("$path/trimmed");
}



sub splitReads ($$) {
    my $wrk    = shift @_;
    my $asm    = shift @_;
    my $bin    = getBinDirectory();
    my $cmd;
    my $path   = "$wrk/3-overlapbasedtrimming";

    return  if (getGlobal("doChimeraDetection") eq 'off');
    return  if (-e "$path/splitted");  #  Splitted?

    make_path($path)  if (! -d $path);

    my $erate  = getGlobal("ovlErrorRate");  #  Was this historically

    $cmd  = "$bin/chimera \\\n";
    $cmd .= "  -G  $wrk/$asm.gkpStore \\\n";
    $cmd .= "  -O  $wrk/$asm.ovlStore \\\n";
    $cmd .= "  -e  $erate \\\n";
    $cmd .= "  -o  $path/$asm.chimera \\\n";
    $cmd .= "  -Ci $path/$asm.latest.clear \\\n"       if (-e "$path/$asm.latest.clear");
    $cmd .= "  -Cm $path/$asm.max.clear \\\n"          if (-e "$path/$asm.max.clear");
    $cmd .= "  -Co $path/$asm.chimera.clear \\\n";
    $cmd .= "  -mininniepair 0 -minoverhanging 0 \\\n" if (getGlobal("doChimeraDetection") eq "aggressive");
    $cmd .= "> $path/$asm.chimera.err 2>&1";

    stopBefore("chimeraDetection", $cmd);

    if (runCommand($path, $cmd)) {
        caFailure("chimera cleaning failed", "$path/$asm.chimera.err");
    }

    caFailure("chimera finished, but no '$asm.chimera.clear' output found", undef)  if (! -e "$path/$asm.chimera.clear");

    unlink("$path/$asm.latest.clear");
    symlink("$path/$asm.chimera.clear", "$path/$asm.latest.clear");

    touch("$path/splitted", "Splitted?  Is that even a word?");
}



sub dumpReads ($$) {
    my $wrk    = shift @_;
    my $asm    = shift @_;
    my $bin    = getBinDirectory();
    my $cmd;
    my $path   = "$wrk/3-overlapbasedtrimming";

    return  if (-e "$wrk/$asm.trimmed.gkp");

    make_path($path)  if (! -d $path);

    $cmd  = "$bin/gatekeeperDump \\\n";
    $cmd .= "  -g gkpStore \\\n";
    $cmd .= "  -dumpfastq $wrk/$asm.trimmed \\\n";
    $cmd .= "  -withlibnames \\\n";
    $cmd .= "> $wrk/$asm.trimmed.err 2>&1\n";

    #$cmd .= "  -Ci $path/$asm.initialTrim.clear \\\n"  if (-e "$path/$asm.initialTrim.clear");
    #$cmd .= "  -Ci $path/$asm.dedupe.clear \\\n"       if (-e "$path/$asm.dedupe.clear");
    #$cmd .= "  -Ci $path/$asm.finalTrim.clear \\\n"    if (-e "$path/$asm.finalTrim.clear");
    #$cmd .= "  -Cm $path/$asm.max.clear \\\n"          if (-e "$path/$asm.max.clear");

    if (runCommand($wrk, $cmd)) {
        caFailure("dumping trimmed reads failed", "$wrk/$asm.trimmed.err");
    }
}
