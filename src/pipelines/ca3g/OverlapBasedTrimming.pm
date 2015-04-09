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



sub dedupeReads ($$$) {
    my $idx    = shift @_;
    my $wrk    = shift @_;
    my $asm    = shift @_;
    my $bin    = getBinDirectory();
    my $cmd;
    my $path   = "$wrk/3-overlapbasedtrimming";

    return  if (-e "$path/deduplicated");

    make_path($path)  if (! -d $path);

    $cmd  = "$bin/deduplicate \\\n";
    $cmd .= "  -G  $wrk/$asm.gkpStore \\\n";
    $cmd .= "  -O  $wrk/$asm.ovlStore \\\n";
    $cmd .= "  -e  " . getGlobal("obtErrorRate") . " \\\n";
    $cmd .= "  -Ci $path/$asm.latest.clear \\\n"  if (-e "$path/$asm.latest.clear");
    $cmd .= "  -Co $path/$asm.$idx.dedupe.clear \\\n";
    $cmd .= "  -o  $path/$asm.$idx.dedupe \\\n";
    $cmd .= ">     $path/$asm.$idx.dedupe.err 2>&1";

    stopBefore("deDuplication", $cmd);

    if (runCommand($path, $cmd)) {
        caFailure("failed to deduplicate the reads", "$path/$asm.$idx.dedupe.err");
    }

    caFailure("dedupe finished, but no '$asm.$idx.dedupe.clear' output found", undef)  if (! -e "$path/$asm.$idx.dedupe.clear");

    unlink("$path/$asm.latest.clear");
    symlink("$path/$asm.$idx.dedupe.clear", "$path/$asm.latest.clear");

    $cmd  = "$bin/gatekeeperDumpFASTQ \\\n";
    $cmd .= "  -G $wrk/$asm.gkpStore \\\n";
    $cmd .= "  -c $path/$asm.$idx.dedupe.clear \\\n";
    $cmd .= "  -o $path/$asm.$idx.dedupe.trimmed \\\n";
    $cmd .= ">    $path/$asm.$idx.dedupe.trimmed.err 2>&1\n";

    if (runCommand($wrk, $cmd)) {
        caFailure("dumping trimmed reads failed", "$wrk/$asm.$idx.dedupe.trimmed.err");
    }

    $cmd  = "$bin/gatekeeperDumpFASTQ \\\n";
    $cmd .= "  -G $wrk/$asm.gkpStore \\\n";
    $cmd .= "  -c $path/$asm.$idx.dedupe.clear \\\n";
    $cmd .= "  -onlydeleted \\\n";
    $cmd .= "  -o $path/$asm.$idx.dedupe.deleted \\\n";
    $cmd .= ">    $path/$asm.$idx.dedupe.deleted.err 2>&1\n";

    if (runCommand($wrk, $cmd)) {
        caFailure("dumping trimmed reads failed", "$wrk/$asm.$idx.dedupe.trimmed.err");
    }

    touch("$path/deduplicated");
}



sub trimReads ($$$) {
    my $idx    = shift @_;
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
    $cmd .= "  -Co $path/$asm.$idx.finalTrim.clear \\\n";
    $cmd .= "  -o  $path/$asm.$idx.finalTrim \\\n";
    $cmd .= ">     $path/$asm.$idx.finalTrim.err 2>&1";

    stopBefore("finalTrimming", $cmd);

    if (runCommand($path, $cmd)) {
        caFailure("failed to compute final trimming", "$path/$asm.$idx.finalTrim.err");
    }

    caFailure("finalTrim finished, but no '$asm.$idx.finalTrim.clear' output found", undef)  if (! -e "$path/$asm.$idx.finalTrim.clear");

    unlink("$path/$asm.latest.clear");
    symlink("$path/$asm.$idx.finalTrim.clear", "$path/$asm.latest.clear");

    $cmd  = "$bin/gatekeeperDumpFASTQ \\\n";
    $cmd .= "  -G $wrk/$asm.gkpStore \\\n";
    $cmd .= "  -c $path/$asm.$idx.finalTrim.clear \\\n";
    $cmd .= "  -o $path/$asm.$idx.finalTrim.trimmed \\\n";
    $cmd .= ">    $path/$asm.$idx.finalTrim.trimmed.err 2>&1\n";

    if (runCommand($wrk, $cmd)) {
        caFailure("dumping trimmed reads failed", "$wrk/$asm.$idx.finalTrim.trimmed.err");
    }

    touch("$path/trimmed");
}



sub splitReads ($$$) {
    my $idx    = shift @_;
    my $wrk    = shift @_;
    my $asm    = shift @_;
    my $bin    = getBinDirectory();
    my $cmd;
    my $path   = "$wrk/3-overlapbasedtrimming";

    return  if (-e "$path/splitted");  #  Splitted?

    make_path($path)  if (! -d $path);

    my $erate  = getGlobal("ovlErrorRate");  #  Was this historically

    #$cmd .= "  -mininniepair 0 -minoverhanging 0 \\\n" if (getGlobal("doChimeraDetection") eq "aggressive");

    $cmd  = "$bin/chimera \\\n";
    $cmd .= "  -G  $wrk/$asm.gkpStore \\\n";
    $cmd .= "  -O  $wrk/$asm.ovlStore \\\n";
    $cmd .= "  -e  $erate \\\n";
    $cmd .= "  -o  $path/$asm.$idx.chimera \\\n";
    $cmd .= "  -Ci $path/$asm.latest.clear \\\n"       if (-e "$path/$asm.latest.clear");
    $cmd .= "  -Cm $path/$asm.max.clear \\\n"          if (-e "$path/$asm.max.clear");
    $cmd .= "  -Co $path/$asm.$idx.chimera.clear \\\n";
    $cmd .= ">     $path/$asm.$idx.chimera.err 2>&1";

    stopBefore("chimeraDetection", $cmd);

    if (runCommand($path, $cmd)) {
        caFailure("chimera cleaning failed", "$path/$asm.$idx.chimera.err");
    }

    caFailure("chimera finished, but no '$asm.$idx.chimera.clear' output found", undef)  if (! -e "$path/$asm.$idx.chimera.clear");

    unlink("$path/$asm.latest.clear");
    symlink("$path/$asm.$idx.chimera.clear", "$path/$asm.latest.clear");

    $cmd  = "$bin/gatekeeperDumpFASTQ \\\n";
    $cmd .= "  -G $wrk/$asm.gkpStore \\\n";
    $cmd .= "  -c $path/$asm.$idx.chimera.clear \\\n";
    $cmd .= "  -o $path/$asm.$idx.chimera.trimmed \\\n";
    $cmd .= ">    $path/$asm.$idx.chimera.trimmed.err 2>&1\n";

    if (runCommand($wrk, $cmd)) {
        caFailure("dumping trimmed reads failed", "$wrk/$asm.$idx.chimera.trimmed.err");
    }

    touch("$path/splitted", "Splitted?  Is that even a word?");
}



sub dumpReads ($$$) {
    my $idx    = shift @_;
    my $wrk    = shift @_;
    my $asm    = shift @_;
    my $bin    = getBinDirectory();
    my $cmd;
    my $path   = "$wrk/3-overlapbasedtrimming";

    return  if (-e "$wrk/$asm.$idx.trimmed.gkp");

    make_path($path)  if (! -d $path);

    $cmd  = "$bin/gatekeeperDumpFASTQ \\\n";
    $cmd .= "  -G $wrk/$asm.gkpStore \\\n";
    $cmd .= "  -c $path/$asm.latest.clear \\\n";
    $cmd .= "  -o $wrk/$asm.trimmed \\\n";
    $cmd .= ">    $wrk/$asm.trimmed.err 2>&1\n";

    if (runCommand($wrk, $cmd)) {
        caFailure("dumping trimmed reads failed", "$wrk/$asm.trimmed.err");
    }
}
