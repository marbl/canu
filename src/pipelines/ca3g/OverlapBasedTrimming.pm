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

    $cmd  = "$bin/trimReads \\\n";
    $cmd .= "  -G  $wrk/$asm.gkpStore \\\n";
    $cmd .= "  -O  $wrk/$asm.ovlStore \\\n";
    $cmd .= "  -e  " . getGlobal("obtErrorRate") . " \\\n";
    $cmd .= "  -Ci $path/$asm.latest.clear \\\n"       if (-e "$path/$asm.latest.clear");
    $cmd .= "  -Cm $path/$asm.max.clear \\\n"          if (-e "$path/$asm.max.clear");      #  Probably not useful anymore
    $cmd .= "  -Co $path/$asm.$idx.trimReads.clear \\\n";
    $cmd .= "  -o  $path/$asm.$idx.trimReads \\\n";
    $cmd .= ">     $path/$asm.$idx.trimReads.err 2>&1";

    stopBefore("trimReads", $cmd);

    if (runCommand($path, $cmd)) {
        caFailure("trimReads failed", "$path/$asm.$idx.trimReads.err");
    }

    caFailure("trimReads finished, but no '$asm.$idx.trimReads.clear' output found", undef)  if (! -e "$path/$asm.$idx.trimReads.clear");

    unlink("$path/$asm.latest.clear");
    symlink("$path/$asm.$idx.trimReads.clear", "$path/$asm.latest.clear");

    $cmd  = "$bin/gatekeeperDumpFASTQ \\\n";
    $cmd .= "  -G $wrk/$asm.gkpStore \\\n";
    $cmd .= "  -c $path/$asm.$idx.trimReads.clear \\\n";
    $cmd .= "  -o $path/$asm.$idx.trimReads.trimmed \\\n";
    $cmd .= ">    $path/$asm.$idx.trimReads.trimmed.err 2>&1\n";

    if (runCommand($wrk, $cmd)) {
        caFailure("dumping trimmed reads failed", "$wrk/$asm.$idx.trimReads.trimmed.err");
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

    $cmd  = "$bin/splitReads \\\n";
    $cmd .= "  -G  $wrk/$asm.gkpStore \\\n";
    $cmd .= "  -O  $wrk/$asm.ovlStore \\\n";
    $cmd .= "  -e  $erate \\\n";
    $cmd .= "  -o  $path/$asm.$idx.splitReads \\\n";
    $cmd .= "  -Ci $path/$asm.latest.clear \\\n"       if (-e "$path/$asm.latest.clear");
    $cmd .= "  -Cm $path/$asm.max.clear \\\n"          if (-e "$path/$asm.max.clear");
    $cmd .= "  -Co $path/$asm.$idx.splitReads.clear \\\n";
    $cmd .= ">     $path/$asm.$idx.splitReads.err 2>&1";

    stopBefore("splitReads", $cmd);

    if (runCommand($path, $cmd)) {
        caFailure("splitReads failed", "$path/$asm.$idx.splitReads.err");
    }

    caFailure("splitReads finished, but no '$asm.$idx.splitReads.clear' output found", undef)  if (! -e "$path/$asm.$idx.splitReads.clear");

    unlink("$path/$asm.latest.clear");
    symlink("$path/$asm.$idx.splitReads.clear", "$path/$asm.latest.clear");

    $cmd  = "$bin/gatekeeperDumpFASTQ \\\n";
    $cmd .= "  -G $wrk/$asm.gkpStore \\\n";
    $cmd .= "  -c $path/$asm.$idx.splitReads.clear \\\n";
    $cmd .= "  -o $path/$asm.$idx.splitReads.trimmed \\\n";
    $cmd .= ">    $path/$asm.$idx.splitReads.trimmed.err 2>&1\n";

    if (runCommand($wrk, $cmd)) {
        caFailure("dumping trimmed reads failed", "$wrk/$asm.$idx.splitReads.trimmed.err");
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
