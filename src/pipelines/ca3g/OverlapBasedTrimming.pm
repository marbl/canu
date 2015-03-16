package ca3g::OverlapBasedTrimming;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(qualTrimReads dedupeReads trimReads splitReads dumpReads);

use strict;

use ca3g::Defaults;
use ca3g::Execution;

#  Assumes overlaps exist, in $asm.ovlStore.
#
#  Prior to calling overlapBasedTrimming, we must have computed overlaps (with -G enabled) and built a
#  ovlStore (with proper obt filtering) and possibly a dupStore.


#  Do an initial overly-permissive quality trimming, intersected with any known vector trimming.
#  This step also applies any clear range computed by MBT.
#  This step must be done before gkpStore is constructed.
#
#  Initial trimming can't be done anymore.  With no clear ranges in the store, we need to chop
#  reads.  To do initial trimming, we need to rewrite fastq, before gkpStore is created.
#
sub qualTrimReads ($$) {
    my $wrk    = shift @_;
    my $asm    = shift @_;
    my $bin    = getBinDirectory();
    my $cmd;

    return  if (-e "$wrk/$asm.gkpStore");
    return  if (-e "$wrk/0-initialtrimming/trimmed.gkp");

    $cmd  = "$bin/initialTrim \\\n";
    $cmd .= "  -log $wrk/0-overlaptrim/$asm.initialTrim.log \\\n";
    $cmd .= "  -frg $wrk/$asm.gkpStore \\\n";
    $cmd .= ">  $wrk/0-overlaptrim/$asm.initialTrim.summary \\\n";
    $cmd .= "2> $wrk/0-overlaptrim/$asm.initialTrim.err ";

    stopBefore("initialTrim", $cmd);

    if (runCommand("$wrk/0-overlaptrim", $cmd)) {
        rename "$wrk/0-overlaptrim/$asm.initialTrim.log", "$wrk/0-overlaptrim/$asm.initialTrim.log.FAILED";
        caFailure("initial trimming failed", "$wrk/0-overlaptrim/$asm.initialTrim.err");
    }
}


sub dedupeReads ($$) {
    my $wrk    = shift @_;
    my $asm    = shift @_;
    my $bin    = getBinDirectory();
    my $cmd;

    return  if (getGlobal("doDeDuplication") == 0);

    return  if (-e "$wrk/3-overlapbasedtrimming/deduplicated");

    caFailure("didn't find '$wrk/$asm.dupStore' with overlaps to use for duplicate finding", undef)  if (! -e "$wrk/$asm.dupStore");

    $cmd  = "$bin/deduplicate \\\n";
    $cmd .= "  -gkp     $wrk/$asm.gkpStore \\\n";
    $cmd .= "  -ovs     $wrk/0-overlaptrim/$asm.obtStore \\\n";
    $cmd .= "  -ovs     $wrk/0-overlaptrim/$asm.dupStore \\\n";
    $cmd .= "  -report  $wrk/0-overlaptrim/$asm.deduplicate.log \\\n";
    $cmd .= "  -summary $wrk/0-overlaptrim/$asm.deduplicate.summary \\\n";
    $cmd .= "> $wrk/0-overlaptrim/$asm.deduplicate.err 2>&1";

    stopBefore("deDuplication", $cmd);

    if (runCommand("$wrk/0-overlaptrim", $cmd)) {
        caFailure("failed to deduplicate the reads", "$wrk/0-overlaptrim/$asm.deduplicate.err");
    }

    touch("$wrk/3-overlapbasedtrimming/deduplicated");
}



sub trimReads ($$) {
    my $wrk    = shift @_;
    my $asm    = shift @_;
    my $bin    = getBinDirectory();
    my $cmd;

    return  if (-e "$wrk/3-overlapbasedtrimming/trimmed");

    #  Previously, we'd pick the error rate used by unitigger.  Now, we don't know unitigger here,
    #  and require an obt specific error rate.

    my $erate  = getGlobal("obtErrorRate");

    $cmd  = "$bin/finalTrim \\\n";
    $cmd .= "  -G $wrk/$asm.gkpStore \\\n";
    $cmd .= "  -O $wrk/0-overlaptrim/$asm.obtStore \\\n";
    $cmd .= "  -e $erate \\\n";
    $cmd .= "  -o $wrk/0-overlaptrim/$asm.finalTrim \\\n";
    $cmd .= "> $wrk/0-overlaptrim/$asm.finalTrim.err 2>&1";

    stopBefore("finalTrimming", $cmd);

    if (runCommand("$wrk/0-overlaptrim", $cmd)) {
        caFailure("failed to compute final trimming", "$wrk/0-overlaptrim/$asm.finalTrim.err");
    }

    touch("$wrk/3-overlapbasedtrimming/trimmed");
}



sub splitReads ($$) {
    my $wrk    = shift @_;
    my $asm    = shift @_;
    my $bin    = getBinDirectory();
    my $cmd;

    return  if (getGlobal("doChimeraDetection") eq 'off');

    return  if (-e "$wrk/3-overlapbasedtrimming/splitted");  #  Splitted?

    my $erate  = getGlobal("ovlErrorRate");  #  Was this historically

    $cmd  = "$bin/chimera \\\n";
    $cmd .= "  -G $wrk/$asm.gkpStore \\\n";
    $cmd .= "  -O $wrk/0-overlaptrim/$asm.obtStore \\\n";
    $cmd .= "  -e $erate \\\n";
    $cmd .= "  -o $wrk/0-overlaptrim/$asm.chimera \\\n";
    $cmd .= "  -mininniepair 0 -minoverhanging 0 \\\n" if (getGlobal("doChimeraDetection") eq "aggressive");
    $cmd .= "> $wrk/0-overlaptrim/$asm.chimera.err 2>&1";

    stopBefore("chimeraDetection", $cmd);

    if (runCommand("$wrk/0-overlaptrim", $cmd)) {
        caFailure("chimera cleaning failed", "$wrk/0-overlaptrim/$asm.chimera.err");
    }

    touch("$wrk/3-overlapbasedtrimming/splitted", "Splitted?  Is that even a word?");
}



sub dumpReads ($$) {
    my $wrk    = shift @_;
    my $asm    = shift @_;
    my $bin    = getBinDirectory();
    my $cmd;

    return  if (-e "$wrk/$asm.trimmed.gkp");

    $cmd  = "$bin/gatekeeperDump \\\n";
    $cmd .= "  -dumpfastq $wrk/$asm.trimmed \\\n";
    $cmd .= "  -withlibnames \\\n";
    $cmd .= "  -g gkpStore \\\n";
    $cmd .= "> $wrk/$asm.trimmed.err 2>&1\n";

    if (runCommand($wrk, $cmd)) {
        caFailure("dumping trimmed reads failed", "$wrk/$asm.trimmed.err");
    }
}
