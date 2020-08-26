
###############################################################################
 #
 #  This file is part of canu, a software program that assembles whole-genome
 #  sequencing reads into contigs.
 #
 #  This software is based on:
 #    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 #    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 #
 #  Except as indicated otherwise, this is a 'United States Government Work',
 #  and is released in the public domain.
 #
 #  File 'README.licenses' in the root directory of this distribution
 #  contains full conditions and disclaimers.
 ##

package canu::OverlapBasedTrimming;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(qualTrimReads dedupeReads trimReads splitReads loadTrimmedReads dumpTrimmedReads);

use strict;
use warnings "all";
no  warnings "uninitialized";

use File::Path 2.08 qw(make_path remove_tree);

use canu::Defaults;
use canu::Execution;

use canu::SequenceStore;
use canu::Report;
use canu::Output;

use canu::Grid_Cloud;


sub trimReads ($) {
    my $asm    = shift @_;
    my $bin    = getBinDirectory();
    my $cmd;
    my $path   = "trimming/3-overlapbasedtrimming";

    goto allDone   if (fileExists("trimming/3-overlapbasedtrimming/$asm.1.trimReads.clear"));

    make_path($path)  if (! -d $path);

    fetchOvlStore($asm, "trimming");

    #  Previously, we'd pick the error rate used by unitigger.  Now, we don't know unitigger here,
    #  and require an obt specific error rate.

    $cmd  = "$bin/trimReads \\\n";
    $cmd .= "  -S  ../../$asm.seqStore \\\n";
    $cmd .= "  -O  ../$asm.ovlStore \\\n";
    $cmd .= "  -Co ./$asm.1.trimReads.clear \\\n";
    $cmd .= "  -e  " . getGlobal("obtErrorRate") . " \\\n";
    $cmd .= "  -minlength " . getGlobal("minReadLength") . " \\\n";
    #$cmd .= "  -Cm ./$asm.max.clear \\\n"          if (-e "./$asm.max.clear");
    $cmd .= "  -ol " . getGlobal("trimReadsOverlap") . " \\\n";
    $cmd .= "  -oc " . getGlobal("trimReadsCoverage") . " \\\n";
    $cmd .= "  -o  ./$asm.1.trimReads \\\n";
    $cmd .= ">     ./$asm.1.trimReads.err 2>&1";

    if (runCommand($path, $cmd)) {
        caFailure("trimReads failed", "$path/$asm.1.trimReads.err");
    }

    caFailure("trimReads finished, but no '$asm.1.trimReads.clear' output found", undef)  if (! -e "$path/$asm.1.trimReads.clear");

    unlink("$path/$asm.1.trimReads.err");

    stashFile("./trimming/3-overlapbasedtrimming/$asm.1.trimReads.clear");

    my $report;

#FORMAT
    open(F, "< trimming/3-overlapbasedtrimming/$asm.1.trimReads.stats") or caExit("can't open 'trimming/3-overlapbasedtrimming/$asm.1.trimReads.stats' for reading: $!", undef);
    while (<F>) {
        $report .= "--  $_";
    }
    close(F);

    addToReport("trimming", $report);


    if (0) {
        $cmd  = "$bin/sqStoreDumpFASTQ \\\n";
        $cmd .= "  -S ../../$asm.seqStore \\\n";
        $cmd .= "  -c ./$asm.1.trimReads.clear \\\n";
        $cmd .= "  -o ./$asm.1.trimReads.trimmed \\\n";
        $cmd .= ">    ./$asm.1.trimReads.trimmed.err 2>&1";

        if (runCommand($path, $cmd)) {
            caFailure("dumping trimmed reads failed", "$path/$asm.1.trimReads.trimmed.err");
        }
    }

  finishStage:
    generateReport($asm);
    resetIteration("obt-trimReads");

  allDone:
}



sub splitReads ($) {
    my $asm    = shift @_;
    my $bin    = getBinDirectory();
    my $cmd;
    my $path   = "trimming/3-overlapbasedtrimming";

    goto allDone   if (fileExists("trimming/3-overlapbasedtrimming/$asm.2.splitReads.clear"));

    make_path($path)  if (! -d $path);

    fetchOvlStore($asm, "trimming");

    fetchFile("./trimming/3-overlapbasedtrimming/$asm.1.trimReads.clear");

    my $erate  = getGlobal("obtErrorRate");  #  Was this historically

    #$cmd .= "  -mininniepair 0 -minoverhanging 0 \\\n" if (getGlobal("doChimeraDetection") eq "aggressive");

    $cmd  = "$bin/splitReads \\\n";
    $cmd .= "  -S  ../../$asm.seqStore \\\n";
    $cmd .= "  -O  ../$asm.ovlStore \\\n";
    $cmd .= "  -Ci ./$asm.1.trimReads.clear \\\n"       if (-e "trimming/3-overlapbasedtrimming/$asm.1.trimReads.clear");
    #$cmd .= "  -Cm ./$asm.max.clear \\\n"               if (-e "trimming/3-overlapbasedtrimming/$asm.max.clear");
    $cmd .= "  -Co ./$asm.2.splitReads.clear \\\n";
    $cmd .= "  -e  $erate \\\n";
    $cmd .= "  -minlength " . getGlobal("minReadLength") . " \\\n";
    $cmd .= "  -o  ./$asm.2.splitReads \\\n";
    $cmd .= ">     ./$asm.2.splitReads.err 2>&1";

    if (runCommand($path, $cmd)) {
        caFailure("splitReads failed", "$path/$asm.2.splitReads.err");
    }

    caFailure("splitReads finished, but no '$asm.2.splitReads.clear' output found", undef)  if (! -e "$path/$asm.2.splitReads.clear");

    unlink("$path/$asm.2.splitReads.err");

    stashFile("./trimming/3-overlapbasedtrimming/$asm.2.splitReads.clear");

    my $report;

#FORMAT
    open(F, "< trimming/3-overlapbasedtrimming/$asm.2.splitReads.stats") or caExit("can't open 'trimming/3-overlapbasedtrimming/$asm.2.splitReads.stats' for reading: $!", undef);
    while (<F>) {
        $report .= "--  $_";
    }
    close(F);

    addToReport("splitting", $report);

    if (0) {
        $cmd  = "$bin/sqStoreDumpFASTQ \\\n";
        $cmd .= "  -S ../../$asm.seqStore \\\n";
        $cmd .= "  -c ./$asm.2.splitReads.clear \\\n";
        $cmd .= "  -o ./$asm.2.splitReads.trimmed \\\n";
        $cmd .= ">    ./$asm.2.splitReads.trimmed.err 2>&1";

        if (runCommand($path, $cmd)) {
            caFailure("dumping trimmed reads failed", "$path/$asm.2.splitReads.trimmed.err");
        }
    }

  finishStage:
    generateReport($asm);
    resetIteration("obt-splitReads");

  allDone:
}



sub loadTrimmedReads ($) {
    my $asm    = shift @_;
    my $bin    = getBinDirectory();
    my $cmd;
    my $path   = "trimming/3-overlapbasedtrimming";
    my $inp;

    goto allDone   if (getNumberOfBasesInStore($asm, "utg") > 0);

    make_path($path)  if (! -d $path);

    fetchFile("./trimming/3-overlapbasedtrimming/$asm.1.trimReads.clear");
    fetchFile("./trimming/3-overlapbasedtrimming/$asm.2.splitReads.clear");

    $inp = "./$asm.1.trimReads.clear"   if (-e "$path/$asm.1.trimReads.clear");
    $inp = "./$asm.2.splitReads.clear"  if (-e "$path/$asm.2.splitReads.clear");

    caFailure("loading trimmed reads failed; no 'clear' input", "trimming/$asm.trimmedReads.err")  if (!defined($inp));

    $cmd  = "$bin/loadTrimmedReads \\\n";
    $cmd .= "  -S ../../$asm.seqStore \\\n";
    $cmd .= "  -c $inp \\\n";
    $cmd .= "> ./$asm.loadTrimmedReads.err 2>&1";

    if (runCommand($path, $cmd)) {
        caFailure("loading clear ranges failed", "$path/$asm.loadTrimmedReads.err");
    }

    unlink("$path/$asm.loadTrimmedReads.err");

    #  Summarize and save results.

    generateReadLengthHistogram("utg", $asm);
    stashSeqStore($asm);

    #  If no trimmed reads were generated, stop now so we
    #  preserve the stores.

    if (getNumberOfReadsInStore($asm, "utg") > 0) {
        print STDERR "--\n";
        print STDERR "-- No trimmed reads generated, overlaps used for trimming saved.\n";
    }
    elsif (getGlobal("saveOverlaps") eq "0") {
        print STDERR "--\n";
        print STDERR "-- Purging overlaps used for trimming.\n";
        remove_tree("trimming/$asm.ovlStore")
    }
    else {
        print STDERR "--\n";
        print STDERR "-- Overlaps used for trimming saved.\n";
    }

  finishStage:
    generateReport($asm);
    resetIteration("obt-dumpReads");

  allDone:
}



sub dumpTrimmedReads ($) {
    my $asm     = shift @_;
    my $bin     = getBinDirectory();
    my $cmd;

    goto allDone   if (fileExists("$asm.trimmedReads.fasta.gz"));
    goto allDone   if (getGlobal("saveReads") == 0);

    $cmd  = "$bin/sqStoreDumpFASTQ \\\n";
    $cmd .= "  -trimmed \\\n";
    $cmd .= "  -S ./$asm.seqStore \\\n";
    $cmd .= "  -o ./$asm.trimmedReads.gz \\\n";
    $cmd .= "  -fasta \\\n";
    $cmd .= "  -trimmed -normal -nolibname \\\n";
    $cmd .= "> ./$asm.trimmedReads.fasta.err 2>&1";

    if (runCommand(".", $cmd)) {
        caExit("failed to output trimmed reads", "./$asm.trimmedReads.fasta.err");
    }

    unlink "./$asm.trimmedReads.fasta.err";

    stashFile("$asm.trimmedReads.fasta.gz");

    print STDERR "--\n";
    print STDERR "-- Trimmed reads saved in '$asm.trimmedReads.fasta.gz'.\n";

  finishStage:
    generateReport($asm);
    resetIteration("cor-dumpTrimmedReads");

  allDone:
    stopAfter("trimming");
}
