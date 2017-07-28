
###############################################################################
 #
 #  This file is part of canu, a software program that assembles whole-genome
 #  sequencing reads into contigs.
 #
 #  This software is based on:
 #    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 #    the 'kmer package' (http://kmer.sourceforge.net)
 #  both originally distributed by Applera Corporation under the GNU General
 #  Public License, version 2.
 #
 #  Canu branched from Celera Assembler at its revision 4587.
 #  Canu branched from the kmer project at its revision 1994.
 #
 #  This file is derived from:
 #
 #    src/pipelines/ca3g/OverlapBasedTrimming.pm
 #
 #  Modifications by:
 #
 #    Brian P. Walenz from 2015-MAR-16 to 2015-AUG-25
 #      are Copyright 2015 Battelle National Biodefense Institute, and
 #      are subject to the BSD 3-Clause License
 #
 #    Brian P. Walenz beginning on 2015-NOV-04
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #    Sergey Koren beginning on 2017-MAR-03
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

package canu::OverlapBasedTrimming;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(qualTrimReads dedupeReads trimReads splitReads dumpReads);

use strict;

use File::Path 2.08 qw(make_path remove_tree);

use canu::Defaults;
use canu::Execution;
use canu::Gatekeeper;
use canu::Report;
use canu::HTML;
use canu::Grid_Cloud;


sub trimReads ($) {
    my $asm    = shift @_;
    my $bin    = getBinDirectory();
    my $cmd;
    my $path   = "trimming/3-overlapbasedtrimming";

    goto allDone   if (skipStage($asm, "obt-trimReads") == 1);
    goto allDone   if (fileExists("trimming/3-overlapbasedtrimming/$asm.1.trimReads.clear"));

    make_path($path)  if (! -d $path);

    fetchStore("./trimming/$asm.ovlStore");

    #  Previously, we'd pick the error rate used by unitigger.  Now, we don't know unitigger here,
    #  and require an obt specific error rate.

    $cmd  = "$bin/trimReads \\\n";
    $cmd .= "  -G  ../$asm.gkpStore \\\n";
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
        $cmd  = "$bin/gatekeeperDumpFASTQ \\\n";
        $cmd .= "  -G ../$asm.gkpStore \\\n";
        $cmd .= "  -c ./$asm.1.trimReads.clear \\\n";
        $cmd .= "  -o ./$asm.1.trimReads.trimmed \\\n";
        $cmd .= ">    ./$asm.1.trimReads.trimmed.err 2>&1";

        if (runCommand($path, $cmd)) {
            caFailure("dumping trimmed reads failed", "$path/$asm.1.trimReads.trimmed.err");
        }
    }

  finishStage:
    emitStage($asm, "obt-trimReads");
    buildHTML($asm, "obt");

  allDone:
}



sub splitReads ($) {
    my $asm    = shift @_;
    my $bin    = getBinDirectory();
    my $cmd;
    my $path   = "trimming/3-overlapbasedtrimming";

    goto allDone   if (skipStage($asm, "obt-splitReads") == 1);
    goto allDone   if (fileExists("trimming/3-overlapbasedtrimming/$asm.2.splitReads.clear"));

    make_path($path)  if (! -d $path);

    fetchStore("./trimming/$asm.ovlStore");
    fetchFile("./trimming/3-overlapbasedtrimming/$asm.1.trimReads.clear");

    my $erate  = getGlobal("obtErrorRate");  #  Was this historically

    #$cmd .= "  -mininniepair 0 -minoverhanging 0 \\\n" if (getGlobal("doChimeraDetection") eq "aggressive");

    $cmd  = "$bin/splitReads \\\n";
    $cmd .= "  -G  ../$asm.gkpStore \\\n";
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
        $cmd  = "$bin/gatekeeperDumpFASTQ \\\n";
        $cmd .= "  -G ../$asm.gkpStore \\\n";
        $cmd .= "  -c ./$asm.2.splitReads.clear \\\n";
        $cmd .= "  -o ./$asm.2.splitReads.trimmed \\\n";
        $cmd .= ">    ./$asm.2.splitReads.trimmed.err 2>&1";

        if (runCommand($path, $cmd)) {
            caFailure("dumping trimmed reads failed", "$path/$asm.2.splitReads.trimmed.err");
        }
    }

  finishStage:
    emitStage($asm, "obt-splitReads");
    buildHTML($asm, "obt");

  allDone:
}



sub dumpReads ($) {
    my $asm    = shift @_;
    my $bin    = getBinDirectory();
    my $cmd;
    my $path   = "trimming/3-overlapbasedtrimming";
    my $inp;

    goto allDone   if (skipStage($asm, "obt-dumpReads") == 1);
    goto allDone   if (sequenceFileExists("$asm.trimmedReads"));

    make_path($path)  if (! -d $path);

    fetchFile("./trimming/3-overlapbasedtrimming/$asm.1.trimReads.clear");
    fetchFile("./trimming/3-overlapbasedtrimming/$asm.2.splitReads.clear");

    $inp = "./3-overlapbasedtrimming/$asm.1.trimReads.clear"   if (-e "$path/$asm.1.trimReads.clear");
    $inp = "./3-overlapbasedtrimming/$asm.2.splitReads.clear"  if (-e "$path/$asm.2.splitReads.clear");

    caFailure("dumping trimmed reads failed; no 'clear' input", "trimming/$asm.trimmedReads.err")  if (!defined($inp));

    $cmd  = "$bin/gatekeeperDumpFASTQ -fasta -nolibname \\\n";
    $cmd .= "  -G ./$asm.gkpStore \\\n";
    $cmd .= "  -c $inp \\\n";
    $cmd .= "  -o ../$asm.trimmedReads.gz \\\n";     #  Adds .fasta
    $cmd .= ">    ../$asm.trimmedReads.err 2>&1";

    if (runCommand("trimming", $cmd)) {
        caFailure("dumping trimmed reads failed", "./$asm.trimmedReads.err");
    }

    unlink("./$asm.trimmedReads.err");

    stashFile("./$asm.trimmedReads.fasta.gz");

    remove_tree("trimming/$asm.ovlStore")   if (getGlobal("saveOverlaps") eq "0");

  finishStage:
    emitStage($asm, "obt-dumpReads");
    buildHTML($asm, "obt");

  allDone:
    print STDERR "--\n";
    print STDERR "-- Trimmed reads saved in 'trimming/$asm.trimmedReads.fasta.gz'\n";

    stopAfter("readTrimming");
}
