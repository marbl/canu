
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
 #    Brian P. Walenz beginning on 2015-MAR-16
 #      are Copyright 2015 Battelle National Biodefense Institute, and
 #      are subject to the BSD 3-Clause License
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

package canu::OverlapBasedTrimming;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(qualTrimReads dedupeReads trimReads splitReads dumpReads);

use strict;

use File::Path qw(make_path remove_tree);

use canu::Defaults;
use canu::Execution;


sub trimReads ($$) {
    my $WRK    = shift @_;
    my $wrk    = "$WRK/trimming";
    my $asm    = shift @_;
    my $bin    = getBinDirectory();
    my $cmd;
    my $path   = "$wrk/3-overlapbasedtrimming";

    return  if (skipStage($WRK, $asm, "obt-trimReads") == 1);
    return  if (-e "$path/trimmed");

    make_path($path)  if (! -d $path);

    #  Previously, we'd pick the error rate used by unitigger.  Now, we don't know unitigger here,
    #  and require an obt specific error rate.

    $cmd  = "$bin/trimReads \\\n";
    $cmd .= "  -G  $wrk/$asm.gkpStore \\\n";
    $cmd .= "  -O  $wrk/$asm.ovlStore \\\n";
    $cmd .= "  -Co $path/$asm.1.trimReads.clear \\\n";
    $cmd .= "  -e  " . getGlobal("obtErrorRate") . " \\\n";
    $cmd .= "  -minlength " . getGlobal("minReadLength") . " \\\n";
    #$cmd .= "  -Cm $path/$asm.max.clear \\\n"          if (-e "$path/$asm.max.clear");
    $cmd .= "  -ol " . getGlobal("trimReadsOverlap") . " \\\n";
    $cmd .= "  -oc " . getGlobal("trimReadsCoverage") . " \\\n";
    $cmd .= "  -o  $path/$asm.1.trimReads \\\n";
    $cmd .= ">     $path/$asm.1.trimReads.err 2>&1";

    stopBefore("trimReads", $cmd);

    if (runCommand($path, $cmd)) {
        caFailure("trimReads failed", "$path/$asm.1.trimReads.err");
    }

    caFailure("trimReads finished, but no '$asm.1.trimReads.clear' output found", undef)  if (! -e "$path/$asm.1.trimReads.clear");

    $cmd  = "$bin/gatekeeperDumpFASTQ \\\n";
    $cmd .= "  -G $wrk/$asm.gkpStore \\\n";
    $cmd .= "  -c $path/$asm.1.trimReads.clear \\\n";
    $cmd .= "  -o $path/$asm.1.trimReads.trimmed \\\n";
    $cmd .= ">    $path/$asm.1.trimReads.trimmed.err 2>&1";

    if (runCommand($path, $cmd)) {
        caFailure("dumping trimmed reads failed", "$wrk/$asm.1.trimReads.trimmed.err");
    }

    touch("$path/trimmed");
    emitStage($WRK, $asm, "obt-trimReads");
}



sub splitReads ($$) {
    my $WRK    = shift @_;
    my $wrk    = "$WRK/trimming";
    my $asm    = shift @_;
    my $bin    = getBinDirectory();
    my $cmd;
    my $path   = "$wrk/3-overlapbasedtrimming";

    return  if (skipStage($WRK, $asm, "obt-splitReads") == 1);
    return  if (-e "$path/splitted");  #  Splitted?

    make_path($path)  if (! -d $path);

    my $erate  = getGlobal("obtErrorRate");  #  Was this historically

    #$cmd .= "  -mininniepair 0 -minoverhanging 0 \\\n" if (getGlobal("doChimeraDetection") eq "aggressive");

    $cmd  = "$bin/splitReads \\\n";
    $cmd .= "  -G  $wrk/$asm.gkpStore \\\n";
    $cmd .= "  -O  $wrk/$asm.ovlStore \\\n";
    $cmd .= "  -Ci $path/$asm.1.trimReads.clear \\\n"       if (-e "$path/$asm.1.trimReads.clear");
    #$cmd .= "  -Cm $path/$asm.max.clear \\\n"               if (-e "$path/$asm.max.clear");
    $cmd .= "  -Co $path/$asm.2.splitReads.clear \\\n";
    $cmd .= "  -e  $erate \\\n";
    $cmd .= "  -minlength " . getGlobal("minReadLength") . " \\\n";
    $cmd .= "  -o  $path/$asm.2.splitReads \\\n";
    $cmd .= ">     $path/$asm.2.splitReads.err 2>&1";

    stopBefore("splitReads", $cmd);

    if (runCommand($path, $cmd)) {
        caFailure("splitReads failed", "$path/$asm.2.splitReads.err");
    }

    caFailure("splitReads finished, but no '$asm.2.splitReads.clear' output found", undef)  if (! -e "$path/$asm.2.splitReads.clear");

    $cmd  = "$bin/gatekeeperDumpFASTQ \\\n";
    $cmd .= "  -G $wrk/$asm.gkpStore \\\n";
    $cmd .= "  -c $path/$asm.2.splitReads.clear \\\n";
    $cmd .= "  -o $path/$asm.2.splitReads.trimmed \\\n";
    $cmd .= ">    $path/$asm.2.splitReads.trimmed.err 2>&1";

    if (runCommand($path, $cmd)) {
        caFailure("dumping trimmed reads failed", "$wrk/$asm.2.splitReads.trimmed.err");
    }

    touch("$path/splitted", "Splitted?  Is that even a word?");
    emitStage($WRK, $asm, "obt-splitReads");
}



sub dumpReads ($$) {
    my $WRK    = shift @_;
    my $wrk    = "$WRK/trimming";
    my $asm    = shift @_;
    my $bin    = getBinDirectory();
    my $cmd;
    my $path   = "$wrk/3-overlapbasedtrimming";
    my $inp;

    return  if (skipStage($WRK, $asm, "obt-dumpReads") == 1);
    return  if (-e "$wrk/$asm.trimmedReads.fastq");

    make_path($path)  if (! -d $path);

    $inp = "$path/$asm.1.trimReads.clear"   if (-e "$path/$asm.1.trimReads.clear");
    $inp = "$path/$asm.2.splitReads.clear"  if (-e "$path/$asm.2.splitReads.clear");

    caFailure("dumping trimmed reads failed; no 'clear' input", "$wrk/$asm.trimmedReads.err")  if (!defined($inp));

    $cmd  = "$bin/gatekeeperDumpFASTQ -nolibname \\\n";
    $cmd .= "  -G $wrk/$asm.gkpStore \\\n";
    $cmd .= "  -c $inp \\\n";
    $cmd .= "  -o $wrk/$asm.trimmedReads \\\n";
    $cmd .= ">    $wrk/$asm.trimmedReads.err 2>&1";

    if (runCommand($wrk, $cmd)) {
        caFailure("dumping trimmed reads failed", "$wrk/$asm.trimmedReads.err");
    }

    #  Need gatekeeperDumpFASTQ to also write a gkp input file
    #touch("$wrk/$asm.trimmedReads.gkp");

    print STDERR "--\n";
    print STDERR "-- Wrote trimmed reads into '$wrk/$asm.trimmedReads.fastq'\n";
    print STDERR "--\n";

    emitStage($WRK, $asm, "obt-dumpReads");
}
