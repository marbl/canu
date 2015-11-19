
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
 #    src/pipelines/ca3g/Unitig.pm
 #
 #  Modifications by:
 #
 #    Brian P. Walenz beginning on 2015-FEB-27
 #      are Copyright 2015 Battelle National Biodefense Institute, and
 #      are subject to the BSD 3-Clause License
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

package canu::Unitig;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(unitig reportUnitigSizes);

use strict;

use File::Path qw(make_path remove_tree);

use canu::Defaults;
use canu::Execution;
use canu::Gatekeeper;
use canu::HTML;


sub bogart ($$$) {
    my $wrk  = shift @_;
    my $asm  = shift @_;
    my $per  = shift @_;
    my $bin  = getBinDirectory();
    my $cmd;

    #getAllowedResources("", "bat");

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



sub reportUnitigSizes ($$$$) {
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

    my $V = substr("000000" . $version, -3);
    my $N = "$wrk/$asm.tigStore/seqDB.v$V.sizes.txt";

    if (! -e $N) {
        $cmd  = "$bin/tgStoreDump \\\n";
        $cmd .= "  -G $wrk/$asm.gkpStore \\\n";
        $cmd .= "  -T $wrk/$asm.tigStore $version \\\n";
        $cmd .= "  -sizes -s " . getGlobal("genomeSize") . " \\\n";
        $cmd .= "> $N";

        if (runCommand($wrk, $cmd)) {
            caExit("failed to generate unitig sizes", undef);
        }
    }

    open(F, "< $N") or caExit("failed to open '$N' for reading: $!\n", undef);
    while (<F>) {
        $singbases = $1  if (m/lenSingleton\s+sum\s+(\d+)/);
        $singnum   = $1  if (m/lenSingleton\s+num\s+(\d+)/);
        $asmbases  = $1  if (m/lenAssembled\s+sum\s+(\d+)/);
        $asmnum    = $1  if (m/lenAssembled\s+num\s+(\d+)/);

        $asmsizes .= "--   $_"  if (m/lenAssembled\s+(n\d+)\s+siz/);
    }
    close(F);

    print STDERR "-- Found, in version $version, $label:\n";
    print STDERR "--   unitigs:     $asmnum sequences, total length $asmbases bp.\n";
    print STDERR "--   singletons:  $singnum sequences, total length $singbases bp.\n";
    print STDERR "--\n";
    print STDERR "$asmsizes";
    print STDERR "--\n";
}



sub unitig ($$) {
    my $WRK     = shift @_;           #  Root work directory (the -d option to canu)
    my $wrk     = "$WRK/unitigging";  #  Local work directory
    my $asm     = shift @_;

    goto allDone    if (skipStage($wrk, $asm, "unitig") == 1);
    goto allDone    if (-d "$wrk/$asm.tigStore");

    make_path("$wrk/4-unitigger")  if (! -d "$wrk/4-unitigger");

    #  How many reads per partition?  This will change - it'll move to be after unitigs are constructed.

    my $perPart = int(getNumberOfReadsInStore($wrk, $asm) / getGlobal("cnsPartitions"));
    my $minPart = getGlobal("cnsPartitionMin");

    $perPart = ($perPart < $minPart) ? ($perPart) : ($minPart);

    my $cmd;

    if      (getGlobal("unitigger") eq "bogart") {
        $cmd = bogart($wrk, $asm, $perPart);

    } else {
        caFailure("unknown unitigger '" . getGlobal("unitigger") . "'", undef);
    }

    stopBefore("unitig", $cmd);

    if (runCommand("$wrk/4-unitigger", $cmd)) {
        caExit("failed to unitig", "$wrk/4-unitigger/unitigger.err");
    }

  finishStage:
    emitStage($WRK, $asm, "unitig");
    buildHTML($WRK, $asm, "utg");
    stopAfter("unitig");

  allDone:
    reportUnitigSizes($wrk, $asm, 1, "after unitig construction");
}
