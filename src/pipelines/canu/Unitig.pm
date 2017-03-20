
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
 #    Brian P. Walenz from 2015-FEB-27 to 2015-AUG-26
 #      are Copyright 2015 Battelle National Biodefense Institute, and
 #      are subject to the BSD 3-Clause License
 #
 #    Brian P. Walenz beginning on 2015-NOV-02
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #    Sergey Koren beginning on 2015-NOV-27
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

package canu::Unitig;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(reportUnitigSizes unitig unitigCheck);

use strict;

use File::Path 2.08 qw(make_path remove_tree);
use POSIX qw(ceil);

use canu::Defaults;
use canu::Execution;
use canu::Configure;    #  For displayGenomeSize
use canu::Gatekeeper;
use canu::Report;
use canu::HTML;
use canu::Meryl;
use canu::Grid_Cloud;




sub reportUnitigSizes ($$$) {
    my $asm       = shift @_;
    my $version   = shift @_;
    my $label     = shift @_;

    my $bin       = getBinDirectory();
    my $cmd       = "";

    my $ctgNum    = 0;
    my $ctgBases  = 0;
    my $ctgSizes  = "";

    my $usmNum   = 0;
    my $usmBases = 0;

    my $bubNum   = 0;
    my $bubBases = 0;

    my $rptNum    = 0;
    my $rptBases  = 0;

    my $gs        = getGlobal("genomeSize");

    my $V = substr("000000" . $version, -3);
    my $N = "$asm.ctgStore/seqDB.v$V.sizes.txt";

    fetchFile("unitigging/$N");

    if (! -e "unitigging/$N") {
        fetchStore("unitigging/$asm.gkpStore");

        fetchFile("unitigging/$asm.ctgStore/seqDB.v$V.dat");
        fetchFile("unitigging/$asm.ctgStore/seqDB.v$V.tig");

        $cmd  = "$bin/tgStoreDump \\\n";                     #  Duplicated at the end of unitigger.sh
        $cmd .= "  -G ./$asm.gkpStore \\\n";
        $cmd .= "  -T ./$asm.ctgStore $version \\\n";
        $cmd .= "  -sizes -s " . getGlobal("genomeSize") . " \\\n";
        $cmd .= "> ./$N";

        if (runCommand("unitigging", $cmd)) {
            caExit("failed to generate unitig sizes", undef);
        }

        stashFile("unitigging/$N");
    }

    $ctgSizes .= "--            NG (bp)  LG (contigs)    sum (bp)\n";
    $ctgSizes .= "--         ----------  ------------  ----------\n";

    open(F, "< unitigging/$N") or caExit("failed to open 'unitigging/$N' for reading: $!", undef);
    while (<F>) {
        $rptBases  = $1  if (m/lenSuggestRepeat\s+sum\s+(\d+)/);
        $rptNum    = $1  if (m/lenSuggestRepeat\s+num\s+(\d+)/);

        $usmBases  = $1  if (m/lenUnassembled\s+sum\s+(\d+)/);
        $usmNum    = $1  if (m/lenUnassembled\s+num\s+(\d+)/);

        $bubBases  = $1  if (m/lenBubble\s+sum\s+(\d+)/);
        $bubNum    = $1  if (m/lenBubble\s+num\s+(\d+)/);

        $ctgBases  = $1  if (m/lenContig\s+sum\s+(\d+)/);
        $ctgNum    = $1  if (m/lenContig\s+num\s+(\d+)/);

        if (m/lenContig\s+ng(\d+)\s+(\d+)\s+bp\s+lg(\d+)\s+(\d+)\s+sum\s+(\d+)\s+bp$/) {
            my $n  = substr("          $1", -3);
            my $ng = substr("          $2", -10);
            my $lg = substr("          $4", -6);
            my $ss = substr("          $5", -10);

            $ctgSizes .= "--    $n  $ng        $lg  $ss\n";
        }

        #  User could have changed genome size since this report was dumped.
        if (m/lenContig.*genomeSize\s+(\d+)/) {
            $gs = $1;
        }
    }
    close(F);

    my $report;

    $report .= "-- Found, in version $version, $label:\n";
    $report .= "--   contigs:      $ctgNum sequences, total length $ctgBases bp (including $rptNum repeats of total length $rptBases bp).\n";
    $report .= "--   bubbles:      $bubNum sequences, total length $bubBases bp.\n";
    $report .= "--   unassembled:  $usmNum sequences, total length $usmBases bp.\n";
    $report .= "--\n";
    $report .= "-- Contig sizes based on genome size ", displayGenomeSize($gs), "bp:\n";
    $report .= "--\n";
    $report .= $ctgSizes;
    $report .= "--\n";

    #  Hmmm.  Report wants a small tag, but $label is a bit verbose.

    addToReport("contigs",   $report)    if ($label eq "after unitig construction");
    addToReport("consensus", $report)    if ($label eq "after consensus generation");
}



sub unitig ($) {
    my $asm     = shift @_;
    my $path    = "unitigging/4-unitigger";

    goto allDone    if (skipStage($asm, "unitig") == 1);
    goto allDone    if (fileExists("unitigging/4-unitigger/unitigger.sh"));
    goto allDone    if (fileExists("unitigging/$asm.utgStore/seqDB.v001.tig") &&
                        fileExists("unitigging/$asm.ctgStore/seqDB.v001.tig"));

    make_path($path)  if (! -d $path);

    #  How many reads per partition?  This will change - it'll move to be after unitigs are constructed.

    my $overlapLength = getGlobal("minOverlapLength");

    #  Dump a script to run the unitigger.

    open(F, "> $path/unitigger.sh") or caExit("can't open '$path/unitigger.sh' for writing: $!\n", undef);

    print F "#!" . getGlobal("shell") . "\n";
    print F "\n";
    print F getBinDirectoryShellCode();
    print F "\n";
    print F setWorkDirectoryShellCode($path);
    print F "\n";
    print F fetchStoreShellCode("unitigging/$asm.gkpStore", $path, "");
    print F fetchStoreShellCode("unitigging/$asm.ovlStore", $path, "");
    print F "\n";
    print F fetchFileShellCode("unitigging/$asm.ovlStore", "evalues", "");
    print F "\n";
    print F getJobIDShellCode();
    print F "\n";
    print F "if [ -e unitigging/$asm.ctgStore/seqDB.v001.tig -a -e unitigging/$asm.utgStore/seqDB.v001.tig ] ; then\n";
    print F "  exit 0\n";
    print F "fi\n";
    print F "\n";

    if      (getGlobal("unitigger") eq "bogart") {
        print F "\$bin/bogart \\\n";
        print F " -G ../$asm.gkpStore \\\n";
        print F " -O ../$asm.ovlStore \\\n";
        print F " -o ./$asm \\\n";
        print F " -gs "             . getGlobal("genomeSize")         . " \\\n";
        print F " -eg "             . getGlobal("utgErrorRate")       . " \\\n";
        print F " -eM "             . getGlobal("utgErrorRate")       . " \\\n";
        print F " -el "             . $overlapLength                  . " \\\n";
        print F " -dg "             . getGlobal("utgGraphDeviation")  . " \\\n";
        print F " -db "             . getGlobal("utgGraphDeviation")  . " \\\n";
        print F " -dr "             . getGlobal("utgRepeatDeviation") . " \\\n";
        print F " -ca "             . getGlobal("utgRepeatConfusedBP"). " \\\n";
        print F " -cp "             . "200"                           . " \\\n";
        print F " -threads "        . getGlobal("batThreads")         . " \\\n"   if (defined(getGlobal("batThreads")));
        print F " -M "              . getGlobal("batMemory")          . " \\\n"   if (defined(getGlobal("batMemory")));
        print F " -unassembled "    . getGlobal("contigFilter")       . " \\\n"   if (defined(getGlobal("contigFilter")));
        print F " "                 . getGlobal("batOptions")         . " \\\n"   if (defined(getGlobal("batOptions")));
        print F " > ./unitigger.err 2>&1 \\\n";
        print F "&& \\\n";
        print F "mv ./$asm.ctgStore ../$asm.ctgStore \\\n";
        print F "&& \\\n";
        print F "mv ./$asm.utgStore ../$asm.utgStore\n";
    } else {
        caFailure("unknown unitigger '" . getGlobal("unitigger") . "'", undef);
    }

    print F "\n";
    print F stashFileShellCode("unitigging/4-unitigger", "$asm.unitigs.gfa", "");
    print F stashFileShellCode("unitigging/4-unitigger", "$asm.contigs.gfa", "");
    print F "\n";
    print F stashFileShellCode("unitigging/$asm.ctgStore", "seqDB.v001.dat", "");
    print F stashFileShellCode("unitigging/$asm.ctgStore", "seqDB.v001.tig", "");
    print F "\n";
    print F stashFileShellCode("unitigging/$asm.utgStore", "seqDB.v001.dat", "");
    print F stashFileShellCode("unitigging/$asm.utgStore", "seqDB.v001.tig", "");
    print F "\n";
    print F "\$bin/tgStoreDump \\\n";                    #  Duplicated in reportUnitigSizes()
    print F "  -G ../$asm.gkpStore \\\n";                 #  Done here so we don't need another
    print F "  -T ../$asm.ctgStore 1 \\\n";               #  pull of gkpStore and ctgStore
    print F "  -sizes -s " . getGlobal("genomeSize") . " \\\n";
    print F "> ../$asm.ctgStore/seqDB.v001.sizes.txt";
    print F "\n";
    print F stashFileShellCode("unitigging/$asm.ctgStore", "seqDB.v001.sizes.txt", "");
    print F "\n";
    print F "exit 0\n";

    close(F);

    stashFile("$path/unitigger.sh");

  finishStage:
    emitStage($asm, "unitig");
    buildHTML($asm, "utg");

  allDone:
}




sub unitigCheck ($) {
    my $asm     = shift @_;
    my $attempt = getGlobal("canuIteration");
    my $path    = "unitigging/4-unitigger";

    goto allDone      if (skipStage($asm, "unitigCheck", $attempt) == 1);
    goto allDone      if (fileExists("$path/unitigger.success"));
    goto finishStage  if (fileExists("unitigging/$asm.utgStore/seqDB.v001.tig") &&
                          fileExists("unitigging/$asm.ctgStore/seqDB.v001.tig"));

    fetchFile("$path/unitigger.sh");

    #  Since there is only one job, if we get here, we're not done.  Any other 'check' function
    #  shows how to process multiple jobs.  This only checks for the existence of the final outputs.
    #  (meryl is the same)

    #  If not the first attempt, report the jobs that failed, and that we're recomputing.

    if ($attempt > 1) {
        print STDERR "--\n";
        print STDERR "-- Unitigger failed.\n";
        print STDERR "--\n";
    }

    #  If too many attempts, give up.

    if ($attempt > getGlobal("canuIterationMax")) {
        caExit("failed to generate unitigs.  Made " . ($attempt-1) . " attempts, jobs still failed", undef);
    }

    #  Otherwise, run some jobs.

    emitStage($asm, "unitigCheck", $attempt);
    buildHTML($asm, "utg");

    submitOrRunParallelJob($asm, "bat", $path, "unitigger", (1));
    return;

  finishStage:
    print STDERR "-- Unitigger finished successfully.\n";

    make_path($path);   #  With object storage, we might not have this directory!

    open(F, "> $path/unitigger.success") or caExit("can't open '$path/unitigger.success' for writing: $!", undef);
    close(F);

    stashFile("$path/unitigger.success");

    reportUnitigSizes($asm, 1, "after unitig construction");

    emitStage($asm, "unitigCheck");
    buildHTML($asm, "utg");

  allDone:
    stopAfter("unitig");
}
