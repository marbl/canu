
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
use canu::HTML;
use canu::Meryl;




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

    if (! -e "unitigging/$N") {
        $cmd  = "$bin/tgStoreDump \\\n";
        $cmd .= "  -G ./$asm.gkpStore \\\n";
        $cmd .= "  -T ./$asm.ctgStore $version \\\n";
        $cmd .= "  -sizes -s " . getGlobal("genomeSize") . " \\\n";
        $cmd .= "> ./$N";

        if (runCommand("unitigging", $cmd)) {
            caExit("failed to generate unitig sizes", undef);
        }
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

    print STDERR "-- Found, in version $version, $label:\n";
    print STDERR "--   contigs:      $ctgNum sequences, total length $ctgBases bp (including $rptNum repeats of total length $rptBases bp).\n";
    print STDERR "--   bubbles:      $bubNum sequences, total length $bubBases bp.\n";
    print STDERR "--   unassembled:  $usmNum sequences, total length $usmBases bp.\n";
    print STDERR "--\n";
    print STDERR "-- Contig sizes based on genome size ", displayGenomeSize($gs), "bp:\n";
    print STDERR "--\n";
    print STDERR "$ctgSizes";
    print STDERR "--\n";
}



sub unitig ($) {
    my $asm     = shift @_;

    goto allDone    if (skipStage($asm, "unitig") == 1);
    goto allDone    if ((-d "unitigging/$asm.ctgStore") && (-d "unitigging/$asm.utgStore"));

    make_path("unitigging/4-unitigger")  if (! -d "unitigging/4-unitigger");

    #  How many reads per partition?  This will change - it'll move to be after unitigs are constructed.

    my $overlapLength = getGlobal("minOverlapLength");

    #  Dump a script to run the unitigger.

    open(F, "> unitigging/4-unitigger/unitigger.sh") or caExit("can't open 'unitigging/4-unitigger/unitigger.sh' for writing: $!\n", undef);

    print F "#!" . getGlobal("shell") . "\n";
    print F "\n";
    print F "if [ -e unitigging/$asm.ctgStore/seqDB.v001.tig -a -e unitigging/$asm.utgStore/seqDB.v001.tig ] ; then\n";
    print F "  exit 0\n";
    print F "fi\n";
    print F "\n";
    print F getBinDirectoryShellCode();
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
    print F "exit 0\n";

    close(F);

  finishStage:
    emitStage($asm, "unitig");
    buildHTML($asm, "utg");

  allDone:
}




sub unitigCheck ($) {
    my $asm     = shift @_;
    my $attempt = getGlobal("canuIteration");
    my $path    = "unitigging/4-unitigger";

    goto allDone  if (skipStage($asm, "unitigCheck", $attempt) == 1);
    goto allDone  if ((-e "unitigging/$asm.ctgStore/seqDB.v001.tig") && (-e "unitigging/$asm.utgStore/seqDB.v001.tig"));

    #  Since there is only one job, if we get here, we're not done.  Any other 'check' function
    #  shows how to process multiple jobs.  This only checks for the existence of the final outputs.

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

    print STDERR "-- Unitigger attempt $attempt begins.\n";

    emitStage($asm, "unitigCheck", $attempt);
    buildHTML($asm, "utg");

    submitOrRunParallelJob($asm, "bat", $path, "unitigger", (1));
    return;

    #  If onGrid, the submitOrRun() has submitted parallel jobs to the grid, and resubmitted the
    #  executive.  #  The parallel version never gets here

  finishStage:
    print STDERR "-- Unitigger finished successfully.\n";

    reportUnitigSizes($asm, 1, "after unitig construction");

    setGlobal("canuIteration", 1);
    emitStage($asm, "unitigCheck");
    buildHTML($asm, "utg");

  allDone:
    stopAfter("unitig");
}
