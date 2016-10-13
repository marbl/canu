
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
 #    src/pipelines/ca3g/Output.pm
 #
 #  Modifications by:
 #
 #    Brian P. Walenz from 2015-MAR-16 to 2015-AUG-25
 #      are Copyright 2015 Battelle National Biodefense Institute, and
 #      are subject to the BSD 3-Clause License
 #
 #    Brian P. Walenz beginning on 2015-NOV-02
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #    Sergey Koren beginning on 2015-DEC-02
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

package canu::Output;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(outputLayout outputGraph outputSequence outputSummary);

use strict;

use File::Copy;

use canu::Defaults;
use canu::Execution;
use canu::HTML;


#  In this module, the inputs are read from $wrk (the local work directory), and the output are
#  written to $WRK (the root work directory).  This is a change from CA8 where outputs
#  were written to the 9-terminator directory.


sub outputLayout ($$) {
    my $WRK     = shift @_;           #  Root work directory (the -d option to canu)
    my $wrk     = "$WRK/unitigging";  #  Local work directory
    my $asm     = shift @_;
    my $bin     = getBinDirectory();
    my $cmd;

    goto allDone   if (skipStage($WRK, $asm, "outputLayout") == 1);
    goto allDone   if ((-e "$WRK/$asm.contigs.layout") && (-e "$WRK/$asm.unitigs.layout"));

    if (! -e "$WRK/$asm.contigs.layout") {
        $cmd  = "$bin/tgStoreDump \\\n";
        $cmd .= "  -G $wrk/$asm.gkpStore \\\n";
        $cmd .= "  -T $wrk/$asm.ctgStore 2 \\\n";
        $cmd .= "  -o $WRK/$asm.contigs \\\n";
        $cmd .= "  -layout \\\n";
        $cmd .= "> $WRK/$asm.contigs.layout.err 2>&1";

        if (runCommand($wrk, $cmd)) {
            caExit("failed to output contig layouts", "$WRK/$asm.contigs.layout.err");
        }

        unlink "$WRK/$asm.contigs.layout.err";
    }

    if (! -e "$WRK/$asm.unitigs.layout") {
        $cmd  = "$bin/tgStoreDump \\\n";
        $cmd .= "  -G $wrk/$asm.gkpStore \\\n";
        $cmd .= "  -T $wrk/$asm.utgStore 2 \\\n";
        $cmd .= "  -o $WRK/$asm.unitigs \\\n";
        $cmd .= "  -layout \\\n";
        $cmd .= "> $WRK/$asm.unitigs.layout.err 2>&1";

        if (runCommand($wrk, $cmd)) {
            caExit("failed to output unitig layouts", "$WRK/$asm.unitigs.layout.err");
        }

        unlink "$WRK/$asm.unitigs.layout.err";
    }

  finishStage:
    emitStage($WRK, $asm, "outputLayout");
    buildHTML($WRK, $asm, "utg");

  allDone:
    print STDERR "-- Contig layouts saved in '$WRK/$asm.contigs.layout'.\n";
    print STDERR "-- Unitig layouts saved in '$WRK/$asm.unitigs.layout'.\n";
}




sub outputGraph ($$) {
    my $WRK     = shift @_;           #  Root work directory (the -d option to canu)
    my $wrk     = "$WRK/unitigging";  #  Local work directory
    my $asm     = shift @_;
    my $bin     = getBinDirectory();
    my $cmd;

    my $type = "fasta";  #  Should probably be an option.

    goto allDone   if (skipStage($WRK, $asm, "outputGraph") == 1);
    goto allDone   if (-e "$WRK/$asm.unitigs.gfa");

    if ((! -e "$WRK/$asm.unitigs.gfa") &&
        (  -e "$wrk/4-unitigger/$asm.unitigs.gfa")) {
        copy("$wrk/4-unitigger/$asm.unitigs.gfa", "$WRK/$asm.unitigs.gfa");
    }

  finishStage:
    emitStage($WRK, $asm, "outputGraph");
    buildHTML($WRK, $asm, "utg");

  allDone:
    print STDERR "-- Unitig graph saved in '$WRK/$asm.unitigs.gfa'.\n";
}




sub outputSequence ($$) {
    my $WRK     = shift @_;           #  Root work directory (the -d option to canu)
    my $wrk     = "$WRK/unitigging";  #  Local work directory
    my $asm     = shift @_;
    my $bin     = getBinDirectory();
    my $cmd;

    my $type = "fasta";  #  Should probably be an option.

    goto allDone   if (skipStage($WRK, $asm, "outputSequence") == 1);
    goto allDone   if (-e "$WRK/$asm.contigs.$type");

    foreach my $tt ("unassembled", "bubbles", "contigs") {
        if (! -e "$WRK/$asm.$tt.$type") {
            $cmd  = "$bin/tgStoreDump \\\n";
            $cmd .= "  -G $wrk/$asm.gkpStore \\\n";
            $cmd .= "  -T $wrk/$asm.ctgStore 2 \\\n";
            $cmd .= "  -consensus -$type \\\n";
            $cmd .= "  -$tt \\\n";
            $cmd .= "> $WRK/$asm.$tt.$type\n";
            $cmd .= "2> $WRK/$asm.$tt.err";

            if (runCommand($WRK, $cmd)) {
                caExit("failed to output $tt consensus sequences", "$WRK/$asm.$tt.err");
            }

            unlink "$WRK/$asm.$tt.err";
        }
    }

    if (! -e "$WRK/$asm.unitigs.$type") {
        $cmd  = "$bin/tgStoreDump \\\n";
        $cmd .= "  -G $wrk/$asm.gkpStore \\\n";
        $cmd .= "  -T $wrk/$asm.ctgStore 2 \\\n";
        $cmd .= "  -consensus -$type \\\n";
        $cmd .= "  -contigs \\\n";
        $cmd .= "> $WRK/$asm.unitigs.$type\n";
        $cmd .= "2> $WRK/$asm.unitigs.err";

        if (runCommand($WRK, $cmd)) {
            caExit("failed to output unitig consensus sequences", "$WRK/$asm.unitigs.err");
        }

        unlink "$WRK/$asm.unitigs.err";
    }

  finishStage:
    emitStage($WRK, $asm, "outputSequence");
    buildHTML($WRK, $asm, "utg");

  allDone:
    print STDERR "-- Sequences saved:\n";
    print STDERR "--   Contigs       -> '$WRK/$asm.contig.$type'\n";
    print STDERR "--   Bubbles       -> '$WRK/$asm.bubble.$type'  (DEPRECATED)\n";
    print STDERR "--   Unassembled   -> '$WRK/$asm.unassembled.$type'\n";
    print STDERR "--   Unitigs       -> '$WRK/$asm.unitigs.$type'\n";
}



sub outputSummary ($$) {
    my $WRK     = shift @_;           #  Root work directory (the -d option to canu)
    my $wrk     = "$WRK/unitigging";  #  Local work directory
    my $asm     = shift @_;
    my $bin     = getBinDirectory();
    my $cmd;

    goto allDone   if (skipStage($WRK, $asm, "outputSummary") == 1);
    goto allDone   if (-e "$WRK/unitiggging.html");

  finishStage:
    emitStage($WRK, $asm, "outputSummary");
    buildHTML($WRK, $asm, "utg");

  allDone:
    print STDERR "-- Summary saved in '$WRK/unitigging.html'.\n";
}
