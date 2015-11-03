
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
 #    Brian P. Walenz beginning on 2015-MAR-16
 #      are Copyright 2015 Battelle National Biodefense Institute, and
 #      are subject to the BSD 3-Clause License
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

package canu::Output;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(outputLayout outputGraph outputSequence);

use strict;

use canu::Defaults;
use canu::Execution;


sub outputLayout ($$) {
    my $wrk    = shift @_;
    my $asm    = shift @_;
    my $bin    = getBinDirectory();
    my $cmd;

    goto allDone   if (skipStage($wrk, $asm, "outputLayout") == 1);
    goto allDone   if (-e "$wrk/$asm.layout");

    if (-e "$wrk/$asm.tigStore/seqDB.v002.tig") {
        $cmd  = "$bin/tgStoreDump \\\n";
        $cmd .= "  -G $wrk/$asm.gkpStore \\\n";
        $cmd .= "  -T $wrk/$asm.tigStore 2 \\\n";
        $cmd .= "  -o $wrk/$asm \\\n";
        $cmd .= "  -layout \\\n";
        $cmd .= "> $wrk/$asm.layout.err 2>&1\n";

        if (runCommand($wrk, $cmd)) {
            caExit("failed to output layouts", "$wrk/$asm.layout.err");
        }

    } else {
        open(O, "> $wrk/$asm.layout") or caExit("can't open '$wrk/$asm.layout' for writing: $!", undef);
        open(F, "< $wrk/5-consensus/cnsjob.files") or caExit("can't open '$wrk/5-consensus/cnsjob.files' for reading: $!", undef);
        while (<F>) {
            my $prefix = $1  if (m/^(.*).cns/);

            if (-e "$prefix.layout") {
                open(L, "< $prefix.layout") or caExit("can't open '$prefix.layout' for reading: $!", undef);
                while (<L>) {
                    next  if (m/^cns\s/);
                    next  if (m/^qlt\s/);

                    print O $_;
                }
                close(L);
            } else {
                caExit("can't open '$prefix.layout' for reading: $!", undef);
            }
        }
        close(F);
        close(O);
    }

  finishStage:
    emitStage($wrk, $asm, "outputLayout");

  allDone:
    print STDERR "--  Unitig layouts saved in '$wrk/$asm.layout'.\n";
}




sub outputGraph ($$) {
    my $wrk    = shift @_;
    my $asm    = shift @_;
    my $bin    = getBinDirectory();
    my $cmd;

    goto allDone   if (skipStage($wrk, $asm, "outputGraph") == 1);
    goto allDone     if (-e "$wrk/$asm.graph");

    #
    #  Stuff here.
    #

    if (runCommand($wrk, $cmd)) {
        caExit("failed to output consensus", "$wrk/$asm.graph.err");
    }


  finishStage:
    emitStage($wrk, $asm, "outputGraph");

  allDone:
    print STDERR "-- Unitig graph saved in (nowhere, yet).\n";
}




sub outputSequence ($$) {
    my $wrk    = shift @_;
    my $asm    = shift @_;
    my $bin    = getBinDirectory();
    my $cmd;

    goto allDone   if (skipStage($wrk, $asm, "outputSequence") == 1);
    goto allDone   if (-e "$wrk/$asm.fastq");

    if (-e "$wrk/$asm.tigStore/seqDB.v002.tig") {
        $cmd  = "$bin/tgStoreDump \\\n";
        $cmd .= "  -G $wrk/$asm.gkpStore \\\n";
        $cmd .= "  -T $wrk/$asm.tigStore 2 \\\n";
        $cmd .= "  -o $wrk/$asm \\\n";
        $cmd .= "  -consensus -fasta \\\n";
        $cmd .= "2> $wrk/$asm.consensus.err\n";

        if (runCommand($wrk, $cmd)) {
            caExit("failed to output consensus", "$wrk/$asm.consensus.fastq.err");
        }

    } else {
        open(O, "> $wrk/$asm.fastq") or caExit("can't open '$wrk/$asm.fastq' for writing: $!", undef);
        open(F, "< $wrk/5-consensus/cnsjob.files") or caExit("can't open '$wrk/5-consensus/cnsjob.files' for reading: $!", undef);
        while (<F>) {
            my $prefix = $1  if (m/^(.*).cns/);

            if (-e "$prefix.fastq") {
                print STDERR "Copying sequencess from $prefix.fastq\n";

                open(L, "< $prefix.fastq") or caExit("can't open '$prefix.fastq' for reading: $!", undef);
                while (<L>) {
                    print O $_;
                }
                close(L);
            } else {
                caExit("can't open '$prefix.fastq' for reading: $!", undef);
            }
        }
        close(F);
        close(O);
    }

  finishStage:
    emitStage($wrk, $asm, "outputSequence");

  allDone:
    print STDERR "--  Unitig sequences saved in '$wrk/$asm.consensus.fastq'.\n";
}
