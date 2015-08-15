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

    goto stopAfter   if (skipStage($wrk, $asm, "outputLayout") == 1);
    goto allDone     if (-e "$wrk/$asm.layout");

    if (-e "$wrk/$asm.tigStore/seqDB.v002.tig") {
        $cmd  = "$bin/tgStoreDump \\\n";
        $cmd .= "  -G $wrk/$asm.gkpStore \\\n";
        $cmd .= "  -T $wrk/$asm.tigStore 2 \\\n";
        $cmd .= "  -U -d layout \\\n";
        $cmd .= ">  $wrk/$asm.layout \\\n";
        $cmd .= "2> $wrk/$asm.layout.err\n";

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

  stopBefore:
    #stopBefore("outputLayout", $cmd);

    if (runCommand($wrk, $cmd)) {
        caExit("failed to output layouts", "$wrk/$asm.layout.err");
    }

  allDone:
    emitStage($wrk, $asm, "outputLayout");
  stopAfter:

    print STDERR "--\n";
    print STDERR "-- wrote unitig layouts to '$wrk/$asm.layout'.\n";
    print STDERR "--\n";
}




sub outputGraph ($$) {
    my $wrk    = shift @_;
    my $asm    = shift @_;
    my $bin    = getBinDirectory();
    my $cmd;

    goto stopAfter   if (skipStage($wrk, $asm, "outputGraph") == 1);
    goto allDone     if (-e "$wrk/$asm.graph");

    #
    #  Stuff here.
    #

  stopBefore:
    #stopBefore("outputSequence", $cmd);

    if (runCommand($wrk, $cmd)) {
        caExit("failed to output consensus", "$wrk/$asm.graph.err");
    }

  allDone:
    emitStage($wrk, $asm, "outputGraph");
  stopAfter:

    print STDERR "--\n";
    print STDERR "-- wrote unitig graph to (nowhere, yet).\n";
    print STDERR "--\n";
}




sub outputSequence ($$) {
    my $wrk    = shift @_;
    my $asm    = shift @_;
    my $bin    = getBinDirectory();
    my $cmd;

    goto stopAfter   if (skipStage($wrk, $asm, "outputSequence") == 1);
    goto allDone     if (-e "$wrk/$asm.fastq");

    if (-e "$wrk/$asm.tigStore/seqDB.v002.tig") {
        $cmd  = "$bin/tgStoreDump \\\n";
        $cmd .= "  -G $wrk/$asm.gkpStore \\\n";
        $cmd .= "  -T $wrk/$asm.tigStore 2 \\\n";
        $cmd .= "  -U -d consensus \\\n";
        $cmd .= ">  $wrk/$asm.consensus.fasta \\\n";
        $cmd .= "2> $wrk/$asm.consensus.fasta.err\n";

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

  stopBefore:
    #stopBefore("outputSequence", $cmd);

    if (runCommand($wrk, $cmd)) {
        caExit("failed to output consensus", "$wrk/$asm.consensus.fasta.err");
    }

  allDone:
    emitStage($wrk, $asm, "outputSequence");
  stopAfter:

    print STDERR "--\n";
    print STDERR "-- wrote unitig sequence to $wrk/$asm.consensus.f''.\n";
    print STDERR "--\n";
}
