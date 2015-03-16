package ca3g::Output;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(outputLayout outputSequence);

use strict;

use ca3g::Defaults;
use ca3g::Execution;

sub outputLayout ($$) {
    my $wrk    = shift @_;
    my $asm    = shift @_;
    my $bin    = getBinDirectory();
    my $cmd;

    $cmd  = "$bin/tigStore \\\n";
    $cmd .= "  -g $wrk/$asm.gkpStore \\\n";
    $cmd .= "  -t $wrk/$asm.tigStore 2 \\\n";
    $cmd .= "  -U -d layout \\\n";
    $cmd .= ">  $wrk/$asm.layout \\\n";
    $cmd .= "2> $wrk/$asm.layout.err\n";

    if (runCommand($wrk, $cmd)) {
        caFailure("failed to output layouts", "$wrk/$asm.layout.err");
    }
}


sub outputSequence ($$) {
    my $wrk    = shift @_;
    my $asm    = shift @_;
    my $bin    = getBinDirectory();
    my $cmd;

    $cmd  = "$bin/tigStore \\\n";
    $cmd .= "  -g $wrk/$asm.gkpStore \\\n";
    $cmd .= "  -t $wrk/$asm.tigStore 2 \\\n";
    $cmd .= "  -U -d consensus \\\n";
    $cmd .= ">  $wrk/$asm.consensus.fasta \\\n";
    $cmd .= "2> $wrk/$asm.consensus.fasta.err\n";

    if (runCommand("$wrk/0-overlaptrim", $cmd)) {
        caFailure("failed to output consensus", "$wrk/$asm.consensus.fasta.err");
    }
}
