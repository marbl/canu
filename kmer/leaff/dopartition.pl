#!/bin/perl

# Takes a partition output from leaff, and splits the original
# sequence into those partitions.

my $genomeseq = "h.fasta";
my $genomeout = "h-part";

open(P, "/home/bwalenz/src/genomics/leaff/leaff -F $genomeseq --partition 16 |");
my $numParts = int(<P>);

while (<P>) {
    my @p = split '\s+', $_;

    $id = $p[0];
    shift @p;
    if ($id =~ m/^(\d+)\]/) {
        $id = $1;
    } else {
        print STDERR "Failed to parse id from $id\n";
    }

    my @q;
    foreach my $x (@p) {
        if ($x =~ m/^(\d+)\(\d+\)$/) {
            push @q, $1;
        } else {
            print STDERR "$x\n";
        }
    }

    @q = sort { $a <=> $b } @q;

    open(Z, "| /home/bwalenz/src/genomics/leaff/leaff -F $genomeseq -A - > $genomeout-$id.fasta");
    foreach my $x (@q) {
        print Z "-s $x\n";
    }
    close(Z);

}
close(P);
