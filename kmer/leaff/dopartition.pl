#!/bin/perl

# Takes a partition output from leaff, and splits the original
# sequence into those partitions.

if (scalar(@ARGV) != 3) {
    die "usage: $0 genome.fasta partition prefix\n";
}

my $genomeseq = shift @ARGV;
my $partition = shift @ARGV;
my $genomeout = shift @ARGV;

my $leaff = "/bioinfo/assembly/walenz/src/genomics/osf1/bin/leaff";

open(P, "$leaff -F $genomeseq --partition $partition |");
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

    open(Z, "| $leaff -F $genomeseq -A - > $genomeout-$id.fasta");
    foreach my $x (@q) {
        print Z "-s $x\n";
    }
    close(Z);

}
close(P);
