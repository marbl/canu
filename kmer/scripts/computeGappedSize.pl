#!/bin/perl

#
#  Reads a genome breakpoint file, outputs the number of bases as N and non-N
#

my $n = 0;
my $b = 0;

my $lb = 0;
my $lp = 0;

#my $brkfile = "/IR/devel/genomics/GENOMES/B34C/B34C.breakpoints";
#my $deffile = "/IR/devel/genomics/GENOMES/B34C/B34C.fasta.deflines";

my $brkfile = "/IR/devel/genomics/GENOMES/MB32A/MB32A.breakpoints";
my $deffile = "/IR/devel/genomics/GENOMES/MB32A/MB32A.fasta.deflines";

open(F, "< $deffile") or die "can't get defflie $deffile\n";
my @defs = <F>;
close(F);
chomp @defs;

my @Ns;
my @Bs;

open(F, "< $brkfile") or die "can't open brkfile $brkfile\n";
while (<F>) {
    my @vals = split '\s+', $_;

    if (($lb eq "N") || ($lb eq "n")) {
        $n += $vals[2] - $lp;
    } else {
        $b += $vals[2] - $lp;
    }

    $lb = $vals[0];
    $lp = $vals[2];

    if ($vals[0] eq ".") {
        push @Ns, $n;
        push @Bs, $b;
            
        $b  = 0;
        $n  = 0;
        $lb = 0;
        $lp = 0;
    }
}
close(F);

if (scalar(@Ns) != scalar(@defs)) {
    print STDERR "Got " . scalar(@Ns) . " Ns and " . scalar(@defs) . " deflines?\n";
    exit(1);
}

while (scalar(@Ns) > 0) {
    my $d = shift @defs;
    my $n = shift @Ns;
    my $b = shift @Bs;

    my $g = $n + $b;

    print STDOUT "$b\t$g\t$d\n";
}
