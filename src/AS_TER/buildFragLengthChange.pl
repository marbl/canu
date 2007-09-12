#!/usr/bin/perl

use strict;

my $pfx = "gambiae";
my $dir = "/bioinfo/assembly/walenz/anopheles/g4";
my $wgs = "/bioinfo/assembly/walenz/wgs403/Linux-amd64/bin";

$dir = "/bioinfo/assembly/walenz/anopheles/g4/OLD";

my %clearLength;
my %lengths;

print STDERR "Reading $pfx.gkpStore.\n";
open(F, "$wgs/gatekeeper -dumpfragments -tabular -clear ecr2 $dir/$pfx.gkpStore |") or die;
while (<F>) {
    my @v = split '\s+', $_;
    $clearLength{$v[0]} = $v[12] - $v[11];
}
close(F);

print STDERR "Reading 9-terminator/$pfx.posmap.frgscf.\n";
open(F, "< $dir/9-terminator/$pfx.posmap.frgctg") or die;
while (<F>) {
    my @v = split '\s+', $_;
    my $l = $v[3] - $v[2];
    if (defined($clearLength{$v[0]})) {
        $l = $clearLength{$v[0]} - $l;
        $lengths{$l}++;

    } else {
        print "Got frag $v[0] in assembly not in clear?\n";
    }
}
close(F);

print STDERR "Output.\n";
foreach my $l (sort {$a <=> $b} keys %lengths) {
    print "$l\t$lengths{$l}\n";
}

#  then gnuplot
#  set logscale y
#  plot [-200:200][1:] "frag-clearsize-asmsize.out" using 1:2 with lines
