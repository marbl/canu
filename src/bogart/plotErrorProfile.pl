#!/usr/local/bin/perl

use strict;

my $nn = shift @ARGV;
my $pn = shift @ARGV;

die "usage: $0 <prefix> <plot-number>\n"  if (!defined($pn));

my $name = "$nn." . substr("00000000$pn", -8) . ".profile";

print "Plotting '$nn' '$pn' - '$name'\n";

my $lastX = 0;
my $lastY = 0;

open(O, "> $name.dat") or die;
open(F, "< $name")     or die;
while (<F>) {
    if (m/^(\d+)\s+(\d+)\s+(\d+.\d+)\s+\+-\s+(\d+.\d+)\s\(\d+\s+overlaps\)/) {
        print O "$lastX\t$lastY\n";
        print O "$1\t$3\n";
        print O "\n";

        print O "$1\t$3\n";
        print O "$2\t$3\n";
        print O "\n";

        $lastX = $2;
        $lastY = $3;
    }
}
close(F);

print O "$lastX\t$lastY\n";
print O "$lastX\t0\n";
print O "\n";

close(O);

open(O, "> $name.gp");
print O "plot '$name.dat' with lines\n";
print O "pause -1\n";
close(O);

system("gnuplot $name.gp");
