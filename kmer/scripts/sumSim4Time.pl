#!/usr/local/bin/perl

$| = 1;

my $p = "0000";
my $c = 0;
my $s = 0;
my $u = 0;
my $C = 0;
my $S = 0;
my $U = 0;

open(F, "ls *stats |");

while (<F>) {
    chomp;

    open(Z, "< $_");
    while (<Z>) {
        if (m/^clockTime:\s+(\d+.\d+)$/) {
            $c += $1;
        }
        if (m/^userTime:\s+(\d+.\d+)$/) {
            $u += $1;
        }
        if (m/^systemTime:\s+(\d+.\d+)$/) {
            $s += $1;
        }
    }
    close(Z);
}

close(F);

$C = int(1000.0*$c/3600) / 1000.0;
$S = int(1000.0*$s/3600) / 1000.0;
$U = int(1000.0*$u/3600) / 1000.0;

$c = int(1000.0*$c) / 1000.0;
$s = int(1000.0*$s) / 1000.0;
$u = int(1000.0*$u) / 1000.0;

print "clk,sys,usr: sec: $c, $s, $u  hrs: $C, $S, $U\n";
