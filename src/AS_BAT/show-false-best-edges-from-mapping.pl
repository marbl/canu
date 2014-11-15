#!/usr/bin/perl

use strict;

my $asm = shift @ARGV;

#  Reads:
#    $asm//best.edges
#    filtered-overlaps.true.ova
#    filtered-overlaps.false.ova
#
#  Reports which best.edges are false

my %true;
my %false;

open(F, "< filtered-overlaps.true.ova") or die;
while (<F>) {
    s/^\s+//;
    s/\s+$//;

    my @v = split '\s+', $_;

    $true{"$v[0]-$v[1]"}++;
    $true{"$v[1]-$v[0]"}++;
}
close(F);


open(F, "< filtered-overlaps.false.ova") or die;
while (<F>) {
    s/^\s+//;
    s/\s+$//;

    my @v = split '\s+', $_;

    $false{"$v[0]-$v[1]"}++;
    $false{"$v[1]-$v[0]"}++;
}
close(F);


    my $true5  = 0;
    my $true3  = 0;
    my $false5 = 0;
    my $false3 = 0;
    my $novel5 = 0;
    my $novel3 = 0;

open(F, "< $asm/best.edges") or die;
while (<F>) {
    s/^\s+//;
    s/\s+$//;

    my @v = split '\s+', $_;

    my $p5 = "$v[0]-$v[2]";
    my $p3 = "$v[0]-$v[4]";

    if      (exists($true{$p5})) {
        $true5++;
    } elsif (exists($false{$p5})) {
        print STDERR "FALSE 5 $_\n";
        $false5++;
    } else {
        print STDERR "NOVEL 5 $_\n";
        $novel5++;
    }

    if      (exists($true{$p3})) {
        $true3++;
    } elsif (exists($false{$p3})) {
        print STDERR "FALSE 3 $_\n";
        $false3++;
    } else {
        print STDERR "NOVEL 3 $_\n";
        $novel3++;
    }
}
close(F);

print "true5   $true5\n";
print "true3   $true3\n";
print "false5  $false5\n";
print "false3  $false3\n";
print "novel5  $novel5\n";
print "novel3  $novel3\n";
