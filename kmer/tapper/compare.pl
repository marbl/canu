#!/usr/bin/perl
use strict;

#  Hack to compare tapperconvert to corona output.
#
#  corona:
#  926_28_374      T0233011320231302110223300      G2002221011312002112001121      1       0       1       1       1       -1809022        -1810770        AAA
#  926_33_329      T3000003020011330112020200      G2200202231211010022113031      0       3       3       1       1       -152955 -154703 AAA
#  926_34_440      T0121010212320132021200211      G0312100032003101301113001      2       1       3       1       1       -1618712        -1620852        AAA
#  926_38_533      T0011012010023302310021321      G0200103332011100011013103      0       4       4       1       1       -251624 -253847 AAA
#  926_42_329      T2331321031230210321102001      G2330012312202211022311020      1       3       4       1       1       -1093994        -1098892        AAA
#
#   0 beadId
#   1 F3 sequence
#   2 R3 sequence
#   3 num F3 mismatches
#   4 num R3 mismatches
#   5 total mismatches
#   6 F3 reference
#   7 R3 reference
#   8 F3 position
#   9 R3 position
#  10 category
#
#
#  tapperconvert:
#  M 12345_926_28_197 0 633481 f 0 0/0/3 54321_926_28_197 0 631038 f 0 0/0/2
#  M 12345_926_28_374 0 1808998 r 0 0/0/1 54321_926_28_374 0 1810746 r 0 0/0/0
#  M 12345_926_29_486 0 129939 f 0 0/0/2 54321_926_29_486 0 127944 f 0 0/0/4
#  M 12345_926_33_329 0 152931 r 0 0/0/0 54321_926_33_329 0 154679 r 0 0/0/3
#  M 12345_926_34_440 0 1618688 r 0 1/2/0 54321_926_34_440 0 1620828 r 0 0/0/1
#  M 12345_926_38_533 0 251600 r 0 0/0/0 54321_926_38_533 0 253823 r 0 0/0/4

my $tinput  = shift @ARGV;
my $terrors = shift @ARGV;
my $cinput  = "pgingivali.F3_R3.mates";
my $cerrors = 3;

if (!defined($tinput) || !defined($terrors)) {
    die "usage: $0 tapper-input-prefix num-errors\n";
}

print STDERR "Reading tangles.\n";
my %tangled;
open(TT, "./tapperconvert -dumpt $tinput |") or die;
while (<TT>) {
    my @v = split '\s+', $_;
    if ($v[1] =~ m/^\d+_(\d+_\d+_\d+)$/) {
        $v[1] = $1;
    }
    $tangled{$v[1]}++;
}
close(TT);

print STDERR "Reading tapper mate for counts.\n";
my %tcounts;
open(TT, "./tapperconvert -dumpm $tinput |") or die;
while (<TT>) {
    my @v = split '\s+', $_;
    if ($v[1] =~ m/^\d+_(\d+_\d+_\d+)$/) {
        $v[1] = $1;
    }
    $tcounts{$v[1]}++;
}
close(TT);

print STDERR "Processing.\n";

open(FC, "< $cinput") or die;
open(FT, "./tapperconvert -dumpm $tinput |") or die;

open(GC, "> compare.pl.corona.out") or die;
open(GT, "> compare.pl.tapper.out") or die;

my $same = 0;
my $qual = 0;
my $diffmultiple = 0;
my $diff = 0;

my $onlyc = 0;
my $onlycerror = 0;
my $onlyctooshort = 0;
my $onlyctoolong = 0;
my $onlyctangled = 0;

my $onlyt = 0;
my $onlyterror = 0;
my $onlyttooshort = 0;
my $onlyttoolong = 0;
my $onlytdupl = 0;

my $cid  = undef;
my $cstr = undef;
my @c;

my $tid  = undef;
my $tstr = undef;
my @t;

while (!eof(FC) && !eof(FT)) {

    if (!defined($cid)) {
        $_ = <FC>;
        while (m/^#/) {
            $_ = <FC>;
        }
        my @v = split '\s+', $_;

        if ($v[10] ne "AAA") {
            goto again;
        }

        my $ori = "f";
        if ($v[8] < 0) {
            $v[8] = -int($v[8]);
            $v[9] = -int($v[9]);
            $ori  = "r";
        }

        my $dist = $v[9] - $v[8];
        if ($dist < 0) {
            $dist = -$dist;
        }

        $cid  = $v[0];
        $cstr = "$v[0] $v[3] $v[4] $ori $v[8] $v[9] $dist";
        $c[0] = $v[0];
        $c[1] = $v[3];
        $c[2] = $v[4];
        $c[3] = $ori;
        $c[4] = $v[8];
        $c[5] = $v[9];
        $c[6] = $dist;

        {
            my @xxx = split '_', $cid;
            $xxx[0] = substr("00000$xxx[0]", -5);
            $xxx[1] = substr("00000$xxx[1]", -5);
            $xxx[2] = substr("00000$xxx[2]", -5);

            $cid = "$xxx[0]$xxx[1]$xxx[2]";
        }

    }

    if (!defined($tid)) {
        $_ = <FT>;
        my @v = split '\s+', $_;

        if ($v[1] =~ m/^\d+_(\d+_\d+_\d+)$/) {
            $v[1] = $1;
        }
        if ($v[6] =~ m!\d+/(\d+)/(\d+)$!) {
            $v[6] = $1 + $2;
        }
        if ($v[12] =~ m!\d+/(\d+)/(\d+)$!) {
            $v[12] = $1 + $2;
        }

        #  Correct for reverse?  Why?
        if ($v[4] eq "r") {
            $v[3] += 24;
            $v[9] += 24;
        }

        my $dist = $v[9] - $v[3];
        if ($dist < 0) {
            $dist = -$dist;
        }

        $tid  = $v[1];
        $tstr = "$v[1] $v[6] $v[12] $v[4] $v[3] $v[9] $dist";
        $t[0] = $v[1];
        $t[1] = $v[6];
        $t[2] = $v[12];
        $t[3] = $v[4];
        $t[4] = $v[3];
        $t[5] = $v[9];
        $t[6] = $dist;

        {
            my @xxx = split '\D+', $tid;
            $xxx[0] = substr("00000$xxx[0]", -5);
            $xxx[1] = substr("00000$xxx[1]", -5);
            $xxx[2] = substr("00000$xxx[2]", -5);

            $tid = "$xxx[0]$xxx[1]$xxx[2]";
        }
    }

    if ($cid eq $tid) {
        print GC "$cstr\n";
        print GT "$tstr\n";

        if      ($cstr eq $tstr) {
            $same++;
        } elsif (($c[3] == $t[3]) && ($c[4] == $t[4]) && ($c[5] == $t[5]) && ($c[6] == $t[6])) {
            $qual++;
        } else {
            #print STDERR "DIFF $cstr == $tstr\n";
            if ($tcounts{$t[0]} > 1) {
                $diffmultiple++;
            } else {
                $diff++;
            }
        }

        undef $cid;
        undef $cstr;
        undef @c;

        undef $tid;
        undef $tstr;
        undef @t;
    } elsif ($cid lt $tid) {
        print GC "$cstr\n";

        if (($c[1] > $terrors) || ($c[2] > $terrors)) {
            $onlycerror++;
        } elsif ($c[6] < 1400) {
            $onlyctooshort++;
        } elsif ($c[6] > 2600) {
            $onlyctoolong++;
        } elsif (exists($tangled{$c[0]})) {
            #print STDERR "TANGLED $cstr\n";
            $onlyctangled++;
        } else {
            #print STDERR "MISSED $cstr\n";
            $onlyc++;
        }

        undef $cid;
        undef $cstr;
        undef @c;
    } else {
        print GT "$tstr\n";

        if (($t[1] > $cerrors) || ($t[2] > $cerrors)) {
            $onlyterror++;
        } elsif ($t[6] < 1400) {
            $onlyttooshort++;
        } elsif ($t[6] > 2600) {
            $onlyttoolong++;
        } elsif ($tcounts{$t[0]} > 1) {
            $onlytdupl++;
        } else {
            $onlyt++;
        }

        undef $tid;
        undef $tstr;
        undef @t;
    }

again:
}

print STDERR "same $same  qual $qual  diff $diff diffmultiple $diffmultiple\n";
print STDERR "onlyc $onlyc err $onlycerror short $onlyctooshort long $onlyctoolong TANGLED $onlyctangled\n";
print STDERR "onlyt $onlyt err $onlyterror short $onlyttooshort long $onlyttoolong DUPLICATE $onlytdupl\n";

close(FC);
close(FT);

close(GC);
close(GT);
