#!/usr/bin/perl

use strict;

my $atacFile    = undef;
my $minLen      = 10000;
my $maxChr      = 24;
my $reference   = 0;

while (scalar(@ARGV)) {
    my $arg = shift @ARGV;

    if      ($arg eq "-m") {
        $minLen = shift @ARGV;
    } elsif ($arg eq "-c") {
        $maxChr = shift @ARGV;
    } elsif ($arg eq "-a") {
        $atacFile = shift @ARGV;
    } elsif ($arg eq "-A") {
        $reference = 0;
    } elsif ($arg eq "-B") {
        $reference = 1;
    } else {
        die "Invalid option '$arg'\n";
    }
}

if (!defined($atacFile)) {
    print STDERR "usage: $0 [-m minlen] [-c maxchr] [-A | -B] -a <some.atac.ckpLast>\n";
    print STDERR "  -m m      Include matches larger than 'm'.  Default: 10000\n";
    print STDERR "  -c c      Include chromosomes in the reference below 'c'.  Default: 24 (1-22+X+Y)\n";
    print STDERR "  -A | -B   The reference genome is sequence A (B).\n";
    print STDERR "  -a x      Process matches from atac file 'x'\n";
    exit(1);
}

#if (! -e "$atacFile.gaps.sorted") {
    open(F, "< $atacFile");
    open(G, "| sort -k1n -k2n > $atacFile.gaps.sorted");
    while (<F>) {
        if (m/^M\sr\s/) {
            my @vals = split '\s+', $_;
            (undef, $vals[4]) = split ':', $vals[4];
            (undef, $vals[8]) = split ':', $vals[8];

            if (($reference == 0) && ($vals[4] < $maxChr) && ($vals[6] > $minLen)) {
                print G "$vals[4] $vals[5] $vals[6] - $vals[8] $vals[9] $vals[10] - $vals[11]\n";
            }
            if (($reference == 1) && ($vals[8] < $maxChr) && ($vals[10] > $minLen)) {
                print G "$vals[8] $vals[9] $vals[10] - $vals[4] $vals[5] $vals[6] - $vals[11]\n";
            }
        }
    }
    close(F);
    close(G);
#}


my $lastscf = -1;
my @chr;
my @pos;
my @len;
my @scf;

open(F, "< $atacFile.gaps.sorted");
while (<F>) {
    my @vals = split '\s+', $_;
    push @chr, $vals[0];
    push @pos, $vals[1];
    push @len, $vals[2];
    push @scf, $vals[4];
}
close(F);

push @scf, -1;
push @scf, -1;

my $num = scalar(@chr) - 1;

#  We compute stats on the distance between elements i and j
#
my $i   = 0;
my $j   = 1;

print "GAPS\n";

while ($i < $num) {

    $j = $i+1;

    if ($chr[$i] != $chr[$j]) {
        $i++;
        next;
    }

  again:

    #  Move j ahead if it's an interleaved scaffold
    #
    #  If our current scaffold is interleaved by someone else,
    #  skip that someone else.  Yes, the end result of this
    #  is to have i and j point to the same scaffold.
    #
    if (($chr[$i] == $chr[$j]) &&
        ($scf[$i] == $scf[$j+1]) &&
        ($len[$j] < 5000)) {
        $j++;
    }

    #  Move j ahead if it's the same scaffold
    #
    if (($chr[$i] == $chr[$j]) &&
        ($scf[$i] == $scf[$j])) {
        $i = $j;
        $j++;
        goto again;
    }



    #  Report, if begin and end are on the same chromosome.
    if ($chr[$i] == $chr[$j]) {
        my $aend = $pos[$i] + $len[$i];
        my $bsta = $pos[$j];
        my $gapl = $bsta - $aend;

        if ($gapl > 0) {
            print "GAP: $gapl -- $i ($chr[$i] $pos[$i] $len[$i] $scf[$i]) -- $j ($chr[$j] $pos[$j] $len[$j] $scf[$j])\n";
        }
    }

    $i = $j;
}

