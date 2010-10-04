#!/usr/bin/perl

use strict;

my $TOOSHORT = 500;

if (! -e "tig.fasta") {
    die "Need tig.fasta.\n";
}

if (! -e "tig.coords") {
    die "Need tig.coords.\n";
}

#  Load the lengh of each sequence.

my %length;

if (-e "tig.fasta") {
    open(F, "< tig.fasta") or die;
    while (<F>) {
        if (m/^>(utg\d+)\s+len=(\d+)$/) {
            $length{$1} = $2;
        }
    }
    close(F);
} else {
    print STDERR "No tig.fasta file found, cannot .....\n";
}

#  Pass one, load the coords.  Analyze anything with one match, counting
#  the amount of coverage.

my %nucmer;
my %nMatch;

my $matches = 0;
my $matchesSmall = 0;

if (-e "tig.coords") {
    open(F, "< tig.coords") or die;
    $_ = <F>;
    $_ = <F>;
    $_ = <F>;
    $_ = <F>;
    $_ = <F>;
    while (<F>) {
        s/^\s+//;
        s/\s$+//;
        my @v = split '\s+', $_;
        my $str = "$v[3]\t$v[4]\t$v[0]\t$v[1]";

        my $len = ($v[3] < $v[4]) ? $v[4] - $v[3] : $v[3] - $v[4];

        #  Filter out small matches due to repeats
        if ($len < $TOOSHORT) {
            $matchesSmall++;
            next;
        } else {
            $matches++;
        }

        if (exists($nucmer{$v[12]})) {
            $nucmer{$v[12]} .= "\n" . $str;
            $nMatch{$v[12]}++;
        } else {
            $nucmer{$v[12]} = $str;
            $nMatch{$v[12]} = 1;
        }
    }
    close(F);
} else {
    die "No tig.coords??\n";
}

print STDERR "matches $matches matchesSmall $matchesSmall\n";

my $spur = 0;
my $spurShort = 0;

foreach my $k (keys %nucmer) {
    next if ($nMatch{$k} != 1);

    my @v = split '\s+', $nucmer{$k};

    die if (scalar(@v) != 4);

    my $len = ($v[0] < $v[1]) ? $v[1] - $v[0] : $v[0] - $v[1];

    die if ($len < 0);

    if ($len < 0.95 * $length{$k}) {
        if ($len < $TOOSHORT) {
            $spurShort++;
        } else {
            print STDERR "Possible spur in $k\t- length $length{$k}\t- mapped $len\t- coords $v[0]\t$v[1]\n";
            $spur++;
        }
    }

    delete $nucmer{$k};
    delete $nMatch{$k};
}

print STDERR "spur $spur spurShort $spurShort\n";
print STDERR "matches ", scalar(keys %nucmer), "\n";


foreach my $k (keys %nucmer) {
    die if ($nMatch{$k} < 2);
    my @m = split '\n', $nucmer{$k};

    print STDERR "$k\n";

    foreach my $m (@m) {
        my @v = split '\s+', $m;

        die if (scalar(@v) != 4);

        print STDERR "$v[0]\t$v[1]\t$v[2]\t$v[3]\n";
    }
}
