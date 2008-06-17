#!/usr/bin/perl

use strict;

#  Given an assembly (gkpStore), map the trimmed reads to a reference,
#  report half-matches and chimera.

my $gatekeeper = "wgs/FreeBSD-amd64/bin/gatekeeper";
my $snapper    = "wgs/kmer/FreeBSD-amd64/bin/snapper2";
my $convert    = "wgs/kmer/FreeBSD-amd64/bin/convertToExtent";

my $gkp = shift @ARGV;
my $ref = shift @ARGV;

die "usage: $0 asm/asm.gkpStore ref.fasta\n" if (!defined($gkp) || ! -d $gkp);
die "usage: $0 asm/asm.gkpStore ref.fasta\n" if (!defined($ref) || ! -e $ref);

if (! -e "$gkp/chimeric-read-test.fasta") {
    system("$gatekeeper -dumpfastaseq -clear OBT $gkp > $gkp/chimeric-read-test.fasta");
}
if (! -e "$gkp/chimeric-read-test.sim4db") {
    system("$snapper -queries $gkp/chimeric-read-test.fasta -genomic $ref -minhitlength 0 -minhitcoverage 0 -discardexonlength 40 -minmatchcoverage 25 -minmatchidentity 90 -verbose -noaligns -numthreads 1 > $gkp/chimeric-read-test.sim4db");
}
if (! -e "$gkp/chimeric-read-test.extents") {
    system("$convert < $gkp/chimeric-read-test.sim4db | sort -k1,1 -k4n -k5n > $gkp/chimeric-read-test.extents");
}

my %reads;

open(F, "< $gkp/chimeric-read-test.fasta") or die;
while (<F>) {
    if (m/^>(\S+)\s/) {
        $reads{$1} = 'u';
    }
}
close(F);

print STDERR "found ", scalar(keys %reads), " reads.\n";

open(F, "< $gkp/chimeric-read-test.extents") or die;
$_ = <F>;  #  labels
while (<F>) {
    my @v = split '\s+', $_;
    my $len;

    if ($v[3] < $v[4]) {
        $len = $v[4] - $v[3];
    } else {
        $len = $v[3] - $v[4];
    }

    if ($len + 10 >= $v[1]) {
        $reads{$v[0]} = 'f';
    } else {
        if ($reads{$v[0]} eq 'u') {
            $reads{$v[0]} = 'p';
        }
    }
}
close(F);

my %partial;
my %unmapped;

foreach my $id (keys %reads) {
    if ($reads{$id} eq 'p') {
        $partial{$id} = 1;
    }
    if ($reads{$id} eq 'u') {
        $unmapped{$id} = 1;
    }
}

print STDERR "Found ", scalar(keys %partial), " partially mapped reads.\n";
print STDERR "Found ", scalar(keys %unmapped), " completely unmapped reads.\n";

my %hits;

open(F, "< $gkp/chimeric-read-test.extents") or die;
$_ = <F>;  #  labels
while (<F>) {
    my @v = split '\s+', $_;

    if (exists($partial{$v[0]})) {
        $hits{$v[0]}++;
    }
}
close(F);

open(P, "> $gkp/chimeric-read-test.partial") or die;
open(C, "> $gkp/chimeric-read-test.chimeric") or die;

my @multi;
my $lastid;

open(F, "< $gkp/chimeric-read-test.extents") or die;
$_ = <F>;  #  labels
while (<F>) {
    my @v = split '\s+', $_;

    next if (!exists($partial{$v[0]}));

    if (($lastid ne $v[0]) && (scalar(@multi) > 0)) {
        @multi = sort { $a <=> $b } @multi;

        my ($beg, $end) = split '\s+', $multi[0];
        my $reports = 0;

        foreach my $m (@multi) {
            my ($b, $e) = split '\s+', $m;

            print C "$lastid\t$m\n";

            if ($end < $b) {
                print C "$lastid\t$beg\t$end\n";
                $beg = $b;
                $end = $e;
                $reports++;
            }
            if ($b < $beg) {
                #  if properly sorted, shouldn't occur
                print STDERR "SORT ERROR.\n";
                $beg = $b;
            }
            if ($e > $end) {
                $end = $e;
            }
        }

        print C "$lastid\t$beg\t$end\n";

        if ($reports > 0) {
            print STDERR "$lastid looks chimeric.\n";
        }

        undef @multi;
    }

    if ($hits{$v[0]} == 1) {
        print P $_;
    } else {
        if ($v[3] < $v[4]) {
            push @multi, "$v[3]\t$v[4]\t$v[6]\t$v[7]\t$v[8]\t$v[9]";
        } else {
            push @multi, "$v[4]\t$v[3]\t$v[6]\t$v[7]\t$v[8]\t$v[9]";
        }
    }

    $lastid = $v[0];
}
close(F);

close(C);
close(P);
