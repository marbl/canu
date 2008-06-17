#!/usr/bin/perl

use strict;

#  Compare trimmed reads against the untrimmed mapping to the reference.

my $gatekeeper = "wgs/FreeBSD-amd64/bin/gatekeeper";
my $snapper    = "wgs/kmer/FreeBSD-amd64/bin/snapper2";
my $convert    = "wgs/kmer/FreeBSD-amd64/bin/convertToExtent";

my $gkp = shift @ARGV;
my $ref = shift @ARGV;

die "usage: $0 asm/asm.gkpStore ref.fasta\n" if (!defined($gkp) || ! -d $gkp);
die "usage: $0 asm/asm.gkpStore ref.fasta\n" if (!defined($ref) || ! -e $ref);

if (! -e "$gkp/reads-untrimmed.fasta") {
    system("$gatekeeper -dumpfastaseq -clear OBT $gkp > $gkp/reads-untrimmed.fasta");
}
if (! -e "$gkp/reads-untrimmed.sim4db") {
    system("$snapper -queries $gkp/reads-untrimmed.fasta -genomic $ref -minhitlength 0 -minhitcoverage 0 -discardexonlength 40 -minmatchcoverage 25 -minmatchidentity 90 -verbose -noaligns -numthreads 1 > $gkp/reads-untrimmed.sim4db");
}
if (! -e "$gkp/reads-untrimmed.extents") {
    system("$convert < $gkp/reads-untrimmed.sim4db | sort -k1,1 -k4n -k5n > $gkp/reads-untrimmed.extents");
}

my %mapBeg;
my %mapEnd;

open(F, "< $gkp/reads-untrimmed.extents");
while (<F>) {
    my @v = split '\s+', $_;

    if ($v[3] < $v[4]) {
        $mapBeg{$v[0]} = $v[3];
        $mapEnd{$v[0]} = $v[4];
    } else {
        $mapBeg{$v[0]} = $v[4];
        $mapEnd{$v[0]} = $v[3];
    }
}
close(F);

my $falseDelete    = 0;
my $noMapping      = 0;
my $wrongChimerEnd = 0;
my $trimTooShort   = 0;
my $trimTooLong    = 0;
my $trimPerfect    = 0;

my $tooShort = 0;
my $tooLong  = 0;

open(F, "$gatekeeper -dumpfragments -tabular -clear OBT $gkp |");
$_ = <F>;  # labels
while (<F>) {
    my @v = split '\s+', $_;

    $v[0] = "$v[0],$v[1]";

    if ($v[6]) {
        #  deleted fragment
        if (exists($mapBeg{$v[0]})) {
            $falseDelete++;
            #print "FALSE DELETE $v[0]\n";
        }
    } else {
        my $mapBeg = $mapBeg{$v[0]};
        my $mapEnd = $mapEnd{$v[0]};
        my $mapLen = $mapEnd - $mapBeg;

        my $trimBeg = $v[11];
        my $trimEnd = $v[12];
        my $trimLen = $trimEnd - $trimBeg;

        if (!defined($mapBeg)) {
            $noMapping++;
        } elsif ($trimLen < $mapLen) {
            $trimTooShort++;
            $tooShort += $mapLen - $trimLen;
            #print "$v[0]\tTRIM\t$trimBeg\t$trimEnd\tMAP\t$mapBeg\t$mapEnd\tSHORT\n";
        } elsif (($trimEnd - 20 < $mapBeg) || ($trimBeg > $mapEnd - 20)) {
            $wrongChimerEnd++;
        } elsif ($trimLen > $mapLen) {
            $trimTooLong++;
            $tooLong += $trimLen - $mapLen;
            if ($trimLen - $mapLen > 30) {
                #print "$v[0]\tTRIM\t$trimBeg\t$trimEnd\tMAP\t$mapBeg\t$mapEnd\tLONG\n";
            }
        } else {
            $trimPerfect++;
        }
    }
}
close(F);

print "false delete        $falseDelete\n";
print "unmapped            $noMapping\n";
print "wrong chimer end    $wrongChimerEnd\n";
print "trim perfect        $trimPerfect\n";
print "trim too short      $trimTooShort\n";
print "trim too long       $trimTooLong\n";

$trimTooShort = 1 if ($trimTooShort == 0);
$trimTooLong  = 1 if ($trimTooLong == 0);

print "too short ", $tooShort / $trimTooShort, "\n";
print "too long  ", $tooLong / $trimTooLong, "\n";
