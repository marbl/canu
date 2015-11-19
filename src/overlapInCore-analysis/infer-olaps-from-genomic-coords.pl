#!/usr/bin/env perl

###############################################################################
 #
 #  This file is part of canu, a software program that assembles whole-genome
 #  sequencing reads into contigs.
 #
 #  This software is based on:
 #    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 #    the 'kmer package' (http://kmer.sourceforge.net)
 #  both originally distributed by Applera Corporation under the GNU General
 #  Public License, version 2.
 #
 #  Canu branched from Celera Assembler at its revision 4587.
 #  Canu branched from the kmer project at its revision 1994.
 #
 #  Modifications by:
 #
 #    Brian P. Walenz beginning on 2015-OCT-12
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

use strict;

my $outputPrefix  = $ARGV[0];
my $readsFile     = $ARGV[1];
my $coordsFile    = $ARGV[2];
my $uidMapFile    = $ARGV[3];

my $orientForward = 0;

#
#  Load name to IID map
#

my %NAMEtoIID;
my %UIDtoNAME;

if (defined($uidMapFile)) {
    print STDERR "Read names from '$uidMapFile'\n";

    open(F, "< $uidMapFile");
    while (<F>) {
        my @v = split '\s+', $_;

        #  UID IID NAME  UID IID NAME

        if (scalar(@v) == 3) {
            $NAMEtoIID{$v[2]} = $v[1];
            $UIDtoNAME{$v[0]} = $v[2];
        }

        if (scalar(@v) == 6) {
            $NAMEtoIID{$v[2]} = $v[1];
            $NAMEtoIID{$v[5]} = $v[4];
            $UIDtoNAME{$v[0]} = $v[2];
            $UIDtoNAME{$v[3]} = $v[5];
        }
    }
    close(F);
}

#
#  Load read lengths.  Needed to find reads that do not map full length, and to compute hangs.
#  Assumes reads are dumped from gkpStore and have names 'uid,iid'.

my %readLength;
my %totalAligns;
my %goodAligns;

print STDERR "Read lengths from '$readsFile'\n";

open(F, "< $readsFile") or die "Couldn't read '$readsFile'\n";
while (!eof(F)) {
    my $a = <F>;
    my $b = <F>;
    my $c = <F>;
    my $d = <F>;

    if ($a =~ m/^.(\S+)$/) {
        if (! exists($NAMEtoIID{$1})) {
            #print STDERR "WARNING: read '$1' not in map.\n"
        } else {
            $readLength{$1} = length($b) - 1;
            $totalAligns{$1} = 0;
            $goodAligns{$1} = 0;
        }
    } else {
        chomp;
        die "Failed to parse read name from '$a'\n";
    }
}
close(F);





my @aligns;
my $totalAligns = 0;
my $goodAligns  = 0;

my $readOrder = 0;
my %positions;


open(F, "< $coordsFile");

$_ = <F>;
$_ = <F>;
$_ = <F>;
$_ = <F>;

while (<F>) {
    my @v = split '\s+', $_;

    my $b = $v[0];  #  reference begin
    my $e = $v[1];  #  reference end

    my $r = $v[9];  #  reference name
    my $n = $v[10]; #  read name

    my $f = 1;
    my $L = $v[2];  #  read begin, unused
    my $R = $v[3];  #  read end, unused

    if ($R < $L) {
        $f = 0;
        $L = $v[3];
        $R = $v[2];
    }

    die "L=$L R=$R\n" if ($R < $L);
    die "b=$b e=$e\n" if ($e < $b);

    next                                           if (!exists($readLength{$n}));
    die "Didn't find read length for read '$n'\n"  if (!exists($readLength{$n}));

    $totalAligns++;
    $totalAligns{$n}++;

    next  if (10 < $L);
    next  if (10 < $readLength{$n} - $R);

    $goodAligns++;
    $goodAligns{$n}++;

    my $scale = ($R - $L) / ($e - $b);  #  scale from reference coords to read coords

    #  The aligns line starts with reference name.  This lets us group aligns per reference sequence.

    push @aligns, "$r\0$b\0$e\0$f\0$scale\0$n";
}


#  Generate olaps.  We sort by reference name, then step through all aligns to that reference and
#  sort by begin coordinate.

my $olaps = 0;

open(OVA, "> $outputPrefix.ova");

@aligns = sort @aligns;

while (scalar(@aligns) > 0) {
    my @coords;

    my @v = split '\0', $aligns[0];
    my $r = $v[0];

    #  Find all coords for this reference sequence

    while ($r eq $v[0]) {
        push @coords, "$v[1]\0$v[2]\0$v[3]\0$v[4]\0$v[5]";

        shift @aligns;
        @v = split '\0', $aligns[0];
    }

    #  Sort those coords by start position

    @coords = sort { $a <=> $b } @coords;

    #  Process each one

    while (scalar(@coords) > 0) {
        my ($rb, $re, $rf, $rs, $rn) = split '\0', $coords[0];

        die "Reversed $coords[0]\n" if ($re < $rb);

        next if (!exists($NAMEtoIID{$rn}));  #  Just skip reads we don't care about.
        die  if (!exists($NAMEtoIID{$rn}));  #  For development, fail.

        #  Save this reads position

        if (exists($positions{$rn})) {
            print STDERR "WARNING: read $rn already has a position!\n";
        }
        $positions{$rn} = $readOrder++;

        #  Relabel the read name to an IID.

        $rn = $NAMEtoIID{$rn};

        shift @coords;

        foreach my $cc (@coords) {
            my ($cb, $ce, $cf, $cs, $cn) = split '\0', $cc;
            my $olap   = $re - $cb;

            last if ($olap < 40);

            my $ahang;
            my $bhang;

            if ($rf == 1) {
                #  Read is forward
                $ahang = $cb - $rb;
                $bhang = $ce - $re;
            } else {
                #  Read is reverse; hangs are swapped and negative from forward case
                $ahang = $re - $ce;
                $bhang = $rb - $cb;
            }

            next if (!exists($NAMEtoIID{$cn}));
            die  if (!exists($NAMEtoIID{$cn}));

            $cn = $NAMEtoIID{$cn};

            #  Scale the hangs.  Hangs are computed using genomic positions, but deletions from the read
            #  cause this hang to be too large.

            if ($ahang < 0) {
                $ahang = $ahang * $cs;  #  hang is composed of the C read
            } else {
                $ahang = $ahang * $rs;  #  hang is composed of the R read
            }

            if ($bhang < 0) {
                $bhang = $bhang * $rs;  #  hang is composed of the R read
            } else {
                $bhang = $bhang * $cs;  #  hang is composed of the C read
            }

            printf(OVA "%9s %9s %s %5d %5d %.2f %.2f\n",
                   $rn,
                   $cn,
                   ($rf == $cf) ? "N" : "I",
                   $ahang,
                   $bhang,
                   0.0,
                   0.0);

            $olaps++;
        }
    }
}

close(OVA);

#  Output reads that aligned, didn't align.

open(R, "< $readsFile");
open(M, "> $outputPrefix.mapped.fastq");
open(P, "> $outputPrefix.partial.fastq");
open(F, "> $outputPrefix.failed.fastq");

my $mappedReads  = 0;
my $partialReads = 0;
my $failedReads  = 0;

my @mappedReads;
my @orientReads;

while (!eof(R)) {
    my $a = <R>;
    my $b = <R>;
    my $c = <R>;
    my $d = <R>;
    my $n;

    if ($a =~ m/\@(\S+)\s*/) {
        $n = $1;
    } else {
        chomp;
        die "Failed to parse read name from '$a'\n";
    }

    #die "Failed to find UIDtoNAME for UID '$n'\n"  if (!exists($UIDtoNAME{$n}));
    #
    #$a = "\@$UIDtoNAME{$n}\n";

    if      ($goodAligns{$n} > 0) {
        $mappedReads++;
        print M "$a$b$c$d";

        die "No position for '$n'\n" if (!exists($positions{$n}));

        push @mappedReads, "$positions{$n}\0$a$b$c$d";
        push @orientReads, "$positions{$n}\0$a$b$c$d";

    } elsif ($totalAligns{$n} > 0) {
        $partialReads++;
        print P "$a$b$c$d";

    } else {
        $failedReads++;
        print F "$a$b$c$d";
    }
}

close(F);
close(P);
close(M);
close(R);

@mappedReads = sort { $a <=> $b } @mappedReads;
@orientReads = sort { $a <=> $b } @orientReads;

open(A, "> $outputPrefix.mapped.ordered.fastq");
foreach my $a (@mappedReads) {
    my ($idx,$dat) = split '\0', $a;
    print A "$dat";
}
close(A);

open(A, "> $outputPrefix.mapped.ordered.oriented.fastq");
foreach my $a (@orientReads) {
    my ($idx,$dat) = split '\0', $a;
    print A "$dat";
}
close(A);


print STDERR "RESULT $outputPrefix ", scalar(keys %readLength), " input $mappedReads mapped $partialReads partial $failedReads unaligned reads, $totalAligns aligns, $goodAligns good aligns, $olaps overlaps.\n";
