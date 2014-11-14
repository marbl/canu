#!/usr/bin/perl

use strict;

#  tail -f cgw.out | egrep Examine\|isQual.\*gap.\*weight

my @p;  #  previous numbers
my @r;  #  current numbers

$p[0]  = $r[0]  = 0;
$p[1]  = $r[1]  = 0;
$p[2]  = $r[2]  = 0;
$p[3]  = $r[3]  = 0;
$p[4]  = $r[4]  = 0;
$p[5]  = $r[5]  = 0;
$p[6]  = $r[6]  = 0;
$p[7]  = $r[7]  = 0;
$p[8]  = $r[8]  = 0;
$p[9]  = $r[9]  = 0;
$p[10] = $r[10] = 0;
$p[11] = $r[11] = 0;
$p[12] = $r[12] = 0;
$p[13] = $r[13] = 0;
$p[14] = $r[14] = 0;
$p[15] = $r[15] = 0;

my $nTests       = 0;
my $nInterleaves = 0;
my $nMerged      = 0;

my $nContigs     = 0;

my $iter         = 0;

open(C, "> cgw.out.cumulative.stats");
open(D, "> cgw.out.difference.stats");

while (!eof(STDIN)) {

    #                                                                     scf1     len1           scf2     len2             gap                weight  orient
    while ($_ !~ m/^isQualityScaffoldMergingEdge\(\)--\sMerge\sscaffolds\s(\d+)\s\((.*bp)\)\sand\s(\d+)\s\((.*bp)\)\:\sgap\s(.*bp.*bp)\sweight\s(\d+)\s(.*)\sedge/) {

        if (m/Insert\sscaffold/) {
            $nMerged++;
        }

        if (m/^MergeScaffoldsAggressive.*iter\s(\d+)\s/) {
            $iter = $1;

            for (my $i=0; $i<15; $i++) {
                $p[$i] = $r[$i] - $p[$i];
            }

            print C "$iter\t$nTests\t$nInterleaves\t$nMerged\t$nContigs\t$r[0]\t$r[1]\t$r[2]\t$r[3]\t$r[4]\t$r[5]\t$r[6]\t$r[7]\t$r[8]\t$r[9]\t$r[10]\t$r[11]\t$r[12]\t$r[13]\t$r[14]\n";
            print D "$iter\t$nTests\t$nInterleaves\t$nMerged\t$nContigs\t$p[0]\t$p[1]\t$p[2]\t$p[3]\t$p[4]\t$p[5]\t$p[6]\t$p[7]\t$p[8]\t$p[9]\t$p[10]\t$p[11]\t$p[12]\t$p[13]\t$p[14]\n";

            @p = @r;
        }

        if (m/^CreateAContigInScaffold.*new\scontig/) {
            $nContigs++;
        }

        $_ = <STDIN>;  chomp $_;
    }            

    $nTests++;

    my ($scf1, $len1, $scf2, $len2, $gapSize, $weight, $orient);

    #print "$_\n";
    if ($_ =~ m/^isQualityScaffoldMergingEdge\(\)--\sMerge\sscaffolds\s(\d+)\s\((.*bp)\)\sand\s(\d+)\s\((.*bp)\)\:\sgap\s(.*bp.*bp)\sweight\s(\d+)\s(.*)\sedge/) {
        $scf1    = $1;
        $len1    = $2;
        $scf2    = $3;
        $len2    = $4;
        $gapSize = $5;
        $weight  = $6;
        $orient  = $7;
    } else {
        die "Failed on '$_'\n";
    }


    my ($happy1, $sad1);

    $_ = <STDIN>;  chomp $_;  #  First scaffold instrument

    while (m/^instrument/) {
        $_ = <STDIN>;  chomp $_;  #  Ignore lines about setting up the instrumenter
    }

    #print "$_\n";
    if ($_ =~ m/isQualityScaffoldMergingEdge\(\)--\s+scaffold\s(\d+)\sinstrumenter\shappy\s(.*)\sgap\s(.*)\smisorient\sclose\s(.*)\scorrect\s(.*)\sfar\s(.*)\soriented\sclose\s(.*)\sfar\s(.*)\smissing\s(.*)\sexternal\s(.*)/) {
         $happy1 = $1;
         $sad1   = $3 + $4 + $5 + $6 + $7 + $8;
    } else {
        die "Failed on '$_'\n";
    }


    my ($happy2, $sad2);

    #print "$_\n";
    $_ = <STDIN>;  chomp $_;  #  Second scaffold instrument
    if ($_ =~ m/isQualityScaffoldMergingEdge\(\)--\s+scaffold\s(\d+)\sinstrumenter\shappy\s(.*)\sgap\s(.*)\smisorient\sclose\s(.*)\scorrect\s(.*)\sfar\s(.*)\soriented\sclose\s(.*)\sfar\s(.*)\smissing\s(.*)\sexternal\s(.*)/) {
        $happy2 = $1;
        $sad2   = $3 + $4 + $5 + $6 + $7 + $8;
    } else {
        die "Failed on '$_'\n";
    }


    my ($happy, $gap, $misclose, $misfar, $mis, $close, $far, $missing, $external);
    my ($happyN, $sadN);

    #print "$_\n";
    $_ = <STDIN>;  chomp $_;  #  Merged scaffold instrument
    if ($_ =~ m/isQualityScaffoldMergingEdge\(\)--\s+scaffold\s\(new\)\sinstrumenter\shappy\s(.*)\sgap\s(.*)\smisorient\sclose\s(.*)\scorrect\s(.*)\sfar\s(.*)\soriented\sclose\s(.*)\sfar\s(.*)\smissing\s(.*)\sexternal\s(.*)/) {
        $happy    = $1;
        $gap      = $2;
        $misclose = $3;
        $misfar   = $4;
        $mis      = $5;
        $close    = $6;
        $far      = $7;
        $missing  = $8;
        $external = $9;

        $happyN = $1;
        $sadN   = $3 + $4 + $5 + $6 + $7 + $8;
    } else {
        die "Failed on '$_'\n";
    }


    my $happyB = $happy1 + $happy2;
    my $sadB   = $sad1 + $sad2;

    my $interleave;

    $_ = <STDIN>;  chomp $_;  #  Before/After mates
    $_ = <STDIN>;  chomp $_;  #  Happy enough to merge?
    $_ = <STDIN>;  chomp $_;  #  pass/fail stats
    $_ = <STDIN>;  chomp $_;  #  Interleave result (optional)

    if ($_ =~ m/ExamineSEdgeForUsability_Interleaved\(\)--\s(.*)/) {
        $nInterleaves++;

        $interleave = $1;

        #print  "$scf1 ($len1) $scf2 ($len2) $gapSize $weight -- + $happy - $misclose $mis $misfar - $close $far - $missing -- $interleave\n";
        #print  "$scf1 ($len1) $scf2 ($len2) gap $gapSize weight $weight -- $happyB / $sadB -- $happyN / $sadN -- $interleave\n";
    } else {
        #print  "$scf1 ($len1) $scf2 ($len2) $gapSize $weight -- + $happy - $misclose $mis $misfar - $close $far - $missing\n";
        #print  "$scf1 ($len1) $scf2 ($len2) gap $gapSize weight $weight -- $happyB / $sadB -- $happyN / $sadN\n";
    }

    #  The classifications, straight from CIScaffoldT_Merge_Interleaved.c

    if      ($interleave eq undef) {
        $r[0]++;
    } elsif($interleave =~ m/Didn't expect end contigs to overlap, didn't find it, and will join with existing edge/) {
        $r[1]++;
    } elsif ($interleave =~ m/Expected end contigs to overlap, found overlap, will merge/) {
        $r[2]++;
    } elsif ($interleave =~ m/Expected end contigs to overlap, didn't find it, but will abut/) {
        $r[3]++;
    } elsif ($interleave =~ m/Expected end contigs to overlap, didn't find it, and abutting not allowed, will not merge/) {
        $r[4]++;
    } elsif ($interleave =~ m/Expected end contigs to overlap, didn't find it, will not merge/) {
        $r[5]++;
    } elsif ($interleave =~ m/Interleaving succeeded with contig overlaps, but failed to make adjustments; will not merge/) {
        $r[6]++;
    } elsif ($interleave =~ m/Interleaving succeeded with contig overlaps; will merge/) {
        $r[7]++;
    } elsif ($interleave =~ m/Interleaving succeeded without contig overlaps, but failed to make adjustments; will not merge/) {
        $r[8]++;
    } elsif ($interleave =~ m/Interleaving succeeded without contig overlaps; will merge/) {
        $r[9]++;
    } elsif ($interleave =~ m/Interleaving failed, abut not allowed, will not merge/) {
        $r[10]++;
    } elsif ($interleave =~ m/Interleaving failed, will abut/) {
        $r[11]++;
    } elsif ($interleave =~ m/Interleaving failed, didn't expect overlap, will join with existing edge/) {
        $r[12]++;
    } elsif ($interleave =~ m/Interleaving failed, will not merge/) {
        $r[13]++;
    } else {
        $r[14]++;
    }
}

close(C);
close(D);
