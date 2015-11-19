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

my $alignedReads = shift @ARGV;

die "Can't find blasr-coords aligned reads in '$alignedReads'\n"  if (! -e "$alignedReads");

my $b1last = 0;
my $e1last = 0;
my $n1last = "";
my $f1last = "fwd";

my $b2last = 0;
my $e2last = 0;
my $n2last = "";
my $f2last = "fwd";

my $nBubble = 0;
my $nSpan   = 0;
my $nGap    = 0;
my $nGapEv  = 0;
my $nAbut   = 0;
my $nOvl    = 0;
my $nOvlEv  = 0;

#  Load the name to iid map.

my %NAMtoIID;
my %NAMtoUID;
my %UIDtoIID;

my %IIDtoDEG;
my %IIDtoUTG;

open(F, "< test.gkpStore.fastqUIDmap") or die "Failed to open 'test.gkpStore.fastqUIDmap'\n";
while (<F>) {
    my @v = split '\s+', $_;

    if (scalar(@v) == 3) {
        $NAMtoIID{$v[2]} = $v[1];
        $NAMtoUID{$v[2]} = $v[0];
        $UIDtoIID{$v[0]} = $v[1];
    } else {
        $NAMtoIID{$v[2]} = $v[1];
        $NAMtoUID{$v[2]} = $v[0];
        $UIDtoIID{$v[0]} = $v[1];

        $NAMtoIID{$v[5]} = $v[4];
        $NAMtoUID{$v[5]} = $v[3];
        $UIDtoIID{$v[3]} = $v[4];
    }
}
close(F);

print STDERR "Loaded ", scalar(keys %NAMtoIID), " read names.\n";

open(F, "< 9-terminator/test.posmap.frgdeg") or die "Failed to open '9-terminator/test.posmap.frgdeg'\n";
while (<F>) {
    my @v = split '\s+', $_;
    my $i = $UIDtoIID{$v[0]};

    die "$_" if (!defined($i));

    $IIDtoDEG{$i} = $v[1];
}
close(F);

print STDERR "Loaded ", scalar(keys %IIDtoDEG), " degenerate unitig reads.\n";

open(F, "< 9-terminator/test.posmap.frgutg") or die "Failed to open '9-terminator/test.posmap.frgdeg'\n";
while (<F>) {
    my @v = split '\s+', $_;
    my $i = $UIDtoIID{$v[0]};

    die "$_" if (!defined($i));

    $IIDtoUTG{$i} = $v[1];
}
close(F);

print STDERR "Loaded ", scalar(keys %IIDtoUTG), " unitig reads.\n";


#  UTG.coords is standard nucmer output.

open(F, "< 9-terminator/UTG.coords") or die "Failed to open '9-terminator/UTG.coords'\n";

$_ = <F>;
$_ = <F>;
$_ = <F>;
$_ = <F>;
$_ = <F>;

while (<F>) {
    s/^\s+//;
    s/\s+$//;

    my @v = split '\s+', $_;

    my $b1 = $v[0];  #  Assembly coords
    my $e1 = $v[1];

    my $b2 = $v[3];  #  Read coords
    my $e2 = $v[4];

    my $f1 = ($b1 < $e1) ? "fwd" : "rev";
    my $f2 = ($b2 < $e2) ? "fwd" : "rev";

    ($b1,$e1) = ($e1,$b1)  if ($b1 > $e1);  #  Doesn't seem to occur.
    ($b2,$e2) = ($e2,$b2)  if ($b2 > $e2);  #  Definitely does occur.

    #my $l1 = $v[6];
    #my $l2 = $v[7];

    my $id = $v[9];

    my $n1 = $v[11];
    my $n2 = $v[12];

    if (($b1last <= $b1) && ($e1 <= $e1last)) {
        $nBubble++;
        #print STDERR "BUBBLE $n2 ($b1,$e1) in $n2last ($b1last,$e1last)\n";
        next;
    }

    $nSpan++;

    goto bail  if ($n1last ne $n1);

    #print STDERR "last $b1last,$e1last curr $b1,$e1\n";

    die "coords not sorted by increasing assembly start!\n"  if ($b1last > $b1);

    $nGap++    if ($e1last <  $b1);
    $nAbut++   if ($e1last == $b1);
    $nOvl++    if ($e1last  > $b1);


    my $minOvlSpan = 150;


    if      ($e1last < $b1) {
        my $ge   = $b1;
        my $gb   = $e1last;
        my $gap  = $ge - $gb;
        my $span = 0;

        open(C, "< $alignedReads") or die "Failed to open '$alignedReads'\n";
        while (<C>) {
            chomp;

            my @v = split '\s+', $_;

            next if ($v[9] ne $n1);

            #  STOP changing the read name, blasr!

            if ($v[10] =~ m/^(.*\d+_\d+)\/\d+_\d+/) {
                $v[10] = $1;
            }

            my $iid = $NAMtoIID{$v[10]};
            my $utg = $IIDtoUTG{$iid};
            my $deg = $IIDtoDEG{$iid};
            my $ann;

            #next if (defined($utg) || defined($deg));

            if      (defined($utg)) {
                $ann = "utg $utg";
            } elsif (defined($deg)) {
                $ann = "deg $deg";
            } else {
                $ann = "singleton       ";
            }

            if (($v[0] + $minOvlSpan < $ge) && ($gb < $v[1] - $minOvlSpan)) {
                #print "$ann iid $iid -- $_\n";
                $span++;
            }
        }
        close(C);

        if ($span > 0) {
            $nGapEv++;
            print "GAP $n1 $n2last ($b1last,$e1last,$f2last) -- $gap -- $n2 ($b1,$e1,$f2) -- SPAN $span\n";
            print "\n";
            print "\n";
        }
    }


    if ($e1last == $b1) {
    }


    if      ($e1last > $b1) {
        my $ge   = $e1last;
        my $gb   = $b1;
        my $ovl  = $ge - $gb;
        my $span = 0;

        open(C, "< $alignedReads") or die "Failed to open '$alignedReads'\n";
        while (<C>) {
            chomp;

            my @v = split '\s+', $_;

            next if ($v[9] ne $n1);

            if ($v[10] =~ m/^(.*\d+_\d+)\/\d+_\d+/) {
                $v[10] = $1;
            }

            my $iid = $NAMtoIID{$v[10]};
            my $utg = $IIDtoUTG{$iid};
            my $deg = $IIDtoDEG{$iid};
            my $ann;

            #next if (defined($utg) || defined($deg));

            if      (defined($utg)) {
                $ann = "utg $utg";
            } elsif (defined($deg)) {
                $ann = "deg $deg";
            } else {
                $ann = "singleton       ";
            }

            if (($v[0] + $minOvlSpan < $ge) && ($gb < $v[1] - $minOvlSpan)) {
                print "$ann iid $iid -- $_\n";
                $span++;
            }
        }
        close(C);

        if ($span > 0) {
            $nOvlEv++;
            print "OVL $n1 $n2last ($b1last,$e1last,$f2last) -- $ovl -- $n2 ($b1,$e1,$f2) -- SPAN $span\n";
            print "\n";
            print "\n";
        }
    }


  bail:
    ($b1last,$e1last,$n1last,$f1last) = ($b1,$e1,$n1,$f1);
    ($b2last,$e2last,$n2last,$f2last) = ($b2,$e2,$n2,$f2);
}
close(F);

print STDERR "nBubble  $nBubble\n";
print STDERR "nSpan    $nSpan\n";

print STDERR "nGap     $nGap\n";
print STDERR "nGap     $nGapEv with evidence\n";
print STDERR "nAbut    $nAbut\n";
print STDERR "nOvl     $nOvl\n";
print STDERR "nOvl     $nOvlEv with evidence\n";
