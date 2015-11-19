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

my $prefix = shift @ARGV;
my @reads;

my @tigStore;

#  If no read supplied, scan the tigStore for all unitigs that have a single read spanning
#  This read should be first.

print STDERR "Searching in assembly prefix '$prefix'.\n";

print STDERR "Loading tigStore.\n";
open(F, "tigStore -g $prefix.gkpStore -t $prefix.tigStore 1 -U -d layout |");
while (!eof(F)) {
    my $utg  = <F>;
    my $len  = <F>;
    my $cns  = <F>;
    my $qlt  = <F>;
    my $d1   = <F>;  #  cov stat
    my $d2   = <F>;  #  microhet
    my $d3   = <F>;  #  status
    my $d4   = <F>;  #  unique_rept
    my $d5   = <F>;  #  status
    my $nfrg = <F>;  #  num frags
    my $d7   = <F>;  #  num unitigs

    my $maxPos = 0;

    if ($utg =~ m/unitig\s+(\d+)$/) {
        $utg = $1;
    } else {
        die "Out of sync on '$utg'\n";
    }

    if ($nfrg =~ m/data.num_frags\s+(\d+)$/) {
        $nfrg = $1;
    } else {
        die "Out of sync on '$nfrg'\n";
    }

    my @frgs;

    for (my $nnnn=0; $nnnn < $nfrg; $nnnn++) {
        $_ = <F>;
        chomp;

        if (m/^FRG\stype\sR\sident\s+(\d+)\s+.*position\s+(\d+)\s+(\d+)$/) {
            push @frgs, "FRG\t$1\t$2\t$3";

            $maxPos = ($maxPos < $2) ? $2 : $maxPos;
            $maxPos = ($maxPos < $3) ? $3 : $maxPos;
        } else {
            die "Nope: '$_'\n";
        }
    }

    push @tigStore, "UTG\t$utg\t$nfrg\t$maxPos";
    push @tigStore, @frgs;
}
close(F);




if (defined($ARGV[0])) {
    push @reads, $ARGV[0];

} else {
    print STDERR "SCANNING tigStore for reads that span a unitig.\n";

    my $chk;
    my $utg;
    my $nfrg = 0;
    my $tfrg = 0;
    my $maxPos;

    foreach my $l (@tigStore) {
        if ($nfrg == 0) {
            ($chk, $utg, $nfrg, $maxPos) = split '\s+', $l;
            die if ($chk ne "UTG");

            next;
        }

        $nfrg--;
        $tfrg++;

        my ($chk, $iid, $bgn, $end) = split '\s+', $l;
        die if ($chk ne "FRG");

        if ($bgn > $end) {
            ($bgn, $end) = ($end, $bgn);
        }

        if (($bgn == 0) && ($end == $maxPos)) {
            push @reads, $iid;
        }
    }

    print STDERR "Searching for placements for ", scalar(@reads), " reads.\n";
}



foreach my $read (@reads) {
    my %ovlRead;
    my %ovlType;

    #print STDERR "Searching for placements for read $read in unitigs.\n";
    print "\n";

    #  Find overlapping reads.

    open(F, "overlapStore -b $read -e $read -d $prefix.ovlStore |");
    while (<F>) {
        s/^\s+//;
        s/\s+$//;

        my @v = split '\s+', $_;

        $ovlRead{$v[1]} = 1;
        $ovlType{$v[1]} = '';

        if      (($v[3] < 0) && ($v[4] < 0)) {
            $ovlType{$v[1]} = '5';
        } elsif (($v[3] > 0) && ($v[4] > 0)) {
            $ovlType{$v[1]} = '3';
        } elsif (($v[3] < 0) && ($v[4] > 0)) {
            $ovlType{$v[1]} = 'c';  #  Contained
        } else {
            $ovlType{$v[1]} = 'C';  #  Container
        }
    }
    close(F);

    #print STDERR "Searching for placements for read $read in unitigs, using ", scalar(keys %ovlRead), " overlapping reads.\n";

    #  Scan unitigs.  For any unitig with more than zero overlapping reads present, report.

    my $chk;
    my $utg;
    my $nfrg = 0;
    my $cfrg = 0;
    my $maxPos;

    my $nmat;
    my %ntyp;
    my $self;

    foreach my $l (@tigStore) {
        if ($cfrg == 0) {
            ($chk, $utg, $nfrg, $maxPos) = split '\s+', $l;
            die if ($chk ne "UTG");

            $cfrg      = $nfrg;

            $nmat      = 0;
            $ntyp{'5'} = 0;
            $ntyp{'3'} = 0;
            $ntyp{'C'} = 0;
            $ntyp{'c'} = 0;
            $self      = 0;

            next;
        }

        $cfrg--;

        my ($chk, $iid, $bgn, $end) = split '\s+', $l;
        die if ($chk ne "FRG");

        if ($bgn > $end) {
            ($bgn, $end) = ($end, $bgn);
        }

        if (exists($ovlRead{$iid})) {
            $nmat++;
            $ntyp{$ovlType{$iid}}++;

            #print "FRG\t$iid\tat\t$bgn\t$end\t$ovlType{$iid}\n";
        }

        if ($iid == $read) {
            $self++;
        }

        if (($cfrg == 0) && ($self + $nmat > 0)) {
            if ($self) {
                print "read $read\tUNITIG $utg\tof length $maxPos with 5'=$ntyp{'5'} 3'=$ntyp{'3'} contained=$ntyp{'c'} container=$ntyp{'C'} overlapping reads out of $nfrg.\n";
            } else {
                print "read $read\tunitig $utg\tof length $maxPos with 5'=$ntyp{'5'} 3'=$ntyp{'3'} contained=$ntyp{'c'} container=$ntyp{'C'} overlapping reads out of $nfrg.\n";
            }
        }
    }
    close(F);
}
