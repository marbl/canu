#!/usr/bin/perl

use strict;

#  Converts a FRGv1 or FRGv2 to fastq format.
#    o Preserves library name and insert size, but NO LIBRARY FEATURES.
#    o Preserves clear ranges.  Assumes ALL fragments in the input are in the same library.
#    o QVs are reencoded to the 'sanger' standard.
#    o Fragments are held in memory until the LKG record arrives; if no LKG record arrives,
#      the fragment is output as unmated.
#
#  DOES NOT CURRENTLY READ FRGv1 FRAGMENTS.
#

my %fg;

my $outPrefix = shift @ARGV;
my $outAname  = "$outPrefix.1.fastq";
my $outBname  = "$outPrefix.2.fastq";
my $outIname  = "$outPrefix.i.fastq";
my $outSname  = "$outPrefix.s.fastq";

open(outA, "> $outAname");
open(outB, "> $outBname");
open(outI, "> $outIname");
open(outS, "> $outSname");

my $v1 = 1;  #  LKG act: typ: fg1: fg2:
my $v2 = 0;  #  LKG act: frg: frg:

my $libName = undef;
my $libMean = 0.0;
my $libSDev = 0.0;


while (<STDIN>) {
    chomp;

    if (m/^ver:1$/) {
        $v1 = 1;
        $v2 = 0;
        next;
    }
    if (m/^ver:2$/) {
        $v1 = 0;
        $v2 = 1;
        next;
    }

    #  FRGv1 DISTANCE record OR FRGv2 LIBRARY record.

    if ((m/^{DST$/) ||
        (m/^{LIB$/)) {
        while ($_ ne "}") {
            $libName = $1  if (m/^acc:(.*)$/);
            $libMean = $1  if (m/^mea:(.*)$/);
            $libSDev = $1  if (m/^std:(.*)$/);
            
            $_ = <STDIN>;  chomp;
        }
    }


    if (m/^{FRG/) {
        my $uid = undef;
        my $seq = undef;
        my $qlt = undef;
        my $clr = "";

        if ($v1 == 1) {
            print STDERR "v1 fragments not supported.\n";
        }

        if ($v2 == 1) {
            $_      = <STDIN>;  chomp;  # act:
            my $acc = <STDIN>;  chomp;  # acc:
            if ($acc =~ m/^acc:(.*)$/) {
                $uid = $1;
            }

            my $rnd = <STDIN>;  chomp;  # rnd:
            if ($rnd =~ m/^rnd:0$/) {
                $clr .= " rnd=f";
            }

            $_      = <STDIN>;  chomp;  # sta:
            $_      = <STDIN>;  chomp;  # lib:
            $_      = <STDIN>;  chomp;  # pla:
            $_      = <STDIN>;  chomp;  # loc:
            $_      = <STDIN>;  chomp;  # src:

            do {
                $_  = <STDIN>;  chomp;
            } while ($_ ne ".");

            $_      = <STDIN>;  chomp;  # seq:
            $_      = "";
            do {
                $seq .= $_;
                $_  = <STDIN>;  chomp;
            } while ($_ ne ".");

            $_      = <STDIN>;  chomp;  # qlt:
            $_      = "";
            do {
                $qlt .= $_;
                $_  = <STDIN>;  chomp;
            } while ($_ ne ".");
                
            do {
                $_  = <STDIN>;  chomp;
            } while ($_ ne ".");

            $_      = <STDIN>;  chomp;  # clear?

            #  Order here is dictated by AS_MSG_pmesg2.c

            if (m/con:(\d+),(\d+)/) {
                $clr .= " tnt=$1,$2";
                $_    = <STDIN>;  chomp;
            }

            if (m/clv:(\d+),(\d+)/) {
                $clr .= " clv=$1,$2";
                $_    = <STDIN>;  chomp;
            }

            if (m/clq:(\d+),(\d+)/) {
                #$clr .= " max=$1,$2";  Obsolete FRGv2 clear range
                $_    = <STDIN>;  chomp;
            }

            if (m/clm:(\d+),(\d+)/) {
                $clr .= " max=$1,$2";
                $_    = <STDIN>;  chomp;
            }

            if (m/clr:(\d+),(\d+)/) {
                $clr .= " clr=$1,$2";
                $_    = <STDIN>;  chomp;
            }
        }

        #  Rebuild QVs.  The CA QV spec (AS_PER_encodeSequenceQuality.c) indicates that the max QV
        #  is 61.  The base is '0' (ascii 48).  The maximum is thus (ascii 48 + 61 = 109) 'm'.
        #
        #  In the conversion below we convert well past this maximum, all the way to 'z'.

        $qlt =~ tr/[0-z]/[!-k]/;

        $fg{$uid} = "\@$uid$clr\n$seq\n+\n$qlt\n";
    }



    if (m/^{LKG$/) {
        my $l1 = undef;
        my $l2 = undef;

        if ($v1 == 1) {
            $_  = <STDIN>;  #  act:
            $_  = <STDIN>;  #  typ:
            $l1 = <STDIN>;  #  fg1;
            $l2 = <STDIN>;  #  fg2;
        }
        if ($v2 == 1) {
            $_  = <STDIN>;  #  act:
            $l1 = <STDIN>;  #  frg;
            $l2 = <STDIN>;  #  frg;
        }

        $l1 = $1  if ($l1 =~ m/^...:(.*)$/);
        $l2 = $1  if ($l2 =~ m/^...:(.*)$/);
            
        die "No frg '$l1'\n"  if (!exists($fg{$l1}));
        die "No frg '$l2'\n"  if (!exists($fg{$l2}));

        print outA $fg{$l1};  print outI $fg{$l1};  delete $fg{$l1};
        print outB $fg{$l2};  print outI $fg{$l2};  delete $fg{$l2};
    }
}

foreach my $id (keys %fg) {
    print outS $fg{$id};
}

close(outA);
close(outB);
close(outI);
close(outS);

my $fqca;
$fqca .= "fastqToCA \\\n";
$fqca .= "  -insertsize $libMean $libSDev \\\n";
$fqca .= "  -libraryname $libName \\\n";
$fqca .= "  -technology sanger \\\n";
$fqca .= "  -type sanger \\\n";
$fqca .= "  -innie \\\n";
$fqca .= "  -mates $outAname,$outBname \\\n";
$fqca .= "  -reads $outSname \\\n";
$fqca .= "> $outPrefix.frg\n";

print STDERR "$fqca";
system($fqca);

exit(0);
