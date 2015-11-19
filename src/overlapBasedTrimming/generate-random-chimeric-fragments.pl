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

se strict;

#
###########################################################################
#
# This file is part of Celera Assembler, a software program that
# assembles whole-genome shotgun reads into contigs and scaffolds.
# Copyright (C) 2009, J. Craig Venter Instititue.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received (LICENSE.txt) a copy of the GNU General Public
# License along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
###########################################################################
#
#  Generates random chimera and spur fragments from a supplied reference
#  genome.
#
#  User unfriendly.
#
#  A spur:     [junk][sequence]
#  A chimer:   [junk][sequence][junk][sequence][junk]

my $genome     = shift @ARGV;
my $readLength = 400;

my $leaff = "/work/wgs/kmer/leaff/leaff";
my $cvtto = "perl /work/wgs/src/AS_MSG/convert-fasta-to-v2.pl";

die "Reference sequence '$genome' not found.\n" if (! -e $genome);

#  generate some random reads from the genome, both forward and reverse.

system("$leaff --errors $readLength 10000 1 0.01 $genome > sample3.fasta");
system("$leaff -f sample3.fasta -R -C -W > sample1.fasta");
system("$leaff --errors $readLength 10000 1 0.01 $genome > sample2.fasta");

#  generate random small bits of sequence from those samples; leaff can't do this nicely.

my @sequences;

open(F, "< sample1.fasta");
while (!eof(F)) {
    my $h = <F>;  chomp $h;
    my $s = <F>;  chomp $s;
    my $r = rand();
    push @sequences, "$r\0seq\0$s";
}
close(F);

open(F, "< sample2.fasta");
while (!eof(F)) {
    my $h = <F>;  chomp $h;
    my $s = <F>;  chomp $s;
    my $r = rand();
    push @sequences, "$r\0rev\0$s";
}
close(F);

@sequences = sort { $a <=> $b } @sequences;

open(F, "> chimer.fasta");
open(R, "$leaff -G 100000 $readLength $readLength |");

my $iid    = 0;
my $def    = "";
my $seq    = "";
my $seqLen = 0;

while (scalar(@sequences) > 0) {
    my $np = int(rand() * 3 + 1);
    my $gl;

    $iid++;
    $def    = "$iid";
    $seq    = "";
    $seqLen = 0;

    $gl = int(rand() * 200 - 50);  #  between -50 and 150, junk on the front
    if ($gl > 0) {
        my $h = <R>;
        my $s = <R>;

        $def    .= "-jnk-$gl";
        $seq    .= substr($s, 0, $gl);
        $seqLen += $gl;
    }

    for (my $pp=0; $pp < $np; $pp++) {
        my $S = pop @sequences;
        my ($r, $h, $s) = split '\0', $S;

        my $l = int(rand() * 200 + 50);  #  between 50 and 250.

        $def    .= "-$h-$l";
        $seq    .= substr($s, 0, $l);
        $seqLen += $l;

        if ($pp + 1 == $np) {
            $gl = 200;
        } else {
            $gl = 100;
        }

        $gl = int(rand() * $gl - 50);  #  between -50 and 50, junk in the middle
        if ($gl > 0) {
            my $h = <R>;
            my $s = <R>;

            $def    .= "-jnk-$gl";
            $seq    .= substr($s, 0, $gl);
            $seqLen += $gl;
        }
    }

    print F ">$def\n$seq\n";
}

close(F);

system("tr 'ACGT' 'NNNN' < chimer.fasta > chimer.qv");
system("perl /work/wgs/src/AS_MSG/convert-fasta-to-v2.pl -l CHIMER -s chimer.fasta -q chimer.qv > chimer.frg");

unlink "sample1.fasta";
unlink "sample2.fasta";
unlink "sample3.fasta";
unlink "sample1.fastaidx";
unlink "sample2.fastaidx";
unlink "sample3.fastaidx";
unlink "chimer.fasta";
unlink "chimer.qv";
