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
 #  This file is derived from:
 #
 #    src/AS_GKP/fastqSimulate-checkCoverage.pl
 #
 #  Modifications by:
 #
 #    Brian P. Walenz on 2014-FEB-20
 #      are Copyright 2014 J. Craig Venter Institute, and
 #      are subject to the GNU General Public License version 2
 #
 #    Brian P. Walenz on 2015-JAN-13
 #      are Copyright 2015 Battelle National Biodefense Institute, and
 #      are subject to the BSD 3-Clause License
 #
 #    Brian P. Walenz beginning on 2015-OCT-12
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

use strict;

system("leaff    -G 1 2000 2000                 > base.fasta");
system("leaff -H -G 1  500  500 | tr ACGT NNNN >> base.fasta");
system("leaff -H -G 1 2000 2000                >> base.fasta");

my (@rcov1, @ccov1);
my (@rcov2, @ccov2);
my (@rcov3, @ccov3);

if (! -e "cov") {
    my $max = -s "base.fasta";
    my $maxiter = 100;

    for (my $iter=0; $iter<$maxiter; $iter++) {
        system("fastqSimulate -em 0 -ei 0 -ed 0 -f base.fasta -o t1 -l 100 -x 100                     -pe 1000 0 > /dev/null 2>&1");
        system("fastqSimulate -em 0 -ei 0 -ed 0 -f base.fasta -o t2 -l 100 -x 100 -allowgaps          -pe 1000 0 > /dev/null 2>&1");
        system("fastqSimulate -em 0 -ei 0 -ed 0 -f base.fasta -o t3 -l 100 -x 100 -allowgaps -allowns -pe 1000 0 > /dev/null 2>&1");

        system("rm *.1.fastq *.2.fastq *.c.fastq");

        print STDERR "Summarizing for iter $iter\n";

        open(IN, "< t1.i.fastq");  summarize(\@rcov1, \@ccov1);  close(IN);
        open(IN, "< t2.i.fastq");  summarize(\@rcov2, \@ccov2);  close(IN);
        open(IN, "< t3.i.fastq");  summarize(\@rcov3, \@ccov3);  close(IN);
    }


    open(COV, "> cov");

    for (my $x=0; $x<$max; $x++) {
        $rcov1[$x] += 0;  $rcov1[$x] /= $maxiter;
        $ccov1[$x] += 0;  $ccov1[$x] /= $maxiter;

        $rcov2[$x] += 0;  $rcov2[$x] /= $maxiter;
        $ccov2[$x] += 0;  $ccov2[$x] /= $maxiter;

        $rcov3[$x] += 0;  $rcov3[$x] /= $maxiter;
        $ccov3[$x] += 0;  $ccov3[$x] /= $maxiter;

        print COV "$x\t$rcov1[$x]\t$ccov1[$x]\t$rcov2[$x]\t$ccov2[$x]\t$rcov3[$x]\t$ccov3[$x]\n";
    }

    close(COV);
}


open(GP, "> cov.gp");
print GP "set terminal png size 1280,800\n";
print GP "set output 'cov.png'\n";
print GP "plot 'cov' using 1:2 with lines lt 1 lw 1 title 'read cov, no Ns', \\\n";
print GP "     'cov' using 1:4 with lines lt 2 lw 1 title 'read cov, clone span Ns', \\\n";
print GP "     'cov' using 1:6 with lines lt 3 lw 1 title 'read cov, both span Ns', \\\n";
print GP "     'cov' using 1:3 with lines lt 1 lw 2 title 'clone cov, no Ns', \\\n";
print GP "     'cov' using 1:5 with lines lt 2 lw 2 title 'clone cov, clone span Ns', \\\n";
print GP "     'cov' using 1:7 with lines lt 3 lw 2 title 'clone cov, both span Ns'\n";
close(GP);

system("gnuplot cov.gp");






sub summarize(@@) {
    my ($rcov, $ccov) = @_;

    while (!eof(IN)) {
        my $a = <IN>;  chomp $a;
        my $b = <IN>;  chomp $b;  my $l = length($b);
        my $c = <IN>;
        my $d = <IN>;

        my $bgn;
        my $end;
        my $spn;

        if ($a =~ m/^\@PE_\d+_\d+\@(\d+)-(\d+)\/1/) {
            $bgn = $1;
            $end = $1 + $l;
            $spn = $2;
        }

        if ($a =~ m/^\@PE_\d+_\d+\@(\d+)-(\d+)\/2/) {
            $bgn = $2 - $l;
            $end = $2;
            $spn = $bgn;  #  No clone span for the second read
        }

        die if (!defined($bgn));

        for (my $x=$bgn; $x<$end; $x++) {
            ${$rcov}[$x]++;
        }

        for (my $x=$bgn; $x<$spn; $x++) {
            ${$ccov}[$x]++;
        }
    }
}
