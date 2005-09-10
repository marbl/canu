#!/usr/local/bin/perl

#
#  bzip2 -dc /raid/WORK/EMpaper/run1-nofiltering/all-scored-hits.sorted.bz2 | perl dumpScores.pl | perl plotScoresSingly.pl
#

use strict;
$| = 1;

my $tmppath = "/tmp";

#  First line was blank.  Don't know why.
#my $junk = <STDIN>;

while (!eof(STDIN)) {
    $_ = <STDIN>;

    my ($estid, $maxscore, $numhits, @vals) = split '\s+', $_;

    open(A, "> $tmppath/hits-$estid.dat");
    open(B, "> $tmppath/iden-$estid.dat");
    open(C, "> $tmppath/covr-$estid.dat");
    foreach my $h (@vals) {
        my ($a, $b, $i, $c) = split ',', $h;

        $a /= $maxscore;
        print A "$a\n";

        $i /= 100.0;
        print B "$i\n";

        $c /= 100.0;
        print C "$c\n";
    }
    close(C);
    close(B);
    close(A);

    my $output = substr("0000000000$estid", -6, 6);
    my $direct = substr("0000000000$estid", -6, 3);

    print "$output\r";

    system("mkdir $direct") if (! -d "$direct");

    open(O, "> $tmppath/plot-$estid.gpl");
    print O "set terminal pbm color\n";
    print O "set output\n";
    print O "set pointsize 0.5\n";
    print O "set xtics 10\n";
    #print O "set size 1.5,1.5\n";
    print O "plot [-5:200][0.0:1.2] ";
    print O "              0.95 notitle lt 0, ";
    print O "              0.80 notitle lt 0, ";
    print O "              0.50 notitle lt 0, ";
    print O "              \"$tmppath/hits-$estid.dat\" using 1 notitle with linespoints 1, ";
    print O "              \"$tmppath/covr-$estid.dat\" using 1 notitle with points 3, ";
    print O "              \"$tmppath/iden-$estid.dat\" using 1 notitle with points 2\n";
    close(O);

    my $cmd = "";

    $cmd  = "gnuplot $tmppath/plot-$estid.gpl | ppmtogif -quiet > $direct/$output.gif";
    $cmd .= " && rm -f";
    $cmd .= "  $tmppath/plot-$estid.gpl";
    $cmd .= "  $tmppath/hits-$estid.dat";
    $cmd .= "  $tmppath/iden-$estid.dat";
    $cmd .= "  $tmppath/covr-$estid.dat";
    system("$cmd");
}
close(STDIN);
