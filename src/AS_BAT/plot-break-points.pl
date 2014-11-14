#!/usr/bin/perl

use strict;

my %nucmer;
my %length;

if (! -e "tig.fasta") {
}

if (! -e "tig.coords") {
}

if (-e "tig.fasta") {
    open(F, "< tig.fasta") or die;
    while (<F>) {
        if (m/^>(utg\d+)\s+len=(\d+)$/) {
            $length{$1} = $2;
        }
    }
    close(F);
} else {
    print STDERR "No tig.fasta file found, cannot .....\n";
}

if (-e "tig.coords") {
    open(F, "< tig.coords") or die;
    $_ = <F>;
    $_ = <F>;
    $_ = <F>;
    $_ = <F>;
    $_ = <F>;
    while (<F>) {
        s/^\s+//;
        s/\s$+//;
        my @v = split '\s+', $_;
        my $str = "$v[3]\t$v[4]";
        if (exists($nucmer{$v[12]})) {
            $nucmer{$v[12]} .= "\t" . $str;
        } else {
            $nucmer{$v[12]} = $str;
            #print "$v[12] -- $str\n";
        }
    }
    close(F);
} else {
    print STDERR "No tig.coords file found, not plotting true break points.\n";
}


open(F, "find coverageplot -name '*.badCoverage' -print | sort |");
my @files = <F>;
chomp @files;
close(F);



#undef @files;
#push @files, "utg000002841.badCoverage";
#push @files, "utg000008002.badCoverage";
#push @files, "utg000008093.badCoverage";

foreach my $file (@files) {
    my $num = 0;

    if ($file =~ m/utg0+(\d+).badCoverage/) {
        $num = $1;
    } else {
        print STDERR "No match '$file'\n";
        exit(1);
    }

    open(GP, "> $file.gnuplot");

    print GP "set terminal png size 1280,800\n";
    print GP "set output \"$file.png\"\n";

    my $cnt = 1;

    open(F, "cat unitigger.*.intersectionBreaking.log | ") or die;
    while (<F>) {
        if (m/BREAK unitig (\d+) at position (\d+),(\d+) from inSize (\d+) inFrags (\d+)./) {
            if ($1 == $num) {
                my $s = int(-$4 / 100);  #  Negative: length / 100
                my $f = int( $5 / 10);   #  Positive: num frags / 10
                #print GP "set object rect from $2,$s to $3,$f fc lt 2\n";
                print GP "set object rect from $2,$s to $3,$f fc lt $cnt\n";
                $cnt++;
            }
        }
        if (m/BREAK unitig (\d+) at position (\d+),(\d+) from MATES./) {
            if ($1 == $num) {
                my $s = -0.5;
                my $f =  0.5;
                print GP "set object rect from $2,$s to $3,$f fc lt 0\n";
                $cnt++;
            }
        }
    }
    close(F);

    if (exists($nucmer{"utg$1"})) {
        my $positionB = 5;
        my $positionT = 5.5;
        my $nucmer   = $nucmer{"utg$1"};
        my @L        = split '\t', $nucmer;

        #print "utg$1 -- $nucmer\n";

        delete $nucmer{"utg$1"};

        while (scalar(@L) > 1) {
            my $b = shift @L;
            my $e = shift @L;

            print GP "set object rect from $b,$positionB to $e,$positionT fc lt 0\n";
            $positionB++;
            $positionT++;
        }
    }


    print GP "plot \\\n";
    #print GP "     \"$file\" using 1:3 with lines title \"fwdBad\", \\\n";
    #print GP "     \"$file\" using 1:4 with lines title \"revBad\", \\\n";
    print GP "     \"$file\" using 1:5 with lines title \"badExtFwd\", \\\n";
    print GP "     \"$file\" using 1:6 with lines title \"badExtRev\", \\\n";
    print GP "     \"$file\" using 1:7 with lines title \"badCompressed\", \\\n";
    print GP "     \"$file\" using 1:8 with lines title \"badStretched\", \\\n";
    print GP "     \"$file\" using 1:9 with lines title \"badNormal\", \\\n";
    print GP "     \"$file\" using 1:10 with lines title \"badAnti\", \\\n";
    print GP "     \"$file\" using 1:11 with lines title \"badOuttie\", \\\n";
    print GP "     \"$file\" using 1:2 with lines title \"Good\"\n";
    print GP "\n";

    close(GP);

    print STDERR "$file\n";
    system("gnuplot < $file.gnuplot > /dev/null 2>&1");
}


foreach my $k (keys %nucmer) {
    my $nucmer   = $nucmer{$k};
    my @L        = split '\t', $nucmer;

    my $c        = 0;
    my $t        = 0;

    while (scalar(@L) > 0) {
        my $b = shift @L;
        my $e = shift @L;
        my $l = ($b < $e) ? ($e - $b) : ($b - $e);

        $t++;
        $c++ if (($l < $length{$k} - 10) && ($l > 150));
    }

    if ($c > 1) {
        print STDERR "$k with $c partial out of $t total regions.  length $length{$k}.\n";
    }
}
