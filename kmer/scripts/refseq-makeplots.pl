#!/bin/perl

#  Given a set of mappings of refseq to various genomes, plot the
#  number mapped at various percent identities over all coverages.
#  One plot for each identity level, that shows the number mapped at
#  all coverage levels.

use strict;

use lib "/work/assembly/walenzbp/projects/scripts";
use libBri;

my @datasets = ( "refseq-ncbi31/polishes-good",
                 "refseq-r27/polishes-good",
                 "refseq-ucsc-20000717/polishes-good",
                 "refseq-ucsc-20000905/polishes-good",
                 "refseq-ucsc-20001007/polishes-good",
                 "refseq-vanilla/polishes-good");

my %polishes;

foreach my $f (@datasets) {
    $polishes{$f} = ();

    if (! -e "$f.quality") {
        print STDERR "reading $f\n";

        open(F, "< $f");
        while (!eof(F)) {
            my %p = &libBri::readPolish(*F);
            push @{$polishes{$f}}, "$p{'percentID'} $p{'coverage'} $p{'estID'}";
            #print STDERR "$p{'percentID'} $p{'coverage'} $p{'estID'}\n";
        }
        close(F);

        open(F, "> $f.quality");
        $, = "\n";
        $\ = "\n";
        print F @{$polishes{$f}};
        $, = "";
        $\ = "";
        close(F);
    } else {
        open(F, "< $f.quality");
        @{$polishes{$f}} = <F>;
        close(F);
    }
}






#
#  Figure out which refseq are "old" or "new"
#
my %idmap;
my @refseqOld;
my @refseqNew;
my $id = 0;

my $refseqTotalAll = 0;
my $refseqTotalNew = 0;
my $refseqTotalOld = 0;

open(F, "grep '>' /dev5/walenz/TRANSCRIPTS/refseq_hum_mRNA_raw_all.fasta |");
while (<F>) {
    chomp;
    if (m/^>CRA\|(\d+)\s/) {
        $idmap{$1} = $id;
        $refseqOld[$id] = 0;
        $refseqNew[$id] = 1;
        $id++;
        $refseqTotalAll++;
        $refseqTotalNew++;
    } else {
        print STDERR "defline error in refseq_hum_mRNA_raw_all; '$_'\n";
    }
}
close(F);

open(F, "grep '>' /dev5/walenz/TRANSCRIPTS/refseq_hum_20010209.fasta |");
while (<F>) {
    chomp;
    if (m/^>CRA\|(\d+)\s/) {
        if (!defined($idmap{$1})) {
            $idmap{$1} = $id++;
        } else {
            $refseqTotalNew--;
        }
        $refseqOld[$idmap{$1}] = 1;
        $refseqNew[$idmap{$1}] = 0;
        $refseqTotalOld++;
    } else {
        print STDERR "defline error in refseq_hum_20010209; '$_'\n";
    }
}
close(F);

print STDERR "all: $refseqTotalAll (from file)\n";
print STDERR "old: $refseqTotalOld (from file)\n";
print STDERR "new: $refseqTotalNew (from file)\n";

my $j;
my $k;
$j = 0;
foreach my $x (@refseqOld) {
    $j++ if ($x == 1);
    $k++;
}
print STDERR "old: $j out of $k\n";

$j = 0;
$k = 0;
foreach my $x (@refseqNew) {
    $j++ if ($x == 1);
    $k++;
}
print STDERR "new: $j out of $k\n";



sub plotMappedAtVariousCoverage {
    my $style    = shift @_;

    for (my $identity=95; $identity<100; $identity+=5) {
        my $datafile = "refseq-coverage-$style-$identity";

        if (! -e "$datafile.dat") {
            open(F, "> $datafile.dat");
            for (my $coverage=50; $coverage <= 100; $coverage++) {
                print F "$coverage $identity  ";
                foreach my $f (@datasets) {
                    my $count = 0;
                    my %thing;
                    undef %thing;
                    foreach my $t (@{$polishes{$f}}) {
                        my ($i, $c, $id) = split '\s+', $t;


                        if (!defined($refseqNew[$id])) {
                            print STDERR "undefined new id - $id\n";
                        }
                        if (!defined($refseqOld[$id])) {
                            print STDERR "undefined old id - $id\n";
                        }




                        #  Count all refseq
                        if ($style eq "all") {
                            $thing{$id} = 1 if (($i >= $identity) && ($c >= $coverage));
                        }

                        #  Count only old refseq
                        if ($style eq "old") {
                            $thing{$id} = 1 if (($refseqOld[$id]) && ($i >= $identity) && ($c >= $coverage));
                        }

                        #  Count only new refseq
                        if ($style eq "new") {
                            $thing{$id} = 1 if (($refseqNew[$id]) && ($i >= $identity) && ($c >= $coverage));
                        }

                    }
                    $count = scalar(keys %thing);
                    print F "$count ";
                }
                print F "\n";
            }
            close(F);
        }

        open(F, "> $datafile.gnuplot");
        print F "set terminal postscript color solid\n";
        print F "set output \"$datafile.ps\"\n";

        print F "set title \"all refseq at $identity percent identity\n" if ($style eq "all");
        print F "set title \"old refseq at $identity percent identity\n" if ($style eq "old");
        print F "set title \"new refseq at $identity percent identity\n" if ($style eq "new");

        print F "set xlabel \"percent coverage\"\n";
        print F "set ylabel \"number of refseq\"\n";
        print F "plot [50:][0:19000] \\\n";
        print F "     $refseqTotalAll title \"Total ($refseqTotalAll)\", \\\n" if ($style eq "all");
        print F "     $refseqTotalNew title \"Total ($refseqTotalNew)\", \\\n" if ($style eq "new");
        print F "     $refseqTotalOld title \"Total ($refseqTotalOld)\", \\\n" if ($style eq "old");
        print F "     \"$datafile.dat\" using 1:3 title \"Mapped to ncbi31\" with lines, \\\n";
        print F "     \"$datafile.dat\" using 1:4 title \"Mapped to r27\" with lines, \\\n";
        print F "     \"$datafile.dat\" using 1:8 title \"Mapped to vanilla\" with lines, \\\n";
        print F "     \"$datafile.dat\" using 1:5 title \"Mapped to ucsc-20000717\" with lines, \\\n";
        print F "     \"$datafile.dat\" using 1:6 title \"Mapped to ucsc-20000905\" with lines, \\\n";
        print F "     \"$datafile.dat\" using 1:7 title \"Mapped to ucsc-20001007\" with lines\n";
        print F "set terminal pbm color small\n";
        print F "set output \"$datafile.ppm\"\n";
        print F "replot\n";
        close(F);

        system("gnuplot < $datafile.gnuplot");
        system("ppmtogif -quiet < $datafile.ppm > $datafile.gif");

        unlink("$datafile.gnuplot");
    }
}




plotMappedAtVariousCoverage("all");
plotMappedAtVariousCoverage("old");
plotMappedAtVariousCoverage("new");



exit;




#
#  Make the plots of genomeA vs genomeB; common, onlyA, onlyB, neither
#
#
#
#


my @datasets2 = @datasets;
foreach my $d1 (@datasets) {
    shift @datasets2;
    foreach my $d2 (@datasets2) {

        my $name1;
        if ($d1 =~ m/refseq-(.*)\/polishes-good/) {
            $name1 = $1;
        } else {
            print STDERR "Error parsing $d1 for name\n";
        }

        my $name2;
        if ($d2 =~ m/refseq-(.*)\/polishes-good/) {
            $name2 = $1;
        } else {
            print STDERR "Error parsing $d2 for name\n";
        }

        for (my $identity=95; $identity<=100; $identity++) {
            my $datafile = "refseq-pair-$identity-$name1-$name2";

            if (! -e "$datafile.dat") {
                print STDERR "Computing $datafile.dat\n";

                open(F, "> $datafile.dat");
                for (my $coverage=50; $coverage <= 100; $coverage++) {
                    print F "$coverage $identity  ";

                    my %thingd1;
                    my %thingd2;
                    my %thing;

                    foreach my $t (@{$polishes{$d1}}) {
                        my ($i, $c, $id) = split '\s+', $t;
                        if (($i >= $identity) && ($c >= $coverage)) {
                            $thingd1{$id} = 1;
                            $thing{$id} = 1;
                        }
                    }

                    foreach my $t (@{$polishes{$d2}}) {
                        my ($i, $c, $id) = split '\s+', $t;
                        if (($i >= $identity) && ($c >= $coverage)) {
                            $thingd2{$id} = 1;
                            $thing{$id} = 1;
                        }
                    }

                    my $both   = 0;
                    my $onlyd1 = 0;
                    my $onlyd2 = 0;

                    foreach my $t (keys %thing) {
                        my $a = defined($thingd1{$t});
                        my $b = defined($thingd2{$t});

                        if ($a && $b) {
                            $both++;
                        } elsif ($a) {
                            $onlyd1++;
                        } elsif ($b) {
                            $onlyd2++;
                        }
                    }

                    print F "$both $onlyd1 $onlyd2\n";
                }
                close(F);
            }

            open(F, "> $datafile.gnuplot");
            print F "set terminal postscript solid\n";
            print F "set output \"$datafile.ps\"\n";
            print F "set title \"refseq on $name1 and $name2 at >= $identity percent identity\n";
            print F "set xlabel \"percent coverage\"\n";
            print F "set ylabel \"number of refseq\"\n";
            print F "plot [50:][0:19000] \"$datafile.dat\" using 1:3 title \"Mapped to both\" with lines, \\\n";
            print F "     \"$datafile.dat\" using 1:4 title \"Mapped only to $name1\" with lines, \\\n";
            print F "     \"$datafile.dat\" using 1:5 title \"Mapped only to $name2\" with lines \n";
            print F "set terminal pbm small color\n";
            print F "set output \"$datafile.ppm\"\n";
            print F "replot\n";
            close(F);

            system("gnuplot < $datafile.gnuplot");
            system("ppmtogif -quiet < $datafile.ppm > $datafile.gif");

            unlink("$datafile.gnuplot");


            my $datafile95  = "refseq-pair-95-$name1-$name2";
            my $datafile96  = "refseq-pair-96-$name1-$name2";
            my $datafile97  = "refseq-pair-97-$name1-$name2";
            my $datafile98  = "refseq-pair-98-$name1-$name2";
            my $datafile99  = "refseq-pair-99-$name1-$name2";
            my $datafile100 = "refseq-pair-100-$name1-$name2";
            my $datafile    = "refseq-pair-$name1-$name2";

            open(F, "> $datafile.gnuplot");
            print F "set terminal postscript solid\n";
            print F "set output \"$datafile.ps\"\n";
            print F "set title \"refseq on $name1 and $name2 at >= $identity percent identity\n";
            print F "set xlabel \"percent coverage\"\n";
            print F "set ylabel \"number of refseq\"\n";
            print F "plot [50:][0:19000] \\\n";
            print F "\"$datafile95.dat\" using 1:3 title \"95 both\" with lines, \\\n";
            print F "\"$datafile95.dat\" using 1:4 title \"95 $name1\" with lines, \\\n";
            print F "\"$datafile95.dat\" using 1:5 title \"95 $name2\" with lines, \\\n";

            print F "\"$datafile96.dat\" using 1:3 title \"96 both\" with lines, \\\n";
            print F "\"$datafile96.dat\" using 1:4 title \"96 $name1\" with lines, \\\n";
            print F "\"$datafile96.dat\" using 1:5 title \"96 $name2\" with lines, \\\n";

            print F "\"$datafile97.dat\" using 1:3 title \"97 both\" with lines, \\\n";
            print F "\"$datafile97.dat\" using 1:4 title \"97 $name1\" with lines, \\\n";
            print F "\"$datafile97.dat\" using 1:5 title \"97 $name2\" with lines, \\\n";

            print F "\"$datafile98.dat\" using 1:3 title \"98 both\" with lines, \\\n";
            print F "\"$datafile98.dat\" using 1:4 title \"98 $name1\" with lines, \\\n";
            print F "\"$datafile98.dat\" using 1:5 title \"98 $name2\" with lines, \\\n";

            print F "\"$datafile99.dat\" using 1:3 title \"99 both\" with lines, \\\n";
            print F "\"$datafile99.dat\" using 1:4 title \"99 $name1\" with lines, \\\n";
            print F "\"$datafile99.dat\" using 1:5 title \"99 $name2\" with lines, \\\n";

            print F "\"$datafile100.dat\" using 1:3 title \"100 both\" with lines, \\\n";
            print F "\"$datafile100.dat\" using 1:4 title \"100 $name1\" with lines, \\\n";
            print F "\"$datafile100.dat\" using 1:5 title \"100 $name2\" with lines \n";

            print F "set terminal pbm small color\n";
            print F "set output \"$datafile.ppm\"\n";
            print F "replot\n";
            close(F);
            
            system("gnuplot < $datafile.gnuplot");
            system("ppmtogif -quiet < $datafile.ppm > $datafile.gif");

            unlink("$datafile.gnuplot");
        }
    }
}



