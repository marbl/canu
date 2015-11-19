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
 #    src/AS_GKP/fastqSimulate-perfectSep.pl
 #
 #  Modifications by:
 #
 #    Brian P. Walenz from 2011-DEC-28 to 2013-AUG-01
 #      are Copyright 2011,2013 J. Craig Venter Institute, and
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

#  Build MP, PE and 454 reads for a reference sequence.  The MP reads are filtered into
#  a perfectly classified set with two libraries PEperfect and MPperfect.  Chimers
#  can be rescued or deleted.

my $prefix     = shift @ARGV;
my $reference  = shift @ARGV;

my $error      = "0.01";

my $BBcoverage = "08";
my $PEcoverage = "25";  my $PEinsert = "0500";  my $PEstddev = "050";
my $MPcoverage = "10";  my $MPinsert = "3000";  my $MPstddev = "300";  my $MPenrich = "0.8";

my $BBreadLen = "400";
my $PEreadLen = "150";
my $MPreadLen = "150";

my $name;

#
#  BUILD 454
#

$name = "$prefix.${BBreadLen}bp.fragment.${BBcoverage}x";

if (! -e "$name.s.fastq") {
    system("fastqSimulate -f $reference -o $name -l $BBreadLen -x $BBcoverage -e $error -se");
    unlink("$name.frg");
}
if (! -e "$name.frg") {
    system("fastqToCA -libraryname BB -type sanger -reads $name.s.fastq > $name.frg");
}

#
#  BUILD PE
#

$name = "$prefix.${PEreadLen}bp.${PEinsert}bpPE.${PEcoverage}x";

if (! -e "$name.i.fastq") {
    system("fastqSimulate -f $reference -o $name -l $PEreadLen -x $PEcoverage -e $error -pe $PEinsert $PEstddev");
    unlink("$name.frg");
}
if (! -e "$name.frg") {
    system("fastqToCA -libraryname BB -insertsize $PEinsert $PEstddev -innie -type sanger -reads $name.1.fastq,$name.2.fastq > $name.frg");
}

#
#  BUILD MP
#

$name = "$prefix.${MPreadLen}bp.${MPinsert}bpMP.${MPcoverage}x";

if (! -e "$name.i.fastq") {
    print STDERR "fastqSimulate -f $reference -o $name -l $MPreadLen -x $MPcoverage -e $error -mp $MPinsert $MPstddev $PEinsert $PEstddev $MPenrich\n";
    system("fastqSimulate -f $reference -o $name -l $MPreadLen -x $MPcoverage -e $error -mp $MPinsert $MPstddev $PEinsert $PEstddev $MPenrich");
    unlink("$name.frg");
}
if (! -e "$name.frg") {
    system("fastqToCA -libraryname MP -insertsize $MPinsert $MPstddev -outtie -type sanger -reads $name.1.fastq,$name.2.fastq > $name.frg");
}

#
#  FILTER MP
#

open(tMP1, "> $name.tMP.1.fastq");
open(tMP2, "> $name.tMP.2.fastq");
open(tMPi, "> $name.tMP.i.fastq");

open(fPE1, "> $name.fPE.1.fastq");
open(fPE2, "> $name.fPE.2.fastq");
open(fPEi, "> $name.fPE.i.fastq");

open(aPE1, "> $name.aPE.1.fastq");  #  A-junction reads that are PE
open(aPE2, "> $name.aPE.2.fastq");
open(aPEi, "> $name.aPE.i.fastq");
open(aMP1, "> $name.aMP.1.fastq");  #  A-junction reads that are MP
open(aMP2, "> $name.aMP.2.fastq");
open(aMPi, "> $name.aMP.i.fastq");

open(bPE1, "> $name.bPE.1.fastq");  #  B-junction reads that are PE
open(bPE2, "> $name.bPE.2.fastq");
open(bPEi, "> $name.bPE.i.fastq");
open(bMP1, "> $name.bMP.1.fastq");  #  B-junction reads that are MP
open(bMP2, "> $name.bMP.2.fastq");
open(bMPi, "> $name.bMP.i.fastq");

open(MP1, "< $name.1.fastq");
open(MP2, "< $name.2.fastq");

my ($a1, $b1, $c1, $d1);
my ($a2, $b2, $c2, $d2);

while (!eof(MP1) && !eof(MP2)) {
    $a1 = <MP1>;  chomp $a1;  $a2 = <MP2>;  chomp $a2;  #  Painful.  The reverse()s below also reverse the newline.
    $b1 = <MP1>;  chomp $b1;  $b2 = <MP2>;  chomp $b2;
    $c1 = <MP1>;  chomp $c1;  $c2 = <MP2>;  chomp $c2;
    $d1 = <MP1>;  chomp $d1;  $d2 = <MP2>;  chomp $d2;

    if ($a1 =~ m/^.(...)_/) {
        my $typ = $1;

        if      ($typ eq "tMP") {
            $b1 = reverse($b1);  $b1 =~ tr/ACGTacgt/TGCAtgca/;  $d1 = reverse($d1);
            $b2 = reverse($b2);  $b2 =~ tr/ACGTacgt/TGCAtgca/;  $d2 = reverse($d2);

            print tMP1 "$a1\n$b1\n$c1\n$d1\n";
            print tMP2 "$a2\n$b2\n$c2\n$d2\n";
            print tMPi "$a1\n$b1\n$c1\n$d1\n";
            print tMPi "$a2\n$b2\n$c2\n$d2\n";

        } elsif ($typ eq "fPE") {
            print fPE1 "$a1\n$b1\n$c1\n$d1\n";
            print fPE2 "$a2\n$b2\n$c2\n$d2\n";
            print fPEi "$a1\n$b1\n$c1\n$d1\n";
            print fPEi "$a2\n$b2\n$c2\n$d2\n";

        } elsif ($typ eq "aMP") {
            if ($a1 =~ m/^...._\d+_\d+.\d+-\d+_(\d+)\/(\d+)\/\d+\/\d$/) {
                my $len = $2;
                my $pos = $1;

                if ($pos < $MPreadLen / 2) {
                    #  Saving the 5' end, makes MP.
                    $b1 = substr($b1, 0, $pos);  $b1 = reverse($b1);  $b1 =~ tr/ACGTacgt/TGCAtgca/;
                    $d1 = substr($d1, 0, $pos);  $d1 = reverse($d1);

                    $b2 = reverse($b2);  $b2 =~ tr/ACGTacgt/TGCAtgca/;
                    $d2 = reverse($d2);

                    print aMP1 "$a1\n$b1\n$c1\n$d1\n";
                    print aMP2 "$a2\n$b2\n$c2\n$d2\n";
                    print aMPi "$a1\n$b1\n$c1\n$d1\n";
                    print aMPi "$a2\n$b2\n$c2\n$d2\n";

                } else {
                    #  Saving the 3' end, makes PE.
                    $b1 = substr($b1, $pos);
                    $d1 = substr($d1, $pos);

                    print aPE1 "$a1\n$b1\n$c1\n$d1\n";
                    print aPE2 "$a2\n$b2\n$c2\n$d2\n";
                    print aPEi "$a1\n$b1\n$c1\n$d1\n";
                    print aPEi "$a2\n$b2\n$c2\n$d2\n";
                }

            } else {
                chomp $a1;
                die "No aMP '$a1'\n";
            }

        } elsif ($typ eq "bMP") {
            if ($a1 =~ m/^...._\d+_\d+.\d+-\d+_(\d+)\/(\d+)\/\d+\/\d$/) {
                my $pos = $1;
                my $len = $2;

                if ($pos < $MPreadLen / 2) {
                    #  Saving the 5' end, makes MP.
                    $b1 = reverse($b1);  $b1 =~ tr/ACGTacgt/TGCAtgca/;
                    $d1 = reverse($d1);

                    $b2 = substr($b2, 0, $pos);  $b2 = reverse($b2);  $b2 =~ tr/ACGTacgt/TGCAtgca/;
                    $d2 = substr($d2, 0, $pos);  $d2 = reverse($d2);

                    print bMP1 "$a1\n$b1\n$c1\n$d1\n";
                    print bMP2 "$a2\n$b2\n$c2\n$d2\n";
                    print bMPi "$a1\n$b1\n$c1\n$d1\n";
                    print bMPi "$a2\n$b2\n$c2\n$d2\n";

                } else {
                    #  Saving the 3' end, makes PE.
                    $b2 = substr($b2, $pos);
                    $d2 = substr($d2, $pos);

                    print bPE1 "$a1\n$b1\n$c1\n$d1\n";
                    print bPE2 "$a2\n$b2\n$c2\n$d2\n";
                    print bPEi "$a1\n$b1\n$c1\n$d1\n";
                    print bPEi "$a2\n$b2\n$c2\n$d2\n";
                }

            } else {
                chomp $a1;
                die "No bMP '$a1'\n";
            }


        } elsif ($typ eq "cMP") {
            chomp;
            print STDERR "$a1\n";

        } else {
            die "Unknown type '$typ'\n";
        }
    } else {
        chomp $a1;
        die "no match: '$a1'\n";
    }
}

close(tMP1);
close(tMP2);
close(tMPi);

close(fPE1);
close(fPE2);
close(fPEi);

close(aMP1);
close(aMP2);
close(aMPi);
close(aPE1);
close(aPE2);
close(aPEi);

close(bMP1);
close(bMP2);
close(bMPi);
close(bPE1);
close(bPE2);
close(bPEi);


system("fastqToCA -libraryname tMP -insertsize $MPinsert $MPstddev -innie -type sanger -mates $name.tMP.1.fastq,$name.tMP.2.fastq > $name.tMP.frg");
system("fastqToCA -libraryname fPE -insertsize $PEinsert $PEstddev -innie -type sanger -mates $name.fPE.1.fastq,$name.fPE.2.fastq > $name.fPE.frg");
system("fastqToCA -libraryname aMP -insertsize $MPinsert $MPstddev -innie -type sanger -mates $name.aMP.1.fastq,$name.aMP.2.fastq > $name.aMP.frg");
system("fastqToCA -libraryname aPE -insertsize $PEinsert $PEstddev -innie -type sanger -mates $name.aPE.1.fastq,$name.aPE.2.fastq > $name.aPE.frg");
system("fastqToCA -libraryname bMP -insertsize $MPinsert $MPstddev -innie -type sanger -mates $name.bMP.1.fastq,$name.bMP.2.fastq > $name.bMP.frg");
system("fastqToCA -libraryname bPE -insertsize $PEinsert $PEstddev -innie -type sanger -mates $name.bPE.1.fastq,$name.bPE.2.fastq > $name.bPE.frg");


system("perl /work/FRAGS/map-illumina-pairs.pl plot-tMP JUNK.150bp.3000bpMP.10x.tMP.1.fastq JUNK.150bp.3000bpMP.10x.tMP.2.fastq Yersinia_pestis_KIM10_AE009952.fasta");
system("perl /work/FRAGS/map-illumina-pairs.pl plot-fPE JUNK.150bp.3000bpMP.10x.fPE.1.fastq JUNK.150bp.3000bpMP.10x.fPE.2.fastq Yersinia_pestis_KIM10_AE009952.fasta");
system("perl /work/FRAGS/map-illumina-pairs.pl plot-aMP JUNK.150bp.3000bpMP.10x aMP.1.fastq JUNK.150bp.3000bpMP.10x.aMP.2.fastq Yersinia_pestis_KIM10_AE009952.fasta");
system("perl /work/FRAGS/map-illumina-pairs.pl plot-aPE JUNK.150bp.3000bpMP.10x.aPE.1.fastq JUNK.150bp.3000bpMP.10x.aPE.2.fastq Yersinia_pestis_KIM10_AE009952.fasta");
system("perl /work/FRAGS/map-illumina-pairs.pl plot-bMP JUNK.150bp.3000bpMP.10x.bMP.1.fastq JUNK.150bp.3000bpMP.10x.bMP.2.fastq Yersinia_pestis_KIM10_AE009952.fasta");
system("perl /work/FRAGS/map-illumina-pairs.pl plot-bPE JUNK.150bp.3000bpMP.10x.bPE.1.fastq JUNK.150bp.3000bpMP.10x.bPE.2.fastq Yersinia_pestis_KIM10_AE009952.fasta");
