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

#  Reads true overlaps from the first file, separates the second file into
#  true/false based on iid.  Hangs are NOT checked.

my $trueOverlaps = "BL/test.ovlStore";
my $testOverlaps = "CAordered/test.ovlStore";

my $filtTrue  = "filtered.true";
my $filtFalse = "filtered.false";

#
#  Load read lengths, make sure both assemblies are the same
#

my @readLengths;
my %contained;

open(F, "gatekeeper -dumpfragments -tabular BL/test.gkpStore |");
while (<F>) {
    my @v = split '\s+', $_;

    $readLengths[$v[1]] = $v[9];
}
close(F);

open(F, "gatekeeper -dumpfragments -tabular CAordered/test.gkpStore |");
while (<F>) {
    my @v = split '\s+', $_;

    die "IID $v[1] BL $readLengths[$v[1]] != CA $v[9]\n"   if ($readLengths[$v[1]] != $v[9]);
}
close(F);

open(F, "< CAordered/4-unitigger/best.contains");
while (<F>) {
    my @v = split '\s+', $_;
    $contained{$v[0]}++;
}
close(F);



#
#  Load true overlaps
#

my %true;           #  If set, id pair is a true overlap
my %trueComputed;   #  If set, id pair is a true overlap, and was computed

{
    my $last;       #  local storage
    my @lenIdent;   #  local storage of length-vs-ident for true overlaps.  ident is not known.

    open(F, "overlapStore -d $trueOverlaps |") or die "Failed to open '$trueOverlaps' for reading.\n";
    while (<F>) {
        s/^\s+//;
        s/\s+$//;

        my ($aIID, $bIID, $orient, $aHang, $bHang, $origErate, $corrErate) = split '\s+', $_;

        if ($last != $aIID) {
            open(O, "> iid$last.len-vs-ident.dat");
            print O @lenIdent;
            close(O);

            undef @lenIdent;
        }

        my $oLen = computeOverlapLength($aIID, $bIID, $aHang, $bHang);

        $true{"$aIID-$bIID"} = $oLen;
        $true{"$bIID-$aIID"} = $oLen;

        next if (exists($contained{$aIID}));
        next if (exists($contained{$bIID}));

        push @lenIdent, "$oLen\t$12.0\n";

        $last = $aIID;
    }
    close(F);
}

print STDERR "true:  ", scalar(keys %true) / 2, "\n";




#
#  Process computed overlaps
#

my %trueH;
my %falseH;
my %keysH;

my @lenIdent5T;
my @lenIdent5F;

my @lenIdent3T;
my @lenIdent3F;

my $longestF = 0;
my $longestT = 0;


{
my $last  = 0;

open(LT, "> all.true.length-vs-ident.dat");
open(LF, "> all.false.length-vs-ident.dat");

open(TO, "> all.true.filtered-overlaps.ova");
open(FO, "> all.false.filtered-overlaps.ova");

open(F, "overlapStore -d $testOverlaps |") or die "Failed to open '$testOverlaps' for reading.\n";
while (<F>) {
    s/^\s+//;
    s/\s+$//;

    my ($aIID, $bIID, $orient, $aHang, $bHang, $origErate, $corrErate) = split '\s+', $_;

    if ($last != $aIID) {
        if (scalar(keys %falseH) > 0) {

            #  error rate histogram
            if (1) {
                open(O, "> iid$last.eratehistogram.dat");
                foreach my $ii (sort { $a <=> $b } keys %keysH) {
                    $trueH{$ii}  += 0.0;
                    $falseH{$ii} += 0.0;
                    print O "$ii\t$trueH{$ii}\t$falseH{$ii}\n";
                }
                close(O);

                open(O, "| gnuplot > /dev/null 2> /dev/null");
                print O "set terminal 'png'\n";
                print O "set output   'iid$last.eratehistogram.png'\n";
                print O "plot 'iid$last.eratehistogram.dat' using 1:3 with lines title 'FALSE',";
                print O "     'iid$last.eratehistogram.dat' using 1:2 with lines title 'TRUE'\n";
                close(O);
            }


            #  length vs identitiy
            if (1) { #0.5 * $longestT < $longestF) {
                #print STDERR "IID $last  true: len=$longestT   false: len=$longestF\n";

                open(O, "> iid$last.len-vs-ident.5.true.dat");
                print O @lenIdent5T;
                close(O);

                open(O, "> iid$last.len-vs-ident.5.false.dat");
                print O @lenIdent5F;
                close(O);

                open(O, "> iid$last.len-vs-ident.3.true.dat");
                print O @lenIdent3T;
                close(O);

                open(O, "> iid$last.len-vs-ident.3.false.dat");
                print O @lenIdent3F;
                close(O);

                open(O, "| gnuplot > /dev/null 2> /dev/null");
                print O "set terminal 'png'\n";
                print O "set output   'iid$last.len-vs-ident.png'\n";
                print O "plot [10:35] [0:15000]";
                print O "   'iid$last.len-vs-ident.5.false.dat' using 2:1 title 'FALSE 5' pt 1 lc 1,";
                print O "   'iid$last.len-vs-ident.5.true.dat'  using 2:1 title 'TRUE 5'  pt 1 lc 2,";
                print O "   'iid$last.len-vs-ident.3.false.dat' using 2:1 title 'FALSE 3' pt 2 lc 1,";
                print O "   'iid$last.len-vs-ident.3.true.dat'  using 2:1 title 'TRUE 3'  pt 2 lc 2,";
                print O "   'iid$last.len-vs-ident.dat'         using 2:1 title 'MISSED'  pt 3 lc 3\n";

                print O "set output   'iid$last.5.len-vs-ident.png'\n";
                print O "plot [10:35] [0:15000]";
                print O "   'iid$last.len-vs-ident.5.false.dat' using 2:1 title 'FALSE 5' pt 1 lc 1,";
                print O "   'iid$last.len-vs-ident.5.true.dat'  using 2:1 title 'TRUE 5'  pt 1 lc 2\n";

                print O "set output   'iid$last.3.len-vs-ident.png'\n";
                print O "plot [10:35] [0:15000]";
                print O "   'iid$last.len-vs-ident.3.false.dat' using 2:1 title 'FALSE 3' pt 2 lc 1,";
                print O "   'iid$last.len-vs-ident.3.true.dat'  using 2:1 title 'TRUE 3'  pt 2 lc 2\n";
                close(O);
            }
        }

        undef %trueH;
        undef %falseH;

        undef %keysH;

        undef @lenIdent5T;
        undef @lenIdent5F;

        undef @lenIdent3T;
        undef @lenIdent3F;

        $longestT = 0;
        $longestF = 0;
    }

    if (exists($true{"$aIID-$bIID"})) {
        print TO "$_\n";
    } else {
        print FO "$_\n";
    }

    $trueComputed{"$aIID-$bIID"}++;

    next if (exists($contained{$aIID}));
    next if (exists($contained{$bIID}));

    my $oLen = computeOverlapLength($aIID, $bIID, $aHang, $bHang);

    my $intErate = int($origErate);

    my $is5 = ($aHang < 0);

    #  Computing the number of matches instead of actual length makes every overlap shorter
    #  than truth, and they all get flagged in output
    #
    #$oLen -= int($oLen * $origErate / 100);

    if (exists($true{"$aIID-$bIID"})) {
        $trueH{$intErate}++;
        $keysH{$intErate}++;

        print LT         "$oLen\t$origErate\n";

        if ($is5) {
            push @lenIdent5T, "$oLen\t$origErate\n";
        } else {
            push @lenIdent3T, "$oLen\t$origErate\n";
        }

        if ($longestT < $oLen) {
            $longestT = $oLen;
        }

    } else {
        $falseH{$intErate}++;
        $keysH{$intErate}++;

        print LF         "$oLen\t$origErate\n";

        if ($is5) {
            push @lenIdent5F, "$oLen\t$origErate\n";
        } else {
            push @lenIdent3F, "$oLen\t$origErate\n";
        }

        if ($longestF < $oLen) {
            $longestF = $oLen;
        }
    }

    $last = $aIID;
}
close(F);

close(LT);
close(LF);

close(TO);
close(FO);
}


print "convertOverlap -ovl < all.true.filtered-overlaps.ova > all.true.filtered-overlaps.ovb\n";
system("convertOverlap -ovl < all.true.filtered-overlaps.ova > all.true.filtered-overlaps.ovb");

print "convertOverlap -ovl < all.false.filtered-overlaps.ova > all.false.filtered-overlaps.ovb\n";
system("convertOverlap -ovl < all.false.filtered-overlaps.ova > all.false.filtered-overlaps.ovb");

print "rm -rf all.true.filtered-overlaps.ovlStore all.true.filtered-overlaps+false.ovlStore\n";
system("rm -rf all.true.filtered-overlaps.ovlStore all.true.filtered-overlaps+false.ovlStore");

print "overlapStoreBuild -g CAordered/test.gkpStore -o all.true.filtered-overlaps.ovlStore -F 1 all.true.filtered-overlaps.ovb\n";
system("overlapStoreBuild -g CAordered/test.gkpStore -o all.true.filtered-overlaps.ovlStore -F 1 all.true.filtered-overlaps.ovb");

print "overlapStoreBuild -g CAordered/test.gkpStore -o all.true.filtered-overlaps+false.ovlStore -F 1 all.true.filtered-overlaps.ovb all.false.filtered-overlaps.ovb\n";
system("overlapStoreBuild -g CAordered/test.gkpStore -o all.true.filtered-overlaps+false.ovlStore -F 1 all.true.filtered-overlaps.ovb all.false.filtered-overlaps.ovb");








my %bestTrueEdge;
my %bestTrueContain;

open(F, "< BL/4-unitigger/best.edges") or die;
while (<F>) {
    my @v = split '\s+', $_;
    $bestTrueEdge{"$v[0]-$v[2]"}++;
    $bestTrueEdge{"$v[2]-$v[0]"}++;

    $bestTrueEdge{"$v[0]-$v[4]"}++;
    $bestTrueEdge{"$v[4]-$v[0]"}++;
}
close(F);

open(F, "< BL/4-unitigger/best.contains") or die;
while (<F>) {
    my @v = split '\s+', $_;
    $bestTrueContain{"$v[0]-$v[3]"}++;
}
close(F);



my %bestTestEdge;
my %bestTestContain;

open(F, "< CAordered/4-unitigger/best.edges") or die;
while (<F>) {
    my @v = split '\s+', $_;
    $bestTestEdge{"$v[0]-$v[2]"}++;
    $bestTestEdge{"$v[2]-$v[0]"}++;

    $bestTestEdge{"$v[0]-$v[4]"}++;
    $bestTestEdge{"$v[4]-$v[0]"}++;
}
close(F);

open(F, "< CAordered/4-unitigger/best.contains") or die;
while (<F>) {
    my @v = split '\s+', $_;
    $bestTestContain{"$v[0]-$v[3]"}++;
}
close(F);


#
#  Output true overlaps that we missed
#

open(F, "overlapStore -d $trueOverlaps |") or die "Failed to open '$trueOverlaps' for reading.\n";

open(O, "> true-overlaps.missed.other.ova");
open(E, "> true-overlaps.missed.edge.ova");
open(C, "> true-overlaps.missed.contain.ova");

while (<F>) {
    s/^\s+//;
    s/\s+$//;

    my ($aIID, $bIID, $orient, $aHang, $bHang, $origErate, $corrErate) = split '\s+', $_;

    next  if (exists($trueComputed{"$aIID-$bIID"}));

    if (exists($bestTrueEdge{"$aIID-$bIID"})) {
        print E "$_\n";

    } elsif (exists($bestTrueContain{"$aIID-$bIID"})) {
        print C "$_\n";

    } else {
        print O "$_\n";
    }
}
close(O);
close(F);


#
#  Classify computed overlaps into true/false for each type (best edge, best contain, the others).
#

open(F, "overlapStore -d $testOverlaps |") or die "Failed to open '$trueOverlaps' for reading.\n";

open(BET, "> computed-overlaps.bestE.true.dat");
open(BEN, "> computed-overlaps.bestE.near.dat");
open(BEF, "> computed-overlaps.bestE.false.dat");

open(BCT, "> computed-overlaps.bestC.true.dat");
open(BCN, "> computed-overlaps.bestC.near.dat");
open(BCF, "> computed-overlaps.bestC.false.dat");

open(OTT, "> computed-overlaps.other.true.dat");
open(OTF, "> computed-overlaps.other.false.dat");

open(BETo, "> computed-overlaps.bestE.true.ova");
open(BENo, "> computed-overlaps.bestE.near.ova");
open(BEFo, "> computed-overlaps.bestE.false.ova");

open(BCTo, "> computed-overlaps.bestC.true.ova");
open(BCNo, "> computed-overlaps.bestC.near.ova");
open(BCFo, "> computed-overlaps.bestC.false.ova");

open(OTTo, "> computed-overlaps.other.true.ova");
open(OTFo, "> computed-overlaps.other.false.ova");

while (<F>) {
    s/^\s+//;
    s/\s+$//;

    my ($aIID, $bIID, $orient, $aHang, $bHang, $origErate, $corrErate) = split '\s+', $_;

    my $pair = "$aIID-$bIID";
    my $len  = computeOverlapLength($aIID, $bIID, $aHang, $bHang);

    $origErate = int(10 * $origErate) / 10;

    my $diff = ($aIID < $bIID) ? $bIID - $aIID : $aIID - $bIID;

    die if ($diff < 0);

    if (exists($bestTestEdge{$pair})) {
        if (exists($bestTrueEdge{$pair})) {
            print BET  "$origErate\t$len\n";
            print BETo "$_\n";
        } elsif (exists($true{$pair})) {
            print BEN  "$origErate\t$len\n";
            print BENo "$_\n";
        } else {
            print BEF  "$origErate\t$len\n";
            print BEFo "$_\n";
        }

    } elsif (exists($bestTestContain{$pair})) {
        if (exists($bestTrueContain{$pair})) {
            print BCT  "$origErate\t$len\n";
            print BCTo "$_\n";
        } elsif (exists($true{$pair})) {
            print BCN  "$origErate\t$len\n";
            print BCNo "$_\n";
        } else {
            print BCF  "$origErate\t$len\n";
            print BCFo "$_\n";
        }

    } else {
        if (exists($true{$pair})) {
            print OTT  "$origErate\t$len\n";
            print OTTo "$_\n";
        } else {
            print OTF  "$origErate\t$len\n";
            print OTFo "$_\n";
        }
    }
}

close(OTF);  close(OTFo);
close(OTT);  close(OTTo);

close(BCF);  close(BCFo);
close(BCN);  close(BCNo);
close(BCT);  close(BCTo);

close(BEF);  close(BEFo);
close(BEN);  close(BENo);
close(BET);  close(BETo);

system("awk '{ print \$1 }' computed-overlaps.bestE.true.dat | sort -n | uniq -c > computed-overlaps.bestE.true.eratehist");
system("awk '{ print \$1 }' computed-overlaps.bestE.near.dat | sort -n | uniq -c > computed-overlaps.bestE.near.eratehist");
system("awk '{ print \$1 }' computed-overlaps.bestE.false.dat | sort -n | uniq -c > computed-overlaps.bestE.false.eratehist");

system("awk '{ print \$1 }' computed-overlaps.bestC.true.dat | sort -n | uniq -c > computed-overlaps.bestC.true.eratehist");
system("awk '{ print \$1 }' computed-overlaps.bestC.near.dat | sort -n | uniq -c > computed-overlaps.bestC.near.eratehist");
system("awk '{ print \$1 }' computed-overlaps.bestC.false.dat | sort -n | uniq -c > computed-overlaps.bestC.false.eratehist");

system("awk '{ print \$1 }' computed-overlaps.other.true.dat | sort -n | uniq -c > computed-overlaps.other.true.eratehist");
system("awk '{ print \$1 }' computed-overlaps.other.false.dat | sort -n | uniq -c > computed-overlaps.other.false.eratehist");



open(F, "> tmp.gp");

print F "set terminal png\n";
print F "\n";
print F "set output 'computed-overlaps.bestE.eratehist.png'\n";
print F "set title 'Computed Overlaps, best edges, error rate histogram'\n";
print F "plot 'computed-overlaps.bestE.true.eratehist'  using 2:1 title 'true'  ps 1.00 lc 2, \\\n";
print F "     'computed-overlaps.bestE.near.eratehist'  using 2:1 title 'near'  ps 1.00 lc 3, \\\n";
print F "     'computed-overlaps.bestE.false.eratehist' using 2:1 title 'false' ps 1.00 lc 1\n";
print F "\n";
print F "set title 'Computed Overlaps, best contains, error rate histogram'\n";
print F "set output 'computed-overlaps.bestC.eratehist.png'\n";
print F "plot 'computed-overlaps.bestC.true.eratehist'  using 2:1 title 'true'  ps 1.00 lc 2, \\\n";
print F "     'computed-overlaps.bestC.near.eratehist'  using 2:1 title 'near'  ps 1.00 lc 3, \\\n";
print F "     'computed-overlaps.bestC.false.eratehist' using 2:1 title 'false' ps 1.00 lc 1\n";
print F "\n";
print F "set title 'Computed Overlaps, non-best, error rate histogram'\n";
print F "set output 'computed-overlaps.other.eratehist.png'\n";
print F "plot 'computed-overlaps.other.true.eratehist'  using 2:1 title 'true'  ps 1.00 lc 2, \\\n";
print F "     'computed-overlaps.other.false.eratehist' using 2:1 title 'false' ps 1.00 lc 1\n";
print F "\n";
print F "\n";
print F "set title 'Computed Overlaps, best edges non-best, error-vs-length'\n";
print F "set output 'computed-overlaps.bestE.length-vs-erate.png'\n";
print F "plot 'computed-overlaps.bestE.true.dat'  using 1:2 title 'true'  ps 0.50 lc 2, \\\n";
print F "     'computed-overlaps.bestE.near.dat'  using 1:2 title 'near'  ps 0.25 lc 3, \\\n";
print F "     'computed-overlaps.bestE.false.dat' using 1:2 title 'false' ps 0.25 lc 1\n";
print F "\n";
print F "set title 'Computed Overlaps, non-best, error-vs-length'\n";
print F "set output 'computed-overlaps.bestC.length-vs-erate.png'\n";
print F "plot 'computed-overlaps.bestC.true.dat'  using 1:2 title 'true'  ps 0.50 lc 2, \\\n";
print F "     'computed-overlaps.bestC.near.dat'  using 1:2 title 'near'  ps 0.25 lc 3, \\\n";
print F "     'computed-overlaps.bestC.false.dat' using 1:2 title 'false' ps 0.25 lc 1\n";
print F "\n";
print F "set title 'Computed Overlaps, non-best, error-vs-length'\n";
print F "set output 'computed-overlaps.other.length-vs-erate.png'\n";
print F "plot 'computed-overlaps.other.true.dat'  using 1:2 title 'true'  ps 0.50 lc 2, \\\n";
print F "     'computed-overlaps.other.false.dat' using 1:2 title 'false' ps 0.25 lc 1\n";
close(F);

system("gnuplot tmp.gp && rm tmp.gp");








sub computeOverlapLength ($$$$) {
    my ($aIID, $bIID, $aHang, $bHang) = @_;

    my $aLen = $readLengths[$aIID];
    my $bLen = $readLengths[$bIID];

    my $aOvl = 0;
    my $bOvl = 0;

    #  Swiped from bogart

    if ($aHang < 0) {
      #  bHang < 0      ?     ----------  :     ----
      #                 ?  ----------     :  ----------
      #
      $aOvl = ($bHang < 0) ? ($aLen + $bHang) : ($aLen);
      $bOvl = ($bHang < 0) ? ($bLen + $aHang) : ($bLen + $aHang - $bHang);
    } else {
      #  bHang < 0      ?  ----------              :  ----------
      #                 ?     ----                 :     ----------
      #
      $aOvl = ($bHang < 0) ? ($aLen - $aHang + $bHang) : ($aLen - $aHang);
      $bOvl = ($bHang < 0) ? ($bLen)                   : ($bLen - $bHang);
    }

    return(($aOvl + $bOvl) / 2);
}
