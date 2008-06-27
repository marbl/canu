#!/usr/local/bin/perl

#  Refactored from Aaron Halpern's originals (asm_parse.pl and
#  frg_onto_scf_map.pl, checked in in AS_TER, but deleted).
#
#  Expect it to take:
#    10 minutes for dros
#    90 minutes for human

use strict;

print STDERR " ******************************************************************\n";
print STDERR " **                                                              **\n";
print STDERR " **  WARNING!  This is OBSOLETE code.  Use buildPosMap instead.  **\n";
print STDERR " **                                                              **\n";
print STDERR " ******************************************************************\n";

#  SPECIAL for huref4!  Look for minAsmGap.

my $minAsmGap = 20;
my $prefix    = shift;

if (!defined($prefix)) {
    die "usage: $0 prefix < prefix.asm\n";
}

open(CTGLEN, "> $prefix.ctglen");
open(SCFLEN, "> $prefix.scflen");

open(FRGCTG, "> $prefix.frgctg");
open(UTGCTG, "> $prefix.utgctg");
open(VARCTG, "> $prefix.varctg");

open(FRGSCF, "> $prefix.frgscf");
open(UTGSCF, "> $prefix.utgscf");
open(CTGSCF, "> $prefix.ctgscf");
open(VARSCF, "> $prefix.varscf");

open(FRGSURR, "> $prefix.frgsurr");


sub readMultiLineDot {
    $_ = <STDIN>;                #  Read and discard the tag
    my $save = $/;  $/ = ".\n";  #  Prepare to read the whole thing
    my $src = <STDIN>;           #  Read it
    $/ = $save;                  #  Reset our end of line marker
    $src =~ s/\s+//g;            #  Replace spaces and newlines
    chop $src;                   #  Get rid of the trailing .
    return($src);
}


sub skipRecord {
    my $save = $/;  $/ = "}\n";  #  Prepare to read the whole thing
    my $src = <STDIN>;           #  Read it
    $/ = $save;                  #  Reset our end of line marker
}

my %contigLength;
my %surrFrags;

my $numCCO = 0;
my $numSCF = 0;
my $numDSC = 0;

while (!eof(STDIN)) {
    my $tag = <STDIN>;  chomp $tag;

    if (($tag eq "{MDI") ||
        ($tag eq "{AFG") ||
        ($tag eq "{AMP") ||
        ($tag eq "{ULK") ||
        ($tag eq "{CLK") ||
        ($tag eq "{SLK")) {
        skipRecord();
        next;
    } elsif ($tag eq "{UTG") {
        my $acc = <STDIN>;
        my $src = readMultiLineDot();
        my $cov = <STDIN>;
        my $sta = <STDIN>;
        my $abp = <STDIN>;
        my $bbp = <STDIN>;
        my $len = <STDIN>;
        my $cns = readMultiLineDot();
        my $qlt = readMultiLineDot();
        my $for = <STDIN>;
        my $nfr = <STDIN>;

        if ($acc =~ m/acc:\((\S+),(\S+)\)/) {
            $acc= $1;
        } else {
            die "Failed to find acc in $acc\n";
        }

        my $tag = <STDIN>;
        while ($tag =~ m/\{MPS/) {
            if ($sta =~ m/^sta:S/) {
                my $typ = <STDIN>;  chomp $typ;
                my $mid = <STDIN>;  chomp $mid;
                my $src = readMultiLineDot();
                my $pos = <STDIN>;  chomp $pos;
                my $dln = <STDIN>;  chomp $dln;
                my $del = <STDIN>;  chomp $del;
                my $jnk = <STDIN>;  chomp $jnk;  #  closing bracket

                #  If there are del's read them too -- we could read
                #  in $dln integers, but we could also just skip the
                #  rest of the message!  Some del's are multiple
                #  lines!
                #
                if ($jnk ne "}") {
                    skipRecord();
                }

                if ($mid =~ m/mid:(\S+)$/) {
                    $mid = $1;
                    $surrFrags{ $acc }{ $mid } = 1;
                } else {
                    die "Failed to find mid in $mid\n";
                }

                $tag = <STDIN>;  chomp $tag;

                if (($tag ne "{MPS") && ($tag ne "}")) {
                    die "Incorrect end tag $tag in MPS\n";
                }
            } else {
                #  Read and ignore MPS records
                skipRecord();
                $tag = <STDIN>;  chomp $tag;
            }
        }
        next;

    } elsif ($tag eq "{CCO") {
        $numCCO++;

        my $acc = <STDIN>;  chomp $acc;
        my $pla = <STDIN>;  chomp $pla;
        my $len = <STDIN>;  chomp $len;
        my $cns = readMultiLineDot();
        my $qlt = readMultiLineDot();
        my $for = <STDIN>;  chomp $for;
        my $npc = <STDIN>;  chomp $npc;
        my $nou = <STDIN>;  chomp $nou;
        my $nvr = <STDIN>;  chomp $nvr;

        #  Process what we have so far.

        if ($acc =~ m/acc:\((\S+),(\S+)\)/) {
            $acc= $1;
        } else {
            die "Failed to find acc in $acc\n";
        }

        my @ctgCoords;

        $ctgCoords[0]=0;
        my $I = 0;
        my $L = 0;

        foreach my $c (split //, $cns) {
            $I++;
            $L++ if ($c ne "-");
            $ctgCoords[$I] = $L;
        }

        print CTGLEN "$acc $L\n";
        $contigLength{$acc} = $L;

        #  Handle VAR's

        my $tag = <STDIN>;
        while ($tag =~ m/\{VAR/) {

            if (0) {
                my $pos = <STDIN>;  chomp $pos;
                my $nrd = <STDIN>;  chomp $nrd;
                my $nba = <STDIN>;  chomp $nba;
                my $nal = <STDIN>;  chomp $nal;
                my $rat = <STDIN>;  chomp $rat;
                my $win = <STDIN>;  chomp $win;
                my $len = <STDIN>;  chomp $len;
                my $var = readMultiLineDot();
                my $jnk = <STDIN>;  chomp $jnk;  #  closing bracket
            }

            my $pos = <STDIN>;  chomp $pos;
            my $nrd = <STDIN>;  chomp $nrd;
            my $nca = <STDIN>;  chomp $nca;
            my $anc = <STDIN>;  chomp $anc;
            my $len = <STDIN>;  chomp $len;
            my $nra = readMultiLineDot();
            my $wgt = readMultiLineDot();
            my $var = readMultiLineDot();
            my $rid = readMultiLineDot();
            my $jnk = <STDIN>;  chomp $jnk;  #  closing bracket

            if ($pos =~ m/pos:(\S+),(\S+)$/) {
                my $b = $ctgCoords[$1];
                my $e = $ctgCoords[$2];
                print VARCTG "$var $acc $b $e $nrd $nca $anc $len $nra $wgt\n";
            } else {
                die "Failed to find VAR pos in $pos\n";
            }

            $tag = <STDIN>;  chomp $tag;
        }

        #  Handle MPS's

        while ($tag =~ m/\{MPS/) {
            my $typ = <STDIN>;  chomp $typ;
            my $mid = <STDIN>;  chomp $mid;
            my $src = readMultiLineDot();
            my $pos = <STDIN>;  chomp $pos;
            my $dln = <STDIN>;  chomp $dln;
            my $del = <STDIN>;  chomp $del;
            my $jnk = <STDIN>;  chomp $jnk;  #  closing bracket

            #  If there are del's read them too -- we could read in
            #  $dln integers, but we could also just skip the rest of
            #  the message!  Some del's are multiple lines!
            #
            chomp $jnk;
            if ($jnk ne "}") {
                skipRecord();
            }

            if ($mid =~ m/mid:(\S+)$/) {
                $mid = $1;
            } else {
                die "Failed to find mid in $mid\n";
            }

            if ($pos =~ m/pos:(\d+),(\d+)$/) {
                my $b = $ctgCoords[$1];
                my $e = $ctgCoords[$2];

                # Do not use $b and $e to decide orientation -- it is
                # possible, albeit very unusual, for a fragment to fit
                # entirely in a gap in the consensus seq, in which
                # case b == e; we need $rev to be set right in this
                # case ...
                #
                my $rev = 0;
                if ($1 > $2){
                    ($b, $e) = ($e, $b);
                    $rev = 1;
                }

                print FRGCTG "$mid $acc $b $e $rev\n";
            } else {
                die "Failed to find MPS pos in $pos\n";
            }

            $tag = <STDIN>;  chomp $tag;

            if (($tag ne "{MPS") && ($tag ne "{UPS") && ($tag ne "}")) {
                die "Incorrect end tag $tag in MPS\n";
            }
        }

        #  Handle UPS's

        while ($tag =~ m/\{UPS/) {

            my $typ = <STDIN>;  chomp $typ;
            my $lid = <STDIN>;  chomp $lid;
            my $pos = <STDIN>;  chomp $pos;
            my $dln = <STDIN>;  chomp $dln;
            my $del = <STDIN>;  chomp $del;
            my $jnk = <STDIN>;  chomp $jnk;  #  closing bracket

            #  If there are del's read them too -- we could read in
            #  $dln integers, but we could also just skip the rest of
            #  the message!  Some del's are multiple lines!
            #
            if ($jnk ne "}") {
                skipRecord();
            }

            if ($lid =~ m/lid:(\S+)$/) {
                $lid = $1;
            } else {
                die "Failed to find lid in $lid\n";
            }

            if (exists $surrFrags{ $lid }) {
                foreach my $frag (keys %{$surrFrags{$lid}}) {
                    print FRGSURR "$frag $lid $acc\n";
                }
            }

            #  See similar block above
            #
            if ($pos =~ m/pos:(\d+),(\d+)$/) {
                my $b = $ctgCoords[$1];
                my $e = $ctgCoords[$2];

                my $rev = 0;
                if ($1 > $2){
                    ($b, $e) = ($e, $b);
                    $rev = 1;
                }

                print UTGCTG "$lid $acc $b $e $rev\n";
            } else {
                die "Failed to find UPS pos in $pos\n";
            }

            $tag = <STDIN>;  chomp $tag;

            if (($tag ne "{UPS") && ($tag ne "}")) {
                die "Incorrect end tag $tag in MPS\n";
            }
        }

        die "End of CCO, got '$tag'\n", if (($tag ne "") && ($tag ne "}"));
    } elsif ($tag eq "{DSC"){
        $numDSC++;

        my $acc = <STDIN>;  chomp $acc;
        my $ctg = <STDIN>;  chomp $ctg;
        my $jnk = <STDIN>;  #  the closing bracket

        #  Yes, the acc is unused -- processScaffolds uses the contig
        #  UID in the dregs file, not the scaffold UID.
        #
        if ($acc =~ m/acc:(\S+)$/) {
            $acc = $1;
        } else {
            die "Failed to find acc in DSC '$acc'\n";
        }

        if ($ctg =~ m/ctg:(\S+)$/) {
            $ctg = $1;

            print CTGSCF "$ctg $ctg 0 $contigLength{$ctg} 0\n";
            print SCFLEN "$ctg $contigLength{$ctg}\n";
        } else {
            die "Failed to find ctg in DSC '$ctg'\n";
        }
    } elsif ($tag eq "{SCF"){
        $numSCF++;

        my $acc = <STDIN>;  chomp $acc;
        my $noc = <STDIN>;  chomp $noc;

        if ($acc =~ m/acc:\((\S+),\S+\)$/) {
            $acc = $1;
        } else {
            die "Failed to find acc in SCF '$acc'\n";
        }

        if ($noc =~ m/noc:(\d+)$/) {
            $noc = $1;
        } else {
            die "Failed to find noc in SCF '$noc'\n";
        }

        #  Handle CTP's

        my $numCTPs = 0;
        my $scfLen  = 0;

        my $tag = <STDIN>;
        while ($tag =~ m/\{CTP$/) {
            my $ct1 = <STDIN>;  chomp $ct1;
            my $ct2 = <STDIN>;  chomp $ct2;
            my $mea = <STDIN>;  chomp $mea;
            my $std = <STDIN>;  chomp $std;
            my $ori = <STDIN>;  chomp $ori;
            my $jnk = <STDIN>;  #  Closing bracket

            $numCTPs++;

            if ($ct1 =~ m/ct1:(\S+)$/) {
                $ct1 = $1;
            } else {
                die "Failed to find ct1 in MPS '$ct1'\n";
            }

            if ($ct2 =~ m/ct2:(\S+)$/) {
                $ct2 = $1;
            } else {
                die "Failed to find ct2 in MPS '$ct2'\n";
            }

            if ($mea =~ m/mea:(.*)$/) {

                #  HUREF special!
                #
                #$mea = $1;
                #$mea = 50 if ($mea <= 0.0);
                #$mea = int($1);

                #  NORMAL operation
                $mea = int($1);
                $mea = $minAsmGap if ($mea < $minAsmGap);
            } else {
                die "Failed to find mea in MPS '$mea'\n";
            }

            if ($std =~ m/std:(.*)$/) {
                $std = int($1);
            } else {
                die "Failed to find std in MPS '$std'\n";
            }

            if ($ori =~ m/ori:(.*)$/) {
                $ori = $1;
            } else {
                die "Failed to find ori in MPS '$ori'\n";
            }

            my $ctgB;
            my $ctgE;
            my $ctgI;

            if ($numCTPs == 1) {
                #  If this is the first CTP, print the first contig.
                #
                $ctgB    = $scfLen;
                $scfLen += $contigLength{$ct1};
                $ctgE    = $scfLen;
                $ctgI    = $ct1;

                my $rev = (($ori eq "A") || ($ori eq "O")) ? 1 : 0;
                print CTGSCF "$ctgI $acc $ctgB $ctgE $rev\n";
            }

            if (($numCTPs > 1) || ($noc > 0)) {
                #  If it's not the first, or there's a gap, print the second contig.
                #
                $scfLen += $mea;
                $ctgB    = $scfLen;
                $scfLen += $contigLength{$ct2};
                $ctgE    = $scfLen;
                $ctgI    = $ct2;

                my $rev = (($ori eq "A") || ($ori eq "I")) ? 1 : 0;
                print CTGSCF "$ctgI $acc $ctgB $ctgE $rev\n";
            }

            $tag = <STDIN>;  chomp $tag;
        }

        #  All done with the scaffold!  Save the length
        #
        print SCFLEN "$acc $scfLen\n";
    } elsif ($tag eq "{ADT") {
        my $tag = <STDIN>;
        while ($tag =~ m/\{ADL/) {
            skipRecord();
            $tag = <STDIN>;  chomp $tag;
        }
        skipRecord();
    } else {
        print STDERR "WARNING:  Didn't handle '$tag'\n";
    }

    if (($numCCO + $numSCF + $numDSC) % 10 == 27) {
        print STDERR "CCO:$numCCO SCF:$numSCF DSC:$numDSC\r";
    }
}


close(CTGLEN);
close(SCFLEN);
close(FRGCTG);
close(UTGCTG);
close(VARCTG);
close(CTGSCF);
close(FRGSURR);

undef %surrFrags;

#
#  frag onto scaffold mapping
#

my %ctgScfID;
my %ctgScfLo;
my %ctgScfHi;
my %ctgScfRv;

open(CTGSCF, "< $prefix.ctgscf");
while(<CTGSCF>){
    my @w = split;
    $ctgScfID{$w[0]}  = $w[1];
    $ctgScfLo{$w[0]}  = $w[2];
    $ctgScfHi{$w[0]}  = $w[3];
    $ctgScfRv{$w[0]}  = $w[4];
}
close(CTGSCF);

# could do stuff to get surrogate frags by scaffold, but join works to
#
system "sort -n $prefix.ctgscf |join -1 3 -2 1 -o '1.1 1.2 1.3 2.2' $prefix.frgsurr - > $prefix.surrscf";

open(FRGCTG, "< $prefix.frgctg");
while(<FRGCTG>){
    my @w =split;

    if(defined($ctgScfID{$w[1]})){
        my $frgScfRv = 1;
        my $frgScfLo;
        my $frgScfHi;

        if ($w[4] == $ctgScfRv{$w[1]}) { # orientations match
            $frgScfRv = 0;
        }

        if ($ctgScfRv{$w[1]} == 0) { # ctg fwd in scf
            $frgScfLo = $ctgScfLo{$w[1]} + $w[2];
            $frgScfHi = $ctgScfLo{$w[1]} + $w[3];
        } else {
            $frgScfLo = $ctgScfHi{$w[1]} - $w[3];
            $frgScfHi = $ctgScfHi{$w[1]} - $w[2];
        }

        print FRGSCF "$w[0] $ctgScfID{$w[1]} $frgScfLo $frgScfHi $frgScfRv\n";
    } else {
        print FRGSCF "$w[0] $w[1] $w[2] $w[3] $w[4] CTG\n";
    }
}
close(FRGCTG);
close(FRGSCF);


open(UTGCTG, "< $prefix.utgctg");
while(<UTGCTG>){
    my @w =split;

    if(defined($ctgScfID{$w[1]})){
        my $utgScfRv = 1;
        my $utgScfLo;
        my $utgScfHi;

        if ($w[4] == $ctgScfRv{$w[1]}) { # orientations match
            $utgScfRv = 0;
        }

        if ($ctgScfRv{$w[1]} == 0) { # ctg fwd in scf
            $utgScfLo = $ctgScfLo{$w[1]} + $w[2];
            $utgScfHi = $ctgScfLo{$w[1]} + $w[3];
        } else {
            $utgScfLo = $ctgScfHi{$w[1]} - $w[3];
            $utgScfHi = $ctgScfHi{$w[1]} - $w[2];
        }

        print UTGSCF "$w[0] $ctgScfID{$w[1]} $utgScfLo $utgScfHi $utgScfRv\n";
    } else {
        print UTGSCF "$w[0] $w[1] $w[2] $w[3] $w[4] CTG\n";
    }
}
close(UTGCTG);
close(UTGSCF);


open(VARCTG, "< $prefix.varctg");
while(<VARCTG>){
    my @w =split;

    if(defined($ctgScfID{$w[1]})){
        my $varScfLo;
        my $varScfHi;

        if ($ctgScfRv{$w[1]} == 0) { # ctg fwd in scf
            $varScfLo = $ctgScfLo{$w[1]} + $w[2];
            $varScfHi = $ctgScfLo{$w[1]} + $w[3];
        } else {
            $varScfLo = $ctgScfHi{$w[1]} - $w[3];
            $varScfHi = $ctgScfHi{$w[1]} - $w[2];
        }

        print VARSCF "$w[0] $ctgScfID{$w[1]} $varScfLo $varScfHi $w[4] $w[5] $w[6] $w[7] $w[8] $w[9]\n";
    } else {
        print VARSCF "$w[0] $w[1] $w[2] $w[3] $w[4] $w[5] $w[6] $w[7] $w[8] $w[9] CTG\n";
    }
}
close(VARCTG);
close(VARSCF);
