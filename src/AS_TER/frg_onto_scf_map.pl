#!/usr/local/bin/perl

$ctgOnScf = shift;
$frgOnCtg = shift;

open(CTGONSCF,"$ctgOnScf");
while(<CTGONSCF>){
    @w = split;
    $ctgScf{$w[0]}=$w[1];
    $ctgScfLo{$w[0]} = $w[2];
    $ctgScfHi{$w[0]} = $w[3];
    $ctgScfRev{$w[0]}=$w[4];
}
close(CTGONSCF);

open(FRGONCTG,"$frgOnCtg");
while(<FRGONCTG>){

    @w =split;

#    $frgCtg{$w[0]}=$w[1];
#    $frgCtgLo{$w[0]} = $w[2];
#    $frgCtgHi{$w[0]} = $w[3];
#    $frgCtgRev{$w[0]}=$w[4];

    if(defined($ctgScf{$w[1]})){

	if($w[4] == $ctgScfRev{$w[1]}){ # orientations match
	    $frgScfRev = 0;
	} else {
	    $frgScfRev = 1;
	}
	if($ctgScfRev{$w[1]} == 0) { # ctg fwd in scf
	    $frgScfLo = $ctgScfLo{$w[1]} + $w[2];
	    $frgScfHi = $ctgScfLo{$w[1]} + $w[3];
	} else {
	    $frgScfLo = $ctgScfHi{$w[1]} - $w[3];
	    $frgScfHi = $ctgScfHi{$w[1]} - $w[2];
	}

	print "$w[0] $ctgScf{$w[1]} $frgScfLo $frgScfHi $frgScfRev\n";

    } else {
	
	print "$w[0] $w[1] $w[2] $w[3] $w[4] CTG\n";

    }

}
close(FRGONCTG);

exit(0);
