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
 #    src/AS_BAT/erate-estimate-test-based-on-mapping.pl
 #
 #  Modifications by:
 #
 #    Brian P. Walenz from 2014-NOV-15 to 2015-AUG-07
 #      are Copyright 2014-2015 Battelle National Biodefense Institute, and
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

my $erateFile  = "BO/test.blasr.sam.coords.longestErate";
my $uidMapFile = "CAordered/test.gkpStore.fastqUIDmap";
my $olapStore  = "CAordered/test.ovlStore";

my %erate5;  #  Either the mapped global erate, or the lowest erate for each end
my %erate3;


#  Load the map from IID to UID and name.

my %UIDtoIID;
my %UIDtoNAM;

my %IIDtoUID;
my %IIDtoNAM;

my %NAMtoUID;
my %NAMtoIID;


if (defined($uidMapFile)) {
    print STDERR "Read names from '$uidMapFile'\n";

    open(F, "< $uidMapFile");
    while (<F>) {
        my @v = split '\s+', $_;

        #  UID IID NAM  UID IID NAM

        $UIDtoIID{$v[0]} = $v[1];
        $UIDtoNAM{$v[0]} = $v[2];

        $IIDtoUID{$v[1]} = $v[0];
        $IIDtoNAM{$v[1]} = $v[2];

        $NAMtoUID{$v[2]} = $v[0];
        $NAMtoIID{$v[2]} = $v[1];

        if (scalar(@v) == 6) {
            $UIDtoIID{$v[3]} = $v[4];
            $UIDtoNAM{$v[3]} = $v[5];

            $IIDtoUID{$v[4]} = $v[3];
            $IIDtoNAM{$v[4]} = $v[5];

            $NAMtoUID{$v[5]} = $v[3];
            $NAMtoIID{$v[5]} = $v[4];
        }
    }
    close(F);
}



if (0) {
print STDERR "Load global read erates from '$erateFile'\n";

my $found = 0;
my $lost  = 0;

if (! -e "$erateFile") {
    die "Run erate-coords.pl in BO directory.\n";
}

open(F, "< $erateFile") or die;
while(<F>) {
    my @v = split '\s+', $_;

    #  This should be read NAME, not assmelby UID.
    my $iid = $NAMtoIID{$v[0]};

    if (defined($iid)) {
        $found++;
        $erate5{$iid} = 100 - $v[1];
        $erate3{$iid} = 100 - $v[1];
    } else {
        $lost++;
        #print STDERR "Didn't find name '$v[0]' - dropped from read set?\n";
    }
}
close(F);

print STDERR "Found $found erates, lost $lost.\n";
}




if (1) {
print STDERR "Load overlap erates from '$olapStore'\n";

open(F, "overlapStore -d $olapStore |");
while (<F>) {
    s/^\s+//;
    s/\s+$//;

    #print "$_\n";

    my ($aIID, $bIID, $orient, $aHang, $bHang, $erate, $crate) = split '\s+', $_;

    if (($aHang <= 0) && ($bHang >= 0)) {
        #  A contained in B
        next;
    }

    elsif (($aHang >= 0) && ($bHang <= 0)) {
        #  B contained in A
        next;
    }

    elsif (($aHang <= 0) && ($bHang <= 0)) {
        #  Dovetail off the 5' end
        $erate5{$aIID} = $erate   if ((!exists($erate5{$aIID})) || ($erate < $erate5{$aIID}));
    }

    elsif (($aHang >= 0) && ($bHang >= 0)) {
        #  Dovetail off the 3' end
        $erate3{$aIID} = $erate   if ((!exists($erate3{$aIID})) || ($erate < $erate3{$aIID}));
    }

    else {
        die;
    }
}
}





#  Scan overlaps, try to decide true from false.

print STDERR "Read overlas from '$olapStore'\n";

my @ratioTrue;
my @ratioFalse;

open(F, "overlapStore -d $olapStore |");
while (<F>) {
    s/^\s+//;
    s/\s+$//;

    #print "$_\n";

    my ($aIID, $bIID, $orient, $aHang, $bHang, $erate, $crate) = split '\s+', $_;

    my $aErate = 0;
    my $bErate = 0;

    if (($aHang <= 0) && ($bHang >= 0)) {
        #  A contained in B
        next;
    }

    elsif (($aHang >= 0) && ($bHang <= 0)) {
        #  B contained in A
        next;
    }


    elsif (($aHang <= 0) && ($bHang <= 0) && ($orient eq "N")) {
        #  Dovetail off the 5' end of A to the 3' end of B
        $aErate = $erate5{$aIID};
        $bErate = $erate3{$bIID};
    }

    elsif (($aHang <= 0) && ($bHang <= 0) && ($orient eq "I")) {
        #  Dovetail off the 5' end of A to the 5' end of B
        $aErate = $erate5{$aIID};
        $bErate = $erate5{$bIID};
    }

    elsif (($aHang >= 0) && ($bHang >= 0) && ($orient eq "N")) {
        #  Dovetail off the 3' end of A to the 5' end of B
        $aErate = $erate3{$aIID};
        $bErate = $erate5{$bIID};
    }

    elsif (($aHang >= 0) && ($bHang >= 0) && ($orient eq "I")) {
        #  Dovetail off the 3' end of A to the 3' end of B
        $aErate = $erate3{$aIID};
        $bErate = $erate3{$bIID};
    }

    else {
        die;
    }


    print STDERR "Didn't find erate for IID $aIID.\n  $_\n"  if (!defined($aErate));
    print STDERR "Didn't find erate for IID $bIID.\n  $_\n"  if (!defined($bErate));

    next  if (!defined($aErate));
    next  if (!defined($bErate));

    my $ratio = int(100 * ($erate) / ($aErate + $bErate));

    my $dist = ($aIID < $bIID) ? ($bIID - $aIID) : ($aIID - $bIID);

    die if ($dist < 0);

    if ($dist < 50) {
        $ratioTrue[$ratio]++;
    } else {
        $ratioFalse[$ratio]++;
    }
}
close(F);

open(F, "> r");
for (my $x=0; $x<1000; $x++) {
    next  if (!defined($ratioTrue[$x]) && (!defined($ratioFalse[$x])));

    $ratioTrue[$x] += 0;
    $ratioFalse[$x] += 0;

    my $v = $x / 100;

    print F "$v\t$ratioTrue[$x]\t$ratioFalse[$x]\n";
}
close(F);
