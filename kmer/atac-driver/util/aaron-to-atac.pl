#!/usr/bin/perl

use strict;

#
#  Reads Aaron formatted files, converts into ATAC.
#


#  XXX DOES NOT CURRENTLY MAP IDs -> IIDs


#  Fudge up an atac header.
#

print "! format atac 1.0\n";
print "/assemblyFile1=MERYL/B35LC.fasta\n";
print "/assemblyFile2=MERYL/HUREF2.fasta\n";
print "/assemblyId1=B35LC\n";
print "/assemblyId2=HUREF2\n";


#  Build a map from ID to IID
#
my %hmap;
my %bmap;
my $iid = 0;
open(F, "< GENOMES/HUREF2.deflines") or die "Failed to open GENOMES/HUREF2.deflines\n";
while (<F>) {
    chomp;
    if (m/>(\S+)\s*/) {
        #print "$1 -> $iid\n";
        $hmap{$1} = $iid;
        $iid++;
    } else {
        print "HUREF2 defline mismatch: $_\n";
    }
}
close(F);


$iid = 0;
open(F, "< GENOMES/B35LC.deflines") or die "Failed to open GENOMES/B35LC.deflines\n";
while (<F>) {
    chomp;
    if (m/>chr(\S+)\s*/) {
        #print "$1 -> $iid\n";
        $bmap{$1} = $iid;
        $iid++;
    } else {
        print "B35LC defline mismatch: $_\n";
    }
}
close(F);


#  Count the number of times we see an ID that we don't know
my %skipped1;
my %skipped2;
my $skipped    = 0;
my $skippedLen = 0;

my %used1;
my %used2;
my $used    = 0;
my $usedLen = 0;

my $matchId = 0;

while (<STDIN>) {
    chomp;

    my @vals = split '\s+', $_;

    if (scalar(@vals) != 7) {
        #print STDERR $_;
        next;
    }

    # id1 sta1 end1 : id2 sta2 end2

    my $id1  = $vals[0];
    my $sta1 = $vals[1];
    my $end1 = $vals[2];
    my $len1 = $end1 - $sta1;

    if ($sta1 > $end1) {
        print STDERR "Flipped first sequence?  $_\n";
        next;
    }

    if ($sta1 < 0) {
        print STDERR "Negative?  $_\n";
        next;
    }

    my $id2  = $vals[4];
    my $sta2 = $vals[5];
    my $end2 = $vals[6];
    my $len2 = $end2 - $sta2;

    my $ori  = 1;

    if ($sta2 > $end2) {
        $ori  = -1;
        $sta2 = $vals[6];
        $end2 = $vals[5];
        $len2 = $end2 - $sta2;
    }

    if ($sta2 < 0) {
        print STDERR "Negative?  $_\n";
        next;
    }

    #  first . is match index
    #  second . is parent

    if      (!defined($hmap{$id1})) {
        $skipped1{$id1}++;
        $skipped++;
        $skippedLen += $len1;
    } elsif (!defined($bmap{$id2})) {
        $skipped2{$id2}++;
        $skipped++;
        $skippedLen += $len2;
    } else {
        $used1{$id1}++;
        $used2{$id2}++;
        $used++;
        $usedLen += $len1;

        print "M u $matchId . B35LC:$bmap{$id2} $sta2 $len2 1 HUREF2:$hmap{$id1} $sta1 $len1 $ori\n";

        if ($len1 < 0) {
            print STDERR "$_\n";
            print STDERR "M u $matchId . B35LC:$bmap{$id2} $sta2 $len2 1 HUREF2:$hmap{$id1} $sta1 $len1 $ori\n";
        }

        if ($len2 < 0) {
            print STDERR "$_\n";
            print STDERR "M u $matchId . B35LC:$bmap{$id2} $sta2 $len2 1 HUREF2:$hmap{$id1} $sta1 $len1 $ori\n";
        }

        $matchId++;
    }
}

print STDERR "Skipped ", scalar(keys %skipped1), " sequences in HUREF2.\n";
print STDERR "Skipped ", scalar(keys %skipped2), " sequences in B35LC.\n";

foreach my $k (sort keys %used2) {
    print STDERR "  $k (map=$bmap{$k})-- $used2{$k} matches used\n";
}
foreach my $k (sort keys %skipped2) {
    print STDERR "  $k (map=$bmap{$k}) -- $skipped2{$k} matches skipped\n";
}

print STDERR "Skipped a total of $skipped matches with length $skippedLen.\n";
print STDERR "\n";
print STDERR "Used ", scalar(keys %used1), " sequences in HUREF2.\n";
print STDERR "Used ", scalar(keys %used2), " sequences in B35LC.\n";
print STDERR "Used a total of $used matches with length $usedLen.\n";
