#!/usr/bin/perl

use strict;

#
#  Reads Aaron formatted files, converts into ATAC.
#


#  Build a map from ID to IID
#
my %hmap;
my %bmap;
my $iid = 0;
open(F, "< MERYL/HUREF2.deflines") or die "Failed to open MERYL/HUREF2.deflines\n";
while (<F>) {
    chomp;
    if (m/>(\S+)\s*/) {
        #print "$1 -> $iid\n";
        $hmap{$iid} = $1;
        $iid++;
    } else {
        print "HUREF2 defline mismatch: $_\n";
    }
}
close(F);


$iid = 0;
open(F, "< MERYL/B35LC.deflines") or die "Failed to open MERYL/B35LC.deflines\n";
while (<F>) {
    chomp;
    if (m/>chr(\S+)\s*/) {
        #print "$1 -> $iid\n";
        $bmap{$iid} = $1;
        $iid++;
    } else {
        print "B35LC defline mismatch: $_\n";
    }
}
close(F);


while (<STDIN>) {
    chomp;

    if (m/^M\su\s/) {
        my @vals = split '\s+', $_;

        my $id1  = $vals[4];
        my $sta1 = $vals[5];
        my $len1 = $vals[6];
        my $end1 = $sta1 + $len1;

        my $id2  = $vals[8];
        my $sta2 = $vals[9];
        my $len2 = $vals[10];
        my $end2 = $sta2 + $len2;

        if ($vals[11] < 0) {
            my $t = $sta2;
            $sta2 = $end2;
            $end2 = $t;
        }



        #  Replace iid1 with the id.  If 1 represents a B35 object,
        #  swap the objects.
        #
        if      ($id1 =~ m/B35LC:(\d+)/) {
            #  Shucks, this should be HUREF!  Make it be.
            my $t;
            $t = $id1;   $id1  = $id2;  $id2  = $t;
            $t = $sta1;  $sta1 = $sta2; $sta2 = $t;
            $t = $end1;  $end1 = $end2; $end2 = $t;
        }

        if ($id1 =~ m/HUREF2:(\d+)/) {
            if (defined($hmap{$1})) {
                $id1 = $hmap{$1};
            } else {
                print STDERR "Didn't find $1 in hmap!\n";
                exit(1);
            }
        } else {
            print STDERR "Failed to match B35LC or HUREF2 in id1 = $id1\n";
            exit(1);
        }



        #  Replace iid2 with the id.  This is now always B35LC.
        #
        if      ($id2 =~ m/B35LC:(\d+)/) {
            if (defined($bmap{$1})) {
                $id2 = $bmap{$1};
            } else {
                print STDERR "Didn't find $1 in bmap!\n";
                exit(1);
            }
        } else {
            print STDERR "Failed to match B35LC in id2 = $id2\n";
            exit(1);
        }


        print "$id1 $sta1 $end1 : $id2 $sta2 $end2\n";

        undef $id1;
        undef $sta1;
        undef $len1;
        undef $end1;
        undef $id2;
        undef $sta2;
        undef $len2;
        undef $end2;
    }
}
