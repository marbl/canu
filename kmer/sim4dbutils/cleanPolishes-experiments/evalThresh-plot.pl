open(F, "< evalThresh.pl.out");

while (!eof(F)) {
    $_ = <F>;
    if (m/at least (\d+)bp/) {
        print "$1\t";
        $_ = <F>;
        $_ = <F>;if (m/oneExon:\s+(\d+)/) { print "$1\t"; } else { print STDERR "no 1\n"; }
        $_ = <F>;if (m/allSmall.*:\s+(\d+)/) { print "$1\t"; } else { print STDERR "no 2\n"; }
        $_ = <F>;if (m/good:\s+(\d+)/) { print "$1\t"; } else { print STDERR "no 3\n"; }
        $_ = <F>;if (m/probably\sgood:\s+(\d+)/) { print "$1\t"; } else { print STDERR "no 4\n"; }
        $_ = <F>;if (m/junkExonsLeft:\s+(\d+)/) { print "$1\t"; } else { print STDERR "no 5\n"; }
        $_ = <F>;if (m/junkExonsRight:\s+(\d+)/) { print "$1\t"; } else { print STDERR "no 6\n"; }
        $_ = <F>;if (m/junkExonsBoth:\s+(\d+)/) { print "$1\t"; } else { print STDERR "no 7\n"; }
        $_ = <F>;if (m/intronOnGap:\s+(\d+)/) { print "$1\t"; } else { print STDERR "no 8\n"; }
        $_ = <F>;if (m/total:\s+(\d+)/) { print "$1"; } else { print STDERR "no 9\n"; }
        print "\n";
    }
}

close(F);




