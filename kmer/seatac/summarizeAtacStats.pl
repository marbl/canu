#!/usr/local/bin/perl




#  generates summary of the meryl min phase

printf("\n%-30s %5d %5d %8d\n",
       "meryl min", "user", "sys", "maxRSS");

open(F, "ls min*stats |");
while (<F>) {
    chomp $_;
    my $file = $_;

    my $ut = 0;
    my $st = 0;
    my $mr = 0;
    my $bt = 0;

    my $tm = 0;
    my $dm = 0;
    my $um = 0;
    my $b  = 0;
    my $p  = 0;

    open(G, "< $file");
    while (<G>) {
        if (m/userTime:\s+(\d+)/) {
            $ut = $1;
        }
        if (m/systemTime:\s+(\d+)/) {
            $st = $1;
        }
        if (m/maxrss:\s+(\d+)/) {
            $mr = $1;
        }
    }
    close(G);

    printf("%-30s %5d %5d %8d\n",
           $file, $ut, $st, $mr);
}
close(F);






#  generates summary of the build phase

printf("\n%-30s %5s %5s %5s %8s %10s %10s %10s %8s   %8s\n",
       "seatac build", "user", "sys", "wall", "maxRSS", "totMer", "distinctMer", "uniqueMer", "bktSize", "posnSize");

open(F, "ls *build*out |");
while (<F>) {
    chomp $_;
    my $file = $_;

    my $ut = 0;
    my $st = 0;
    my $mr = 0;
    my $bt = 0;

    my $tm = 0;
    my $dm = 0;
    my $um = 0;
    my $b  = 0;
    my $p  = 0;

    open(G, "< $file");
    while (<G>) {
        if (m/userTime:\s+(\d+)/) {
            $ut = $1;
        }
        if (m/systemTime:\s+(\d+)/) {
            $st = $1;
        }
        if (m/maxrss:\s+(\d+)/) {
            $mr = $1;
        }
        if (m/build:\s+(\d+)/) {
            $bt = $1;
        }


        if (m/Found\s+(\d+)\s+total/) {
            $tm = $1;
        }
        if (m/Found\s+(\d+)\s+distinct/) {
            $dm = $1;
        }
        if (m/Found\s+(\d+)\s+unique/) {
            $um = $1;
        }
        if (m/Allocated\s+(\d+)\s*KB\s+for\s+buckets/) {
            $b = $1;
        }
        if (m/Allocated\s+(\d+)\s*KB\s+for\s+positions/) {
            $p = $1;
        }
    }
    close(G);

    printf("%-30s %5d %5d %5d %8d %10d %10d %10d %8dKB %8dKB\n",
           $file, $ut, $st, $bt, $mr, $tm, $dm, $um, $b, $p);
}
close(F);

#  generates summary of the search phase

printf("\n%-50s %5s %5s %5s %5s %5s %9s %9s %10s\n",
       "seatac search", "user", "sys", "build", "srch", "total", "usr/srch", "usr/totl", "maxRSS");

open(F, "ls *segment*stats |");
while (<F>) {
    chomp $_;
    my $file = $_;

    my $ut = 0;
    my $st = 0;
    my $mr = 0;

    my $btt = 0;
    my $stt = 0;
    my $ttt = 0;

    open(G, "< $file");
    while (<G>) {
        if (m/userTime:\s+(\d+)/) {
            $ut = $1;
        }
        if (m/systemTime:\s+(\d+)/) {
            $st = $1;
        }
        if (m/maxrss:\s+(\d+)/) {
            $mr = $1;
        }
        if (m/build:\s+(\d+)/) {
            $btt = $1;
        }
        if (m/search:\s+(\d+)/) {
            $stt = $1;
        }
        if (m/total:\s+(\d+)/) {
            $ttt = $1;
        }
    }
    close(G);

    printf("%-50s %5d %5d %5d %5d %5d %9.6f %9.6f %10d\n",
           $file, $ut, $st, $btt, $stt, $ttt, $ut / $stt, $ut / $ttt, $mr);
}
close(F);
