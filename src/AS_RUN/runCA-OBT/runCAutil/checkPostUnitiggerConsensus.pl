use strict;

#  Check consensus output, prod the user into fixing the broken ones.

sub checkPostUnitiggerConsensus {

    if (! -e "$wrk/5-consensus/$asm.cgi") {
        my $failedJobs = 0;
        my @CGBfiles;

        #
        #  Check that consensus finished properly
        #
        
        open(CGB, "ls $wrk/5-consensus/*.cgb |") or die;
        while (<CGB>) {
            chomp;

            if (m/^.*(\d\d\d).cgb$/) {
                push @CGBfiles, $1;

                if ((! -e "$wrk/5-consensus/$1.success") ||
                    (! -e "$wrk/5-consensus/${asm}_$1.cgi")) {
                    print STDERR "$wrk/5-consensus/$1 failed -- no .success or no .cgi!\n";
                    $failedJobs++;
                }
            } else {
                print STDERR "WARNING: didn't match CGB '$_'!\n";
            }
        }
        close(CGB);

        if ($failedJobs) {
            print STDERR "$failedJobs failed.  Good luck.\n";
            exit(1);
        }

        #
        #  Consolidate all the output
        #

        foreach my $fid (@CGBfiles) {
            if (runCommand("cat $wrk/5-consensus/${asm}_$fid.cgi >> $wrk/5-consensus/$asm.cgi")) {
                print STDERR "Failed.\n";
                rename "$wrk/5-consensus/$asm.cgi", "$wrk/5-consensus/$asm.cgi.FAILED";
                exit(1);
            }
        }
        close(CGB);
    }
}

1;
