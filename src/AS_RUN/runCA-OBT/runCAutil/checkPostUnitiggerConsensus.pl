use strict;

#  Check consensus output, prod the user into fixing the broken ones.

sub checkPostUnitiggerConsensus (@) {
    my @cgbFiles = @_;

    return if (-e "$wrk/5-consensus/consensus.success");

    if (! -e "$wrk/5-consensus/$asm.cgi") {
        my $failedJobs = 0;
        my @cgbIndices;

        #
        #  Check that consensus finished properly
        #

        foreach my $f (@cgbFiles) {
            if ($f =~ m/^.*(\d\d\d).cgb$/) {
                push @cgbIndices, $1;
            } else {
                die "Didn't match '$f' for CGB filename!\n";
            }
        }

        foreach my $f (@cgbIndices) {
            if ((! -e "$wrk/5-consensus/${asm}_$f.success") ||
                (! -e "$wrk/5-consensus/${asm}_$f.cgi")) {
                print STDERR "$wrk/5-consensus/$f failed -- no .success or no .cgi!\n";
                $failedJobs++;
            }
        }

        die  "$failedJobs failed.  Good luck.\n" if ($failedJobs);

        #
        #  Consolidate all the output
        #

        foreach my $fid (@cgbIndices) {
            if (runCommand("cat $wrk/5-consensus/${asm}_$fid.cgi >> $wrk/5-consensus/$asm.cgi")) {
                rename "$wrk/5-consensus/$asm.cgi", "$wrk/5-consensus/$asm.cgi.FAILED";
                die "cat failed?\n";
            }
        }
    }

    touch ("$wrk/5-consensus/consensus.success");
}

1;
