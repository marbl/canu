use strict;

sub runLSF {
    my $cmd = shift @_;

    open(F, "$cmd |") or die "LSF submission failed:\n  $cmd\n";
    my $res = <F>;
    close(F);

    if ($res =~ m/Job\s+<(\d+)>\s/) {
        $res = $1;
    } else {
        die "ERROR:  LSF possibly failed to start; job number in result unknown: '$res'\n";
    }

    return($res);
}



1;
