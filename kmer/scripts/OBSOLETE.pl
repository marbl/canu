

######################################################################
#
#  Subroutine to read a single sim4 polish, and return the coverage,
#  percent identity and text of the polish.
#
sub readPolishForScore (*) {
    local *READPOLISHFH = shift;
    my $line;
    my $save;

    #  Skip lines until the next match.  If used properly, on a proper
    #  file, this should be silent.  After the loop, we are at the
    #  start of a polish; the line should be "sim4begin".
    #
    $line = <READPOLISHFH>;
    while (defined($line) && ($line !~ m/^sim4begin$/)) {
        chomp $line;
        print STDERR "Skipped: '$line'\n";
        $line = <READPOLISHFH>;
    }
    $save = $line;

    #  Return now if were are out of file
    #
    return(undef, undef, undef) if (eof(READPOLISHFH));

    #  Read the description line
    #
    $line  = <READPOLISHFH>;
    $save .= $line;

    my $estLen;
    my $pA;
    my $pT;
    my $percentID;
    my $matchOrientation;
    my $strand;

    if ($line =~ m/^\d+\[(\d+)-+(\d+)-+(\d+)\]\s+\d+\[\d+-\d+\]\s+\<\d+-\d+-(\d+)-(\w+)-(\w+)\>$/) {
        $estLen           = $1;
        $pA               = $2;
        $pT               = $3;
        $percentID        = $4;
        $matchOrientation = $5;
        $strand           = $6;
    } else {
        print STDERR "expecting description line, got: '$line'\n";
        return(undef, undef, undef);
    }

    my $coverage = 0;
    my $estStart = $estLen;
    my $estEnd   = 0;

    $line  = <READPOLISHFH>;
    $save .= $line;

    while (defined($line) && ($line !~ m/^sim4end$/)) {
        if ($line =~ /^(\d+)-(\d+)/) {
            $coverage += $2 - $1 + 1;
            $estStart  = $1 if ($estStart > $1);
            $estEnd    = $2 if ($estEnd < $2);
        }

        $line  = <READPOLISHFH>;
        $save .= $line;
    }

    if ($matchOrientation eq "complement") {
        ($estEnd, $estStart) = ($estLen - $estStart + 1, $estLen - $estEnd + 1);
    }

    $estLen -= $pA + $pT;
    if ($estLen > 0) {
        $coverage = 100.0 * $coverage / $estLen;
    } else {
        $coverage = 0;
    }

    return($coverage, $percentID, $strand, $save);
}


#
#  Returns an alternate (and hopefully, more generic) absolute path
#  for a given input path.  Input path can be relative or absolute.
#
#  If the path doesn't exist, it returns undef.
#
#  XXX: This doesn't always work:
#
#    /home/walenzbp/home/walenzbp resolves to /home/walenzbp
#    /work/assembly is mounted as /devel/I$10/assembly
#

sub resolvePath {
    my $IN = shift @_;
    my $I = $IN;

    #  Fail if the path doesn't exist
    #
    return undef if (! -e $IN);

    #  Save the inode of the path.  We'll use this later to make sure
    #  that we have found the exact same path.
    #
    my $inode = (stat($IN))[1];

    #print STDERR "Got $inode for $IN\n";

    #  If the path is relative, make it absolute.
    #
    $I = cwd() . "/" . $I if ($I !~ m!^/!);

    #  Squash multiple "/"s into a single one
    #
    $I =~ tr!/!/!s;

    #  Replace "/./" with "/" -- the loop is needed because it doesn't
    #  seem to use replaced things in the next patten.  "/././" is
    #  replaced with "/./" (by replacing the first "/./" with "/", but
    #  then the replace examines only the last "./", which doesn't
    #  match.
    #
    while ($I =~ s!/\./!/!g) {};

    #  Remove any trailing "/"
    #
    $I =~ s!/$!!;

    #  Remove "./" from the start, and "/." from the end.
    #
    $I =~ s!^(\./)*!!;
    $I =~ s!(/\.)*$!!;

    #  Traverse the path, building a new one, but removing ".." and
    #  the previous component.
    #
    my @P;
    foreach my $p (split '/', $I) {
        if ($p eq "..") {
            pop @P;
        } else {
            push @P, $p;
        }
    }

    #  The first thing on @P is empty (or, should be empty).
    #
    shift @P if ($P[0] eq "");

    #  Now, while @P has stuff, find the shortest path that exists.
    #  If you gave me a path that doesn't exist, I'll return an empty
    #  path.
    #
    my $r = "";
    my $p = "";
    while (scalar(@P) > 0) {
        $p = "/" . join '/', @P;
        $r = $p if (-e $p);
        shift @P;
    }

    if (($r ne "") && (-e $r)) {
        return($r);
    } else {
        if (-e $IN) {
            print STDERR "WARNING: Couldn't resolve path '$IN', but it exists!\n";
        } else {
            print STDERR "WARNING: Couldn't resolve path '$IN', and it doesn't exist!\n";
        }
        return($IN);
    }
}
