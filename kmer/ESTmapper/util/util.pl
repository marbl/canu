

######################################################################
#
#  Copies ESTs from $input to $outputYes or $outputNay, depending on
#  if the EST is present in the sim4db output file $polishes
#
sub splitFastABasedOnPolishes {
    my ($input, $polishes, $outputYea, $outputNay) = @_;
    my %est;

    #  Read the polishes from $polishes, and save the EST deflines
    #
    open (F, "< $polishes");
    while (<F>) {
        if (m/^edef=>(.*)$/) {
            my $h = $1;
            $h =~ s/^\s+//;
            $h =~ s/\s+$//;
            $est{$h} = 1;
        }
    }
    close(F);

    #  Copy the input to the two output files
    #
    open(F, "< $input");
    open(Y, "> $outputYea");
    open(N, "> $outputNay");

    while (!eof(F)) {
        my ($h, $s) = fasta::nextFastA(*F);

        if (defined($est{$h})) {
            print Y ">$h\n$s\n";
        } else {
            print N ">$h\n$s\n";
        }
    }
    close(N);
    close(Y);
    close(F);
}



######################################################################
#
#  Copies ESTs from $input to $output, given that the corresponding
#  value in the $counts file is ZERO.
#
sub copyZeroFastA {
    my ($input, $counts, $output) = @_;

    open(F, "< $counts");
    open(I, "< $input");
    open(O, "> $output");
    while (!eof(I)) {
        my ($h, $s) = fasta::nextFastA(*I);
        print O ">$h\n$s\n" if (0 == int <F>);
    }
    close(O);
    close(I);
    close(F);
}


######################################################################
#
#  Copies ESTs from $input to $output, given that their defline doesn't
#  occur in $discard.  All files are multi-fasta.
#
sub subtractFastAfromFastA {
    my ($input, $discard1, $discard2, $output) = @_;
    my %deflines;

    #  Read the deflines in $discard[12], and save them.
    #
    if (-e $discard1) {
        open(F, "< $discard1");
        while (!eof(F)) {
            my ($h, $s) = fasta::nextFastA(*F);
            $deflines{$h} = 1;
        }
        close(F);
    }

    if (-e $discard2) {
        open(F, "< $discard2");
        while (!eof(F)) {
            my ($h, $s) = fasta::nextFastA(*F);
            $deflines{$h} = 1;
        }
        close(F);
    }

    #  Now, copy $input to $output
    #
    open(I, "< $input");
    open(O, "> $output");
    while (!eof(I)) {
        my ($h, $s) = fasta::nextFastA(*I);
        print O ">$h\n$s\n" if (!defined($deflines{$h}));
    }
    close(O);
    close(I);
}



######################################################################
#
#  Generates a report on a set of polishes.
#
#  number of cDNA-scaffold matches
#  number of different cDNA sequences in the set
#  number of different scaffolds in the set
#
sub summarizePolishes {
    my (@files) = @_;

    my %est;
    my %scf;
    my $mat = 0;
    my $ests = 0;
    my $scfs = 0;

    foreach my $infile (@files) {
        open(INPUT, "< $infile");

        while (<INPUT>) {
            if (m/^sim4begin$/) {
                $mat++;
            } elsif (m/^edef=/) {
                $ests++;
                $est{$_} = 1;
            } elsif (m/^ddef=/) {
                $scfs++;
                $scf{$_} = 1;
            }
        }

        close(INPUT);
    }

    if (($ests != $mat) || ($scfs != $mat)) {
        print STDERR "WARNING: summarizePolishes counted\n";
        print STDERR "           $mat matches\n";
        print STDERR "           $ests cDNA deflines\n";
        print STDERR "           $scfs scaffold deflines\n";
        print STDERR "         The number of deflines and the number of matches should be the same!\n";
    }

    return($mat, (scalar (keys %est)), (scalar (keys %scf)));
}


1;
