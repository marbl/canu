

my %LSFvars;
my %LSFjobs;

sub initializeLSF {
    $LSFvars{'queue'} = "";
    $LSFvars{'project'} = "";
    $LSFvars{'select'} = 500;
    $LSFvars{'rusage'} = 500;
    $LSFvars{'output'} = "/dev/null";
    $LSFvars{'error'} = undef;
    $LSFvars{'email'} = undef;
    $LSFvars{'options'} = undef;
}

sub setLSFQueue {
    $LSFvars{'queue'} = shift;
}

sub setLSFProject {
    $LSFvars{'project'} = shift;
}

sub setLSFMemory {
    $LSFvars{'select'} = shift;
    $LSFvars{'rusage'} = shift;
}

sub setLSFOutput {
    $LSFvars{'output'} = shift;
    undef $LSFvars{'email'};
}

sub setLSFError {
    $LSFvars{'error'} = shift;
}

sub setLSFEmail {
    $LSFvars{'email'} = shift;
    undef $LSFvars{'output'};
}

sub setLSFOptions {
    $LSFvars{'options'} = shift;
}

sub submitToLSF {
    my @jobs = @_;

    while (scalar(@JOBS) > 0) {
        my $cmd = "";
        my $res;

        $cmd .= "bsub -q $LSF{'queue'} -P $LSF{'project'} ";
        $cmd .= "-R \"select[physmem>$LSF{'select'}]rusage[physmem=$LSF{'rusage'}]\" ";
        $cmd .= "-e $LSF{'email'} "  if (defined($LSF{'email'}));
        $cmd .= "-o $LSF{'output'} " if (defined($LSF{'output'}));
        $cmd .= "-u $LSF{'error'} "  if (defined($LSF{'error'}));
        $cmd .= shift @jobs;

        $res = `$cmd`;
        chomp $res;

        if ($res =~ m/^Job <(\d+)> is submitted to queue <\S+>.$/) {
            $LSFjobs{$1} = $res;
        } else {
            print STDERR "Job might not be submitted.  Result is '$res'.\n";
        }
    }
}

sub waitForLSF {

    while (checkLSF()) {
        sleep 60;
    }
}

sub checkLSF {
    my %LSFruns;
    my $notDone = 0;

    open(F, "bjobs |");
    while (<F>) {
        chomp;

        if (m/^(\d+)\s/) {
            $LSFruns{int($1)} = $_;
        }
    }
    close(F);

    foreach my $k (keys %LSFjobs) {
        if (defined($LSFruns{$k})) {
            $notDone++;
        }
    }

    return($notDone);
}
