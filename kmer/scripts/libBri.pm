#!/usr/local/bin/perl

#                    Confidential -- Do Not Distribute
#   Copyright (c) 2002 PE Corporation (NY) through the Celera Genomics Group
#                           All Rights Reserved.

package libBri;

use strict;
use POSIX "sys_wait_h";
use Cwd;

$| = 1;

#  Called by "use libBri;"
sub import () {
}


#  Sun Jul 14 17:41:05 EDT 2002
#  Modified nextFastA() and splitFastABasedOnPolishes() to strip whitespace
#  from the start and end of each defline.  Not explicitly tested.

#  Fri Aug  2 15:40:31 EDT 2002
#  Added resolvePath() to strip off the temporary mount points from a path





######################################################################
#
#  Function to stream a multi-fasta file.  Usage:
#
#  my $header;
#  my $sequence;
#
#  open(F, "< $ARGV[0]");
#  while (!eof(F)) {
#    ($header, $sequence) = nextFastA(*F);
#    print ">$header\n";
#    print "$sequence\n";
#  }
#  close(F);
#
sub nextFastA {
    my ($fh) = @_;
    my $h;
    my $s;

    if (!eof($fh)) {
        my $sep = $/;        #  Save the file line separator

        $/ = "\n";
        $h = <$fh>;

        $h =~ s/^\s+//;
        $h =~ s/\s+$//;

        #chomp $h;

        $/ = ">";
        $s = <$fh>;
        $s =~ s/\s//gm;

        #  Fix things -- the '>' is at the start of the header only if this
        #  is the first sequence read from the file.  Likewise, the sequence
        #  will have the '>' at the end if it is not the last sequence read.
        #
        $h =~ s/^>//g;
        $s =~ s/>$//g;

        $/ = $sep;        #  Restore the file line separator
    }

    return($h, $s);
}


######################################################################
#
#  Returns a modified 'raw' string, using the current values for the
#  info line.  DOES NOT rewrite the exons.
#
sub updatePolishInfoLine {
    my %p = @_;
    my @L = split '\n', $p{'raw'};
    my $l;

    shift @L;
    shift @L;

    $l  = "sim4begin\n";
    $l .= "$p{'estID'}\[$p{'estLen'}-$p{'pA'}-$p{'pT'}\] ";
    $l .= "$p{'dbID'}\[$p{'dbLo'}-$p{'dbHi'}\] ";
    $l .= "<$p{'numMatches'}-$p{'numMatchesN'}-$p{'percentID'}-$p{'matchOrientation'}-$p{'strandPrediction'}>\n";

    foreach my $x (@L) {
        $l .= "$x\n";
    }

    return($l);
}


######################################################################
#
#  Subroutine to read a single sim4 polish, and return it as a structure.
#
sub readPolish (*) {
    local *READPOLISHFH = shift;
    my %p;
    my $line;
    my $save;

    #  These are the fields returned.
    #
    $p{'raw'} = undef;

    $p{'estID'} = undef;
    $p{'estDefLine'} = undef;
    $p{'estLen'} = undef;
    $p{'pA'} = undef;
    $p{'pT'} = undef;

    $p{'dbID'} = undef;
    $p{'dbDefLine'} = undef;
    $p{'dbLen'} = undef;
    $p{'dbLo'} = undef;
    $p{'dbHi'} = undef;

    $p{'comment'} = undef;

    $p{'numMatches'} = undef;
    $p{'numMatchesN'} = undef;
    $p{'percentID'} = undef;
    $p{'coverage'} = undef;
    $p{'matchOrientation'} = undef;
    $p{'strandPrediction'} = undef;

    #  An array of references to hashes, one hash for each exon.
    $p{'exons'} = ();


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
    return(%p) if (eof(READPOLISHFH));


    #  Read the description line
    #
    $line  = <READPOLISHFH>;
    $save .= $line;

    if ($line =~ m/^(\d+)\[(\d+)-+(\d+)-+(\d+)\]\s+(\d+)\[(\d+)-(\d+)\]\s+\<(\d+)-(\d+)-(\d+)-(\w+)-(\w+)\>$/) {
        $p{'estID'}            = $1;
        $p{'estLen'}           = $2;
        $p{'pA'}               = $3;
        $p{'pT'}               = $4;
        $p{'dbID'}             = $5;
        $p{'dbLo'}             = $6;
        $p{'dbHi'}             = $7;
        $p{'numMatches'}       = $8;
        $p{'numMatchesN'}      = $9;
        $p{'percentID'}        = $10;
        $p{'matchOrientation'} = $11;
        $p{'strandPrediction'} = $12;
    } else {
        print STDERR "expecting description line, got: '$line'\n";
        return(%p);
    }


    #  Read the two deflines, if they exist.
    #
    $line = <READPOLISHFH>;

    if ($line =~ m/^comment=\s*(.*)\s*$/) {
        $p{'comment'} = $1;
        $save .= $line;
        $line  = <READPOLISHFH>;
    } else {
        #print STDERR "libBri::readPolish()-- WARNING:  Didn't get comment!\n";
        #print STDERR "libBri::readPolish()-- WARNING:  $line";
    }
    if ($line =~ m/^edef=(.*)$/) {
        $p{'estDefLine'} = $1;
        $save .= $line;
        $line  = <READPOLISHFH>;
    } else {
        #print STDERR "libBri::readPolish()-- WARNING:  Didn't get edef!\n";
        #print STDERR "libBri::readPolish()-- WARNING:  $line";
    }

    if ($line =~ m/^ddef=(.*)$/) {
        $p{'dbDefLine'} = $1;
        $save .= $line;
        $line  = <READPOLISHFH>;
    } else {
        #print STDERR "libBri::readPolish()-- WARNING:  Didn't get ddef!\n";
        #print STDERR "libBri::readPolish()-- WARNING:  $line";
    }


    #  Read the exons
    #
    my $exonAlign      = 0;
    my $exonAlignFirst = 1;
    my $exonCoverage   = 0;

    while (defined($line) && ($line !~ m/^sim4end$/)) {

        #  If this match succeeds, we have an exon description.
        #  Otherwise, it's an alignment line.
        #
        if ($line =~ /^(\d+)-(\d+)\s+\((\d+)-(\d+)\)\s+\<(\d+)-(\d+)-(\d+)\>\s+(.*)$/) {
            my $e = {};

            $exonCoverage += $2 - $1 + 1;

            $e->{'cDNAstart'}         = $1;
            $e->{'cDNAend'}           = $2;
            $e->{'GENOMICstart'}      = $3;
            $e->{'GENOMICend'}        = $4;
            $e->{'numMatches'}        = $5;
            $e->{'numMatchesN'}       = $6;
            $e->{'percentID'}         = $7;
            $e->{'intronOrientation'} = $8;

            push @{$p{'exons'}}, $e;
        } else {
            if ($exonAlignFirst) {
                $p{'exons'}[$exonAlign]->{'cDNAalign'} = $line;
                chomp $p{'exons'}[$exonAlign]->{'cDNAalign'};
                $exonAlignFirst = 0;
            } else {
                $p{'exons'}[$exonAlign]->{'GENOMICalign'} = $line;
                chomp $p{'exons'}[$exonAlign]->{'GENOMICalign'};
                $exonAlignFirst = 1;
                $exonAlign++;
            }
        }

        $save .= $line;
        $line  = <READPOLISHFH>;
    }

    $save .= $line;

    if (($p{'pA'} + $p{'pT'}) >= $p{'estLen'}) {
        $p{'coverage'} = 0;
    } else {
        $p{'coverage'} = 100.0 * $exonCoverage / ($p{'estLen'} - $p{'pA'} - $p{'pT'});
    }

    $p{'raw'} = $save;

    return(%p);
}



######################################################################
#
#  Functions for running multiple processes at the same time.
#
my $numberOfProcesses       = 0;
my $numberOfProcessesToWait = 0;
my @processQueue;
my @processesRunning;
my $printProcessCommand = 0;
my $printProcessStatus  = 0;

sub schedulerSetNumberOfProcesses {
    ($numberOfProcesses) = @_;
}

sub schedulerSetNumberOfProcessesToWaitFor {
    ($numberOfProcessesToWait) = @_;
}

sub schedulerSetShowCommands {
    ($printProcessCommand) = @_;
}

sub schedulerSetShowStatus {
    ($printProcessStatus) = @_;
}


#  Submit a task to the scheduler
#
sub schedulerSubmit {
    chomp @_;
    push @processQueue, @_;
}

sub forkProcess {
    my($process) = @_;
    my($pid);

    #  From Programming Perl, page 167
  FORK: {
      if ($pid = fork) {
          # Parent
          #
          return($pid);
     } elsif (defined $pid) {
         # Child
         #
         exec($process);
      } elsif ($! =~ /No more processes/) {
          # EAGIN, supposedly a recoverable fork error
          sleep 1;
          redo FORK;
      } else {
          die "Can't fork: $!\n";
      }
  }
}

sub reapProcess {
    my($pid) = @_;

    if (waitpid($pid, &WNOHANG) > 0) {
        return(1);
    } else {
        return(0);
    }
}

sub schedulerRun {
    my(@newProcesses);

    #  Reap any processes that have finished
    #
    undef @newProcesses;
    foreach my $i (@processesRunning) {
        if (reapProcess($i) == 0) {
            push @newProcesses, $i;
        }
    }
    undef @processesRunning;
    @processesRunning = @newProcesses;

    #  Run processes in any available slots
    #
    while (((scalar @processesRunning) < $numberOfProcesses) &&
           ((scalar @processQueue) > 0)) {
        my $process = shift @processQueue;

        if ($printProcessCommand) {
            print "sched()-- starting '$process'";
        }

        push @processesRunning, forkProcess($process);

        if ($printProcessStatus) {
            my $remain = scalar(@processQueue);
            my $prefix;

            if ($printProcessCommand) {
                $prefix = " -- ";
            } else {
                $prefix = "sched()-- ";
            }

            if ($remain == 0) {
                print "${prefix}No jobs remain in the queue.\n";
            } elsif ($remain == 1) {
                print "${prefix}1 job remains in the queue.\n";
            } else {
                print "${prefix}$remain jobs remain in the queue.\n";
            }
        } elsif ($printProcessCommand) {
            print "\n";
        }
    }
}


#  Wait for all processes in the scheduler to finish.
#
sub schedulerFinishStatusReport {
    my ($remain) = @_;

}

sub schedulerFinish {
    my $child;
    my @newProcesses;
    my $remain;

    $remain = scalar @processQueue;

    #  Run all submitted jobs
    #
    while ($remain > 0) {
        schedulerRun();

        $remain = scalar @processQueue;

        if ($remain > 0) {
            $child = waitpid -1, 0;

            undef @newProcesses;
            foreach my $i (@processesRunning) {
                push @newProcesses, $i if ($child != $i);
            }
            undef @processesRunning;
            @processesRunning = @newProcesses;
        }
    }

    if ($printProcessStatus) {
        print "sched()-- All jubs submitted.  Waiting for completion.\n";
    }

    #  Wait for them to finish, if requested
    #
    while ((scalar @processesRunning) > $numberOfProcessesToWait) {
        if ($printProcessStatus) {
            my $remain = scalar(@processesRunning);

            if ($remain == 0) {
                print "sched()-- No jobs running.\n";
            } elsif ($remain == 1) {
                print "sched()-- 1 job running.\n";
            } else {
                print "sched()-- $remain jobs running.\n";
            }
        }

        waitpid(shift @processesRunning, 0);
    }

    if ($printProcessStatus) {
        print "sched()-- All done!\n";
    }
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
            my ($h, $s) = nextFastA(*F);
            $deflines{$h} = 1;
        }
        close(F);
    }

    if (-e $discard2) {
        open(F, "< $discard2");
        while (!eof(F)) {
            my ($h, $s) = nextFastA(*F);
            $deflines{$h} = 1;
        }
        close(F);
    }

    #  Now, copy $input to $output
    #
    open(I, "< $input");
    open(O, "> $output");
    while (!eof(I)) {
        my ($h, $s) = nextFastA(*I);
        print O ">$h\n$s\n" if (!defined($deflines{$h}));
    }
    close(O);
    close(I);
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
        my ($h, $s) = nextFastA(*I);
        print O ">$h\n$s\n" if (0 == int <F>);
    }
    close(O);
    close(I);
    close(F);
}


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
        my ($h, $s) = nextFastA(*F);

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

######################################################################
#
#  Filters a set of polishes according to coverage and identity requirements.
#  Things that pass are put into $save.
#
#  Percentages should be specified as integers, e.g., 50, not 0.5
#
#sub filterPolishes {
#    my ($polishes, $coverage, $identity, $save) = @_;
#}


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

1;
