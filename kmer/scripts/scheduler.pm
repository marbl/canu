#!/usr/local/bin/perl
#
#  Functions for running multiple processes at the same time.
#

package scheduler;

use strict;
use POSIX "sys_wait_h";

#  Called by "use scheduler;"
sub import () {
}

my $numberOfProcesses = 0;
my @processQueue      = ();

sub schedulerSetNumberOfProcesses {
    $numberOfProcesses = shift @_;
}

sub schedulerSubmit {
    chomp @_;
    push @processQueue, @_;
}

sub forkProcess {
    my $process = shift @_;
    my $pid;

    #  From Programming Perl, page 167
  FORK:
    if ($pid = fork) {
        return($pid);    #  Parent, returns child id
    } elsif (defined $pid) {
        exec($process);  #  Child, runs the process
    } elsif ($! =~ /No more processes/) {
        sleep 1;         # EAGIN, supposedly a recoverable fork error
        redo FORK;
    } else {
        die "Can't fork: $!\n";
    }

    die "scheduler::forkProcess()--  Shouldn't be here.\n";
}

sub schedulerFinish {
    my @processesRunning;
    my @newProcesses;
    my $remain = scalar(@processQueue);

    my $t = localtime();
    my $d = time();

    print STDERR "----------------------------------------START CONCURRENT $t\n";

    while ($remain > 0) {

        #  Reap any processes that have finished

        undef @newProcesses;
        foreach my $i (@processesRunning) {
            if (waitpid($i, &WNOHANG) <= 0) {
                push @newProcesses, $i;
            }
        }
        undef @processesRunning;
        @processesRunning = @newProcesses;

        #  Run processes in any available slots

        while ((scalar(@processesRunning) < $numberOfProcesses) &&
               (scalar(@processQueue) > 0)) {
            my $process = shift @processQueue;
            print STDERR "$process\n";
            push @processesRunning, forkProcess($process);
        }

        $remain = scalar(@processQueue);

        #  If still stuff out there, wait for something to finish.

        if ($remain > 0) {
            my $child = waitpid -1, 0;

            undef @newProcesses;
            foreach my $i (@processesRunning) {
                push @newProcesses, $i if ($child != $i);
            }
            undef @processesRunning;
            @processesRunning = @newProcesses;
        }
    }

    while (scalar(@processesRunning) > 0) {
        waitpid(shift @processesRunning, 0);
    }

    $t = localtime();
    print STDERR "----------------------------------------END CONCURRENT $t (", time() - $d, " seconds)\n";
}

1;
