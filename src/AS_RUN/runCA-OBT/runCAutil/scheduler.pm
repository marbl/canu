#!/usr/local/bin/perl

#                    Confidential -- Do Not Distribute
#   Copyright (c) 2002 PE Corporation (NY) through the Celera Genomics Group
#                           All Rights Reserved.

package scheduler;

use strict;
use POSIX "sys_wait_h";

$| = 1;

#  Called by "use scheduler;"
sub import () {
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


1;
