#!/usr/bin/perl

use strict;
use lib "/bioinfo/assembly/walenz/src/genomics/scripts";
use scheduler;

#  Submit a set of jobs using the smallest number of job arrays.
#  Returns the number of jobs submitted -- zero if no jobs were
#  submitted because all jobs are done.
#
#  It is assumed that your job numbering starts 1, not 0.  If you have
#  jobs that need a parameter 0 <= x < jobid, it's up to your script
#  to subtract 1 from the jobid supplied.
#
#  Here's the block of code I use to get either the SGE_TASK_ID, or if
#  that is not set, then then command line, and then subtract one.
#
#        print F "jobid=\$SGE_TASK_ID\n";
#        print F "if [ x\$jobid = x -o x\$jobid = xundefined ]; then\n";
#        print F "  jobid=\$1\n";
#        print F "fi\n";
#        print F "if [ x\$jobid = x ]; then\n";
#        print F "  echo Error: I need SGE_TASK_ID set, or a job index on the command line.\n";
#        print F "  exit 1\n";
#        print F "fi\n";
#        print F "\n";
#        print F "jobid=`expr \$jobid - 1`\n";
#        print F "jobid=`printf %03d \$jobid`\n";
#
#  The number of #'s in the output template indicates the actual number
#  of digits in it.
#
#  You probably want to use full absolute paths for everything.

#runBatchSGE("/bioinfo/assembly/walenz/bigmerjunk/1-markDistinct",
#            16,
#            "/bioinfo/assembly/walenz/bigmerjunk/1-markDistinct/pass1-###.markbits.success",
#            "/bioinfo/assembly/walenz/bigmerjunk/1-markDistinct/mark.sh",
#            "-N bigM1");
            
sub runBatchSGE ($$$$$) {
    my $path    = shift @_;  #  Directory to run in
    my $nJobs   = shift @_;  #  Number of jobs to run/check
    my $output  = shift @_;  #  Final output template created (e.g., "stuff-#####.success")
    my $script  = shift @_;  #  Script to run
    my $sge     = shift @_;  #  General SGE options

    system("mkdir $path/sgeout") if (! -d "$path/sgeout");

    #  Build our prototype job number.
    my $numberc = ($output =~ tr/\#//);
    my $number0 = "0" x $numberc;
    my $numberh = "#" x $numberc;

    #  Check if we need to submit pieces of the array, or if we can
    #  submit the whole thing.

    my @ap;
    my @sf;
    my $wholeThing = 1;  #  true, just submit the whole batch

    my $seg = "1";
    $seg = substr("$number0$seg", -$numberc);
    while ($seg <= $nJobs) {
        $seg = substr("$number0$seg", -$numberc);

        my $successfile = $output;
        $successfile =~ s/$numberh/$seg/;

        push @sf, $successfile;

        if (-e $successfile) {
            $wholeThing = 0;
        } else {
            push @ap, $seg;
        }

        $seg++;
    }


    #  If we are told to run right now (not submit to the grid) then
    #  run right now.
    #
    if ($sge =~ m/local.*(\d+)/) { 
        &scheduler::schedulerSetNumberOfProcesses($1);
        &scheduler::schedulerSetShowCommands(1);

        foreach my $seg (@ap) {
            my $cmd = "cd $path && $script $seg > $path/sgeout/seg$seg.out 2>&1";
            &scheduler::schedulerSubmit($cmd);
        }

        &scheduler::schedulerFinish();

        #  Return if all jobs all finished successfully.

        my $failed = 0;
        foreach my $successfile (@sf) {
            if (! -e $successfile) {
                print STDERR "Didn't find '$successfile]' -- JOB FAILED?\n";
                $failed++;
            }
        }
        die "Some jobs didn't finish successfully.  Bye.\n" if ($failed);
        return(0);
    }


    my $nSubmitted = 0;

    if ($wholeThing == 1) {
        #  Yippee!  Submit all at once!
        #
        my $cmd;
        $cmd  = "cd $path && ";
        $cmd .= "qsub ";
        $cmd .= "-cwd -j y -o $path/sgeout/seg\\\$TASK_ID.out ";
        $cmd .= "-t 1-$nJobs $sge $script";

        print STDERR "$cmd\n";

        die "SGE submission failed?\n" if (runCommand($cmd));

        $nSubmitted = $nJobs;
    } elsif (scalar(@ap) > 0) {
        #  Dang, we need to submit individually....or we can take five
        #  minutes and figure out ranges to submit.
        #
        my $st;
        my $ed;

        #  To catch the last chunk, we need to know where we fall off
        #  the end...so $it is a global.

        my $it = undef;

        while (scalar(@ap) > 0) {
            $it = shift @ap;

            if (!defined($st)) {
                $st = $it;
                $ed = $it + 1;
                undef $it;
            } elsif ($ed == $it) {
                $ed = $it + 1;
                undef $it;
            } else {
              submitagain:
                #print STDERR "submit $st - $ed\n";
                my $cmd;
                $cmd  = "cd $path && ";
                $cmd .= "qsub ";
                $cmd .= "-cwd -j y -o $path/sgeout/seg\\\$TASK_ID.out ";
                $cmd .= "-t $st-$ed $sge $script";

                print STDERR "$cmd\n";

                die "SGE submission failed?\n" if (runCommand($cmd));

                $nSubmitted = $ed - $st + 1;

                $st = $it;
                $ed = $it + 1;
                undef $it;
            }
        }

        if (defined($it)) {
            $ed = $it;
            goto submitagain;
        }
    } else {
        print STDERR "All segments computed successfully!\n";
    }

    print STDERR "Submitted $nSubmitted jobs.\n" if ($nSubmitted > 0);
    return($nSubmitted);
}

1;
