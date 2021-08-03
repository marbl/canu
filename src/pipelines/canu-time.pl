#!/usr/bin/env perl

###############################################################################
 #
 #  This file is part of canu, a software program that assembles whole-genome
 #  sequencing reads into contigs.
 #
 #  This software is based on:
 #    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 #    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 #
 #  Except as indicated otherwise, this is a 'United States Government Work',
 #  and is released in the public domain.
 #
 #  File 'README.licenses' in the root directory of this distribution
 #  contains full conditions and disclaimers.
 ##

use strict;
use Time::Piece;

my $asmPath;
my $gridEngine;

my %ignore;

my @jobNames;
my @jobIDs;


while (scalar (@ARGV > 0)) {
    my $arg = shift @ARGV;

    if    ($arg eq "-ignore") {
        $ignore{shift @ARGV} = 1;
    }

    elsif ($arg eq "-job") {
        my $jobID = shift @ARGV;

        push @jobNames, "";
        push @jobIDs,   $jobID;
    }

    elsif (!defined($asmPath)) {
        $asmPath = $arg;
    }

    else {
        die "usage: $0 [-ignore <jobID>] [-job <jobID>] <path-to-assembly>\n";
    }
}

if (scalar(@jobIDs) == 0) {
    findCanuJobs($asmPath);
}

else {
    detectGridEngine($asmPath);
}


#  Grab accounting for all these
my $totalUser    = 0;
my $totalSystem  = 0;
my $totalWall    = 0;
my $maxRSS       = 0;
my $maxVM        = 0;
my $totalStart   = 9999999999;
my $totalEnd     = 0;
my $totalElapsed = 0;

printf STDOUT "         -----sum-of-per-job-times---- ------per-job------ ---------------------------per-stage---------------------------\n";
printf STDOUT "jobID    user-time  sys-time wall-time   max-rss    max-vm                 start-time                   end-time elapsed  stage\n";
printf STDOUT "           (hours)   (hours)   (hours)      (MB)      (MB)                                                       (hours)\n";
printf STDOUT "-------- --------- --------- --------- --------- --------- -------------------------- -------------------------- -------  --------------------\n";

while (scalar(@jobIDs) > 0) {
    my ($jobUser, $jobSystem, $jobWall, $jobRSS, $jobVM, $jobStart, $jobEnd, $ignore);

    my $jobName = shift @jobNames;
    my $jobID   = shift @jobIDs;

    #print STDERR "\n";
    #print STDERR "jobName - $jobName\n";
    #print STDERR "jobID   - $jobID\n";

    if      ($gridEngine eq "SGE") {
        ($jobName, $jobUser, $jobSystem, $jobWall, $jobRSS, $jobVM, $jobStart, $jobEnd, $ignore) = getTimeSGE($jobName, $jobID);
    } elsif ($gridEngine eq "SLURM") {
        ($jobName, $jobUser, $jobSystem, $jobWall, $jobRSS, $jobVM, $jobStart, $jobEnd, $ignore) = getTimeSlurm($jobName, $jobID);
    } else {
        die "Unknown grid engine '$gridEngine'.\n";
    }

    next if ($ignore);

    my $st = gmtime($jobStart);
    my $et = gmtime($jobEnd);

    my $jobElapsed = ($jobEnd - $jobStart) / 3600.0;

    $totalUser    += $jobUser;
    $totalSystem  += $jobSystem;
    $totalWall    += $jobWall;
    $maxRSS        = ($maxRSS     < $jobRSS) ? $jobRSS : $maxRSS;
    $maxVM         = ($maxVM      < $jobVM)  ? $jobVM  : $maxVM;
    $totalStart    = ($totalStart < $jobStart) ? $totalStart : $jobStart;
    $totalEnd      = ($totalEnd   > $jobEnd)   ? $totalEnd   : $jobEnd  ;
    $totalElapsed += $jobElapsed;

    printf STDOUT "%-8d %9.2f %9.2f %9.2f %9.3f %9.3f %26s %26s %7.2f  %s\n", $jobID, $jobUser, $jobSystem, $jobWall, $jobRSS, $jobVM, $st, $et, $jobElapsed, $jobName;
}

my $st = gmtime($totalStart);
my $et = gmtime($totalEnd);

my $totalInQueue = ($totalEnd - $totalStart) / 3600.0;

printf STDOUT "-------- --------- --------- --------- --------- --------- -------------------------- -------------------------- -------  --------------------\n";
printf STDOUT "         %9.2f %9.2f %9.2f %9.3f %9.3f %26s %26s %7.2f\n", $totalUser, $totalSystem, $totalWall, $maxRSS, $maxVM, $st, $et, $totalElapsed;
printf STDOUT "\n";
printf STDOUT "Execution time: %s.\n", reportTimeDHM($totalElapsed);
printf STDOUT "Queue time:     %s.\n", reportTimeDHM($totalInQueue);
printf STDOUT "CPU time:       %s.\n", reportTimeDHM($totalUser + $totalSystem);
printf STDOUT "\n";





sub findCanuJobs ($) {
    my $asmPath = shift @_;
    my @files;

    #  Find the list of canu log files.
    open(F, "ls $asmPath/canu-scripts/canu.*.out |");
    @files = <F>;
    chomp @files;
    close(F);

    if (scalar(@files) == 0) {
        die "Didn't find any canu.*.out files in $asmPath/canu-scripts/.\n";
    }


    #  Parse log files, looking for submitted jobs.
    #    -- 'meryl-count.jobSubmit-01.sh' -> job 21767697 task 1.
    foreach my $file (@files) {
        #print STDERR "Scanning '$file'.\n";

        open(F, "< $file");
        while (! eof(F)) {
            $_ = <F>;

            if (m/--\sDetected\sSun\sGrid\sEngine\s/) {
                $gridEngine = "SGE";
            }

            if (m/--\sDetected\sSlurm\s/) {
                $gridEngine = "SLURM";
            }

            if (m/--\s'(.*).jobSubmit-\d\d.sh'\s->\sjob\s(\d+)\s/) {
                if (!exists($ignore{$2})) {
                    push @jobNames, $1;
                    push @jobIDs,   $2;
                }

                #print STDERR " $2 - $1\n";
            }

            #  Look for 'canu-scripts/canu.##.sh', and then parse the next line
            #  for the job ID.  Ugly.
            if (m/canu-scripts/) {
                $_ = <F>;  chomp;

                if ($_ =~ m/^Your\s+job\s+(\d+)\s/) {
                    if (!exists($ignore{$1})) {
                        push @jobNames, "canu";    #  SGE
                        push @jobIDs,   $1;
                    }
                } else {
                    if (!exists($ignore{$_})) {
                        push @jobNames, "canu";    #  Slurm
                        push @jobIDs,   $_;
                    }
                }

                #print STDERR " $_ - canu\n";
            }
        }
        close(F);
    }
}



sub detectGridEngine ($) {
    my $asmPath = shift @_;

    if (exists($ENV{'SGE_CELL'})) {
        $gridEngine = "SGE";
    } else {
        $gridEngine = "SLURM";
    }
}




sub reportTimeDHM ($) {
    my $t = shift @_;

    my $d = int(($t) / 24);
    my $h = int(($t - 24 * $d));
    my $m = int(($t - 24 * $d - $h) * 60);

    return(sprintf "%3d day%s %2d hour%s %2d minute%s", $d, ($d == 1) ? " " : "s", $h, ($h == 1) ? " " : "s", $m, ($m == 1) ? " " : "s");
}


sub getTimeSGE ($) {
    my $jobName = shift @_;
    my $jobID   = shift @_;
    my $jn      = $jobName;
    my $ju      = 0;
    my $js      = 0;
    my $jr      = 0;
    my $jv      = 0;
    my $je      = 0;
    my $ST      = time();
    my $st      = $ST;
    my $ET      = 0;
    my $et      = 0;
    my $ig      = 0;

    open(F, "qacct -j $jobID 2> /dev/null |");
    while (<F>) {
        s/^\s+//;
        s/\s+$//;

        my @v   = split '\s+', $_;
        my $tag = shift @v;
        my $val = join ' ', @v;

        if    (($tag eq "jobname") && ($jobName eq "")) {
            $jn = $val;
        }

        elsif ($tag eq "start_time") {
            my $t = parseDate($val, $st);
            $st = ($st < $t) ? $st : $t;
        }

        elsif ($tag eq "end_time") {
            my $t = parseDate($val, $et);
            $et = ($t < $et) ? $et : $t;
        }

        elsif ($tag eq "ru_wallclock") {
            $je += parseTime($val);
        }

        elsif ($tag eq "ru_utime") {
            $ju += parseTime($val);
        }

        elsif ($tag eq "ru_stime") {
            $js += parseTime($val);
        }

        elsif ($tag eq "ru_maxrss") {
            my $r = parseMemory($val);

            $jr = ($r < $jr) ? $jr : $r;
        }

        elsif ($tag eq "cpu") {
        }

        elsif ($tag eq "maxvmem") {
            my $v = parseMemory($val);

            $jv = ($v < $jv) ? $jv : $v;
        }
    }
    close(F);

    $ig = 1   if ($st == $ST);   #  Ignore this result if either
    $ig = 1   if ($et == $ET);   #  time was missing.

    return($jn, $ju, $js, $je, $jr, $jv, $st, $et, $ig);
}



sub getTimeSlurm ($) {
    my $jobName = shift @_;
    my $jobID   = shift @_;
    my $jn      = $jobName;
    my $ju      = 0;
    my $js      = 0;
    my $jr      = 0;
    my $jv      = 0;
    my $je      = 0;
    my $st      = 9999999999;
    my $et      = 0;
    my $ig      = 0;

    open(F, "sacct -n -o 'JobID%30,UserCPU,SystemCPU,TotalCPU,MaxRSS,MaxVMSize,CPUTimeRAW,Elapsed,Start,End' -j $jobID 2> /dev/null |");
    while (<F>) {
        s/^\s+//;
        s/\s+$//;

        my @v = split '\s+', $_;
        my $t;

        next if ($v[0] !~ m/batch/);

        $jn = $v[0]   if ($jn eq "");

        $ju += parseTime($v[1]);
        $js += parseTime($v[2]);
        $jr  = ($v[4] < $jr) ? $jr : $v[4];
        $jv  = ($v[5] < $jv) ? $jv : $v[5];
        $je += parseTime($v[7]);

        $t = parseDate($v[8], $st);
        $st = ($st < $t) ? $st : $t;

        $t = parseDate($v[9], $et);
        $et = ($t < $et) ? $et : $t;
    }
    close(F);

    $jr /= 1024.0;
    $jv /= 1024.0;

    return($jn, $ju, $js, $je, $jr, $jv, $st, $et, $ig);
}




sub parseTime ($) {
    my $t = $_[0];
    my $s = 0;

    if      ($t =~ m/^(\d+)-(\d+):(\d+):(\d+)$/) {
        $s += $1 * 24 * 60 * 60;
        $s += $2 *      60 * 60;
        $s += $3 *           60;
        $s += $4;

    } elsif ($t =~ m/^(\d+):(\d+):(\d+)$/) {
        $s += $1 *      60 * 60;
        $s += $2 *           60;
        $s += $3;

    } elsif ($t =~ m/^(\d+):(\d+).(\d+)$/) {
        $s += $1 *           60;
        $s += $2;
        $s += $3 / 100;

    } elsif ($t =~ m/(\d+.*\d*)s/) {
        $s += $1;

    } else {
        die "Failed to parse time '$t'\n";
    }

    return($s / 3600.0);
}



#  SGE returns Mon Mar 18 10:32:31 2019
#              %a  %b  %e %H:%M:%S %Y
#  Slurm returns YYYY-MM-DDTHH:MM:SS
#                %Y-%m-%dT%H:%M:%S
sub parseDate ($$) {
    my $d = $_[0];   #  Input date/time string.
    my $f = $_[1];   #  If fail to parse, return this instead.
    my $t;

    if      ($d =~ m/^\d+-\d+-\d/) {
        $t = Time::Piece->strptime($d, "%Y-%m-%dT%H:%M:%S");
    } elsif ($d =~ m/\s\d+:\d+:\d+\s\d\d\d\d$/) {
        $t = Time::Piece->strptime($d, "%a %b %e %H:%M:%S %Y");
    } else {
        return($f);
    }

    #my $T = gmtime($t->epoch);
    #print "$d -> $t -> $T\n";

    return($t->epoch);
}


sub parseMemory ($) {
    my $m = $_[0];

    if    ($m =~ m/^(\d+.*\d*)\s*GB$/) {
        $m = $1 / 1024.0;
    }

    elsif ($m =~ m/^(\d+.*\d*)\s*MB$/) {
        $m = $1;
    }

    elsif ($m =~ m/^(\d+.*\d*)\s*KB$/) {
        $m = $1 / 1024.0;
    }

    elsif ($m =~ m/^(\d+.*\d*)\s*B$/) {
        $m = $1 / 1024.0 / 1024.0;
    }

    else {
        die "Failed to parse memory '$m'\n";
    }

    return($m);
}
