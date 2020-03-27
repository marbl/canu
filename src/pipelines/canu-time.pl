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

my $asmPath = shift @ARGV;
my $gridEngine;

my @files;

#  Find the list of canu log files.
open(F, "ls $asmPath/canu-scripts/canu.*.out |");
@files = <F>;
chomp @files;
close(F);

my @jobNames;
my @jobIDs;

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
            push @jobNames, $1;
            push @jobIDs,   $2;

            #print STDERR " $2 - $1\n";
        }

        #  Look for 'canu-scripts/canu.##.sh', and then parse the next line
        #  for the job ID.  Ugly.
        if (m/canu-scripts/) {
            $_ = <F>;  chomp;

            if ($_ =~ m/^Your\s+job\s+(\d+)\s/) {
                push @jobNames, "canu";    #  SGE
                push @jobIDs,   $1;
            } else {
                push @jobNames, "canu";    #  Slurm
                push @jobIDs,   $_;
            }

            #print STDERR " $_ - canu\n";
        }
    }
    close(F);
}

#  Grab accounting for all these
my $totalUser    = 0;
my $totalSystem  = 0;
my $totalWall    = 0;
my $maxRSS       = 0;
my $maxVM        = 0;
my $totalStart   = 9999999999;
my $totalEnd     = 0;

printf STDOUT "user-time  sys-time wall-time   max-rss    max-vm   elapsed                 start-time                   end-time    jobID stage\n";
printf STDOUT "  (hours)   (hours)   (hours)      (MB)      (MB)   (hours)\n";
printf STDOUT "--------- --------- --------- --------- --------- --------- -------------------------- -------------------------- -------- --------------------\n";

while (scalar(@jobNames) > 0) {
    my ($jobUser, $jobSystem, $jobWall, $jobRSS, $jobVM, $jobStart, $jobEnd);

    my $jobName = shift @jobNames;
    my $jobID   = shift @jobIDs;

    #print STDERR "\n";
    #print STDERR "jobName - $jobName\n";
    #print STDERR "jobID   - $jobID\n";

    if      ($gridEngine eq "SGE") {
        ($jobUser, $jobSystem, $jobWall, $jobRSS, $jobVM, $jobStart, $jobEnd) = getTimeSGE($jobID);
    } elsif ($gridEngine eq "SLURM") {
        ($jobUser, $jobSystem, $jobWall, $jobRSS, $jobVM, $jobStart, $jobEnd) = getTimeSlurm($jobID);
    } else {
        die "Unknown grid engine '$gridEngine'.\n";
    }

    $totalUser    += $jobUser;
    $totalSystem  += $jobSystem;
    $totalWall    += $jobWall;
    $maxRSS        = ($maxRSS     < $jobRSS) ? $jobRSS : $maxRSS;
    $maxVM         = ($maxVM      < $jobVM)  ? $jobVM  : $maxVM;
    $totalStart    = ($totalStart < $jobStart) ? $totalStart : $jobStart;
    $totalEnd      = ($totalEnd   > $jobEnd)   ? $totalEnd   : $jobEnd  ;

    my $st = gmtime($jobStart);
    my $et = gmtime($jobEnd);

    my $jobElapsed = ($jobEnd - $jobStart) / 3600.0;

    printf STDOUT "%9.2f %9.2f %9.2f %9.3f %9.3f %9.2f %26s %26s %8d %s\n", $jobUser, $jobSystem, $jobWall, $jobRSS, $jobVM, $jobElapsed, $st, $et, $jobID, $jobName;
}

my $st = gmtime($totalStart);
my $et = gmtime($totalEnd);

my $totalElapsed = ($totalEnd - $totalStart) / 3600.0;

printf STDOUT "--------- --------- --------- --------- --------- --------- -------------------------- -------------------------- -------- --------------------\n";
printf STDOUT "%9.2f %9.2f %9.2f %9.3f %9.3f %9.2f %26s %26s\n", $totalUser, $totalSystem, $totalWall, $maxRSS, $maxVM, $totalElapsed, $st, $et;




sub getTimeSGE ($) {
    my $jobID = shift @_;
    my $ju    = 0;
    my $js    = 0;
    my $jr    = 0;
    my $jv    = 0;
    my $je    = 0;

    my $st    = time();
    my $et    = 0;

    open(F, "qacct -j $jobID |");
    while (<F>) {
        s/^\s+//;
        s/\s+$//;

        my @v   = split '\s+', $_;
        my $tag = shift @v;
        my $val = join ' ', @v;

        if    ($tag eq "start_time") {
            my $t = parseDate($val);

            $st = ($st < $t) ? $st : $t;
        }

        elsif ($tag eq "end_time") {
            my $t = parseDate($val);

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

    return($ju, $js, $je, $jr, $jv, $st, $et);
}



sub getTimeSlurm ($) {
    my $jobID = shift @_;
    my $ju    = 0;
    my $js    = 0;
    my $jr    = 0;
    my $jv    = 0;
    my $je    = 0;
    my $st    = 9999999999;
    my $et    = 0;

    open(F, "sacct -n -o 'JobID%30,UserCPU,SystemCPU,TotalCPU,MaxRSS,MaxVMSize,CPUTimeRAW,Elapsed,Start,End' -j $jobID |");
    while (<F>) {
        s/^\s+//;
        s/\s+$//;

        my @v = split '\s+', $_;
        my $t;

        next if ($v[0] !~ m/batch/);

        $ju += parseTime($v[1]);
        $js += parseTime($v[2]);
        $jr  = ($v[4] < $jr) ? $jr : $v[4];
        $jv  = ($v[5] < $jv) ? $jv : $v[5];
        $je += parseTime($v[7]);

        $t = parseDate($v[8]);
        $st = ($st < $t) ? $st : $t;

        $t = parseDate($v[9]);
        $et = ($t < $et) ? $et : $t;
    }
    close(F);

    $jr /= 1024.0;
    $jv /= 1024.0;

    return($ju, $js, $je, $jr, $jv, $st, $et);
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
sub parseDate ($) {
    my $d = $_[0];
    my $t;

    if ($d =~ m/^\d+-\d+-\d/) {
        $t = Time::Piece->strptime($d, "%Y-%m-%dT%H:%M:%S");
    } else {
        $t = Time::Piece->strptime($d, "%a %b %e %H:%M:%S %Y");
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
