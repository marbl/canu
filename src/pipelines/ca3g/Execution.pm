package ca3g::Execution;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(touch getInstallDirectory getBinDirectory getBinDirectoryShellCode submitScript buildGridJob submitOrRunParallelJob runCommand stopBefore stopAfter diskSpace);

use strict;
use Config;            #  for @signame

use POSIX ":sys_wait_h";  #  For waitpid(..., &WNOHANG)

#use lib "$FindBin::RealBin";
#use lib "$FindBin::RealBin/ca3g/lib/perl5";

use ca3g::Defaults;
use Filesys::Df;  #  for diskSpace()


#
#  Functions for running multiple processes at the same time.  This is private to the module.
#

my $numberOfProcesses       = 0;     #  Number of jobs concurrently running
my $numberOfProcessesToWait = 0;     #  Number of jobs we can leave running at exit
my @processQueue            = ();
my @processesRunning        = ();
my $printProcessCommand     = 1;     #  Show commands as they run

sub schedulerSetNumberOfProcesses {
    $numberOfProcesses = shift @_;
}

sub schedulerSubmit ($) {
    my $cmd = shift @_;

    chomp $cmd;

    push @processQueue, $cmd;
}

sub schedulerForkProcess ($) {
    my $process = shift @_;
    my $pid;

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

sub schedulerReapProcess ($) {
    my $pid = shift @_;

    if (waitpid($pid, &WNOHANG) > 0) {
        return(1);
    } else {
        return(0);
    }
}

sub schedulerRun () {
    my @newProcesses;

    #  Reap any processes that have finished
    #
    undef @newProcesses;
    foreach my $i (@processesRunning) {
        if (schedulerReapProcess($i) == 0) {
            push @newProcesses, $i;
        }
    }
    undef @processesRunning;
    @processesRunning = @newProcesses;

    #  Run processes in any available slots
    #
    while ((scalar(@processesRunning) < $numberOfProcesses) &&
           (scalar(@processQueue) > 0)) {
        my $process = shift @processQueue;
        print STDERR "$process\n";
        push @processesRunning, schedulerForkProcess($process);
    }
}

sub schedulerFinish ($) {
    my $dir = shift @_;
    my $child;
    my @newProcesses;
    my $remain;

    $remain = scalar(@processQueue);

    if (defined($dir)) {
        my $f = diskSpace($dir);
        print STDERR "----------------------------------------SPACE $f GB\n";
    }

    my $t = localtime();
    my $d = time();
    print STDERR "----------------------------------------START CONCURRENT $t ($remain processes; $numberOfProcesses concurrently)\n";

    #  Run all submitted jobs
    #
    while ($remain > 0) {
        schedulerRun();

        $remain = scalar(@processQueue);

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

    #  Wait for them to finish, if requested
    #
    while (scalar(@processesRunning) > $numberOfProcessesToWait) {
        waitpid(shift @processesRunning, 0);
    }

    $t = localtime();
    print STDERR "----------------------------------------END CONCURRENT $t (", time() - $d, " seconds)\n";

    if (defined($dir)) {
        my $f = diskSpace($dir);
        print STDERR "----------------------------------------SPACE $f GB DANGEROUSLY LOW\n"  if ($f < 1.0);
    }
}






#
#  File Management
#

sub touch ($@) {
    open(F, "> $_[0]") or caFailure("failed to touch file '$_[0]'", undef);
    print F "$_[1]\n"  if (defined($_[1]));
    close(F);
}





#  Decide what bin directory to use.
#
#  When we are running on the grid, the path of this perl script is NOT
#  always the correct architecture.  If the submission host is
#  FreeBSD, but the grid is Linux, the BSD box will submit
#  FreeBSD/bin/ca3g.pl to the grid -- unless it knows in advance,
#  there is no way to pick the correct one.  The grid host then has to
#  have enough smarts to choose the correct binaries, and that is what
#  we're doing here.
#
#  To make it more trouble, shell scripts need to do all this by
#  themselves.
#
sub getInstallDirectory () {
    my @t = split '/', "$FindBin::RealBin";
    pop @t;                         #  bin
    pop @t;                         #  arch, e.g., FreeBSD-amd64
    my $installDir = join '/', @t;  #  path to the assembler

    return($installDir);
}


#  Used inside ca3g to find where binaries are located.  It uses uname to find OS, architecture and
#  system name, then uses that to construct a path to binaries.  If a "pathMap" is defined, this is
#  used to hardcode a path to a system name.
#
sub getBinDirectory () {
    my $installDir = getInstallDirectory();

    my $syst = `uname -s`;    chomp $syst;  #  OS implementation
    my $arch = `uname -m`;    chomp $arch;  #  Hardware platform
    my $name = `uname -n`;    chomp $name;  #  Name of the system

    $arch = "amd64"  if ($arch eq "x86_64");
    $arch = "ppc"    if ($arch eq "Power Macintosh");

    my $path = "$installDir/$syst-$arch/bin";

    while ((! -d $path) && ($installDir ne "")) {
        print STDERR "Failed to find bin directory '$path'.\n";

        my @id = split '/', $installDir;
        pop @id;
        $installDir = join('/', @id);

        $path = "$installDir/$syst-$arch/bin";
    }

    my $pathMap = getGlobal("pathMap");
    if (defined($pathMap)) {
        open(F, "< $pathMap") or caFailure("failed to open pathMap '$pathMap'", undef);
        while (<F>) {
            my ($n, $b) = split '\s+', $_;
            $path = $b if ($name eq $n);
        }
        close(F);
    }

    if (! -d "$path") {
        print STDERR "Failed to find bin directory.\n";
        exit(1);
    }

    #print STDERR "Found binaries in '$path'\n";

    return($path);
}


#  Emits a block of shell code to locate binaries during shell scripts.  See comments on
#  getBinDirectory.
#
sub getBinDirectoryShellCode () {
    my $installDir = getInstallDirectory();
    my $string;

    $string .= "syst=`uname -s`\n";
    $string .= "arch=`uname -m`\n";
    $string .= "name=`uname -n`\n";
    $string .= "\n";
    $string .= "if [ \"\$arch\" = \"x86_64\" ] ; then\n";
    $string .= "  arch=\"amd64\"\n";
    $string .= "fi\n";
    $string .= "if [ \"\$arch\" = \"Power Macintosh\" ] ; then\n";
    $string .= "  arch=\"ppc\"\n";
    $string .= "fi\n";
    $string .= "\n";
    $string .= "bin=\"$installDir/\$syst-\$arch/bin\"\n";

    my $pathMap = getGlobal("pathMap");
    if (defined($pathMap)) {
        open(PM, "< $pathMap") or caFailure("failed to open pathMap '$pathMap'", undef);
        while (<PM>) {
            my ($n, $b) = split '\s+', $_;
            $string .= "if [ \"\$name\" = \"$n\" ] ; then\n";
            $string .= "  bin=\"$b\"\n";
            $string .= "fi\n";
        }
        close(PM);
        $string .= "\n";
    }

    return($string);
}







#  Submit ourself back to the grid.  If the one argument is defined, make us hold on jobs with that
#  name.
#
#  The previous version (CA) would use "gridPropagateHold" to reset holds on existing jobs so that
#  they would also hold on this job.
#
sub submitScript ($$$) {
    my $wrk         = shift @_;
    my $asm         = shift @_;
    my $jobToWaitOn = shift @_;

    #  If not requested to run on the grid, return.

    return if (getGlobal("useGridScript") == 0);

    #  If no job to wait on, and we are already on the grid, do NOT resubmit ourself.
    #
    #  When the user launches ca3g on the head node, a call to submitScript() is made to launch ca3g
    #  under grid control.  That results in a restart of ca3g, and another call to submitScript(),
    #  but this time, the envorinment variable is set, we we can skip the resubmission, and continue
    #  with ca3g execution.

    return if (($jobToWaitOn eq undef) && (exists($ENV{getGlobal("gridEngineJobID")})));

    #  Find the next available output file.

    my $idx = "01";

    while (-e "$wrk/ca3g.$idx.out") {
        $idx++;
    }

    my $output = "$wrk/ca3g.$idx.out";
    my $script = "$wrk/ca3g.$idx.sh";

    #  Make a script for us to submit.

    open(F, "> $script") or caFailure("failed to open '$script' for writing", undef);
    print F "#!" . getGlobal("shell") . "\n";
    print F "#\n";
    print F "if [ \"x\$SGE_ROOT\" != \"x\" ]; then \n";
    print F "   #  Attempt to (re)configure SGE.  For reasons Bri doesn't know,\n";
    print F "   #  jobs submitted to SGE, and running under SGE, fail to read his\n";
    print F "   #  .tcshrc (or .bashrc, limited testing), and so they don't setup\n";
    print F "   #  SGE (or ANY other paths, etc) properly.  For the record,\n";
    print F "   #  interactive SGE logins (qlogin, etc) DO set the environment.\n";
    print F "   \n";
    print F "   . \$SGE_ROOT/\$SGE_CELL/common/settings.sh\n";
    print F "fi\n";
    print F "\n";
    print F "#  On the off chance that there is a pathMap, and the host we\n";
    print F "#  eventually get scheduled on doesn't see other hosts, we decide\n";
    print F "#  at run time where the binary is.\n";

    print F getBinDirectoryShellCode();

    #print F "/usr/bin/env perl \$bin/ca3g " . getCommandLineOptions() . "\n";
    print F "/usr/bin/env perl $0 " . getCommandLineOptions() . "\n";
    close(F);

    system("chmod +x $script");

    #  Construct a submission command line.

    my $jobName              = "c3g_" . $asm . ((defined(getGlobal("gridOptionsJobName"))) ? ("_" . getGlobal("gridOptionsJobName")) : (""));

    my $gridOpts             = getGlobal("gridOptions") . " " . getGlobal("gridOptionsScript");

    #  If the jobToWaitOn is defined, make the script wait for that to complete.  LSF might need to
    #  query jobs in the queue and figure out the job ID (or IDs) for the jobToWaitOn.  Reading LSF
    #  docs online (for bsub.1) claim that we can still use jobToWaitOn.

    if (defined($jobToWaitOn)) {
        (my $hold = getGlobal("gridEngineHoldOption")) =~ s/WAIT_TAG/$jobToWaitOn/;
        $gridOpts .= $hold;
    }

    my $submitCommand        = getGlobal("gridEngineSubmitCommand");
    my $nameOption           = getGlobal("gridEngineNameOption");
    my $outputOption         = getGlobal("gridEngineOutputOption");

    my $qcmd = "$submitCommand $gridOpts $nameOption \"$jobName\" $outputOption $output $script";

    runCommand($wrk, $qcmd) and caFailure("Failed to submit script", undef);

    exit(0);
}






#  Expects
#    job name
#    number of jobs
#    global pattern for option
#
sub buildGridArray ($$$$) {
    my $r = $_[3];

    $r =~ s/ARRAY_NAME/$_[0]/g;        #  Replace ARRAY_NAME with 'job name'
    $r =~ s/ARRAY_JOBS/$_[1]-$_[2]/g;  #  Replace ARRAY_JOBS with 'bgn-end'

    return($r);
}



sub buildOutputName ($$$) {
    my $path   = shift @_;
    my $script = shift @_;
    my $tid    = shift @_;

    my $outName = "$path/$script.$tid.out";

    if ((-e "$path/logs") && ($script =~ m/scripts\/(.*)/)) {
        $outName = "$path/logs/$1.$tid.out";
    }

    return($outName);
}



sub buildGridJob ($$$$$$) {
    my $asm     = shift @_;
    my $jobType = shift @_;
    my $path    = shift @_;
    my $script  = shift @_;
    my $bgnJob  = shift @_;
    my $endJob  = shift @_;

    #  Unpack the job range if needed.

    if ($bgnJob =~ m/^(\d+)-(\d+)$/) {
        $bgnJob = $1;
        $endJob = $2;
    }

    if (!defined($endJob)) {
        $endJob = $bgnJob;
    }

    #  Figure out the command and options needed to run the job.

    my $submitCommand = getGlobal("gridEngineSubmitCommand");
    my $gridOpts      = getGlobal("gridOptions") . " " . getGlobal("gridOptions$jobType");
    my $nameOption    = getGlobal("gridEngineNameOption");

    my $jobNameT = "${jobType}_" . $asm . ((defined(getGlobal("gridOptionsJobName"))) ? ("_" . getGlobal("gridOptionsJobName")) : (""));

    my $jobName       = buildGridArray($jobNameT, $bgnJob, $endJob, getGlobal("gridEngineArrayName"));
    my $arrayOpt      = buildGridArray($jobNameT, $bgnJob, $endJob, getGlobal("gridEngineArrayOption"));

    my $outputOption  = getGlobal("gridEngineOutputOption");
    my $tid           = getGlobal("gridEngineArraySubmitID");
    my $outName       = buildOutputName($path, $script, $tid);

    #  Build the command line.

    my $cmd;
    $cmd  = "  $submitCommand \\\n";
    $cmd .= "    $gridOpts \\\n"  if ($gridOpts ne " ");
    $cmd .= "    $nameOption \"$jobName\" \\\n";
    $cmd .= "    $arrayOpt \\\n";
    $cmd .= "    $outputOption $outName \\\n";
    $cmd .= "    $path/$script.sh\n";

    #  Save it, just because.

    open(F, "> $path/$script.jobSubmit.sh") or die;
    print F $cmd;
    close(F);

    #  Return the command and the job name it will be submitted with.

    return($cmd, $jobName);
}




#  Convert @jobs to a list of ranges, a-b, c, d-e, etc.  These will be directly submitted to the
#  grid, or run one-by-one locally.
#
#  If we're SGE, we can combine everything to one job range: a-b,c,d,e-f.  Except that
#  buildGridJob() doesn't know how to handle that.

sub convertToJobRange (@) {
    my @jobs;

    foreach my $j (@_) {
        if        ($j =~ m/^(\d+)-(\d+)$/) {
            for (my $a=$1; $a<=$2; $a++) {
                push @jobs, $a;
            }

        } elsif ($j =~ m/^(\d+)$/) {
            push @jobs, $1;

        } else {
            caFailure("invalid job format in '$j'", undef);
        }
    }

    my @jobsA = sort { $a <=> $b } @jobs;

    undef @jobs;

    my $st = $jobsA[0];
    my $ed = $jobsA[0];

    shift @jobsA;

    foreach my $j (@jobsA) {
        if ($ed + 1 == $j) {
            $ed = $j;
        } else {
            push @jobs, ($st == $ed) ? "$st" : "$st-$ed";
            $st = $j;
            $ed = $j;
        }
    }

    push @jobs, ($st == $ed) ? "$st" : "$st-$ed";

    return(@jobs);
}





#  Expects
#    job type ("ovl", etc)
#    output directory
#    script name with no directory or .sh
#    number of jobs in the task
#
#  If under grid control, submit grid jobs.  Otherwise, run in parallel locally.
#
sub submitOrRunParallelJob ($$$$$$@) {
    my $wrk          = shift @_;
    my $asm          = shift @_;

    my $jobType      = shift @_;  #  E.g., ovl, cns, ... - populates 'gridOptionsXXX and useGridXXX (former XXXOnGrid)
                             #                         - also becomes the grid job name prefix, so three letters suggested
    my $path         = shift @_;
    my $script       = shift @_;  #  Runs $path/$script.sh > $path/$script.######.out

    my $nParallel    = shift @_;

    #my $holds        = shift @_;
    #my $submitScript = shift @_;

    my @jobs         = convertToJobRange(@_);

    #  The script MUST be executable.

    system("chmod +x \"$path/$script.sh\"");

    #  If 'gridEngineJobID' environment variable exists (SGE: JOB_ID; LSF: LSB_JOBID) then we are
    #  currently running under grid crontrol.  If so, run the grid command to submit more jobs, then
    #  submit ourself back to the grid.  If not, tell the user to run the grid command by hand.

    #  Jobs under grid control, and we submit them

    if (getGlobal("useGrid") && getGlobal("useGrid$jobType") && (exists($ENV{getGlobal("gridEngineJobID")}))) {
        my $cmd;
        my $jobName;

        foreach my $j (@jobs) {
            ($cmd, $jobName) = buildGridJob($asm, $jobType, $path, $script, $j, undef);

            runCommand($path, $cmd) and caFailure("Failed to submit batch jobs", undef);
        }

        submitScript($wrk, $asm, $jobName);

        exit(0);
    }

    #  Jobs under grid control, but the user must submit them

    if (getGlobal("useGrid") && getGlobal("useGrid$jobType") && (! exists($ENV{getGlobal("gridEngineJobID")}))) {
        print STDERR "Please submit the following jobs to the grid for execution:\n";
        print STDERR "\n";

        foreach my $j (@jobs) {
            my ($cmd, $jobName) = buildGridJob($asm, $jobType, $path, $script, $j, undef);

            print $cmd;
        }

        print STDERR "\n";
        print STDERR "When all jobs complete, restart ca3g as before.\n";

        exit(0);
    }

    #  Standard jobs, run locally.

    foreach my $j (@jobs) {
        my $st;
        my $ed;

        if ($j =~ m/^(\d+)-(\d+)$/) {
            $st = $1;
            $ed = $2;
        } else {
            $st = $ed = $j;
        }

        for (my $i=$st; $i<=$ed; $i++) {
            my $outName  = buildOutputName($path, $script, substr("000000" . $i, -6));

            schedulerSubmit("$path/$script.sh $i > $outName 2>&1");
        }
    }

    schedulerSetNumberOfProcesses($nParallel) if (defined($nParallel));
    schedulerFinish($path);
}








#  Utility to run a command and check the exit status, report time used.
#
sub runCommand ($$) {
    my $dir = shift @_;
    my $cmd = shift @_;

    if (! -d $dir) {
        caFailure("Directory '$dir' doesn't exist, can't run command", "");
    }

    if (getGlobal('showNext')) {
        print STDERR "----------------------------------------NEXT-COMMAND\n$cmd\n";
        exit(0);
    }

    if (defined($dir)) {
        my $f = diskSpace($dir);
        print STDERR "----------------------------------------SPACE $f GB\n";
    }

    my $t = localtime();
    my $d = time();
    print STDERR "----------------------------------------START $t\n$cmd\n";

    my $rc = 0xffff & system("cd $dir && $cmd");

    $t = localtime();
    print STDERR "----------------------------------------END $t (", time() - $d, " seconds)\n";

    if (defined($dir)) {
        my $f = diskSpace($dir);
        print STDERR "----------------------------------------SPACE $f GB DANGEROUSLY LOW\n"  if ($f < 1.0);
    }

    #  Pretty much copied from Programming Perl page 230

    return(0) if ($rc == 0);

    #  Bunch of busy work to get the names of signals.  Is it really worth it?!
    #
    my @signame;
    if (defined($Config{sig_name})) {
        my $i = 0;
        foreach my $n (split('\s+', $Config{sig_name})) {
            $signame[$i] = $n;
            $i++;
        }
    }

    my $error = "ERROR: Failed with ";

    if ($rc == 0xff00) {
        $error .= "$!\n";
    } else {
        if ($rc & 0x80) {
            $error .= "coredump from ";
        }

        if ($rc > 0x80) {
            $rc >>= 8;
        }
        $rc &= 127;

        if (defined($signame[$rc])) {
            $error .= "signal $signame[$rc] ($rc)\n";
        } else {
            $error .= "signal $rc\n";
        }
    }

    print STDERR $error;

    return(1);
}



sub stopBefore ($$) {
    my $stopBefore = shift @_;
    my $cmd        = shift @_;

    $stopBefore =~ tr/A-Z/a-z/;

    if ((defined($stopBefore)) &&
        (defined(getGlobal('stopBefore'))) &&
        (getGlobal('stopBefore') eq $stopBefore)) {
        print STDERR "Stop requested before '$stopBefore'.\n";
        print STDERR "Command:\n$cmd\n" if (defined($cmd));
        exit(0);
    }
}



sub stopAfter ($) {
    my $stopAfter = shift @_;

    $stopAfter =~ tr/A-Z/a-z/;

    if ((defined($stopAfter)) &&
        (defined(getGlobal('stopAfter'))) &&
        (getGlobal('stopAfter') eq $stopAfter)) {
        print STDERR "Stop requested after '$stopAfter'.\n";
        exit(0);
    }
}



#sub diskSpace ($) {
#    my $wrk   = shift @_;
#
#    my ($bsize, $frsize, $blocks, $bfree, $bavail, $files, $ffree, $favail, $flag, $namemax) = statvfs($wrk);
#
#    print STDERR "bsize  $bsize\n";
#    print STDERR "frsize  $frsize\n";
#    print STDERR "blocks  $blocks\n";
#    print STDERR "bfree  $bfree\n";
#    print STDERR "bavail  $bavail\n";
#    print STDERR "files  $files\n";
#    print STDERR "ffree  $ffree\n";
#    print STDERR "favail  $favail\n";
#    print STDERR "flag  $flag\n";
#    print STDERR "namemax  $namemax\n";
#
#    my $used = $blocks - $bfree;
#
#    #my $total = int($bsize * $blocks / 1048576);
#    #my $used  = int($bsize * $used   / 1048576);
#    #my $free  = int($bsize * $bfree  / 1048576);
#    #my $avail = int($bsize * $bavail / 1048576);
#
#    my $total = $bsize * $blocks;
#    my $used  = $bsize * $used;
#    my $free  = $bsize * $bfree;
#    my $avail = $bsize * $bavail;
#
#    print STDERR "Disk space: total $total GB, used $used GB, free $free GB, available $avail GB\n";
#
#    return (wantarray) ? ($total, $used, $free, $avail) : $avail;
#}


sub diskSpace ($) {
    my $wrk   = shift @_;
    my $df    = df($wrk, 1024);

    my $total = int(10 * $df->{blocks} / 1048576) / 10;
    my $used  = int(10 * $df->{used}   / 1048576) / 10;
    my $free  = int(10 * $df->{bfree}  / 1048576) / 10;
    my $avail = int(10 * $df->{bavail} / 1048576) / 10;

    #print STDERR "Disk space: total $total GB, used $used GB, free $free GB, available $avail GB\n";

    return (wantarray) ? ($total, $used, $free, $avail) : $avail;
}



1;
