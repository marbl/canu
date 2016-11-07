
###############################################################################
 #
 #  This file is part of canu, a software program that assembles whole-genome
 #  sequencing reads into contigs.
 #
 #  This software is based on:
 #    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 #    the 'kmer package' (http://kmer.sourceforge.net)
 #  both originally distributed by Applera Corporation under the GNU General
 #  Public License, version 2.
 #
 #  Canu branched from Celera Assembler at its revision 4587.
 #  Canu branched from the kmer project at its revision 1994.
 #
 #  This file is derived from:
 #
 #    kmer/ESTmapper/scheduler.pm
 #    kmer/scripts/libBri.pm
 #    kmer/scripts/scheduler.pm
 #    src/pipelines/ca3g/Execution.pm
 #
 #  Modifications by:
 #
 #    Brian P. Walenz from 2003-JAN-03 to 2003-NOV-11
 #      are Copyright 2003 Applera Corporation, and
 #      are subject to the GNU General Public License version 2
 #
 #    Brian P. Walenz on 2004-MAR-22
 #      are Copyright 2004 Brian P. Walenz, and
 #      are subject to the GNU General Public License version 2
 #
 #    Brian P. Walenz from 2006-APR-07 to 2011-DEC-28
 #      are Copyright 2006,2008-2009,2011 J. Craig Venter Institute, and
 #      are subject to the GNU General Public License version 2
 #
 #    Brian P. Walenz from 2015-FEB-27 to 2015-SEP-11
 #      are Copyright 2015 Battelle National Biodefense Institute, and
 #      are subject to the BSD 3-Clause License
 #
 #    Brian P. Walenz beginning on 2015-NOV-03
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #    Sergey Koren beginning on 2015-NOV-25
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

package canu::Execution;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(stopBefore stopAfter skipStage emitStage touch getInstallDirectory getJobIDShellCode getLimitShellCode getBinDirectory getBinDirectoryShellCode submitScript submitOrRunParallelJob runCommand runCommandSilently findCommand findExecutable caExit caFailure);

use strict;
use Config;            #  for @signame
use Cwd qw(getcwd);
use Carp qw(cluck);

use POSIX ":sys_wait_h";  #  For waitpid(..., &WNOHANG)
use List::Util qw(min max);
use File::Path qw(make_path remove_tree);
use File::Spec;

use canu::Defaults;


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
    my @running;

    #  Reap any processes that have finished

    foreach my $i (@processesRunning) {
        push @running, $i  if (schedulerReapProcess($i) == 0);
    }

    @processesRunning = @running;

    #  Run processes in any available slots

    while ((scalar(@processesRunning) < $numberOfProcesses) &&
           (scalar(@processQueue) > 0)) {
        my $process = shift @processQueue;
        print STDERR "    $process\n";
        push @processesRunning, schedulerForkProcess($process);
    }
}

sub schedulerFinish ($) {
    my $dir = shift @_;
    my $child;
    my @newProcesses;
    my $remain;

    $remain = scalar(@processQueue);

    my $startsecs = time();
    my $diskfree  = (defined($dir)) ? (diskSpace($dir)) : (0);

    print STDERR "----------------------------------------\n";
    print STDERR "-- Starting concurrent execution on ", scalar(localtime()), " with $diskfree GB free disk space ($remain processes; $numberOfProcesses concurrently)\n"  if  (defined($dir));
    print STDERR "-- Starting concurrent execution on ", scalar(localtime()), " ($remain processes; $numberOfProcesses concurrently)\n"                                    if (!defined($dir));
    print STDERR "\n";

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

    $diskfree = (defined($dir)) ? (diskSpace($dir)) : (0);

    my $warning = "  !!! WARNING !!!" if ($diskfree < 10);
    my $elapsed = time() - $startsecs;

    $elapsed = "lickety-split"    if ($elapsed eq "0");
    $elapsed = "$elapsed second"  if ($elapsed eq "1");
    $elapsed = "$elapsed seconds" if ($elapsed  >  1);

    print STDERR "\n";
    print STDERR "-- Finished on ", scalar(localtime()), " ($elapsed) with $diskfree GB free disk space$warning\n";
    print STDERR "----------------------------------------\n";
}



#
#  File Management
#

sub touch ($@) {
    open(F, "> $_[0]") or caFailure("failed to touch file '$_[0]'", undef);
    print F "$_[1]\n"  if (defined($_[1]));
    close(F);
}



#
#  State management
#

sub stopBefore ($$) {
    my $stopBefore = shift @_;
    my $cmd        = shift @_;

    $stopBefore =~ tr/A-Z/a-z/;

    if ((defined($stopBefore)) &&
        (defined(getGlobal("stopBefore"))) &&
        (getGlobal("stopBefore") eq $stopBefore)) {
        print STDERR "\n";
        print STDERR "Stop requested before '$stopBefore'.\n";
        print STDERR "\n";
        print STDERR "Command:\n  $cmd\n" if (defined($cmd));
        exit(0);
    }
}

sub stopAfter ($) {
    my $stopAfter = shift @_;

    $stopAfter =~ tr/A-Z/a-z/;

    if ((defined($stopAfter)) &&
        (defined(getGlobal("stopAfter"))) &&
        (getGlobal("stopAfter") eq $stopAfter)) {
        print STDERR "Stop requested after '$stopAfter'.\n";
        exit(0);
    }
}


sub emitStage ($$$@) {
    return;
}


sub skipStage ($$$@) {
    return(0);
}



#  Decide what bin directory to use.
#
#  When we are running on the grid, the path of this perl script is NOT always the correct
#  architecture.  If the submission host is FreeBSD, but the grid is Linux, the BSD box will submit
#  FreeBSD/bin/canu to the grid.  Unless it knows which grid host it will run on in advance, there
#  is no way to pick the correct one.  The grid host then has to have enough smarts to choose the
#  correct binaries, and that is what we're doing here.
#
#  To make it more trouble, shell scripts need to do all this by themselves.
#
sub getInstallDirectory () {
    my $installDir = $FindBin::RealBin;

    if ($installDir =~ m!^(.*)/\w+-\w+/bin$!) {
        $installDir = $1;
    }

    return($installDir);
}


#  Emits a block of shell code to parse the grid task id and offset.
#  Expects zero or one argument, which is interpreted different in grid and non-grid mode.
#    Off grid - the job to run
#    On grid  - an offset to add to SGE_TASK_ID or SLURM_ARRAY_TASK_ID to compute the job to run
#
#  PBSPro refuses to run an array job with one element.  They're submitted as a normal job.  Here,
#  we check if it is running on the grid and if the task ID (aka, array ID) isn't set.  If so, we
#  assume it is job 1.
#
sub getJobIDShellCode () {
    my $string;
    my $taskenv = getGlobal('gridEngineTaskID');

    $string .= "#  Discover the job ID to run, from either a grid environment variable and a\n";
    $string .= "#  command line offset, or directly from the command line.\n";
    $string .= "#\n";
    $string .= "if [ x\$PBS_JOBID != x -a x\$$taskenv = x ]; then\n"   if (uc(getGlobal("gridEngine")) eq "PBSPRO");
    $string .= "  $taskenv=1\n"                                        if (uc(getGlobal("gridEngine")) eq "PBSPRO");
    $string .= "fi\n"                                                  if (uc(getGlobal("gridEngine")) eq "PBSPRO");
    $string .= "if [ x\$$taskenv = x -o x\$$taskenv = xundefined -o x\$$taskenv = x0 ]; then\n";
    $string .= "  baseid=\$1\n";           #  Off grid
    $string .= "  offset=0\n";
    $string .= "else\n";
    $string .= "  baseid=\$$taskenv\n";    #  On Grid
    $string .= "  offset=\$1\n";
    $string .= "fi\n";
    $string .= "if [ x\$offset = x ]; then\n";
    $string .= "  offset=0\n";
    $string .= "fi\n";
    $string .= "if [ x\$baseid = x ]; then\n";
    $string .= "  echo Error: I need $taskenv set, or a job index on the command line.\n";
    $string .= "  exit\n";
    $string .= "fi\n";
    $string .= "jobid=`expr \$baseid + \$offset`\n";
    $string .= "if [ x\$$taskenv = x ]; then\n";
    $string .= "  echo Running job \$jobid based on command line options.\n";
    $string .= "else\n";
    $string .= "  echo Running job \$jobid based on $taskenv=\$$taskenv and offset=\$offset.\n";
    $string .= "fi\n";
}


#  Emits a block of shell code to change shell imposed limit on the number of open files and
#  processes.
#
sub getLimitShellCode ($) {
    my $which = shift @_;
    my $string;

    if ($which eq "processes") {
        $string .= "\n";
        $string .= "max=`ulimit -Hu`\n";
        $string .= "bef=`ulimit -Su`\n";
        $string .= "if [ \$bef -lt \$max ] ; then\n";
        $string .= "  ulimit -Su \$max\n";
        $string .= "  aft=`ulimit -Su`\n";
        $string .= "  echo \"Changed max processes per user from \$bef to \$aft (max \$max).\"\n";
        $string .= "  echo \"\"\n";
        $string .= "else\n";
        $string .= "  echo \"Max processes per user limited to \$bef, no increase possible.\"\n";
        $string .= "  echo \"\"\n";
        $string .= "fi\n";
        $string .= "\n";
    }

    if ($which eq "files") {
        $string .= "\n";
        $string .= "max=`ulimit -Hn`\n";
        $string .= "bef=`ulimit -Sn`\n";
        $string .= "if [ \$bef -lt \$max ] ; then\n";
        $string .= "  ulimit -Sn \$max\n";
        $string .= "  aft=`ulimit -Sn`\n";
        $string .= "  echo \"Changed max open files from \$bef to \$aft (max \$max).\"\n";
        $string .= "  echo \"\"\n";
        $string .= "else\n";
        $string .= "  echo \"Max open files limited to \$bef, no increase possible.\"\n";
        $string .= "  echo \"\"\n";
        $string .= "fi\n";
        $string .= "\n";
    }

    return($string);
}


#  Used inside canu to find where binaries are located.  It uses uname to find OS, architecture and
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
        $path = $installDir;
    }

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
    $string .= "\n";

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

    $string .= "if [ ! -d \"\$bin\" ] ; then\n";
    $string .= "  bin=\"$installDir\"\n";
    $string .= "fi\n";
    $string .= "\n";

    return($string);
}





#  Spend too much effort ensuring that the name is unique in the system.  For 'canu' jobs, we don't
#  care.

sub makeRandomSuffix ($) {
    my $length = shift @_;
    my @chars  = +('0'..'9', 'a'..'k', 'm'..'z', 'A'..'H', 'J'..'N', 'P'..'Z');    #  Remove 'l', 'I' and 'O'
    my $suffix;

    while ($length-- > 0) {
        $suffix .= @chars[int(rand(59))];
    }

    return($suffix);
}


sub makeUniqueJobName ($$) {
    my $jobType = shift @_;
    my $asm     = shift @_;

    #  If a canu job, just return the standard name.  No uniquification needed.

    if ($jobType eq "canu") {
        return("canu_" . $asm . ((defined(getGlobal("gridOptionsJobName"))) ? ("_" . getGlobal("gridOptionsJobName")) : ("")));
    }

    #  For all other jobs, we need to ensure the name is unique.  We do this by adding digits at the end.

    my $jobName = "${jobType}_" . $asm . ((defined(getGlobal("gridOptionsJobName"))) ? ("_" . getGlobal("gridOptionsJobName")) : (""));
    my %jobs;

    #  First, find the list of all jobs that exist.

    if (uc(getGlobal("gridEngine")) eq "SGE") {
        open(F, "qstat -xml |");
        while (<F>) {
            $jobs{$1}++  if (m/^\s*<JB_name>(.*)<\/JB_name>$/);
        }
        close(F);
    }

    if (uc(getGlobal("gridEngine")) eq "PBS") {
    }

    if (uc(getGlobal("gridEngine")) eq "PBSPro") {
    }

    if (uc(getGlobal("gridEngine")) eq "LSF") {
    }

    #  If the jobName doesn't exist, we can use it.

    return($jobName)  if (! exists($jobs{$jobName}));

    #  Otherwise, find a unique random 2-letter suffix.

    my $jobIdx  = makeRandomSuffix(2);

    while (exists($jobs{"${jobName}_$jobIdx"})) {
        $jobIdx = makeRandomSuffix(2);
    }

    #  And return it!  Simple!

    # this was breaking dependencies when multiple jobs were submitted like for a failed consensus run, turn off for now
    return("${jobName}");
    #return("${jobName}_$jobIdx");
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

    return   if (getGlobal("useGrid")       ne "1");      #  If not requested to run on the grid,
    return   if (getGlobal("gridEngine")    eq undef);    #  or can't run on the grid, don't run on the grid.

    #  If no job to wait on, and we are already on the grid, do NOT resubmit ourself.
    #
    #  When the user launches canu on the head node, a call to submitScript() is made to launch canu
    #  under grid control.  That results in a restart of canu, and another call to submitScript(),
    #  but this time, the envorinment variable is set, we we can skip the resubmission, and continue
    #  with canu execution.

    return   if (($jobToWaitOn eq undef) && (exists($ENV{getGlobal("gridEngineJobID")})));

    #  Find the next available output file.

    make_path("$wrk/canu-scripts")  if (! -d "$wrk/canu-scripts");  #  Done in canu.pl, just being paranoid

    my $idx = "01";

    while (-e "$wrk/canu-scripts/canu.$idx.out") {
        $idx++;
    }

    my $output    = "$wrk/canu-scripts/canu.$idx.out";
    my $script    = "$wrk/canu-scripts/canu.$idx.sh";
    my $iteration = getGlobal("canuIteration");

    #  Make a script for us to submit.

    open(F, "> $script") or caFailure("failed to open '$script' for writing", undef);
    print F "#!" . getGlobal("shell") . "\n";
    print F "\n"                                                                        if (getGlobal("gridEngine") eq "SGE");
    print F "#  Attempt to (re)configure SGE.  For unknown reasons, jobs submitted\n"   if (getGlobal("gridEngine") eq "SGE");
    print F "#  to SGE, and running under SGE, fail to read the shell init scripts,\n"  if (getGlobal("gridEngine") eq "SGE");
    print F "#  and so they don't set up SGE (or ANY other paths, etc) properly.\n"     if (getGlobal("gridEngine") eq "SGE");
    print F "#  For the record, interactive logins (qlogin) DO set the environment.\n"  if (getGlobal("gridEngine") eq "SGE");
    print F "\n";
    print F "if [ \"x\$SGE_ROOT\" != \"x\" ]; then \n"                                  if (getGlobal("gridEngine") eq "SGE");
    print F "  . \$SGE_ROOT/\$SGE_CELL/common/settings.sh\n"                            if (getGlobal("gridEngine") eq "SGE");
    print F "fi\n"                                                                      if (getGlobal("gridEngine") eq "SGE");
    print F "\n";
    print F "#  On the off chance that there is a pathMap, and the host we\n";
    print F "#  eventually get scheduled on doesn't see other hosts, we decide\n";
    print F "#  at run time where the binary is.\n";
    print F "\n";
    print F getBinDirectoryShellCode();
    print F "\n";
    print F "/usr/bin/env perl \\\n";
    print F "\$bin/canu " . getCommandLineOptions() . " canuIteration=$iteration\n";
    close(F);

    system("chmod +x $script");

    #  Construct a submission command line.

    my ($jobName, $memOption, $thrOption, $gridOpts);

    $jobName   = makeUniqueJobName("canu", $asm);

    #  The canu.pl script isn't expected to take resources.  We'll default to 4gb and one thread.

    my $mem = 4;
    my $thr = 1;

    #  However, the sequential overlap store is still built from within the canu process.

    if (getGlobal("ovsMethod") eq "sequential") {
        $mem = getGlobal("ovsMemory");
        $mem = $2  if ($mem =~ m/^(\d+)-(\d+)$/);
    }

    $memOption = buildMemoryOption($mem, 1);
    $thrOption = buildThreadOption($thr);

    $gridOpts  = $memOption                           if (defined($memOption));
    $gridOpts .= " "                                  if (defined($gridOpts));
    $gridOpts .= $thrOption                           if (defined($thrOption));
    $gridOpts .= " "                                  if (defined($gridOpts));
    $gridOpts .= getGlobal("gridOptions")             if (defined(getGlobal("gridOptions")));
    $gridOpts .= " "                                  if (defined($gridOpts));
    $gridOpts .= getGlobal("gridOptionsExecutive")    if (defined(getGlobal("gridOptionsExecutive")));

    #  If the jobToWaitOn is defined, make the script wait for that to complete.  LSF might need to
    #  query jobs in the queue and figure out the job ID (or IDs) for the jobToWaitOn.  Reading LSF
    #  docs online (for bsub.1) claim that we can still use jobToWaitOn.

    if (defined($jobToWaitOn)) {
        my $hold = getGlobal("gridEngineHoldOption");

        # most grid engines don't understand job names to hold on, only IDs
        if ((uc(getGlobal("gridEngine")) eq "PBS") ||
            (uc(getGlobal("gridEngine")) eq "PBSPRO") ||
            (uc(getGlobal("gridEngine")) eq "SLURM")){
           my $tcmd = getGlobal("gridEngineNameToJobIDCommand");
           $tcmd =~ s/WAIT_TAG/$jobToWaitOn/g;
           my $propJobCount = `$tcmd |wc -l`;
           chomp $propJobCount;
           if ($propJobCount == 0) {
              $tcmd = getGlobal("gridEngineNameToJobIDCommandNoArray");
              $tcmd =~ s/WAIT_TAG/$jobToWaitOn/g;
              $hold = getGlobal("gridEngineHoldOptionNoArray");
              $propJobCount = `$tcmd |wc -l`;
           }
           if ($propJobCount != 1) {
              print STDERR "Warning: multiple IDs for job $jobToWaitOn got $propJobCount and should have been 1.\n";
           }
           my $jobID = undef;
           open(F,  "$tcmd |awk '{print \$1}' |");
           while (<F>) {
              chomp $_;
              if (defined($jobID)) {
                 $jobID = "$jobID:$_";
              } else {
                 $jobID = $_;
              }
           }
           close(F);
           $hold =~ s/WAIT_TAG/$jobID/g;
        } else {
           $hold =~ s/WAIT_TAG/$jobToWaitOn/;
        }
        $gridOpts .= " " . $hold;
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
    my ($name, $bgn, $end, $opt) = @_;
    my  $off = 0;

    #  In some grids (SGE)   this is the maximum size of an array job.
    #  In some grids (Slurm) this is the maximum index of an array job.
    #
    #  So, here, we just don't let any index be above the value.  Both types will be happy.

    if ($end > getGlobal('gridEngineArrayMaxJobs')) {
        $off  = $bgn - 1;
        $bgn -= $off;
        $end -= $off;
    }

    #  PBSPro requires array jobs to have bgn < end.  When $bgn == $end, we
    #  just remove the array qualifier.  But only if this option is setting
    #  the number of jobs, not if it is setting the name.

    if (uc(getGlobal("gridEngine")) eq "PBSPRO") {
        $opt = ""  if (($bgn == $end) && ($opt =~ m/ARRAY_JOBS/));
        $off = "";
    }

    #  Further, PBS/Torque won't let scripts be passed options unless they
    #  are prefixed with a -F....and PBSPro doesn't need this.

    if (uc(getGlobal("gridEngine")) eq "PBS") {
        $off = "-F \"$off\"";
        $off = "";
    }

    $opt =~ s/ARRAY_NAME/$name/g;        #  Replace ARRAY_NAME with 'job name'
    $opt =~ s/ARRAY_JOBS/$bgn-$end/g;    #  Replace ARRAY_JOBS with 'bgn-end'

    return($opt, $off);
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


sub buildMemoryOption ($$) {
    my $m = shift @_;
    my $t = shift @_;
    my $r;
    my $u = "g";

    if (uc(getGlobal("gridEngine")) eq "SGE") {
        $m /= $t;
    }

    if ((uc(getGlobal("gridEngine")) eq "SLURM") && (getGlobal("gridEngineMemoryOption") =~ m/mem-per-cpu/i)) {
        $m /= $t;
    }

    if (int($m) != $m) {
        $m = int($m * 1024);
        $u = "m";
    }

    if (uc(getGlobal("gridEngine")) eq "LSF") {
        $m = $m / 1024          if (getGlobal("gridEngineMemoryUnits") eq "t");
        $m = $m * 1             if (getGlobal("gridEngineMemoryUnits") eq "g");
        $m = $m * 1024          if (getGlobal("gridEngineMemoryUnits") eq "m");
        $m = $m * 1024 * 1024   if (getGlobal("gridEngineMemoryUnits") eq "k");
        $u = "";
    }

    $r =  getGlobal("gridEngineMemoryOption");
    $r =~ s/MEMORY/${m}${u}/g;

    return($r);
}


sub buildThreadOption ($) {
    my $t = shift @_;
    my $r;

    $r =  getGlobal("gridEngineThreadsOption");
    $r =~ s/THREADS/$t/g;

    return($r);
}


sub buildGridJob ($$$$$$$$) {
    my $asm     = shift @_;
    my $jobType = shift @_;
    my $path    = shift @_;
    my $script  = shift @_;
    my $mem     = shift @_;
    my $thr     = shift @_;
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

    my $submitCommand          = getGlobal("gridEngineSubmitCommand");
    my $nameOption             = getGlobal("gridEngineNameOption");

    my $jobNameT               = makeUniqueJobName($jobType, $asm);

    my ($jobName,  $jobOff)    = buildGridArray($jobNameT, $bgnJob, $endJob, getGlobal("gridEngineArrayName"));
    my ($arrayOpt, $arrayOff)  = buildGridArray($jobNameT, $bgnJob, $endJob, getGlobal("gridEngineArrayOption"));

    my $outputOption           = getGlobal("gridEngineOutputOption");
    my $outName                = buildOutputName($path, $script, getGlobal("gridEngineArraySubmitID"));

    my $memOption              = buildMemoryOption($mem, $thr);
    my $thrOption              = buildThreadOption($thr);

    my $gridOpts;

    $gridOpts  = getGlobal("gridOptions")          if (defined(getGlobal("gridOptions")));
    $gridOpts .= " "                               if (defined($gridOpts));
    $gridOpts .= getGlobal("gridOptions$jobType")  if (defined(getGlobal("gridOptions$jobType")));
    $gridOpts .= " "                               if (defined($gridOpts));
    $gridOpts .= $memOption                        if (defined($memOption));
    $gridOpts .= " "                               if (defined($gridOpts));
    $gridOpts .= $thrOption                        if (defined($thrOption));

    #  Build the command line.

    my $cmd;
    $cmd  = "  $submitCommand \\\n";
    $cmd .= "    $gridOpts \\\n"  if (defined($gridOpts));
    $cmd .= "    $nameOption \"$jobName\" \\\n";
    $cmd .= "    $arrayOpt \\\n";
    $cmd .= "    $outputOption $outName \\\n";
    $cmd .= "    $path/$script.sh $arrayOff\n";

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

    #  Expand the ranges into a simple list of job ids.

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

    #  Sort.

    my @jobsA = sort { $a <=> $b } @jobs;

    undef @jobs;

    #  Merge adjacent ids into a range.

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


    #  In some grids (SGE)   this is the maximum size of an array job.
    #  In some grids (Slurm) this is the maximum index of an array job.
    #
    #  So, here, we make blocks that have at most that many jobs.  When we submit the job, we'll
    #  offset the indices to be 1..Max.

    my $l = getGlobal("gridEngineArrayMaxJobs") - 1;

    if ($l > 0) {
        @jobsA = @jobs;
        undef @jobs;

        foreach my $j (@jobsA) {
            if ($j =~ m/^(\d+)-(\d+)$/) {
                my $b = $1;
                my $e = $2;

                while ($b <= $e) {
                    my $B = ($b + $l < $e) ? ($b + $l) : $e;
                    push @jobs, "$b-$B";
                    $b += $l + 1;
                }
            } else {
                push @jobs, $j
            }
        }

        undef @jobsA;
    }

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
sub submitOrRunParallelJob ($$$$$@) {
    my $wrk          = shift @_;  #  Root of the assembly (NOT wrk/correction or wrk/trimming)
    my $asm          = shift @_;  #  Name of the assembly

    my $jobType      = shift @_;  #  E.g., ovl, cns, ... - populates 'gridOptionsXXX
                                  #                      - also becomes the grid job name prefix, so three letters suggested

    my $path         = shift @_;  #  Location of script to run
    my $script       = shift @_;  #  Runs $path/$script.sh > $path/$script.######.out

    my $mem          = getGlobal("${jobType}Memory");
    my $thr          = getGlobal("${jobType}Threads");

    my @jobs         = convertToJobRange(@_);

    #  The script MUST be executable.

    system("chmod +x \"$path/$script.sh\"");

    #  Report what we're doing.

    #my $t = localtime();
    #print STDERR "----------------------------------------GRIDSTART $t\n";
    #print STDERR "$path/$script.sh with $mem gigabytes memory and $thr threads.\n";

    #  Check stopping rules.

    stopBefore($jobType, "$path/$script.sh");

    #  Break infinite loops.  If the grid jobs keep failing, give up after a few attempts.
    #
    #  submitScript() passes canuIteration on to the next call.
    #  canuIteration is reset to zero if the Check() for any parallel step succeeds.
    #
    #  Assuming grid jobs die on each attempt:
    #    0) canu run from the command line submits iteration 1; canuIteration is NOT incremented
    #       because no parallel jobs have been submitted.
    #    1) Iteration 1 - canu.pl submits jobs, increments the interation count, and submits itself as iteration 2
    #    2) Iteration 2 - canu.pl submits jobs, increments the interation count, and submits itself as iteration 3
    #    3) Iteration 3 - canu.pl fails with the error below
    #
    #  If the jobs succeed in Iteration 2, the canu in iteration 3 will pass the Check(), never call
    #  this function, and continue the pipeline.

    caExit("canu iteration count too high, stopping pipeline (most likely a problem in the grid-based computes)", undef)
        if (getGlobal("canuIteration") > getGlobal("canuIterationMax"));

    setGlobal("canuIteration", getGlobal("canuIteration") + 1);


    #  If 'gridEngineJobID' environment variable exists (SGE: JOB_ID; LSF: LSB_JOBID) then we are
    #  currently running under grid crontrol.  If so, run the grid command to submit more jobs, then
    #  submit ourself back to the grid.  If not, tell the user to run the grid command by hand.

    #  Jobs under grid control, and we submit them

    if (defined(getGlobal("gridEngine")) &&
        (getGlobal("useGrid") eq "1") &&
        (getGlobal("useGrid$jobType") eq "1") &&
        (exists($ENV{getGlobal("gridEngineJobID")}))) {
        my $cmd;
        my $jobName;

        foreach my $j (@jobs) {
            ($cmd, $jobName) = buildGridJob($asm, $jobType, $path, $script, $mem, $thr, $j, undef);

            runCommand($path, $cmd) and caFailure("Failed to submit batch jobs", undef);
        }

        submitScript($wrk, $asm, $jobName);

        #  submitScript() should never return.  If it does, then a parallel step was attempted too many time.

        caExit("Too many attempts to run a parallel stage on the grid.  Stop.", undef);
    }

    #  Jobs under grid control, but the user must submit them

    if (defined(getGlobal("gridEngine")) &&
        (getGlobal("useGrid") ne "0") &&
        (getGlobal("useGrid$jobType") eq "1") &&
        (! exists($ENV{getGlobal("gridEngineJobID")}))) {
        print STDERR "\n";
        print STDERR "Please submit the following jobs to the grid for execution using $mem gigabytes memory and $thr threads:\n";
        print STDERR "\n";

        foreach my $j (@jobs) {
            my ($cmd, $jobName) = buildGridJob($asm, $jobType, $path, $script, $mem, $thr, $j, undef);

            print $cmd;
        }

        print STDERR "\n";
        print STDERR "When all jobs complete, restart canu as before.\n";

        exit(0);
    }

    #  Standard jobs, run locally.

    my $cwd = getcwd();  #  Remember where we are.
    chdir($path);        #  So we can root the jobs in the correct location.

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


    # compute limit based on # of cpus
    my $nCParallel  = getGlobal("${jobType}Concurrency");
    $nCParallel     = int(getGlobal("maxThreads") / $thr)  if ((!defined($nCParallel)) || ($nCParallel == 0));
    $nCParallel     = 1                                    if ((!defined($nCParallel)) || ($nCParallel == 0));

    # compute limit based on physical memory
    my $nMParallel = getGlobal("${jobType}Concurrency");
    $nMParallel    = int(getGlobal("maxMemory") / getGlobal("${jobType}Memory")) if ((!defined($nMParallel)) || ($nMParallel == 0));
    $nMParallel    = 1                                                           if ((!defined($nMParallel)) || ($nMParallel == 0));

    # run min of our limits
    my $nParallel  = $nCParallel < $nMParallel ? $nCParallel : $nMParallel;

    schedulerSetNumberOfProcesses($nParallel);
    schedulerFinish($path);

    chdir($cwd);
}




#  Pretty-ify the command.  If there are no newlines already in it, break
#  before every switch and before file redirects.

sub prettifyCommand ($) {
    my $dis = shift @_;

    if (($dis =~ tr/\n/\n/) == 0) {
        $dis =~ s/\s-/ \\\n  -/g;    #  Replace ' -' with '\n  -' (newline, two spaces, then the dash)
        $dis =~ s/\s>\s/ \\\n> /;    #  Replace ' > ' with '\n> '
        $dis =~ s/\s2>\s/ \\\n2> /;  #  Replace ' 2> ' with '\n2> '
    }

    $dis = "    " . $dis;    #  Indent the command by four spaces.
    $dis =~ s/\n/\n    /g;

    return($dis);
}


sub reportRunError ($) {
    my $rc  = shift @_;

    #  Bunch of busy work to get the names of signals.  Is it really worth it?!

    my @signame;
    if (defined($Config{sig_name})) {
        my $i = 0;
        foreach my $n (split('\s+', $Config{sig_name})) {
            $signame[$i] = $n;
            $i++;
        }
    } else {
        for (my $i=0; $i<127; $i++) {
            $signame[$i] = "signal $i";
        }
    }

    #  The rest is rather straightforward at least.

    print STDERR "ERROR:\n";

    if      ($rc ==  -1) {
        print STDERR "ERROR:  Failed to run the command.  (rc=$rc)\n";
    } elsif ($rc  & 127) {
        print STDERR "ERROR:  Failed with signal $signame[$rc & 127].  (rc=$rc)\n";
    } else {
        print STDERR "ERROR:  Failed with exit code ", $rc >> 8 , ".  (rc=$rc)\n";
    }

    print STDERR "ERROR:\n";
}


#  Utility to run a command and check the exit status, report time used.
#
sub runCommand ($$) {
    my $dir = shift @_;
    my $cmd = shift @_;
    my $dis = prettifyCommand($cmd);

    return  if ($cmd eq "");

    #  Check if the directory exists.

    if (! -d $dir) {
        caFailure("Directory '$dir' doesn't exist, can't run command", "");
    }

    #  If only showing the next command, show it and stop.

    if (getGlobal("showNext")) {
        print STDERR "--NEXT-COMMAND\n";
        print STDERR "$dis\n";
        exit(0);
    }

    #  Log that we're starting, and show the pretty-ified command.

    my $startsecs = time();
    my $diskfree  = (defined($dir)) ? (diskSpace($dir)) : (0);

    print STDERR "----------------------------------------\n";
    print STDERR "-- Starting command on ", scalar(localtime()), " with $diskfree GB free disk space\n"  if  (defined($dir));
    print STDERR "-- Starting command on ", scalar(localtime()), "\n"                                    if (!defined($dir));
    print STDERR "\n";
    print STDERR "$dis\n";

    my $rc = 0xffff & system("cd $dir && $cmd");

    $diskfree = (defined($dir)) ? (diskSpace($dir)) : (0);

    my $warning = "  !!! WARNING !!!" if ($diskfree < 10);
    my $elapsed = time() - $startsecs;

    $elapsed = "lickety-split"    if ($elapsed eq "0");
    $elapsed = "$elapsed second"  if ($elapsed eq "1");
    $elapsed = "$elapsed seconds" if ($elapsed  >  1);

    print STDERR "\n";
    print STDERR "-- Finished on ", scalar(localtime()), " ($elapsed) with $diskfree GB free disk space$warning\n";
    print STDERR "----------------------------------------\n";

    #  Pretty much copied from Programming Perl page 230

    return(0) if ($rc == 0);

    reportRunError($rc);

    return(1);
}



sub runCommandSilently ($$$) {
    my $dir      = shift @_;
    my $cmd      = shift @_;
    my $dis      = prettifyCommand($cmd);
    my $critical = shift @_;

    return(0)   if ($cmd eq "");

    if (! -d $dir) {
        caFailure("Directory '$dir' doesn't exist, can't run command", "");
    }

    my $rc = 0xffff & system("cd $dir && $cmd");

    return(0) if ($rc == 0);         #  No errors, return no error.
    return(1) if ($critical == 0);   #  If not critical, return that it failed, otherwise, report error and fail.

    print STDERR "$dis\n";

    reportRunError($rc);

    return(1);
}



sub findCommand ($) {
    my $cmd  = shift @_;
    my @path = File::Spec->path;

    for my $path (@path) {
        if (-x "$path/$cmd") {
            return("$path/$cmd");
        }
    }

    return(undef);
}



sub findExecutable ($) {
    my $exec = shift @_;

    my $path = `which \"$exec\" 2> /dev/null`;

    $path =~ s/^\s+//;
    $path =~ s/\s+$//;

    return(undef)  if ($path eq "");
    return($path);
}


#  Use caExit() for transient errors, like not opening files, processes that die, etc.
sub caExit ($$) {
    my  $wrk   = getGlobal("onExitDir");
    my  $asm   = getGlobal("onExitNam");
    my  $msg   = shift @_;
    my  $log   = shift @_;

    print STDERR "================================================================================\n";
    print STDERR "Don't panic, but a mostly harmless error occurred and canu failed.\n";
    print STDERR "\n";

    #  Really should pass in $wrk
    if (defined($log)) {
        my  $df = diskSpace($log);

        print STDERR "Disk space available:  $df GB\n";
        print STDERR "\n";
    }

    if (-e $log) {
        print STDERR "Last 50 lines of the relevant log file ($log):\n";
        print STDERR "\n";
        system("tail -n 50 $log");
        print STDERR "\n";
    }

    print STDERR "canu failed with '$msg'.\n";
    print STDERR "\n";

    my $fail = getGlobal('onFailure');
    if (defined($fail)) {
        runCommandSilently($wrk, "$fail $asm", 0);
    }

    exit(1);
}


#  Use caFailure() for errors that definitely will require code changes to fix.
sub caFailure ($$) {
    my  $wrk   = getGlobal("onExitDir");
    my  $asm   = getGlobal("onExitNam");
    my  $msg   = shift @_;
    my  $log   = shift @_;

    print STDERR "================================================================================\n";
    print STDERR "Please panic.  canu failed, and it shouldn't have.\n";
    print STDERR "\n";
    print STDERR "Stack trace:\n";
    print STDERR "\n";
    cluck;
    print STDERR "\n";

    if (-e $log) {
        print STDERR "Last few lines of the relevant log file ($log):\n";
        print STDERR "\n";
        system("tail -n 50 $log");
    }

    print STDERR "\n";
    print STDERR "canu failed with '$msg'.\n";

    my $fail = getGlobal('onFailure');
    if (defined($fail)) {
        runCommandSilently($wrk, "$fail $asm", 0);
    }

    exit(1);
}


1;
