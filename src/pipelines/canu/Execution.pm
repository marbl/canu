
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

package canu::Execution;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(stopAfter
             resetIteration
             touch
             makeExecutable
             getJobIDShellCode
             getLimitShellCode
             getBinDirectory
             getBinDirectoryShellCode
             setWorkDirectory
             setWorkDirectoryShellCode
             submitScript
             submitOrRunParallelJob
             runCommand
             runCommandSilently
             findCommand
             findExecutable
             caExit
             caFailure);

use strict;
use warnings "all";
no  warnings "uninitialized";

use Config;            #  for @signame
use Cwd qw(getcwd);
use Carp qw(longmess);

use POSIX ":sys_wait_h";  #  For waitpid(..., &WNOHANG)
use File::Basename;
use List::Util qw(min max);
use File::Path 2.08 qw(make_path remove_tree);
use File::Spec;

use canu::Defaults;





#  Log that we've finished a task.

sub logFinished ($$) {
    my $dir       = shift @_;
    my $startsecs = shift @_;

    my $diskfree  = diskSpace(".");

    my $warning = "  !!! WARNING !!!" if ($diskfree < 10);
    my $elapsed = time() - $startsecs;
    my $message;

    my @fast;

    push @fast, "lickety-split";
    push @fast, "fast as lightning";
    push @fast, "furiously fast";
    push @fast, "like a bat out of hell";
    push @fast, "in the blink of an eye";

    my @slow;

    push @slow, "fashionably late";
    push @slow, "better late than never";
    push @slow, "like watching paint dry";
    push @slow, "at least I didn't crash";
    push @slow, "it'll be worth it in the end";
    push @slow, "no bitcoins found either";

    my $rf = int(rand(scalar(@fast)));
    my $rs = int(rand(scalar(@slow)));
    my $rp = int(rand(100));

    $message  = "$elapsed seconds" if ($elapsed  > 1);
    $message  = "one second"       if ($elapsed == 1);
    $message  = $fast[$rf]         if ($elapsed  < 1);

    $message .= ", " . $slow[$rs]  if ((($elapsed > 1000)  && ($rp < 1)) ||
                                       (($elapsed > 10000) && ($rp < 50)) ||
                                       (($elapsed > 86400)));

    print STDERR "\n";
    print STDERR "-- Finished on ", scalar(localtime()), " ($message) with $diskfree GB free disk space$warning\n";
    print STDERR "----------------------------------------\n";
}



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

sub schedulerFinish ($$) {
    my $dir = shift @_;
    my $nam = shift @_;
    my $child;
    my @newProcesses;
    my $remain;

    $remain = scalar(@processQueue);

    my $startsecs = time();
    my $diskfree  = (defined($dir)) ? (diskSpace($dir)) : (0);

    print STDERR "----------------------------------------\n";
    print STDERR "-- Starting '$nam' concurrent execution on ", scalar(localtime()), " with $diskfree GB free disk space ($remain processes; $numberOfProcesses concurrently)\n"  if  (defined($dir));
    print STDERR "-- Starting '$nam' concurrent execution on ", scalar(localtime()), " ($remain processes; $numberOfProcesses concurrently)\n"                                    if (!defined($dir));
    print STDERR "\n";
    print STDERR "    cd $dir\n";

    my $cwd = getcwd();  #  Remember where we are.
    chdir($dir);        #  So we can root the jobs in the correct location.

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

    logFinished($dir, $startsecs);

    chdir($cwd);
}



#
#  File Management
#

sub touch ($@) {
    open(F, "> $_[0]") or caFailure("failed to touch file '$_[0]'", undef);
    print F "$_[1]\n"  if (defined($_[1]));
    close(F);
}



sub makeExecutable ($) {
    my $file = shift @_;

    chmod(0755 & ~umask(), $file);
}


#
#  State management
#

sub stopAfter ($) {
    my $stopAfter = shift @_;

    $stopAfter =~ tr/A-Z/a-z/;

    return   if (($stopAfter ne "theend") &&
                 ($stopAfter ne getGlobal("stopAfter")));

    if ($stopAfter ne "theend") {
        print STDERR "--\n";
        print STDERR "--  Stop requested after '$stopAfter'.\n";
    }

    if (defined(getGlobal("onSuccess"))) {
        print STDERR "--\n";
        print STDERR "--  Running user-supplied termination command.\n";

        runCommand(getGlobal("onExitDir"), getGlobal("onSuccess") . " " . getGlobal("onExitNam"));
    }

    print STDERR "--\n";
    print STDERR "-- Bye.\n";

    exit(0);
}


sub resetIteration ($) {
    my $stage = shift @_;

    print STDERR "-- Finished stage '$stage', reset canuIteration.\n"       if (defined($stage));

    setGlobal("canuIteration", 0);
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
#sub getInstallDirectory () {
#    my $installDir = ;
#
#    if ($installDir =~ m!^(.*)/\w+-\w+/bin$!) {
#        $installDir = $1;
#    }
#
#    return($installDir);
#}


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
    $string .= "jobid=`expr -- \$baseid + \$offset`\n";
    $string .= "if [ x\$baseid = x0 ]; then\n";
    $string .= "  echo Error: jobid 0 is invalid\\; I need $taskenv set, or a job index on the command line.\n";
    $string .= "  exit\n";
    $string .= "fi\n";
    $string .= "if [ x\$$taskenv = x ]; then\n";
    $string .= "  echo Running job \$jobid based on command line options.\n";
    $string .= "else\n";
    $string .= "  echo Running job \$jobid based on $taskenv=\$$taskenv and offset=\$offset.\n";
    $string .= "fi\n";
}


#  Emits a block of shell code to change shell imposed limit on the number of open files and
#  processes.
#
sub getLimitShellCode () {
    my $string;

    $string .= "echo \"\"\n";
    $string .= "echo \"Attempting to increase maximum allowed processes and open files.\"";
    $string .= "\n";
    $string .= "max=`ulimit -Hu`\n";
    $string .= "bef=`ulimit -Su`\n";
    $string .= "if [ \$bef -lt \$max ] ; then\n";
    $string .= "  ulimit -Su \$max\n";
    $string .= "  aft=`ulimit -Su`\n";
    $string .= "  echo \"  Changed max processes per user from \$bef to \$aft (max \$max).\"\n";
    $string .= "else\n";
    $string .= "  echo \"  Max processes per user limited to \$bef, no increase possible.\"\n";
    $string .= "fi\n";
    $string .= "\n";
    $string .= "max=`ulimit -Hn`\n";
    $string .= "bef=`ulimit -Sn`\n";
    $string .= "if [ \$bef -lt \$max ] ; then\n";
    $string .= "  ulimit -Sn \$max\n";
    $string .= "  aft=`ulimit -Sn`\n";
    $string .= "  echo \"  Changed max open files from \$bef to \$aft (max \$max).\"\n";
    $string .= "else\n";
    $string .= "  echo \"  Max open files limited to \$bef, no increase possible.\"\n";
    $string .= "fi\n";
    $string .= "\n";
    $string .= "echo \"\"\n";
    $string .= "\n";

    return($string);
}


#  Used inside canu to find where binaries are located.
#
sub getBinDirectory () {
    return($FindBin::RealBin);

    #my $idir = getInstallDirectory();
    #my $path = $idir;
    #
    #$path = "$idir/bin"   if (-d "$idir/bin");
    #
    #return($path);
}


#  Emits a block of shell code to locate binaries during shell scripts.  See comments on
#  getBinDirectory.
#
sub getBinDirectoryShellCode () {
    my $idir = $FindBin::RealBin;
    my $string;

    #  First, run any preExec command that might exist.

    if (defined(getGlobal("preExec"))) {
        $string .= "#  Pre-execution commands.\n";
        $string .= "\n";
        $string .= getGlobal('preExec') . "\n";
        $string .= "\n";
    }

    #  Then, setup and report paths.

    my $javaPath = getGlobal("java");
    my $canu     = "\$bin/" . basename($0);   #  NOTE: $bin decided at script run time

    $string .= "\n";
    $string .= "#  Path to Canu.\n";
    $string .= "\n";
    $string .= "bin=\"$idir\"\n";
    $string .= "\n";
    $string .= "#  Report paths.\n";
    $string .= "\n";
    $string .= "echo \"\"\n";
    $string .= "echo \"Found perl:\"\n";
    $string .= "echo \"  \" `which perl`\n";
    $string .= "echo \"  \" `perl --version | grep version`\n";
    $string .= "echo \"\"\n";
    $string .= "echo \"Found java:\"\n";
    $string .= "echo \"  \" `which $javaPath`\n";
    $string .= "echo \"  \" `$javaPath -showversion 2>&1 | head -n 1`\n";
    $string .= "echo \"\"\n";
    $string .= "echo \"Found canu:\"\n";
    $string .= "echo \"  \" $canu\n";
    $string .= "echo \"  \" `$canu -version`\n";
    $string .= "echo \"\"\n";
    $string .= "\n";
    $string .= "\n";
    $string .= "#  Environment for any object storage.\n";
    $string .= "\n";
    $string .= "export CANU_OBJECT_STORE_CLIENT="    . getGlobal("objectStoreClient")    . "\n";
    $string .= "export CANU_OBJECT_STORE_CLIENT_UA=" . getGlobal("objectStoreClientUA")  . "\n";
    $string .= "export CANU_OBJECT_STORE_CLIENT_DA=" . getGlobal("objectStoreClientDA")  . "\n";
    $string .= "export CANU_OBJECT_STORE_NAMESPACE=" . getGlobal("objectStoreNameSpace") . "\n";
    $string .= "export CANU_OBJECT_STORE_PROJECT="   . getGlobal("objectStoreProject")   . "\n";
    $string .= "\n";
    $string .= "\n";

    return($string);
}




#
#  If running on a cloud system, shell scripts are started in some random location.
#  setWorkDirectory() will create the directory the script is supposed to run in (e.g.,
#  correction/0-mercounts) and move into it.  This will keep the scripts compatible with the way
#  they are run from within canu.pl.
#
#  If you're fine running in 'some random location' do nothing here.
#
#  Note that canu does minimal cleanup.
#

sub setWorkDirectory ($$) {
    my $asm     = shift @_;
    my $rootdir = shift @_;

    #  Set the initial directory based on various rules.
    #
    #  For the canu executive, in grid mode, both setWorkDirectoryShellCode and
    #  this (in that order) are called.  TEST is assuming that all (non-executive)
    #  compute jobs are run as arrays.

    if    ((getGlobal("objectStore") eq "TEST") && (defined($ENV{"JOB_ID"}))) {
        my $jid = $ENV{'JOB_ID'};
        my $tid = $ENV{'SGE_TASK_ID'};   #  'undefined' since this isn't an array job.

        remove_tree("/assembly/objectstore/job-$jid");
        make_path  ("/assembly/objectstore/job-$jid");
        chdir      ("/assembly/objectstore/job-$jid");
    }

    elsif (getGlobal("objectStore") eq "DNANEXUS") {
    }

    elsif (getGlobal("gridEngine") eq "PBSPRO") {
        chdir($ENV{"PBS_O_WORKDIR"})   if (exists($ENV{"PBS_O_WORKDIR"}));
        delete $ENV{"PBS_O_WORKDIR"};
    }

    #  Now move into the assembly directory.

    if (defined($rootdir)) {
        make_path($rootdir)  if (! -d $rootdir);
        chdir($rootdir);
    }

    #  And save some pieces we need when we quit.

    setGlobal("onExitDir", getcwd());
    setGlobal("onExitNam", $asm);
}



sub setWorkDirectoryShellCode ($) {
    my $path = shift @_;
    my $code = "";

    if    (getGlobal("objectStore") eq "TEST") {
        $code .= "if [ z\$SGE_TASK_ID != z ] ; then\n";
        $code .= "  jid=\$JOB_ID\n";
        $code .= "  tid=\$SGE_TASK_ID\n";
        $code .= "  if [ x\$tid != xundefined ] ; then\n";
        $code .= "    rm   -rf /assembly/objectstore/job-\$jid-\$tid/\n";
        $code .= "    mkdir -p /assembly/objectstore/job-\$jid-\$tid/$path\n";
        $code .= "    cd       /assembly/objectstore/job-\$jid-\$tid/$path\n";
        $code .= "  fi\n";
        $code .= "fi\n";
    }

    elsif (getGlobal("objectStore") eq "DNANEXUS") {
        #  You're probably fine running in some random location, but if there is faster disk
        #  available, move there.
    }

    elsif (getGlobal("gridEngine") eq "PBSPRO") {
        $code .= "if [ z\$PBS_O_WORKDIR != z ] ; then\n";
        $code .= "  cd \$PBS_O_WORKDIR\n";
        $code .= "fi\n";
    }

    return($code);
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

    if (uc(getGlobal("gridEngine")) eq "DNANEXUS") {
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
sub submitScript ($$) {
    my $asm         = shift @_;
    my $jobHold     = shift @_;

    return   if (getGlobal("useGrid")       ne "1");      #  If not requested to run on the grid,
    return   if (getGlobal("gridEngine")    eq undef);    #  or can't run on the grid, don't run on the grid.

    #  If no job hold, and we are already on the grid, do NOT resubmit ourself.
    #
    #  When the user launches canu on the head node, a call to submitScript() is made to launch canu
    #  under grid control.  That results in a restart of canu, and another call to submitScript(),
    #  but this time, the envorinment variable is set, we we can skip the resubmission, and continue
    #  with canu execution.

    return   if (($jobHold eq undef) && (exists($ENV{getGlobal("gridEngineJobID")})));

    #  Figure out the name of the script we want to be making, and a place for it to write output.

    my $idx = "01";

    while ((-e "canu-scripts/canu.$idx.sh") ||
           (-e "canu-scripts/canu.$idx.out")) {
        $idx++;
    }

    my $script    = "canu-scripts/canu.$idx.sh";
    my $scriptOut = "canu-scripts/canu.$idx.out";

    #  Make a script for us to submit.

    open(F, "> $script") or caFailure("failed to open '$script' for writing", undef);
    print F "#!" . getGlobal("shell") . "\n";
    print F "\n";
    print F "#  Attempt to (re)configure SGE.  For unknown reasons, jobs submitted\n"   if (getGlobal("gridEngine") eq "SGE");
    print F "#  to SGE, and running under SGE, fail to read the shell init scripts,\n"  if (getGlobal("gridEngine") eq "SGE");
    print F "#  and so they don't set up SGE (or ANY other paths, etc) properly.\n"     if (getGlobal("gridEngine") eq "SGE");
    print F "#  For the record, interactive logins (qlogin) DO set the environment.\n"  if (getGlobal("gridEngine") eq "SGE");
    print F "\n"                                                                        if (getGlobal("gridEngine") eq "SGE");
    print F "if [ \"x\$SGE_ROOT\" != \"x\" -a \\\n"                                     if (getGlobal("gridEngine") eq "SGE");
    print F "     -e  \$SGE_ROOT/\$SGE_CELL/common/settings.sh ]; then\n"               if (getGlobal("gridEngine") eq "SGE");
    print F "  . \$SGE_ROOT/\$SGE_CELL/common/settings.sh\n"                            if (getGlobal("gridEngine") eq "SGE");
    print F "fi\n"                                                                      if (getGlobal("gridEngine") eq "SGE");
    print F "\n";
    print F getBinDirectoryShellCode();
    print F "\n";
    print F setWorkDirectoryShellCode(".");
    print F "\n";
    print F "rm -f canu.out\n";
    print F "ln -s $scriptOut canu.out\n";
    print F "\n";
    print F "/usr/bin/env perl \\\n";
    print F "\$bin/" . basename($0) . " " . getCommandLineOptions() . " canuIteration=" . getGlobal("canuIteration") . "\n";
    close(F);

    makeExecutable("$script");

    #  Construct a submission command line.

    my $jobName   = makeUniqueJobName("canu", $asm);

    #  The canu.pl script isn't expected to take resources.  We'll default to 4gb and one thread.

    my $mem = getGlobal("executiveMemory");
    my $thr = getGlobal("executiveThreads");

    my $resOption = buildResourceOption($mem, $thr);

    my $gridOpts;

    $gridOpts  = $jobHold;
    $gridOpts .= " "                                     if (defined($gridOpts));

    #  LSF ignores all but the first option, so options need to be reversed.
    #  DNAnexus doesn't use threads, and memory is the instance type.

    if    (uc(getGlobal("gridEngine")) eq "LSF") {
        $gridOpts .= getGlobal("gridOptionsExecutive")   if (defined(getGlobal("gridOptionsExecutive")));
        $gridOpts .= " "                                 if (defined($gridOpts));
        $gridOpts .= getGlobal("gridOptions")            if (defined(getGlobal("gridOptions")));
        $gridOpts .= " "                                 if (defined($gridOpts));
        $gridOpts .= $resOption                          if (defined($resOption));
    }

    elsif (uc(getGlobal("gridEngine")) eq "DNANEXUS") {
        $gridOpts .= getGlobal("gridOptions")            if (defined(getGlobal("gridOptions")));
        $gridOpts .= " "                                 if (defined($gridOpts));
        $gridOpts .= getGlobal("gridOptionsExecutive")   if (defined(getGlobal("gridOptionsExecutive")));
    }

    else {
        $gridOpts .= $resOption                          if (defined($resOption));
        $gridOpts .= " "                                 if (defined($gridOpts));
        $gridOpts .= getGlobal("gridOptions")            if (defined(getGlobal("gridOptions")));
        $gridOpts .= " "                                 if (defined($gridOpts));
        $gridOpts .= getGlobal("gridOptionsExecutive")   if (defined(getGlobal("gridOptionsExecutive")));
    }

    my $submitCommand        = getGlobal("gridEngineSubmitCommand");
    my $nameOption           = getGlobal("gridEngineNameOption");
    my $outputOption         = getGlobal("gridEngineOutputOption");

    my $qcmd = "$submitCommand $gridOpts";
    $qcmd   .= " $nameOption '$jobName'"   if defined($nameOption);
    $qcmd   .= " $outputOption $scriptOut" if defined($outputOption);

    #  DNAnexus doesn't submit scripts; all parameters are passed through
    #  '-i' options.  The 'fetch_and_run' magic is in dx-canu/src/canu-job-launcher.sh.
    #  It will download the requested shell script and execute said function
    #  in it.

    if (uc(getGlobal("gridEngine")) eq "DNANEXUS") {
        $qcmd .= " \\\n";
        $qcmd .= " -ioutput_folder:string=\"" . getGlobal("objectStoreNamespace") . "\" \\\n";
        $qcmd .= " -iscript_path:string=\"\" \\\n";
        $qcmd .= " -iscript_name:string=\"canu-executive.sh\" \\\n";
        $qcmd .= " -icanu_iteration:int="     . getGlobal("canuIteration")        .   " \\\n";
        $qcmd .= " -icanu_iteration_max:int=" . getGlobal("canuIterationMax")     .   " \\\n";
        $qcmd .= " fetch_and_run \\\n";
    }

    else {
        $qcmd .= "  $script";
    }

    if (runCommand(getcwd(), $qcmd) == 0) {      #  Exit sucessfully if we've submitted
        exit(0);                                 #  the next part successfully.
    }

    print STDERR "-- Failed to submit Canu executive.  Delay 10 seconds and try again.\n";

    sleep(10);

    if (runCommand(getcwd(), $qcmd) == 0) {
        exit(0);
    }

    print STDERR "-- Failed to submit Canu executive.  Giving up after two tries.\n";

    exit(1);
}



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
    #  New versions of PBS have this behavior too

    if (uc(getGlobal("gridEngine")) eq "PBSPRO" || uc(getGlobal("gridEngine")) eq "PBS") {
        if (($bgn == $end) && ($opt =~ m/ARRAY_JOBS/)) {
            $opt = "";
            $off = $bgn;
        }
    }
    # DNA nexus doesn't have arrays and only supports 1 job, which we use to pass the identifier
    # Set the offset to blank since it is not supported as well
    if (uc(getGlobal("gridEngine")) eq "DNANEXUS" && ($bgn == $end) && ($opt =~ m/ARRAY_JOBS/)) {
           my $jid = $bgn + $off;
           $opt =~ s/ARRAY_JOBS/$jid/g;
           $off = "";
    }

    #  Further, PBS/Torque won't let scripts be passed options unless they
    #  are prefixed with a -F....and PBSPro doesn't need this.

    if (uc(getGlobal("gridEngine")) eq "PBS") {
        $off = "-F \"$off\"";
    }

    $opt =~ s/ARRAY_NAME/$name/g;        #  Replace ARRAY_NAME with 'job name'
    $opt =~ s/ARRAY_JOBS/$bgn-$end/g;    #  Replace ARRAY_JOBS with 'bgn-end'

    return($opt, $off);
}


sub buildOutputName ($$$) {
    my $path   = shift @_;
    my $script = shift @_;
    my $tid    = substr("000000" . (shift @_), -6);
    my $o;

    #  When this function is called, canu.pl is running in the assembly directory.
    #  But, when the script is executed, it is rooted in '$path'.  To get the
    #  'logs' working, we need to check if the directory relative to the assembly root exists,
    #  but set it relative to $path (which is also where $script is relative to).

    $o = "$script.$tid.out";
    $o = "logs/$1.$tid.out"   if ((-e "$path/logs") && ($script =~ m/scripts\/(.*)/));

    return($o);
}


sub buildOutputOption ($$) {
    my $path   = shift @_;
    my $script = shift @_;
    my $tid    = getGlobal("gridEngineArraySubmitID");
    my $opt    = getGlobal("gridEngineOutputOption");

    if (defined($tid) && defined($opt)) {
        my $o;

        $o = "$script.$tid.out";
        $o = "logs/$1.$tid.out"   if ((-e "$path/logs") && ($script =~ m/scripts\/(.*)/));

        return("$opt $o");
    }

    return(undef);
}


sub buildStageOption ($$) {
    my $t = shift @_;
    my $d = shift @_;
    my $r;

    if ($t eq "cor" || $t eq "cormhap" || $t eq "obtmhap" || $t eq "utgmhap") {
        $r =  getGlobal("gridEngineStageOption");
        $r =~ s/DISK_SPACE/${d}/g;
    }

    return($r);
}


sub buildResourceOption ($$) {
    my $m = shift @_;
    my $t = shift @_;
    my $u = "g";

    #  Increase memory slightly if this is a retry.

    if (getGlobal("canuIteration") > 0) {
        $m *= 1.25 ** (getGlobal("canuIteration")-1);
    }

    #  Massage the memory requested into a format the grid is happy with.

    if (getGlobal("gridEngineMemoryPerJob") != "1") {    #  If anything but "1", divide the memory request
        $m /= $t;                                        #  by the number of slots we request.  Default behavior
    }                                                    #  for SGE and Slurm when mem-per-cpu is used.

    if (uc(getGlobal("gridEngine")) eq "LSF") {                                     #  But then reset for LSF,
        # always round up
        $m = int($m / 1024 + 0.5)          if (getGlobal("gridEngineMemoryUnits") =~ m/t/i);   #  because LSF wants to
        $m = int($m * 1 + 0.5)             if (getGlobal("gridEngineMemoryUnits") =~ m/g/i);   #  enforce the units used.
        $m = int($m * 1024 + 0.5)          if (getGlobal("gridEngineMemoryUnits") =~ m/m/i);
        $m = int($m * 1024 * 1024 + 0.5)   if (getGlobal("gridEngineMemoryUnits") =~ m/k/i);
        $u = "";
    }

    if (uc(getGlobal("gridEngine")) eq "DNANEXUS") {
       $m = canu::Grid_DNANexus::getDNANexusInstance($m, $t);
       $u = "";
    }

    if (($u eq "g") &&          #  If we're not an integral number of gigabytes,
        (int($m) != $m)) {      #  switch over to megabytes.
        $m = int($m * 1024);    #    But only if we're still gigabytes!
        $u = "m";               #    In particular, both LSF and DNANEXUS set units to "".
    }

    #  Replace MEMORY and THREADS with actual values.

    my $r = getGlobal("gridEngineResourceOption");

    $r =~ s/MEMORY/${m}${u}/g;
    $r =~ s/THREADS/$t/g;

    return($r);
}


sub purgeGridJobSubmitScripts ($$) {
    my $path    = shift @_;
    my $script  = shift @_;
    my $idx     = "01";

    while (-e "$path/$script.jobSubmit-$idx.sh") {
        unlink "$path/$script.jobSubmit-$idx.sh";
        $idx++;
    }
}


sub buildGridJob ($$$$$$$$$) {
    my $asm     = shift @_;
    my $jobType = shift @_;
    my $path    = shift @_;
    my $script  = shift @_;
    my $mem     = shift @_;
    my $thr     = shift @_;
    my $dsk     = shift @_;
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

    my $outputOption           = buildOutputOption($path, $script);

    my $stageOption            = buildStageOption($jobType, $dsk);
    my $resOption              = buildResourceOption($mem, $thr);
    my $globalOptions          = getGlobal("gridOptions");
    my $jobOptions             = getGlobal("gridOptions$jobType");

    my $opts;

    $opts  = "$stageOption "    if (defined($stageOption));
    $opts .= "$resOption "      if (defined($resOption));
    $opts .= "$globalOptions "  if (defined($globalOptions));
    $opts .= "$jobOptions "     if (defined($jobOptions));
    $opts .= "$outputOption "   if (defined($outputOption));
    $opts =~ s/\s+$//;

    #  Find a unique file name to save the command.

    my $idx = "01";

    while (-e "$path/$script.jobSubmit-$idx.sh") {
        $idx++;
    }

    #  Build and save the command line.  Return the command PREFIX (we'll be adding .sh and .out as
    #  appropriate), and the job name it will be submitted with (which isn't expected to be used).

    open(F, "> $path/$script.jobSubmit-$idx.sh") or die;
    print F "#!/bin/sh\n";
    print F "\n";
    print F "$submitCommand \\\n";
    print F "  $opts \\\n"  if (defined($opts));
    print F "  $nameOption \"$jobName\" \\\n";
    print F "  $arrayOpt \\\n";


    if (uc(getGlobal("gridEngine")) eq "PBSPRO") {    #  PBSpro needs '--' to tell it to
        print F " -- ";                               #  stop parsing the command line.
    }

    #  DNAnexus wants job parameters via options to the submit command;
    #  everyone else wants the script itself.

    if (uc(getGlobal("gridEngine")) eq "DNANEXUS") {
        print F "  -ioutput_folder:string=\"" . getGlobal("objectStoreNamespace") . "\" \\\n";
        print F "  -iscript_path:string=\"$path\" \\\n";
        print F "  -iscript_name:string=\"$script.sh\" \\\n";
        print F "  fetch_and_run \\\n";
    } else {
        print F "  `pwd`/$script.sh $arrayOff \\\n";
    }

    print F "> ./$script.jobSubmit-$idx.out 2>&1\n";
    close(F);

    makeExecutable("$path/$script.jobSubmit-$idx.sh");

    return("$script.jobSubmit-$idx", $jobName);
}




#  Convert @jobs to a list of ranges, 1-4, 5, 10-20, etc.  These will be directly submitted to the
#  grid, or run one-by-one locally.
#
#  If we're SGE, we can combine everything to one job range: 1-4, 5, 10-20.  Except that
#  buildGridJob() doesn't know how to handle that.

sub convertToJobRange (@) {
    my @jobs;

    #  Expand the ranges into a simple list of job ids.

    foreach my $j (@_) {
        if        ($j =~ m/^0*(\d+)-0*(\d+)$/) {
            for (my $a=$1; $a<=$2; $a++) {
                push @jobs, $a;
            }

        } elsif ($j =~ m/^0*(\d+)$/) {
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

    if ($l >= 0) {
        @jobsA = @jobs;
        undef @jobs;

        foreach my $j (@jobsA) {
            if ($j =~ m/^0*(\d+)-0*(\d+)$/) {
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



sub countJobsInRange (@) {
    my @jobs  = @_;
    my $nJobs = 0;

    foreach my $j (@jobs) {
        if ($j =~ m/^(\d+)-(\d+)$/) {
            $nJobs += $2 - $1 + 1;
        } else {
            $nJobs++;
        }
    }

    return($nJobs);
}



#  Expects
#    job type ("ovl", etc)
#    output directory
#    script name with no directory or .sh
#    number of jobs in the task
#
#  If under grid control, submit grid jobs.  Otherwise, run in parallel locally.
#
sub submitOrRunParallelJob ($$$$@) {
    my $asm          = shift @_;  #  Name of the assembly

    my $jobType      = shift @_;  #  E.g., ovl, cns, ... - populates 'gridOptionsXXX
                                  #                      - also becomes the grid job name prefix, so three letters suggested

    my $path         = shift @_;  #  Location of script to run
    my $script       = shift @_;  #  Runs $path/$script.sh > $path/$script.######.out

    my $mem          = getGlobal("${jobType}Memory");
    my $thr          = getGlobal("${jobType}Threads");
    my $dsk          = getGlobal("${jobType}StageSpace");

    my @jobs         = convertToJobRange(@_);
    my $nJobs        = countJobsInRange(@jobs);

    my $runDirectly  = 0;

    #  The script MUST be executable.

    makeExecutable("$path/$script.sh");

    #  If the job can fit in the task running the executive, run it right here.

    if (($nJobs * $mem + 0.5 <= getGlobal("executiveMemory")) &&
        ($nJobs * $thr       <= getGlobal("executiveThreads"))) {
        $runDirectly = 1;
    }

    #  Report what we're doing.

    #my $t = localtime();
    #print STDERR "----------------------------------------GRIDSTART $t\n";
    #print STDERR "$path/$script.sh with $mem gigabytes memory and $thr threads.\n";

    #  Break infinite loops.  If the grid jobs keep failing, give up after a few attempts.
    #
    #  submitScript() passes canuIteration on to the next call.
    #  canuIteration is reset to zero if the Check() for any parallel step succeeds.
    #
    #  Assuming grid jobs die on each attempt:
    #    0) canu run from the command line submits iteration 1; canuIteration is NOT incremented
    #       because no parallel jobs have been submitted.
    #    1) Iteration 1 - canu.pl submits jobs, increments the iteration count, and submits itself as iteration 2
    #    2) Iteration 2 - canu.pl submits jobs, increments the iteration count, and submits itself as iteration 3
    #    3) Iteration 3 - canu.pl fails with the error below
    #
    #  If the jobs succeed in Iteration 2, the canu in iteration 3 will pass the Check(), never call
    #  this function, and continue the pipeline.

    my $iter = getGlobal("canuIteration");
    my $max  = getGlobal("canuIterationMax");

    if ($iter >= $max) {
        caExit("canu iteration count too high, stopping pipeline (most likely a problem in the grid-based computes)", undef);
    } elsif ($iter == 0) {
        $iter = "First";
    } elsif ($iter == 1) {
        $iter = "Second";
    } elsif ($iter == 2) {
        $iter = "Third";
    } elsif ($iter == 3) {
        $iter = "Fourth";
    } elsif ($iter == 4) {
        $iter = "Fifth";
    } else {
        $iter = "${iter}th";
    }

    print STDERR "--\n";
    print STDERR "-- Running jobs.  $iter attempt out of $max.\n";

    setGlobal("canuIteration", getGlobal("canuIteration") + 1);

    #  If 'gridEngineJobID' environment variable exists (SGE: JOB_ID; LSF: LSB_JOBID) then we are
    #  currently running under grid crontrol.  If so, run the grid command to submit more jobs, then
    #  submit ourself back to the grid.  If not, tell the user to run the grid command by hand.

    #  Jobs under grid control, and we submit them

    if (defined(getGlobal("gridEngine")) &&
        (getGlobal("useGrid") eq "1") &&
        (getGlobal("useGrid$jobType") eq "1") &&
        (exists($ENV{getGlobal("gridEngineJobID")})) &&
        ($runDirectly == 0)) {
        my @jobsSubmitted;

        print STDERR "--\n";

        purgeGridJobSubmitScripts($path, $script);

        if (getGlobal("showNext")) {
            print STDERR "--\n";
            print STDERR "-- NEXT COMMANDS\n";
            print STDERR "--\n";
            print STDERR "\n";
            print STDERR prettifyCommand("cd $path"), "\n";
            foreach my $j (@jobs) {
                my ($cmd, $jobName) = buildGridJob($asm, $jobType, $path, $script, $mem, $thr, $dsk, $j, undef);

                print STDERR prettifyCommand("./$cmd.sh") . "\n";
            }
            exit(0);
        }

        foreach my $j (@jobs) {
            my ($cmd, $jobName) = buildGridJob($asm, $jobType, $path, $script, $mem, $thr, $dsk, $j, undef);

            if (runCommandSilently($path, "./$cmd.sh", 0)) {
                print STDERR "-- Failed to submit compute jobs.  Delay 10 seconds and try again.\n";
                sleep(10);

                runCommandSilently($path, "./$cmd.sh", 0) and caFailure("Failed to submit compute jobs", "$path/$cmd.out");
            }

            #  Parse the stdout/stderr from the submit command to find the id of the job
            #  we just submitted.  We'll use this to hold the next iteration until all these
            #  jobs have completed.

            open(F, "< $path/$cmd.out");
            while (<F>) {
                chomp;

                if (uc(getGlobal("gridEngine")) eq "SGE") {
                    #  Your job 148364 ("canu_asm") has been submitted
                    if (m/Your\sjob\s(\d+)\s/) {
                        $jobName = $1;
                    }
                    #  Your job-array 148678.1500-1534:1 ("canu_asm") has been submitted
                    if (m/Your\sjob-array\s(\d+).\d+-\d+:\d\s/) {
                        $jobName = $1;
                    }
                }

                if (uc(getGlobal("gridEngine")) eq "LSF") {
                    #  Job <759810> is submitted to queue <14>.
                    if (m/Job\s<(\d+)>\sis/) {
                        $jobName = "ended($1)";
                    }
                }

                if (uc(getGlobal("gridEngine")) eq "PBS") {
                    #  123456.qm2
                    $jobName = $_;
                }

                if (uc(getGlobal("gridEngine")) eq "PBSPRO") {
                    #  ??
                    $jobName = $_;
                }

                if (uc(getGlobal("gridEngine")) eq "SLURM") {
                    #  BPW has seen Slurm report "ERROR" instead of something
                    #  useful here.  If that is seen, report the error to the
                    #  screen and ignore this job.  We'll redo it on the next
                    #  iteration (unless this is the second iteration, then
                    #  we're screwed either way).
                    if (m/Submitted\sbatch\sjob\s(\d+)/) {
                        $jobName = $1;
                    } elsif (m/ERROR/) {
                        $jobName = undef;
                    } else {
                        $jobName = $_;
                    }
                }

                if (uc(getGlobal("gridEngine")) eq "DNANEXUS") {
                   $jobName = $_;
                }
            }
            close(F);

            if      (!defined($jobName)) {
                print STDERR "-- '$cmd.sh' -> returned an error; job not submitted.\n";
            } elsif ($j =~ m/^\d+$/) {
                print STDERR "-- '$cmd.sh' -> job $jobName task $j.\n";
            } else {
                print STDERR "-- '$cmd.sh' -> job $jobName tasks $j.\n";
            }

            if (defined($jobName)) {
                push @jobsSubmitted, $jobName;
            }
        }

        print STDERR "--\n";

        #  All jobs submitted.  Make an option to hold the executive on those jobs.

        my $jobHold;

        if (uc(getGlobal("gridEngine")) eq "SGE") {
            $jobHold = "-hold_jid " . join ",", @jobsSubmitted;
        }

        if (uc(getGlobal("gridEngine")) eq "LSF") {
            $jobHold = "-w \"" . (join "&&", @jobsSubmitted) . "\"";
        }

        if (uc(getGlobal("gridEngine")) eq "PBS") {
            # new PBS versions dont have 1-task arrays like PBSPro but still have afteranyarray (which doesn't work on a not-array task)
            # so we need to check if we are waiting for a regular job or array
            my $holdType = (join ":", @jobsSubmitted)  =~ m/^(\d+)\[(.*)\]/ ? "afteranyarray" : "afterany";
            $jobHold = "-W depend=$holdType:" . join ":", @jobsSubmitted;
        }

        if (uc(getGlobal("gridEngine")) eq "PBSPRO") {
            $jobHold = "-W depend=afterany:" . join ":", @jobsSubmitted;
        }

        if (uc(getGlobal("gridEngine")) eq "SLURM") {
            $jobHold = "--depend=afterany:" . join ":", @jobsSubmitted;
        }

        if (uc(getGlobal("gridEngine")) eq "DNANEXUS") {
            $jobHold = "--depends-on " . join " ", @jobsSubmitted;
        }

        submitScript($asm, $jobHold);

        #  submitScript() should never return.  If it does, then a parallel step was attempted too many time.

        caExit("Too many attempts to run a parallel stage on the grid.  Stop.", undef);
    }

    #  Jobs under grid control, but the user must submit them

    if (defined(getGlobal("gridEngine")) &&
        (getGlobal("useGrid") ne "0") &&
        (getGlobal("useGrid$jobType") eq "1") &&
        (! exists($ENV{getGlobal("gridEngineJobID")})) &&
        ($runDirectly == 0)) {
        my $cwd = getcwd();
        my $s   = (scalar(@jobs) == 1) ? "" : "s";

        print STDERR "\n";
        print STDERR "Please run the following command$s to submit tasks to the grid for execution.\n";
        print STDERR "Each task will use $mem gigabytes memory and $thr threads.\n";
        print STDERR "\n";
        print STDERR "  cd $cwd/$path\n";

        purgeGridJobSubmitScripts($path, $script);

        foreach my $j (@jobs) {
            my ($cmd, $jobName) = buildGridJob($asm, $jobType, $path, $script, $mem, $thr, $dsk, $j, undef);

            print "  ./$cmd.sh\n";
        }

        print STDERR "\n";
        print STDERR "When all tasks are finished, restart canu as before.  The output of the grid\n";
        print STDERR "submit command$s will be in *jobSubmit*out.\n";
        print STDERR "\n";

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

        if (getGlobal("showNext")) {
            print STDERR "--\n";
            print STDERR "-- NEXT COMMANDS\n";
            print STDERR "--\n";
            print STDERR "\n";
            print STDERR prettifyCommand("cd $path") . "\n";
            for (my $i=$st; $i<=$ed; $i++) {
                print STDERR prettifyCommand("./$script.sh $i") . "\n";
            }
            exit(0);
        }

        for (my $i=$st; $i<=$ed; $i++) {
            schedulerSubmit("./$script.sh $i > ./" . buildOutputName($path, $script, $i) . " 2>&1");
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
    schedulerFinish($path, $jobType);
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

    print STDERR "\n";
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

    return(0)  if ($cmd eq "");

    #  Check if the directory exists.

    if (! -d $dir) {
        caFailure("Directory '$dir' doesn't exist, can't run command", "");
    }

    #  If only showing the next command, show it and stop.

    if (getGlobal("showNext")) {
        print STDERR "--\n";
        print STDERR "-- NEXT COMMAND\n";
        print STDERR "--\n";
        print STDERR "\n";
        print STDERR prettifyCommand("cd $dir") . "\n";
        print STDERR "$dis\n";
        exit(0);
    }

    #  Log that we're starting, and show the pretty-ified command.

    my $cwd = getcwd();        #  Remember where we are.
    chdir($dir);               #  So we can root the jobs in the correct location.

    my $startsecs = time();
    my $diskfree  = diskSpace(".");

    print STDERR "----------------------------------------\n";
    print STDERR "-- Starting command on ", scalar(localtime()), " with $diskfree GB free disk space\n";
    print STDERR "\n";
    print STDERR "    cd $dir\n";
    print STDERR "$dis\n";

    my $rc = 0xffff & system($cmd);

    logFinished(".", $startsecs);

    chdir($cwd);

    #  Pretty much copied from Programming Perl page 230

    return(0) if ($rc == 0);

    reportRunError($rc);

    return(1);
}



#  Duplicated in Grid_Cloud.pm to get around recursive 'use' statements.

sub runCommandSilently ($$$) {
    my $dir      = shift @_;
    my $cmd      = shift @_;
    my $dis      = prettifyCommand($cmd);
    my $critical = shift @_;

    return(0)   if ($cmd eq "");

    my $cwd       = getcwd();  #  Remember where we are.
    chdir($dir);               #  So we can root the jobs in the correct location.

    my $rc = 0xffff & system($cmd);

    chdir($cwd);

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
    my  $asm     = getGlobal("onExitNam");
    my  $msg     = shift @_;
    my  $log     = shift @_;
    my  $version = getGlobal("version");

    $msg = undef   if ($msg eq "");
    $log = undef   if ($log eq "");

    print STDERR "\n";
    print STDERR "ABORT:\n";
    print STDERR "ABORT: $version\n";
    print STDERR "ABORT: Don't panic, but a mostly harmless error occurred and Canu stopped.\n";
    print STDERR "ABORT: Try restarting.  If that doesn't work, ask for help.\n";
    print STDERR "ABORT:\n";
    print STDERR "ABORT:   $msg.\n"     if (defined($msg));
    print STDERR "ABORT:\n"             if (defined($msg));

    if (defined($log) && -e $log) {
        my  $df = diskSpace($log);

        print STDERR "ABORT: Disk space available:  $df GB\n";
        print STDERR "ABORT:\n";
    }

    if (-e $log) {
        print STDERR "ABORT: Last 50 lines of the relevant log file ($log):\n";
        print STDERR "ABORT:\n";

        open(Z, "tail -n 50 $log |");
        while (<Z>) {
            print STDERR "ABORT:   $_";
        }
        close(Z);

        print STDERR "ABORT:\n";
    }

    my $fail = getGlobal('onFailure');
    if (defined($fail)) {
        runCommandSilently(getGlobal("onExitDir"), "$fail $asm", 0);
    }

    exit(1);
}


#  Use caFailure() for errors that definitely will require code changes to fix.
sub caFailure ($$) {
    my  $asm     = getGlobal("onExitNam");
    my  $msg     = shift @_;
    my  $log     = shift @_;
    my  $version = getGlobal("version");
    my  $trace   = longmess("Failed");

    $trace =~ s/\n/\nCRASH: /g;

    print STDERR "\n";
    print STDERR "CRASH:\n";
    print STDERR "CRASH: $version\n";
    print STDERR "CRASH: Please panic, this is abnormal.\n";
    print STDERR "CRASH:\n";
    print STDERR "CRASH:   $msg.\n";
    print STDERR "CRASH:\n";
    print STDERR "CRASH: $trace\n";
    #print STDERR "CRASH:\n";   #  $trace has an extra CRASH: at the end

    if (-e $log) {
        print STDERR "CRASH: Last 50 lines of the relevant log file ($log):\n";
        print STDERR "CRASH:\n";

        open(Z, "tail -n 50 $log |");
        while (<Z>) {
            print STDERR "CRASH: $_";
        }
        close(Z);

        print STDERR "CRASH:\n";
    } else {
        print STDERR "CRASH: No log file supplied.\n";
        print STDERR "CRASH:\n";
    }

    my $fail = getGlobal('onFailure');
    if (defined($fail)) {
        runCommandSilently(getGlobal("onExitDir"), "$fail $asm", 0);
    }

    exit(1);
}


1;
