
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
 #    Brian P. Walenz beginning on 2015-FEB-27
 #      are Copyright 2015 Battelle National Biodefense Institute, and
 #      are subject to the BSD 3-Clause License
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

package canu::Execution;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(stopBefore stopAfter skipStage emitStage touch getInstallDirectory getBinDirectory getBinDirectoryShellCode submitScript submitOrRunParallelJob runCommand runCommandSilently);

use strict;
use Config;            #  for @signame

use POSIX ":sys_wait_h";  #  For waitpid(..., &WNOHANG)

use File::Path qw(make_path remove_tree);

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

    my $startsecs = time();
    my $diskfree  = (defined($dir)) ? (diskSpace($dir)) : (0);

    print STDERR "--\n";
    print STDERR "-- Starting concurrent execution on ", scalar(localtime()), " with $diskfree GB free disk space ($remain processes; $numberOfProcesses concurrently)\n"  if  (defined($dir));
    print STDERR "-- Starting concurrent execution on ", scalar(localtime()), " ($remain processes; $numberOfProcesses concurrently)\n"                                    if (!defined($dir));

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
    print STDERR "--\n";
    print STDERR "-- Finished on ", scalar(localtime()), " (", time() - $startsecs, " seconds) with $diskfree GB free disk space$warning\n";
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
        print STDERR "Stop requested before '$stopBefore'.\n";
        print STDERR "Command:\n$cmd\n" if (defined($cmd));
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


sub lookupStageLabel ($) {
    my $label  = shift @_;
    my %ckp;
    my $index;

    #  For correction

    $index = 100;

    $ckp{'cor-gatekeeper'}                  = $index++;

    $ckp{'cor-meryl'}                       = $index++;

    $ckp{'cor-mhapConfigure'}               = $index++;
    $ckp{'cor-mhapPrecomputeCheck'}         = $index++;  # + attempt
    $ckp{'cor-mhapCheck'}                   = $index++;  # + attempt
    $ckp{'cor-overlapConfigure'}            = $index++;
    $ckp{'cor-overlapCheck'}                = $index++;  # + attempt

    $ckp{'cor-overlapStoreConfigure'}       = $index++;
    $ckp{'cor-overlapStoreBucketizerCheck'} = $index++;
    $ckp{'cor-overlapStoreSorterCheck'}     = $index++;
    $ckp{'cor-createOverlapStore'}          = $index++;

    $ckp{'cor-buildCorrectionLayouts'}      = $index++;
    $ckp{'cor-generateCorrectedReads'}      = $index++;  # + attempt
    $ckp{'cor-dumpCorrectedReads'}          = $index++;

    #  For trimming

    $index = 200;

    $ckp{'obt-gatekeeper'}                  = $index++;

    $ckp{'obt-meryl'}                       = $index++;

    $ckp{'obt-mhapConfigure'}               = $index++;
    $ckp{'obt-mhapPrecomputeCheck'}         = $index++;  # + attempt
    $ckp{'obt-overlapConfigure'}            = $index++;
    $ckp{'obt-overlapCheck'}                = $index++;  # + attempt

    $ckp{'obt-overlapStoreConfigure'}       = $index++;
    $ckp{'obt-overlapStoreBucketizerCheck'} = $index++;
    $ckp{'obt-overlapStoreSorterCheck'}     = $index++;
    $ckp{'obt-createOverlapStore'}          = $index++;

    $ckp{'obt-trimReads'}                   = $index++;
    $ckp{'obt-splitReads'}                  = $index++;
    $ckp{'obt-dumpReads'}                   = $index++;  #  rename this

    #  For assembly

    $index = 300;

    $ckp{'utg-gatekeeper'}                  = $index++;

    $ckp{'utg-meryl'}                       = $index++;

    $ckp{'utg-mhapConfigure'}               = $index++;
    $ckp{'utg-mhapPrecomputeCheck'}         = $index++;  # + attempt
    $ckp{'utg-overlapConfigure'}            = $index++;
    $ckp{'utg-overlapCheck'}                = $index++;  # + attempt

    $ckp{'utg-overlapStoreConfigure'}       = $index++;
    $ckp{'utg-overlapStoreBucketizerCheck'} = $index++;
    $ckp{'utg-overlapStoreSorterCheck'}     = $index++;
    $ckp{'utg-createOverlapStore'}          = $index++;

    $ckp{'overlapFilterDetectConfigure'}    = $index++;
    $ckp{'overlapFilterDetectCheck'}        = $index++;  # + attempt
    $ckp{'overlapFilterConfigure'}          = $index++;
    $ckp{'overlapFilter'}                   = $index++;  # + attempt

    $ckp{'readErrorDetectionConfigure'}     = $index++;
    $ckp{'readErrorDetectionCheck'}         = $index++;  # + attempt
    $ckp{'overlapErrorAdjustmentConfigure'} = $index++;
    $ckp{'overlapErrorAdjustmentCheck'}     = $index++;  # + attempt
    $ckp{'updateOverlapStore'}              = $index++;

    $ckp{'unitig'}                          = $index++;

    $ckp{'consensusConfigure'}              = $index++;
    $ckp{'consensusCheck'}                  = $index++;  # + attempt
    $ckp{'consensusLoad'}                   = $index++;
    $ckp{'consensusFilter'}                 = $index++;

    $ckp{'outputLayout'}                    = $index++;
    $ckp{'outputGraph'}                     = $index++;
    $ckp{'outputSequence'}                  = $index++;

    caFailure("invalid checkpoint label '$label'", undef)  if (!defined($ckp{$label}));

    return($ckp{$label});
}



#  Returns true if we should skip this stage.  No used, but left in for possible use in cleaning up things.  Signals that we're
#  at the start of some stage, and we could clean up earlier stages.
#
sub skipStage ($$$@) {
    my $wrk         = shift @_;
    my $asm         = shift @_;
    my $stage       = shift @_;
    my $attempt     = shift @_;

    my $ckpstage    = "";
    my $ckpattempt  = undef;

    #  DISABLED.
    return(0);

    if (! -e "$wrk/$asm.stage") {
        #  No checkpoint file exists, must compute!
        print STDERR "No $wrk/$asm.stage file, compute it all!\n";
        return(0);
    }

    open(F, "< $wrk/$asm.stage") or caFailure("failed to open '$wrk/$asm.stage' for reading", undef);
    while (<F>) {
        if (m/canu\s+at\s+stage\s+(\S*)\s+\(#\d+\)\sattempt\s+(\d+)$/) {
            $ckpstage   = $1;
            $ckpattempt = $2;
        }
        if (m/canu\s+at\s+stage\s+(\S*)\s+\(#\d+\)$/) {
            $ckpstage = $1;
        }
    }
    close(F);

    caFailure("didn't find stage in '$wrk/$asm.stage'", undef)  if ($ckpstage eq "");

    #  Don't skip it.  The stage to run is after the stage in the checkpoint.
    if (lookupStageLabel($ckpstage) < lookupStageLabel($stage)) {
        #print STDERR "PURGE at stage $stage (ignored attempt $attempt)\n";
        purgeRecomputable($wrk, $asm, $stage);
        return(0);
    };

    #  Don't skip it.  The stage to run is the stage in the checkpoint, but the attempt we're trying
    #  is after.
    if ((lookupStageLabel($stage) == lookupStageLabel($ckpstage) &&
         (defined($attempt)) &&
         (defined($ckpattempt)) &&
         ($ckpattempt < $attempt))) {
        #print STDERR "PURGE at stage $stage attempt $attempt\n";
        purgeRecomputable($wrk, $asm, $stage);
        return(0);
    };

    #print STDERR "skipStage()-- target $stage/" . lookupStageLabel($stage);
    #print STDERR " attempt $attempt"     if (defined($attempt));
    #print STDERR " -- checkpoint $ckpstage/" . lookupStageLabel($ckpstage);
    #print STDERR " attempt $ckpattempt"  if (defined($ckpattempt));
    #print STDERR " -- from '$wrk/$asm.stage'\n";

    #  Skip it.  But first, purge any files we'll never need again.

    print STDERR "PURGE extraneous before stage $stage attempt $attempt\n";
    purgeExtraneous($wrk, $asm, $stage);

    return(1);
}


#  Same as skipStage(), left in for future use cleaning up.  Signals that we're done with a stage.
#
sub emitStage ($$$@) {

    return;

    my $wrk         = shift @_;
    my $asm         = shift @_;
    my $stage       = shift @_;
    my $attempt     = shift @_;
    my $time        = localtime();

    my $label       = lookupStageLabel($stage);
    my $label1      = $label - 1;
    my $attempt     = (defined($attempt)) ? " attempt $attempt" : "";
    my $ATTEMPT     = (defined($attempt)) ? " ATTEMPT $attempt" : "";

    open(F, ">> $wrk/$asm.stage") or caFailure("failed to open '$wrk/$asm.stage' for appending\n", undef);
    print F "$time -- canu at stage $stage (#$label)$attempt\n";
    close(F);

    #print "----------------------------------------STAGE $stage (#$label)$ATTEMPT FINISHED.\n";

    make_path("$wrk/$asm.stage.fileLists")  if (! -d "$wrk/$asm.stage.fileLists");

    #  Find all files created since the last checkpoint, or accessed since the last.  Linux claims either -anewer or -newera will work
    #
    #  Format is: access-time -- modification-time -- status-change-time -- filename

    #  Problem - label-1 doesn't always exist because we occasionally reset the label to 200, 300.  Fixed by not doing that.

    while (($label1 > 100) && (! -e "$wrk/$asm.stage.fileLists/stage.$label1.created")) {
        $label1--;
    }

    if (-e "$wrk/$asm.stage.fileLists/stage.$label1.created") {
        runCommandSilently($wrk, "find . -type f -and -newer  $wrk/$asm.stage.fileLists/stage.$label1.created -print > $wrk/$asm.stage.fileLists/stage.$label.created");
        runCommandSilently($wrk, "find . -type f -and -anewer $wrk/$asm.stage.fileLists/stage.$label1.created -print > $wrk/$asm.stage.fileLists/stage.$label.accessed");
    } else {
        runCommandSilently($wrk, "find . -type f -print > $wrk/$asm.stage.fileLists/stage.$label.created");
    }
}


#
#my %firstAccessed;
#my %lastAccessed;
#
#
#sub readAccessed ($$$$) {
#    my $wrk   = shift @_;
#    my $asm   = shift @_;
#    my $stage = shift @_;
#    my $file  = shift @_;
#
#    open(G, "< $file") or die "Failed to open '$file' for reading: $!\n";
#    while (<G>) {
#        chomp;
#
#        s!^\.\/!!;
#        s!$asm\.!PREFIX.!;
#        #s!\d\d\d\d\d\d\.out!TASKID.out!;
#        s!\d\d\d\d\d\d\.!TASKID.!;
#
#        next  if (m/PREFIX.stage.fileLists/);
#        next  if (m/PREFIX.stage/);
#        next  if (m/runCA-logs/);
#
#        $firstAccessed{$_} = $stage  if (!exists($firstAccessed{$_}));
#        $lastAccessed{$_}  = $stage;
#    }
#    close(G);
#}
#
#
##  Remove files we won't ever need again.
#sub purgeExtraneous ($$$$) {
#    my $wrk         = shift @_;
#    my $asm         = shift @_;
#    my $stage       = shift @_;
#
#    open(F, "ls $wrk/$asm.stage.fileLists/ |");
#    while (<F>) {
#        my $file = $_;  chomp $file;
#
#        if ($file =~ m/stage.(\d+).accessed/) {
#            readAccessed($wrk, $asm, $1, "$wrk/test.stage.fileLists/$file");
#        }
#    }
#    close(F);
#
#    foreach my $f (keys %lastAccessed) {
#        next  if ($lastAccessed{$f} < $stage);
#
#        #print STDERR "AT stage $stage REMOVE extraneous file $f last needed in stage $lastAccessed{$f}\n";
#    }
#}
#
#
#
##  Remove files we'll be computing again.  This is based on knowing the current stage we're at,
##  then examining the sge.fileLists directory for future stages and removing those files.
##
##  If any are detected, we emit a new stage (to find files between the last stage and the last stop),
##  then delete those files.
##
#sub purgeRecomputable ($$$) {
#    my $wrk         = shift @_;
#    my $asm         = shift @_;
#    my $stage       = shift @_;
#
#    open(F, "ls $wrk/$asm.stage.fileLists/ |");
#    while (<F>) {
#        my $file = $_;  chomp $file;
#
#        if ($file =~ m/stage.(\d+).accessed/) {
#            readAccessed($wrk, $asm, $1, "$wrk/test.stage.fileLists/$file");
#        }
#    }
#    close(F);
#
#    foreach my $f (keys %firstAccessed) {
#        next  if ($stage < $firstAccessed{$f});
#
#        print STDERR "AT stage $stage REMOVE premature file $f first needed in stage $firstAccessed{$f}\n";
#    }
#}
#


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
    my @t = split '/', "$FindBin::RealBin";
    pop @t;                         #  bin
    pop @t;                         #  arch, e.g., FreeBSD-amd64
    my $installDir = join '/', @t;  #  path to the assembler

    return($installDir);
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

    #  If not requested to run on the grid, or can't run on the grid, fail.

    return   if (getGlobal("useGrid")       == 0);
    return   if (getGlobal("useGridMaster") == 0);
    return   if (getGlobal("gridEngine")    eq undef);

    #  If no job to wait on, and we are already on the grid, do NOT resubmit ourself.
    #
    #  When the user launches canu on the head node, a call to submitScript() is made to launch canu
    #  under grid control.  That results in a restart of canu, and another call to submitScript(),
    #  but this time, the envorinment variable is set, we we can skip the resubmission, and continue
    #  with canu execution.

    return   if (($jobToWaitOn eq undef) && (exists($ENV{getGlobal("gridEngineJobID")})));

    #  Find the next available output file.

    my $idx = "01";

    while (-e "$wrk/canu.$idx.out") {
        $idx++;
    }

    my $output    = "$wrk/canu.$idx.out";
    my $script    = "$wrk/canu.$idx.sh";
    my $iteration = getGlobal("canuIteration") + 1;

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
    print F "\n";
    print F getBinDirectoryShellCode();
    print F "\n";
    print F "/usr/bin/env perl \\\n";
    print F "\$bin/canu " . getCommandLineOptions() . " canuIteration=$iteration\n";
    close(F);

    system("chmod +x $script");

    #  Construct a submission command line.

    my ($jobName, $memOption, $thrOption, $gridOpts);

    $jobName   = "canu_" . $asm . ((defined(getGlobal("gridOptionsJobName"))) ? ("_" . getGlobal("gridOptionsJobName")) : (""));

    $memOption = buildMemoryOption(getGlobal("masterMemory"), getGlobal("masterThreads"));
    $thrOption = buildThreadOption(getGlobal("masterThreads"));

    $gridOpts  = getGlobal("gridOptions")          if (defined(getGlobal("gridOptions")));
    $gridOpts .= " "                               if (defined($gridOpts));
    $gridOpts  = getGlobal("gridOptionsMaster")    if (defined(getGlobal("gridOptionsMaster")));
    $gridOpts .= " "                               if (defined($gridOpts));
    $gridOpts .= $memOption                        if (defined($memOption));
    $gridOpts .= " "                               if (defined($gridOpts));
    $gridOpts .= $thrOption                        if (defined($thrOption));

    #  If the jobToWaitOn is defined, make the script wait for that to complete.  LSF might need to
    #  query jobs in the queue and figure out the job ID (or IDs) for the jobToWaitOn.  Reading LSF
    #  docs online (for bsub.1) claim that we can still use jobToWaitOn.

    if (defined($jobToWaitOn)) {
        (my $hold = getGlobal("gridEngineHoldOption")) =~ s/WAIT_TAG/$jobToWaitOn/;
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


sub buildMemoryOption ($$) {
    my $m = shift @_;
    my $t = shift @_;
    my $r;

    if (getGlobal("gridEngine") eq "SGE") {
        $m /= $t;
    }

    $r =  getGlobal("gridEngineMemoryOption");
    $r =~ s/MEMORY/${m}g/;

    return($r);
}


sub buildThreadOption ($) {
    my $t = shift @_;
    my $r;

    $r =  getGlobal("gridEngineThreadsOption");
    $r =~ s/THREADS/$t/;

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

    my $submitCommand = getGlobal("gridEngineSubmitCommand");
    my $nameOption    = getGlobal("gridEngineNameOption");

    my $jobNameT = "${jobType}_" . $asm . ((defined(getGlobal("gridOptionsJobName"))) ? ("_" . getGlobal("gridOptionsJobName")) : (""));

    my $jobName       = buildGridArray($jobNameT, $bgnJob, $endJob, getGlobal("gridEngineArrayName"));
    my $arrayOpt      = buildGridArray($jobNameT, $bgnJob, $endJob, getGlobal("gridEngineArrayOption"));

    my $outputOption  = getGlobal("gridEngineOutputOption");
    my $tid           = getGlobal("gridEngineArraySubmitID");
    my $outName       = buildOutputName($path, $script, $tid);

    my $memOption     = buildMemoryOption($mem, $thr);
    my $thrOption     = buildThreadOption($thr);

    my $gridOpts;

    $gridOpts  = getGlobal("gridOptions")          if (defined(getGlobal("gridOptions")));
    $gridOpts .= " "                               if (defined($gridOpts));
    $gridOpts  = getGlobal("gridOptions$jobType")  if (defined(getGlobal("gridOptions$jobType")));
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
sub submitOrRunParallelJob ($$$$$@) {
    my $wrk          = shift @_;  #  Root of the assembly (NOT wrk/correction or wrk/trimming)
    my $asm          = shift @_;  #  Name of the assembly

    my $jobType      = shift @_;  #  E.g., ovl, cns, ... - populates 'gridOptionsXXX and useGridXXX (former XXXOnGrid)
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

    #  Break infinite loops.  If the grid jobs keep failing, give up after a few attempts.
    #
    #  If useGridMaster = 0, this has no impact; canuIteration is incremented in submitScript() and
    #  the initial value is zero.  It is reset during the Check() operation on each parallel step.
    #
    #  submitScript() will increment canuIteration on each call.  It is reset to zero if the Check()
    #  for any parallel step succeeds.
    #
    #  Assuming grid jobs die on each attempt:
    #    Iteration 0 - the one run on the command line; submits iteration 1
    #    Iteration 1 - run on the grid, submits parallel jobs and iteration 2
    #    Iteration 2 - run on the grid, submits parallel jobs and iteration 3
    #    Iteration 3 - run on the grid, fails
    #
    #  If the jobs succeed in Iteration 2, the canu in iteration 3 will pass the Check(), and
    #  continue the pipeline.

    caExit("canu iteration count too high, stopping pipeline (most likely a problem in the grid-based computes)", undef)
        if (getGlobal("canuIteration") > getGlobal("canuIterationMax"));

    #  If 'gridEngineJobID' environment variable exists (SGE: JOB_ID; LSF: LSB_JOBID) then we are
    #  currently running under grid crontrol.  If so, run the grid command to submit more jobs, then
    #  submit ourself back to the grid.  If not, tell the user to run the grid command by hand.

    #  Jobs under grid control, and we submit them

    if (getGlobal("gridEngine") && getGlobal("useGrid") && getGlobal("useGrid$jobType") && (exists($ENV{getGlobal("gridEngineJobID")}))) {
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

    if (getGlobal("gridEngine") && getGlobal("useGrid") && getGlobal("useGrid$jobType") && (! exists($ENV{getGlobal("gridEngineJobID")}))) {
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

    my $nParallel  = getGlobal("${jobType}Concurrency");
    $nParallel     = int(getNumberOfCPUs() / $thr)  if ((!defined($nParallel)) || ($nParallel == 0));

    schedulerSetNumberOfProcesses($nParallel);
    schedulerFinish($path);
}





#  Utility to run a command and check the exit status, report time used.
#
sub runCommand ($$) {
    my $dir = shift @_;
    my $cmd = shift @_;
    my $dis = $cmd;

    return  if ($cmd eq "");

    #  Pretty-ify the command.  If there are no newlines already in it, break
    #  before every switch and before file redirects.

    if (($dis =~ tr/\n/\n/) == 0) {
        $dis =~ s/\s-/ \\\n  -/g;
        $dis =~ s/\s>\s/ \\\n> /;
        $dis =~ s/\s2>\s/ \\\n2> /;
    }

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

    print STDERR "--\n";
    print STDERR "-- Starting command on ", scalar(localtime()), " with $diskfree GB free disk space\n"  if  (defined($dir));
    print STDERR "-- Starting command on ", scalar(localtime()), "\n"                                    if (!defined($dir));
    print STDERR "--\n";
    print STDERR "$dis\n";

    my $rc = 0xffff & system("cd $dir && $cmd");

    $diskfree = (defined($dir)) ? (diskSpace($dir)) : (0);

    my $warning = "  !!! WARNING !!!" if ($diskfree < 10);
    print STDERR "--\n";
    print STDERR "-- Finished on ", scalar(localtime()), " (", time() - $startsecs, " seconds) with $diskfree GB free disk space$warning\n";

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



sub runCommandSilently ($$) {
    my $dir = shift @_;
    my $cmd = shift @_;

    return  if ($cmd eq "");

    if (! -d $dir) {
        caFailure("Directory '$dir' doesn't exist, can't run command", "");
    }

    my $rc = 0xffff & system("cd $dir && $cmd");

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



1;
