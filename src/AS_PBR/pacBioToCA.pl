#!/usr/bin/perl

#   Copyright (C) 2011, Battelle National Biodefense Institute (BNBI);
#   all rights reserved. Authored by: Sergey Koren
#   
#   This Software was prepared for the Department of Homeland Security
#   (DHS) by the Battelle National Biodefense Institute, LLC (BNBI) as
#   part of contract HSHQDC-07-C-00020 to manage and operate the National
#   Biodefense Analysis and Countermeasures Center (NBACC), a Federally
#   Funded Research and Development Center.
#   
#   Redistribution and use in source and binary forms, with or without
#   modification, are permitted provided that the following conditions are
#   met:
#   
#   * Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
#   
#   * Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in the
#     documentation and/or other materials provided with the distribution.
#   
#   * Neither the name of the Battelle National Biodefense Institute nor
#     the names of its contributors may be used to endorse or promote
#     products derived from this software without specific prior written
#     permission.
#   
#   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#   HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
#   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
#   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
#   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
###########################################################################
#
#  Read in fragments from fastq-format sequence and quality files,
#  correct to pacbio fragments.
#

use strict;

use Config;  #  for @signame
use FindBin;
use Cwd;
use Carp;
use FileHandle;

use POSIX qw(ceil floor sys_wait_h);

my $TAB_SIZE = 8;

my %global;



my @nonCAOptions = ("maxCorrectionRounds", "genomeSize", "shortReads", "longReads", "libraryname", "specFile", "length", "coverageCutoff", "maxCoverage", "maxGap", "maxUncorrectedGap", "samFOFN", "blasr", "bowtie", "threads", "repeats", "fastqFile", "partitions", "submitToGrid", "sgeCorrection", "consensusConcurrency", "cleanup");

my $commandLineOptions = "";

sub getGlobal ($) {
    my $var = shift @_;
    if (!exists($global{$var})) {
    	return undef;
    }
    return($global{$var});
}

sub setGlobal ($$) {
    my $var = shift @_;
    my $val = shift @_;

    #  If no value, set the field to undefined, the default for many of the options.

    $val = undef  if ($val eq "");

    $global{$var} = $val;
    $val =~ s/\\\"/\"/g;
    $val =~ s/\"/\\\"/g;
    $val =~ s/\\\$/\$/g;
    $val =~ s/\$/\\\$/g;
    
    # some vars should be global like specFile
    if ($var eq "specFile" || $var eq "fastqFile" || $var eq "samFOFN") {
       $val = makeAbsolute($val);
    }
    $commandLineOptions .= " \"$var=$val\" ";
}

sub setDefaults() {
	# grid options, duplicate of runCA
    $global{"gridSubmitCommand"}		   = "qsub";
    $global{"gridHoldOption"}		       = "-hold_jid \"WAIT_TAG\""; # for lsf it is -w "done("WAIT_TAG")"
    $global{"gridSyncOption"}			   = "-sync y"; # for lsf it is -K
    $global{"gridNameOption"}			   = "-cwd -N";         # for lsf it is -J
    $global{"gridArrayOption"}			   = "-t ARRAY_JOBS";	# for lsf, empty ("")
    $global{"gridArrayName"}			   = "ARRAY_NAME";		# for lsf, it is ARRAY_NAME[ARRAY_JOBS]
    $global{"gridOutputOption"}			   = "-j y -o";         # for lsf, it is -o
    $global{"gridPropagateCommand"}		   = "qalter -hold_jid \"WAIT_TAG\""; # for lsf it is bmodify -w "done(WAIT_TAG)"
    $global{"gridNameToJobIDCommand"}      = undef;             # for lsf it is bjobs -J "WAIT_TAG" | grep -v JOBID    
    $global{"gridTaskID"}				   = "SGE_TASK_ID";     # for lsf it is LSB_JOBINDEX
    $global{"gridArraySubmitID"}           = "\$TASK_ID";       # for lsf it is %I
    
    $global{"shell"}                       = "/bin/sh";
	
    $global{"utgErrorRate"} = 0.25;
    $global{"utgErrorLimit"} = 6.5;
    $global{"ovlErrorRate"} = 0.25;
    $global{"cnsErrorRate"} = 0.25;
    $global{"cgwErrorRate"} = 0.25;
	
    $global{"shortReads"} = undef;
    $global{"longReads"} = undef;
    $global{"libraryname"} = undef;
    $global{"specFile"} = undef;
    $global{"length"} = 500;
	
    $global{"coverageCutoff"} = 0;
    $global{"maxCoverage"} = 50;
    $global{"genomeSize"} = 0;
    $global{"maxUncorrectedGap"} = 0;
    $global{"maxCorrectionRounds"} = 0;
    $global{"repeats"} = "";
    $global{"samFOFN"} = undef;
    $global{"blasr"} = undef;
    $global{"bowtie"} = undef;
	
    $global{"threads"} = 1;
    $global{"partitions"} = 1;
    $global{"submitToGrid"} = 0;
    $global{"sge"} = undef;
    $global{"sgeCorrection"} = undef;
    $global{"consensusConcurrency"} = 8;

    $global{"doOverlapBasedTrimming"} = 0;  # should this be defaulted to true or false?

    $global{"cleanup"} = 1;
		
    $global{"fastqFile"} = undef;
}

sub updateSpecFile($) {
	my $outFile = shift @_;
    open(W, "> $outFile") or die("Couldn't open '$outFile'", undef);

	my $print = 1;
	foreach my $key (keys %global) {
    	$print = 1;        
    	foreach my $nonCaOption (@nonCAOptions) {
    		if (index($key, $nonCaOption) == 0) {
    			$print = 0;
    			last;
    		}
    	}
    	if ($print == 1) {
    		print W "$key = " . getGlobal($key) . "\n";
    	}
    }
    close(W);
}

sub setParametersFromFile ($@) {
    my $specFile  = shift @_;
    my @fragFiles = @_;

    #  Client should be ensuring that the file exists before calling this function.
    die "specFile '$specFile' not found.\n"  if (! -e "$specFile");

    print STDERR "\n";
    print STDERR "###\n";
    print STDERR "###  Reading options from '$specFile'\n";
    print STDERR "###\n";
    print STDERR "\n";

    open(F, "< $specFile") or die("Couldn't open '$specFile'", undef);

    while (<F>) {
        s/^\s+//;
        s/\s+$//;

        next if (m/^\s*\#/);
        next if (m/^\s*$/);

        if (-e $_) {
            my $xx = $_;
            $xx = "$ENV{'PWD'}/$xx" if ($xx !~ m!^/!);
            if (-e $xx) {
                push @fragFiles, $xx;
            } else {
                setGlobal("help", getGlobal("help") . "File not found '$_' after appending absolute path.\n");
            }
        } elsif (m/\s*(\w*)\s*=([^#]*)#*.*$/) {
            my ($var, $val) = ($1, $2);
            $var =~ s/^\s+//; $var =~ s/\s+$//;
            $val =~ s/^\s+//; $val =~ s/\s+$//;
            undef $val if ($val eq "undef");
            if ($var eq "maxGap") {
               $var = "maxUncorrectedGap"; 
            }
            setGlobal($var, $val);
        } else {
            setGlobal("help", getGlobal("help") . "File not found or unknown specFile option line '$_'.\n");
        }
    }
    close(F);

    return(@fragFiles);
}

sub setParametersFromCommandLine(@) {
    my @specOpts = @_;

    if (scalar(@specOpts) > 0) {
        print STDERR "\n";
        print STDERR "###\n";
        print STDERR "###  Reading options from the command line.\n";
        print STDERR "###\n";
        print STDERR "\n";
    }

    foreach my $s (@specOpts) {
        if ($s =~ m/\s*(\w*)\s*=(.*)/) {
            my ($var, $val) = ($1, $2);
            $var =~ s/^\s+//; $var =~ s/\s+$//;
            $val =~ s/^\s+//; $val =~ s/\s+$//;
            setGlobal($var, $val);
        } else {
            setGlobal("help", getGlobal("help") . "Misformed command line option '$s'.\n");
        }
    }
}

sub makeAbsolute ($) {
    my $val = shift @_;
    if (defined($val) && ($val !~ m!^/!)) {
        $val = "$ENV{'PWD'}/$val";
    }
    return $val;
}

sub getInstallDirectory () {
    my @t = split '/', "$FindBin::RealBin";
    pop @t;                         #  bin
    pop @t;                         #  arch, e.g., FreeBSD-amd64
    my $installDir = join '/', @t;  #  path to the assembler

    return($installDir);
}

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

    return($path);
}

sub getBinDirectoryShellCode () {
    my $installDir = getInstallDirectory();
    my $string;

    $string  = "\n";
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

    return($string);
}

sub runCommand($$) {
   my $dir = shift;
   my $cmd = shift;

   my $t = localtime();
   my $d = time();
   print STDERR "----------------------------------------START $t\n$cmd\n";

   my $rc = 0xffff & system("cd $dir && $cmd");

   $t = localtime();
   print STDERR "----------------------------------------END $t (", time() - $d, " seconds)\n";

   return(0) if ($rc == 0);

   die "Failed to execute $cmd\n";
}

################################################################################

#  Functions for running multiple processes at the same time.

my $numberOfProcesses       = 0;     #  Number of jobs concurrently running
my $numberOfProcessesToWait = 0;     #  Number of jobs we can leave running at exit
my @processQueue            = ();
my @processesRunning        = ();
my $printProcessCommand     = 1;     #  Show commands as they run

sub schedulerSetNumberOfProcesses {
    $numberOfProcesses = shift @_;
}

sub schedulerSubmit {
    chomp @_;
    push @processQueue, @_;
}

sub schedulerForkProcess {
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

sub schedulerReapProcess {
    my $pid = shift @_;

    if (waitpid($pid, &WNOHANG) > 0) {
        return(1);
    } else {
        return(0);
    }
}

sub schedulerRun {
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

sub schedulerFinish {
    my $child;
    my @newProcesses;
    my $remain;

    my $t = localtime();
    my $d = time();
    print STDERR "----------------------------------------START CONCURRENT $t\n";

    $remain = scalar(@processQueue);

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
}

################################################################################
#
#  Functions for running jobs on grid

sub runningOnGrid () {
	my $taskID = getGlobal("gridTaskID");
    return(defined($ENV{$taskID}));
}

sub findNextScriptOutputFile ($$) {
    my $wrk = shift @_;
    my $prefix = shift @_;
    my $idx = "00";
    while (-e "$wrk/$prefix.$idx.sh") {
        $idx++;
    }
    return("$wrk/$prefix.$idx");
}

sub buildGridArray($$$) {
	my $name = shift @_;
	my $maxLimit = shift @_;
	my $globalValue = shift @_;
	
	my $arrayJobName = getGlobal($globalValue);
	$arrayJobName =~ s/ARRAY_NAME/$name/g;
	$arrayJobName =~ s/ARRAY_JOBS/1-$maxLimit/g;
	
	return $arrayJobName;
}

sub getGridArrayName($$) {
	my $name = shift @_;
	my $maxLimit = shift @_;	
	return buildGridArray($name, $maxLimit, "gridArrayName");
}

sub getGridArrayOption($$) {
	my $name = shift @_;
	my $maxLimit = shift @_;	
	return buildGridArray($name, $maxLimit, "gridArrayOption");	
}

sub writeScriptHeader($$$$$) {
    my $wrk = shift @_;
    my $asm = shift @_;
    my $script = shift @_;
    my $executable = shift @_;
    my $options = shift @_;

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

    print F "/usr/bin/env perl \$bin/$executable $options\n";
    close(F);

    system("chmod +x $script");
}

sub submit($$$$$) {
    my $wrk = shift @_;
    my $script = shift @_;
    my $output = shift @_;
    my $prefix = shift @_;
    my $waitTag = shift @_;

    my $sge         = getGlobal("sge");
    my $sgeScript   = getGlobal("sgeScript");
    my $sgePropHold = getGlobal("sgePropagateHold");

    my $submitCommand        = getGlobal("gridSubmitCommand");
    my $holdOption           = getGlobal("gridHoldOption");
    my $nameOption           = getGlobal("gridNameOption");
    my $outputOption         = getGlobal("gridOutputOption");
    my $holdPropagateCommand = getGlobal("gridPropagateCommand");

    my $jobName = $prefix;
    if (defined($waitTag)) {
       my $hold = $holdOption;
       $hold =~ s/WAIT_TAG/$waitTag/g;
       $waitTag = $hold;
    }

    my $qcmd = "$submitCommand $sge $sgeScript $nameOption \"$jobName\" $waitTag $outputOption $output $script";
    runCommand($wrk, $qcmd) and caFailure("Failed to submit script.\n");

    if (defined($sgePropHold)) {
		if (defined($holdPropagateCommand)) {
           my $translateCmd = getGlobal("gridNameToJobIDCommand");
           
           # translate hold option to job id if necessar
           if (defined($translateCmd) && $translateCmd ne "") {
              my $tcmd = $translateCmd;
              $tcmd =~ s/WAIT_TAG/$sgePropHold/g;
              my $propJobCount = `$tcmd |wc -l`;
              chomp $propJobCount;
              if ($propJobCount != 1) {
                 print STDERR "Warning: multiple IDs for job $sgePropHold got $propJobCount and should have been 1.\n";
              }
              #my $jobID = `$tcmd |head -n 1 |awk '{print \$1}'`;
              #chomp $jobID;
              #print STDERR "Translated job ID $sgePropHold to be job $jobID\n";
              #$sgePropHold = $jobID;
              open(PROPS, "$tcmd |awk '{print \$1}' | ") or die("Couldn't get list of jobs that need to hold", undef);
              
              # now we can get the job we are holding for
              $tcmd = $translateCmd;
              $tcmd =~ s/WAIT_TAG/$jobName/g;              
              my $holdJobCount = `$tcmd |wc -l`;
              chomp $propJobCount;
              if ($propJobCount != 1) {
                 print STDERR "Warning: multiple IDs for job $jobName got $propJobCount and should have been 1.\n";
              }
              #$jobID = `$tcmd |head -n 1 |awk '{print \$1}'`;
              #chomp $jobID;
              #print STDERR "Translated job ID $sgePropHold to be job $jobID\n";
              #$jobName = $jobID;
              open(HOLDS, "$tcmd |awk '{print \$1}' | ") or die("Couldn't get list of jobs that should be held for", undef);
              
              # loop over all jobs and all sge hold commands to modify the jobs. We have no way to know which is the right one unfortunately               
              while (my $prop = <PROPS>) {
                 while (my $hold = <HOLDS>) {
                  my $hcmd = $holdPropagateCommand;
                  $hcmd =~ s/WAIT_TAG/$hold/g;
                  my $acmd = "$hcmd $prop";
                  print STDERR "Propagating hold to $prop to wait for job $hold\n";
                  system($acmd) and print STDERR "WARNING: Failed to reset hold_jid trigger on '$prop'.\n";                  
                 }
              }
              close(HOLDS);
              close(PROPS);                
           } else {
              $sgePropHold = "\"$sgePropHold\"";
              $holdPropagateCommand =~ s/WAIT_TAG/$jobName/g;
              my $acmd = "$holdPropagateCommand $sgePropHold";
print "The command I am running to alter is $acmd\n";
              system($acmd) and print STDERR "WARNING: Failed to reset hold_jid trigger on '$sgePropHold'.\n";
           }
		} else {
			print STDERR "WARNING: Failed to reset hold '$sgePropHold', not supported on current grid environment.\n";
		}
    }
    
    if (defined($sgePropHold)) {
       if (defined($holdPropagateCommand) && defined($holdOption)) {
          $holdOption =~ s/WAIT_TAG//g;
          my $acmd = "$holdPropagateCommand $holdOption \"$prefix\" \"$sgePropHold\"";
          system($acmd) and print STDERR "WARNING: Failed to reset hold_jid trigger on '$sgePropHold'.\n";
       } else {
          print STDERR "WARNING: Failed to reset hold '$sgePropHold', not supported on current grid environment.\n";
       }
    }
}

sub submitScript ($$$) {
    my $wrk = shift @_;
    my $asm = shift @_;
    my $waitTag = shift @_;
    my $libraryname = "temp" . getGlobal("libraryname");

    $wrk =~ s/$libraryname//g;
    my $sgeName     = getGlobal("sgeName");
    $sgeName = "_$sgeName"              if (defined($sgeName));

    return if (getGlobal("scriptOnGrid") == 0);
    my $output = findNextScriptOutputFile("$wrk/$libraryname", "runPBcR.sge.out");
    my $script = "$output.sh";
    writeScriptHeader($wrk, $asm, $script, "pacBioToCA", $commandLineOptions);
    submit($wrk, $script, $output, "pBcR_$asm$sgeName", $waitTag);

    exit(0);
}

sub submitBatchJobs($$$$) {
   my $wrk = shift @_;
   my $asm = shift @_;
   my $SGE = shift @_;
   my $TAG = shift @_;

   if (getGlobal("scriptOnGrid")) {
       runCommand($wrk, $SGE) and die("Failed to submit batch jobs.");
       submitScript($wrk, $asm, $TAG);
   } else {
       print "Please execute:\n$SGE\n";
   }
}

sub submitRunCA($$$$) {
   my $wrk = shift @_;
   my $asm = shift @_;
   my $options = shift @_;
   my $TAG = shift @_;
   my $holdPropagateCommand    = getGlobal("gridPropagateCommand");
   my $holdOption	           = getGlobal("gridHoldOption");
   my $syncOption              = getGlobal("gridSyncOption");
   my $sge                     = getGlobal("sge");

   if (defined($holdPropagateCommand) && $holdPropagateCommand ne "") {
      # do nothing, we can use the hold option to propagate
      $syncOption = "";
   } else {
      print "Warning: grid environment does not support hold propagate, running in sync mode\n";
   }

   if (getGlobal("scriptOnGrid")) {
       my $output = findNextScriptOutputFile($wrk, "runPBcR.runCA.sge.out");
       my $script = "$output.sh";
       writeScriptHeader($wrk, $asm, $script, "runCA", $options);
       submit($wrk, "$syncOption $script", $output, $TAG, undef);
       submitScript($wrk, $asm, $TAG);
   } else {
       print "Please execute:\nrunCA $options\n";
   }
}

sub submitRunPBcR($$$$) {
   my $wrk = shift @_;
   my $asm = shift @_;
   my $options = shift @_;
   my $TAG = shift @_;
   my $holdPropagateCommand    = getGlobal("gridPropagateCommand");
   my $holdOption                  = getGlobal("gridHoldOption");
   my $syncOption              = getGlobal("gridSyncOption");
   my $sge                     = getGlobal("sge");

   if (defined($holdPropagateCommand) && $holdPropagateCommand ne "") {
      # do nothing, we can use the hold option to propagate
      $syncOption = "";
   } else {
      print "Warning: grid environment does not support hold propagate, running in sync mode\n";
   }

   if (getGlobal("scriptOnGrid")) {
       my $output = findNextScriptOutputFile($wrk, "runPBcR.sge.out");
       my $script = "$output.sh";
       writeScriptHeader($wrk, $asm, $script, "pacBioToCA", $options);
       submit($wrk, "$syncOption $script", $output, $TAG, undef);
       submitScript($wrk, $asm, $TAG);
   } else {
       print "Please execute:\npacBioToCA $options\n";
   }
}

################################################################################

my $HISTOGRAM_CUTOFF = 0.985; #seems like updated blasr parameters map more

sub pickMappingThreshold($$) {
    my $histogram = shift @_;
    my $conservative = shift @_;
    my $threshold = 0;
    
    # pick threshold as 2sd (95% of the bases)
    my @hist;
    my $total = 0;
    my $sum = 0;
    my $mean = 0;
    my $variance = 0;

    open(F, "< $histogram") or die("Couldn't open '$histogram'", undef);
    while (<F>) {
       s/^\s+//;
       s/\s+$//;

       next if (m/^\s*\#/);
       next if (m/^\s*$/);
    
       $hist[$_]++;
       $total++;
       my $delta = $_ - $mean;
       $mean += $delta / $total;
       $variance += $delta * ($_ - $mean);
    }
    $variance /= $total;
    my $sd = sqrt($variance);
    close(F);
    
    for (my $i = 0; $i <= $#hist; $i++) {
       $sum += $hist[$i];
       if ($sum / $total > $HISTOGRAM_CUTOFF) {
          $threshold = $i - 1;
          last;
       }
    }
    $threshold = ($threshold < ($mean + 3*$sd) ? floor($mean + 3*$sd) : $threshold);
    if (defined($conservative) && $conservative == 1) {
       $threshold = ceil($mean * 2);
    }
    print STDERR "Picked mapping cutoff of $threshold\n";
    
    return $threshold;
}

my $MIN_FILES_WITHOUT_PARTITIONS = 20;
my $REPEAT_MULTIPLIER = 10;
my $MAX_CORRECTION_COVERAGE = 100;
my $MAX_SPLIT_PERCENTAGE = 0.30;
my $MAX_BLASR_LIMIT = 4294967296;

setDefaults();

my @fragFiles;
my @specOpts;
my @cmdArgs;

my $srcstr;

{
    local $, = " ";
    $srcstr = "$0 @ARGV";
}

# process the spec file and the fragment files
while (scalar(@ARGV) > 0) {
    my $arg = shift @ARGV;
    if      ($arg eq "-s") {
        setGlobal("specFile", shift @ARGV);
    
    } elsif (($arg =~ /\.frg$|frg\.gz$|frg\.bz2$/i)) {
       if (-e $arg) {
          push @fragFiles, makeAbsolute($arg);
       } else {
          print "Invalid file $arg could not be found.\n";
       }

    } elsif ($arg =~ m/=/) {
        push @specOpts, $arg;
    
    } else {
        push @cmdArgs, $arg;
    }
}

# initialize parameters from the spec files
my $err = 0;
if (-e getGlobal("specFile")) {
   @fragFiles = setParametersFromFile(getGlobal("specFile"), @fragFiles);
} else {
	print STDERR "Error: spec file " . getGlobal("specFile") . " does not exist. Double-check your paths and try again.\n";
   $err++;
}
setParametersFromCommandLine(@specOpts);

# finally, command-line parameters take precedence
while (scalar(@cmdArgs) > 0) {
    my $arg = shift @cmdArgs;

    if ($arg eq "-shortReads") {
        setGlobal("shortReads", 1);

    } elsif ($arg eq "-longReads") {
        setGlobal("longReads", 1);

    } elsif ($arg eq "-coverageCutoff") {
        setGlobal("coverageCutoff", shift @cmdArgs);
        if (getGlobal("coverageCutoff") < 0) { setGlobal("coverageCutoff", 0); }

    } elsif ($arg eq "-maxCoverage") {
        setGlobal("maxCoverage", shift @cmdArgs);
        if (getGlobal("maxCoverage") < 0) { setGlobal("maxCoverage", 0); }

    } elsif ($arg eq "-genomeSize") {
       setGlobal("genomeSize", shift @cmdArgs);
       if (getGlobal("genomeSize") < 0) { setGlobal("genomeSize", 0); }

    } elsif ($arg eq "-maxGap") {
        setGlobal("maxUncorrectedGap", shift @cmdArgs);
        if (getGlobal("maxUncorrectedGap") < 0) { setGlobal("maxUncorrectedGap", 0); }

    } elsif ($arg eq "-length") {
        setGlobal("length", shift @cmdArgs);

    } elsif ($arg eq "-repeats") {
       setGlobal("repeats", shift @cmdArgs);

    } elsif ($arg eq "-fastq") {
       setGlobal("fastqFile", shift @cmdArgs);

    } elsif ($arg eq "-t" || $arg eq "-threads") {
       setGlobal("threads", shift @cmdArgs);
       if (getGlobal("threads") <= 0) { setGlobal("threads", 1); }

    } elsif ($arg eq "-l" || $arg eq "-libraryname") {
       setGlobal("libraryname", shift @cmdArgs);

    } elsif ($arg eq "-partitions") {
       setGlobal("partitions", shift @cmdArgs);
    
    } elsif ($arg eq "-sge") {
       setGlobal("sge", shift @cmdArgs);

    } elsif ($arg eq "-sgeCorrection") {
       setGlobal("sgeCorrection", shift @cmdArgs);

    } elsif ($arg eq "-noclean") {
       setGlobal("cleanup", 0);
    
    } else {
       print STDERR "Unknown parameter " + $err + "\n";
       $err++;
    }   
 }
    
if (($err) || ((scalar(@fragFiles) == 0) && (!defined(getGlobal("longReads")) || getGlobal("longReads") == 0)) || (!defined(getGlobal("fastqFile"))) || (!defined(getGlobal("specFile"))) || (!defined(getGlobal("libraryname")))) {
   print STDERR "No frag files specified. Please specify at least one frag file containing high-accuracy data for correction.\n" if (scalar(@fragFiles) == 0);
   print STDERR "No fastq file specified. Please specify a fastq file containing PacBio sequences for correction using the -fastq parameter.\n" if (!defined(getGlobal("fastqFile")));
   print STDERR "No spec file defined. Please specify a spec file using the -s option.\n" if (!defined(getGlobal("specFile")));
   print STDERR "No library name provided. Please specify a name using the -library option.\n" if (!defined(getGlobal("libraryname")));  

    print STDERR "usage: $0 [options] -s spec.file -fastq fastqfile <frg>\n";
    print STDERR "  -length                  Minimum length of PacBio sequences to correct/output.\n";
    print STDERR "  -partitions              Number of partitions for consensus\n";
    print STDERR "  -sge                     Submit consensus jobs to the grid\n";
    print STDERR "  -sgeCorrection           Parameters for the correction step for the grid. This should match the threads specified below, for example by using -pe threaded\n";
    print STDERR "  -libraryname libraryname Name of the library; freeformat text.\n";
    print STDERR "  -threads threads         Number of threads to use for correction.\n";
    print STDERR "  -shortReads              Use if the sequences for correction are 100bp or shorter.\n";

    print STDERR "  -genomeSize		     Specify the approximate genome size. This overwrites the coverage parameter above and will be used to compute estimated coverage instead of automatically estimating.\n\n";
    print STDERR "  -maxCoverage		 Maximum coverage of PacBio sequences to correct. Only the longest sequences adding up to this coverage will be corrected. Requires genomeSize or coverage parameter to be specified\n";

    print STDERR "\nAdvanced options (EXPERT):\n";
    print STDERR "  -coverageCutoff	     Specify the pacBio coverage (integer) used to separate repeat copies instead of automatically estimating.\n";
    print STDERR "  -maxGap		     The maximum uncorrected PacBio gap that will be allowed. When there is no short-read coverage for a region, by default the pipeline will split a PacBio sequence. This option allows a number of PacBio sequences without short-read coverage to remain. For example, specifying 50, will mean 50bp can have no short-read coverage without splitting the PacBio sequence. Warning: this will allow more sequences that went through the SMRTportal to not be fixed.\n";
    exit(1);
}

# get grid options
my $submitCommand 	= getGlobal("gridSubmitCommand");
my $nameOption 		= getGlobal("gridNameOption");
my $outputOption 	= getGlobal("gridOutputOption");

#check for valid parameters for requested partitions and threads
my $limit = 1024;
my $FILES_PER_PARTITION = (getGlobal("maxUncorrectedGap") > 0 ? 3 : 1);
my $FILES_PER_THREAD = 3;

my $sysLimit = `ulimit -Sn`;
chomp($sysLimit);
if (defined($sysLimit)) {
   $limit = $sysLimit;
}
if ($limit - $MIN_FILES_WITHOUT_PARTITIONS <= $FILES_PER_PARTITION*getGlobal("partitions") || $limit - $MIN_FILES_WITHOUT_PARTITIONS <= getGlobal("threads") * $FILES_PER_THREAD) {
	my $maxPartitions = floor(($limit - $MIN_FILES_WITHOUT_PARTITIONS) / $FILES_PER_PARTITION);
	my $maxThreads = (floor($limit - $MIN_FILES_WITHOUT_PARTITIONS) / $FILES_PER_THREAD);
   setGlobal("partitions", ($maxPartitions < getGlobal("partitions") ? $maxPartitions : getGlobal("partitions")));
   setGlobal("threads", ($maxThreads < getGlobal("threads") ? $maxThreads : getGlobal("threads")));
   print STDERR "Warning: file handle limit of $limit prevents using requested partitions. Reset partitions to " . getGlobal("partitions") . " and " . getGlobal("threads") . " threads. If you want more partitions, reset the limit using ulimit -Sn and try again.\n";
}
if (getGlobal("threads") > getGlobal("partitions")) {
   setGlobal("threads", getGlobal("partitions") - 1);
   print STDERR "Warning: number of partitions should be > # threads. Adjusted threads to be ". getGlobal("threads") . "\n";
}

print STDERR "Running with " . getGlobal("threads") . " threads and " . getGlobal("partitions") . " partitions\n";

my $CA = getBinDirectory();
my $AMOS = "$CA/../../../AMOS/bin/";
my $BLASR = "$CA/../../../smrtanalysis/analysis/bin/";
my $BOWTIE = "$CA/../../../bowtie2/";
my $wrk = makeAbsolute("");
my $asm = "asm";
my $scriptParams = getGlobal("sgeScript");

if (defined($scriptParams)) {
   if (!defined(getGlobal("sgeCorrection"))) {
   	setGlobal("sgeCorrection", $scriptParams);
   }
}

my $useGrid = getGlobal("useGrid");
if (defined($useGrid)) {
	setGlobal("submitToGrid", $useGrid);
}
elsif (defined(getGlobal("sge"))) {
   setGlobal("submitToGrid", 1);
}

my $caCNS  = getGlobal("cnsConcurrency");
if (defined($caCNS)) {
   setGlobal("consensusConcurrency", $caCNS);
}

if (! -e "$AMOS/bank-transact") {
   # try to use path
   my $amosPath = `which bank-transact`;
   chomp $amosPath;
   my @t = split '/', "$amosPath";
   pop @t;                      #  bank-transact 
   $AMOS = join '/', @t;  #  path to the assembler

   # if we really can't find it just give up
   if (! -e "$AMOS/bank-transact") {
      die "AMOS binaries: bank-transact not found in $AMOS\n";
   }
}

if (! -e "$BLASR/blasr") {
   # try to use path
   my $amosPath = `which blasr`;
   chomp $amosPath;
   my @t = split '/', "$amosPath";
   pop @t;                      #  blasr 
   $BLASR = join '/', @t;  #  path to the assembler

   # if we really can't find it just give up
   if (! -e "$BLASR/blasr") {
      die "BLASR binaries: blasr not found in $BLASR\n" if defined(getGlobal("blasr"));
   }
}

my $blasrVersion = 1.2;
if (-e "$BLASR/blasr") {
   $blasrVersion = `$BLASR/blasr -version |tail -n 1`;
   chomp $blasrVersion;
}

if (! -e "$BOWTIE/bowtie2-build") {
   # try to use path
   my $amosPath = `which bowtie2-build`;
   chomp $amosPath;
   my @t = split '/', "$amosPath";
   pop @t;                      #  bowtie2-build 
   $BOWTIE = join '/', @t;  #  path to the assembler

   # if we really can't find it just give up
   if (! -e "$BOWTIE/bowtie2-build") {
      die "Bowtie2 binaries: bowtie not found in $BOWTIE\n" if defined(getGlobal("bowtie"));
   }
}

print STDERR "********* Starting correction...\n CA: $CA\nAMOS:$AMOS\nSMRTportal:$BLASR ($blasrVersion)\nBowtie:$BOWTIE\n";
print STDERR "******** Configuration Summary ********\n";
foreach my $key (keys %global) {
	if (length($key) < $TAB_SIZE) {
		print STDERR "$key\t\t\t=";
	} elsif (length($key) < 2*$TAB_SIZE) {
		print STDERR "$key\t\t="
	} else {
		print STDERR "$key\t=";
	}
	print STDERR "\t" . $global{$key} . "\n";
}
print STDERR "****************\n";

# get global variable values
my $specFile = getGlobal("specFile");
my $libraryname = getGlobal("libraryname");
my $fastqFile = getGlobal("fastqFile");

my $ovlErrorRate = getGlobal("utgErrorRate");
my $ovlErrorLimit = getGlobal("utgErrorLimit");

my $maxUncorrectedGap = getGlobal("maxUncorrectedGap");
my $coverage = getGlobal("coverageCutoff");
my $genomeSize = getGlobal("genomeSize");
my $shortReads = getGlobal("shortReads");
my $longReads = getGlobal("longReads");
my $length = getGlobal("length");
my $repeats = getGlobal("repeats");

my $threads = getGlobal("threads");
my $partitions = getGlobal("partitions");

my $submitToGrid = getGlobal("submitToGrid");
my $sge = getGlobal("sge");
my $sgeOvl = getGlobal("sgeOverlap");
my $sgeCorrection = getGlobal("sgeCorrection");
my $sgeConsensus = getGlobal("sgeConsensus");
my $sgeTaskID = getGlobal("gridTaskID");

my $consensusConcurrency = getGlobal("consensusConcurrency");

my $cleanup = getGlobal("cleanup");

my $cmd = "";

if (! -e "$wrk/temp$libraryname") {
   runCommand("$wrk", "mkdir temp$libraryname");
   # generate the ca spec file, since we support additional options in the spec file, we need to strip those out before passing it to ca
   updateSpecFile("$wrk/temp$libraryname/$libraryname.spec");
}
$specFile = "$wrk/temp$libraryname/$libraryname.spec";

if (! -e "$wrk/temp$libraryname/$libraryname.frg") {
   runCommand($wrk, "$CA/fastqToCA -libraryname $libraryname -type sanger -innie -technology pacbio-long -reads " . makeAbsolute($fastqFile) . " > $wrk/temp$libraryname/$libraryname.frg"); 
}

# now that were ready, add the frg file info to the command line args
my $commandLineOptionsNoFrgs = "-s $specFile $commandLineOptions";
$commandLineOptions = "-s $specFile $commandLineOptions @fragFiles";

# and we're off
submitScript($wrk, $asm, undef) if (!runningOnGrid());
my $sgeName = (defined(getGlobal("sgeName")) ? "_" . getGlobal("sgeName") : "");
if (! -d "$wrk/temp$libraryname/$asm.gkpStore") {
   $cmd = "$CA/runCA ";
   $cmd .= " -s $specFile -p $asm -d temp$libraryname ";
   $cmd .= " stopAfter=initialStoreBuilding ";
   $cmd .= " sgeName=\"" . getGlobal("sgeName") . "\" " if defined(getGlobal("sgeName"));
   $cmd .= " @fragFiles $wrk/temp$libraryname/$libraryname.frg";
   runCommand($wrk, $cmd);
}

# make assumption that we correct using all libraries preceeding pacbio
# figure out what number of libs we have and what lib is pacbio
my $numLib = `$CA/gatekeeper -dumpinfo temp$libraryname/$asm.gkpStore | grep LIB |awk '{print \$1}'`;
chomp($numLib);

my $minCorrectLib = 0;
my $maxCorrectLib = 0;
my $libToCorrect = 0;
for (my $i = 1; $i <= $numLib; $i++) {
   if (system("$CA/gatekeeper -isfeatureset $i doConsensusCorrection temp$libraryname/$asm.gkpStore") == 0) {
   	  if ($libToCorrect != 0) {
   	  	die("Error: only one PacBio library can be corrected. Both libraries $libToCorrect and $i are set to be corrected. Please double-check your input files and try again", undef);
   	  }
      $libToCorrect = $i;
    } else {
      if ($minCorrectLib == 0) { $minCorrectLib = $i; }
      $maxCorrectLib = $i;
   }
}

# check that we were able to find the libraries for correction as expected
if ($libToCorrect <= $minCorrectLib) {
	die("Error: The PacBio library $libToCorrect must be the last library loaded but it preceedes $minCorrectLib. Please double-check your input files and try again.", undef);
}
if (defined(getGlobal("longReads")) && getGlobal("longReads") == 1) {
   $maxCorrectLib = $libToCorrect;
   $minCorrectLib = $libToCorrect if ($minCorrectLib == 0);
} elsif (($libToCorrect == 0 || ($minCorrectLib == 0 && $maxCorrectLib == 0))) {
   die ("Error: unable to find a library to correct. Please double-check your input files and try again.", undef);
}
print STDERR "Will be correcting PacBio library $libToCorrect with librarie[s] $minCorrectLib - $maxCorrectLib\n";
my $totalBP = 0;
my $totalNumToCorrect = 0;
# compute the number of bases in the gateeeker to be corrected
open(F, "$CA/gatekeeper -dumpinfo temp$libraryname/$asm.gkpStore |") or die("Couldn't open gatekeeper store", undef);

while (<F>) {
   s/^\s+//;
   s/\s+$//;

   my @array = split '\s+';
   if ($#array == 8 && $array[0] == $libToCorrect) { 
      $totalBP = $array[6];
   }
}
close(F); 

if (! -e "$wrk/temp$libraryname/$asm.toerase.out") {
   # here is where we filter for specified length as well as max of longest X of coverage for correction
   # use the genome size/coverage, if available to subset the sequences
   if ($genomeSize != 0 && getGlobal("maxCoverage") != 0) {
      $totalBP = $genomeSize * getGlobal("maxCoverage"); 
   }
   runCommand($wrk, "$CA/gatekeeper -dumpfragments -invert -tabular -longestovermin $libToCorrect $length -longestlength $libToCorrect $totalBP temp$libraryname/$asm.gkpStore |awk '{if (!(match(\$1, \"UID\") != 0 && length(\$1) == " . length("UID") . ")) { print \"frg uid \"\$1\" isdeleted 1\"; } }' > $wrk/temp$libraryname/$asm.toerase.uid");
   runCommand($wrk, "$CA/gatekeeper --edit $wrk/temp$libraryname/$asm.toerase.uid $wrk/temp$libraryname/$asm.gkpStore > $wrk/temp$libraryname/$asm.toerase.out 2> $wrk/temp$libraryname/$asm.toerase.err");
}

# compute the number of bases left after our filtering gateeeker to be corrected
my $totalCorrectingWith = 0;
open(F, "$CA/gatekeeper -dumpinfo temp$libraryname/$asm.gkpStore |") or die("Couldn't open gatekeeper store", undef);
while (<F>) {
   s/^\s+//;
   s/\s+$//;

   my @array = split '\s+';
   if ($#array == 8) {
   	  if ($array[0] == $libToCorrect) { 
      	$totalBP = $array[6];
      	$totalNumToCorrect = $array[3] + $array[4];
   	  } 
   	  if ($minCorrectLib <= $array[0] && $array[0] <= $maxCorrectLib) {
   	  	$totalCorrectingWith += $array[6];
   	  }
   }
}
close(F);

# check that we have good data
if ($genomeSize != 0) {
    print STDERR "Running with " . ($totalBP / $genomeSize) . "X (for genome size $genomeSize) of $libraryname sequences ($totalBP bp).\n";
    print STDERR "Correcting with " . floor($totalCorrectingWith / $genomeSize) . "X sequences ($totalCorrectingWith bp).\n";
} else {
    print STDERR "Running with $totalBP bp for $libraryname.\n";
    print STDERR "Correcting with $totalCorrectingWith bp.\n";
}
    
if ($totalBP == 0) {
	print STDERR "Error: All $libraryname sequences were eliminated. Please check the length threshold of $length and your input file $fastqFile.\n";
    runCommand("$wrk", "rm -rf temp$libraryname");
	die;
}
if ($totalCorrectingWith == 0) {
	print STDERR "Error: No high-accuracy sequences for correction. Please check your input FRG files " . join(", ", @fragFiles) . "\n";
	runCommand("$wrk", "rm -rf temp$libraryname");
	die;
}
if ($genomeSize != 0 && floor($totalCorrectingWith / $genomeSize) > $MAX_CORRECTION_COVERAGE) {
	print STDERR "Warning: input a total of " . floor($totalCorrectingWith / $genomeSize) . " of high-accuracy coverage for correction. For best performance, at most $MAX_CORRECTION_COVERAGE is recommended.\n";
	# could randomly subsample here
}

if (!defined(getGlobal("bowtie"))) {
   # run correction up thorough meryl
   if (! -e "$wrk/temp$libraryname/0-mercounts/$asm.nmers.ovl.fasta") {
      $cmd  = "$CA/runCA ";
      $cmd .=    "-s $specFile ";
      $cmd .=    "-p $asm -d temp$libraryname ";
      $cmd .=    "ovlHashLibrary=$libToCorrect ";
      $cmd .=    "ovlRefLibrary=$minCorrectLib-$maxCorrectLib ";
      $cmd .=    "ovlCheckLibrary=1 ";
      $cmd .=    "obtHashLibrary=$minCorrectLib-$maxCorrectLib ";
      $cmd .=    "obtRefLibrary=$minCorrectLib-$maxCorrectLib ";
      $cmd .=    "obtCheckLibrary=0 ";
      $cmd .=    "sgeName=\"" . getGlobal("sgeName") . "\" " if defined(getGlobal("sgeName"));
      $cmd .=    "doOverlapBasedTrimming=0 stopAfter=meryl";
      runCommand($wrk, $cmd);
   }
   
   # set the meryl threshold based on the genome size and coverage (if specified)
   if (!defined(getGlobal("ovlMerThreshold")) || getGlobal("ovlMerThreshold") == "auto") {
      # no threshold specified, check if the chosen one is OK
      my $autoSetThreshold = `cat temp$libraryname/0-mercounts/*estMerThresh.out`;
      chomp($autoSetThreshold);
      setGlobal("ovlMerThreshold", $autoSetThreshold);
      
      my $corrCov = (defined($genomeSize) && $genomeSize != 0 ? floor($totalCorrectingWith / $genomeSize) : 0);
      my $pacCov = (defined($genomeSize) && $genomeSize != 0 ? floor($totalBP / $genomeSize) : 0);
      my $maxCov = $pacCov + (defined($longReads) && $longReads == 1 ? 0 : $corrCov);
      
      if ($autoSetThreshold < ($maxCov * $REPEAT_MULTIPLIER)) {
         setGlobal("ovlMerThreshold", $maxCov * $REPEAT_MULTIPLIER);
         unlink("$wrk/temp$libraryname/0-mercounts/$asm.nmers.ovl.fasta");
         print STDERR "Resetting from auto threshold of $autoSetThreshold to be " . ($maxCov * 10) . "\n";
      }
   }
}

if (! -d "$wrk/temp$libraryname/$asm.ovlStore") {
   # run trimming if needed first
   if (getGlobal("doOverlapBasedTrimming") != 0 && ! -e "$wrk/temp$libraryname/0-overlaptrim/overlaptrim.success") {
      if (-e "$wrk/temp$libraryname/0-overlaptrim-overlap/overlap.sh") {
         die("Overlap trimming failed. Remove the temporary directory temp$libraryname and try again.", undef);   
      }
      $cmd  =    "-s $specFile ";
      $cmd .=    "-p $asm -d . ";
      $cmd .=    "ovlHashLibrary=$libToCorrect ";
      $cmd .=    "ovlRefLibrary=$minCorrectLib-$maxCorrectLib ";
      $cmd .=    "ovlCheckLibrary=1 ";
      $cmd .=    "obtHashLibrary=$minCorrectLib-$maxCorrectLib ";
      $cmd .=    "obtRefLibrary=$minCorrectLib-$maxCorrectLib ";
      $cmd .=    "obtCheckLibrary=0 ";
      $cmd .=    "sgeName=\"" . getGlobal("sgeName") . "\" " if defined(getGlobal("sgeName"));
      $cmd .=    "sgePropagateHold=\"pBcR_$asm$sgeName\" "; 
      $cmd .=    "stopAfter=overlapBasedTrimming";
      if ($submitToGrid == 1) {
         submitRunCA("$wrk/temp$libraryname", $asm, $cmd, "runCA_obt_$asm$sgeName");
      } else {
         runCommand("$wrk/temp$libraryname", "$CA/runCA $cmd");
      }
   }
   
   # now run the correction
   my $ovlThreshold = getGlobal("ovlMerThreshold");

   if (defined(getGlobal("blasr")) && $totalBP > $MAX_BLASR_LIMIT) {
       print STDERR "Warning: Blasr is currently limited to $MAX_BLASR_LIMIT. Switching to default aligner.\n";
       setGlobal("blasr", undef);
   }
   if (defined($longReads) && $longReads == 1) {
      if ($totalBP > $MAX_BLASR_LIMIT) {
          print STDERR "Warning: Cannot use blasr due to reference size limitation. Turning off long read option.\n";
          $longReads = undef;
      } else {
         print STDERR "Warning: Blasr is required to align long reads to long reads. Switching blasr ON.\n";
         setGlobal("blasr", defined(getGlobal("blasr")) ? getGlobal("blasr") . " -minSubreadLength 500 -minReadLength 500 -maxLCPLength 15 -minPctIdentity 70.0" : "-minSubreadLength 500 -minReadLength 500 -advanceHalf -noRefineAlign -ignoreQuality -minMatch 10 -maxLCPLength 15 -minPctIdentity 70.0");
      }
   }

   # when we were asked to not use CA's overlapper, first perform common options
   if (defined(getGlobal("samFOFN")) || defined(getGlobal("blasr")) || defined(getGlobal("bowtie"))) {
      # check for sam tools
      if (!-e "$CA/convertSamToCA") { 
         die("Error: request to use Blasr, Bowtie, or SAM files requires CA to be built with SAMTOOLS. Please rebuild CA and try again");
      }
      # create directories
      runCommand($wrk, "mkdir $wrk/temp$libraryname/1-overlapper") if (! -d "$wrk/temp$libraryname/1-overlapper");

      # first partition the data so we know how many jobs we have
      my $cmd;

      my $ovlHashBlockLength = getGlobal("ovlHashBlockLength");
      my $ovlHashBlockSize   = $totalNumToCorrect;
      my $ovlRefBlockSize    = getGlobal("ovlRefBlockSize");
      my $ovlRefBlockLength  = getGlobal("ovlRefBlockLength");

      if (($ovlRefBlockSize > 0) && ($ovlRefBlockLength > 0)) {
          die("can't set both ovlRefBlockSize and ovlRefBlockLength", undef);
      }

      $cmd  = "$CA/overlap_partition \\\n";
      $cmd .= " -g  $wrk/temp$libraryname/$asm.gkpStore \\\n";
      #$cmd .= " -bl $totalBP \\\n"; # for now we cant partition the hash with other tools
      $cmd .= " -bs $ovlHashBlockSize \\\n";
      $cmd .= " -rs $ovlRefBlockSize \\\n" if (defined(getGlobal("ovlRefBlockSize")));
      $cmd .= " -rl $ovlRefBlockLength \\\n"  if (defined(getGlobal("ovlRefBlockLength")));
      $cmd .= " -H $libToCorrect \\\n";
      $cmd .= " -R $minCorrectLib-$maxCorrectLib \\\n";
      $cmd .= " -C \\\n";
      $cmd .= " -o  $wrk/temp$libraryname/1-overlapper";
      runCommand($wrk, $cmd);
      open(F, "< $wrk/temp$libraryname/1-overlapper/ovljob") or die("failed partition for overlapper: no ovljob file found", undef);
      my @job = <F>;
      close(F);
    
      my $blasrThreshold = 0;
      my $ovlThreads = getGlobal("ovlThreads");   
      my $numOvlJobs = 0;
      my $blasrOpts = getGlobal("blasr");
      my $bowtieOpts = getGlobal("bowtie");
      my $suffix = "sam";

      # SMRTportal 1.4 changed parameters so update our blasr options
      if ($blasrVersion >= 1.4) {
         $blasrOpts =~ s/-ignoreQuality//g;
      }
   
      if (defined(getGlobal("samFOFN"))) {
         my $fofn = getGlobal("samFOFN");
         $numOvlJobs = `wc -l $fofn`;
         chomp $numOvlJobs;
      } else {
         $numOvlJobs = $#job + 1;
      }
      if (-e "$wrk/temp$libraryname/1-overlapper/overlap.sh") {
         goto checkforerror;
      }
      
      # dump gatekeeper store and print eid to iid, add to it the gkpStore fastq file into eidToIID and lengths
      runCommand("$wrk/temp$libraryname", "$CA/gatekeeper -dumpfragments -tabular $asm.gkpStore |awk '{print \$1\"\\t\"\$2}' > $asm.eidToIID");
      runCommand("$wrk/temp$libraryname", "$CA/gatekeeper -dumpfragments -tabular $asm.gkpStore |awk '{print \$2\"\\t\"\$10}' > $asm.iidToLen");
   
      # now do specific steps (either import SAMs or run blasr)
      if (defined(getGlobal("samFOFN"))) {
         my $fofn = getGlobal("samFOFN");

         # get the external UIDs from the input file rather than gatekeeper since these won't match our externally-generated SAM
         runCommand("$wrk/temp$libraryname", "cat $asm.gkpStore.fastqUIDmap | awk '{if (NF > 3) { print \$NF\"\\t\"\$(NF-1); } print \$3\"\\t\"\$2; }' >> $asm.eidToIID");

         # convert sam to ovl and put it into the overlap directory
         open(F, "< $fofn") or die("Couldn't open '$fofn'", undef);
         while (<F>) {
            s/^\s+//;
            s/\s+$//;
      
            next if (m/^\s*\#/);
            next if (m/^\s*$/);
            if (! -e $_) {
               print STDERR "Warning: could not open SAM file $_, skipping\n";
               next;
            } elsif ($_ =~ /\.sam$/) {
                $suffix = "sam";
            } elsif ($_ =~ /\.sam\.gz$/) {
               $suffix = "sam.gz";
            } elsif ($_ =~ /\.sam\.bz2$/) {
               $suffix = "sam.bz2";
            }
            
            $numOvlJobs++;
            runCommand("$wrk/temp$libraryname/1-overlapper/001", "unlink $wrk/temp$libraryname/1-overlapper/001/$numOvlJobs.$suffix") if (-e "$wrk/temp$libraryname/1-overlapper/001/$numOvlJobs.$suffix");
            runCommand("$wrk/temp$libraryname/1-overlapper/001", "cp $_ $wrk/temp$libraryname/1-overlapper/001/$numOvlJobs.$suffix");
         }
         close(F); 
      } else {
         runCommand("$wrk/temp$libraryname/1-overlapper", "$CA/gatekeeper -dumpfasta long_reads -randomsubset $libToCorrect 1 $wrk/temp$libraryname/$asm.gkpStore");
         for (my $i = $minCorrectLib; $i <= $maxCorrectLib; $i++) {
           runCommand("$wrk/temp$libraryname/1-overlapper", "$CA/gatekeeper -dumpfasta " . $i . "_lib ". (defined($longReads) && $longReads == 1 && $i == $libToCorrect ? " -allreads -allbases" : "") . " -randomsubset $i 1 $wrk/temp$libraryname/$asm.gkpStore");
           runCommand("$wrk/temp$libraryname/1-overlapper", "$CA/gatekeeper -dumpfasta " . $i . "_subset " . (defined($longReads) && $longReads == 1 && $i == $libToCorrect ? " -allreads -allbases" : "") . " -randomsubset $i 0.01 $wrk/temp$libraryname/$asm.gkpStore");
         }
         runCommand("$wrk/temp$libraryname/1-overlapper", "cat `ls [0-9]*_lib.fasta` > correct_reads.fasta");
         runCommand("$wrk/temp$libraryname/1-overlapper", "cat `ls [0-9]*_subset.fasta` > correct_subset.fasta");
         runCommand("$wrk/temp$libraryname/1-overlapper", "rm `ls [0-9]*_subset.fasta*`");
         runCommand("$wrk/temp$libraryname/1-overlapper", "rm `ls [0-9]*_lib.fasta*`");
       
         if (defined(getGlobal("blasr"))) {
            # run with best limit set to kmer with a small (0.5-1% of the data)   
            my $ovlThreads = getGlobal("ovlThreads");
            runCommand("$wrk/temp$libraryname/1-overlapper", "$BLASR/sawriter long.sa long_reads.fasta");
            if (-e "$BLASR/../../etc/setup.sh") {
               runCommand("$wrk/temp$libraryname/1-overlapper", "source $BLASR/../../etc/setup.sh && $BLASR/blasr -sa long.sa correct_subset.fasta long_reads.fasta -nproc $ovlThreads -bestn $ovlThreshold $blasrOpts -sam -out subset.sam");
            } else {
               runCommand("$wrk/temp$libraryname/1-overlapper", "$BLASR/blasr -sa long.sa correct_subset.fasta long_reads.fasta -nproc $ovlThreads -bestn $ovlThreshold $blasrOpts -sam -out subset.sam");
            }
            runCommand("$wrk/temp$libraryname/1-overlapper", "cat subset.sam |grep -v \"@\" | awk '{print \$1}' |sort -T . |uniq -c |awk '{print \$1}' > subset.hist");
         
            # pick threshold as 2sd (95% of the bases)
            $blasrThreshold = pickMappingThreshold("$wrk/temp$libraryname/1-overlapper/subset.hist", ($blasrOpts =~ m/maxLCPLength/));
            runCommand("$wrk/temp$libraryname/1-overlapper", "unlink subset.sam");
            runCommand("$wrk/temp$libraryname/1-overlapper", "unlink correct_subset.fasta");
         } elsif (defined(getGlobal("bowtie"))) {     
         # run with best limit set to all with a small (0.5-1% of the data)   
            my $ovlThreads = getGlobal("ovlThreads");
            runCommand("$wrk/temp$libraryname/1-overlapper", "$BOWTIE/bowtie2-build -f long_reads.fasta long.sa");
            runCommand("$wrk/temp$libraryname/1-overlapper", "$BOWTIE/bowtie2 -x long.sa -f correct_subset.fasta -p $ovlThreads -a $bowtieOpts -S subset.sam");
            runCommand("$wrk/temp$libraryname/1-overlapper", "cat subset.sam |grep -v \"@\" | awk '{print \$1}' |sort -T . |uniq -c |awk '{print \$1}' > subset.hist");
         
            # pick threshold as 2sd (95% of the bases)
            $blasrThreshold = pickMappingThreshold("$wrk/temp$libraryname/1-overlapper/subset.hist", 0);
            runCommand("$wrk/temp$libraryname/1-overlapper", "unlink subset.sam");
            runCommand("$wrk/temp$libraryname/1-overlapper", "unlink correct_subset.fasta");
         }
      }
   
      # now that we have our parameters, create a run job
      open F, "> $wrk/temp$libraryname/1-overlapper/overlap.sh" or die ("can't open '$wrk/temp$libraryname/1-overlapper/overlap.sh'");
      print F "#!" . getGlobal("shell") ."\n";
      print F getBinDirectoryShellCode();
      print F "\n";
      print F "jobid=\$$sgeTaskID\n";
      print F "if [ x\$jobid = x -o x\$jobid = xundefined -o x\$jobid = x0 ]; then\n";
      print F "jobid=\$1\n";
      print F "fi\n";
      print F "\n";
      print F "if test x\$jobid = x; then\n";
      print F "  echo Error: I need $sgeTaskID set, or a job index on the command line\n";
      print F "  exit 1\n";
      print F "fi\n";
      print F "\n";
      print F "bat=`head -n \$jobid $wrk/temp$libraryname/1-overlapper/ovlbat | tail -n 1`\n";
      print F "job=`head -n \$jobid $wrk/temp$libraryname/1-overlapper/ovljob | tail -n 1`\n";
      print F "opt=`head -n \$jobid $wrk/temp$libraryname/1-overlapper/ovlopt | tail -n 1`\n";
      print F "\n";
      print F "if [ ! -d $wrk/temp$libraryname/1-overlapper/\$bat ]; then\n";
      print F "  mkdir $wrk/temp$libraryname/1-overlapper/\$bat\n";
      print F "fi\n";
      print F "\n";
      print F "if [ -e $wrk/temp$libraryname/1-overlapper/\$bat/\$job.ovb ]; then\n";
      print F "  echo Job previously completed successfully.\n";
      print F "  exit\n";
      print F "fi\n";
      print F "\n";
      print F "if [ x\$bat = x ]; then\n";
      print F "  echo Error: Job index out of range.\n";
      print F "  exit 1\n";
      print F "fi\n";
      print F "\n";
      print F "options=(\$opt)\n";
      print F "params=(\$(echo \${options[3]}  | tr \"-\" \"\\n\")) \n";
      print F "zeroJobId=\$((\$jobid - 1))\n";
      print F "numOvlJobs=$numOvlJobs\n";
      print F "start=\${params[0]} \n";
      print F "start=\$((\$start - 1))\n";
      print F "end=\${params[1]} \n";
      print F "total=\$((\$end - \$start))\n";
      print F "echo \"Running partition \$job with options \$opt start \$start end \$end total \$total zero job \$zeroJobId and stride \$numOvlJobs\"\n";
      print F "\n";
      if (defined(getGlobal("blasr"))) {
         if (-e "$BLASR/../../etc/setup.sh") {
            print F "   source $BLASR/../../etc/setup.sh \n";
         }
         print F "   $BLASR/blasr \\\n";
         print F "          -start \$zeroJobId -stride \$numOvlJobs \\\n";
         print F "          -sa $wrk/temp$libraryname/1-overlapper/long.sa $wrk/temp$libraryname/1-overlapper/correct_reads.fasta \\\n";
         print F "          $wrk/temp$libraryname/1-overlapper/long_reads.fasta \\\n";
         print F "          -nproc $ovlThreads -bestn $blasrThreshold \\\n";
         print F "          $blasrOpts \\\n";
         print F "          -sam -out $wrk/temp$libraryname/1-overlapper/\$bat/\$job.sam \\\n";
         print F "          > $wrk/temp$libraryname/1-overlapper/\$jobid.out 2>&1 \\\n";
         print F "            && touch $wrk/temp$libraryname/1-overlapper/\$jobid.blasr.success\n";
         print F "   if test -e $wrk/temp$libraryname/1-overlapper/\$jobid.blasr.success ; then\n";
         print F "      echo Blasr completed.\n";
         print F "   else\n";
         print F "      echo Blasr failed.\n";
         print F "      tail $wrk/temp$libraryname/1-overlapper/\$jobid.out && exit\n";
         print F "   fi\n";
      } elsif (defined(getGlobal("bowtie"))) {
         print F "   $BOWTIE/bowtie2 \\\n";
         print F "          -s \$start -u \$total \\\n";
         print F "          -x $wrk/temp$libraryname/1-overlapper/long.sa -f $wrk/temp$libraryname/1-overlapper/correct_reads.fasta \\\n";
         print F "          -p $ovlThreads -k $blasrThreshold \\\n";
         print F "          $bowtieOpts \\\n";
         print F "          -S $wrk/temp$libraryname/1-overlapper/\$bat/\$job.sam \\\n";
         print F "          > $wrk/temp$libraryname/1-overlapper/\$jobid.out 2>&1 \\\n";
         print F "            && touch $wrk/temp$libraryname/1-overlapper/\$jobid.blasr.success\n";
         print F "   if test -e $wrk/temp$libraryname/1-overlapper/\$jobid.blasr.success ; then\n";
         print F "      echo Bowtie completed.\n";
         print F "   else\n";
         print F "      echo Bowtie failed.\n";
         print F "      tail $wrk/temp$libraryname/1-overlapper/\$jobid.out && exit\n";
         print F "   fi\n";
      }
      print F "   \$bin/convertSamToCA \\\n";
      print F "        $wrk/temp$libraryname/1-overlapper/\$bat/\$job.$suffix $wrk/temp$libraryname/$asm.eidToIID $wrk/temp$libraryname/$asm.iidToLen \\\n";
      print F "        > $wrk/temp$libraryname/1-overlapper/\$bat/\$job.ovls 2> $wrk/temp$libraryname/1-overlapper/\$jobid.java.err\\\n";
      print F "            && touch $wrk/temp$libraryname/1-overlapper/\$jobid.java.success\n";
      print F "   if test -e $wrk/temp$libraryname/1-overlapper/\$jobid.java.success ; then\n";
      print F "      echo SamToCA conversion completed.\n";
      print F "   else\n";
      print F "      echo SamToCA conversion failed.\n";
      print F "      tail $wrk/temp$libraryname/1-overlapper/\$jobid.java.err && exit\n";
      print F "   fi\n";
      print F "   \$bin/convertOverlap -ovl < $wrk/temp$libraryname/1-overlapper/\$bat/\$job.ovls > $wrk/temp$libraryname/1-overlapper/\$bat/\$job.ovb\n";
      close(F);
      chmod 0755, "$wrk/temp$libraryname/1-overlapper/overlap.sh";
      
      if ($submitToGrid == 1) {
          my $sgeName = "pBcR_ovl_$asm$sgeName";
      	  my $jobName = getGridArrayName($sgeName, $numOvlJobs);
      	  my $arrayOpt = getGridArrayOption($sgeName, $numOvlJobs);
         submitBatchJobs($wrk, $asm, "$submitCommand $sge $sgeOvl $nameOption \"$jobName\" $arrayOpt $outputOption /dev/null $wrk/temp$libraryname/1-overlapper/overlap.sh", $jobName);
      } else {
         for (my $i = 1; $i <= $numOvlJobs; $i++) {
            schedulerSubmit("$wrk/temp$libraryname/1-overlapper/overlap.sh $i");
         }
         schedulerSetNumberOfProcesses(getGlobal("ovlConcurrency"));
         schedulerFinish();
      } 
      
  checkforerror:
      my $failedJobs = 0;
      my $failureMessage = "";

      open(B, "< $wrk/temp$libraryname/1-overlapper/ovlbat") or die("failed to open '$wrk/temp$libraryname/1-overlapper/ovlbat'", undef);
      open(J, "< $wrk/temp$libraryname/1-overlapper/ovljob") or die("failed to open '$wrk/temp$libraryname/1-overlapper/ovljob'", undef);

      while (!eof(B) && !eof(J)) {
         my $b = <B>;  chomp $b;
         my $j = <J>;  chomp $j;

         if ((! -e "$wrk/temp$libraryname/1-overlapper/$b/$j.ovb.gz") &&
             (! -e "$wrk/temp$libraryname/1-overlapper/$b/$j.ovb")) {
            $failureMessage .= "ERROR:  Overlap job $wrk/temp$libraryname/1-overlapper/$b/$j FAILED.\n";
            $failedJobs++;
         } else {
            runCommand("$wrk/temp$libraryname/1-overlapper/$b/", "rm $j.sam");
            runCommand("$wrk/temp$libraryname/1-overlapper/$b/", "rm $j.ovls");
         }
      }
      if (!eof(B) || !eof(J)) {
         print STDERR "Partitioning error; '$wrk/temp$libraryname/1-overlapper/ovlbat' and '$wrk/temp$libraryname/1-overlapper/ovljob' have extra lines.\n";
      }
      close(B);
      close(J);
      
      my $errorType = "sam conversion";
      if (defined(getGlobal("blasr"))) { $errorType = "blasr mapping"; }
      if (defined(getGlobal("bowtie"))) { $errorType = "bowtie mapping"; }
      $failureMessage .= "\n$failedJobs $errorType jobs failed.";
      if ($failedJobs > 0) {
          die "$failureMessage\n";
      }  
   }

   $cmd  =    "-s $specFile ";
   $cmd .=    "-p $asm -d . ";
   $cmd .=    "ovlMerThreshold=$ovlThreshold ";
   $cmd .=    "ovlHashLibrary=$libToCorrect ";
   $cmd .=    "ovlRefLibrary=$minCorrectLib-$maxCorrectLib ";
   $cmd .=    "ovlCheckLibrary=1 ";
   $cmd .=    "obtHashLibrary=$minCorrectLib-$maxCorrectLib ";
   $cmd .=    "obtRefLibrary=$minCorrectLib-$maxCorrectLib ";
   $cmd .=    "obtCheckLibrary=0 ";
   $cmd .=    "sgeName=\"" . getGlobal("sgeName") . "\" " if defined(getGlobal("sgeName"));
   $cmd .=    "sgePropagateHold=\"pBcR_$asm$sgeName\" "; 
   $cmd .=    "stopAfter=overlapper";
   if ($submitToGrid == 1) {
      submitRunCA("$wrk/temp$libraryname", $asm, $cmd, "runCA_ovl_$asm$sgeName");      
   } else {
      runCommand("$wrk/temp$libraryname", "$CA/runCA $cmd");
   }
}

if (! -e "$wrk/temp$libraryname/$asm.layout.success") {
   open F, "> $wrk/temp$libraryname/runCorrection.sh" or die ("can't open '$wrk/temp$libraryname/runCorrection.sh'");
   print F "#!" . getGlobal("shell") ."\n";
   print F getBinDirectoryShellCode();
   print F "\n";
   print F " if test -e $wrk/temp$libraryname/$asm.layout.success; then\n";
   print F "    echo Job previously completed successfully.\n";
   print F " else\n";
   print F "   \$bin/correctPacBio \\\n";
   if (defined($longReads) && $longReads == 1) {
   print F "      -L \\\n";
   }
   print F "      -C $coverage \\\n";
   print F "      -M $maxUncorrectedGap \\\n";
   print F "      -t $threads \\\n";
   print F "       -p $partitions \\\n";
   print F "       -o $asm \\\n";
   print F "       -l $length \\\n";
   print F "        $repeats \\\n";
   print F "        -O $wrk/temp$libraryname/$asm.ovlStore \\\n";
   print F "        -G $wrk/temp$libraryname/$asm.gkpStore \\\n";
   print F "        -e $ovlErrorRate -c $ovlErrorRate  -E $ovlErrorLimit > $wrk/temp$libraryname/$asm.layout.err 2> $wrk/temp$libraryname/$asm.layout.err && touch $wrk/temp$libraryname/$asm.layout.success\n";
   print F " fi\n";
   close(F);
   chmod 0755, "$wrk/temp$libraryname/runCorrection.sh";

   if ($submitToGrid == 1) {    
      my $sgeName = "pBcR_correct_$asm$sgeName";
      submitBatchJobs("$wrk/temp$libraryname", $asm, "$submitCommand $sge $sgeCorrection $nameOption \"$sgeName\" $outputOption /dev/null $wrk/temp$libraryname/runCorrection.sh", $sgeName);
   } else {
      runCommand("$wrk/temp$libraryname", "$wrk/temp$libraryname/runCorrection.sh");
   }
}

if (! -e "$wrk/temp$libraryname/runPartition.sh") {
   open F, "> $wrk/temp$libraryname/runPartition.sh" or die ("can't open '$wrk/temp$libraryname/runPartition.sh'");
   print F "#!" . getGlobal("shell") ."\n";
   print F getBinDirectoryShellCode();
   print F "\n";
   print F "jobid=\$$sgeTaskID\n";
   print F "if [ x\$jobid = x -o x\$jobid = xundefined -o x\$jobid = x0 ]; then\n";
   print F "jobid=\$1\n";
   print F "fi\n";
   print F "\n";
   print F "if test x\$jobid = x; then\n";
   print F "  echo Error: I need $sgeTaskID set, or a job index on the command line\n";
   print F "  exit 1\n";
   print F "fi\n";
   print F "\n";
   print F "if test -e $wrk/temp$libraryname/\$jobid.success ; then\n";
   print F "   echo Job previously completed successfully.\n";
   print F "else\n";
   print F "   numLays=`cat $wrk/temp$libraryname/$asm" . ".\$jobid.lay |grep \"{LAY\" |wc -l`\n";
   print F "   if test \$numLays = 0 ; then\n";
   print F "      touch $wrk/temp$libraryname/\$jobid.fasta\n";
   print F "      touch $wrk/temp$libraryname/\$jobid.qual\n";
   print F "      touch $wrk/temp$libraryname/\$jobid.success\n";
   print F "   else\n";
   print F "      rm -rf $wrk/temp$libraryname/$asm" . ".bnk_partition\$jobid.bnk\n";
   print F "      $AMOS/bank-transact -b $wrk/temp$libraryname/$asm" . ".bnk_partition\$jobid.bnk -m $wrk/temp$libraryname/$asm.\$jobid" . ".lay -c > $wrk/temp$libraryname/bank-transact.\$jobid.err 2>&1\n";
   if (defined($shortReads) && $shortReads == 1) {
   print F "      $AMOS/make-consensus -e 0.03 -w 5 -x $wrk/temp$libraryname/\$jobid.excluded -B -b $wrk/temp$libraryname/" . $asm . ".bnk_partition\$jobid.bnk > $wrk/temp$libraryname/\$jobid.out 2>&1 && touch $wrk/temp$libraryname/\$jobid.success\n";
   } elsif (defined($longReads) && $longReads == 1) {
   print F "      $AMOS/make-consensus -e 0.30 -w 50 -o 5 -x $wrk/temp$libraryname/\$jobid.excluded -B -b $wrk/temp$libraryname/" . $asm . ".bnk_partition\$jobid.bnk > $wrk/temp$libraryname/\$jobid.out 2>&1 && touch $wrk/temp$libraryname/\$jobid.success\n";
   } else {
   print F "      $AMOS/make-consensus -x $wrk/temp$libraryname/\$jobid.excluded -B -b $wrk/temp$libraryname/" . $asm . ".bnk_partition\$jobid.bnk > $wrk/temp$libraryname/\$jobid.out 2>&1 && touch $wrk/temp$libraryname/\$jobid.success\n";
   }
   print F "      if test -e $wrk/temp$libraryname/\$jobid.success ; then\n";
   print F "         $AMOS/bank2fasta -e -q $wrk/temp$libraryname/\$jobid.qual -b $wrk/temp$libraryname/" . $asm . ".bnk_partition\$jobid.bnk > $wrk/temp$libraryname/\$jobid.fasta\n";
   print F "      else\n";
   print F "         rm -rf $wrk/temp$libraryname/$asm" . ".bnk_partition\$jobid.bnk\n";
   print F "         $AMOS/bank-transact -b $wrk/temp$libraryname/$asm" . ".bnk_partition\$jobid.bnk -m $wrk/temp$libraryname/$asm.\$jobid" . ".lay -c > $wrk/temp$libraryname/bank-transact.\$jobid.err 2>&1\n";
   if (defined($shortReads) && $shortReads == 1) {
   print F "         $AMOS/make-consensus -e 0.03 -w 5 -B -b $wrk/temp$libraryname/" . $asm . ".bnk_partition\$jobid.bnk > $wrk/temp$libraryname/\$jobid.out 2>&1 && touch $wrk/temp$libraryname/\$jobid.success\n";
   } elsif (defined($longReads) && $longReads == 1) {
   print F "         $AMOS/make-consensus -e 0.30 -w 50 -o 5 -B -b $wrk/temp$libraryname/" . $asm . ".bnk_partition\$jobid.bnk > $wrk/temp$libraryname/\$jobid.out 2>&1 && touch $wrk/temp$libraryname/\$jobid.success\n";
   } else {
   print F "         $AMOS/make-consensus -B -b $wrk/temp$libraryname/" . $asm . ".bnk_partition\$jobid.bnk > $wrk/temp$libraryname/\$jobid.out 2>&1 && touch $wrk/temp$libraryname/\$jobid.success\n";
   }
   print F "         $AMOS/bank2fasta -e -q $wrk/temp$libraryname/\$jobid.qual -b $wrk/temp$libraryname/" . $asm . ".bnk_partition\$jobid.bnk > $wrk/temp$libraryname/\$jobid.fasta\n";
   print F "      fi\n";
   print F "      if test -e $wrk/temp$libraryname/\$jobid.success ; then\n";
   print F "         rm -rf $wrk/temp$libraryname/$asm" . ".bnk_partition\$jobid.bnk\n";
   print F "         rm $wrk/temp$libraryname/bank-transact.\$jobid.err\n";
   print F "         rm $wrk/temp$libraryname/\$jobid.excluded\n";
   print F "         rm $wrk/temp$libraryname/\$jobid.out\n";
   print F "      fi\n";
   print F "   fi\n";
   print F "fi\n";
   close(F);

   chmod 0755, "$wrk/temp$libraryname/runPartition.sh";

   if ($submitToGrid == 1) {
      my $sgeName = "pBcR_cns_$asm$sgeName";
   	  my $jobName = getGridArrayName($sgeName, $partitions);
   	  my $arrayOpt = getGridArrayOption($sgeName, $partitions);
      submitBatchJobs("$wrk/temp$libraryname", $asm, "$submitCommand $sge $sgeConsensus $nameOption \"$jobName\" $arrayOpt $outputOption /dev/null $wrk/temp$libraryname/runPartition.sh", $jobName);
   } else {
      for (my $i = 1; $i <=$partitions; $i++) {
         schedulerSubmit("$wrk/temp$libraryname/runPartition.sh $i");
      }
      schedulerSetNumberOfProcesses($consensusConcurrency);
      schedulerFinish();
   } 
}

for (my $i = 1; $i <= $partitions; $i++) {
  if (! -e "$wrk/temp$libraryname/$i.success") {
    die "Failed to run correction job $i. Remove $wrk/temp$libraryname/runPartition.sh to try again.\n";
  }
}
runCommand("$wrk/temp$libraryname", "rm -rf core.*");

# combine partitions
runCommand("$wrk/temp$libraryname", "cat `ls $asm.[0-9]*.log | sort -T . -rnk1` > corrected.log");
runCommand("$wrk/temp$libraryname", "cat `ls [0-9]*.fasta |sort -T . -rnk1` > corrected.fasta");
runCommand("$wrk/temp$libraryname", "cat `ls [0-9]*.qual |sort -T . -rnk1` > corrected.qual");

# check how much we split the sequences up
my %splitUIDs = ();
runCommand("$wrk", "$CA/convert-fasta-to-v2.pl -pacbio -s $wrk/temp$libraryname/corrected.fasta -q $wrk/temp$libraryname/corrected.qual -l unsplit_$libraryname > $wrk/temp$libraryname/$libraryname.frg");
runCommand("$wrk/temp$libraryname", "cat corrected.log |awk '{print \$1}' |sort -T . |uniq -c |awk '{if (\$1 > 1) print \$2}' > $asm.split.uid");

# read in the split UIDs so we can process them
open(F, "< $wrk/temp$libraryname/$asm.split.uid") or die("Couldn't open '$wrk/temp$libraryname/$asm.split.uid'", undef);
while (<F>) {
   s/^\s+//;
   s/\s+$//;
   chomp $_;

   $splitUIDs{$_} = 1;
}
close(F);

open (F, "< $wrk/temp$libraryname/corrected.log") or die ("Couldn't open '$wrk/temp$libraryname/corrected.log", undef);
open (O, "> $wrk/temp$libraryname/$asm.split.allEdit") or die ("Couldn't open '$wrk/temp$libraryname/$asm.split.allEdit", undef);

my $totalSplitBases = 0;
while (<F>) {
   s/^\s+//;
   s/\s+$//;
   chomp $_;

   my ($uid, $newname, $subread, $start, $end) = split '\s+', $_;
   my $len = $end - $start;
   if (defined($splitUIDs{$uid}) && $splitUIDs{$uid} == 1) {
      print O "$newname\n";

      $totalSplitBases += $len;
   }
}
close(F);
close(O);

print STDOUT "Total split bases is $totalSplitBases vs $totalBP so ratio is " . ($totalSplitBases/$totalBP) . "\n";
if (($totalSplitBases / $totalBP) > $MAX_SPLIT_PERCENTAGE && defined(getGlobal("maxCorrectionRounds")) && getGlobal("maxCorrectionRounds") > 1) { 
   if (! -e "$wrk/temp$libraryname/temproundTwo/corrected.fasta") {
      runCommand("$wrk", "$CA/gatekeeper -dumpfastq $wrk/temp$libraryname/$asm.split -iid $wrk/temp$libraryname/$asm.split.uid $wrk/temp$libraryname/$asm.gkpStore > $wrk/temp$libraryname/$asm.split.out 2> $wrk/temp$libraryname/$asm.split.err");
      runCommand("$wrk", "$CA/gatekeeper -T -o $wrk/temp$libraryname/$asm.corrected.gkpStore $wrk/temp$libraryname/$libraryname.frg > $wrk/temp$libraryname/$asm.corrected.out 2> $wrk/temp$libraryname/$asm.corrected.err");
      runCommand("$wrk", "$CA/gatekeeper -dumpfragments -tabular $wrk/temp$libraryname/$asm.corrected.gkpStore |grep -v UID > $wrk/temp$libraryname/corrected.uids");

      my %gkpSeqs = ();
      open(F, "< $wrk/temp$libraryname/corrected.uids") or die("Couldn't open '$wrk/temp$libraryname/corrected.uids'", undef);
      while (<F>) {
         s/^\s+//;
         s/\s+$//;
         chomp $_;
         my ($uid, $iid, $mate, $mateUID, $libUID, $libIID, $deleted, $random, $orient, $length, $clrBgn, $clrEnd) = split '\s+', $_;
         $gkpSeqs{$uid} = 1;
      }
      close(F);

      open (F, "< $wrk/temp$libraryname/$asm.split.allEdit") or die ("Couldn't open '$wrk/temp$libraryname/$asm.split.allEdit'", undef);
      open (O, "> $wrk/temp$libraryname/$asm.split.edit") or die ("Couldn't open '$wrk/temp$libraryname/$asm.split.edit", undef);
      while (<F>) {
         s/^\s+//;
         s/\s+$//;
         chomp $_;

         my ($uid) = split '\s+', $_;
         if (defined($gkpSeqs{$uid}) && $gkpSeqs{$uid} == 1) {
            print O "frg uid $uid isdeleted 1\n";
         }
      }
      close(F);
      close(O);

      runCommand("$wrk", "$CA/gatekeeper --edit $wrk/temp$libraryname/$asm.split.edit $wrk/temp$libraryname/$asm.corrected.gkpStore > $wrk/temp$libraryname/$asm.split.out 2> $wrk/temp$libraryname/$asm.split.err");
      runCommand("$wrk", "$CA/gatekeeper -dumpfrg $wrk/temp$libraryname/$asm.corrected.gkpStore > $wrk/temp$libraryname/$asm.correctedUnsplit.frg");
      runCommand("$wrk", "$CA/gatekeeper -dumpfasta $wrk/temp$libraryname/corrected.roundOne $wrk/temp$libraryname/$asm.corrected.gkpStore");
      runCommand("$wrk", "mv $wrk/temp$libraryname/corrected.roundOne.fasta.qual $wrk/temp$libraryname/corrected.roundOne.qual");
   
      # clean up intermediate files
      runCommand("$wrk", "rm -rf $wrk/temp$libraryname/$libraryname.frg");
      runCommand("$wrk", "rm -rf $wrk/temp$libraryname/asm.corrected.*");
      runCommand("$wrk", "rm -rf $wrk/temp$libraryname/$asm.split.[12].fastq");
      runCommand("$wrk", "rm -rf $wrk/temp$libraryname/$asm.split.paired.fastq");
   
      #run second round
      if ($submitToGrid == 1) {
         submitRunPBcR("$wrk/temp$libraryname", $asm, "$commandLineOptions \"fastqFile=$wrk/temp$libraryname/$asm.split.unmated.fastq\" \"maxCorrectionRounds=1\" \"libraryname=roundTwo\" sgeName=pBcR_roundTwo_$asm$sgeName sgePropagateHold=\"pBcR_$asm$sgeName\" $wrk/temp$libraryname/$asm.correctedUnsplit.frg", "pbcr_roundTwo_$asm$sgeName");
      } else {
         runCommand("$wrk/temp$libraryname", "$CA/pacBioToCA $commandLineOptions \"fastqFile=$wrk/temp$libraryname/$asm.split.unmated.fastq\" \"maxCorrectionRounds=1\" \"libraryname=roundTwo\" $wrk/temp$libraryname/$asm.correctedUnsplit.frg");
     } 
   }
   runCommand("$wrk/temp$libraryname", "cat corrected.roundOne.fasta temproundTwo/corrected.fasta |awk '{if (match(\$1, \"^>\") && match(\$1, \",\")) print substr(\$1, 1, index(\$1, \",\")-1); else print \$0}' > $asm.merged.fasta");
   runCommand("$wrk/temp$libraryname", "cat corrected.roundOne.qual temproundTwo/corrected.qual |awk '{if (match(\$1, \"^>\") && match(\$1, \",\")) print substr(\$1, 1, index(\$1, \",\")-1); else print \$0}' > $asm.merged.qual");
   runCommand("$wrk/temp$libraryname", "mv $asm.merged.fasta corrected.fasta");
   runCommand("$wrk/temp$libraryname", "mv $asm.merged.qual corrected.qual");
   
   #runCommand("$wrk/temp$libraryname", "rm -rf $wrk/temp$libraryname/temproundTwo");
}

#save output from temporary directory
runCommand("$wrk/temp$libraryname", "cp corrected.fasta $wrk/$libraryname.fasta");
runCommand("$wrk/temp$libraryname", "cp corrected.qual  $wrk/$libraryname.qual");
runCommand("$wrk/temp$libraryname", "cp $asm.layout.err  $wrk/$libraryname.correction.err");
runCommand("$wrk/temp$libraryname", "cp $asm.layout.hist  $wrk/$libraryname.correction.hist");

if (defined($longReads) && $longReads == 1) {
   runCommand("$wrk", "$CA/trimFastqByQVWindow $wrk/$libraryname.fasta $wrk/$libraryname.qual > $wrk/$libraryname.fastq");
   runCommand("$wrk", "$CA/fastqToCA -libraryname $libraryname -technology pacbio -type sanger -reads $wrk/$libraryname.fastq > $wrk/$libraryname.frg");
} else {
   # should do trimming here
   runCommand("$wrk", "$CA/convert-fasta-to-v2.pl -pacbio -s $wrk/temp$libraryname/corrected.fasta -q $wrk/temp$libraryname/corrected.qual -l $libraryname > $wrk/$libraryname.frg");
}

my $numOutput = `ls $wrk/temp$libraryname/*sge.out* |wc -l`;
chomp $numOutput;
if ($numOutput != 0) {
   runCommand("$wrk/temp$libraryname", "cat `ls *sge.out*` > $wrk/temp$libraryname/corrected.err");
   runCommand("$wrk/temp$libraryname", "cp corrected.err  $wrk/$libraryname.err");
}

# generate summary reports
my %nameTranslation;
open(F, "< $wrk/temp$libraryname/$asm.gkpStore.fastqUIDmap") or die("Couldn't open '$wrk/temp$libraryname/$asm.gkpStore.fastqUIDmap'", undef);

while (<F>) {
   s/^\s+//;
   s/\s+$//;

   next if (m/^\s*\#/);
   next if (m/^\s*$/);
   my @splitLine=split /\s+/;
   $nameTranslation{$splitLine[1]} = $splitLine[2];
}
close(F);

my $maxCorrected = 0;
my $totalCorrected = 0;
my $totalCorrectedBP = 0;
open(F, "< $wrk/temp$libraryname/corrected.log") or die("Couldn't open '$wrk/temp$libraryname/corrected.log'", undef);
open(W, "> $wrk/$libraryname.log") or die("Couldn't open '$wrk/$libraryname.log'", undef);

while (<F>) {
   s/^\s+//;
   s/\s+$//;

   next if (m/^\s*\#/);
   next if (m/^\s*$/);
   my @splitLine=split /\s+/;
   my $len = $splitLine[4] - $splitLine[3];
   my $name = (defined($nameTranslation{$splitLine[0]}) ? $nameTranslation{$splitLine[0]} : $splitLine[0]);
   
   print W "INPUT_NAME\tOUTPUT_NAME\tSUBREAD\tSTART\tEND\tLENGTH\n"if ($totalCorrected == 0);
   print W "$name\t$splitLine[1]\t$splitLine[2]\t$splitLine[3]\t$splitLine[4]\t$len\n";
   
   $maxCorrected = $len if ($len > $maxCorrected);
   $totalCorrected++;
   $totalCorrectedBP += $len;
}
close(F);
close(W);
         
print STDERR "******** Correction Summary ********\n";
printf "Total corrected bp (before consensus): $totalCorrectedBP (%.3f).\n", ($totalCorrectedBP/$totalBP)*100;
printf "Longest sequence: $maxCorrected bp, mean: %.3f bp.\n", ($totalCorrectedBP/$totalCorrected);

# finally clean up the assembly directory
if ($cleanup == 1) {
   runCommand("$wrk", "rm -rf temp$libraryname/[0-9]*.*");
   runCommand("$wrk", "rm -rf temp$libraryname/$asm.[0-9]*.log");
   runCommand("$wrk", "rm -rf temp$libraryname/$asm.gkpStore");
   runCommand("$wrk", "rm -rf temp$libraryname/$asm.ovlStore");
}

# trim here?

print STDERR "********* Finished correcting $totalBP bp (using $totalCorrectingWith bp).\n";
