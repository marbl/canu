#!/usr/bin/env perl

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

my @nonCAOptions = ("QV", "genomeSize", "shortReads", "longReads", "ovlMemory", "assemble", "assembleCoverage","mhapPrecompute", "mhap", "secret", "localStaging", "pbcns", "bankPath","libraryname", "specFile", "length", "coverageCutoff", "maxCoverage", "falconcns", "maxGap", "maxUncorrectedGap", "samFOFN", "blasr", "bowtie", "threads", "repeats", "fastqFile", "partitions", "submitToGrid", "sgeCorrection", "consensusConcurrency", "cleanup","asmOBT", "asmOvlErrorRate","asmUtgErrorRate","javaPath", "pythonPath", "asmCns", "asmMerSize", "asmCnsErrorRate","asmCgwErrorRate","asmObtErrorRate","asmObtErrorLimit");

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
    if ($var eq "utgErrorRate") {
       setGlobal("cgwErrorRate", $val);
       setGlobal("cnsErrorRate", $val);
       setGlobal("ovlErrorRate", $val);
    }
    if (($var eq "gridEngine") && ($val eq "SGE")) {
        setGlobal("gridSubmitCommand",      "qsub");
        setGlobal("gridHoldOption",         "-hold_jid \"WAIT_TAG\"");
        setGlobal("gridHoldOptionNoArray",  undef);
        setGlobal("gridSyncOption",         "-sync y");
        setGlobal("gridNameOption",         "-cwd -N");
        setGlobal("gridArrayOption",        "-t ARRAY_JOBS");
        setGlobal("gridArrayName",          "ARRAY_NAME");
        setGlobal("gridOutputOption",       "-j y -o");
        setGlobal("gridPropagateCommand",   "qalter -hold_jid \"WAIT_TAG\"");
        setGlobal("gridNameToJobIDCommand", undef);
        setGlobal("gridNameToJobIDCommandNoArray", undef);
        setGlobal("gridTaskID",             "SGE_TASK_ID");
        setGlobal("gridArraySubmitID",      "\\\$TASK_ID");
        setGlobal("gridJobID",              "JOB_ID");
    }

    if (($var eq "gridEngine") && ($val eq "LSF")) {
        setGlobal("gridSubmitCommand",      "bsub");
        setGlobal("gridHoldOption",         "-w \"numended\(\"WAIT_TAG\", \*\)\"");
        setGlobal("gridHoldOptionNoArray",  "-w \"done\(\"WAIT_TAG\"\)\"");
        setGlobal("gridSyncOption",         "-K");
        setGlobal("gridNameOption",         "-J");
        setGlobal("gridArrayOption",        "");
        setGlobal("gridArrayName",          "ARRAY_NAME\[ARRAY_JOBS\]");
        setGlobal("gridOutputOption",       "-o");
        setGlobal("gridPropagateCommand",   "bmodify -w \"done\(\"WAIT_TAG\"\)\"");
        setGlobal("gridNameToJobIDCommand", "bjobs -A -J \"WAIT_TAG\" | grep -v JOBID");
        setGlobal("gridNameToJobIDCommandNoArray", "bjobs -J \"WAIT_TAG\" | grep -v JOBID");
        setGlobal("gridTaskID",             "LSB_JOBINDEX");
        setGlobal("gridArraySubmitID",      "%I");
        setGlobal("gridJobID",              "LSB_JOBID");
    }

    # legacy support
    if ($var eq "secret") {
       setGlobal("mhap", $val);
    }

    if ($var eq "falconcns" && ($val eq 1)) {
       setGlobal("pbcns", 0);
    }

    if ($var eq "pbcns" && ($val eq 1)) {
       setGlobal("falconcns", 0);
    }

    if ($var eq "threads") {
       setGlobal("merylThreads",$val);
       setGlobal("ovlThreads", $val);
       setGlobal("consensusConcurrency", $val);
    }

    $commandLineOptions .= " \"$var=$val\" ";
}

sub initializeLocalSystem() {
    if (defined(getGlobal("submitToGrid")) && getGlobal("submitToGrid") == 1) {
       if (!defined(getGlobal("threads"))) {
          $global{"threads"} = 16;
       }
       if (!defined(getGlobal("ovlMemory"))) {
          $global{"ovlMemory"} = 32;
       }
       if (!defined(getGlobal("merylMemory"))) {
          $global{"merylMemory"} = 32000;
       }
       if (!defined(getGlobal("merylThreads"))) {
          $global{"merylThreads"} = 16;
       }
       if (!defined(getGlobal("ovlStoreMemory"))) {
          $global{"ovlStoreMemory"} = 32000;
       }
       if (!defined(getGlobal("ovlThreads"))) {
          $global{"ovlThreads"} = 16;
       }
       if (!defined(getGlobal("ovlConcurrency"))) {
           $global{"ovlConcurrency"} = 1;
       }
       if (!defined(getGlobal("consensusConcurrency"))) {
          $global{"consensusConcurrency"} = 16;
       }

       return;
    }

    # local run, get system info
    my $numCPU = 16;
    my $mem = 32 * 1024 * 1024 * 1024;

    if ( -e "/proc/cpuinfo" && -e "/proc/meminfo") {
       $numCPU = `cat /proc/cpuinfo |grep processor |wc -l`;
       chomp $numCPU;
       $mem = `cat /proc/meminfo |grep MemTotal: |awk '{printf(\"\%d\\n\", \$2*1024)}'`;
       chomp $mem;
    } else {
       my $sysInfo = `which sysctl`;
       chomp $sysInfo;
       if ( $sysInfo ne "") {
          $numCPU = `$sysInfo -n hw.ncpu`;
          chomp $numCPU;
          $mem = `$sysInfo -n hw.memsize`;
          chomp $mem;
       }
    }
    if (!defined(getGlobal("threads"))) {
       $global{"threads"} = $numCPU;
    }
    if (!defined(getGlobal("ovlMemory"))) {
       $global{"ovlMemory"} = int $mem / 1024 / 1024 / 1024;
    }
    if (!defined(getGlobal("merylMemory"))) {
       $global{"merylMemory"} = int $mem / 1024 / 1024;
    }
    if (!defined(getGlobal("merylThreads"))) {
       $global{"merylThreads"} = $numCPU;
    }
    if (!defined(getGlobal("ovlStoreMemory"))) {
       $global{"ovlStoreMemory"} = int $mem / 1024 / 1024;
    }
    if (!defined(getGlobal("ovlThreads"))) {
       $global{"ovlThreads"} = $numCPU;
    }
    if (!defined(getGlobal("ovlConcurrency"))) {
       $global{"ovlConcurrency"} = 1;
    }
    if (!defined(getGlobal("consensusConcurrency")) && !defined(getGlobal("cnsConcurrency"))) {
       $global{"consensusConcurrency"} = $numCPU;
       $global{"cnsConcurrency"} = $numCPU;
    }
}

sub setDefaults() {
    # grid submission options, duplicate of runCA
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
    $global{"gridJobID"}		   = "JOB_ID";
    $global{"gridArraySubmitID"}           = "\$TASK_ID";       # for lsf it is %I
    
    $global{"shell"}                       = "/bin/bash";
	
    $global{"utgErrorRate"} = 0.25;
    $global{"utgErrorLimit"} = 6.5;
    $global{"ovlErrorRate"} = 0.25;
    $global{"cnsErrorRate"} = 0.25;
    $global{"cgwErrorRate"} = 0.25;
	

    $global{"ovlHashBlockLength"} = 1000000000;
    $global{"ovlRefBlockLength"} = 100000000000;
    $global{"ovlRefBlockSize"} = 0;
    $global{"shortReads"} = undef;
    $global{"longReads"} = undef;
    $global{"QV"} = 54.5;
    $global{"libraryname"} = undef;
    $global{"falconcns"} = 0;
    $global{"pbcns"} = 1;
    $global{"bankPath"} = undef;
    $global{"specFile"} = undef;
    $global{"length"} = 500;
	
    # algorithm defaults
    $global{"frgMinLen"} = 64;
    $global{"coverageCutoff"} = 0;
    $global{"maxCoverage"} = 40;
    $global{"genomeSize"} = 0;
    $global{"maxUncorrectedGap"} = 1500;
    $global{"repeats"} = "";
    $global{"samFOFN"} = undef;
    $global{"blasr"} = undef;
    $global{"bowtie"} = undef;
    $global{"mhap"} = undef;
    $global{"mhapPrecompute"} = 1;
    $global{"merSize"} = 16;
    $global{"doOverlapBasedTrimming"} = 0;  # should this be defaulted to true or false?

    # now the assembly parameters
    $global{"assemble"} = 1;
    $global{"asmCns"} = "pbutgcns";
    $global{"cnsMinFrags"} = 1000;
    $global{"asmMerSize"} = 22;
    $global{"assembleCoverage"} = 25;
    $global{"unitigger"} = "bogart";
    $global{"asmOvlErrorRate"} = 0.03;
    $global{"asmUtgErrorRate"} = 0.013;
    $global{"asmCgwErrorRate"} = 0.10;
    $global{"asmCnsErrorRate"} = 0.10;
    $global{"asmOBT"} = 1;
    $global{"utgGraphErrorLimit"} = 0;
    $global{"utgGraphErrorRate"} = 0.013;
    $global{"utgMergeErrorLimit"} = 0;
    $global{"utgMergeErrorRate"} = 0.02;
    $global{"asmObtErrorRate"} = 0.03;
    $global{"asmObtErrorLimit"} =4.5;
    $global{"doFragmentCorrection"}=1;
    $global{"doUnitigSplitting"}=0;

    # default thread jobs/grid params based on system
    $global{"partitions"} = 200;
    $global{"submitToGrid"} = 0;
    $global{"sge"} = undef;
    $global{"sgeCorrection"} = "-pe threads 16 -l mem=2GB";
    $global{"sgeOverlap"} = "-pe threads 16 -l mem=2GB";
    $global{"sgeConsensus"} = "-pe threads 16";
    $global{"sgeScript"} = "-pe threads 1";
    $global{"useGrid"} = 0;
    $global{"scriptOnGrid"} = 0;

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
	my $taskID = getGlobal("gridJobID");
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
       if (getGlobal("gridEngine") eq "LSF"){
          my $tcmd = getGlobal("gridNameToJobIDCommand");
          $tcmd =~ s/WAIT_TAG/$waitTag/g;
          my $propJobCount = `$tcmd |wc -l`;
          chomp $propJobCount;
          if ($propJobCount == 0) {
             $tcmd = getGlobal("gridNameToJobIDCommandNoArray");
             $tcmd =~ s/WAIT_TAG/$waitTag/g;
             $hold = getGlobal("gridHoldOptionNoArray");
             $propJobCount = `$tcmd |wc -l`;
          }
          if ($propJobCount != 1) {
             print STDERR "Warning: multiple IDs for job $sgePropHold got $propJobCount and should have been 1.\n";
          }
          my $jobID = `$tcmd |tail -n 1 |awk '{print \$1}'`;
          chomp $jobID;
          $hold =~ s/WAIT_TAG/$jobID/g;
       } else{
          $hold =~ s/WAIT_TAG/$waitTag/g;
       }
       $waitTag = $hold;
    }

    my $qcmd = "$submitCommand $sge $sgeScript $nameOption \"$jobName\" $waitTag $outputOption $output $script";
    runCommand($wrk, $qcmd) and caFailure("Failed to submit script.\n");

    if (defined($sgePropHold)) {
	if (defined($holdPropagateCommand)) {
           my $translateCmd = getGlobal("gridNameToJobIDCommandNoArray");
           
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
                  chomp $hold;
                  chomp $prop;
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
              system($acmd) and print STDERR "WARNING: Failed to reset hold_jid trigger on '$sgePropHold'.\n";
           }
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
    writeScriptHeader($wrk, $asm, $script, "PBcR", $commandLineOptions);
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

sub submitRunCAHelper($$$$$) {
   my $wrk = shift @_;
   my $asm = shift @_;
   my $options = shift @_;
   my $TAG = shift @_;
   my $submitSelf = shift @_;
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
       if ($submitSelf) { 
          submitScript($wrk, $asm, $TAG);
       }
   } else {
       print "Please execute:\nrunCA $options\n";
   }
}

sub submitRunCA($$$$) {
   my $wrk = shift @_;
   my $asm = shift @_;
   my $options = shift @_;
   my $TAG = shift @_;
   submitRunCAHelper($wrk, $asm, $options, $TAG, 1);
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
       writeScriptHeader($wrk, $asm, $script, "PBcR", $options);
       submit($wrk, "$syncOption $script", $output, $TAG, undef);
       submitScript($wrk, $asm, $TAG);
   } else {
       print "Please execute:\nPBcR $options\n";
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
my $MAX_TO_PRECOMPUTE = 2000000000; 
my $MIN_SELF_CORRECTION = 50;
my $MAX_SPLIT_PERCENTAGE = 0.30;
my $FALCON_ERATE_ADJUST = 1.0;
my $MIN_COVERAGE_TO_ASM = 10;

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
} elsif (defined(getGlobal("specFile"))) {
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

    } elsif ($arg eq "-pbCNS") {
        setGlobal("falconcns", 1);

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
       my $thread = shift @cmdArgs;
       if ($thread <= 0) { $thread = 1; }

       setGlobal("threads", $thread);
    } elsif ($arg eq "-l" || $arg eq "-libraryname") {
       my $name  = shift @cmdArgs;
       if ($name =~ /^[0-9]/) {
          print STDERR "Warning: library name must start with a character, not a number\n";
          $name = "lib$name";
       }
       setGlobal("libraryname", $name);

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

if ((scalar(@fragFiles) == 0) && (!defined(getGlobal("longReads")) || getGlobal("longReads") == 0)) {
      print STDERR "Warning: no frag files specified, assuming self-correction of pacbio sequences.\n";
      setGlobal("longReads", 1);
}
  
if (($err) || (!defined(getGlobal("fastqFile"))) || (!defined(getGlobal("specFile"))) || (!defined(getGlobal("libraryname")))) {
    print STDERR "usage: $0 [options] -libraryname <name> -s <specfile> -fastq <fastqfile> [optional frg files]\n";
    print STDERR "  -length <int>                Minimum length of PacBio sequences to correct/output.\n";
    print STDERR "  -partitions <int>            Number of partitions for consensus\n";
    print STDERR "  -libraryname <string>        Name of the library; freeformat text. Must be unique from any library names in the FRG files used for correction\n";
    print STDERR "  -threads <int>               Number of threads to use for correction. Defaults to available cores on the local system\n";
    print STDERR "  -shortReads                  Use if the sequences for correction are 100bp or shorter.\n";

    print STDERR "  -genomeSize	<int>          Specify the approximate genome size. This will be used to compute the maximum number of bases to correct\n";
    print STDERR "  -maxCoverage <int>           Maximum coverage of PacBio sequences to correct. Only the longest sequences adding up to this coverage will be corrected. Requires genomeSize to be specified. Defaults to 40X\n";
    print STDERR " \"localStaging=<string>\"    Specify a local path (such as /scratch) to use for caching overlap computation. Will speed up grid-based computation by avoiding disk contention. Only use when running on a grid system (SGE or LSF), not a single machine\n";
    print STDERR "\nAdvanced options (EXPERT):\n";
    print STDERR "  -maxGap <int>                The maximum uncorrected PacBio gap that will be allowed. When there is no short-read coverage for a region, by default the pipeline will split a PacBio sequence. This option will attempt to use other PacBio sequences to patch the gap and avoid splitting the read. Sequences where the gaps have no support will still be broken. For example, specifying 50, will mean any gap 50bp or smaller can have no short-read coverage (but has other PacBio sequence support) without splitting the PacBio sequence. Warning: this can allow more sequences that went through the SMRTbell to not be fixed.\n";
    print STDERR "  -coverageCutoff                Specify the pacBio coverage (integer) used to separate repeat copies instead of automatically estimating.\n";
    print STDERR "  \"blasr=<string>\"               Use blasr for overlap computation instead of CA's built-in overlapper. This parameter specifies the blasr parameters to use.\n";
    print STDERR "  \"bowtie=<string>\"              Use bowtie2 for overlap computation instead of CA's built-in overlapper. This parameter specifies the bowtie 2 parameters to use.\n";
    print STDERR "  samFOFN=<file>                 Skip overlap computation. Use the provided file of file names of SAM files as the overlaps instead. Any valid sam files should be accepted.\n";

    print STDERR "\nComplete documentation at http://wgs-assembler.sourceforge.net/\n\n";

    print STDERR "No fastq file specified. Please specify a fastq file containing PacBio sequences for correction using the -fastq parameter.\n" if (!defined(getGlobal("fastqFile")));
    print STDERR "No spec file defined. Please specify a spec file using the -s option.\n" if (!defined(getGlobal("specFile")));
    print STDERR "No library name provided. Please specify a name using the -library option.\n" if (!defined(getGlobal("libraryname")));
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

# finally configure for local system if we are not on grid
initializeLocalSystem();

print STDERR "Running with " . getGlobal("threads") . " threads and " . getGlobal("partitions") . " partitions\n";

# try to find the tools we optionally support
my $CA = getBinDirectory();
my $JELLYFISH = getBinDirectory();
my $AMOS = "$CA/../../../AMOS/bin/";
my $BLASR = "$CA/../../../smrtanalysis/current/analysis/bin/";
my $FALCON = "$CA/../../../FALCON-0.1.2/bin/";
my $BOWTIE = "$CA/../../../bowtie2/";
my $MHAP_OVL = "$CA/../lib/java/mhap-0.1-ob.jar";
my $JELLYFISH = "$CA/../../../jellyfish/bin/";
my $wrk = makeAbsolute("");
my $asm = "asm";
my $scriptParams = getGlobal("sgeScript");
my $javaPath = "java";
if (defined(getGlobal("javaPath"))) {
   $javaPath = getGlobal("javaPath") . "/java";
}
my $pythonPath = "python";
if (defined(getGlobal("pythonPath"))) {
   $pythonPath = getGlobal("pythonPath") . "/python";
}
if (defined($scriptParams)) {
   if (!defined(getGlobal("sgeCorrection"))) {
   	setGlobal("sgeCorrection", $scriptParams);
   }
}

# check for previously existing file and do not overwrite
if ( -e "$wrk/" . getGlobal("libraryname") . ".frg") {
   print STDERR "Error: requested to output " . getGlobal("libraryname") . ".frg but file already exists. Will not overwrite.\n";
goto assemble;
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

if ( ! -e $MHAP_OVL ) {
  setGlobal("mhap", undef);
}

if (! -e "$JELLYFISH/jellyfish") {
   if (-e "$CA/jellyfish") {
      $JELLYFISH = $CA;
   } else {
      # try to use path
      my $amosPath = `which jellyfish`;
      chomp $amosPath;
      my @t = split '/', "$amosPath";
      pop @t;                      #  blasr
      $JELLYFISH = join '/', @t;  #  path to the assembler
   }
}

if (! -e "$BLASR/blasr") {
   if (-e "$CA/blasr") {
      $BLASR = $CA;
   } else {
      # try to use path
      my $amosPath = `which blasr`;
      chomp $amosPath;
      my @t = split '/', "$amosPath";
      pop @t;                      #  blasr 
      $BLASR = join '/', @t;  #  path to the assembler
   }
   # if we really can't find it just give up
   if (! -e "$BLASR/blasr") {
      die "BLASR binaries: blasr not found in $ENV{PATH}. Please download it from http://pacificbiosciences.github.com/DevNet/ and add it to your path.\n" if (defined(getGlobal("blasr")) || (defined(getGlobal("longReads")) && getGlobal("longReads") == 1 && ! -e $MHAP_OVL));
   }

   # check for consensus too
   # make sure we have the pb consensus module available if it was requested
   if (! -e "$BLASR/blasr" || ! -e "$BLASR/pbdagcon") {
      print STDERR "Warning: requested PBDAGON but either BLASR or pbdagcon executables were not found.\n";
      setGlobal("pbcns", 0);
   }
}

if (! -e "$FALCON/falcon_sense") {
   if (-e "$CA/falcon_sense") {
      $FALCON = $CA;
   } else {
      # try to use path
      my $amosPath = `which falcon_sense`;
      chomp $amosPath;
      my @t = split '/', "$amosPath";
      pop @t;                      #  falcon_sense
      $FALCON = join '/', @t;  #  path to the assembler
   }
   # if we really can't find it just give up
   if (! -e "$FALCON/falcon_sense") {
      if (-e "$BLASR/pbdagcon" && defined(getGlobal("longReads")) && getGlobal("longReads") == 1) {
         setGlobal("pbcns", 1);
      } else {
         setGlobal("pbcns", 0);
      }
      setGlobal("falconcns", 0);
   }
}

my $blasrVersion = 1.2;
if (-e "$BLASR/blasr") {
   $blasrVersion = `$BLASR/blasr -version |tail -n 1`;
   chomp $blasrVersion;
}

if (! -e "$AMOS/bank-transact" && (!defined(getGlobal("pbcns")) || getGlobal("pbcns") == 0)) {
   # check ca path
   if ( -e "$CA/bank-transact") {
      $AMOS = $CA;
   } else {
      # try to use path
      my $amosPath = `which bank-transact`;
      chomp $amosPath;
      my @t = split '/', "$amosPath";
      pop @t;                      #  bank-transact
      $AMOS = join '/', @t;  #  path to the assembler
   }
   # if we really can't find it just give up
   if (! -e "$AMOS/bank-transact") {
      die "AMOS binaries: bank-transact not found in $ENV{PATH}. Please download it from http://amos.sf.net and add it to your path.\n" if (! -e "$BLASR/pbdagcon" && ! -e "$FALCON/falcon_sense");
   }
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
my $ovlMemory = getGlobal("ovlMemory");

my $maxUncorrectedGap = getGlobal("maxUncorrectedGap");
my $coverage = getGlobal("coverageCutoff");
my $genomeSize = getGlobal("genomeSize");
my $shortReads = getGlobal("shortReads");
my $longReads = getGlobal("longReads");
my $length = getGlobal("length");
my $repeats = getGlobal("repeats");
my $QV = getGlobal("QV");

my $threads = getGlobal("threads");
my $partitions = getGlobal("partitions");

my $submitToGrid = getGlobal("submitToGrid");
my $sge = getGlobal("sge");
my $sgeOvl = getGlobal("sgeOverlap");
my $sgeCorrection = getGlobal("sgeCorrection");
my $sgeConsensus = getGlobal("sgeConsensus");
my $sgeTaskID = getGlobal("gridTaskID");

my $bankPath = getGlobal("bankPath");

my $consensusConcurrency = getGlobal("consensusConcurrency");

my $cleanup = getGlobal("cleanup");

my $cmd = "";

if (defined(getGlobal("longReads")) && getGlobal("longReads") == 1) {
    setGlobal("frgMinLen", 200);
}

if (! -e "$wrk/temp$libraryname") {
   runCommand("$wrk", "mkdir temp$libraryname");
   # generate the ca spec file, since we support additional options in the spec file, we need to strip those out before passing it to ca
   updateSpecFile("$wrk/temp$libraryname/$libraryname.spec");
}
$specFile = "$wrk/temp$libraryname/$libraryname.spec";

if (! -e "$wrk/temp$libraryname/$libraryname.frg") {
   runCommand($wrk, "$CA/fastqToCA -libraryname $libraryname -type sanger -technology none -feature doConsensusCorrection 1 -reads " . makeAbsolute($fastqFile) . " > $wrk/temp$libraryname/$libraryname.frg"); 
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
my $numLib = `$CA/gatekeeper -dumpinfo $wrk/temp$libraryname/$asm.gkpStore | grep LIB |awk '{print \$1}'`;
chomp($numLib);

my $minCorrectLib = 0;
my $maxCorrectLib = 0;
my $libToCorrect = 0;
for (my $i = 1; $i <= $numLib; $i++) {
   if (system("$CA/gatekeeper -isfeatureset $i doConsensusCorrection $wrk/temp$libraryname/$asm.gkpStore") == 0) {
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
   $ovlErrorRate = 0.35;
   $ovlErrorLimit = 6.5;
} elsif (($libToCorrect == 0 || ($minCorrectLib == 0 && $maxCorrectLib == 0))) {
   die ("Error: unable to find a library to correct. Please double-check your input files and try again.", undef);
}
print STDERR "Will be correcting PacBio library $libToCorrect with librarie[s] $minCorrectLib - $maxCorrectLib\n";
my $totalBP = 0;
my $totalNumToCorrect = 0;
my $minCorrectID = 0;
my $maxCorrectID = 0;
# compute the number of bases in the gateeeker to be corrected
open(F, "$CA/gatekeeper -dumpinfo $wrk/temp$libraryname/$asm.gkpStore |") or die("Couldn't open gatekeeper store", undef);

while (<F>) {
   s/^\s+//;
   s/\s+$//;

   my @array = split '\s+';
   if ($#array == 8 && $array[0] == $libToCorrect) { 
      $totalBP = $array[6];
      $minCorrectID = $array[1];
      $maxCorrectID = $array[2];
   }
}
close(F); 

my $totalInputBP = 0;
open(F, "$CA/gatekeeper -dumpinfo $wrk/temp$libraryname/$asm.gkpStore |") or die("Couldn't open gatekeeper store", undef);
while (<F>) {
   s/^\s+//;
   s/\s+$//;

   my @array = split '\s+';
   if ($#array == 8) {
      if ($array[0] == $libToCorrect) {
         $totalInputBP = $array[6];
      }
   }
}
close(F);

my $seedLength = 0;
if (! -e "$wrk/temp$libraryname/$asm.toerase.out") {
   # here is where we filter for specified length as well as max of longest X of coverage for correction
   # use the genome size/coverage, if available to subset the sequences
   if ($genomeSize != 0 && getGlobal("maxCoverage") != 0) {
      $totalBP = $genomeSize * getGlobal("maxCoverage");
   }
   runCommand($wrk, "$CA/gatekeeper -dumpfragments -invert -tabular -longestovermin $libToCorrect $length -longestlength $libToCorrect $totalBP $wrk/temp$libraryname/$asm.gkpStore 2> $wrk/temp$libraryname/$asm.seedlength |awk '{if (!(match(\$1, \"UID\") != 0 && length(\$1) == " . length("UID") . ")) { print \"frg uid \"\$1\" isdeleted 1\"; } }' > $wrk/temp$libraryname/$asm.toerase.uid");
   runCommand($wrk, "$CA/gatekeeper --edit $wrk/temp$libraryname/$asm.toerase.uid $wrk/temp$libraryname/$asm.gkpStore > $wrk/temp$libraryname/$asm.toerase.out 2> $wrk/temp$libraryname/$asm.toerase.err");
}
open(F, "cat $wrk/temp$libraryname/$asm.seedlength |awk '{print \$NF}' |") or die ("Couldn't open seed length file", undef);
$seedLength = <F>;
chomp $seedLength;
close(F);

# compute the number of bases left after our filtering gateeeker to be corrected
my $totalCorrectingWith = 0;
my $totalNumCorrectingWith = 0;
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
         $totalNumCorrectingWith += $array[3] + $array[4];
      }
   }
}
close(F);

if (defined(getGlobal("longReads")) && getGlobal("longReads") != 0) {
   $totalCorrectingWith = $totalInputBP;
}

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
if ($genomeSize != 0 && floor($totalCorrectingWith / $genomeSize) > $MAX_CORRECTION_COVERAGE && !(defined(getGlobal("longReads")) || getGlobal("longReads") == 0)) {
	print STDERR "Warning: input a total of " . floor($totalCorrectingWith / $genomeSize) . " of high-accuracy coverage for correction. For best performance, at most $MAX_CORRECTION_COVERAGE is recommended.\n";
	# could randomly subsample here
}
if ($genomeSize != 0 && defined(getGlobal("longReads")) && getGlobal("longReads") != 0 && ($totalCorrectingWith / $genomeSize) < $MIN_SELF_CORRECTION) {
   print STDERR "Warning: performing self-correction with a total of " .floor($totalCorrectingWith / $genomeSize) . ". For best performance, at least $MIN_SELF_CORRECTION is recommended.\n";
}

if (defined($longReads) && $longReads == 1) {
   my $merSize = getGlobal("merSize");
   if ( -e $MHAP_OVL && !defined(getGlobal("blasr"))) {
      if (!defined(getGlobal("mhap"))) {
         print STDERR "Warning: enabling MHAP overlapper to align long reads to long reads.\n";
         setGlobal("mhap", "-k $merSize --num-hashes 512 --num-min-matches 3 --threshold 0.04");
      }
      if ($totalCorrectingWith > $MAX_TO_PRECOMPUTE) {
         setGlobal("mhapPrecompute", undef);
      }
      # check java version and stop if not found or too old
      if (defined(getGlobal("javaPath")) && ! -e $javaPath) {
         die "Error: java is required to use MHAP and is not found in $javaPath specified.\n";
      } elsif (!(defined(getGlobal("javaPath")))) {
         my $java = `which java`;
         chomp $java;
         my @t = split '/', "$java";
         pop @t;                      #  java
         my $jPath = join '/', @t;
         if (! -e "$jPath/java") {
            die "Error: java is required to use MHAP and is not found in $ENV{PATH}.\n";
         }
      }
      open(F, "$javaPath -version 2>&1 |") or die("Couldn't find java in $ENV{PATH} or $javaPath. Specify javaPath=<path to java> or add it to your path.", undef);
      while (<F>) {
         s/^\s+//;
         s/\s+$//;

         if (m/java version/) {
            my @vArray = split '\s+';
            my $v = $vArray[2];
            $v =~ s/\"//g;
            @vArray = split '\.', $v;
            if ($vArray[0] == 1 && $vArray[1] < 7) {
               die "Error: java version 1.7 or newer is required for MHAP and found version $v. Please update and try again.\n";
            }
         }
      }
      close(F);
   } elsif (-e "$BLASR/blasr")  {
      print STDERR "Warning: Blasr is required to align long reads to long reads. Switching blasr ON.\n";
      setGlobal("blasr", defined(getGlobal("blasr")) ? getGlobal("blasr") . " -maxLCPLength 16" : "-minReadLength 200 -maxScore -1000 -maxLCPLength 16");
   } else {
      die "Aligning long reads requires either MHAP or BLASR and neither was available. Please either add SMRTportal to your path or download MHAP.";
  }
  if ((!defined(getGlobal("pbcns")) || getGlobal("pbcns") == 0) && -e "$FALCON/falcon_sense") {
     setGlobal("falconcns", 1);
  }
} else {
   # always default to blasr alignments when available
   if (-e "$BLASR/blasr") {
      setGlobal("blasr", defined(getGlobal("blasr")) ? getGlobal("blasr") : "-noRefineAlign -advanceHalf -noSplitSubreads -minMatch 10 -minPctIdentity 70");
   }
   setGlobal("falconcns", 0);
   if (-e "$BLASR/pbdagcon") {
      setGlobal("pbcns", 1);
   } else {
      die "AMOS or PBDAGCON is required to call consensus and neither was available. Please add either SMRTportal or AMOS to your path and try again\n" if (! -e "$AMOS/bank-transact");
      setGlobal("pbcns", 0);
  }
}

if (-e $MHAP_OVL && defined(getGlobal("mhap"))) {
   $MHAP_OVL="\$bin/../lib/java/mhap-0.1-ob.jar";
}

my $cutoffSpecified = 0;

if (defined(getGlobal("blasr"))) {
   if (getGlobal("blasr") =~ m/bestn/) {
      $cutoffSpecified = 1;
      print STDERR "Warning: cutoff manually specified, skip meryl run\n";
   }
}
 
my $ignore = "";
if (!defined(getGlobal("bowtie")) && $cutoffSpecified == 0 && !(defined(getGlobal("mhap")) && -e "$JELLYFISH/jellyfish")) {
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

my $inMer = getGlobal("merSize");
my $inFile = 0;
if ( -e "$wrk/*.$inMer.ignore" ) {
    my $inFile = `ls $wrk/*.$inMer.ignore |wc -l |awk '{print \$1}'`;
}
chomp $inFile;
if ($inFile == 0 && !-e "$wrk/temp$libraryname/$asm.ignore" &&  -e "$JELLYFISH/jellyfish" && defined(getGlobal("mhap"))) {
   my $f = makeAbsolute($fastqFile);
   if ($f =~ /\.bz2$/i || $f =~ /\.zip$/i || $f =~ /\.gz$/i) {
      die "Input fastq file $f is compressed. This is currently not supported, please uncompress the file and try again.";
   }
   $cmd  = "$JELLYFISH/jellyfish count ";
   $cmd .= " -m " . getGlobal("merSize") . " -s 120000000 -t " . getGlobal("merylThreads") . " -o $wrk/temp$libraryname/$asm.mers " . makeAbsolute($fastqFile);
   runCommand("$wrk/temp$libraryname", $cmd);
   if ( -e "$wrk/temp$libraryname/$asm.mers_0") {
      $cmd = "$JELLYFISH/jellyfish merge -s 120000000 -o $wrk/temp$libraryname/$asm.mers $wrk/temp$libraryname/$asm.mers_*";
      runCommand("$wrk/temp$libraryname", $cmd);
   } 

   $cmd =  "$JELLYFISH/jellyfish histo -t " . getGlobal("merylThreads") . " -f $wrk/temp$libraryname/$asm.mers > $wrk/temp$libraryname/$asm.hist";
   runCommand("$wrk/temp$libraryname", $cmd);
   my $total = 0;
   my $runningSum = 0;
   my $sum = `cat $wrk/temp$libraryname/$asm.hist |awk '{SUM+=\$NF; } END {print SUM}'`;
   chomp $sum;
   my $cut = `cat $wrk/temp$libraryname/$asm.hist |awk -v TOTAL=$sum '{SUM+=\$NF; if (SUM/TOTAL > 0.99) { print \$1; } }' |head -n 1`;
   chomp $cut;
   $cut++;

   $cmd = "$JELLYFISH/jellyfish dump -c -t -L $cut $wrk/temp$libraryname/$asm.mers |awk -v TOTAL=$sum '{printf(\"\%s\\t\%0.10f\\t\%d\\t\%d\\n\", \$1, \$2/TOTAL, \$2, TOTAL)}' |sort -T . -rnk2> $wrk/temp$libraryname/$asm.ignore";
   runCommand("$wrk/temp$libraryname", $cmd);
   runCommand("$wrk/temp$libraryname", "rm $wrk/temp$libraryname/$asm.mers*");
} elsif (!-e "$wrk/temp$libraryname/$asm.ignore" &&  defined(getGlobal("mhap"))) {
   runCommand("$wrk/temp$libraryname", "cp $wrk/*.$inMer.ignore $wrk/temp$libraryname/$asm.ignore");
}
if ( -e "$wrk/temp$libraryname/$asm.ignore") {
   $ignore = " -f $wrk/temp$libraryname/$asm.ignore";
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

   # when we were asked to not use CA's overlapper, first perform common options
   if (defined(getGlobal("samFOFN")) || defined(getGlobal("blasr")) || defined(getGlobal("bowtie")) || defined(getGlobal("mhap"))) {
      # check for sam tools
      if (!-e "$CA/convertSamToCA") { 
         die("Error: request to use Blasr, Bowtie, or SAM files requires CA to be built with SAMTOOLS. Please rebuild CA and try again");
      }
      # create directories
      runCommand($wrk, "mkdir $wrk/temp$libraryname/1-overlapper") if (! -d "$wrk/temp$libraryname/1-overlapper");

      # first partition the data so we know how many jobs we have
      my $cmd;

      my $ovlHashBlockSize   = $totalNumToCorrect;
      my $ovlHashBlockLength = undef;
      my $maxBatch = 1000;

      if (defined(getGlobal("blasr"))) {
         $ovlHashBlockLength = getGlobal("ovlHashBlockLength");
         $ovlHashBlockSize   = getGlobal("ovlHashBlockSize");
      }
      my $ovlRefBlockSize    = getGlobal("ovlRefBlockSize");
      my $ovlRefBlockLength  = getGlobal("ovlRefBlockLength");
      if (defined(getGlobal("mhap"))) { 
         setGlobal("blasr", undef);
         setGlobal("bowtie", undef);
         if (!defined($ovlRefBlockSize) || $ovlRefBlockSize == 0) {
            my $memFactor = getGlobal("ovlMemory") / 16;
            my @mhapSplit = split(/\s+/, getGlobal("mhap"));
            my $numHashes = 1256;
            for my $pos (0 .. $#mhapSplit) {
               if ($mhapSplit[$pos] =~ m/num-hashes/) {
                  $numHashes = int($mhapSplit[$pos+1]);
                  last;
                }
            }
            my $loadFactor = 5000;
            if ($numHashes < 786) { $loadFactor = $loadFactor * 2; }
            $ovlRefBlockSize = floor($memFactor * $loadFactor);
         }
         $ovlRefBlockLength = undef; 
         $ovlHashBlockLength = undef;
         $ovlHashBlockSize = $ovlRefBlockSize;
      } 

      if (($ovlRefBlockSize > 0) && ($ovlRefBlockLength > 0)) {
          die("can't set both ovlRefBlockSize and ovlRefBlockLength", undef);
      }

      my @job = ();
      if (!defined(getGlobal("mhap"))) {
         $cmd  = "$CA/overlap_partition \\\n";
         $cmd .= " -g  $wrk/temp$libraryname/$asm.gkpStore \\\n";
         $cmd .= " -bs $ovlHashBlockSize \\\n" if (defined($ovlHashBlockSize));
         $cmd .= " -bl $ovlHashBlockLength \\\n" if (defined($ovlHashBlockLength));
         $cmd .= " -rs $ovlRefBlockSize \\\n" if (defined($ovlRefBlockSize));
         $cmd .= " -rl $ovlRefBlockLength \\\n"  if (defined($ovlRefBlockLength));
         $cmd .= " -H $libToCorrect \\\n";
         $cmd .= " -R $minCorrectLib-$maxCorrectLib \\\n";
         $cmd .= " -C \\\n";
         $cmd .= " -o  $wrk/temp$libraryname/1-overlapper";
         runCommand($wrk, $cmd);
         open(F, "< $wrk/temp$libraryname/1-overlapper/ovlopt") or die("failed partition for overlapper: no ovlopt file found", undef);
         @job = <F>;
         close(F);
      } else {
         open(BAT, "> $wrk/temp$libraryname/1-overlapper/ovlbat") or die ("Couldn't open $wrk/temp$libraryname/1-overlapper/ovlbat", undef);
         open(JOB, "> $wrk/temp$libraryname/1-overlapper/ovljob") or die ("Couldn't open $wrk/temp$libraryname/1-overlapper/ovljob", undef);
         open(OPT, "> $wrk/temp$libraryname/1-overlapper/ovlopt") or die ("Couldn't open $wrk/temp$libraryname/1-overlapper/ovlopt", undef);
         
         my $numJobs = 1;
         my $numInBatch = 1;
         my $numBatch = 1;
         my $maxRefsToCompare = floor($totalNumCorrectingWith / $ovlRefBlockSize / 2);
         if ($maxRefsToCompare == 0) {
           $maxRefsToCompare = floor($totalNumCorrectingWith / $ovlRefBlockSize);
         }
         my $maxToCompare = $maxRefsToCompare * $ovlRefBlockSize;
         for (my $i = 1; $i <= $totalNumCorrectingWith; $i+=$ovlRefBlockSize) {
            if ($numInBatch > $maxBatch) {
               $numBatch++;
               $numInBatch=1;
            }
            my $max = $i+$ovlHashBlockSize-1;
            my $next = $max+1;
            if ($max > $totalNumCorrectingWith) { $max = $totalNumCorrectingWith; }

            my $limit = $next;
            if (($totalNumCorrectingWith - $next + 1) > $maxToCompare) {
               # split into two jobs for better load balancing
               push(@job, "-h $i-$max -r $i-$max");
               printf(BAT "%03d\n", $numBatch);
               printf(JOB "%06d\n", $numJobs);
               $numJobs++;
               $numInBatch++;
               $limit = $next + $maxToCompare - 1;
               if ($next > $totalNumCorrectingWith) { $next = $totalNumCorrectingWith+1; }
               if ($limit > $totalNumCorrectingWith) { $limit = $totalNumCorrectingWith+1; }
               print OPT "-h $i-$max -r $next-$limit\n";
               $limit++;
            }
            push(@job, "-h $i-$max -r $i-$max"); 
            printf(BAT "%03d\n", $numBatch);
            printf(JOB "%06d\n", $numJobs);
            if ($next > $totalNumCorrectingWith) { $next = $totalNumCorrectingWith+1; }
            if ($limit > $totalNumCorrectingWith) { $limit = $totalNumCorrectingWith+1; }
            print OPT "-h $i-$max -r $limit-$totalNumCorrectingWith\n";
            $numJobs++;
            $numInBatch++;
         }
         close(BAT);
         close(JOB);
         close(OPT);
      }
    
      my $blasrThreshold = 0;
      my $ovlThreads = getGlobal("ovlThreads");   
      my $numOvlJobs = 0;
      my $numIndicies = 1;
      my $blasrOpts = getGlobal("blasr");
      my $bowtieOpts = getGlobal("bowtie");
      my $suffix = "sam";

      $blasrVersion="1.3.1.116174";
      # SMRTportal 1.4 changed parameters so update our blasr options
      if ($blasrVersion >= 1.3.1.116174) {
         $blasrOpts =~ s/-ignoreQuality//g;
      }
   
      if (defined(getGlobal("samFOFN"))) {
         my $fofn = getGlobal("samFOFN");
         $numOvlJobs = 0;
      } else {
         $numOvlJobs = $#job + 1;
      }

      my $wrkDir = $wrk;
      if (defined(getGlobal("localStaging"))) {
         $wrkDir = getGlobal("localStaging") . "/\$USER/";
      }

      if (!-e "$wrk/temp$libraryname/1-overlapper/ovlpreplocal.sh" && defined(getGlobal("localStaging"))) {
         open F, "> $wrk/temp$libraryname/1-overlapper/ovlpreplocal.sh" or die ("can't open '$wrk/temp$libraryname/1-overlapper/ovlpreplocal.sh'");
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
         print F "\n";
         my $wrkDir = $wrkDir = getGlobal("localStaging") . "/\$USER/";
         print F "rm -rf $wrkDir/temp$libraryname/1-overlapper/\n";
         close(F);
         chmod 0755, "$wrk/temp$libraryname/1-overlapper/ovlpreplocal.sh";

         if ($submitToGrid == 1) {
            my $sgeName = "pBcR_ovlpreplocal_$asm$sgeName";
            my $jobName = getGridArrayName($sgeName, $numOvlJobs);
            my $arrayOpt = getGridArrayOption($sgeName, $numOvlJobs);
           submitBatchJobs($wrk, $asm, "$submitCommand $sge $sgeOvl $nameOption \"$jobName\" $arrayOpt $outputOption /dev/null $wrk/temp$libraryname/1-overlapper/ovlpreplocal.sh", $jobName);
         } else {
            for (my $i = 1; $i <= $numOvlJobs; $i++) {
               schedulerSubmit("$wrk/temp$libraryname/1-overlapper/ovlpreplocal.sh $i");
            }
            schedulerSetNumberOfProcesses(getGlobal("ovlConcurrency"));
            schedulerFinish();
         }
      }

      if (-e "$wrk/temp$libraryname/1-overlapper/ovlpindex") {
         open(F, "< $wrk/temp$libraryname/1-overlapper/ovlpindex") or die("failed partition for overlapper: no ovlopt file found", undef);
         my @indicies = <F>;
         close(F);
         $numIndicies += $#indicies + 1 ;
      }

      if (-e "$wrk/temp$libraryname/$asm.ovlStore.BUILDING") {
          die "Error: overlap store building failed. See $wrk/temp$libraryname/$asm.ovlStore.err for more details";
      }
      if (-e "$wrk/temp$libraryname/1-overlapper/overlap.sh") {
         goto checkforerror;
      }
      if (-e "$wrk/temp$libraryname/1-overlapper/ovlprep.sh") {
         goto checkovlprep;
      }
      
      # dump gatekeeper store and print eid to iid, add to it the gkpStore fastq file into eidToIID and lengths
      runCommand("$wrk/temp$libraryname", "$CA/gatekeeper -dumpfragments -tabular $asm.gkpStore |awk '{print \$1\"\\t\"\$2}' > $asm.eidToIID");
      runCommand("$wrk/temp$libraryname", "$CA/gatekeeper -dumpfragments -tabular $asm.gkpStore |awk '{print \$2\"\\t\"\$10}' > $asm.iidToLen");
       open F, "> $wrk/temp$libraryname/1-overlapper/ovlindex" or die ("can't open '$wrk/temp$libraryname/1-overlapper/ovlindex");

   
      # now do specific steps (either import SAMs or run blasr)
      if (defined(getGlobal("samFOFN"))) {
         my $fofn = getGlobal("samFOFN");

         # get the external UIDs from the input file rather than gatekeeper since these won't match our externally-generated SAM
         runCommand("$wrk/temp$libraryname", "cat $asm.gkpStore.fastqUIDmap | awk '{if (NF > 3) { print \$NF\"\\t\"\$(NF-1); } print \$3\"\\t\"\$2; }' >> $asm.eidToIID");

         # convert sam to ovl and put it into the overlap directory
         my $batch = 1;

         open(BAT, "> $wrk/temp$libraryname/1-overlapper/ovlbat") or die ("Couldn't open $wrk/temp$libraryname/1-overlapper/ovlbat", undef);
         open(JOB, "> $wrk/temp$libraryname/1-overlapper/ovljob") or die ("Couldn't open $wrk/temp$libraryname/1-overlapper/ovljob", undef);
         open(OPT, "> $wrk/temp$libraryname/1-overlapper/ovlopt") or die ("Couldn't open $wrk/temp$libraryname/1-overlapper/ovlopt", undef);
         open(SA, "> $wrk/temp$libraryname/1-overlapper/ovlindex") or die ("Couldn't open $wrk/temp$libraryname/1-overlapper/ovlindex", undef);
         $numIndicies++;

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
            if ($numOvlJobs > $maxBatch) {
               $batch++; 
            }
            my $batchName = sprintf("%03d", $batch);
            my $jobName = sprintf("%06d", $numOvlJobs);
            print BAT "$batchName\n";
            print OPT "-r 1-$totalNumToCorrect -l 1-$totalNumToCorrect\n";
            print JOB "$jobName\n";
            print SA "1\n";
            runCommand("$wrk/temp$libraryname/1-overlapper", "mkdir -p $batchName");
            runCommand("$wrk/temp$libraryname/1-overlapper/$batchName", "unlink $wrk/temp$libraryname/1-overlapper/$batchName/$jobName.$suffix") if (-e "$wrk/temp$libraryname/1-overlapper/$batchName/$jobName.$suffix");
            runCommand("$wrk/temp$libraryname/1-overlapper/$batchName", "cp $_ $wrk/temp$libraryname/1-overlapper/$batchName/$jobName.$suffix");
         }
         close(F);
         close(BAT);
         close(OPT);
         close(JOB);
         close(SA);
      } else {
         if (!defined(getGlobal("mhap"))) {
            for (my $i = $minCorrectLib; $i <= $maxCorrectLib; $i++) {
              runCommand("$wrk/temp$libraryname/1-overlapper", "$CA/gatekeeper -dumpfasta " . $i . "_subset " . (defined($longReads) && $longReads == 1 && $i == $libToCorrect ? " -allreads -allbases" : "") . " -randomsubset $i 0.01 $wrk/temp$libraryname/$asm.gkpStore");
            }
            runCommand("$wrk/temp$libraryname/1-overlapper", "cat `ls [0-9]*_subset.fasta` > correct_subset.fasta");
            runCommand("$wrk/temp$libraryname/1-overlapper", "rm `ls [0-9]*_subset.fasta*`");
         }

         # read in ovl options and partition data from the gatekeeper store
         my $i = 1;
         my $part = 1;
         my %longIndex = {};
         my %correctIndex = {};
         my $numPrepJobs = 0;

         open(PREP, "> $wrk/temp$libraryname/1-overlapper/ovlprep") or die ("Couldn't open $wrk/temp$libraryname/1-overlapper/ovlprep", undef);
         open(PINDEX, "> $wrk/temp$libraryname/1-overlapper/ovlpindex") or die ("Couldn't open $wrk/temp$libraryname/1-overlapper/ovlpindex", undef);

         foreach my $ovlJob (@job) {
            $ovlJob =~ s/^\s+//;
            $ovlJob =~ s/\s+$//;
            my @values=split(/\s+/, $ovlJob);
            my $indexIds = "";
            my $correctIds = "";

            if (!defined($longIndex{$values[1]})) {
               # first the hash sequences
               my @startend=split(/-/, $values[1]);
               $indexIds = "-b $startend[0] -e $startend[1]";
               if (!defined(getGlobal("mhap"))) { 
                  print PREP "-randomsubset $libToCorrect 1 $indexIds\n";
                  print PINDEX "long_reads_part\t$numIndicies\n";
                  $numPrepJobs++;
               }
               $numIndicies++;
               $longIndex{$values[1]} = 1;
            }

            if (!defined($correctIndex{$values[3]})) {
               # now the correction sequences
               my @startend=split(/-/, $values[3]);
               $correctIds="-b $startend[0] -e $startend[1]";
               my $dumpAll = 0;
               my $start = $startend[0] > $minCorrectID ? $startend[0] : $minCorrectID;
               my $end = $startend[1] < $maxCorrectID ? $startend[1] : $maxCorrectID;
               if ($end - $start >= 0) {
                  $dumpAll = 1;
               }
               if (defined($longReads) && $longReads == 1 && $dumpAll) {
                  if (!defined(getGlobal("mhap"))) {
                     print PREP " -allreads -allbases $correctIds\n";
                  } else {
                     print PREP " -allreads -allbases $correctIds\n";
                  }
               } else {
                  print PREP " -randomsubset 0 1 $correctIds\n";
               }
               print PINDEX "correct_reads_part\t$i\n";
               $numPrepJobs++;
               $i++;
               $correctIndex{$values[3]} = 1;
            }

            print F $numIndicies-1 . "\n";
            $part++;
         }
         close(PREP);
         close(PINDEX);
         close(F);

         open F, "> $wrk/temp$libraryname/1-overlapper/ovlprep.sh" or die ("can't open '$wrk/temp$libraryname/1-overlapper/ovlprep.sh'");
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
         print F "job=`head -n \$jobid $wrk/temp$libraryname/1-overlapper/ovlprep | tail -n 1`\n";
         print F "index=`head -n \$jobid $wrk/temp$libraryname/1-overlapper/ovlpindex |tail -n 1`\n";
         print F "name=`echo \"\$index\" |awk '{print \$1}'`\n";
         print F "id=`echo \"\$index\" |awk '{print \$2}'`\n";
         print F "\n";
         if (defined(getGlobal("localStaging"))) {
            print F "mkdir -p $wrkDir/temp$libraryname/1-overlapper/\n";
         }
         print F "   \$bin/gatekeeper -dumpfasta $wrkDir/temp$libraryname/1-overlapper/\$name\$id \$job $wrk/temp$libraryname/$asm.gkpStore\n";
         if (defined(getGlobal("blasr"))) {
            print F "if [ \$name == \"long_reads_part\" ]; then\n";
            print F "   /usr/bin/time $BLASR/sawriter $wrkDir/temp$libraryname/1-overlapper/\$name\$id.sa $wrkDir/temp$libraryname/1-overlapper/\$name\$id.fasta > $wrk/temp$libraryname/1-overlapper/\$jobid.hash.err 2>&1 \n";
           print F "fi\n";
         }
         if (!defined(getGlobal("mhap"))) {
            print F "\$bin/gatekeeper -dumpfragments -tabular \$job $wrk/temp$libraryname/$asm.gkpStore |awk '{print \$1\"\\t\"\$2}' > $wrk/temp$libraryname/1-overlapper/\$name\$id.eidToIID\n";
            print F "\$bin/gatekeeper -dumpfragments -tabular \$job $wrk/temp$libraryname/$asm.gkpStore |awk '{print \$2\"\\t\"\$10}' > $wrk/temp$libraryname/1-overlapper/\$name\$id.iidToLen\n";
         } elsif (!defined(getGlobal("localStaging")) && defined(getGlobal("mhapPrecompute"))) {
            print F "/usr/bin/time $javaPath -XX:+UseG1GC -server -Xmx" . $ovlMemory . "g -jar $MHAP_OVL FastAlignMain " . getGlobal("mhap") . " --min-store-length " . ($seedLength-1) . " --num-threads $ovlThreads" . " $ignore -p $wrkDir/temp$libraryname/1-overlapper/\$name\$id.fasta -q $wrkDir/temp$libraryname/1-overlapper/ > $wrk/temp$libraryname/1-overlapper/\$jobid.hash.err 2>&1\n";
            print F "rm -f $wrkDir/temp$libraryname/1-overlapper/\$name\$id.fasta\n";
         }
         # remove qual files they are not currently used
         print F "rm -f $wrkDir/temp$libraryname/1-overlapper/\$name\$id.fasta.q*\n";

         if (defined(getGlobal("localStaging"))) {
            if (!defined(getGlobal("mhap"))) {
               print F "mv $wrkDir/temp$libraryname/1-overlapper/\$name\$id.fasta $wrk/temp$libraryname/1-overlapper/\n";
               if (defined(getGlobal("blasr"))) {
                  print F "mv $wrkDir/temp$libraryname/1-overlapper/\$name\$id.sa $wrk/temp$libraryname/1-overlapper/\n";
               }
            } else {
               print F "if [ -e \"$wrkDir/temp$libraryname/1-overlapper/\$name\$id.dat\" ]; then\n";
               print F "mv $wrkDir/temp$libraryname/1-overlapper/\$name\$id.dat $wrk/temp$libraryname/1-overlapper/\n";
               print F "else\n";
               print F "mv $wrkDir/temp$libraryname/1-overlapper/\$name\$id.fasta $wrk/temp$libraryname/1-overlapper/\n";
               print F "fi\n";
            }
         }

         close(F);
         chmod 0755, "$wrk/temp$libraryname/1-overlapper/ovlprep.sh";

         if ($submitToGrid == 1) {
            my $sgeName = "pBcR_ovlprep_$asm$sgeName";
            my $jobName = getGridArrayName($sgeName, $numPrepJobs);
            my $arrayOpt = getGridArrayOption($sgeName, $numPrepJobs);
           submitBatchJobs($wrk, $asm, "$submitCommand $sge $sgeOvl $nameOption \"$jobName\" $arrayOpt $outputOption /dev/null $wrk/temp$libraryname/1-overlapper/ovlprep.sh", $jobName);
         } else {
            for (my $i = 1; $i <= $numPrepJobs; $i++) {
               schedulerSubmit("$wrk/temp$libraryname/1-overlapper/ovlprep.sh $i");
            }
            schedulerSetNumberOfProcesses(getGlobal("ovlConcurrency"));
            schedulerFinish();
         }
  checkovlprep:
         my $failedJobs = 0;
         my $failureMessage = "";
         open(B, "< $wrk/temp$libraryname/1-overlapper/ovlpindex") or die("failed to open '$wrk/temp$libraryname/1-overlapper/ovlpindex'", undef);
         while (!eof(B)) {
            my $b = <B>;  chomp $b;
            my ($name, $id) = split '\s+', $b;
            my $fileName = "$wrk/temp$libraryname/1-overlapper/$name$id";

            if (defined(getGlobal("mhap"))) { 
               if ((! -e "$fileName.dat") && ! -e "$fileName.fasta") { 
                  $failureMessage .= "ERROR:  Overlap prep job $wrk/temp$libraryname/1-overlapper/$b FAILED.\n";
                  $failedJobs++;
               }
            } else {
               if (defined(getGlobal("blasr")) && $fileName =~ m/long_reads/) {
                  if ((! -e "$fileName.sa")) {  
                     $failureMessage .= "ERROR:  Overlap prep job $wrk/temp$libraryname/1-overlapper/$b FAILED.\n";
                     $failedJobs++;
                  }
               }
               if (! -e "$fileName.fasta") {  
                  $failureMessage .= "ERROR:  Overlap prep job $wrk/temp$libraryname/1-overlapper/$b FAILED.\n";
                  $failedJobs++;
               }
            }
         }
         close(B);
         $failureMessage .= "\n$failedJobs overlap partitioning jobs failed.";
         if ($failedJobs > 0) {
             die "$failureMessage\n";
         }


         if (defined(getGlobal("blasr")) && $cutoffSpecified == 0) {
            runCommand("$wrk/temp$libraryname/1-overlapper", "ln -s long_reads_part1.fasta long_reads.fasta");
            # run with best limit set to kmer with a small (0.5-1% of the data)   
            my $ovlThreads = getGlobal("ovlThreads");
            runCommand("$wrk/temp$libraryname/1-overlapper", "ln -s long_reads_part1.sa long.sa");

            open F, "> $wrk/temp$libraryname/1-overlapper/pickNBest.sh" or die ("can't open '$wrk/temp$libraryname/1-overlapper/pickNBest.sh");
            print F "#!" . getGlobal("shell") ."\n";
            print F getBinDirectoryShellCode();
            print F "\n";

            if (-e "$BLASR/../../etc/setup.sh") {
               print F "source $BLASR/../../etc/setup.sh && /usr/bin/time $BLASR/blasr -sa long.sa correct_subset.fasta long_reads.fasta -nproc $ovlThreads -bestn $ovlThreshold $blasrOpts -sam -out subset.sam";
            } else {
               print F "/usr/bin/time $BLASR/blasr -sa long.sa correct_subset.fasta long_reads.fasta -nproc $ovlThreads -bestn $ovlThreshold $blasrOpts -sam -out subset.sam";
            }
            close(F);
            chmod 0755, "$wrk/temp$libraryname/1-overlapper/pickNBest.sh";
            runCommand("$wrk/temp$libraryname/1-overlapper", "$wrk/temp$libraryname/1-overlapper/pickNBest.sh");

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
      close(F);
   
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
      print F "sa=`head -n \$jobid $wrk/temp$libraryname/1-overlapper/ovlindex | tail -n 1`\n";
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
      print F "hashParams=(\$(echo \${options[1]}  | tr \"-\" \"\\n\")) \n";
      print F "hashStart=\${hashParams[0]} \n";
      print F "hashStart=\$((\$hashStart - 1))\n";
      print F "hashEnd=\${hashParams[1]} \n";
      print F "params=(\$(echo \${options[3]}  | tr \"-\" \"\\n\")) \n";
      print F "numTotalJobs=$numOvlJobs\n";
      print F "numSA=" , ($numIndicies-1) , "\n";
      print F "numOvlJobs=\$((\$numTotalJobs/\$numSA))\n";
      print F "zeroJobId=\$((\$((\$jobid - 1)) % \$numOvlJobs))\n";
      print F "oneJobId=\$((\$zeroJobId+1))\n";
      print F "start=\${params[0]} \n";
      print F "start=\$((\$start - 1))\n";
      print F "end=\${params[1]} \n";
      print F "total=\$((\$end - \$start))\n"; 
      print F "refOffset=\$((\$start - \$hashEnd))\n";
      print F "echo \"Running partition \$job with options \$opt start \$start end \$end total \$total zero job \$zeroJobId and stride \$numOvlJobs\"\n";
      print F "\n";

      my $wrkDir = $wrk;
      my $suffix = "dat";
      if (defined(getGlobal("mhap"))) {
         if (defined(getGlobal("localStaging"))) {
            $wrkDir = getGlobal("localStaging") . "/\$USER/";
            $suffix = "fasta";
            print F "mkdir -p $wrkDir/temp$libraryname/1-overlapper/\n";
         } elsif (!defined(getGlobal("mhapPrecompute"))) {
            $suffix = "fasta";
         }
      } else { 
         $suffix = "sam";
      }
      if (defined(getGlobal("blasr"))) {
         if (-e "$BLASR/../../etc/setup.sh") {
            print F "   source $BLASR/../../etc/setup.sh \n";
         }
         if (defined(getGlobal("localStaging"))) {
         print F " cp $wrk/temp$libraryname/1-overlapper/long_reads_part\$sa.sa $wrkDir\n";
         print F " cp $wrk/temp$libraryname/1-overlapper/long_reads_part\$sa.fasta $wrkDir\n";
         print F " cp $wrk/temp$libraryname/1-overlapper/correct_reads_part\$oneJobId.fasta $wrkDir\n"
         }
         print F "   /usr/bin/time $BLASR/blasr \\\n";
         print F "          -sa $wrkDir/temp$libraryname/1-overlapper/long_reads_part\$sa.sa $wrkDir/temp$libraryname/1-overlapper/correct_reads_part\$oneJobId.fasta \\\n";
         print F "          $wrkDir/temp$libraryname/1-overlapper/long_reads_part\$sa.fasta \\\n";
         print F "          -nproc $ovlThreads \\\n";
         if ($cutoffSpecified == 0) {
         print F "          -bestn $blasrThreshold \\\n";
         }
         print F "          $blasrOpts \\\n";
         print F "          -sam -out $wrkDir/temp$libraryname/1-overlapper/\$bat/\$job.sam \\\n";
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
         print F "          -x $wrk/temp$libraryname/1-overlapper/long.sa -f $wrk/temp$libraryname/1-overlapper/correct_reads_part\$oneJobId.fasta \\\n";
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
      } elsif (defined(getGlobal("mhap"))) {
         my $javaCmd = "/usr/bin/time $javaPath -server -Xmx" . $ovlMemory . "g -jar $MHAP_OVL FastAlignMain " . getGlobal("mhap") . " --min-store-length " . ($seedLength-1) . " --num-threads $ovlThreads" . " $ignore ";
         my $convertCmd = "awk -v REF_OFFSET=\$refOffset -v OFFSET=\$hashStart '{if (\$5 == 0 && \$9 == 0) { ORI=\"N\"; } if (\$5 == 0 && \$9 == 1) { ORI=\"I\"; } if (\$8 <= \$12 && (\$7-\$6) / \$8 > 0.9) print \$1+OFFSET+REF_OFFSET\"\\t\"\$2+OFFSET\"\\t\"ORI\"\\t\"(-1*\$10)\"\\t\"(\$12-\$11)\"\\t\"\$3/5\"\\t\"\$3/5; else if (\$8 > \$12 && (\$11 - \$10) / \$12 > 0.9) print \$1+OFFSET+REF_OFFSET\"\\t\"\$2+OFFSET\"\\t\"ORI\"\\t\"\$6\"\\t\"(-1*(\$8-\$7))\"\\t\"\$3/5\"\\t\"\$3/5; }' | \$bin/convertOverlap -ovl > $wrk/temp$libraryname/1-overlapper/\$bat/\$job.ovb";
         print F " startIndex=`perl -w -e \"use POSIX; print floor(\$start/$ovlRefBlockSize)\"`\n";
         print F " endIndex=`perl -w -e \"use POSIX; print ceil(\$end/$ovlRefBlockSize)\"`\n";
         if (defined(getGlobal("localStaging"))) {
            print F " if [ ! -e $wrkDir/temp$libraryname/1-overlapper/correct_reads_part\$sa.fasta ]; then\n";
            print F " echo \"Copying file $wrk/temp$libraryname/1-overlapper/correct_reads_part\$sa.fasta\"\n";
            print F " cp $wrk/temp$libraryname/1-overlapper/correct_reads_part\$sa.fasta $wrkDir/temp$libraryname/1-overlapper/correct_reads_part\$sa.fasta\n";
            print F "fi\n";
         }
         print F " rm -rf $wrkDir/temp$libraryname/1-overlapper/stream_\$jobid\n";
         print F " mkdir -p $wrkDir/temp$libraryname/1-overlapper/stream_\$jobid\n";
         print F " for i in \$(seq \$((\$startIndex+1)) \$endIndex); do\n";
         print F "    leadingI=`printf %06d \$i`\n";
         if (defined(getGlobal("localStaging"))) {
            print F "    cp $wrk/temp$libraryname/1-overlapper/correct_reads_part\$i.fasta $wrkDir/temp$libraryname/1-overlapper/stream_\$jobid/correct_reads_part\$leadingI.fasta\n";
         } elsif (!defined(getGlobal("mhapPrecompute"))) {
            print F "    ln -s $wrk/temp$libraryname/1-overlapper/correct_reads_part\$i.fasta $wrkDir/temp$libraryname/1-overlapper/stream_\$jobid/correct_reads_part\$leadingI.fasta\n";
         } else {
            print F "    ln -s $wrk/temp$libraryname/1-overlapper/correct_reads_part\$i.dat $wrkDir/temp$libraryname/1-overlapper/stream_\$jobid/correct_reads_part\$leadingI.dat\n";
         }
         print F " done\n";
         print F " if [ \$sa -eq \$numSA ]; then\n";
         print F "    $javaCmd -s $wrkDir/temp$libraryname/1-overlapper/correct_reads_part\$sa.$suffix 2> $wrk/temp$libraryname/1-overlapper/\$jobid.err | $convertCmd\n";
         print F " else\n";
         print F "    if [ \$start -eq \$hashEnd ]; then\n";
         print F "       $javaCmd -s $wrkDir/temp$libraryname/1-overlapper/correct_reads_part\$sa.$suffix -q $wrkDir/temp$libraryname/1-overlapper/stream_\$jobid 2> $wrk/temp$libraryname/1-overlapper/\$jobid.err | $convertCmd\n";
         print F "    else\n";
         print F "       $javaCmd --no-self -s $wrkDir/temp$libraryname/1-overlapper/correct_reads_part\$sa.$suffix -q $wrkDir/temp$libraryname/1-overlapper/stream_\$jobid 2> $wrk/temp$libraryname/1-overlapper/\$jobid.err | $convertCmd\n";
         print F "    fi\n";
         print F " fi\n";
         print F " rm -rf $wrkDir/temp$libraryname/1-overlapper/stream_\$jobid\n";
      }
      if (!defined(getGlobal("mhap"))) {
      print F "   if test -e $wrk/temp$libraryname/1-overlapper/long_reads_part\$sa.eidToIID ; then\n";
      print F "      \$bin/convertSamToCA \\\n";
      print F "           $wrkDir/temp$libraryname/1-overlapper/\$bat/\$job.$suffix $wrk/temp$libraryname/1-overlapper/long_reads_part\$sa.eidToIID,$wrk/temp$libraryname/1-overlapper/correct_reads_part\$oneJobId.eidToIID $wrk/temp$libraryname/1-overlapper/long_reads_part\$sa.iidToLen,$wrk/temp$libraryname/1-overlapper/correct_reads_part\$oneJobId.iidToLen \\\n";
      print F "           > $wrk/temp$libraryname/1-overlapper/\$bat/\$job.ovls 2> $wrk/temp$libraryname/1-overlapper/\$jobid.java.err\\\n";
      print F "               && touch $wrk/temp$libraryname/1-overlapper/\$jobid.java.success\n";
      print F "   else\n";
      print F "      \$bin/convertSamToCA \\\n";
      print F "           $wrkDir/temp$libraryname/1-overlapper/\$bat/\$job.$suffix $wrk/temp$libraryname/$asm.eidToIID $wrk/temp$libraryname/$asm.iidToLen \\\n";
      print F "           > $wrk/temp$libraryname/1-overlapper/\$bat/\$job.ovls 2> $wrk/temp$libraryname/1-overlapper/\$jobid.java.err\\\n";
      print F "               && touch $wrk/temp$libraryname/1-overlapper/\$jobid.java.success\n";
      print F "   fi\n";
      print F "   if test -e $wrk/temp$libraryname/1-overlapper/\$jobid.java.success ; then\n";
      print F "      echo SamToCA conversion completed.\n";
      print F "   else\n";
      print F "      echo SamToCA conversion failed.\n";
      print F "      tail $wrk/temp$libraryname/1-overlapper/\$jobid.java.err && exit\n";
      print F "   fi\n";
      print F "   \$bin/convertOverlap -ovl -i $wrk/temp$libraryname/1-overlapper/\$bat/\$job.ovls -o $wrk/temp$libraryname/1-overlapper/\$bat/\$job.ovb\n";
      print F "   rm -f $wrk/temp$libraryname/1-overlapper/\$bat/\$job.sam\n";
      }
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
            runCommand("$wrk/temp$libraryname/1-overlapper/$b/", "rm -f $j.ovls");
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
      runCommand("$wrk/temp$libraryname/1-overlapper/", "rm -f *.dat");
      # cleanup local storage if needed
      if (!-e "$wrk/temp$libraryname/1-overlapper/ovlclean.sh" && defined(getGlobal("localStaging"))) {
         open F, "> $wrk/temp$libraryname/1-overlapper/ovlclean.sh" or die ("can't open '$wrk/temp$libraryname/1-overlapper/ovlclean.sh'");
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
         print F "\n";
         my $wrkDir = $wrkDir = getGlobal("localStaging") . "/\$USER/";
         print F "rm -rf $wrkDir/temp$libraryname/1-overlapper/\n";
         close(F);
         chmod 0755, "$wrk/temp$libraryname/1-overlapper/ovlclean.sh";

         if ($submitToGrid == 1) {
            my $sgeName = "pBcR_ovlclean_$asm$sgeName";
            my $jobName = getGridArrayName($sgeName, $numOvlJobs);
            my $arrayOpt = getGridArrayOption($sgeName, $numOvlJobs);
           submitBatchJobs($wrk, $asm, "$submitCommand $sge $sgeOvl $nameOption \"$jobName\" $arrayOpt $outputOption /dev/null $wrk/temp$libraryname/1-overlapper/ovlclean.sh", $jobName);
         } else {
            for (my $i = 1; $i <= $numOvlJobs; $i++) {
               schedulerSubmit("$wrk/temp$libraryname/1-overlapper/ovlclean.sh $i");
            }
            schedulerSetNumberOfProcesses(getGlobal("ovlConcurrency"));
            schedulerFinish();
         }
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
   if ( -e "$wrk/temp$libraryname/$asm.layout.err") {
      # avoid infinite loops
      die "Error: found output from runCorrection. Stopping to avoid infinite loops. To try again, remove $asm.layout.*\n";
   }
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

    runCommand("$wrk/temp$libraryname", "rm -rf $asm.paired.ovlStore");
}

if (! -e "$wrk/temp$libraryname/runPartition.sh") {
   # convert QV to coverage for PBCNS
   # we use an ad-hoc conversion where QV 60 = coverage 8 (max) and QV 50 or below = coverage 1
   my $coverage = floor($QV - 50);
   if ($coverage < 1) { $coverage = 1; }
   if ($coverage > 8) { $coverage = 8; }

   my $cnsType = "-consensus pbdagcon";
   # priority is pbcns then falconcns
   if (defined(getGlobal("pbcns")) && getGlobal("pbcns") == 1) {
     setGlobal("falconcns", 0);
   }
   if (defined(getGlobal("falconcns")) && getGlobal("falconcns") == 1) {
     setGlobal("pbcns", 0);
     $cnsType = " -consensus falcon -pythonPath $FALCON";
     $coverage += 2;
   }
   my $falconcns = getGlobal("falconcns");
   my $pbcns = getGlobal("pbcns");

   my $threads = 1;
   if ((defined($pbcns) && $pbcns == 1) || (defined($falconcns) && $falconcns == 1)) {
      if ($submitToGrid == 1) {
         $threads = $consensusConcurrency;
      } else {
         $threads = 8;
         if ($consensusConcurrency < $threads) {
            $threads = $consensusConcurrency;
            $consensusConcurrency = 1;
         } else {
            $consensusConcurrency = int(int($consensusConcurrency) / int($threads));
         }
      }
   }

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
   my $wrkDir = $wrk;
   print F "if test -e $wrk/temp$libraryname/\$jobid.success ; then\n";
   print F "   echo \"Job previously completed successfully.\"\n";
   print F "else\n";
   print F "if test -e $wrk/temp$libraryname/\$jobid.lay.success ; then\n";
   print F "   echo \"Layout previously generated.\"\n";
   print F "else\n";
   print F "   \$bin/outputLayout \\\n";
   if (defined($longReads) && $longReads == 1) {
   print F "      -L \\\n";
   }
   print F "      -e $ovlErrorRate -M $maxUncorrectedGap \\\n";
   print F "      -i $wrk/temp$libraryname/$asm \\\n";
   print F "      -o $wrkDir/temp$libraryname/$asm \\\n";
   print F "      -p \$jobid \\\n";
   print F "      -l $length \\\n";
   print F "      $repeats \\\n";
   if (defined($pbcns) && $pbcns == 1) {
   print F "      -P \\\n";
   }
   if (defined($falconcns) && $falconcns == 1) {
   print F "     -F \\\n";
   }
   print F "      -G $wrk/temp$libraryname/$asm.gkpStore \\\n";
   if (defined($falconcns) && $falconcns == 1) {
   print F " 2> $wrk/temp$libraryname/\$jobid.lay.err | $FALCON/falcon_sense --max_n_read 500 --trim --min_idt 0.50 --output_multi --local_match_count_threshold 3 --min_cov $coverage --n_core $threads > $wrk/temp$libraryname/\$jobid.fasta  2> $wrk/temp$libraryname/\$jobid.cns.err\n";
   print F " if [ \$? -ge 0 ]; then\n";
   print F "   touch $wrk/temp$libraryname/\$jobid.success\n";
   print F " fi\n";
   } elsif (defined($pbcns) && $pbcns == 1) {
   print F " 2> $wrk/temp$libraryname/\$jobid.lay.err | \$bin/convertToPBCNS $cnsType -path $BLASR -output $wrk/temp$libraryname/\$jobid.fasta -prefix $wrkDir/temp$libraryname/\$jobid.tmp -length ";
   if (defined($longReads) && $longReads == 1) {
      print F "$length";   } else {      print F "50";
   }
   print F " -coverage $coverage -threads $threads > $wrk/temp$libraryname/\$jobid.err 2>&1 && touch $wrk/temp$libraryname/\$jobid.success\n";
# | $FALCON/pbdagcon -a -t 0 -m " . ((defined($longReads) && $longReads == 1)  ? $length : 50) . " -j $threads -c $coverage -b 2000  > $wrk/temp$libraryname/\$jobid.fasta  2> $wrk/temp$libraryname/\$jobid.cns.err\n";
   print F " if [ \$? -ge 0 ]; then\n";
   print F "   touch $wrk/temp$libraryname/\$jobid.success\n";
   print F " fi\n";
   } else {
   print F "      > $wrk/temp$libraryname/\$jobid.lay.err 2>&1 && touch $wrk/temp$libraryname/\$jobid.lay.success\n";
   print F "   if ! test -e $wrk/temp$libraryname/\$jobid.lay.success ; then\n";
   print F "      echo \"Error generating layout for job \$jobid\"\n";
   print F "      exit 1\n";
   print F "   fi\n";
   print F "   fi\n";
   print F "   \n";
   print F "   numLays=`cat $wrkDir/temp$libraryname/$asm" . ".\$jobid.lay |grep \"{LAY\" |wc -l`\n";
   print F "   if test \$numLays = 0 ; then\n";
   print F "      touch $wrk/temp$libraryname/\$jobid.fasta\n";
   print F "      touch $wrk/temp$libraryname/\$jobid.qual\n";
   print F "      touch $wrk/temp$libraryname/\$jobid.success\n";
   print F "   else\n";
   if (!defined $bankPath) {
   print F "      bankPath=\"$wrk/temp$libraryname\"\n";
   } elsif  (uc($bankPath) eq "SHARED") {
   print F "      bankPath=\"/dev/shm/temp$libraryname\"\n";
   print F "      laySize=`ls -la $wrk/temp$libraryname/$asm.\$jobid.lay |awk '{print \$5}'`\n";
   print F "      tmpFSSize=`df |grep /dev/shm |awk '{print \$(NF-2)}'`\n";
   if ($submitToGrid) {
   print F "      let requiredSize=\$laySize*2\n";
   } else {
   print F "      let requiredSize=\$laySize*2*$consensusConcurrency\n";
   }
   print F "      if test x\$tmpFSSize = x; then\n";
   print F "         tmpFSSize=0\n";
   print F "      fi\n";
   print F "      if [ \$requiredSize -lt \$tmpFSSize ]; then\n";
   print F "         mkdir -p /dev/shm/temp$libraryname\n";
   print F "      else\n";
   print F "         bankPath=\"$wrk/temp$libraryname\"\n";
   print F "      fi\n";
   } else {
   print F "      bankPath=\"$bankPath/temp$libraryname\"\n";
   print F "      mkdir -p \"$bankPath/temp$libraryname\"\n";
   }
   print F "      rm -rf \$bankPath/$asm" . ".bnk_partition\$jobid.bnk\n";
   print F "      $AMOS/bank-transact -b \$bankPath/$asm" . ".bnk_partition\$jobid.bnk -m $wrkDir/temp$libraryname/$asm.\$jobid" . ".lay -c > $wrk/temp$libraryname/bank-transact.\$jobid.err 2>&1\n";
   if (defined($shortReads) && $shortReads == 1) {
   print F "      $AMOS/make-consensus -e 0.03 -w 5 -x $wrk/temp$libraryname/\$jobid.excluded -B -b \$bankPath/" . $asm . ".bnk_partition\$jobid.bnk > $wrk/temp$libraryname/\$jobid.out 2>&1 && touch $wrk/temp$libraryname/\$jobid.success && rm $wrkDir/temp$libraryname/$asm.\$jobid.lay\n";
   } elsif (defined($longReads) && $longReads == 1) {
   print F "      $AMOS/make-consensus -e 0.30 -w 50 -o 5 -x $wrk/temp$libraryname/\$jobid.excluded -B -b \$bankPath/" . $asm . ".bnk_partition\$jobid.bnk > $wrk/temp$libraryname/\$jobid.out 2>&1 && touch $wrk/temp$libraryname/\$jobid.success && rm $wrkDir/temp$libraryname/$asm.\$jobid.lay\n";
   } else {
   print F "      $AMOS/make-consensus -x $wrk/temp$libraryname/\$jobid.excluded -B -b \$bankPath/" . $asm . ".bnk_partition\$jobid.bnk > $wrk/temp$libraryname/\$jobid.out 2>&1 && touch $wrk/temp$libraryname/\$jobid.success && rm $wrkDir/temp$libraryname/$asm.\$jobid.lay\n";
   }
   print F "      if test -e $wrk/temp$libraryname/\$jobid.success ; then\n";
   print F "         $AMOS/bank2fasta -e -q $wrk/temp$libraryname/\$jobid.qual -b \$bankPath/" . $asm . ".bnk_partition\$jobid.bnk > $wrk/temp$libraryname/\$jobid.fasta\n";
   print F "      else\n";
   print F "         rm -rf \$bankPath/$asm" . ".bnk_partition\$jobid.bnk\n";
   print F "         $AMOS/bank-transact -b \$bankPath/$asm" . ".bnk_partition\$jobid.bnk -m $wrkDir/temp$libraryname/$asm.\$jobid" . ".lay -c > $wrk/temp$libraryname/bank-transact.\$jobid.err 2>&1\n";
   if (defined($shortReads) && $shortReads == 1) {
   print F "         $AMOS/make-consensus -e 0.03 -w 5 -B -b \$bankPath/" . $asm . ".bnk_partition\$jobid.bnk > $wrk/temp$libraryname/\$jobid.out 2>&1 && touch $wrk/temp$libraryname/\$jobid.success && rm $wrkDir/temp$libraryname/$asm.\$jobid.lay\n";
   } elsif (defined($longReads) && $longReads == 1) {
   print F "         $AMOS/make-consensus -e 0.30 -w 50 -o 5 -B -b \$bankPath/" . $asm . ".bnk_partition\$jobid.bnk > $wrk/temp$libraryname/\$jobid.out 2>&1 && touch $wrk/temp$libraryname/\$jobid.success && rm $wrkDir/temp$libraryname/$asm.\$jobid.lay\n";
   } else {
   print F "         $AMOS/make-consensus -B -b \$bankPath/" . $asm . ".bnk_partition\$jobid.bnk > $wrk/temp$libraryname/\$jobid.out 2>&1 && touch $wrk/temp$libraryname/\$jobid.success && rm $wrkDir/temp$libraryname/$asm.\$jobid.lay\n";
   }
   print F "         $AMOS/bank2fasta -e -q $wrk/temp$libraryname/\$jobid.qual -b \$bankPath/" . $asm . ".bnk_partition\$jobid.bnk > $wrk/temp$libraryname/\$jobid.fasta\n";
   print F "      fi\n";
   print F "      if test -e $wrk/temp$libraryname/\$jobid.success ; then\n";
   print F "         $CA/trimFastqByQVWindow $wrk/temp$libraryname/\$jobid.fasta $wrk/temp$libraryname/\$jobid.qual $wrk/temp$libraryname/\$jobid.trim.fasta $wrk/temp$libraryname/\$jobid.trim.qual $QV $length > $wrk/temp$libraryname/\$jobid.trim.fastq\n";
   print F "         rm -rf \$bankPath/$asm" . ".bnk_partition\$jobid.bnk\n";
   print F "         rm $wrk/temp$libraryname/bank-transact.\$jobid.err\n";
   print F "         rm $wrk/temp$libraryname/\$jobid.excluded\n";
   print F "         rm $wrk/temp$libraryname/\$jobid.out\n";
   print F "      fi\n";
   }
   # cleanup intermediate files if we're successfull
   if ($cleanup == 1) {
   print F "      if test -e $wrk/temp$libraryname/\$jobid.success ; then\n";
   print F "         rm -f $wrk/temp$libraryname/$asm.\$jobid.shortmap*\n";
   print F "         rm -f $wrk/temp$libraryname/$asm.\$jobid.rank\n";
   print F "         rm -f $wrk/temp$libraryname/$asm.\$jobid.olaps\n";
   print F "      fi\n";
   }
   if ((defined($pbcns) && $pbcns == 1) || (defined($falconcns) && $falconcns == 1)) {
   print F "      if test -e $wrk/temp$libraryname/\$jobid.success ; then\n";
   print F "         cat $wrk/temp$libraryname/\$jobid.fasta | awk '{if(\$0~/>/){print;}else{l=length(\$0);q=\"\";while(l--){q=q \" 60\"}printf(\"%s\\n\",q)}}' > $wrk/temp$libraryname/\$jobid.qual\n";
   print F "         cat $wrk/temp$libraryname/\$jobid.fasta | awk '{if(\$0~/>/){sub(/>/,\"@\",\$0);print;}else{l=length(\$0);q=\"\";while(l--){q=q \"]\"}printf(\"%s\\n+\\n%s\\n\",\$0,q)}}' > $wrk/temp$libraryname/\$jobid.trim.fastq\n";
   print F "         ln -s $wrk/temp$libraryname/\$jobid.fasta $wrk/temp$libraryname/\$jobid.trim.fasta\n";
   print F "         ln -s $wrk/temp$libraryname/\$jobid.qual $wrk/temp$libraryname/\$jobid.trim.qual\n";
   print F "      fi\n";
   }
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

# check how much we split the sequences up
my %splitUIDs = ();
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

#save output from temporary directory
runCommand("$wrk/temp$libraryname", "cp $asm.layout.err  $wrk/$libraryname.correction.err");
if (-e "$wrk/temp$libraryname/$asm.layout.hist") { 
   runCommand("$wrk/temp$libraryname", "cp $asm.layout.hist  $wrk/$libraryname.correction.hist");
}
runCommand("$wrk/temp$libraryname", "cat `ls [0-9]*.fasta |grep trim |sort -T . -rnk1` > $wrk/$libraryname.fasta");
runCommand("$wrk/temp$libraryname", "cat `ls [0-9]*.qual |grep trim | sort -T . -rnk1` > $wrk/$libraryname.qual");
runCommand("$wrk/temp$libraryname", "cat `ls [0-9]*.fastq |grep trim | sort -T . -rnk1` > $wrk/$libraryname.fastq");
runCommand("$wrk", "$CA/fastqToCA -libraryname $libraryname -technology pacbio-corrected -type sanger -reads $wrk/$libraryname.fastq > $wrk/$libraryname.frg");

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
   runCommand("$wrk", "rm -rf temp$libraryname/$asm.[0-9]*.lay");
}

print STDERR "********* Finished correcting $totalBP bp (using $totalCorrectingWith bp).\n";

# now assemble 
assemble:
if (defined(getGlobal("assemble")) && getGlobal("assemble") == 1) {
   my $libraryname = getGlobal("libraryname");

   if (! -d "$wrk/$libraryname" && ! -e "$wrk/$libraryname.longest25.frg") { 
      my $specFile = "$wrk/temp$libraryname/$libraryname.spec";
      my $genomeSize = getGlobal("genomeSize");
      if (!defined($genomeSize) || $genomeSize == 0) { # guestimate genome size
         $genomeSize = 5000000;
      }
      print STDERR "********* Assembling corrected sequences.\n";

      my $asmCoverage = getGlobal("assembleCoverage");
      my $totalToAsm = int $asmCoverage * $genomeSize;
      runCommand("$wrk", "$CA/gatekeeper -T -F -o $libraryname.gkpStore $libraryname.frg");
      my $average = 0;
      my $total = 0;
      open(F, "$CA/gatekeeper -dumpfragments -tabular $libraryname.gkpStore |") or die("Couldn't open gatekeeper store", undef);
      while (<F>) {
         s/^\s+//;
         s/\s+$//;

         my @array = split '\s+';
         if ($array[9] >= $length) {
            $average += $array[9];
            $total++;
         }
     }
     close(F);
     my $totalBP = $average;
     $average /= $total;

     # should use a formula, maybe 40% rounded to next 500
     my $frgLen = 3000;
     if ($average < $frgLen) {
        $frgLen = $average / 2;
     }

     if (defined(getGlobal("falconcns")) && getGlobal("falconcns") == 1) {
        setGlobal("asmUtgErrorRate", getGlobal("asmUtgErrorRate") * $FALCON_ERATE_ADJUST);
        setGlobal("utgGraphErrorRate", getGlobal("utgGraphErrorRate") * $FALCON_ERATE_ADJUST);
     }

     my $ovlLen = 80;
     print STDERR "Assembling with average $average and using ovl is $ovlLen\n";
     runCommand("$wrk", "$CA/gatekeeper -dumpfrg -longestlength 0 $totalToAsm $libraryname.gkpStore > $libraryname.longest$asmCoverage.frg");
     runCommand("$wrk", "$CA/gatekeeper -dumpfasta $wrk/$libraryname.longest$asmCoverage -longestlength 0 $totalToAsm $libraryname.gkpStore");
     runCommand("$wrk", "cat $wrk/$libraryname.longest$asmCoverage.fasta |awk '{if (match(\$1, \">\")) { print \$1; } else { print \$0; } }' > tmp.fasta");
     runCommand("$wrk", "mv tmp.fasta $wrk/$libraryname.longest$asmCoverage.fasta");
     runCommand("$wrk", "rm -rf $libraryname.gkpStore*");

     # don't assemble if we don't have enough data
     if ($totalBP < $MIN_COVERAGE_TO_ASM * $genomeSize) {
        print STDERR "Error: after correction only " . ($totalBP / $genomeSize) . "X for genome $genomeSize. Not performing automated assembly\n";
     } else {
        my $ovlPartCmd = "";
        my $ovlRefBlockSize    = getGlobal("ovlRefBlockSize");
        my $ovlRefBlockLength  = getGlobal("ovlRefBlockLength");
        if ($ovlRefBlockLength > 0) {
           $ovlPartCmd = "ovlRefBlockLength=$ovlRefBlockLength ovlRefBlockSize=0";
        } elsif ($ovlRefBlockSize > 0) {
           $ovlPartCmd = "ovlRefBlockSize=$ovlRefBlockSize ovlRefBlockLength=0";
        }
        $cmd  =    "-s $specFile ";
        $cmd .=    "-p $asm -d $libraryname $ovlPartCmd ";
        $cmd .=    "useGrid=" .getGlobal("useGrid") . " scriptOnGrid=" . getGlobal("scriptOnGrid") . " unitigger=" . getGlobal("unitigger") . " ";
        $cmd .=    "ovlErrorRate=" . getGlobal("asmOvlErrorRate") . " utgErrorRate=" . getGlobal("asmUtgErrorRate") . " cgwErrorRate=" . getGlobal("asmCgwErrorRate") . " cnsErrorRate=" . getGlobal("asmCnsErrorRate") . " ";
        $cmd .=    "utgGraphErrorLimit=" . getGlobal("utgGraphErrorLimit") . " utgGraphErrorRate=" . getGlobal("utgGraphErrorRate") . " utgMergeErrorLimit=" . getGlobal("utgMergeErrorLimit") . " utgMergeErrorRate=" . getGlobal("utgMergeErrorRate") . " ";
        $cmd .=    "frgCorrBatchSize=100000 doOverlapBasedTrimming=" . getGlobal("asmOBT") . " obtErrorRate=" . getGlobal("asmObtErrorRate") . " obtErrorLimit=" . getGlobal("asmObtErrorLimit") . " frgMinLen=$frgLen ovlMinLen=$ovlLen ";
        $cmd .=    "consensus=" . getGlobal("asmCns") . " merSize=" . getGlobal("asmMerSize") . " cnsMaxCoverage=1 cnsReuseUnitigs=1 ";
        $cmd .=    "sgeName=\"" . getGlobal("sgeName") . "\" " if defined(getGlobal("sgeName"));
        $cmd .=    "sgePropagateHold=\"pBcR_$asm$sgeName\" ";
        $cmd .=    " $libraryname.longest$asmCoverage.frg ";
        if ($submitToGrid == 1) {
           submitRunCAHelper("$wrk", $asm, $cmd, "runCA_asm_$asm$sgeName", 0);
        } else {
           runCommand("$wrk", "$CA/runCA $cmd");
        }
     }
  }
}
