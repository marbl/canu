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
my @nonCAOptions = ("genomeSize", "shortReads", "libraryname", "specFile", "length", "coverage", "maxUncorrectedGap", "threads", "repeats", "fastqFile", "partitions", "sumitToGrid", "sgeCorrection", "consensusConcurrency", "cleanup");

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
}

sub setDefaults() {
	$global{"shortReads"} = undef;
	$global{"libraryname"} = undef;
	$global{"specFile"} = undef;
	$global{"length"} = 500;
	$global{"coverage"} = 0;
        $global{"genomeSize"} = 0;
	$global{"maxUncorrectedGap"} = 0;
	$global{"threads"} = 1;
	$global{"repeats"} = "";
	$global{"fastqFile"} = undef;
	$global{"partitions"} = 1;
	$global{"sge"} = undef;
	$global{"submitToGrid"} = 0;
	$global{"sgeCorrection"} = undef;
	$global{"consensusConcurrency"} = 8;
	$global{"cleanup"} = 1;
}

sub updateSpecFile($$) {
	my $specFile = shift @_;
	my $outFile = shift @_;
    open(F, "< $specFile") or die("Couldn't open '$specFile'", undef);
    open(W, "> $outFile") or die("Couldn't open '$outFile'", undef);

	my $print = 1;
    while (<F>) {
		s/^\s+//;
        s/\s+$//;
        
    	# skip comments
		next if (m/^\s*\#/);
        next if (m/^\s*$/);

    	$print = 1;        
    	foreach my $key (@nonCAOptions) {
    		if (index($_, $key) != -1) {
    			$print = 0;
    			last;
    		}
    	}
    	if ($print == 1) {
    		print W "$_\n";
    	}
    }
    close(F);
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

sub getBinDirectory () {
    my $installDir;

    ###
    ### CODE DUPLICATION WITH getBinDirectoryShellCode
    ###

    #  Assume the current binary path is the path to the global CA
    #  install directory.

    #  CODE DUPLICATION!!!

    my @t = split '/', "$FindBin::RealBin";
    pop @t;                      #  bin
    pop @t;                      #  arch, e.g., FreeBSD-amd64
    my $installDir = join '/', @t;  #  path to the assembler
    #  CODE DUPLICATION!!!

    #  Guess what platform we are currently running on.

    my $syst = `uname -s`;    chomp $syst;  #  OS implementation
    my $arch = `uname -m`;    chomp $arch;  #  Hardware platform
    my $name = `uname -n`;    chomp $name;  #  Name of the system

    $arch = "amd64"  if ($arch eq "x86_64");
    $arch = "ppc"    if ($arch eq "Power Macintosh");

    my $path = "$installDir/$syst-$arch/bin";
    return($path);
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
my $MIN_FILES_WITHOUT_PARTITIONS = 20;
my $REPEAT_MULTIPLIER = 10;

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
          push @fragFiles, $arg;
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
   $err++;
}
setParametersFromCommandLine(@specOpts);

# finally, command-line parameters take precedence
while (scalar(@cmdArgs) > 0) {
    my $arg = shift @cmdArgs;

    if ($arg eq "-shortReads") {
        setGlobal("shortReads", 1);

    } elsif ($arg eq "-coverage") {
        setGlobal("coverage", shift @cmdArgs);
        if (getGlobal("coverage") < 0) { setGlobal("coverage", 0); }

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
    
if (($err) || (scalar(@fragFiles) == 0) || (!defined(getGlobal("fastqFile"))) || (!defined(getGlobal("specFile"))) || (!defined(getGlobal("libraryname")))) {
   print STDERR "No frag files specified. Please specify at least one frag file containing high-accuracy data for correction.\n" if (scalar(@fragFiles) == 0);
   print STDERR "No fastq file specified. Please specify a fastq file containing PacBio sequences for correction using the -fastq parameter.\n" if (!defined(getGlobal("fastqFile")));
   print STDERR "No spec file defined. Please specify a spec file using the -s option.\n" if (!defined(getGlobal("specFile")));
   print STDERR "No library name provided. Please specify a name using the -library option.\n" if (!defined(getGlobal("libraryname")));  

    print STDERR "usage: $0 [options] -s spec.file -fastq fastqfile <frg>\n";
    print STDERR "  -length                  Minimum length to keep.\n";
    print STDERR "  -partitions              Number of partitions for consensus\n";
    print STDERR "  -sge                     Submit consensus jobs to the grid\n";
    print STDERR "  -sgeCorrection           Parameters for the correction step for the grid. This should match the threads specified below, for example by using -pe threaded\n";
    print STDERR "  -libraryname libraryname Name of the library; freeformat text.\n";
    print STDERR "  -threads threads         Number of threads to use for correction.\n";
    print STDERR "  -shortReads              Use if the sequences for correction are 100bp or shorter.\n";

    print STDERR "  -coverage		     Specify the pacBio coverage (integer) instead of automatically estimating.\n";
    print STDERR "  -genomeSize		     Specify the approximate genome size. This overwrites the coverage parameter above and will be used to compute estimated coverage instead of automatically estimating.\n\n";

    print STDERR "\nAdvanced options (EXPERT):\n";
    print STDERR "  -maxGap		     The maximum uncorrected PacBio gap that will be allowed. When there is no short-read coverage for a region, by default the pipeline will split a PacBio sequence. This option allows a number of PacBio sequences without short-read coverage to remain. For example, specifying 50, will mean 50bp can have no short-read coverage without splitting the PacBio sequence. Warning: this will allow more sequences that went through the SMRTportal to not be fixed.\n";
    exit(1);
}

#check for valid parameters for requested partitions and threads
$MIN_FILES_WITHOUT_PARTITIONS += getGlobal("threads");
my $limit = 1024;
my $sysLimit = `ulimit -Sn`;
chomp($sysLimit);
if (defined($sysLimit)) {
   $limit = $sysLimit;
}
if ($limit - $MIN_FILES_WITHOUT_PARTITIONS <= getGlobal("partitions")) {
   setGlobal("partitions", $limit - $MIN_FILES_WITHOUT_PARTITIONS);
   if (getGlobal("threads") > getGlobal("partitions")) { setGlobal("threads", getGlobal("partitions") - 1); }
   print STDERR "Warning: file handle limit of $limit prevents using requested partitions. Reset partitions to " . getGlobal("partitions") . ". If you want more partitions, reset the limit and try again.\n";
}
if (getGlobal("partitions") <= getGlobal("threads")) {
   setGlobal("partitions", getGlobal("threads") + 1);
   print STDERR "Warning: number of partitions should be > # threads. Adjusted partitions to be ". getGlobal("partitions") . "\n";
}

print STDOUT "Running with " . getGlobal("threads") . " threads and " . getGlobal("partitions") . " partitions\n";

my $CA = getBinDirectory();
my $AMOS = "$CA/../../../AMOS/bin/";
my $wrk = makeAbsolute("");
my $asm = "asm";
my $caSGE  = getGlobal("sge");

if (defined($caSGE)) {
   $caSGE = "\"" .$caSGE . " -sync y\" sgePropagateHold=corAsm";
} else {
   $caSGE = "sge=\"" . " -sync y\" sgePropagateHold=corAsm";
}
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
print STDERR "********* Starting correction...\n CA: $CA\nAMOS:$AMOS \n";
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
my $coverage = getGlobal("coverage");
my $genomeSize = getGlobal("genomeSize");
my $shortReads = getGlobal("shortReads");
my $length = getGlobal("length");
my $repeats = getGlobal("repeats");

my $threads = getGlobal("threads");
my $partitions = getGlobal("partitions");

my $submitToGrid = getGlobal("submitToGrid");
my $sge = getGlobal("sge");
my $sgeCorrection = getGlobal("sgeCorrection");
my $consensusConcurrency = getGlobal("consensusConcurrency");

my $cleanup = getGlobal("cleanup");

my $cmd = "";

if ($coverage != 0 && $genomeSize != 0) {
   print STDERR "Warning, both coverage and genome size is set. Coverage $coverage will be overwritten by genome size\n";
}

if (! -e "temp$libraryname") {
   runCommand("$wrk", "mkdir temp$libraryname");
}
# generate the ca spec file, since we support additional options in the spec file, we need to strip those out before passing it to ca
updateSpecFile($specFile, "$wrk/temp$libraryname/$libraryname.spec");
$specFile = "$wrk/temp$libraryname/$libraryname.spec";

runCommand($wrk, "$CA/fastqToCA -libraryname PacBio -type sanger -innie -technology pacbio-long -reads " . makeAbsolute($fastqFile) . " > $wrk/temp$libraryname/$libraryname.frg"); 
runCommand($wrk, "$CA/runCA -s $specFile -p $asm -d temp$libraryname $caSGE stopAfter=initialStoreBuilding @fragFiles $wrk/temp$libraryname/$libraryname.frg");

# make assumption that we correct using all libraries preceeding pacbio
# figure out what number of libs we have and what lib is pacbio
my $numLib = `$CA/gatekeeper -dumpinfo temp$libraryname/$asm.gkpStore | grep LIB |awk '{print \$1}'`;
chomp($numLib);

my $minCorrectLib = 0;
my $maxCorrectLib = 0;
my $libToCorrect = 0;
for (my $i = 1; $i <= $numLib; $i++) {
   if (system("$CA/gatekeeper -isfeatureset $i doConsensusCorrection temp$libraryname/$asm.gkpStore") == 0) {
      $libToCorrect = $i;
    } else {
      if ($minCorrectLib == 0) { $minCorrectLib = $i; }
      $maxCorrectLib = $i;
   }
}

# run correction up thorough meryl
runCommand($wrk, "$CA/runCA -s $specFile -p $asm -d temp$libraryname $caSGE stopAfter=meryl");

# set the meryl threshold based on the genome size and coverage (if specified)
if ($genomeSize != 0) {
   my $totalBP = 0;

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

   if ($totalBP > 0) {
      $coverage = ceil($totalBP / $genomeSize);
      setGlobal("coverage", $coverage);
   }
}
# set the meryl threshold if necessary
if (getGlobal("ovlMerThreshold") == "auto") {
   # no threshold specified, check if the chosen one is OK
   my $autoSetThreshold = `cat temp$libraryname/0-mercounts/*estMerThresh.out`;
   chomp($autoSetThreshold); 

   if ($autoSetThreshold < ($coverage * $REPEAT_MULTIPLIER)) {
      setGlobal("ovlMerThreshold", $coverage * $REPEAT_MULTIPLIER);
      unlink("temp$libraryname/0-mercounts/$asm.nmers.ovl.fasta");
      print STDERR "Resetting from auto threshold of $autoSetThreshold to be " . ($coverage * 10) . "\n";
   }
}

# now run the correction
my $ovlThreshold = getGlobal("ovlMerThreshold");
$cmd  = "$CA/runCA ";
$cmd .=    "-s $specFile ";
$cmd .=    "-p $asm -d temp$libraryname ";
$cmd .=    "ovlMerThreshold=$ovlThreshold ";
$cmd .=    "ovlHashLibrary=$libToCorrect ";
$cmd .=    "ovlRefLibrary=$minCorrectLib-$maxCorrectLib ";
$cmd .=    "ovlCheckLibrary=1 ";
$cmd .=    "obtHashLibrary=$minCorrectLib-$maxCorrectLib ";
$cmd .=    "obtRefLibrary=$minCorrectLib-$maxCorrectLib ";
$cmd .=    "obtCheckLibrary=0 ";
$cmd .=   "$caSGE stopAfter=overlapper";
runCommand($wrk, $cmd);

if (! -e "$wrk/temp$libraryname/$asm.layout.success") {
   open F, "> $wrk/temp$libraryname/runCorrection.sh" or die ("can't open '$wrk/temp$libraryname/runCorrection.sh'");
   print F "#!" . "/bin/sh\n";
   print F "\n";
   print F " if test -e $wrk/temp$libraryname/$asm.layout.success; then\n";
   print F "    echo Job previously completed successfully.\n";
   print F " else\n";
   print F "   $CA/correctPacBio \\\n";
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
      runCommand("$wrk/temp$libraryname", "qsub $sge $sgeCorrection -sync y -cwd -N correct_$asm -j y -o /dev/null $wrk/temp$libraryname/runCorrection.sh");
   } else {
      runCommand("$wrk/temp$libraryname", "$wrk/temp$libraryname/runCorrection.sh");
   }
}

if (! -e "$wrk/temp$libraryname/runPartition.sh") {
   open F, "> $wrk/temp$libraryname/runPartition.sh" or die ("can't open '$wrk/temp$libraryname/runPartition.sh'");
   print F "#!" . "/bin/sh\n";
   print F "\n";
   print F "jobid=\$SGE_TASK_ID\n";
   print F "if test x\$jobid = x -o x\$jobid = xundefined; then\n";
   print F "jobid=\$1\n";
   print F "fi\n";
   print F "\n";
   print F "if test x\$jobid = x; then\n";
   print F "  echo Error: I need SGE_TASK_ID set, or a job index on the command line\n";
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
   print F "      $AMOS/bank-transact -b $wrk/temp$libraryname/$asm" . ".bnk_partition\$jobid.bnk -m $wrk/temp$libraryname/$asm.\$jobid" . ".lay -c > $wrk/temp$libraryname/bank-transact.\$jobid.err 2>&1\n";
   if (defined($shortReads) && $shortReads == 1) {
   print F "      $AMOS/make-consensus -e 0.03 -w 5 -x $wrk/temp$libraryname/\$jobid.excluded -B -b $wrk/temp$libraryname/" . $asm . ".bnk_partition\$jobid.bnk > $wrk/temp$libraryname/\$jobid.out 2>&1 && touch $wrk/temp$libraryname/\$jobid.success\n";
   } else {
   print F "      $AMOS/make-consensus -x $wrk/temp$libraryname/\$jobid.excluded -B -b $wrk/temp$libraryname/" . $asm . ".bnk_partition\$jobid.bnk > $wrk/temp$libraryname/\$jobid.out 2>&1 && touch $wrk/temp$libraryname/\$jobid.success\n";
   }
   print F "      if test -e $wrk/temp$libraryname/\$jobid.success ; then\n";
   print F "         $AMOS/bank2fasta -e -q $wrk/temp$libraryname/\$jobid.qual -b $wrk/temp$libraryname/" . $asm . ".bnk_partition\$jobid.bnk > $wrk/temp$libraryname/\$jobid.fasta\n";
   print F "      else\n";
   print F "         rm -rf $wrk/temp$libraryname/$asm" . ".bnk_partition\$jobid.bnk\n";
   print F "         $AMOS/bank-transact -b $wrk/temp$libraryname/$asm" . ".bnk_partition\$jobid.bnk -m $wrk/temp$libraryname/$asm.\$jobid" . ".lay -c > $wrk/temp$libraryname/bank-transact.\$jobid.err 2>&1\n";
   print F "         $AMOS/make-consensus -B -b $wrk/temp$libraryname/" . $asm . ".bnk_partition\$jobid.bnk > $wrk/temp$libraryname/\$jobid.out 2>&1 && touch $wrk/temp$libraryname/\$jobid.success\n";
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
      runCommand("$wrk/temp$libraryname", "qsub $sge -sync y -cwd -N utg_$asm -t 1-$partitions -j y -o /dev/null $wrk/temp$libraryname/runPartition.sh");
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

runCommand("$wrk/temp$libraryname", "cat `ls [1234567890]*.fasta |sort -rnk1` > corrected.fasta");
runCommand("$wrk/temp$libraryname", "cat `ls [1234567890]*.qual |sort -rnk1` > corrected.qual");
runCommand("$wrk", "$CA/convert-fasta-to-v2.pl -pacbio -s $wrk/temp$libraryname/corrected.fasta -q $wrk/temp$libraryname/corrected.qual -l $libraryname > $wrk/$libraryname.frg");
runCommand("$wrk/temp$libraryname", "cp corrected.fasta $wrk/$libraryname.fasta");
runCommand("$wrk/temp$libraryname", "cp corrected.qual  $wrk/$libraryname.qual");

# finally clean up the assembly directory
if ($cleanup == 1) {
   runCommand("$wrk", "rm -rf temp$libraryname");
}
