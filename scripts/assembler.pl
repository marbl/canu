#!/usr/bin/perl -w
#
###########################################################################
#
# This file is part of Celera Assembler, a software program that 
# assembles whole-genome shotgun reads into contigs and scaffolds.
# Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received (LICENSE.txt) a copy of the GNU General Public 
# License along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
###########################################################################
#

print "The Celera Whole Genome Shotgun Assembler.\n";
#######################################################################
#
# Script for running the Celera Whole Genome Shotgun Assembler.
#
# CODE STRUCTURE
# + The main() subroutine is located at the bottom.
# + The usage() subroutine runs if script is invoked without parameters.
# + Utility subroutines are located in a separate perl Package.
#
# YOU CAN STOP & RESTART THE SCRIPT
# + You can re-run this script, picking up where a previous run stopped. 
# All steps are numbered. The output (stdout) contains the number 
# of each step. The command line accepts optional start and stop numbers.
# + Numbers advance by 100s for main components and by 1s otherwise.
# For instance, Unitigger is number 500, so individual
# unitigging steps are numbered 501, 502, 503, etc.
# + Every execution step is embedded in a test call to shouldExecute().
# This increments the step counter and tests for within command-line params.
# This facilitates restarting the script. Please add any other tests 
# (like "if LSF") after the call to shouldExecute(), so your test doesn't
# short-circuit the step counter increment.
#
# THE SCRIPT LOGS ACTIVITY & ERRORS
# + The ordinary output (stdout) should be very helpful. For every step,
# the output contains the step number, an English description,
# and the actual command line. 
# + In the event of an error, the output (stdout) should
# contain everything that can be trapped in perl, including 
# stack trace and child process exit status. 
# + Error messages (stderr) from child processes get redirected 
# to a <prefix>.stderr.log file.
#
# HISTORICAL NOTE
# By the way, this script does NOT depend on the environment variables
# $AS_ROOT or $AS_BIN.
#
# Blame: Jason Miller, Clark Mobarry, Randall Bolanos, Karin Remington.
#######################################################################


#----------------------------------------------------------------
# ----------------------- START: LIBRARIES ----------------------
#----------------------------------------------------------------
use strict;   # force variable declartion
use Getopt::Long;   # command line options longer than 1 letter
use diagnostics;    # print stack trace on error
use sigtrap;  # signal handler prints stack trace too
use English;   # full length variable names for perl globals
use Cwd;    # Current working directory

# Add the $AS_ROOT/scripts directory to the library search path @INC.
# Other ways you could help perl find the module are:
#    1) Use command line options -I<dir> -M<module>
#    2) Add those command line options to the first line of this file.
#    3) Add the scripts directory to your PERL5LIB environment variable.
#    4) Hard code "use lib 'myroot/scripts/';" into this script.
#    5) Link to the script and the module from your current directory.
#    6) Copy the script and the module to your current directory.
use FindBin;   # Find directory where this script lives...
use lib "$FindBin::Bin";   # ... and add it to the @INC path.

# Import the module of utilities written for this script.
# The module Assembler lives in $AS_ROOT/scripts/Assembler.pm.
use Assembler;   # Package of utilities written for this script.
#----------------------------------------------------------------
# ----------------------- END: LIBRARIES ------------------------
#----------------------------------------------------------------



#----------------------------------------------------------------
# ----------------------- START: GLOBALS ------------------------
#----------------------------------------------------------------
# !!!!!
# PLEASE KEEP THIS SCRIPT MOSTLY FREE OF GLOBAL VARIABLES!
# PLEASE USE setGlobal() INSTEAD.
# !!!!!
# Temporary: Global parameter switches used by all assembler components.
my ($OUTPUT_MODE,$FILE_CLOBBER);
#----------------------------------------------------------------
# ----------------------- END: GLOBALS --------------------------
#----------------------------------------------------------------



#######################################################################
# Start assembler subroutines.
#######################################################################

#----------------------------------------------------------------
# Run the Gatekeeper.
# This validates inputs while building the Gatekeeper store.
#----------------------------------------------------------------
sub run_gatekeeper {
    my ($nextStep) = @_;
    &Assembler::setNextStep($nextStep);
    my $commandLine;
    my $descriptionLine;
    my $prefix = &Assembler::getGlobal("prefix");

    my $useCountMessages = &Assembler::getGlobal("use countmessages");

    if (&Assembler::shouldExecute() && $useCountMessages) {
	$descriptionLine = "Count fragments in this batch.";
        $commandLine="countmessages < $prefix.frg > $prefix.cnt";
	&Assembler::runLocal($descriptionLine,$commandLine);
    }
    if (&Assembler::shouldExecute() && $useCountMessages) {
	# This step is informational.
	$descriptionLine = "Display fragment count.";
	$commandLine = "cat $prefix.cnt";
	&Assembler::runCommand($descriptionLine,$commandLine);
    }

    # Temporary quick fix.
    # Later, change this to use clobber mode like everything else?
    my $OLD_CREATE_MODE="-f";

    my $GATEKEEPER_CMD=
    "gatekeeper -X -C -N -Q $OUTPUT_MODE $OLD_CREATE_MODE";
    
    if (&Assembler::shouldExecute()) {
	$descriptionLine = "Gatekeeper will verify inputs. " .
            "The gatekeeper store is $prefix.gkpStore, a directory.";
        $commandLine="$GATEKEEPER_CMD $prefix.gkpStore $prefix.frg";
	&Assembler::runLocal($descriptionLine,$commandLine);
    }
}
    
#----------------------------------------------------------------
# Screener. Mask fragments matching ubiquitous repeats and contaminant.
# This step was discontinued in 2001.
# Parameter 1: First number for numbering each step (for the logs).
#----------------------------------------------------------------
sub run_screener {
    my ($nextStep) = @_;
    &Assembler::setNextStep($nextStep);
    my $commandLine;
    my $descriptionLine;
    my $prefix = &Assembler::getGlobal("prefix");
    my ($repeat_lib) = &Assembler::getGlobal("repeat library");

    my $URCSCREENER_CMD=
	"urc_screener -r $FILE_CLOBBER -s $OUTPUT_MODE $repeat_lib";
    
    my $useScreener = &Assembler::getGlobal("use screener");

    if (&Assembler::shouldExecute() && $useScreener) {
	$descriptionLine = "Screener. This will mask repeats. " .
	    "The screen library is $repeat_lib.";
	$commandLine="$URCSCREENER_CMD $prefix.inp";
	&Assembler::runLocal($descriptionLine,$commandLine);
    }
}
    
#----------------------------------------------------------------
# Run Overlapper.
# Build a collection of k-mer matches between fragments.
# Parameter 1: First number for numbering each step (for the logs).
#----------------------------------------------------------------
sub run_overlapper {
    my ($nextStep) = @_;
    &Assembler::setNextStep($nextStep);
    my $commandLine;
    my $descriptionLine;
    my $prefix = &Assembler::getGlobal("prefix");
     
    # Temporary quick fix.
    # This program cannot handle the options -c -f,
    # which should be universal for all assembler components.
#    my $OLD_CREATE_MODE="-f";

    # Whether to run the fragment correction algorithm.
    my $useFragmentCorrection = &Assembler::getGlobal("use fragment correction");

    # Populator runs here <done manually, add to this script later>
    # Populates frag store

#    my $OVERLAP_CMD1= "overlap -n $OUTPUT_MODE $OLD_CREATE_MODE";
    my $OVERLAP_CMD1= "PopulateFragStore $OUTPUT_MODE ";
     
    my $OVERLAP_CMD2=
	"overlap $OUTPUT_MODE -w -h 1- -r 1- "; 

    my ($useScreener , $INPUT_FILENAME);
    $useScreener = &Assembler::getGlobal("use screener");
    if ($useScreener) {
	# If screener ran, read its output file.
	$INPUT_FILENAME = "$prefix.urc" ;
    } else {
	# If screener didn't run, read its input file.
	$INPUT_FILENAME = "$prefix.inp" ;
    }	

    # This step uses overlapper to populate the frag store.
    # On large data sets, we use the PopulateFragStore program
    # just so we can make incremental backups.
    # Later, fix this script to run that.

    if (&Assembler::shouldExecute()) {
	$descriptionLine = "Overlapper makes a fragment store ";
        $commandLine="$OVERLAP_CMD1 -c -f -o $prefix.frgStore $INPUT_FILENAME";
#        $commandLine="$OVERLAP_CMD1 $prefix.frgStore $INPUT_FILENAME";
	&Assembler::runLocal($descriptionLine,$commandLine);
    }    

    # A separate overlapping step for the whole range of fragments.

    if (&Assembler::shouldExecute()) {
	$descriptionLine = "Overlapper compares fragments, all against all. " .
	    "The fragment store is in the $prefix.frgStore directory.";
        $commandLine="$OVERLAP_CMD2 -o $prefix.ovl $prefix.frgStore";
	&Assembler::runLocal($descriptionLine,$commandLine);
    }    

    # Actually, the following commands are for LSF only!
    # Fix this later.
    # NO! We always want to make an overlap store. -- CMM.
    
    if (&Assembler::shouldExecute()) {
	$descriptionLine = "Make overlap store.";
        $commandLine="grow-olap-store " .
	    "-cf -o $prefix.ovlStore $prefix.ovl";
	&Assembler::runLocal($descriptionLine,$commandLine);
    }

    # if fragment correction and lsf
    # then run correct-script here
    # <done manually now, add this to script later>
}

#----------------------------------------------------------------
# Run Unitigger.
# Build maximal contigs containing no contradicted elements.
# Parameter 1: First number for numbering each step (for the logs).
#----------------------------------------------------------------
sub run_unitigger {
    my ($nextStep) = @_;
    &Assembler::setNextStep($nextStep);
    my $commandLine;
    my $descriptionLine;
    my $prefix = &Assembler::getGlobal("prefix");

    my $smooth_spurs = 1;
    my $smooth_bubbles = 1;

    my $num_crappies = 0;
    my $fgbStoreHold = "$prefix.tmp.fgbStore";

    my $overlapErrorThreshold = "";
    my $useFragmentCorrection = 
	&Assembler::getGlobal("use fragment correction");
    if ($useFragmentCorrection) {
	# -e means ignore overlaps with error rates above this value.
	$overlapErrorThreshold = "-e 0.015";
    }

    # -x sets dovetail double-sided degree threshold.
    # -z sets containment degree threshold.
    my $degreeThresholds = "-x 1 -z 10";

    my $ReaperOptions = "$degreeThresholds " .
	"$overlapErrorThreshold " .
	    " -d 0 -M 0 ";
    my $FGBOptions = " -d 1 -M 1 ";    

    my $FGB_LOGGING = " 1>> $prefix.fgb.log ";
    my $ANAL_OPT = "-A 2 ";  
    my $U_UNITIG_ASTAT = "-j " . 
	&Assembler::getGlobal("U-unitig A-statistic");
    
    my $REAPER_CMD="fgb $OUTPUT_MODE $ANAL_OPT $FILE_CLOBBER ".
	"$ReaperOptions ";
    my $FGB_CMD =  "fgb $OUTPUT_MODE $ANAL_OPT $FILE_CLOBBER ".
	"$FGBOptions $FGB_LOGGING";
    my $CGB_CMD =  "cgb $OUTPUT_MODE $ANAL_OPT";

    # Convert a newer data format to an older format.
    # This step is for compatibility with older code.
    # Note FragStore needs changes before fgb can read it.
    # Known bug: looses coordinates from frag source field.

    if (&Assembler::shouldExecute()) {
	$descriptionLine = "Convert FragStore to OFG messages."; 
        $commandLine = "make_OFG_from_FragStore " .
	    "$prefix.frgStore > $prefix.ofg";
	&Assembler::runLocal($descriptionLine,$commandLine);
    }    

    if (&Assembler::shouldExecute()) {
	$descriptionLine = "Reaper. Filters out many overlaps, ".
	    "and labels some fragments spurs or contained."; 
        $commandLine = "$REAPER_CMD " .
	    "-I $prefix.ovlStore -Q 0 -c -o $prefix.fgbStore $prefix.ofg";
	&Assembler::runLocal($descriptionLine,$commandLine);
    }    
    
    if (&Assembler::shouldExecute()) {
	$descriptionLine = "FGB, after reaper."; 
        $commandLine="$FGB_CMD " .
	    " -a -i $prefix.fgbStore " ;
	&Assembler::runLocal($descriptionLine,$commandLine);
    }    

    if (&Assembler::shouldExecute() && $smooth_spurs) {
	$descriptionLine = "CGB. Spur smoothing."; 
        $commandLine = "$CGB_CMD -s -b 0 $U_UNITIG_ASTAT  " . 
	    " -U 0 $prefix.frgStore $prefix.fgbStore ";
	&Assembler::runLocal($descriptionLine,$commandLine);
    }    

    # This next bit is ugly. Should find a better way to count.
    if (&Assembler::shouldExecute()) {
	print "Start steps that clean up crappies.\n";
	$descriptionLine = "Count crappies.";
	$commandLine = 
	    "gawk '{ if(\$0 !~ /^#/) print \$0;}' " .
		"$prefix.cgb_crappies | wc -l";
	$num_crappies = 0 +
	    &Assembler::readChildProcess
		($descriptionLine,$commandLine);
	print "There are $num_crappies crappies.\n";
    }
    if (&Assembler::shouldExecute() && $num_crappies>0) {
	$descriptionLine = "Repair breakers."; 
        $commandLine = "repair_breakers -S $prefix.cgb_crappies " . 
	    "-F $prefix.frgStore -O $prefix.breaker_overlaps.ovl " .
		"$prefix.cgb";
	&Assembler::runLocal($descriptionLine,$commandLine);
    }
    if (&Assembler::shouldExecute() && $num_crappies>0) {
	$descriptionLine = "FGB, after repair breakers."; 
        $commandLine="$FGB_CMD " .
	    " -a -i $prefix.fgbStore " .
		    "$prefix.breaker_overlaps.ovl ";
	&Assembler::runLocal($descriptionLine,$commandLine);
    }    
    if (&Assembler::shouldExecute() && $num_crappies>0) {
	print "Finished steps that clean up crappies.\n";
    }

    if (&Assembler::shouldExecute() && $smooth_bubbles) {
	print "Start steps that smooth bubbles.\n";
	$descriptionLine = "CGB, smooth bubbles.";
        $commandLine = "$CGB_CMD $U_UNITIG_ASTAT -b 0 " .
	    "-U 1 $prefix.frgStore $prefix.fgbStore ";
	&Assembler::runLocal($descriptionLine,$commandLine);
    }
    if (&Assembler::shouldExecute() && $smooth_bubbles) {
	$descriptionLine = "FGB, after smoothing bubbles.";
#        $commandLine = "$FGB_CMD " . 
#	    "-a -i $prefix.fgbStore " .
# A work around for a bug in append mode ...            
        $commandLine = "$FGB_CMD " . 
            "$degreeThresholds " . "$overlapErrorThreshold " .
            "-c -f -o $prefix.fgbStore " .
            "-I $prefix.ovlStore " .
            "$prefix.ofg " .
	      (($num_crappies > 0)?"$prefix.breaker_overlaps.ovl":"").
            "$prefix.bubble_edges.ovl ";
	    
	&Assembler::runLocal($descriptionLine,$commandLine);
    }
    if (&Assembler::shouldExecute() && $smooth_bubbles) {
        print "Finished steps that smooth bubbles.\n";
    }
            
    if (&Assembler::shouldExecute()) {
	$descriptionLine = "CGB, final run.";
        $commandLine = "$CGB_CMD -b 0 $U_UNITIG_ASTAT " .
	    "-U 0 $prefix.frgStore $prefix.fgbStore";
	&Assembler::runLocal($descriptionLine,$commandLine);
    }

}

#----------------------------------------------------------------
# Run Consensus after Unitigger.
# Parameter 1: First number for numbering each step (for the logs).
#----------------------------------------------------------------
sub run_unitig_consensus {
    my ($nextStep) = @_;
    &Assembler::setNextStep($nextStep);
    my $commandLine;
    my $descriptionLine;
    my $prefix = &Assembler::getGlobal("prefix");

    my $CONSENSUS_CMD = "consensus $OUTPUT_MODE";
    my $CONSENSUS_IN_PARALLEL = 
	&Assembler::getGlobal("use LSF for consensus");

    if (&Assembler::shouldExecute() && !($CONSENSUS_IN_PARALLEL)) {
	$descriptionLine = "Consensus on unitigs. Non-parallel.";
	$commandLine = "$CONSENSUS_CMD -U $prefix.frgStore $prefix.cgb";
	&Assembler::runLocal($descriptionLine,$commandLine);
    }

    if (&Assembler::shouldExecute() && $CONSENSUS_IN_PARALLEL) {
	$descriptionLine = "Make partition file for ". 
	    "unitig consensus in parallel.";
        $commandLine = "make_partitionFile $prefix";
	&Assembler::runLocal($descriptionLine,$commandLine);
    }
    if (&Assembler::shouldExecute() && $CONSENSUS_IN_PARALLEL) {
	$descriptionLine = "Make partition FragStore " . 
	    "for unitig consensus in parallel.";
        $commandLine = "partitionFragStore $prefix.partFile " . 
	    "$prefix.frgStore $prefix.frgStore_part";
	&Assembler::runLocal($descriptionLine,$commandLine);
    }
    if (&Assembler::shouldExecute() && $CONSENSUS_IN_PARALLEL) {
	$descriptionLine = "Consensus on unitigs, in parallel.";
        $commandLine = "submit_consensus assembly $prefix";
	&Assembler::runLocal($descriptionLine,$commandLine);
    }
}


#----------------------------------------------------------------
# Run CGW (Scaffolder).
# Parameter 1: First number for numbering each step (for the logs).
#----------------------------------------------------------------
sub run_scaffolder {
    my ($nextStep) = @_;
    &Assembler::setNextStep($nextStep);
    my $commandLine;
    my $descriptionLine;
    my $prefix = &Assembler::getGlobal("prefix");

    my $D_UNITIG_ASTAT = "-j " . 
	&Assembler::getGlobal("D-unitig A-statistic");
    my $U_UNITIG_ASTAT = "-k " . 
	&Assembler::getGlobal("U-unitig A-statistic");

    my $CGW_CMD;
    $CGW_CMD="cgw -c $D_UNITIG_ASTAT $U_UNITIG_ASTAT " .
	"-r 4 -s 2 -w 1 -T $OUTPUT_MODE";

    if (&Assembler::shouldExecute()) {
	$descriptionLine = "CGW, the Scaffolder.";
        $commandLine = "$CGW_CMD -f $prefix.frgStore " . 
	    "-g $prefix.gkpStore " .
		"-o $prefix $prefix.cgi";
	&Assembler::runLocal($descriptionLine,$commandLine);
    }
}

#----------------------------------------------------------------
# Run Consensus after CGW.
# Parameter 1: First number for numbering each step (for the logs).
#----------------------------------------------------------------
sub run_scaffold_consensus {
    my ($nextStep) = @_;
    &Assembler::setNextStep($nextStep);
    my $commandLine;
    my $descriptionLine;
    my $prefix = &Assembler::getGlobal("prefix");

    my $CONSENSUS_CMD="consensus $OUTPUT_MODE";

    if (&Assembler::shouldExecute()) {
	$descriptionLine = "Concatenate the cgw files.";
        $commandLine="cat $prefix.cgw $prefix.cgw_contigs ". 
	    "$prefix.cgw_scaffolds >$prefix.cgw_total";
	&Assembler::runCommand($descriptionLine,$commandLine);
    }

    if (&Assembler::shouldExecute()) {
	$descriptionLine = "Consensus on scaffolds.";
        $commandLine="$CONSENSUS_CMD $prefix.frgStore $prefix.cgw_total";
	&Assembler::runLocal($descriptionLine,$commandLine);
    }

    if (&Assembler::shouldExecute()) {
	$descriptionLine = "Removing concatenated cgw files.";
        $commandLine="rm -f $prefix.cgw_total";
	&Assembler::runCommand($descriptionLine,$commandLine);
    }

}

#----------------------------------------------------------------
# Run Terminator.
# Parameter 1: First number for numbering each step (for the logs).
#----------------------------------------------------------------
sub run_terminator {
    my ($nextStep) = @_;
    &Assembler::setNextStep($nextStep);
    my $commandLine;
    my $descriptionLine;
    my $prefix = &Assembler::getGlobal("prefix");

    my $TERMINATOR_CMD;
    $TERMINATOR_CMD="terminator $OUTPUT_MODE";
    # Terminator does not support the $FILE_CLOBBER option!
    # The -f option specifies the frgStore instead.
    # Terminator clobbers the output file (prefix.asm) always.

    if (&Assembler::shouldExecute()) {
	$descriptionLine = "Terminator. Output the assembly results.";
        $commandLine="$TERMINATOR_CMD -g $prefix.gkpStore " . 
	    "-f $prefix.frgStore " .
		"-i $prefix.cns " .
		    "-o $prefix.asm " . 
			"-m $prefix.map";
	&Assembler::runLocal($descriptionLine,$commandLine);
    }
}
        
#######################################################################
# End assembler subroutines.
#######################################################################


#----------------------------------------------------------------
# Construct a local bin directory.
# Keep executables local to guard against recompilations.
# Later, make this an optional step regardless of -start & -end options.
# Parameter 1: First number for numbering each step (for the logs).
#----------------------------------------------------------------
sub make_bin_dir {
    my ($nextStep) = @_;
    &Assembler::setNextStep($nextStep);
    my $commandLine;
    my $descriptionLine;
    my ($as_root_dir) = &Assembler::getGlobal("AS_ROOT dir");
    my ($local_bin) = &Assembler::getGlobal("local bin");

    if (&Assembler::shouldExecute()) {
	$descriptionLine = "Remove bin directory.";
        $commandLine="rm -rf $local_bin";
	&Assembler::runCommand($descriptionLine,$commandLine);
    }
    if (&Assembler::shouldExecute()) { 
	$descriptionLine = "Make bin directory.";
        $commandLine="mkdir $local_bin";
	&Assembler::runCommand($descriptionLine,$commandLine);
    }
    if (&Assembler::shouldExecute()) { 
	$descriptionLine = "Populate bin directory with executables.";
        $commandLine="cp $as_root_dir/bin/* $local_bin";
	&Assembler::runCommand($descriptionLine,$commandLine);
    }
    if (&Assembler::shouldExecute()) { 
	$descriptionLine = "Populate bin directory with scripts.";
	$commandLine = "cp " .
	    "$as_root_dir/scripts/make_partitionFile " .
	    "$as_root_dir/scripts/submit_consensus " .
	    "$local_bin";
	&Assembler::runCommand($descriptionLine,$commandLine);
    }
}


    
#----------------------------------------------------------------
# Usage.
# Print the usage message if command line is incorrect.
#----------------------------------------------------------------
sub usage {
    my ($progName) = @_;
    print "Usage: $progName \n";
    print " Required arguments:\n";
    print "   -root=<string> # path to directory above bin\n";
    print "   -prefix=<string> # first name of data files\n";
    print " Optional arguments:\n";
    print "   -start=<integer> # number of first step to execute\n";
    print "   -end=<integer> # number of last step to execute\n";
    print "   -echo # Echo commands but do not execute them\n";
}



#----------------------------------------------------------------
# Main.
#----------------------------------------------------------------
sub main {    
    # Process command line ----------------------------
    my $as_root_dir = "";
    my $prefix = "";
    my $startAt = 0; # Which step to execute first (zero for all).
    my $endAt = -1;  # Which step to execute last (negative for all).
    my $echoOnly = 0;   # Echo commands but don't execute them.

    &GetOptions(
                "root=s", => \$as_root_dir, 
                "prefix=s", => \$prefix,
                "start=i", => \$startAt,
                "end=i", => \$endAt,
		"echo", => \$echoOnly);
    if(
       ($as_root_dir eq "") ||
       ($prefix eq "") ) {
        usage($PROGRAM_NAME);
        exit 1;
    }

    # Set up ----------------------------
    # This is a crude attempt to reduce global variables.
    # For now, hash global values to readable strings.
    # This allows us to log all global values.
    # Later, read from an optional properties file.
    &Assembler::setGlobal("AS_ROOT dir",$as_root_dir);
    &Assembler::setGlobal("prefix",$prefix);
    &Assembler::setGlobal("local bin", "./bin");
    &Assembler::setGlobal("repeat library","$as_root_dir/lib/$prefix.lib");
    &Assembler::setGlobal("D-unitig A-statistic", 1);
    &Assembler::setGlobal("U-unitig A-statistic", 5);
    &Assembler::setGlobal("use LSF for consensus", 0);
    &Assembler::setGlobal("use screener", 0);
    &Assembler::setGlobal("use fragment correction", 0);
    &Assembler::setGlobal("use countmessages", 0);
    
    # Handle the -start and -end options for this script.
    &Assembler::setStartAndEnd($startAt,$endAt);

    # Handle the -echo option for this script.
    &Assembler::setEchoMode($echoOnly);

    # Establish an error log. Truncate it if it already exists.
    my ($stderrRedirect) = "$prefix.stderr.log";
    &Assembler::setGlobal("stderr redirected to", $stderrRedirect);
    if (! $echoOnly && $startAt==0) {
	# So we can always append from now on.
        truncate $stderrRedirect, 0;  
    }

    # The -P option uses ('prototype') ascii text files.
    # The -B option uses binary files.
    # These should be universal to all components. Are they?
    &Assembler::setGlobal("output mode", "-P");
    $OUTPUT_MODE = &Assembler::getGlobal("output mode");

    # The -f option is used to clobber an existant store with the same name.
    # The -c option is used to create an empty store.
    # These should be universal to all components. Are they?
    &Assembler::setGlobal("clobber mode", "-f");  # should be "-c -f" ?
    $FILE_CLOBBER = &Assembler::getGlobal("clobber mode");

    # Start logging ----------------------------
    my $CurrentWorkingDirectory = cwd(); 
    print "\n";
    print "Current directory: $CurrentWorkingDirectory\n"; 
    print "Operating system: $OSNAME   Process ID: $PROCESS_ID\n";
    print "Debug flags: $DEBUGGING\n";
    print "Script: $PROGRAM_NAME\n";
    print "Interpreter: $EXECUTABLE_NAME   Version: $PERL_VERSION\n";
    print "Perl library search path: @INC\n";
    print "Time of launch: " . localtime($BASETIME) . "\n";

    print "\nGlobal settings:\n";
    &Assembler::printGlobals();
    print "\n";

    # Run the assembler  ---------------------------
    print "########## START OF LOG ##########\n";
    eval {
	# Parameters are arbitrary widely spaced increasing numbers.
	# They provide capability to restart script at arbitrary point.
	# Incremental numbers are printed to the logs on stdout and stderr.

	make_bin_dir(100);
	run_gatekeeper(200);
	run_screener(300);
	run_overlapper(400);
	run_unitigger(500);
	run_unitig_consensus(600);
	run_scaffolder(700);
	run_scaffold_consensus(800);
	run_terminator(900);
    }; # semicolon required! 

    if ($@) {
	# Error handling for eval() block above.
	print "\n";
	warn $@; 
	print "\nTime of termination: " . localtime() . "\n";
	print "-------- PREMATURE END OF LOG ----------\n";
	exit 1;
    }

    # Clean up --------------------------------
        # Don't yet clean the directory of working files 
	# *.inp *.ovl *.cgb *cgw;
    print "\nTime of completion: " . localtime() . "\n";
    print "########## SUCCESSFUL END OF LOG ##########\n";
}
# ------- Run this program.
main();
exit 0;

#######################################################################
