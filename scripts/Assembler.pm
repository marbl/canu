
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
package Assembler;

use English;   # full length variable names for perl globals
#use strict;   # force variable declartion
use Carp;   # carp() and confess() error reporting
use diagnostics;    # print stack trace on error
use sigtrap;   # signal handler prints stack trace too

my $TRUE = 1;
my $FALSE = 0;

sub import () {
    # This gets called by the perl directive "use Assembler;".
    print "[Importing perl module Assembler.pm]\n";
}

#######################################################################
# 
# Utility functions for Celera's WGS assembly scripts.
#
# Include this package in your scripts with this line:
# use Assembler;
#
# Main blame: Jason Miller 
#######################################################################



#######################################################################
# Global Variables Utility.
#
# If you use this mechanism instead of global variables,
# then your code will be cleaner,
# and you can log all global settings with one function call.
#######################################################################
my (%globalVars); 
#----------------------------------------------------------------
# Set a global.
# See getGlobal(), printGlobals().
# Implemented using the global hash %globalVars.
#----------------------------------------------------------------
sub setGlobal {
    my ($key, $val) = @_;
    $globalVars{$key} = $val;
}
#----------------------------------------------------------------
# Return the value of any global setting.
# See setGlobal().
#----------------------------------------------------------------
sub getGlobal {
    my ($key) = @_;
    my ($val);
    $val = $globalVars{$key};
    die ("No global variable defined for '$key'.") unless defined $val;
    return $val;
}
#----------------------------------------------------------------
# Print the values of all the global settings.
# See setGlobal().
#----------------------------------------------------------------
sub printGlobals {
    my (@keys, @values, $key, $val);
    @keys = sort keys %globalVars;
    foreach $key (@keys) {
	$val = $globalVars{$key};
	print " $key = $val\n";
    }
}


#######################################################################
# System Commands Utility.
#
# If you use these mechanisms instead of system() or `cmd`,
# then you will get automatic logging and full error reporting.
#######################################################################
my $echoOnly;
#----------------------------------------------------------------
# Set Echo Mode.
# Tell module to echo commands but not execute them.
# Parameter 1: 1=> echo only, else echo and execute
#----------------------------------------------------------------
sub setEchoMode {
    my ($param) = @_;
    if ($param==1) {
	print "Echo mode is on! Commands will echo but not execute.\n";
	$echoOnly = 1;
    } else {
	print "Echo mode is off! Commands will echo and execute.\n";
	$echoOnly = 0;
    }
}
#----------------------------------------------------------------
# Execute a unix command.
# Log the activity fully.
# Parameter 1: Some description of this step (for the log).
# Parameter 2: The command to execute.
#----------------------------------------------------------------
sub runCommand {
    my ($description, $commandLine) = @_;
    my ($pattern) = "#-#-#";
    my ($stderrRedirect) = &getGlobal("stderr redirected to");
    my ($errorLoggingCmd) = "2>>$stderrRedirect";
    my ($fullCommand) = "$commandLine $errorLoggingCmd";
    my ($sysReturn) = 0;
    # Improve helpfulness of stderr by marking milestones.
    if (! $echoOnly) {
        ` echo "\n$fullCommand\n" >>$stderrRedirect `;
    }
    # Improve helpfulness of stdout by marking milestones.
    print "$pattern START OF COMMAND $pattern\n";
    print "$pattern $description\n";
    print "Time is " . localtime() . "\n";
    print "Executing> $fullCommand\n";
    &timerUtility($FALSE);
    if (! $echoOnly) {
        $sysReturn = 0xffff & system ($fullCommand);
    }
    &timerUtility($TRUE);
    if ($sysReturn != 0) {
	# Print detailed error message on one grep-able line.
	# Adapted from example shown in Programming Perl, 2nd edition.
	# Adapted from example shown under 'system' command
	# in Programming Perl, 2nd edition.
	print "$pattern Failure! ";
	printf ("Exit status was %#04x (hex). ",$sysReturn);
	if ($sysReturn == 0xff00) {
	    print "Exist status indicates OS-level failure. ";
	}
	elsif ($sysReturn > 0x80) {
	    $sysReturn >>= 8;
	    print "Exit status indicates problem was not due to a signal. ";
	    print "Non-signal exit value would be $sysReturn. ";
	}
	else {
	    if ($sysReturn & 0x80) {
		$sysReturn &= ~0x80;
		print "Exit status indicates coredump. ";
	    }
	    print "Exist status indicates signal $sysReturn. ";
	}
	print "Perl's child error value is $CHILD_ERROR. ";
	print "Perl's OS error value is $OS_ERROR. ";
	print "\n";
	die;
    }
    print "$pattern END OF COMMAND $pattern \n\n";
}
#----------------------------------------------------------------
# Quick fix for executing backtick commands.
# Later, modify runCommand() to do this too.
# Parameters: same as runCommand.
#----------------------------------------------------------------
sub readChildProcess {
    my ($desc, $cmd) = @_;
    my ($returnValue) = 0;
    print "Reading from child process.\n";
    print "Executing> $cmd\n";
    if (! $echoOnly) {
	my $answer = eval `$cmd` ;
        if (! defined $answer) {
	    print "Eval error = $EVAL_ERROR. Child error = $CHILD_ERROR\n";
	    die $CHILD_ERROR;
	}
	$returnValue = $answer;
    }
    return $returnValue;
}
#----------------------------------------------------------------
# Run a unix command from the local bin.
# Same as runCommand(), but prepends the local path.
# Depends on global setting of local bin.
# Parameters: see runCommand.
#----------------------------------------------------------------
sub runLocal {
    my ($desc, $cmd) = @_;
    my ($local_bin) = &getGlobal("local bin");
    &runCommand ($desc, "$local_bin/$cmd");
}





#######################################################################
# Interrupted Runs Utility.
#
# If you invoke this test before each sequential operation,
# then your script can be stopped and restared anywhere,
# for instance based on command-line options.
# For an example, see assembler.pl.
#######################################################################
my ($nextStep,$startAt,$endAt) ;
#----------------------------------------------------------------
# Set the range of steps to be executed.
# For example, (501,502) would execute only step #501.
# Parameter 1: Inclusive. Use 0 to start at the beginning.
# Parameter 2: Exclusive. Use -1 to run to the end.
#----------------------------------------------------------------
sub setStartAndEnd {
    my ($param1,$param2) = @_;
    $startAt = $param1;
    $endAt = $param2;
    print "Process will start at step $startAt and end at step $endAt. ".
	"(0 means beginning, -1 means end.)\n";
}
#----------------------------------------------------------------
# Set the number of the next step to be executed.
# See setStartAndEnd.
# Parameter 1: The next step to be executed.
#----------------------------------------------------------------
sub setNextStep {
    my ($param1) = @_;
    $nextStep = $param1;
}
#----------------------------------------------------------------
# Determine whether to execute this step
# based on globals "step", "start at" and "end at".
# Also, increment the step counter.
# See setStartAndEnd.
#----------------------------------------------------------------
sub shouldExecute {
    my ($doIt);
    $doIt = 
	($nextStep >= $startAt) &&
	($endAt<0 || $nextStep < $endAt) ;
    if ($doIt) {
	print "For restarting this script, next step is #$nextStep...\n";
    }
    $nextStep++;
    return $doIt;
}




#######################################################################
# Timer Utility.
#
#######################################################################
#----------------------------------------------------------------
# Measure elapsed time.
# This is crude: only one timer can run, and global variables are used.
# Later, write a timer object.
# Parameter 1: FALSE=> start the clock, TRUE=> print elapsed time.
#----------------------------------------------------------------
my ($timerStart, $timerEnd);  
sub timerUtility {
    my ($outputTime) = @_;
    my ($elapsedTime);
    if (! $outputTime) {
	$timerStart = time;
    } else {
	$timerEnd = time;
	$elapsedTime = $timerEnd - $timerStart;
	print "Elapsed time $elapsedTime sec.\n", 
    }
}



1;
