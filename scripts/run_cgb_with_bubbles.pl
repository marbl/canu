#!/bin/perl
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

#######################################################################
#
# run_cgb_with_bubbles.pl
#
# This script runs the cgb phase of Assembler Grande with bubbles.
# It assumes that fgb has already been run to create a graph store.
# The script runs cgb twice with an intervening run of fgb.
#
# Usage: run_cgb_with_bubbles.pl -C "cgb opts" -F "fgb opts" <frag store> \
#           <graph_store>
#
# Options:
#
#   -C <string>      Command line options for cgb.  Should not include
#                    the stores, which are appended automatically by the
#                    script.  The -U and -P options are also handled 
#                    automatically by the script and are not needed here.
#   
#   -F <string>      Command line options for fgb.  As with cgb, the store
#                    names and -P option are handled by the script, so
#                    this option is usually unnecessary.
#
#   -P <0 or 1>      (default = 1)  Boolean flag to enable ASCII output.
#
#   -T               Run in test mode.  Echos actions that would be
#                    taken, but does not execute them.
#
#   -r <1, 2>        Restart the script after a failed run.  1 = restart
#                    with FGB, 2 = restart with the final CGB stage.
#
#   Example:
#
#   $ ls
#   genome.fgbStore/     genome.frgStore/
#   $ run_cgb_with_bubbles.pl -C "-S -A 1 -b 4" genome.frgStore genome.fgbStore
#
#######################################################################

require "getopts.pl";

$options = "hTr:C:F:P:";		# Option vector

# $asbin = ".";  # path to executable
$asbin = $ENV{AS_BIN};
print "AS_BIN = $asbin\n";

if ($opt_h || ($#ARGV < 1)) {
    print("Usage: run_cgb_with_bubbles.pl -C 'cgb opts' -F 'fgb opts' <frag store> \ \n");
    print("           <graph_store>\n\n");
    
    print(" Options:\n\n");
    
    print("   -C <string>      Command line options for cgb.  Should not include\n");
    print("                    the stores, which are appended automatically by the\n");
    print("                    script.  The -U and -P options are also handled \n");
    print("                    automatically by the script and are not needed here.\n\n");
    
    print("   -F <string>      Command line options for fgb.  As with cgb, the store\n");
    print("                    names and -P option are handled by the script, so\n");
    print("                    this option is usually unnecessary.\n\n");
    
    print("   -P <0 or 1>      (default = 1)  Boolean flag to enable ASCII output.\n\n");
    
    print("   -T               Run in test mode.  Echos actions that would be\n");
    print("                    taken, but does not execute them.\n\n");
    print("   -r <1, 2>        Restart the script after a failed run.  1 = restart\n");
    print("                    with FGB, 2 = restart with the final CGB stage.\n");
    exit(0);
}

#
# Process arguments and set variables
#

&Getopts($options);

$frg_store = $ARGV[$#ARGV - 1];
$fgb_store = $ARGV[$#ARGV];

$do_phase_1 = (($opt_r ne "1") && ($opt_r ne "2"));
$do_phase_2 = ($opt_r ne "2");

# Set up flags
$cgb_flags_pre = $opt_C . " -U 1"; # Options for bubble-edge generating run
$fgb_flags = $opt_F . " -f";
$cgb_flags_post = $opt_C . " -U 0"; # Options for post bubble-edge run

if ($opt_P ne "0") {
    $cgb_flags_pre = $cgb_flags_pre . " -P";
    $fgb_flags = $fgb_flags . " -P";
    $cgb_flags_post = $cgb_flags_post . " -P";
}

# Find file prefix
if ($fgb_store =~ /([\S]+)\.fgbStore/) {
    $prefix = $1;
}
elsif (($fgb_store =~ /([\S]+)\.Graph/)) {
    $prefix = $1;
}
else {
    die("Could not find prefix from $fgb_store.  Badger Dan for better test.");
}

# Set up file names

$bubble_ovls = $prefix . ".pre.bubble_edges.ovl";
$pre_frg_store = $frg_store;
$pre_fgb_store = $prefix . ".pre.fgbStore";
$post_fgb_store = $prefix . ".post.fgbStore";
$post_frg_store = $frg_store;
$log_file = "run_cgb_with_bubbles.log";

# Set up command strings

$cgb_pre = "$asbin/cgb $cgb_flags_pre $pre_frg_store $pre_fgb_store";
$fgb = "$asbin/fgb $fgb_flags -i $pre_fgb_store -o $post_fgb_store $bubble_ovls";
$cgb_post = "$asbin/cgb $cgb_flags_post $post_frg_store $post_fgb_store"; 

#
# Prepare the directory and test for file existence
#

print "====================================================================\n";
print "Preparing to run ...\n\n";

if ($opt_T) {
    print "# ln -sf $fgb_store $pre_fgb_store\n";
}
else {
    print "Log file will be $log_file \n";
    print "Creating link to fragment store ... \n";
    $cmd_line="ln -s -f $fgb_store  $pre_fgb_store";
    print "Executing: $cmd_line \n";
    (system($cmd_line) == 0) || die("FAILED");
    print "Done\n";
}

if ($opt_T) {
    print "# Open $log_file for output.\n";
}
else {
    print "Opening log file ... ";
    if (-f $log_file) {
	open(LOG_FILE, ">> $log_file") || die("FAILED");
	print LOG_FILE "\n\nWARNING: Script Restart =============================================\n\n";
	print "Done\n";
    }
    else {
	open(LOG_FILE, "> $log_file") || die("FAILED");
	print "Done\n";
    }
}

#
# Run CGB to generate bubble edges
#

if ($do_phase_1) {
    print "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n";
    print "Bubble edge generation CGB phase ...\n\n";
    
    if (!$opt_T) {
	print LOG_FILE "\n!!!!!\n!!!!! BUBBLE EDGE GENERATION PHASE\n!!!!!\n";
	print LOG_FILE "!!!!! Cmd line: $cgb_pre\n!!!!!\n";
	close LOG_FILE;
        $cmd_line="$cgb_pre  >>  $log_file  2>&1";
	print "Executing: $cmd_line \n";
	(system($cmd_line) == 0) || die("FAILED.  Aborting.");
	print "Done.\n";
    }
    else {
	print "# ${cmd_line}\n\n";
    }
    
}
else {
    print "Skipping bubble edge generation phase.\n";
}

#
# Run FGB to incorporate bubble edges
#

if ($do_phase_2) {
    print "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n";
    print "Running FGB to build store with bubble edges ...\n\n";
    
    if (!$opt_T) {
	open(LOG_FILE, ">> $log_file") || die("Weird.  Couldn't re-open log.");
	print LOG_FILE "\n!!!!!\n!!!!! BUBBLE EDGE INCORPORATION\n!!!!!\n";
	print LOG_FILE "!!!!! Cmd line: $fgb\n!!!!!\n";
	close LOG_FILE;
        $cmd_line="$fgb >> $log_file 2>&1";
        print "Executing: $cmd_line \n";
	(system($cmd_line) == 0) || die("FAILED.  Aborting.");
	print "Done.\n";
    }
    else {
	print "# $fgb >> $log_file 2>&1\n\n";
    }
}
else {
    print "Skipping graph store generation phase.\n";
}

#
# Run CGB to generate final .cgb file with bubbles popped
#

print "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n";
print "Final CGB phase ...\n\n";

if (!$opt_T) {
    open(LOG_FILE, ">> $log_file") || die("Weird.  Couldn't re-open log.");
    print LOG_FILE "\n!!!!!\n!!!!! BUBBLE EDGE GENERATION PHASE\n!!!!!\n";
    print LOG_FILE "!!!!! Cmd line: $cgb_post\n!!!!!\n";
    close LOG_FILE;
    $cmd_line="$cgb_post >> $log_file 2>&1";
    print "Executing: $cmd_line \n";
    (system($cmd_line) == 0) || die("FAILED.  Aborting.");
    print "Done.\n\n";

    print "Final CGB file is $prefix.post.cgb\n";
}
else {
    print "# $cgb_post >> $log_file 2>&1\n\n";
    print "Final CGB file would be $prefix.post.cgb\n";
}



