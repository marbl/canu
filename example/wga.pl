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

###########################################################

use Getopt::Long;   # command line options longer than 1 letter
my $commandLine;

#This are the assembler script  arguments

$prefix="a006";
$AS_BIN="../Linux64/bin";

#The update  variable allows cgw to to run with different values for j and K
# -i <thresh>]  Set max coverage stat for microhet determination of non-uniqueness (default -1)         
# -j <thresh>]  Set min coverage stat for definite uniqueness         
# -k <thresh>]  Set max coverage stat for possible uniqueness  


$UTG_Options_with_Frag_corr = "-c -P -A 0  -n 700 -m 3500  -d 1 -x 1 -z 5 -j 5 -U 1  "; 

############################## CGW options  ################################

$CGWOptions = "-c -j 1 -k 5 -r 4 -s 2 -w 0 -T -P";

############################################################################

#$ENV{'PATH'} .= ":$AS_BIN";

############################################################################


############################################################################
#                      Gatekeeper                                           
#                                                                           
############################################################################


print "Execute: gatekeeper\n";

if (! -e "${prefix}.inp") {
#print "Time of launch the gatekeeper : " . localtime() . "\n";

$commandLine = "$AS_BIN/gatekeeper -X -C -N -Q -P -f ${prefix}.gkpStore ${prefix}.frg";

print "Execute: $commandLine\n";
$systemReturn = system($commandLine) && die "Failed: $commandLine";
print "systemReturn = $systemReturn\n";
print "Finishing the gkp run : " . localtime() . "\n";
}

############################################################################
#                      Populator                                            
#                                                                           
############################################################################
print "Execute: populator\n";


if (! -e "done.populator") {
$commandLine = "$AS_BIN/PopulateFragStore -f -c -o ${prefix}.frgStore ${prefix}.inp";
print "Execute: $commandLine\n";
$systemReturn = system($commandLine) && die "Failed: $commandLine";
print "systemReturn = $systemReturn\n";
print "Finishing the populator: " . localtime() . "\n";
system("touch done.populator");
}


print "Time of launch the ofglist: " . localtime() . "\n";

$commandLine = "ls -1 *.ovl > ${prefix}.ofglist";
print "Execute: $commandLine\n";
$systemReturn = system($commandLine) && die "Failed: $commandLine";
print "systemReturn = $systemReturn\n";
print "Time finished the ofglist: " . localtime() . "\n";


############################################################################
#                      Make Range File                                      
#                                                                           
############################################################################

print "Making the range file: \n";

if (! -e "${prefix}.range") {
$commandLine = ("time ${AS_BIN}/make_range_file ${prefix}.ovl  ${prefix}.range");
print "Execute: $commandLine\n";
$systemReturn = system($commandLine) && die "Failed: $commandLine";
print "systemReturn = $systemReturn\n";
print "Finishing makeng the range file" . localtime() . "\n";
}


####################make the ovl script#####################################
# The overlapper need different steps: Meryl, create the ovelap script, create the frag correction script 
# and create the olap store.  LSF NOT  active if the -Q option is present                                 
############################################################################

############################################################################
#                      Meryl                                                
#                                                                           
############################################################################

print "Execute: Meryl\n";

if (! -e "${prefix}.counts.fasta") {
$commandLine = ("time $AS_BIN/meryl -s ${prefix}.frgStore -m 22 -n 100  -o ${prefix}.counts.fasta");
print "Execute: $commandLine\n";
$systemReturn = system($commandLine) && die "Failed: $commandLine";
print "systemReturn = $systemReturn\n";
print "Finishing makeng the range file" . localtime() . "\n";
}


############################################################################
#           Creation and execution of the overlap                           
#                                                                           
############################################################################

if (! -e "done.overlap") {
print "Create and Execute: Overlap\n";

$commandLine = "time $AS_BIN/overlap -M 1GB -h 1-700 -r  1-700  -w -k ${prefix}.counts.fasta -o  batch.r1-700h1-700.ovl ${prefix}.frgStore ";

print "Execute: $commandLine\n";
$systemReturn = system($commandLine) && die "Failed: $commandLine";
print "systemReturn = $systemReturn\n";

system("touch done.overlap");
}


############################################################################
#                Moving the ofg files to OVL directory                      
#                                                                           
############################################################################


print "Create the Overlap Store\n";

print "Time of launch the ovl list: " . localtime() . "\n";

$commandLine = "ls -1 *.ovl > OVL_LIST";
print "Execute: $commandLine\n";
$systemReturn = system($commandLine) && die "Failed: $commandLine";
print "systemReturn = $systemReturn\n";
print "Time finished the ovl list: " . localtime() . "\n";


############################################################################
#                                               created the ovl Store       
#                                                                           
############################################################################
print "Time of lunch the grow-olap-store: " . localtime() . "\n";

$commandLine = "time $AS_BIN/grow-olap-store -cfS -o ${prefix}.ovlStore *.ovl";
print "Execute: $commandLine\n";
$systemReturn = system($commandLine) && die "Failed: $commandLine";
print "systemReturn = $systemReturn\n";
print "Time finished the grow-olap-store: " . localtime() . "\n";


############################################################################
#                                          Fragment correction              
#                                                                           
############################################################################

print "Time of lunch correct frags " . localtime() . "\n";
$commandLine = "time $AS_BIN/correct-frags -k 9 -S ${prefix}.ovlStore -x 1 -o ${prefix}.corr ${prefix}.frgStore 1 700";
print "Execute: $commandLine\n";
$systemReturn = system($commandLine) && die "Failed: $commandLine";
print "systemReturn = $systemReturn\n";
print "Time finished correct frags : " . localtime() . "\n";

print "Time of lunch correct olaps: " . localtime() . "\n";
$commandLine = "time $AS_BIN/correct-olaps -S ${prefix}.ovlStore -e corrolap1695912.erate ${prefix}.frgStore ${prefix}.corr 1 700";
print "Execute: $commandLine\n";
$systemReturn = system($commandLine) && die "Failed: $commandLine";
print "systemReturn = $systemReturn\n";
print "Time finished correct olaps : " . localtime() . "\n";
 

print "Time of lunch the update erates: " . localtime() . "\n";
$commandLine = "time $AS_BIN/update-erates ${prefix}.ovlStore corrolap1695912.erate";
print "Execute: $commandLine\n";
$systemReturn = system($commandLine) && die "Failed: $commandLine";
print "systemReturn = $systemReturn\n";
print "Time to finished update-erates: " . localtime() . "\n";



############################################################################
#                         unitigger                                         
#                                                                           
############################################################################

print "Execute: Unitigger phase\n";

if (! -e "${prefix}.cgb"){
print "Time of launching the UNITIGGER : " . localtime() . "\n";
$commandLine = "time $AS_BIN/unitigger -c -P -A 0  -n 700 -m 3000  -d 1 -x 1 -z 5 -j 5 -U 1 -F ${prefix}.frgStore -L ${prefix}.ofglist -f -o ${prefix}.fgbStore -I ${prefix}.ovlStore 1> ${prefix}.utg.stdout  2> ${prefix}.utg.stderr";
print "Execute: $commandLine\n";
$systemReturn = system($commandLine) && die "Failed: $commandLine";
print "systemReturn = $systemReturn\n";
print "Finishing the unitigger: " . localtime() . "\n";
}


############################################################################
#              post-unitigger consensus
#                                                                           
############################################################################
print "consensus phase\n";

if (! -e "${prefix}.cgi"){
$commandLine = "time $AS_BIN/consensus -P -U ${prefix}.frgStore ${prefix}.cgb 2>>${prefix}.stderr.log";
print "Execute: $commandLine\n";
$systemReturn = system($commandLine) && die "Failed: $commandLine";
print "systemReturn = $systemReturn\n";
}


############################################################################
#              cgw                                                          
#                                                                           
############################################################################

print "Execute: cgw_phase\n";

if (! -e "${prefix}.cgw"){
$commandLine = "time $AS_BIN/cgw -c -j 1 -k 5 -r 4 -s 2 -w 0 -T -P -f ${prefix}.frgStore -g ${prefix}.gkpStore -o ${prefix} *.cgi 2> ${prefix}.stderr.log";
print "Execute: $commandLine\n";
$systemReturn = system($commandLine) && die "Failed: $commandLine";
print "systemReturn = $systemReturn\n";
}

if (! -e "${prefix}.cgw_total"){
$commandLine = "time cat ${prefix}.cgw ${prefix}.cgw_contigs ${prefix}.cgw_scaffolds > ${prefix}.cgw_total";
print "Execute: $commandLine\n";
$systemReturn = system($commandLine) && die "Failed: $commandLine";
print "systemReturn = $systemReturn\n";
}


############################################################################
#              Consensus                                                    
#                                                                           
############################################################################

if (! -e "${prefix}.cns"){
$commandLine = "time $AS_BIN/consensus -P ${prefix}.frgStore ${prefix}.cgw_total 2>>${prefix}.stderr.log";
print "Execute: $commandLine\n";
$systemReturn = system($commandLine) && die "Failed: $commandLine";
print "systemReturn = $systemReturn\n";

}


############################################################################
#              Terminator                                                   
#                                                                           
############################################################################

$commandLine = "time $AS_BIN/terminator -P -g ${prefix}.gkpStore -f ${prefix}.frgStore -i ${prefix}.cns -o ${prefix}.asm -m ${prefix}.map 2>>${prefix}.stderr.log";
print "Execute: $commandLine\n";
$systemReturn = system($commandLine) && die "Failed: $commandLine";
print "systemReturn = $systemReturn\n";


print "The assembler end  at  : " . localtime() . "\n"
