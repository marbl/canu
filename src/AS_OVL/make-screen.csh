#!/usr/bin/csh
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

# Makes a multifasta file of sequences that have high
# overlap degree

# USAGE:  make-screen.csh  <tag>  <fragstore>  <lofrag>  <hifrag>  <olapdegree>

# Use  <tag>  as prefix of filex produced.  Compute overlaps of
# fragments  <lofrag> .. <hifrag>  from  <fragstore>  (numbers are IIDs).
# Create screen library from fragment sections that have >= <olapdegree>
# overlaps.
# Files created and kept are:
#   <tag>.log        Log of standard output and standard error
#   <tag>.ovlStore   Overlap store of overlaps used
#   <tag>.screen     Final screen library
#   <tag>.degr       Overlap degrees of ends of fragments

if  ($#argv < 5)  then
  echo "USAGE:  make-screen.csh  <tag>  <fragstore>  <lofrag>  <hifrag>  <olapdegree>"
  exit -1
endif

echo "Starting"

echo 'Using binaries in $AS_BIN =' $AS_BIN |& tee ${1}.log

echo "Compute overlaps" |& tee -a ${1}.log
$AS_BIN/overlap -s -w -h${3}-${4} -r${3}-${4} -o ${1}.$$.ovl ${2} |& tee -a ${1}.log
if  ($status != 0)  then
  echo "overlap command failed"
  exit -1
endif

echo "Create overlap store" |& tee -a ${1}.log
$AS_BIN/grow-olap-store -cfS -o ${1}.ovlStore ${1}.$$.ovl |& tee -a ${1}.log
if  ($status != 0)  then
  echo "grow-olap-store command failed"
  exit -1
endif

echo "Compute overlap degree on ends of each fragment" |& tee -a ${1}.log
$AS_BIN/get-degrees -S ${1}.ovlStore > ${1}.degr
if  ($status != 0)  then
  echo "get-degrees command failed"
  exit -1
endif

echo "Delete files ${1}.$$.ovl" |& tee -a ${1}.log
rm ${1}.$$.ovl
if  ($status != 0)  then
  echo "rm command failed"
  exit -1
endif

echo "Create raw screen sequences" |& tee -a ${1}.log
$AS_BIN/auto-screen -S -d $5 ${1}.degr ${2} ${1}.ovlStore
mv screen.fasta ${1}.$$.raw.screen
if  ($status != 0)  then
  echo "mv command failed"
  exit -1
endif

echo "Remove duplicate sequences" |& tee -a ${1}.log
$AS_BIN/removedups < ${1}.$$.raw.screen | $AS_BIN/remove-dup-screen -L3 > ${1}.screen
if  ($status != 0)  then
  echo "removedump or remove-dup-screen command failed"
  exit -1
endif

echo "Delete file ${1}.$$.raw.screen" |& tee -a ${1}.log
rm ${1}.$$.raw.screen
if  ($status != 0)  then
  echo "rm command failed"
  exit -1
endif
