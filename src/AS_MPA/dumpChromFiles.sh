#!/usr/local/bin/bash
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
# $Id: dumpChromFiles.sh,v 1.2 2004-09-23 20:25:24 mcschatz Exp $
#

function MoveUnmapped
{
  mv [1-9][0-9][0-9].txt unmapped
  mv 0[3-9][0-9].txt unmapped
  mv 02[4-9].txt unmapped
}

AS=${1}

if [ -z ${DATA_DIR} ]; then
  if [ ${OS} == "AIX" ] || [ ${OS} == "OSF1" ]; then
    export DATA_DIR=/prod/IR01/dewim/mps/human
  else
    export DATA_DIR=/home/dewim/celera/sandbox/cds/IR/COMPASS/data/human
  fi
fi

if [ -z ${AS} ] || [ ${AS} == "bell-style" ] ; then
  echo "Please identify the assembly (B33A, VAN, HG06, ...)"
  return
fi

if [ ! -d ${DATA_DIR}/${AS} ]; then
  mkdir ${DATA_DIR}/${AS}
fi

if [ ! -d ${DATA_DIR}/${AS}/intraChromosome ]; then
  mkdir ${DATA_DIR}/${AS}/intraChromosome
fi

if [ ! -d ${DATA_DIR}/${AS}/intraChromosome/unmapped ]; then
  mkdir ${DATA_DIR}/${AS}/intraChromosome/unmapped
fi

if [ ! -d ${DATA_DIR}/${AS}/interChromosome ]; then
  mkdir ${DATA_DIR}/${AS}/interChromosome
fi

if [ ! -d ${DATA_DIR}/${AS}/interChromosome/unmapped ]; then
  mkdir ${DATA_DIR}/${AS}/interChromosome/unmapped
fi

echo "Dumping intraChromosome mate pairs"
cd ${DATA_DIR}/${AS}/intraChromosome
dumpMappedMatePairs -s ${DATA_DIR}/stores/human_vanilla.asmStore -m ${DATA_DIR}/stores/${AS}.mapStore

for file in `ls *.mmps.txt`; do
  chr=${file%%.*}
  echo "Converting ${chr} records"
  gawk '{print $1, $3, $4, $5, $6, $7}' ${file} > ${chr}.txt
  rm ${file}
  chmod -w ${chr}.txt
done
MoveUnmapped

echo "Dumping interChromosome mate pairs"
cd ${DATA_DIR}/${AS}/interChromosome
dumpMappedElsewheres -s ${DATA_DIR}/stores/human_vanilla.asmStore -m ${DATA_DIR}/stores/${AS}.mapStore

for file in `ls *.elsewhere.txt`; do
  chr=${file%%.*}
  mv ${file} ${chr}.txt
  chmod -w ${chr}.txt
done
MoveUnmapped

