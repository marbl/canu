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
# $Id: fixOffsets.sh,v 1.2 2004-09-23 20:25:24 mcschatz Exp $
#

if [ -z ${DATA_DIR} ]; then
  if [ ${OS} == "AIX" ] || [ ${OS} == "OSF1" ]; then
    export DATA_DIR=/prod/IR01/dewim/mps/human
  else
    export DATA_DIR=/home/dewim/celera/sandbox/cds/IR/COMPASS/data/human
  fi
fi

function ProcessIntraFiles
{
  for file in `ls [0-9][0-9][0-9].txt`; do
    chr=${file%%.*}
    gawk '{print $1, $2, $3, $4, $5-1, $6-1}' ${chr}.txt > ${chr}new.txt
    #chmod +w ${file}
    #mv ${chr}new.txt ${file}
    #chmod -w ${file}
  done
}

function ProcessInterFiles
{
  for file in `ls [0-9][0-9][0-9].txt`; do
    chr=${file%%.*}
    gawk '{print $1, $2, $3-1, $4, $5, $6, $7-1, $8, $9}' ${chr}.txt > ${chr}new.txt
    #chmod +w ${file}
    #mv ${chr}new.txt ${file}
    #chmod -w ${file}
  done
}

Assemblies=(B28 B33A B34 CSAB HG06 VAN WGAB)

for AS in "${Assemblies[@]}"; do
  echo "Fixing offsets in ${AS}"

  echo "  working on intraChromosome files"
  cd ${DATA_DIR}/${AS}/intraChromosome
  ProcessIntraFiles
  cd unmapped
  ProcessIntraFiles

  echo "  working on interChromosome files"
  cd ${DATA_DIR}/${AS}/interChromosome
  ProcessInterFiles
  cd unmapped
  ProcessInterFiles
  
done