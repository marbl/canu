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
# $Id: explainGaps.sh,v 1.4 2005-03-22 19:48:58 jason_miller Exp $
#

AS=VAN

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

# get intra-scaffold gaps not spanned by 2k/10ks
for file in `ls [0-9][0-9][0-9].txt`; do
  chr=${chr%%.*}
  echo ${chr}
  sed 's/[,:]/ /g' ${AS}.${chr}.2k10kGapClones.txt | gawk '{if(NF==2)print $0}' > ${AS}.${chr}.gapsWO2k10kClones.txt
done

# get intervals spanned by 50ks & bacs
for file in `ls [0-9][0-9][0-9].txt`; do
  chr=${file%%.*}
  fgrep -f ${DATA_DIR}/libs/50kLibs.txt ${file} | gawk '{print $5, $6-$5}' > ${AS}.${chr}.50kCloneIntervals.txt
  fgrep -f ${DATA_DIR}/libs/bacEndLibs.txt ${file} | gawk '{print $5, $6-$5}' > ${AS}.${chr}.bacIntervals.txt
done

# intersect left over gaps w 50ks
for file in `ls [0-9][0-9][0-9].txt`; do
  
done

