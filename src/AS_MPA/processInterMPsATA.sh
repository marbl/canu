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
# $Id: processInterMPsATA.sh,v 1.3 2005-03-22 19:05:58 jason_miller Exp $
#

AS=${1}
mapping=${2}

if [ -z ${DATA_DIR} ]; then
  if [ ${OS} == "AIX" ] || [ ${OS} == "OSF1" ]; then
    export DATA_DIR=/prod/IR01/dewim/mps/human
  else
    export DATA_DIR=/home/dewim/celera/sandbox/cds/IR/COMPASS/data/human
  fi
fi

if [ -z ${AS} ] || [ ${AS} == "bell-style" ]; then
  echo "Please identify the assembly (B33A, VAN, HG06, ...)"
  return
fi

#if [ -z ${mapping} ] || [ ${AS} == "bell-style" ]; then
#  echo "Please provide second parameter - mapped or unmapped"
#  return
#fi

for file in `ls [0-9][0-9][0-9].txt | sort -n`; do
  chrom=${file%%.*}
  echo "    working on ${chrom}"
  processInterCG -l ${DATA_DIR}/libs/humanLibs.txt -e ${file} -a ${AS} -c ${chrom} -L 2> ${chrom}.err
done

echo "  summarizing inter-chromosome results"
filterIntersWithSatisfieds.sh ${AS} ${mapping}