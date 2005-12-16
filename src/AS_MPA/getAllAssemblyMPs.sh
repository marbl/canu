#!/usr/bin/env bash
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
# $Id: getAllAssemblyMPs.sh,v 1.5 2005-12-16 22:13:07 catmandew Exp $
#

if [ -z ${DATA_DIR} ]; then
  if [ ${OS} == "AIX" ] || [ ${OS} == "OSF1" ]; then
    export DATA_DIR=/prod/IR01/dewim/mps/human
  else
    export DATA_DIR=/home/dewim/celera/sandbox/cds/IR/COMPASS/data/human
  fi
fi

ASTORE=${DATA_DIR}/stores/human_vanilla.asmStore

ASSEMBLY=${1}

if [ -z ${ASSEMBLY} ] || [ ${ASSEMBLY} == "bell-style" ]; then
  echo "Please specify an assembly name"
  return;
fi
MSTORE=${DATA_DIR}/stores/${ASSEMBLY}.mapStore

cd ${DATA_DIR}/${ASSEMBLY}/intraChromosome
chmod -R +w *
cd unmapped
dumpMappedMatePairs -s ${ASTORE} -m ${MSTORE}

for file in `ls *.mmps.txt`; do
  chr=${file%%.*}
  gawk '{print $1, $3, $4, $5, $6, $7}' ${file} > ${chr}.txt
done

rm *.mmps.txt

mv 00?.txt ..
mv 01?.txt ..
mv 02[0-3].txt ..

cd ${DATA_DIR}/${ASSEMBLY}/interChromosome
chmod -R +w *
cd unmapped
dumpMappedElsewheres -s ${ASTORE} -m ${MSTORE}

for file in `ls *.elsewhere.txt`; do
  chr=${file%%.*}
  mv ${file} ${chr}.txt
done

mv 00?.txt ..
mv 01?.txt ..
mv 02[0-3].txt ..
