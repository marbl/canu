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
# $Id: buildMapStoreFromMatches.sh,v 1.2 2004-09-23 20:25:23 mcschatz Exp $
#

function DeleteFile
{
 if [ -f $1 ]; then
   echo "WARNING: deleting file ${1}"
   rm -f $1
 fi
}


DATA_DIR=/prod/IR01/dewim/mps/human
FUIDS=${DATA_DIR}/stores/matedFragUIDs.txt
MS=mapSummary.txt
ASTORE=human_vanilla.asmStore


# input file is named all.regions
# parameter is assembly name
ASSEMBLY=${1}

if [ -z ${ASSEMBLY} ] || [ ${ASSEMBLY} == "bell-style" ]; then
  echo "Please specify an assembly name"
  return;
fi

MC=${ASSEMBLY}_matchCounts.cgm
UMF=${ASSEMBLY}_uniquelyMappedFragments.txt
MSI=${DATA_DIR}/${ASSEMBLY}/rawData/${ASSEMBLY}_mapStoreInput.txt
CUIDS=${DATA_DIR}/${ASSEMBLY}/rawData/${ASSEMBLY}_UIDs.txt

cd ${DATA_DIR}/${ASSEMBLY}/rawData

echo "Compiling match counts"
DeleteFile ${MC}
matchCounts < all.regions > ${MC}

echo "Identifying uniquely matching fragments"
DeleteFile ${UMF}
DeleteFile ${MS}
getUniqueMatchesFromCGM ${MC} ${FUIDS} all.regions > ${UMF} 2> ${MS}

echo "Building map store input & adjusting for off-by-one difference"
DeleteFile ${MSI}
gawk '{print $1, $3, $4-1, $5-1}' ${UMF} > ${MSI}

echo "Building map store"
cd ${DATA_DIR}/stores
buildMapStore -s ${ASTORE} -m ${ASSEMBLY}.mapStore -c ${CUIDS} -f ${MSI}
