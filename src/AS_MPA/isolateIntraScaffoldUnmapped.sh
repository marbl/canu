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
# $Id: isolateIntraScaffoldUnmapped.sh,v 1.1.1.1 2004-04-14 13:52:03 catmandew Exp $
#

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

cd ${DATA_DIR}/${AS}/intraChromosome/unmapped

for file in `ls [0-9][0-9][0-9].txt`; do
  chr=${file%%.*}
  gawk '{print $5, $6-$5, $0}' ${file} > process.txt
  mv ${file} ${chr}.all.txt

  gawk '{if($4>$3){print $3, $4-$3}else{print $4, $3-$4}}' ${DATA_DIR}/${AS}/scaffoldCoordinates/${chr}.out > ${chr}.scaffoldIntervals.txt

  getIntervalIntersections process.txt ${chr}.scaffoldIntervals.txt -s |gawk '{if(NF>1)print $1}' |sed 's/[,:]/ /g' |gawk '{for(i=3;i<=NF;i++)printf("%s ", $i);printf("\n")}' > ${file}

done