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
# $Id: getCloneLengths.sh,v 1.3 2005-03-22 19:05:53 jason_miller Exp $
#

if [ -z ${DATA_DIR} ]; then
  if [ ${OS} == "AIX" ] || [ ${OS} == "OSF1" ]; then
    export DATA_DIR=/prod/IR01/dewim/mps/human
  else
    export DATA_DIR=/home/dewim/celera/sandbox/cds/IR/COMPASS/data/human
  fi
fi

cd ${DATA_DIR}/B34/intraChromosome

if [ -f ${DATA_DIR}/libs/${lib}NewLengths.txt ]; then
  rm ${DATA_DIR}/libs/${lib}NewLengths.txt
fi

for file in `ls [0-9][0-9][0-9].txt`; do
  echo "Processing clones in ${file}"
  for lib in `gawk '{print $1}' ${DATA_DIR}/libs/humanLibs.txt`; do
    gawk -v l=${lib} '{if($1=="I" && $4==l)print $6-$5}' ${file} >> ${DATA_DIR}/libs/${lib}NewLengths.txt
  done
done