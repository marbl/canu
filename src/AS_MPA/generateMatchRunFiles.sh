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
# $Id: generateMatchRunFiles.sh,v 1.2 2004-09-23 20:25:24 mcschatz Exp $
#

function GenerateMatchRunFilesInDir
{
  for file in `ls [0-9][0-9][0-9].txt`; do
    chr=${file%%.*}
    for i in "${UnsatIntervals2[@]}"; do
      if [ -f ${1}.${chr}.${i}.ata ]; then
        gawk '{if($1=="F" && NF==11)print $6, $7}' ${1}.${chr}.${i}.ata > ${1}.${chr}.${i}.breakpoints.txt
      fi
    done
    if [ -f ${1}.${chr}.compressed.intervals.txt ]; then
      cp ${1}.${chr}.compressed.intervals.txt ${1}.${chr}.compressed.breakpoints.txt
    fi
    if [ ! -f ${1}.${chr}.satisfied.clones.txt ]; then
      gawk '{print $5, $6-$5}' ${1}.${chr}.satisfied.raw > ${1}.${chr}.satisfied.clones.txt
    fi
  done
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

UnsatIntervals2=(inversion stretched transposition)

currdir=`pwd`

cd ${DATA_DIR}/${AS}/intraChromosome
GenerateMatchRunFilesInDir ${AS}

cd unmapped
GenerateMatchRunFilesInDir ${AS}

cd ${currdir}