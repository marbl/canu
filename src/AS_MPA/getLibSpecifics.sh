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
# $Id: getLibSpecifics.sh,v 1.4 2005-03-22 19:48:58 jason_miller Exp $
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

Type=(stretched compressed inversion transposition)

declare -i lcount=0
declare -i total=0

libFile=${DATA_DIR}/libs/humanLibs.txt

for t in "${Type[@]}"; do
  fn=${AS}_${t}_type.txt
  echo "${t}" > ${fn}
  for lib in `sort -k 2n ${libFile} | gawk '{print $1}'`; do
    mean=`egrep ${lib} ${libFile} | gawk '{print $2}'`
    lcount=0
    total=0
    for file in `ls [0-9][0-9][0-9].txt`; do
      chr=${file%%.*}
      lcount=`grep -c ${lib} ${AS}.${chr}.${t}.ata`
      total=${lcount}+${total}
    done
    echo "${lib} ${mean} ${total}" >> ${fn}
  done
done

paste -d ' ' ${AS}_*_type.txt | gawk '{if(NR==1){printf("Lib , mean");for(i=1;i<=NF;i++){printf(" , %s", $i)}}else{printf("%s , %s", $1, $2); for(i=3;i<=NF;i+=3)printf(" , %d", $i)}printf("\n")}' > ${AS}_libSummary.csv