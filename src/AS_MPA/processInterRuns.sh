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
# $Id: processInterRuns.sh,v 1.1.1.1 2004-04-14 13:52:06 catmandew Exp $
#

if [ -z ${DATA_DIR} ]; then
  if [ ${OS} == "AIX" ] || [ ${OS} == "OSF1" ]; then
    export DATA_DIR=/prod/IR01/dewim/mps/human
  else
    export DATA_DIR=/home/dewim/celera/sandbox/cds/IR/COMPASS/data/human
  fi
fi

B_AS=B34
V_AS=VAN
B_DIR=${DATA_DIR}/${B_AS}/intraChromosome
V_DIR=${DATA_DIR}/${V_AS}/intraChromosome
O_DIR=${DATA_DIR}/one2one/interRun

UnsatIntervals=(compressed stretched inversion transposition)

function ProcessDirFiles
{
  for file in `ls [0-9][0-9][0-9].txt`; do
    chr=${file%%.*}
    gawk -v c=${chr} '{if($1==c)print $2, $3, NR}' ${currdir}/${1}.all.txt > ${currdir}/${chr}.${1}.txt
    getIntervalIntersections ${currdir}/${chr}.${1}.txt ${1}.${chr}.satisfied.clones.txt -s > ${currdir}/${chr}.${1}.satisfied.spanned.txt
    for t in "${UnsatIntervals[@]}"; do
      getIntervalIntersections ${currdir}/${chr}.${1}.txt ${1}.${chr}.${t}.breakpoints.txt > ${currdir}/${chr}.${1}.${t}.bpIntersections.txt
    done
  done
}

function ProcessAssemblyFiles
{
  currdir=`pwd`
  
  cd ${2}
  echo "mapped intervals"
  ProcessDirFiles ${1} ${2}
  
  cd unmapped
  echo "unmapped intervals"
  ProcessDirFiles ${1} ${2}

  cd ${currdir}
  
  echo "cat'ing, cut'ing, sort'ing, gawk'ing"
  cat *.${1}.satisfied.spanned.txt | cut -f 3- -d',' | sort -n | gawk '{print NF-1}' > ${1}.satisfied.spanned.sorted.txt
  for t in "${UnsatIntervals[@]}"; do
    cat *.${1}.${t}.bpIntersections.txt | cut -f 3- -d',' | sort -n | gawk '{print NF-1}' > ${1}.${t}.bpIntersections.sorted.txt
  done

  echo "paste'ing"
  paste -d ' ' ${1}.satisfied.spanned.sorted.txt ${1}.*.bpIntersections.sorted.txt | gawk '{printf("%d: %s\n", NR, $0)}' > ${1}.marked.txt


  gawk '{ \
      b=$3+$4+$5+$6; \
      if(b==0) \
      { \
        if($2==0) \
        { \
          printf("%s no evidence\n", $1) \
        } \
        else \
        { \
          printf("%s good\n", $1) \
        } \
      } \
      else \
      { \
        if($3==b) \
        { \
          printf("%s compressed %d (%d)\n", $1, b, $2); \
        } \
        else if($4==b) \
        { \
          printf("%s inversion %d (%d)\n", $1, b, $2); \
        } \
        else if($5==b) \
        { \
          printf("%s stretched %d (%d)\n", $1, b, $2); \
        } \
        else if($6==b || b==$3+$5+$6) \
        { \
          printf("%s transposition %d (%d)\n", $1, b, $2); \
        } \
        else \
        { \
          printf("%s mixed %d (%d)\n", $1, b, $2); \
        } \
      } \
    }' ${1}.marked.txt > ${1}.results.txt

  if [ -f ${1}.summary.txt ]; then
    rm ${1}.summary.txt
  fi

  Categories=(good evidence compressed inversion stretched transposition mixed)
  declare -i r=0
  for c in "${Categories[@]}"; do
    r=`grep -c ${c} ${1}.results.txt`
    echo "${r} ${c}" >> ${1}.summary.txt
  done
  
  # clean up
  #rm ${1}.all.txt
  #rm *.${1}.*.spanned.txt
  #rm ${1}.*.spanned.sorted.txt
}

cd ${O_DIR}/${B_AS}
sed 's/[:-]/ /g' ${B_AS}.notinrun_noends.txt | awk '{print $1, $2, $4}' > ${B_AS}.all.txt
ProcessAssemblyFiles ${B_AS} ${B_DIR}

cd ${O_DIR}/${V_AS}
sed 's/[:-]/ /g' ${V_AS}.notinrun_noends.txt | awk '{print $1, $2, $4}' > ${V_AS}.all.txt
ProcessAssemblyFiles ${V_AS} ${V_DIR}

