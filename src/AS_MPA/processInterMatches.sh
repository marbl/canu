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
# $Id: processInterMatches.sh,v 1.3 2005-03-22 19:05:58 jason_miller Exp $
#

delta=${1}
gapNs=${2}

if [ -z ${delta} ] || [ ${delta} == "bell-style" ]; then
  delta=500
fi

if [ -z ${gapNs} ]; then
  gapNs=20
fi

UnsatIntervals=(compressed stretched inversion transposition)
UnsatIntervals2=(inversion stretched transposition)

function ProcessAssemblyFiles
{
  currdir=`pwd`
  
  cd ${2}
  
  # extract the relevant assembly's intervals from Ms.atac
  # extract from the 'all' file to separate chromosome files for mapped
  # and get spanning clone type
  echo "mapped intervals"
  for file in `ls [0-9][0-9][0-9].txt`; do
    chr=${file%%.*}
    gawk -v c=${chr} '{if($1==c)print $2, $3, NR}' ${currdir}/${1}.all.txt > ${currdir}/${chr}.${1}.txt
    getIntervalIntersections ${currdir}/${chr}.${1}.txt ${1}.${chr}.satisfied.clones.txt -i -q > ${currdir}/${chr}.${1}.satisfied.notSpanned.txt
    for t in "${UnsatIntervals[@]}"; do
      getIntervalIntersections ${currdir}/${chr}.${1}.txt ${1}.${chr}.${t}.breakpoints.txt > ${currdir}/${chr}.${1}.${t}.bpIntersections.txt
    done
  done

  cd unmapped

  # extract from the 'all' file to separate chromosome files for unmapped
  # and get spanning clone type
  echo "unmapped intervals"
  for file in `ls [0-9][0-9][0-9].txt`; do
    chr=${file%%.*}
    gawk -v c=${chr} '{if($1==c)print $2, $3, NR}' ${currdir}/${1}.all.txt > ${currdir}/${chr}.${1}.txt
    getIntervalIntersections ${currdir}/${chr}.${1}.txt ${1}.${chr}.satisfied.clones.txt -i -q > ${currdir}/${chr}.${1}.satisfied.notSpanned.txt
    for t in "${UnsatIntervals[@]}"; do
      getIntervalIntersections ${currdir}/${chr}.${1}.txt ${1}.${chr}.${t}.breakpoints.txt > ${currdir}/${chr}.${1}.${t}.bpIntersections.txt
    done
  done

  cd ${currdir}

  echo "cat'ing, cut'ing, sort'ing, gawk'ing"
  cat *.${1}.satisfied.notSpanned.txt | cut -f 3- -d',' | sort -n | gawk '{print NF-1}' > ${1}.satisfied.notSpanned.sorted.txt
  for t in "${UnsatIntervals[@]}"; do
    cat *.${1}.${t}.bpIntersections.txt | cut -f 3- -d',' | sort -n | gawk '{print NF-1}' > ${1}.${t}.bpIntersections.sorted.txt
  done

  echo "paste'ing"
  paste -d ' ' ${1}.satisfied.notSpanned.sorted.txt ${1}.*.bpIntersections.sorted.txt | gawk '{printf("%d: %s\n", NR, $0)}' > ${1}.marked.txt
  
  # clean up
  #rm ${1}.all.txt
  #rm *.${1}.*.notSpanned.txt
  #rm ${1}.*.notSpanned.sorted.txt
}


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
O_DIR=${DATA_DIR}/one2one/interMatch

IN_FILE=${O_DIR}/pairid.ian.1.atac
IR_FILE=${O_DIR}/IRs.atac
NI_FILE=${O_DIR}/numberedIRs.atac

cd ${O_DIR}

echo "Filtering matches. MinDelta=${delta}, Gapped=${gapNs} contiguous Ns"
sed 's/[=]/ /g' ${IN_FILE} | \
  gawk -v d=${delta} -v g=${gapNs} \
    '{ \
      if($21<g && ($7-$11>d || $11-$7>d)) \
      { \
        for(i=1;i<=13;i++)printf("%s ", $i); \
        for(i=14;i<=NF;i+=2)printf("%s=%s ", $i, $(i+1)); \
        printf("\n") \
      } \
     }' > ${IR_FILE}
nl ${IR_FILE} > ${NI_FILE}

echo "Separating ${V_AS} intervals"
sed 's/:/ /g' ${IR_FILE} |gawk '{print $11, $12, $13}' > ${V_AS}.all.txt
echo "Separating ${B_AS} intervals"
sed 's/:/ /g' ${IR_FILE} |gawk '{print $6, $7, $8}' > ${B_AS}.all.txt

ProcessAssemblyFiles ${V_AS} ${V_DIR}
ProcessAssemblyFiles ${B_AS} ${B_DIR}

# MNum B34sat B34comp B34inv B34str B34trans VANcomp VANinv VANstr VANtrans
paste -d ' ' B34.marked.txt VAN.marked.txt |sed 's/:/ /g' | awk '{print $1, $2, $3, $4, $5, $6, $8, $9, $10, $11, $12}' > B34VAN.marked.txt

awk '{b=$3+$4+$5+$6;v=$8+$9+$10+$11; \
      if(b==0&&v==0) \
      { \
        ans="both okay"; \
      } \
      else if(b>0 && v==0) \
      { \
        ans="B34 bad"; \
      } \
      else if(b==0 && v > 0) \
      { \
        ans="VAN bad" \
      } \
      else \
      { \
        ans="both bad" \
      } \
      print $1, ans \
    }' B34VAN.marked.txt > B34VAN.results.txt

sed 's/ /,/g' B34VAN.marked.txt > B34VAN.marked.csv

paste -d ' ' B34VAN.results.txt ${IR_FILE} | awk -v a=${B_AS} '{if($2==a){for(i=4;i<=NF;i++){printf("%s ", $i)}printf("\n")}}' > ${B_AS}_BadMs.atac

paste -d ' ' B34VAN.results.txt ${IR_FILE} | awk -v a=${V_AS} '{if($2==a){for(i=4;i<=NF;i++){printf("%s ", $i)}printf("\n")}}' > ${V_AS}_BadMs.atac

paste -d ' ' B34VAN.marked.txt ${IR_FILE} | sed 's/:/ /g' | awk '{for(i=1;i<=11;i++){printf("%s , ",$i)}printf("%d , %d , %d , %d , %d , %d\n", $17+1, $18, $19, $22+1, $23, $24)}' > B34VAN.tally.csv
