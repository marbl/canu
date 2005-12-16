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
# $Id: analyzeGaps.sh,v 1.5 2005-12-16 22:13:07 catmandew Exp $
#

if [ -z ${DATA_DIR} ]; then
  if [ ${OS} == "AIX" ] || [ ${OS} == "OSF1" ]; then
    export DATA_DIR=/prod/IR01/dewim/mps/human
  else
    export DATA_DIR=/home/dewim/celera/sandbox/cds/IR/COMPASS/data/human
  fi
fi


AS=${1}

if [ -z ${AS} ] || [ ${AS} == "bell-style" ]; then
  echo "Please identify assembly name"
  return 
fi

Intervals=(stretched compressed)
Libs=(2k 10k 50k bacEnd)

# want to find:
#  stretched intervals that span gaps
#  compressed intervals that span gaps
#  intra-scaffold gaps spanned by satisfied 2k/10k clones
#  inter-scaffold gaps spanned by satisfied BACs

# do stretched & compressed intervals that span gaps
echo "Identifying compressed/stretched clones that span 20bp gaps"
for i in "${Intervals[@]}"; do
  if [ -f ${AS}_${i}_NGaps.txt ]; then
    rm ${AS}_${i}_NGaps.txt
  fi
  for file in `ls [0-9][0-9][0-9].txt`; do
    c=${file%%.*}
    egrep weight ${AS}.${c}.${i}.ata |gawk '{print $6, $7}' > ${AS}.${c}.${i}.intervals.txt
    getIntervalIntersections ${AS}.${c}.${i}.intervals.txt ${DATA_DIR}/${AS}/gaps/${c}.allGaps.txt > ${AS}.${c}.${i}.NGaps.txt
#    egrep "\(20\)" ${AS}.${c}.${i}.NGaps.txt |wc |gawk -v c=${c} -v i=${i} '{print $1, i, "intervals spanning 20bp gaps in chrom", c}' >> ${AS}_${i}_NGaps.txt
    sed 's/,/ /g' ${AS}.${c}.${i}.NGaps.txt | gawk -v c=${c} -v i=${i} 'BEGIN{a=0}{for(j=3;j<=NF;j+=2){if($(j+1)-$(j)==20){a++;break}}}END{print a, i, "intervals spanning 20bp gaps in chrom", c}' >> ${AS}_${i}_NGaps.txt
  done
done

LibSets=(2k10k 50k bacEnd)

# get the gaps & satisfied clones
if [ -f ${DATA_DIR}/${AS}/gaps/000.intraScaffoldGaps.txt ]; then
  echo "Identifying intra-scaffold gaps spanned by 2k/10k clones"
  for file in `ls [0-9][0-9][0-9].txt`; do
    c=${file%%.*}

    # get satisfied clone intervals for other reasons
    gawk '{print $5, $6-$5}' ${AS}.${c}.satisfied.raw > ${AS}.${c}.satisfied.clones.txt

    for lib in "${LibSets[@]}"; do
      # get satisfied intervals
      egrep -f ${DATA_DIR}/libs/${lib}Libs.txt ${AS}.${c}.satisfied.raw | gawk '{print $5,$6-$5}' > ${AS}.${c}.${lib}.intervals.txt

      # do the intersections
      getIntervalIntersections ${DATA_DIR}/${AS}/gaps/${c}.intraScaffoldGaps.txt ${AS}.${c}.${lib}.intervals.txt -s > ${AS}.${c}.${lib}.intraSG.txt

      # generate a 2-column count file from each of the intersection results
      awk '{print $1, NF-1}' ${AS}.${c}.${lib}.intraSG.txt > ${AS}.${c}.${lib}.intraSGC.txt
    done

    # paste the files together & get rid of every other column
    paste -d ' ' ${AS}.${c}.*.intraSGC.txt | gawk '{printf("%s ", $1);for(i=2;i<=NF;i+=2){printf("%s ", $i)}printf("\n");}' > ${AS}.${c}.intraSGSummary.txt
  done
fi

# do full analysis of intra-scaffold gaps
# spanned by satisfied 2ks/10ks, 50ks, bac ends
# intersecting compressed, inversion, stretched, transposition





# do inter-scaffold gaps spanned by BAC end clones
if [ -f ${DATA_DIR}/${AS}/gaps/000.interScaffoldGaps.txt ]; then
  echo "Identifying inter-scaffold gaps spanned by BACs"
  for file in `ls [0-9][0-9][0-9].txt`; do
    c=${file%%.*}
    getIntervalIntersections ${DATA_DIR}/${AS}/gaps/${c}.interScaffoldGaps.txt ${AS}.${c}.bacEnd.intervals.txt > ${AS}.${c}.bacEnd.interSG.txt
  done
  
  # summarize bac ends
  cat ${AS}.*.bacEnd.interSG.txt | gawk '{a++;if(NF>1)b++;}END{printf("%d gaps, %d spanned by BAC end clones\n", a, b)}' > bacEnd.interSGSummary.txt
fi
