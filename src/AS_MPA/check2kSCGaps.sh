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
# $Id: check2kSCGaps.sh,v 1.3 2005-03-22 19:05:47 jason_miller Exp $
#


function ProcessDirFiles
{
  for file in `ls [0-9][0-9][0-9].txt`; do
    chr=${file%%.*}

    echo "Processing ${chr}"

    for type in "${Types[@]}"; do
    
      sed 's/=/ /g' ${AS}.${chr}.${type}.ata | gawk '{if($13=="/weight"){start=$6;len=$7;test=1}else{if(test==1){if($18=="/libUID" && ($19=="19866798939412"||$19=="9000001400016"||$19=="19866798939410"||$19=="19866798939413"||$19=="19866799876498"||$19=="19866798939411"||$19=="3000068808047")){print start, len;test=0}}}}' > ${AS}.${chr}.2k10k.${type}.intervals.txt
      getIntervalIntersections ${DATA_DIR}/${AS}/gaps/${chr}.intraLongScaffoldGaps.txt ${AS}.${chr}.2k10k.${type}.intervals.txt | gawk -v c=${chr} -v t=${type} 'BEGIN{a=0}{if(NF>1)a++}END{print c, t, a, "long scaffold gaps spanned"}'
      getIntervalIntersections ${DATA_DIR}/${AS}/gaps/${chr}.intraShortScaffoldGaps.txt ${AS}.${chr}.2k10k.${type}.intervals.txt | gawk -v c=${chr} -v t=${type} 'BEGIN{a=0}{if(NF>1)a++}END{print c, t, a, "short scaffold gaps spanned"}'
      
    done
    
    sed 's/=/ /g' ${AS}.${chr}.compressed.ata ${AS}.${chr}.stretched.ata | gawk '{if($13=="/weight"){start=$6;len=$7;test=1}else{if(test==1){if($18=="/libUID" && ($19=="19866798939412"||$19=="9000001400016"||$19=="19866798939410"||$19=="19866798939413"||$19=="19866799876498"||$19=="19866798939411"||$19=="3000068808047")){print start, len;test=0}}}}' > ${AS}.${chr}.2k10k.sc.intervals.txt
      getIntervalIntersections ${DATA_DIR}/${AS}/gaps/${chr}.intraLongScaffoldGaps.txt ${AS}.${chr}.2k10k.sc.intervals.txt | gawk -v c=${chr} -v t="stretched/compressed" 'BEGIN{a=0}{if(NF>1)a++}END{print c, t, a, "long scaffold gaps spanned"}'
      getIntervalIntersections ${DATA_DIR}/${AS}/gaps/${chr}.intraShortScaffoldGaps.txt ${AS}.${chr}.2k10k.sc.intervals.txt | gawk -v c=${chr} -v t="stretched/compressed" 'BEGIN{a=0}{if(NF>1)a++}END{print c, t, a, "short scaffold gaps spanned"}'
      
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
  exit
fi

Lib2ks=(19866798939412 9000001400016 19866798939410)
Lib10ks=(19866798939413 19866799876498 19866798939411 3000068808047)
Types=(compressed stretched)

cd ${DATA_DIR}/${AS}/intraChromosome
ProcessDirFiles

cd ${DATA_DIR}/${AS}/intraChromosome/unmapped
ProcessDirFiles
