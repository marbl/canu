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
# $Id: filterIntersWithSatisfieds.sh,v 1.1.1.1 2004-04-14 13:52:00 catmandew Exp $
#


# satfile, chromosome, hi or lo chromosome, and inter-dir
function DoIntersections
{
  gawk -v c=${2} '{if($1==c)print $2, $3, NR}' ${4}/${3}.txt > ${4}/${2}.${3}.txt
  getIntervalIntersections ${4}/${2}.${3}.txt ${1} -i -q | gawk '{print $1, NF-1}' | sed 's/[,:]/ /g' |gawk '{print $3, $4}' > ${4}/${2}.${3}I.txt
}


# parameter is inter-dir
function DoAllIntersections
{
  for file in `ls [0-9][0-9][0-9].txt`; do
    chr=${file%%.*}
    satFile=${AS}.${chr}.satisfied.clones.txt
  
    DoIntersections ${satFile} ${chr} hiChrom ${1}
    DoIntersections ${satFile} ${chr} loChrom ${1}
  done
}


# make sure all required satisfied clone files exist
function MakeSatCloneFiles
{
  for file in `ls [0-9][0-9][0-9].txt`; do
    chr=${file%%.*}
    satFile=${1}.${chr}.satisfied.clones.txt

    if [ ! -f ${satFile} ]; then
      gawk '{print $5, $6-$5}' ${1}.${chr}.satisfied.raw > ${satFile}
    fi
  done
}


if [ -z ${DATA_DIR} ]; then
  if [ ${OS} == "AIX" ] || [ ${OS} == "OSF1" ]; then
    export DATA_DIR=/prod/IR01/dewim/mps/human
  else
    export DATA_DIR=/home/dewim/celera/sandbox/cds/IR/COMPASS/data/human
  fi
fi

AS=${1}
mapping=${2}

if [ -z ${AS} ] || [ ${AS} == "bell-style" ]; then
  echo "Please identify assembly name"
  return 
fi

# work in current directory
currDir=`pwd`

# see if we're in an unmapped dir
declare -i unmapped=0
if [ -f 024.txt ]; then
  unmapped=1
fi

# cut hi chromosome intervals into hiChrom.txt
cat ${AS}.???.interChromosome.ata | sed 's/:/ /g' | gawk '{if($1=="M")print $6, $7, $8}' > hiChrom.txt

# cut lo chromosome intervals into loChrom.txt
cat ${AS}.???.interChromosome.ata | sed 's/:/ /g' | gawk '{if($1=="M")print $11, $12, $13}' > loChrom.txt

# change to mapped intrachromosome dir
cd ${DATA_DIR}/${AS}/intraChromosome
MakeSatCloneFiles ${AS}

# intersect with all satisfieds, mapped chromosome by chromosome
DoAllIntersections ${currDir}

# change to mapped intrachromosome dir
if [ ${unmapped} -eq 1 ]; then
  cd ${DATA_DIR}/${AS}/intraChromosome/unmapped
  MakeSatCloneFiles ${AS}

  cd ${DATA_DIR}/${AS}/intraChromosome
  DoAllIntersections ${currDir}

  cd ${DATA_DIR}/${AS}/intraChromosome/unmapped
  DoAllIntersections ${currDir}
fi

# go to the inter-dir & continue with processing
cd ${currDir}

cat [0-9][0-9][0-9].hiChromI.txt | sort -n > all.hiChromI.txt
cat [0-9][0-9][0-9].loChromI.txt | sort -n > all.loChromI.txt

paste all.hiChromI.txt all.loChromI.txt > pasted.hi_loChromI.txt

paste pasted.hi_loChromI.txt hiChrom.txt > deleteme.txt
paste deleteme.txt loChrom.txt > hi_loChromI.txt

processHiLoInterChromIntervalFiles.sh ${AS}
