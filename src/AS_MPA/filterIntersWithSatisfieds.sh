#!/bin/bash
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
# $Id: filterIntersWithSatisfieds.sh,v 1.5 2005-09-21 20:13:07 catmandew Exp $
#


# satfile, chromosome, hi or lo chromosome
function DoIntersections
{
  gawk -v c=${3} '{if($1==c)print $2, $3, NR}' ${1}.${4}.txt > ${1}.${3}.${4}.txt
  getIntervalIntersections ${1}.${3}.${4}.txt ${2} -i -q | gawk '{print $1, NF-1}' | sed 's/[,:]/ /g' |gawk '{print $3, $4}' > ${1}.${3}.${4}I.txt 2> /dev/null
}


function DoAllIntersections
{
  for file in `ls | egrep "${1}_(.+)_inter.txt"`; do
    chr=${file%_*}
    chr=${chr##*_}
    satFile=${1}.${chr}.satisfied.clones.txt
  
    DoIntersections ${1} ${satFile} ${chr} hiChrom
    DoIntersections ${1} ${satFile} ${chr} loChrom
  done
}


# make sure all required satisfied clone files exist
function MakeSatCloneFiles
{
  for file in `ls | egrep "${1}_(.+)_inter.txt"`; do
    chrom=${file%_*}
    chrom=${chrom##*_}
    satFile=${1}.${chrom}.satisfied.clones.txt
    inFile=${1}.${chrom}.satisfied.raw

    if [ ! -f ${satFile} ] && [ -f ${inFile} ] ; then
      gawk '{print $5, $6-$5}' ${inFile} > ${satFile}
    fi
  done
}


AS=${1}
mapping=${2}

if [ -z ${AS} ] || [ ${AS} == "bell-style" ]; then
  echo "Please identify assembly name"
  return 
fi

outFile="${AS}.hiChrom.txt"
# cut hi chromosome intervals into hiChrom.txt
if [ -f ${outFile} ] ; then
  rm -f ${outFile}
fi
for file in `ls | egrep "${AS}\.(.+)\.interChromosome.ata"`; do
  cat ${file} | \
    sed 's/:/ /g' | \
    gawk '{if($1=="M")print $6, $7, $8}' >> ${outFile}
done

outFile="${AS}.loChrom.txt"
# cut lo chromosome intervals into loChrom.txt
if [ -f ${outFile} ] ; then
  rm -f ${outFile}
fi
for file in `ls | egrep "${AS}\.(.+)\.interChromosome.ata"`; do
  cat ${file} | \
    sed 's/:/ /g' | \
    gawk '{if($1=="M")print $11, $12, $13}' >> ${outFile}
done

# change to mapped intrachromosome dir
MakeSatCloneFiles ${AS}

# intersect with all satisfieds, mapped chromosome by chromosome
DoAllIntersections ${AS}

# change to mapped intrachromosome dir
if [ "${mapping}" == "unmapped" ]; then
  cd unmapped
  MakeSatCloneFiles ${AS}
  cd ..

  DoAllIntersections ${AS}

  cd unmapped
  DoAllIntersections ${AS}
  cd ..
fi

outFile="${AS}.all.hiChromI.txt"
if [ -f ${outFile} ] ; then
  rm -f ${outFile}
fi
for file in `ls | egrep "${AS}\.(.+)\.hiChromI.txt"`; do
  cat ${file}
done | sort -n > ${outFile}

outFile="${AS}.all.loChromI.txt"
if [ -f ${outFile} ] ; then
  rm -f ${outFile}
fi
for file in `ls | egrep "${AS}\.(.+)\.loChromI.txt"`; do
  cat ${file}
done | sort -n > ${outFile}

paste ${AS}.all.hiChromI.txt ${AS}.all.loChromI.txt > ${AS}.pasted.hi_loChromI.txt

paste ${AS}.pasted.hi_loChromI.txt ${AS}.hiChrom.txt > deleteme.txt
paste deleteme.txt ${AS}.loChrom.txt > ${AS}.hi_loChromI.txt
rm -f deleteme.txt

processHiLoInterChromIntervalFiles.sh ${AS}
