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
# $Id: explainIntervalsInterChrom.sh,v 1.4 2005-03-22 19:48:58 jason_miller Exp $
#

# params: 1=assembly, 2=input file, 3=tempdir
function ProcessFiles
{
  echo "In ProcessFiles"
  
  # separate into one file per chromosome
  for file in `ls [0-9][0-9][0-9].txt`; do
    chr=${file%%.*}
    IvFile=${3}/${chr}.intervals.txt
    
    echo "Working on ${chr}"

    # separate this chromosome's intervals from the rest
    echo "Separating chrom ${chr} from ${2} into ${IvFile}"
    gawk -v c=${chr} '{if($1==c)print $2, $3, NR, $0}' ${2} > ${IvFile}

    # do interchrom intervals
    fn=${chr}.interChromosomeIntervals.txt
    echo "Intersecting with ${fn}"
    DoIntersections ${IvFile} ${fn} ${3}
    
    # put together a total file
    outFile=${3}/${chr}_all.csv
    echo "BP = breakpoint intervals at ends of mis-assembled intervals" > ${outFile}
    echo "Interval(start length), InterChromIntervals" >> ${outFile}
    numFields=`gawk '{if(m<NF)m=NF}END{print m}' ${IvFile}`
    paste -d ' ' ${IvFile} ${3}/${chr}.interChromosomeIntervals.txt.counts | gawk -v f=${numFields} '{for(i=1;i<=f;i++){printf("%s ", $i)}printf(", ");for(i=f+1;i<=NF;i++)printf("%s, ", $i);printf("\n")}' >> ${outFile}
  done
}

function DoIntersections
{
  getIntervalIntersections ${1} ${2} > ${3}/${2}.i
  gawk '{print $1, NF-1}' ${3}/${2}.i > ${3}/${2}.ic
  gawk '{print $2}' ${3}/${2}.ic > ${3}/${2}.counts
}

# parameters are 1: assembly, 2: input file, 3: workSubDir

# assume input is a file of intervals with fields (w/o commas)
# chromosome#, start, length

if [ -z ${DATA_DIR} ]; then
  if [ ${OS} == "AIX" ] || [ ${OS} == "OSF1" ]; then
    export DATA_DIR=/prod/IR01/dewim/mps/human
  else
    export DATA_DIR=/home/dewim/celera/sandbox/cds/IR/COMPASS/data/human
  fi
fi

if [ -z ${1} ] || [ ${1} == "bell-style" ]; then
  echo "Please identify assembly name as 1st param"
  return
fi

if [ -z ${2} ] || [ ! -f ${2} ]; then
  echo "Please identify input filename as 2nd param"
  return
fi

if [ -z ${3} ]; then
  echo "Please identify working subdirectory as 3rd param"
  return
fi

currdir=`pwd`

AS=${1}
IN_FILE=${currdir}/${2}
TEMP_DIR=${3}

echo "Input file is ${IN_FILE}"

# run on assembly intraChromosome & unmapped

# where X == 2k, 10k 50k, bacEnd, breakpoint intervals, bad intervals
# X intersecting but not spanning & not spanned by
# X spanning
# X spanned by

# want table with 3 * 6 columns per interval
# loop over all clone library sizes

cd ${DATA_DIR}/${AS}/interChromosome
mkdir ${TEMP_DIR}
ProcessFiles ${AS} ${IN_FILE} ${TEMP_DIR}

cd ${DATA_DIR}/${AS}/interChromosome/unmapped
mkdir ${TEMP_DIR}
ProcessFiles ${AS} ${IN_FILE} ${TEMP_DIR}

# combine files in mapped & unmapped
cd ${DATA_DIR}/${AS}/interChromosome/${TEMP_DIR}
cat *_all.csv ../unmapped/${TEMP_DIR}/*_all.csv | grep -v -i intervals | sort -k 3n > allNoHeaderSorted.csv

cd ${currdir}
