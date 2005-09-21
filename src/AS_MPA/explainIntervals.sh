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
# $Id: explainIntervals.sh,v 1.5 2005-09-21 20:13:07 catmandew Exp $
#

# params: 1=assembly, 2=input file, 3=tempdir
function ProcessFiles
{
  echo "In ProcessFiles"
  
  Libs=(2k 10k 50k bacEnd)
  Unsat1=(compressed)
  Unsat2=(stretched inversion transposition)

  # separate into one file per chromosome
  for file in `ls [0-9][0-9][0-9].txt`; do
    chr=${file%%.*}
    IvFile=${3}/${chr}.intervals.txt
    
    echo "Working on ${chr}"

    # separate this chromosome's intervals from the rest
    echo "Separating chrom ${chr} from ${2} into ${IvFile}"
    gawk -v c=${chr} '{if($1==c)print $2, $3, NR, $0}' ${2} > ${IvFile}

    # loop over all clone library sizes
    for lib in "${Libs[@]}"; do
    
      # do 5sigmas
      #fn=${1}.${chr}.${lib}.5sigma.intervals.txt
      # make sure file of satisfied intervals of 'lib' mates is present
      #if [ ! -f ${fn} ]; then
      #  for libUID in `cat ${DATA_DIR}/libs/${lib}Libs.txt`; do
      #    mean=`grep ${libUID} ${DATA_DIR}/libs/humanLibs.txt | gawk '{print $2}'`
      #    stddev=`grep ${libUID} ${DATA_DIR}/libs/humanLibs.txt | gawk '{print $3}'`
      #    gawk -v l=${libUID} -v m=${mean} -v s=${stddev} 'BEGIN{min=m-5*s;max=m+5*s}{if($1=="I"&&$4==l&&$6-$5>=min&&$6-$5<=max)print $5, $6-$5}' ${chr}.txt >> ${fn}
      #  done
      #fi
      #echo "Intersecting with ${fn}"
      #DoIntersections ${IvFile} ${fn} ${3}
      
      # do satisfied intervals
      fn=${1}.${chr}.${lib}.intervals.txt
      # make sure file of satisfied intervals of 'lib' mates is present
      if [ ! -f ${fn} ] && [ -f ${1}.${chr}.satisfied.raw ] ; then
        fgrep -f ${DATA_DIR}/libs/${lib}Libs.txt ${1}.${chr}.satisfied.raw | gawk '{print $5, $6-$5}' > ${fn}
      fi
      echo "Intersecting with ${fn}"
      DoIntersections ${IvFile} ${fn} ${3}
    done

    # get confirmed stretched/comrpessed with a 2k or 10k mate pair in them
    fn=${1}.${chr}.2k10k.stretchedCompressed.intervals.txt
    if [ ! -f ${fn} ]; then
      sed 's/[=]/ /g' ${1}.${chr}.stretched.ata ${1}.${chr}.compressed.ata | \
        gawk \
          '{ \
            if($13=="/weight") \
            { \
              start=$6; \
              len=$7; \
              consider=1; \
            } \
            else \
            { \
              if(consider==1 && \
                 ($19==19866798939412 || \
                  $19==9000001400016 || \
                  $19==19866798939413 || \
                  $19==19866799876498 || \
                  $19==19866798939410 || \
                  $19==19866798939411 || \
                  $19==3000068808047)) \
              { \
                print start, len; \
                consider=0; \
              } \
            } \
          }' > ${fn}
    fi
    echo "Intersecting with ${fn}"
    DoIntersections ${IvFile} ${fn} ${3}

    # make sure files of unsatisfied intervals are present
    for type in "${Unsat1[@]}"; do
      fn=${1}.${chr}.${type}.intervals.txt
      if [ ! -f ${fn} ] && [ -f ${1}.${chr}.${type}.ata ] ; then
        grep weight ${1}.${chr}.${type}.ata | gawk '{print $6, $7}' > ${fn}
      fi
      echo "Intersecting with ${fn}"
      DoIntersections ${IvFile} ${fn} ${3}
    done
    for type in "${Unsat2[@]}"; do
      # do bad intervals
      fn=${1}.${chr}.${type}.intervals.txt
      if [ ! -f ${fn} ]  && [ -f ${1}.${chr}.${type}.ata ] ; then
        grep weight ${1}.${chr}.${type}.ata | gawk '{print $6, $7}' > ${fn}
      fi
      echo "Intersecting with ${fn}"
      DoIntersections ${IvFile} ${fn} ${3}

      # do breakpoint intervals
      fn=${1}.${chr}.${type}.breakpoints.txt
      if [ ! -f ${fn} ]  && [ -f ${1}.${chr}.${type}.ata ] ; then
        gawk '{if($1=="F" && NF==11)print $6, $7}' ${1}.${chr}.${type}.ata> ${fn}
      fi
      echo "Intersecting with ${fn}"
      DoIntersections ${IvFile} ${fn} ${3}
    done
    
    # put together a total file
    outFile=${3}/${chr}_all.csv
    echo "- = intersection of any type" > ${outFile}
    echo "i = intersects but doesn't span" >> ${outFile}
    echo "s = spans" >> ${outFile}
    echo "q = intersects but isn't spanned by" >> ${outFile}
    echo "5S = 5 stddevs is satisfied" >> ${outFile}
    echo "I = mis-assembled interval. Note transposition intervals may be huge & span good sequence." >> ${outFile}
    echo "BP = breakpoint intervals at ends of mis-assembled intervals" >> ${outFile}
    echo "Interval(start length), 10k, 10ki, 10ks, 10kq, 2k, 2ki, 2ks, 2kq, 2k10kConfirmedCompressedOrStretched, 2k10kConfirmedCompressedOrStretchedi, 2k10kConfirmedCompressedOrStretcheds, 2k10kConfirmedCompressedOrStretchedq, 50k, 50ki, 50ks, 50kq, BacEnd, BacEndi, BacEnds, BacEndq, Compressed, Compressedi, Compresseds, Compressedq, InversionBP, InversionBPi, InversionBPs, InversionBPq, InversionI, InversionIi, InversionIs, InversionIq, StretchedBP, StretchedBPi, StretchedBPs, StretchedBPq, StretchedI, StretchedIi, StretchedIs, StretchedIq, TranspositionBP, TranspositionBPi, TranspositionBPs, TranspositionBPq, TranspositionI, TranspositionIi, TranspositionIs, TranspositionIq" >> ${outFile}
    numFields=`gawk '{if(m<NF)m=NF}END{print m}' ${IvFile}`
    paste -d ' ' ${IvFile} ${3}/${1}.${chr}.*.intervals.txt.counts ${3}/${1}.${chr}.*.breakpoints.txt.counts | gawk -v f=${numFields} '{for(i=1;i<=f;i++){printf("%s ", $i)}for(i=f+1;i<=NF;i++)printf(", %s", $i);printf("\n")}' >> ${outFile}
  done
}

function DoIntersections
{
  getIntervalIntersections ${1} ${2} > ${3}/${2}.n
  gawk '{print $1, NF-1}' ${3}/${2}.n > ${3}/${2}.nc
  getIntervalIntersections ${1} ${2} -i > ${3}/${2}.i
  gawk '{print $1, NF-1}' ${3}/${2}.i > ${3}/${2}.ic
  getIntervalIntersections ${1} ${2} -s > ${3}/${2}.s
  gawk '{print $1, NF-1}' ${3}/${2}.s > ${3}/${2}.sc
  getIntervalIntersections ${1} ${2} -q > ${3}/${2}.q
  gawk '{print $1, NF-1}' ${3}/${2}.q > ${3}/${2}.qc
  paste -d ' ' ${3}/${2}.nc ${3}/${2}.ic ${3}/${2}.sc ${3}/${2}.qc | gawk '{print $2, $4, $6, $8}' > ${3}/${2}.counts
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
  exit
fi

if [ -z ${2} ] || [ ! -f ${2} ]; then
  echo "Please identify input filename as 2nd param"
  exit
fi

if [ -z ${3} ]; then
  echo "Please identify working subdirectory as 3rd param"
  exit
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

cd ${DATA_DIR}/${AS}/intraChromosome
mkdir ${TEMP_DIR}
ProcessFiles ${AS} ${IN_FILE} ${TEMP_DIR}

cd ${DATA_DIR}/${AS}/intraChromosome/unmapped
mkdir ${TEMP_DIR}
ProcessFiles ${AS} ${IN_FILE} ${TEMP_DIR}

cd ${currdir}
