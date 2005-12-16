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
# $Id: makeMockInterVANFiles.sh,v 1.5 2005-12-16 22:13:07 catmandew Exp $
#

if [ -z ${DATA_DIR} ]; then
  if [ ${OS} == "AIX" ] || [ ${OS} == "OSF1" ]; then
    export DATA_DIR=/prod/IR01/dewim/mps/human
  else
    export DATA_DIR=/home/dewim/celera/sandbox/cds/IR/COMPASS/data/human
  fi
fi

cd ${DATA_DIR}/vScaffolds/interScaffold
inFile=allInterScaffoldMPs.txt

# first, adjust for 13 split scaffolds
for ss in `gawk '{print $1}' ${DATA_DIR}/VANSplitPoints.txt`; do
  oldS=`grep ${ss} ${DATA_DIR}/VANSplitPoints.txt | gawk '{print $1}'`
  newLeftS=`grep ${ss} ${DATA_DIR}/VANSplitPoints.txt | gawk '{print $16}'`
  newRightS=`grep ${ss} ${DATA_DIR}/VANSplitPoints.txt | gawk '{print $19}'`
  coord=`grep ${ss} ${DATA_DIR}/VANSplitPoints.txt | gawk '{print $13}'`

  echo "Adjusting ${ss} into ${newLeftS} and ${newRightS} at ${coord}"

  gawk -v c=${coord} -v o=${oldS} -v l=${newLeftS} -v r=${newRightS} \
  '{ \
    if($2==o) \
    { \
      if($3<c) \
      { \
        $2=l; \
      } \
      else \
      { \
        $3-=c; \
        $2=r; \
      } \
    } \
    if($6==o) \
    { \
      if($7<c) \
      { \
        $6=l; \
      } \
      else \
      { \
        $7-=c; \
        $6=r; \
      } \
    } \
    print $0; \
  }' ${inFile} > ${newLeftS}.txtA
  
  inFile=${newLeftS}.txtA
done

cp ${inFile} newAllInterScaffoldMPs.txt;
inFile=${DATA_DIR}/vScaffolds/interScaffold/newAllInterScaffoldMPs.txt;

cd ${DATA_DIR}/scfCoordsVAN
# loop over all chromosome mapping files
for file in `ls chr*.out`; do

  # get the chromosome number
  chr=${file%%.*}
  chr=${chr#*r}
  declare -i c=0
  if [ ${chr} == "U" ]; then
    c=24;
  else
    c=${chr}-1
  fi
  echo "${chr} -> ${c}"

  # loop over all scaffolds in the mapping file
  for scf in `gawk '{print $1}' ${file}`; do
  
    destFile=${DATA_DIR}/vScaffolds/interScaffold/${scf}.txtB
    
    # if the scaffold was in the assembly, process it...
    left=`grep ${scf} ${file} | gawk '{print $3}'`
    right=`grep ${scf} ${file} | gawk '{print $4}'`
    orient=`grep ${scf} ${file} | gawk '{print $5}'`

    if [ ${orient} == "-1" ]; then
      gawk -v c=${c} -v s=${scf} -v r=${right} \
      '{ \
        if($2==s) \
        { \
          $2=c; \
          $3=r-$3; \
          if($4=="A_B") \
          { \
            $4="B_A"; \
          } \
          else \
          { \
            $4="A_B"; \
          } \
        } \
        if($6==s) \
        { \
          $6=c; \
          $7=r-$7; \
          if($8=="A_B") \
          { \
            $8="B_A"; \
          } \
          else \
          { \
            $8="A_B"; \
          } \
        } \
        print $0; \
      }' ${inFile} > ${destFile}
    else
      gawk -v c=${c} -v s=${scf} -v l=${left} \
      '{ \
        if($2==s) \
        { \
          $2=c; \
          $3=l+$3; \
        } \
        if($6==s) \
        { \
          $6=c; \
          $7=l+$7; \
        } \
        print $0; \
      }' ${inFile} > ${destFile}
    fi

    rm ${inFile}
    inFile=${destFile}
  done
done

cp ${inFile} ${DATA_DIR}/VAN_SCF/interChromosome/newAllInterChromosomeMPs.txt
inFile=${DATA_DIR}/VAN_SCF/interChromosome/newAllInterChromosomeMPs.txt

# create chromosome files
declare -i curr=0
declare -i last=24

cd ${DATA_DIR}/VAN_SCF/interChromosome/
while [ ${curr} -le ${last} ]; do
  gawk -v c=${curr} '{if($2==c){print $0}else if($6==c){print $5, $6, $7, $8, $1, $2, $3, $4, $9}}' ${inFile} > ${curr}.txt
  curr=${curr}+1
done

cd ${DATA_DIR}/VAN_SCF/interChromosome/
for file in `ls ?.txt`; do
  mv ${file} 0${file}
done

for file in `ls ??.txt`; do
  mv ${file} 0${file}
done

