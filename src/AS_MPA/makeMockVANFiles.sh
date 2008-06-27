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
# $Id: makeMockVANFiles.sh,v 1.6 2008-06-27 06:29:17 brianwalenz Exp $
#

if [ -z ${DATA_DIR} ]; then
  if [ ${OS} == "AIX" ] || [ ${OS} == "OSF1" ]; then
    export DATA_DIR=/prod/IR01/dewim/mps/human
  else
    export DATA_DIR=/home/dewim/celera/sandbox/cds/IR/COMPASS/data/human
  fi
fi

cd ${DATA_DIR}/vScaffolds/intraScaffold

# create intra-scaffold files for split scaffolds
for ss in `gawk '{print $1}' ${DATA_DIR}/VANSplitPoints.txt`; do
  newLeftS=`grep ${ss} ${DATA_DIR}/VANSplitPoints.txt | gawk '{print $16}'`
  newRightS=`grep ${ss} ${DATA_DIR}/VANSplitPoints.txt | gawk '{print $19}'`
  coord=`grep ${ss} ${DATA_DIR}/VANSplitPoints.txt | gawk '{print $13}'`

  echo "Splitting ${ss} into ${newLeftS} and ${newRightS} at ${coord}"

  inFile=${ss}.txt
  outLeft=${newLeftS}.txt
  outRight=${newRightS}.txt
  outTween=${newLeftS}.elsewheres.txt

  gawk -v c=${coord} '{if($5<c&&$6<c){print $0}}' ${inFile} > ${outLeft}
  gawk -v c=${coord} '{if($5>=c&&$6>=c){$5-=c;$6-=c;print $0}}' ${inFile} > ${outRight}
  gawk -v c=${coord} -v l=${newLeftS} -v r=${newRightS} \
    '{ \
      if($1=="I") \
      { \
        lo="A_B";ro="B_A"; \
      } \
      else if($1=="O") \
      { \
        lo="B_A";ro="A_B"; \
      }else if($1=="A") \
      { \
        lo="B_A";ro="B_A"; \
      }else \
      { \
        lo="A_B"; ro="A_B"; \
      } \
      if($5<c && $6>=c) \
      { \
        print $2, l, $5, lo, $3, r, $6-c, ro, $4;
      }else if($5>=c && $6<c) \
      { \
        print $2, r, $5-c, lo, $3, l, $6, ro, $4;
      } \
    }' ${inFile} > ${outTween}
done


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
  destFile=${DATA_DIR}/VAN_SCF/intraChromosome/${c}.txt

  if [ -f ${destFile} ]; then
    rm ${destFile}
  fi

  # loop over all scaffolds in the mapping file
  for scf in `gawk '{print $1}' ${file}`; do
    # if the scaffold was in the assembly, process it...
    left=`grep ${scf} ${file} | gawk '{print $3}'`
    right=`grep ${scf} ${file} | gawk '{print $4}'`
    if [ -f ${DATA_DIR}/vScaffolds/intraScaffold/${scf}.txt ]; then

      # if the scaffold is reversed, complement coordinates & swap frags
      orient=`grep ${scf} ${file} | gawk '{print $5}'`
      if [ ${orient} == "-1" ]; then
        gawk -v r=${right} '{print $1, $3, $2, $4, r-$6, r-$5}' ${DATA_DIR}/vScaffolds/intraScaffold/${scf}.txt >> ${destFile}
      else
        gawk -v l=${left} '{print $1, $2, $3, $4, l+$5, l+$6}' ${DATA_DIR}/vScaffolds/intraScaffold/${scf}.txt >> ${destFile}
      fi
    else
      declare -i length=${right}-${left}
      echo "Scaffold ${scf} ( ${length} bp) absent from original assembly!"
    fi

  done
done

cd ${DATA_DIR}/VAN_SCF/intraChromosome/
for file in `ls ?.txt`; do
  mv ${file} 0${file}
done

for file in `ls ??.txt`; do
  mv ${file} 0${file}
done

