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
# $Id: findBacDeletions.sh,v 1.1.1.1 2004-04-14 13:52:00 catmandew Exp $
#

function ProcessAssemblyFiles
{
  # extract the relevant assembly's intervals from Ms.atac
  sed 's/:/ /g' Ms.atac |gawk -v f=${3} '{print $(f), $(f+1), $(f+2)}' > ${1}.all.txt

  # extract from the 'all' file to separate chromosome files for mapped
  # and get spanning clone type
  for file in `ls ${2}/[0-9][0-9][0-9].txt`; do
    chr=${file%%.*}
    chr=${chr##*/}
    gawk -v c=${chr} '{if($1==c)print $2, $3, NR}' ${1}.all.txt > ${chr}.${1}.txt
    getIntervalIntersections ${chr}.${1}.txt ${2}/${1}.${chr}.${4}.${5}.txt -s > ${chr}.${1}.spanned.txt
  done
  
  # extract from the 'all' file to separate chromosome files for unmapped
  # and get spanning clone type
  for file in `ls ${2}/unmapped/[0-9][0-9][0-9].txt`; do
    chr=${file%%.*}
    chr=${chr##*/}
    gawk -v c=${chr} '{if($1==c)print $2, $3, NR}' ${1}.all.txt > ${chr}.${1}.txt
    getIntervalIntersections ${chr}.${1}.txt ${2}/unmapped/${1}.${chr}.${4}.${5}.txt -s > ${chr}.${1}.spanned.txt
  done

  # get interval # & any spanning clone intervals
  cat *.${1}.spanned.txt |cut -f 3- -d ',' |sort -n > ${1}.spanned.sorted.txt

  # convert to interval # & a count of spanning clone intervals
  gawk '{print $1, NF-1}' ${1}.spanned.sorted.txt > ${1}.marked.txt

  # clean up
  #rm ${1}.all.txt
  #rm *.${1}.spanned.txt
  #rm ${1}.spanned.sorted.txt
}

if [ -z ${DATA_DIR} ]; then
  if [ ${OS} == "AIX" ] || [ ${OS} == "OSF1" ]; then
    export DATA_DIR=/prod/IR01/dewim/mps/human
  else
    export DATA_DIR=/home/dewim/celera/sandbox/cds/IR/COMPASS/data/human
  fi
fi

declare -i MinVanSatClones=1
declare -i LastVChrom=24
declare -i LastBChrom=40

B_AS=B34
V_AS=VAN
BD_DIR=${DATA_DIR}/bacDeletions/${V_AS}vs${B_AS}wGaps
V_DIR=${DATA_DIR}/${V_AS}/intraChromosome
B_DIR=${DATA_DIR}/${B_AS}/intraChromosome

cd ${BD_DIR}

# egrep "vanNs=0" pairid.D.atac | gawk '{if($1=="M")print $0}' > Ms.atac
gawk '{if($1=="M")print $0}' pairid.D.atac > Ms.atac
nl Ms.atac > numberedMs.atac

# parameters are assembly, directory, first field in .atac, clone type, type
ProcessAssemblyFiles ${V_AS} ${V_DIR} 11 satisfied clones
ProcessAssemblyFiles ${B_AS} ${B_DIR}  6 compressed intervals

# rm Ms.atac

paste -d ' ' ${V_AS}.marked.txt ${B_AS}.marked.txt | sed 's/://g' > joined.marked.txt   
paste -d ' ' joined.marked.txt numberedMs.atac | \
  gawk -v g=${MinVanSatClones} \
    '{ \
       if($2>=g && $4>0) \
       { \
         for(i=6;i<=NF;i++){printf("%s ", $i)} \
         printf("\n") \
       } \
     }' > potentialBacDeletionMs.txt

gawk -v vas=${V_AS} -v bas=${B_AS} \
  '{ \
     a++; \
     b+=$7; \
     v+=$11; \
   }END{printf("%d BAC deletions. %d bp in %s. %d bp in %s.\n", a, v, vas, b, bas) \
   }' potentialBacDeletionMs.txt > bacDeletionSummary.txt
