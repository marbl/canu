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
# $Id: doRunsSpanningTranspositions.sh,v 1.6 2008-06-27 06:29:17 brianwalenz Exp $
#

if [ -z ${DATA_DIR} ]; then
  if [ ${OS} == "AIX" ] || [ ${OS} == "OSF1" ]; then
    export DATA_DIR=/prod/IR01/dewim/mps/human
  else
    export DATA_DIR=/home/dewim/celera/sandbox/cds/IR/COMPASS/data/human
  fi
fi



function ProcessAssemblyDir
{
  for file in `ls [0-9][0-9][0-9].txt`; do
    blah=${file%%.*}
    gawk '{if($2=="tl" || $2=="tr")print $6, $7}' ${1}.${blah}.transposition.ata > ${blah}.transpositionBPs.txt
    gawk -v c=${blah} '{if($1==c)print $2, $3}' ${2} > ${blah}.common.runs.txt

    getIntervalIntersections ${blah}.transpositionBPs.txt ${blah}.common.runs.txt -s | \
      gawk -v c=${blah} 'BEGIN{a=0}{if(NF>1)a++}END{print c, a}' >> ${OutFile}
  done
}

function ProcessAssembly
{
  # do chromosomes 0 to 23
  cd ${DATA_DIR}/${1}/intraChromosome
  ProcessAssemblyDir ${1} ${2}

  cd ${DATA_DIR}/${1}/intraChromosome/unmapped
  ProcessAssemblyDir ${1} ${2}
}


Assemblies=(VAN B34)

currdir=`pwd`

for Ass in "${Assemblies[@]}"; do
  echo "Running ${Ass}"

  OutFile=${DATA_DIR}/${Ass}.CommonRunTranspBP.txt
  if [ -f ${OutFile} ]; then
    rm ${OutFile}
  fi

  ProcessAssembly ${Ass} ${DATA_DIR}/atac/${Ass}CommonRuns.txt
done

cd ${currdir}