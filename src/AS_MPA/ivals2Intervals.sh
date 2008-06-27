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
# $Id: ivals2Intervals.sh,v 1.6 2008-06-27 06:29:17 brianwalenz Exp $
#


function SeparateIntoChromFiles
{
  declare -i maxChrom=`gawk '{if($1>m)m=$1}END{print m}' ${1}`
  declare -i curr=0
  while [ ${curr} -le ${maxChrom} ]; do
    gawk -v c=${curr} '{if($1==c){for(i=2;i<=NF;i++){printf("%s ", $i)}printf("\n")}}' ${1} > ${curr}.${2}.txt
    curr=${curr}+1
  done
  for file in `ls ?.${2}.txt`; do
    mv ${file} 0${file}
  done
  for file in `ls ??.${2}.txt`; do
    mv ${file} 0${file}
  done
}

# create 4 parallel interval files
echo "Creating parallel interval files"
sed 's/:/ /g' ${1} | gawk '{if($6==$20){print $6, $21, $7-$21, NR}else{print $6, -1, -1, NR}}' > B34Left.txt

sed 's/:/ /g' ${1} | gawk '{if($6==$24){print $6, $9, $25-$9, NR}else{print $6, -1, -1, NR}}' > B34Right.txt

sed 's/:/ /g' ${1} | gawk '{if($12==$30){print $12, $31, $13-$31, NR}else{print $12, -1, -1, NR}}' > VANLeft.txt

sed 's/:/ /g' ${1} | gawk '{if($12==$34){print $12, $15, $35-$15, NR}else{print $12, -1, -1, NR}}' > VANRight.txt

# separate into chrom files
echo "Separating B34Left into chromosomes"
SeparateIntoChromFiles B34Left.txt BL
echo "Separating B34Right into chromosomes"
SeparateIntoChromFiles B34Right.txt BR

echo "Separating VANLeft into chromosomes"
SeparateIntoChromFiles VANLeft.txt VL
echo "Separating VANRight into chromosomes"
SeparateIntoChromFiles VANRight.txt VR
