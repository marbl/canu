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
# $Id: getInterFromResults.sh,v 1.4 2005-03-22 19:48:58 jason_miller Exp $
#

chroms=(0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23)

echo "chromosome , intervals , bps"
for chrom in "${chroms[@]}"; do
  egrep " mps in " ???.out | gawk -v c=${chrom} '{if($11==c){a++;b+=$13-$12}else if($5==c){a++;b+=$7-$6}}END{printf("%d , %d , %d\n", c+1, a, b)}'
done

egrep " mps in " ???.out | gawk '{print $5, $6, $7-$6, $11, $12, $13-$12; print $11, $12, $13-$12, $5, $6, $7-$6}' > interChromosomeIntervals.txt

declare -i firstChrom=`gawk 'BEGIN{a=5000}{if($1<a)a=$1}END{print a}' interChromosomeIntervals.txt`

declare -i lastChrom=`gawk 'BEGIN{a=-1}{if($1>a)a=$1}END{print a}' interChromosomeIntervals.txt`

while [ ${firstChrom} -le ${lastChrom} ]; do
  gawk -v c=${firstChrom} '{if($1==c){for(i=2;i<=NF;i++){printf("%s ", $i)}printf("\n")}}' interChromosomeIntervals.txt > ${firstChrom}.interChromosomeIntervals.txt
  firstChrom=${firstChrom}+1
done

for file in `ls ?.interChromosomeIntervals.txt`; do
  mv ${file} 0${file}
done

for file in `ls ??.interChromosomeIntervals.txt`; do
  mv ${file} 0${file}
done
