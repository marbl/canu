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
# $Id: processIntraMPs.sh,v 1.7 2008-06-27 06:29:17 brianwalenz Exp $
#

AS=${1}
LIBFILE=${2}

if [ -z ${AS} ] || [ ${AS} == "bell-style" ] ; then
  echo "Please identify the assembly name as 1st parameter"
  echo "Please identify the library stats file as 2nd parameter"
  return
fi

if [ -z ${LIBFILE} ] || [ ! -f ${LIBFILE} ] ; then
  echo "Please identify the assembly name as 1st parameter"
  echo "Please identify the library stats file as 2nd parameter"
  return
fi

rm *.err
if [ -f stretchedInnies.txt ]; then
  rm compressedInnies.txt stretchedInnies.txt transpInnie*.txt
fi

# process each file
for file in `ls | egrep "${AS}_(.+)_intra.txt"`; do
  chrom=${file%_*}
  chrom=${chrom##*_}
  echo -n -e "    working on ${chrom}                   \r"
  processIntra -l ${LIBFILE} -m ${file} 2> ${AS}.${chrom}.intra.err
done
echo -e "\n\n"

# # summarize results
# echo "  summarizing intra-chromosome results"
getIntraResults.sh ${AS} ${LIBFILE}

# # get Ns in stretched & compressed intervals
# echo "  analyzing sequence gaps"
# analyzeGaps.sh ${AS}

# # identify double-counted mate pairs in transpositions
echo "  identifying double-counted compressed & stretched in transpositions"
getDoubleCountedInTranpositions.sh ${AS} > ${AS}_doubleCountedInTranspositions.txt

