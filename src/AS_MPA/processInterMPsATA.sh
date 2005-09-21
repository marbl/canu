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
# $Id: processInterMPsATA.sh,v 1.5 2005-09-21 20:13:07 catmandew Exp $
#

AS=${1}
LIBFILE=${2}

if [ -z ${AS} ] || [ ${AS} == "bell-style" ] ; then
  echo "Please identify the assembly (B33A, VAN, HG06, ...) as parameter 1"
  echo "Optional 3rd parameter: 'unmapped'"
  return
fi

if [ -z ${LIBFILE} ] || [ ! -f ${LIBFILE} ] ; then
  echo "Please identify the assembly (B33A, VAN, HG06, ...) as parameter 1"
  echo "Please identify the library stats file as parameter 2"
  echo "Optional 3rd parameter: 'unmapped'"
  return
fi

for file in `ls | egrep "${AS}_(.+)_inter.txt"`; do
  chrom=${file%_*}
  chrom=${chrom##*_}
  echo -n -e "    working on ${chrom}                   \r"
  processInterCG -l ${LIBFILE} -e ${file} 2> ${AS}.${chrom}.err
done

echo -e "\n\n  summarizing inter-chromosome results"
filterIntersWithSatisfieds.sh ${AS} ${mapping}