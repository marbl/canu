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
# $Id: processHiLoInterChromIntervalFiles.sh,v 1.2 2004-09-23 20:25:24 mcschatz Exp $
#


if [ -z ${DATA_DIR} ]; then
  if [ ${OS} == "AIX" ] || [ ${OS} == "OSF1" ]; then
    export DATA_DIR=/prod/IR01/dewim/mps/human
  else
    export DATA_DIR=/home/dewim/celera/sandbox/cds/IR/COMPASS/data/human
  fi
fi

AS=${1}
mapping=${2}

gawk '{ \
  l=$5+1;
  r=$8+1;
  if(l>a) a=l; \
  if(r>a) a=r; \
  if($2<2) \
  { \
    if($4<2) \
    { \
      unkI[l]++; \
      unkBP[l] += $7; \
      unkI[r]++; \
      unkBP[r] += $10; \
    } \
    else
    { \
      movingI[l]++; \
      movingBP[l] += $7; \
      pullingI[r]++; \
      pulledBP[r] += $7; \
    } \
  } \
  else if($4<2) \
  { \
    pullingI[l]++; \
    pulledBP[l] += $10; \
    movingI[r]++; \
    movingBP[r] += $10; \
  } \
  else \
  { \
    rptI[l]++; \
    rptBP[l] += $7; \
    rptI[r]++; \
    rptBP[r] += $10; \
  } \
}END{ \
  print "Chrom , pullingI, pulledBP, movingI, movingBP, rptI, rptBP, unkI, unkBP"; \
  for(i=1; i<= a; i++) \
  { \
    printf("%d , %d , %d , %d , %d , %d , %d , %d , %d\n", \
           i, \
           pullingI[i], pulledBP[i], \
           movingI[i], movingBP[i], \
           rptI[i], rptBP[i], \
           unkI[i], unkBP[i]) \
  } \
}' hi_loChromI.txt > ${AS}_${mapping}InterSummary.csv
