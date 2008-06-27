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
# $Id: convertICRuns.sh,v 1.6 2008-06-27 06:29:17 brianwalenz Exp $
#


currdir=`pwd`
OutFile=B34_VAN_EndPointPairs.txt
OutFile2=B34_VAN_IntervalPairs.txt

grep -v "#" ${1} | sed 's/:/ /g' | \
  gawk \
  '{ \
     if(NR%2==1) \
     { \
       r1=$4; \
       bc1=$6; bs1=$7; bl1=$8; bo1=$9; \
       vc1=$11; vs1=$12; vl1=$13; vo1=$14; \
     } \
     else \
     { \
       r2=$4; \
       bc2=$6; bs2=$7; bl2=$8; bo2=$9; \
       vc2=$11; vs2=$12; vl2=$13; vo2=$14; \
       if(vs1<vs2) \
       { \
         print bc1, bs1+bl1, 0, vc1, vs1+vl1, 0, bc2, bs2, 0, vc2, vs2, 0, r1, r2; \
       } \
       else \
       { \
         print bc2, bs2+bl2, 0, vc2, vs2+bl2, 0, bc1, bs1, 0, vc1, vs1, 0, r2, r1; \
       } \
     } \
   }' > ${OutFile}

gawk '{print $1, $2, $3, $13}' ${OutFile} > B34LeftPoint.txt

gawk '{print $7, $8, $9, $14}' ${OutFile} > B34RightPoint.txt

gawk '{print $4, $5, $6, $13}' ${OutFile} > VANLeftPoint.txt

gawk '{print $10, $11, $12, $14}' ${OutFile} > VANRightPoint.txt

grep -v "#" ${1} | sed 's/:/ /g' | \
  gawk \
  '{ \
     if(NR%2==1) \
     { \
       r1=$4; \
       bc1=$6; bs1=$7; bl1=$8; bo1=$9; \
       vc1=$11; vs1=$12; vl1=$13; vo1=$14; \
     } \
     else \
     { \
       r2=$4; \
       bc2=$6; bs2=$7; bl2=$8; bo2=$9; \
       vc2=$11; vs2=$12; vl2=$13; vo2=$14; \
       if(vs1<vs2) \
       { \
         print bc1, bs1, bl1, vc1, vs1, vl1, bc2, bs2, bl2, vc2, vs2, vl2, r1, r2; \
       } \
       else \
       { \
         print bc2, bs2, bl2, vc2, vs2, vl2, bc1, bs1, bl1, vc1, vs1, vl1, r2, r1; \
       } \
     } \
   }' > ${OutFile2}


gawk '{print $1, $2, $3, $13}' ${OutFile2} > B34LeftInterval.txt

gawk '{print $7, $8, $9, $14}' ${OutFile2} > B34RightInterval.txt

gawk '{print $4, $5, $6, $13}' ${OutFile2} > VANLeftInterval.txt

gawk '{print $10, $11, $12, $14}' ${OutFile2} > VANRightInterval.txt

gawk '{print $4, $5+$6, $11-($5+$6), $13, $14}' ${OutFile2} > VANInterval.txt
