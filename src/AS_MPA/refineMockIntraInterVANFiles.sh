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
# $Id: refineMockIntraInterVANFiles.sh,v 1.4 2005-03-22 19:48:58 jason_miller Exp $
#

if [ -z ${DATA_DIR} ]; then
  if [ ${OS} == "AIX" ] || [ ${OS} == "OSF1" ]; then
    export DATA_DIR=/prod/IR01/dewim/mps/human
  else
    export DATA_DIR=/home/dewim/celera/sandbox/cds/IR/COMPASS/data/human
  fi
fi


cd ${DATA_DIR}/VAN_SCF/interChromosome
# loop over all chromosome mapping files
for file in `ls [0-9][0-9][0-9].txt`; do

  # get the chromosome number
  chr=${file%%.*}

  echo "${chr}"

  # gawk '{if($2!=$6)print $0}' ${file} > ${chr}.inter.txt
  gawk \
    '{ \
       if($2==$6) \
       { \
         if($3<$7) \
         { \
           lUID=$1; \
           lc=$3;
           rUID=$5; \
           rc=$7;
           if($4=="A_B") \
           { \
             if($8=="B_A") \
             { \
               o="I"; \
             } \
             else
             { \
               o="N";
             } \
           } \
           else \
           { \
             if($8=="B_A") \
             { \
               o="A"; \
             } \
             else
             { \
               o="O";
             } \
           } \
         } \
         else \
         { \
           lUID=$5; \
           lc=$7;
           rUID=$1; \
           rc=$3;
           if($8=="A_B") \
           { \
             if($4=="B_A") \
             { \
               o="I"; \
             } \
             else
             { \
               o="N";
             } \
           } \
           else \
           { \
             if($4=="B_A") \
             { \
               o="A"; \
             } \
             else
             { \
               o="O";
             } \
           } \
         } \
         print o, lUID, rUID, $9, lc, rc;
       } \
     }' ${file} > ../intraChromosome/${chr}.intra.txt
  
done