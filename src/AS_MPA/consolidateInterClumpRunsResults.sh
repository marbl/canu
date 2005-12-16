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
# $Id: consolidateInterClumpRunsResults.sh,v 1.5 2005-12-16 22:13:07 catmandew Exp $
#

if [ -z ${DATA_DIR} ]; then
  if [ ${OS} == "AIX" ] || [ ${OS} == "OSF1" ]; then
    export DATA_DIR=/prod/IR01/dewim/mps/human
  else
    export DATA_DIR=/home/dewim/celera/sandbox/cds/IR/COMPASS/data/human
  fi
fi

echo "B34 Left Breakpoint, 10k, 2k, 50k, BAC End, CompressedBP, InversionBP, StretchedBP, TranspositionBP, B34 Left Interval, InterChromosome, B34 Right Breakpoint, 10k, 2k, 50k, BAC End, CompressedBP, InversionBP, StretchedBP, TranspositionBP, B34 Right Interval, InterChromosome, VAN Left Breakpoint, 10k, 2k, 50k, BAC End, CompressedBP, InversionBP, StretchedBP, TranspositionBP, VAN Left Interval, InterChromosome, VAN Right Breakpoint, 10k, 2k, 50k, BAC End, CompressedBP, InversionBP, StretchedBP, TranspositionBP, VAN Right Interval, InterChromosome, VAN Inter-clump Interval, 10k, 2k, 50k, BAC End, "
paste -d ' ' ${DATA_DIR}/B34/intraChromosome/B34LeftPoint/allNoHeaderSorted.csv \
  ${DATA_DIR}/B34/interChromosome/B34LeftInterval/allNoHeaderSorted.csv \
  ${DATA_DIR}/B34/intraChromosome/B34RightPoint/allNoHeaderSorted.csv \
  ${DATA_DIR}/B34/interChromosome/B34RightInterval/allNoHeaderSorted.csv \
  ${DATA_DIR}/VAN/intraChromosome/VANLeftPoint/allNoHeaderSorted.csv \
  ${DATA_DIR}/VAN/interChromosome/VANLeftInterval/allNoHeaderSorted.csv \
  ${DATA_DIR}/VAN/intraChromosome/VANRightPoint/allNoHeaderSorted.csv \
  ${DATA_DIR}/VAN/interChromosome/VANRightInterval/allNoHeaderSorted.csv \
  ${DATA_DIR}/VAN/intraChromosome/VANInterval/allNoHeaderSorted.csv
