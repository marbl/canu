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
# $Id: processInterClumpRunsFile.sh,v 1.2 2004-09-23 20:25:24 mcschatz Exp $
#

if [ -z ${DATA_DIR} ]; then
  if [ ${OS} == "AIX" ] || [ ${OS} == "OSF1" ]; then
    export DATA_DIR=/prod/IR01/dewim/mps/human
  else
    export DATA_DIR=/home/dewim/celera/sandbox/cds/IR/COMPASS/data/human
  fi
fi

explainBreakpoints.sh B34 B34LeftPoint.txt B34LeftPoint
explainBreakpoints.sh B34 B34RightPoint.txt B34RightPoint
explainBreakpoints.sh VAN VANLeftPoint.txt VANLeftPoint
explainBreakpoints.sh VAN VANRightPoint.txt VANRightPoint

explainIntervalsInterChrom.sh B34 B34LeftInterval.txt B34LeftInterval
explainIntervalsInterChrom.sh B34 B34RightInterval.txt B34RightInterval
explainIntervalsInterChrom.sh VAN VANLeftInterval.txt VANLeftInterval
explainIntervalsInterChrom.sh VAN VANRightInterval.txt VANRightInterval
explainIntervalsInterChrom2.sh VAN VANInterval.txt VANInterval

consolidateInterClumpRunsResults.sh > B34VANresults.csv
