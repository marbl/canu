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
# $Id: processAssemblyMPs.sh,v 1.4 2005-03-22 19:48:58 jason_miller Exp $
#

######################################################################
# 'main' function

AS=${1}

if [ -z ${DATA_DIR} ]; then
  if [ ${OS} == "AIX" ] || [ ${OS} == "OSF1" ]; then
    export DATA_DIR=/prod/IR01/dewim/mps/human
  else
    export DATA_DIR=/home/dewim/celera/sandbox/cds/IR/COMPASS/data/human
  fi
fi

if [ -z ${AS} ] || [ ${AS} == "bell-style" ] ; then
  echo "Please identify the assembly (B33A, VAN, HG06, ...)"
  return
fi

echo "  analyzing mapped intra-chromosome mate pairs"
cd ${DATA_DIR}/${AS}/intraChromosome
processIntraMPs.sh ${AS}

echo "  analyzing unmapped intra-chromosome mate pairs"
cd ${DATA_DIR}/${AS}/intraChromosome/unmapped
processIntraMPs.sh ${AS}

echo "  analyzing mapped inter-chromosome mate pairs"
cd ${DATA_DIR}/${AS}/interChromosome
processInterMPsATA.sh ${AS} mapped
  
echo "  analyzing unmapped inter-chromosome mate pairs"
cd ${DATA_DIR}/${AS}/interChromosome/unmapped
processInterMPsATA.sh ${AS} unmapped
######################################################################
