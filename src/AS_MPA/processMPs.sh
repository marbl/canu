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
# $Id: processMPs.sh,v 1.3 2005-03-22 19:05:58 jason_miller Exp $
#

if [ ${OS} == "AIX" ] || [ ${OS} == "OSF1" ]; then
  export BASE_DIR=/prod/IR01/dewim/mps/
else
  export BASE_DIR=/home/dewim/celera/sandbox/cds/IR/COMPASS/data
fi
export DATA_DIR=${BASE_DIR}/human

Assemblies=(HG05 HG06 CSAB WGAB VAN VAN_asm B28 R26 R27 B33A B34)

# loop over all assemblies
for assembly in "${Assemblies[@]}"; do

  echo "Processing assembly ${assembly}"
  processAssemblyMPs.sh ${assembly}
done
