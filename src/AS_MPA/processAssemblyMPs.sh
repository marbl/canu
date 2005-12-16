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
# $Id: processAssemblyMPs.sh,v 1.6 2005-12-16 22:13:07 catmandew Exp $
#

######################################################################
# 'main' function

AS=${1}
LIBFILE=${2}

if [ -z ${AS} ] || [ ${AS} == "bell-style" ] ; then
  echo "Please identify the assembly as 1st parameter"
  echo "Please identify the clone library stats file as 2nd parameter"
  return
fi

if [ -z ${LIBFILE} ] || [ ! -f ${LIBFILE} ] ; then
  echo "Please identify the assembly as 1st parameter"
  echo "Please identify the clone library stats file as 2nd parameter"
  return
fi

echo "  analyzing mapped intra-chromosome mate pairs"
processIntraMPs.sh ${AS} ${LIBFILE}

if [ -d "unmapped" ] ; then
  echo "  analyzing unmapped intra-chromosome mate pairs"
  cd unmapped
  processIntraMPs.sh ${AS} ${LIBFILE}
  cd ..
fi

echo "  analyzing mapped inter-chromosome mate pairs"
processInterMPsATA.sh ${AS} ${LIBFILE}
  
if [ -d "unmapped" ] ; then
  echo "  analyzing unmapped inter-chromosome mate pairs"
  cd unmapped
  processInterMPsATA.sh ${AS} ${LIBFILE} unmapped
fi
######################################################################
