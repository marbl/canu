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
# $Id: getCloneLengths.sh,v 1.6 2005-09-21 20:13:07 catmandew Exp $
#

prefix=${1}

if [ -z ${prefix} ] || [ ${prefix} == "bell-style" ] ; then
  echo "Please specify the prefix of the intra & lib files"
  return
fi

for lib in `gawk '{print $1}' ${prefix}Libs.txt`; do
  fn=${lib}Lengths.txt
  if [ -f $fn ] ; then
    rm -f $fn;
  fi
  echo "Creating file ${fn}"
  
  for file in `ls ${prefix}_*_intra.txt`; do
    gawk -v l=${lib} '{if($1=="I" && $4==l){printf("%u ", $6-$5); if($2<$3){printf("%u %u\n", $2, $3)}else{printf("%u %u\n", $3, $2)}}}' ${file} >> ${lib}Lengths.txt
  done

  cut -f 1 -d ' ' ${lib}Lengths.txt | sort -n |uniq -c | gawk '{a+=$1;printf("%u %u\n", $2, a)}' > ${lib}LengthsCumulative.gp
  
done