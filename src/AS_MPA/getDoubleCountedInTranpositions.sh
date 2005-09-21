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
# $Id: getDoubleCountedInTranpositions.sh,v 1.5 2005-09-21 20:13:07 catmandew Exp $
#

if [ -z ${1} ] || [ ${1} == "bell-style" ]; then
  echo "Please enter an assembly name"
  return
fi

AS=${1}

InnieATACode=mi

egrep ${InnieATACode} ${AS}.*.transposition.ata |sed 's/[-=]/ /g' |gawk '{print $4, $17, $20}' |gawk 'BEGIN{last=-1}{if($1!=last){if(last!=-1)printf("\n");last=$1;printf("%d %s %s", $1, $2, $3)}else{printf(" %s %s", $2, $3)}}' > transpInnies.txt

gawk '{for(n=2;n<=NF;n++){print $n}}' transpInnies.txt > transpInnieUIDs.txt

egrep ${InnieATACode} ${AS}.*.compressed.ata |sed 's/[-=]/ /g' |gawk '{print $4, $17, $20}' |gawk 'BEGIN{last=-1}{if($1!=last){if(last!=-1)printf("\n");last=$1;printf("%d %s %s", $1, $2, $3)}else{printf(" %s %s", $2, $3)}}' > compressedInnies.txt

egrep ${InnieATACode} ${AS}.*.stretched.ata |sed 's/[-=]/ /g' |gawk '{print $4, $17, $20}' |gawk 'BEGIN{last=-1}{if($1!=last){if(last!=-1)printf("\n");last=$1;printf("%d %s %s", $1, $2, $3)}else{printf(" %s %s", $2, $3)}}' > stretchedInnies.txt

s=`fgrep -f transpInnieUIDs.txt stretchedInnies.txt |sort -n |uniq |wc -l|gawk '{print $1}'`

echo "${s} stretched intervals also included in transpositions"

c=`fgrep -f transpInnieUIDs.txt compressedInnies.txt |sort -n |uniq |wc -l|gawk '{print $1}'`

echo "${c} compressed intervals also included in transpositions"
