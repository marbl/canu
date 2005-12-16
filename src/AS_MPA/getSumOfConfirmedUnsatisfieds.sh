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
# $Id: getSumOfConfirmedUnsatisfieds.sh,v 1.5 2005-12-16 22:13:07 catmandew Exp $
#

function DoIt
{
  confUns=`cat ${1}/intraChromosome/${1}.0[0-2]?.stretched.ata \
               ${1}/intraChromosome/${1}.0[0-2]?.compressed.ata \
               ${1}/intraChromosome/${1}.0[0-2]?.normal.ata \
               ${1}/intraChromosome/${1}.0[0-2]?.antinormal.ata \
               ${1}/intraChromosome/${1}.0[0-2]?.outtie.ata \
               ${1}/interChromosome/${1}.0[0-2]?.interChromosome.ata | \
             grep weight | \
               sed 's/=/ /'g | \
                 gawk '{a+=$(NF)}END{print a}'`
                 
  confMappable=`cat ${1}/interChromosome/unmapped/${1}.[0-9][0-9][0-9].interChromosome.ata | \
                  grep weight | \
                    sed 's/=/ /'g | \
                      gawk '{a+=$(NF)}END{print a}'`
  
  mapped=`cat ${1}/intraChromosome/[0-9][0-9][0-9].txt \
              ${1}/interChromosome/[0-9][0-9][0-9].txt | \
            gawk 'END{print NR}'`
            
  unmapped=`cat ${1}/intraChromosome/unmapped/[0-9][0-9][0-9].${2} \
                ${1}/interChromosome/unmapped/[0-9][0-9][0-9].txt | \
              gawk 'END{print NR}'`
              
  gawk -v a=${1} \
       -v s=${confUns} \
       -v t=${confMappable} \
       -v m=${mapped} \
       -v u=${unmapped} \
    'END{printf("%s: %d confirmed unsatisfied of %d matepairs in mapped sequence = %.3f%%. %d mappable unmapped matepairs of %d matepairs total = %.3f%%.\n", a, s, m, 100.0 * s/m, t, u+m, 100.0 * t/(u+m))}' ~/.bashrc
}


CeleraAssemblies=(WGAB CSAB VAN R26 R27)
PublicAssemblies=(HG05 HG06 B28 B33A B34)

for ass in "${CeleraAssemblies[@]}"; do
  DoIt ${ass} "all.txt"
done

for ass in "${PublicAssemblies[@]}"; do
  DoIt ${ass} "txt"    
done


