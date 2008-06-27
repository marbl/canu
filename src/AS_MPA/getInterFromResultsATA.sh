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
# $Id: getInterFromResultsATA.sh,v 1.6 2008-06-27 06:29:17 brianwalenz Exp $
#

chroms=(0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23)

echo "chromosome , mis-intervals , mis-bps , trans-intervals , trans-bps"
for chrom in "${chroms[@]}"; do
  cat *.[0-9][0-9][0-9].interChromosome.ata | sed 's/:/ /g' | \
    gawk -v c=${chrom} \
      'BEGIN{from=1;to=2;} \
      { \
        if($2=="xa") \
        { \
          if($17=="/mixedOrientations=0") \
          { \
            mixed=0; \
          } \
          else \
          { \
            mixed=1; \
          } \
        } \
        else if($2=="xb") \
        { \
          if(mixed==1) \
          { \
            if($6==c) \
            { \
              mixedI[from]++; \
              mixedBP[from] += $8; \
            } \
            else if($11==c) \
            { \
              mixedI[to]++; \
              mixedBP[to] += $8; \
            } \
          } \
          else \
          { \
            if($6==c) \
            { \
              pureI[from]++; \
              pureBP[from] += $8; \
            } \
            else if($11==c) \
            { \
              pureI[to]++; \
              pureBP[to] += $8; \
            } \
          } \
        } \
        else if($2=="xc") \
        { \
          if(mixed==1) \
          { \
            if($6==c) \
            { \
              mixedI[to]++; \
              mixedBP[to] += $8; \
            } \
            else if($11==c) \
            { \
              mixedI[from]++; \
              mixedBP[from] += $8; \
            } \
          } \
          else \
          { \
            if($6==c) \
            { \
              pureI[to]++; \
              pureBP[to] += $8; \
            } \
            else if($11==c) \
            { \
              pureI[from]++; \
              pureBP[from] += $8; \
            } \
          } \
        } \
      } \
      END{printf("%d , %d , %d , %d , %d\n", c+1, \
                 (pureI[from] + pureI[to]) / 2, \
                 (pureBP[from] + pureBP[to]) / 2, \
                 (mixedI[from] + mixedI[to]) / 2, \
                 (mixedBP[from] + mixedBP[to]) / 2)}'
done
