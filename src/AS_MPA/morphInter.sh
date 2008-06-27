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
# $Id: morphInter.sh,v 1.6 2008-06-27 06:29:17 brianwalenz Exp $
#

cat *.elsewhere.txt | \
  gawk ' \
  { \
    if($2<$6) \
    { \
      if($4=="A_B") \
      { \
        if($8=="B_A") \
        { \
          printf("I") \
        } \
        else \
        { \
          printf("N") \
        } \
      } \
      else \
      { \
        if($8=="B_A") \
        { \
          printf("A") \
        } \
        else \
        { \
          printf("O") \
        } \
      } \
      printf(" %s %s %s %.f %.f\n", \
             $1, $5, $9, $2 * 1000000000 + $3, $6 * 1000000000 + $7) \
    } \
  } \
  ' > e.txt
