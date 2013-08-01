
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received (LICENSE.txt) a copy of the GNU General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/
static char *rcsid = "$Id: AS_UTL_interval.c,v 1.9 2009-07-30 10:42:56 brianwalenz Exp $";
#include "AS_global.H"
#include "AS_UTL_interval.H"

/* For two unordered intervals (a_bgn, a_end) (b_bgn, b_end), this function
   returns the length of their overlap, or 0 if they don't overlap.
   Only overlaps longer than minimum_overlap are considered.
*/
int32 IntervalsOverlap( int32 a_bgn, int32 a_end,
                              int32 b_bgn, int32 b_end,
                              int32 minimum_overlap){
  int32 a_min, b_min, a_max, b_max;
  int32 themin, themax;
  int32 overlap;

  /* For a linear genome, the fragments truly overlap if the length
     of a line segment that covers both of the fragments is less than
     the sum of the lengths of both fragments. */

  a_min = MIN(a_bgn,a_end);
  a_max = MAX(a_bgn,a_end);
  b_min = MIN(b_bgn,b_end);
  b_max = MAX(b_bgn,b_end);
  themin = MIN(a_min,b_min);
  themax = MAX(a_max,b_max);
  overlap = ((a_max-a_min)+(b_max-b_min)) - (themax-themin);

  if(overlap > minimum_overlap)
    return overlap;
  // else
    return 0;

}





