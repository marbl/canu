
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
static char CM_ID[] = "$Id: AS_UTL_interval.c,v 1.5 2006-11-14 17:52:18 eliv Exp $";
#include "AS_global.h"
#include "AS_UTL_interval.h"

/* For two unordered intervals (a_bgn, a_end) (b_bgn, b_end), this function
   returns the length of their overlap, or 0 if they don't overlap.
   Only overlaps longer than minimum_overlap are considered.
*/
CDS_COORD_t IntervalsOverlap( CDS_COORD_t a_bgn, CDS_COORD_t a_end,
                              CDS_COORD_t b_bgn, CDS_COORD_t b_end,
                              CDS_COORD_t minimum_overlap){
  CDS_COORD_t a_min, b_min, a_max, b_max;
  CDS_COORD_t themin, themax;
  CDS_COORD_t overlap;

  /* For a linear genome, the fragments truly overlap if the length
     of a line segment that covers both of the fragments is less than
     the sum of the lengths of both fragments. */

  a_min = min(a_bgn,a_end);
  a_max = MAX(a_bgn,a_end);
  b_min = min(b_bgn,b_end);
  b_max = MAX(b_bgn,b_end);
  themin = min(a_min,b_min);
  themax = MAX(a_max,b_max);
  overlap = ((a_max-a_min)+(b_max-b_min)) - (themax-themin);

  if(overlap > minimum_overlap)
    return overlap;
  // else
    return 0;

}





