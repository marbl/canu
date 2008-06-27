
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
/********************************************************************/
/* A math utilities package for the Celera Assembler
 *
 *     Clark M. Mobarry
 *     April 1999
 *
 * $Id: math_AS.h,v 1.5 2008-06-27 06:29:21 brianwalenz Exp $
 */

#ifndef AS_UTL_MATH_H
#define AS_UTL_MATH_H

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <limits.h>

static int ceil_log2(const size_t x) {
  /*
   * Return the smallest integer K such that (x <= (1 << K)).
   * The name refers to the ceiling of the logarithm of the argument.
   * This function is useful for memory allocation where the
   * request is rounded up to the smallest power of two that
   * is larger than or equal to x.
   */

  size_t y;
  int k;
  /* Even though the logarithm of zero is not defined... */
  assert(x > 0);
  y = x-1; /* This is so that the half open interval
	      ( ((1L<<k-1)+1), (1L<<k) ] returns the same
	      value of k. */
  /* Remember that there are 8 bits per byte. */
  for(k=0; k < CHAR_BIT*sizeof(size_t) ; k++) {
    if(y == 0) break;
    y = y >> 1;
  }
  //assert(x <= (1L << k));
  //assert(((k == 0)&&(x<2)) || ((1L << k-1) < x));
  return k;
}

#ifdef NEVER
#include <stdio.h>
void main(void) {
    size_t x;
    for(x=1;x < 400; x++) {
      printf("x,ceil_log2(x)=" F_SIZE_T ",%d\n",x,ceil_log2(x));
    }
  }
#endif

#endif // AS_UTL_MATH_H

