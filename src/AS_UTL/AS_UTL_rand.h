
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
#ifndef AS_UTL_RAND_H
#define AS_UTL_RAND_H
/* A normally distributed random number generator with zero mean and unit variance (stddev)
   From "Numerical Recipes in C", pg 217

   Uses drand48 for a uniform random deviate, so srand48(seed) should be called to initialize.

*/

/** For drand48 **/
#include <stdlib.h>
#include "AS_global.h"

double GaussRandomNormalized_AS(void);

static void InitRandom_AS(long seed){
  srand48(seed);
}

static double GaussRandom_AS(double mean, double stddev){
  return GaussRandomNormalized_AS() * stddev + mean;
}


// Return a random number distributed:
//  uniform = TRUE    uniformly on [min,max]
//          = FALSE   normally on mean = (max + min)/2  std = (max - min)/6
//
int GetRand_AS(int min, int max, int uniform);
// returns a uniformly distributed random number from the interval
// [min,max]
double GetDrand_AS(double min, double max);

#endif
