
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
/* 	$Id: AS_CGW_dataTypes.h,v 1.1.1.1 2004-04-14 13:50:12 catmandew Exp $	 */
#ifndef AS_CGW_DATATYPES_H
#define AS_CGW_DATATYPES_H

#undef DEBUG
#include <assert.h>
#include "AS_global.h"
#include "AS_UTL_Var.h"
#ifdef NEVER
#include "AS_UTL_histo.h"
#endif
#include "math.h"

/*** Constants ***/
#define CHECK_CONNECTIVITY TRUE

#define CGW_MIN_READS_IN_UNIQUE  2

#define CGW_FUDGE_FACTOR (0.026)
#define CGW_DP_ERATE .10
#define CGW_DP_THRESH 1e-6
#define CGW_DP_MINLEN 20
#define CGW_DP_TRY_HARDER_MINLEN 20
#define CGW_DP_DESPERATION_MINLEN 10

// Due to the FBAC fragments, we get some pathologically short U-Unitigs
// Set the following threshhold to eliminate the really short ones
#define CGW_MIN_DISCRIMINATOR_UNIQUE_LENGTH 1000


#define NO_END 0
#define A_END 1
#define B_END 2
#define ALL_END (A_END | B_END)

#include "PrimitiveVA.h"

/* We've moved to a linear model of variance as a function of length.
   We chose the FUDGE_FACTOR as follows:
  3 * sqrt(variance(600bp fragment)) = 2% of 600bp   3-sigma = 2% of 600bp
  variance(600bp fragment) = (.02 * 600)^2 /3 = FUDGE_FACTOR * 600
*/

static double ComputeFudgeVariance(double length){
  //  double variance = length * CGW_FUDGE_FACTOR/3.0;
  //  variance *= variance;
  // return variance;

  return fabs(length) * CGW_FUDGE_FACTOR;

}

#define ALIGNER_FUNC Overlap *(*)(char *, char *, int, int, int, double, \
                                 double, int, CompareOptions)

// statistics from building sedges
typedef struct
{
  cds_int32 edgesAttempted;
  cds_int32 edgesSucceeded;
  cds_int32 edgesInternal;
  cds_int32 guidesAttempted;
  cds_int32 guidesSucceeded;
  cds_int32 guidesInternal;
} SEdgeBuildStats;

#endif
