
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
/* 	$Id: AS_CGW_dataTypes.h,v 1.12 2007-08-26 10:11:00 brianwalenz Exp $	 */
#ifndef AS_CGW_DATATYPES_H
#define AS_CGW_DATATYPES_H

#include <assert.h>
#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "math.h"

#define CHECK_CONNECTIVITY TRUE

#define CGW_MIN_READS_IN_UNIQUE  2

//  We've moved to a linear model of variance as a function of length.
//  We chose the FUDGE_FACTOR as follows:
//
//  3 * sqrt(variance(600bp fragment)) = 2% of 600bp
//  3-sigma = 2% of 600bp
//  variance(600bp fragment) = (.02 * 600)^2 /3 = FUDGE_FACTOR * 600
//
#define CGW_FUDGE_FACTOR           (0.026)
#define ComputeFudgeVariance(L)    (fabs(L) * CGW_FUDGE_FACTOR)

#define CGW_DP_THRESH 1e-6
#define CGW_DP_MINLEN 20
#define CGW_DP_DESPERATION_MINLEN 10

// Due to the FBAC fragments, we get some pathologically short U-Unitigs
// Set the following threshhold to eliminate the really short ones
#define CGW_MIN_DISCRIMINATOR_UNIQUE_LENGTH 1000

#define NO_END 0
#define A_END 1
#define B_END 2
#define ALL_END (A_END | B_END)

#endif
