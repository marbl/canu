
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
/*********************************************************************
   CVS_ID:  $Id: ConsistencyChecksREZ.h,v 1.3 2005-03-22 19:07:33 jason_miller Exp $
 *********************************************************************/
#ifndef CONSISTENCY_CHECKS_H
#define CONSISTENCY_CHECKS_H

#include "DataTypesREZ.h"

// The following defines make it easy to adapt changes
// in the typing of variance and mean

#define AS_REZ_MEANTYPE int
#define AS_REZ_VARTYPE  float

// The following define specifies the percentage of bases
// two 3*stdDev intervals must overlap in order to pass 
// the consistency check. This is not yet implemented

//#define OLAP_CONSERVATIVE
// if OLAP_CONSERVATIVE is defined the parameters
// AS_REZ_MAX_REL_ERROR and AS_REZ_MIN_OVERLAP
// are defined more stringent

#define AS_REZ_ERROR_RATE 0.06

#ifdef OLAP_CONSERVATIVE
#define AS_REZ_MIN_OVERLAP 30
#define AS_REZ_MAX_REL_ERROR 0.05
#else
#define AS_REZ_MAX_REL_ERROR 2.0
#define AS_REZ_MIN_OVERLAP 30
#endif
// AS_REZ_MAX_REL_ERROR gives the maxmal relative error that a guessed overlap
// can have relativ to the computed overlap
// AS_REZ_MIN_OVERLAP gives the minium length which we consider to be an overlap


#define AS_REZ_SIMTEST 1000

int check_consistency(Scaffold_Fill_t *, 
		      int,		      
		      int);

void combine_two_distrib(LengthT,
			 LengthT,
			 LengthT *);

OverlapStatusREZ check_overlap(Gap_Chunk_t,
			       Gap_Chunk_t, 
			       int,
			       ChunkOverlapCheckT *);

int Is_Edge_Consistent(CIEdgeT *,
		       const Gap_Chunk_t *,
		       const Gap_Chunk_t *);

int Is_Edge_Orientation_Consistent(CIEdgeT *,
				   const Gap_Chunk_t *,
				   const Gap_Chunk_t *);

#endif








