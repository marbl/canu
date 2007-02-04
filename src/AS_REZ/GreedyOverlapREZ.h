
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

#ifndef GREEDY_OVERLAP_H
#define GREEDY_OVERLAP_H

#include "AS_global.h"
#include "AS_MSG_pmesg.h"
#include "ScaffoldGraph_CGW.h"

/* The main quality function. It assigns a value between 0 (good) and
 * 1 (bad) to an overlap between two fragments. NOTE that it is also
 * used to compute the overlap between ANY two sequences with quality
 * values. It returns the number of matches
 */
int
compute_bayesian_quality(InternalFragMesg* IFG1,
                         InternalFragMesg* IFG2,
                         OverlapMesg* olap,
                         int qt,
                         int* l,
                         FILE* fp);

#endif








