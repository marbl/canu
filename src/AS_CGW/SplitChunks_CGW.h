
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
/* $Id: SplitChunks_CGW.h,v 1.1.1.1 2004-04-14 13:51:05 catmandew Exp $ */

#ifndef SPLITCHUNKS_H
#define SPLITCHUNKS_H

#include "AS_UTL_Var.h"
#include "ScaffoldGraph_CGW.h"

// number of contigs off A/B end with links before/after chimeric interval
#define CHIMERA_SET_THRESHOLD 1

typedef struct
{
  CDS_CID_t id;
  SeqInterval interval;
} SplitInterval;

VA_DEF(SplitInterval)
  
int SplitInputUnitigs(ScaffoldGraphT * graph);
VA_TYPE(SplitInterval) * DetectChimericChunksInGraph(ScaffoldGraphT * graph);

#endif
