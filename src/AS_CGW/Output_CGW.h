
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

#ifndef OUTPUT_CGW_H
#define OUTPUT_CGW_H

static const char *rcsid_OUTPUT_CGW_H = "$Id: Output_CGW.h,v 1.9 2012-10-30 16:48:26 brianwalenz Exp $";

#include "Globals_CGW.h"
#include "GraphCGW_T.h"
#include "AS_MSG_pmesg.h"

void MarkContigEdges(void);

UnitigStatus  finalUnitigStatus(NodeCGW_T *unitig);
ContigStatus  finalContigStatus(NodeCGW_T *contig);

UnitigType    finalUnitigType(NodeCGW_T *unitig);

void OutputUnitigsFromMultiAligns(void);
void OutputContigsFromMultiAligns(int32 outputFragsPerPartition);

#endif

