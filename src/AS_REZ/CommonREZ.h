
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
CVS_ID:  $Id: CommonREZ.h,v 1.2 2004-09-23 20:25:27 mcschatz Exp $
 *********************************************************************/

#ifndef COMMONREZ_H
#define COMMONREZ_H

#include "DataTypesREZ.h"
#include "dpc_CNS.h"


#define  MAKE_CAM_FILE         0
  // if  1  will create celamy file  rez.cam
#define  SIMULATED_DATA        CHECK_CELSIM_COORDS
  // if  1  assume celsim coordinates and info are available


#define  NUM_COLOURS           14

#define  UNIQUE_COLOUR         1
#define  INSERTED_COLOUR       2
#define  NO_CONNECT_COLOUR     3
#define  CONNECT_COLOUR        4
#define  CONSISTENT_COLOUR     5
#define  PLACED_COLOUR         6
#define  MISPLACED_COLOUR      7
#define  REJECT_COLOUR         8
#define  SCAFFOLD_COLOUR       9
#define  LO_COVERSTAT_COLOUR  10
#define  STONE_COLOUR         11
#define  RU_RR_COLOUR         12
#define  UR_COLOUR            13

extern FILE  * Cam_File;


// Flags for Find_Olap_Path edge mask

#define  USE_TANDEM_OLAPS          0
#define  SKIP_TANDEM_OLAPS         1
#define  SKIP_CONTAINMENT_OLAPS    2


char *  CGB_Type_As_String
    (unsigned int t);

float  CIEdge_Quality
    (CIEdgeT *);

int  Find_Olap_Path
    (ChunkInstanceT * from, int chunk_end, ChunkInstanceT * to,
     int num_targets, Target_Info_t target [], double bound,
     int * first, int * max_hits, int * max_first, LengthT * to_position,
     unsigned int edge_mask);

void  Free_Fill_Array
    (Scaffold_Fill_t * fill_chunks);

void  Force_Increasing_Variances
    (void);

int IsSurrogate(ChunkInstanceT * chunk);

int  Is_Unique
    (ChunkInstanceT *);

#if  0
Overlap *  OverlapSequences
    (char * seq1, char * seq2, ChunkOrientationType orientation, 
     int min_ahang, int max_ahang, double erate, double thresh, int minlen,
     CompareOptions what,
     Overlap *(*COMPARE_FUNC)(COMPARE_ARGS));
#endif

void  Print_Fill_Info
    (FILE *, Scaffold_Fill_t *);

Scaffold_Fill_t *  Scan_Gaps
    (void);

#endif
