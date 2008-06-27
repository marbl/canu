
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
#ifndef INTERLEAVED_MERGING_H
#define INTERLEAVED_MERGING_H

#include "AS_global.h"
#include "AS_UTL_Var.h"

#include "InputDataTypes_CGW.h"
#include "ScaffoldGraph_CGW.h"
#include "GraphCGW_T.h"

#include "CA_ALN_scafcomp.h"

#include "Instrument_CGW.h"

VA_DEF(Segment);
VA_DEF(Local_Overlap);
VA_DEF(Scaffold_Gap);
VA_DEF(Scaffold_Tig);

typedef struct
{
  int index;
  CDS_CID_t id;
  double length;
  double minCoord;
  double maxCoord;
  NodeOrient orient;
} ContigElement;

VA_DEF(ContigElement);

typedef struct
{
  VA_TYPE(Scaffold_Gap) * gapPool;
  VA_TYPE(Scaffold_Tig) * tigPool;
} ScaffoldPools;

typedef struct
{
  CDS_COORD_t bandBeg;             // max AHang - only for scaffold A
  CDS_COORD_t bandEnd;             // min AHang - only for scaffold A
  Scaffold * scaffold;
  ScaffoldPools * pools;
  VA_TYPE(ContigElement) * contigs;
  VA_TYPE(ContigElement) * edgeContigs;
  NodeOrient orient;
} ScaffoldStuff;

typedef struct
{
  Segment * segmentList;
  int numSegs;
  int varWin;
  ScaffoldStuff * scaffoldA;
  ScaffoldStuff * scaffoldB;
  int best;
  Segment * lastSegment;
  ScaffoldInstrumenter * scaffInst;       // used in instrumenting
} ScaffoldAlignmentInterface;


VA_DEF(PtrT);

typedef struct
{
  ScaffoldAlignmentInterface * sai;
  int contigNow;
  int checkForTinyScaffolds;
  int checkAbutting;
  float minSatisfied;     // 0 to 1 - .98 or so?
  float maxDelta;         // 0 to 1 - .01 or so?
  VA_TYPE(PtrT) * MIs;
  ChunkOverlapperT * badSEdges;
} InterleavingSpec;


void DeleteScaffoldAlignmentInterface(ScaffoldAlignmentInterface * sai);
ScaffoldAlignmentInterface * CreateScaffoldAlignmentInterface(void);
int PopulateScaffoldAlignmentInterface(CIScaffoldT * scaffoldA,
                                       CIScaffoldT * scaffoldB,
                                       SEdgeT * sEdge,
                                       ScaffoldAlignmentInterface * sai);
SEdgeT * MakeScaffoldAlignmentAdjustments(CIScaffoldT * scaffoldA,
                                          CIScaffoldT * scaffoldB,
                                          SEdgeT * sEdge,
                                          ScaffoldAlignmentInterface * sai);
Segment* DuplicateSegmentList(Segment * segmentList);
int GetNumSegmentsInList(Segment * segmentList);


#endif // INTERLEAVED_MERGING_H
