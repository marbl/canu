
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

static const char *rcsid_INTERLEAVED_MERGING_H = "$Id: InterleavedMerging.h,v 1.14 2011-12-29 09:26:03 brianwalenz Exp $";

#include "AS_global.h"
#include "AS_UTL_Var.h"

#include "InputDataTypes_CGW.h"
#include "ScaffoldGraph_CGW.h"
#include "GraphCGW_T.h"

#include "CA_ALN_scafcomp.h"

#include "Instrument_CGW.h"

VA_DEF(Segment)
VA_DEF(Local_Overlap)
VA_DEF(Scaffold_Gap)
VA_DEF(Scaffold_Tig)
VA_DEF(PtrT)

typedef struct {
  int index;
  CDS_CID_t id;
  double length;
  double minCoord;
  double maxCoord;
  SequenceOrient orient;
} ContigElement;

VA_DEF(ContigElement)
VA_DEF(MateInstrumenterP)

typedef struct {
  VA_TYPE(Scaffold_Gap) * gapPool;
  VA_TYPE(Scaffold_Tig) * tigPool;
} ScaffoldPools;

typedef struct {
  int32 bandBeg;             // max AHang - only for scaffold A
  int32 bandEnd;             // min AHang - only for scaffold A
  Scaffold * scaffold;
  ScaffoldPools * pools;
  VA_TYPE(ContigElement) * contigs;
  VA_TYPE(ContigElement) * edgeContigs;
  SequenceOrient orient;
} ScaffoldStuff;

typedef struct {
  Segment * segmentList;
  int numSegs;
  int varWin;
  ScaffoldStuff * scaffoldA;
  ScaffoldStuff * scaffoldB;
  int best;
  Segment * lastSegment;
  ScaffoldInstrumenter * scaffInst;       // used in instrumenting
} ScaffoldAlignmentInterface;

typedef struct {
  ScaffoldAlignmentInterface * sai;
  int contigNow;
  int checkForTinyScaffolds;
  int checkAbutting;
  double minSatisfied;     // 0 to 1 - .98 or so?
  double maxDelta;         // 0 to 1 - .01 or so?
  VA_TYPE(MateInstrumenterP) *MIs;
  ChunkOverlapperT * badSEdges;
} InterleavingSpec;



ScaffoldAlignmentInterface *
CreateScaffoldAlignmentInterface(void);


int
PopulateScaffoldAlignmentInterface(CIScaffoldT * scaffoldA,
                                   CIScaffoldT * scaffoldB,
                                   SEdgeT * sEdge,
                                   ScaffoldAlignmentInterface * sai);


SEdgeT *
MakeScaffoldAlignmentAdjustments(CIScaffoldT * scaffoldA,
                                 CIScaffoldT * scaffoldB,
                                 SEdgeT * sEdge,
                                 ScaffoldAlignmentInterface * sai);


void
DeleteScaffoldAlignmentInterface(ScaffoldAlignmentInterface * sai);


static
int
GetNumSegmentsInList(Segment *segmentList) {
  int numSegments = 0;

  while (segmentList != NULL) {
    segmentList = segmentList->next;
    numSegments++;
  }

  return(numSegments);
}


#endif // INTERLEAVED_MERGING_H
