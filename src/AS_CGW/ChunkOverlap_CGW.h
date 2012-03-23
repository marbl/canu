
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

#ifndef CHUNKOVERLAP_CGW_H
#define CHUNKOVERLAP_CGW_H

static const char *rcsid_CHUNKOVERLAP_CGW_H = "$Id: ChunkOverlap_CGW.h,v 1.12 2012-03-23 06:45:50 brianwalenz Exp $";

#include "AS_global.h"
#include "AS_UTL_Hash.h"
#include "AS_ALN_aligners.h"

#include "GraphCGW_T.h"  //  For GraphCGW_T, EdgeCGW_T



//  This structure comprises the 'symbol' that is the key for database
//  lookups
typedef struct {
  CDS_CID_t  cidA;
  CDS_CID_t  cidB;
  PairOrient orientation;
} ChunkOverlapSpecT;



//  This is the value stored in the symbol table.  Note that it's
//  first field is the key for the symbol table.  Thus we kill two
//  birds with one stone.  These structures are allocated from a heap
//
typedef struct {
  ChunkOverlapSpecT spec;

  int32   minOverlap;
  int32   maxOverlap;
  int32   cgbMinOverlap;
  int32   cgbMaxOverlap;
  double  errorRate;

  // This is what we found
  uint32       computed:1;
  uint32       fromCGB:1;
  uint32       AContainsB:1;
  uint32       BContainsA:1;
  uint32       suspicious:1;
  uint32       unused:26;

  int32   overlap;  // The overlaplength if there is an overlap, 0 otherwise
  int32   ahg;
  int32   bhg;
  double  quality;
  int32   min_offset;
  int32   max_offset;
} ChunkOverlapCheckT;



typedef struct {
  HashTable_AS *hashTable;
  Heap_AS      *ChunkOverlaps;  //  Heap of ChunkOverlapCheckT
} ChunkOverlapperT;




ChunkOverlapperT *
CreateChunkOverlapper(void);

void
DestroyChunkOverlapper(ChunkOverlapperT *chunkOverlapper);

int
InsertChunkOverlap(ChunkOverlapperT *chunkOverlapper,
                   ChunkOverlapCheckT *olap);

void
SaveChunkOverlapperToStream(ChunkOverlapperT *chunkOverlapper, FILE *stream);

ChunkOverlapperT *
LoadChunkOverlapperFromStream(FILE *stream);

int
InitCanonicalOverlapSpec(CDS_CID_t cidA, CDS_CID_t cidB,
                         PairOrient orientation,
                         ChunkOverlapSpecT *spec);

void
CreateChunkOverlapFromEdge(GraphCGW_T *graph,
                           EdgeCGW_T *edge);

void
FillChunkOverlapWithEdge(EdgeCGW_T *edge, ChunkOverlapCheckT *olap);

ChunkOverlapCheckT *
LookupCanonicalOverlap(ChunkOverlapperT *chunkOverlapper,
                       ChunkOverlapSpecT *spec);

int
LookupOverlap(GraphCGW_T *graph,
              CDS_CID_t cidA,
              CDS_CID_t cidB,
              PairOrient orientation,
              ChunkOverlapCheckT *olap);

CDS_CID_t
InsertComputedOverlapEdge(GraphCGW_T *graph,
                          ChunkOverlapCheckT *olap);

void
CollectChunkOverlap(GraphCGW_T *graph,
                    CDS_CID_t cidA, CDS_CID_t cidB,
                    PairOrient orientation,
                    double   meanOverlap, double   deltaOverlap,
                    double   quality, int bayesian,
                    int fromCGB,
                    int verbose);

ALNoverlap *
OverlapSequences(char *seq1, char *seq2,
                 PairOrient orientation,
                 int32 min_ahang, int32 max_ahang,
                 double erate, double thresh, int32 minlen,
                 uint32 tryLocal = FALSE);

ChunkOverlapCheckT
OverlapChunks(GraphCGW_T *graph,
              CDS_CID_t cidA, CDS_CID_t cidB,
              PairOrient orientation,
              int32 minOverlap,
              int32 maxOverlap,
              double errorRate,
              int insertGraphEdges);

ALNoverlap *
OverlapContigs(NodeCGW_T *contig1, NodeCGW_T *contig2,
               PairOrient *overlapOrientation,
               int32 minAhang, int32 maxAhang,
               int computeAhang,
               uint32 tryLocal = FALSE,
               uint32 tryRev   = FALSE);


void
ComputeOverlaps(GraphCGW_T *graph, int addEdgeMates,
                int recomputeCGBOverlaps);


#endif  //  CHUNKOVERLAP_CGW_H
