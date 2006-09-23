
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
#ifndef AS_CGW_EDGEDIAGNOSTICS_H
#define AS_CGW_EDGEDIAGNOSTICS_H

#include <math.h>

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_UTL_Hash.h"
#include "InputDataTypes_CGW.h"


// used for either contigs or scaffolds
typedef struct
{
  CDS_CID_t  keyID;    // contigID for contigs, scaffoldID for scaffolds
  CDS_CID_t  secondID; // scaffoldID for contigs, else unused
  FragOrient orient;
} OrientHolder;

VA_DEF(OrientHolder);


typedef struct
{
  VA_TYPE(OrientHolder) * array;
  HashTable_AS * ht;
} OrientProcessor;

typedef struct
{
  OrientProcessor * contigs;
  OrientProcessor * scaffolds;
} ContigOrientChecker;


void GetFragment5pAndOrientationInChunk(ScaffoldGraphT * graph,
                                        CIFragT * frag,
                                        LengthT * offset5p,
                                        FragOrient * orient,
                                        ChunkInstanceT * chunk);
void ComputeFrag5pToChunkEndFromOffset(LengthT * fragOffset5p,
                                       ChunkOrientType whichEnd,
                                       LengthT * distFromEnd,
                                       ChunkInstanceT * chunk);
void ComputeFrag5pToChunkEnd(ScaffoldGraphT * graph,
                             CIFragT * frag,
                             ChunkOrientType whichEnd,
                             LengthT * distFromEnd,
                             ChunkInstanceT * chunk);
void ComputeFragToChunkEndForEdge(ScaffoldGraphT * graph,
                                  CIFragT * frag,
                                  ChunkOrientType * endFromFrag,
                                  LengthT * distFromEnd,
                                  ChunkInstanceT * chunk);
void PopulateChunkEdgeBasics(ScaffoldGraphT * graph,
                             CIFragT * fragA,
                             ChunkInstanceT * chunkA,
                             CIFragT * fragB,
                             ChunkInstanceT * chunkB,
                             DistT * dist,
                             EdgeCGW_T * edge);
void ValidateAllContigEdges(ScaffoldGraphT * graph, int fixBadOnes);


void DestroyContigOrientChecker(ContigOrientChecker * coc);
ContigOrientChecker * CreateContigOrientChecker(void);
void ResetContigOrientChecker(ContigOrientChecker * coc);
int AddScaffoldToContigOrientChecker(ScaffoldGraphT * graph,
                                     CIScaffoldT * scaffold,
                                     ContigOrientChecker * coc);
int CompareNewOrientations(ScaffoldGraphT * graph,
                           CIScaffoldT * scaffold,
                           ContigOrientChecker * coc);

int AddAllScaffoldsToContigOrientChecker(ScaffoldGraphT * graph,
                                         ContigOrientChecker * coc);
int CheckAllContigOrientationsInAllScaffolds(ScaffoldGraphT * graph,
                                             ContigOrientChecker * coc,
                                             int populateHashTables);

void PrintScaffoldConnectivity(ScaffoldGraphT * graph,
                               CIScaffoldT * scaffold,
                               CDS_CID_t otherScaffoldID);
void PrintScaffoldConnectivityForIID(ScaffoldGraphT * graph,
                                     CDS_CID_t scaffoldID);
void PrintConnectivityBetweenScaffolds(ScaffoldGraphT * graph,
                                       CDS_CID_t scaffoldID1,
                                       CDS_CID_t scaffoldID2);
void DetectRepetitiveContigs(ScaffoldGraphT * graph);
#endif
