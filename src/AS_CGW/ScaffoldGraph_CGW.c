
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
static char *rcsid = "$Id: ScaffoldGraph_CGW.c,v 1.43 2009-09-12 22:35:57 brianwalenz Exp $";

//#define DEBUG 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_UTL_fileIO.h"
#include "UtilsREZ.h"
#include "AS_CGW_dataTypes.h"
#include "Globals_CGW.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"
#include "RepeatRez.h"
#include "CommonREZ.h"
#include "Stats_CGW.h"

ScaffoldGraphT *ScaffoldGraph = NULL;
tSequenceDB    *SequenceDB    = NULL;

void ClearChunkInstance(ChunkInstanceT *ci){
  memset(ci, 0, sizeof(ChunkInstanceT));
  ci->id = ci->scaffoldID = NULLINDEX;
  ci->info.CI.contigID = NULLINDEX;
  ci->flags.all = 0;
}

void
LoadScaffoldGraphFromCheckpoint(char   *name,
                                int32   checkPointNum,
                                int     readWrite){
  char ckpfile[FILENAME_MAX];
  char tmgfile[FILENAME_MAX];

  sprintf(ckpfile, "%s.ckp.%d", name, checkPointNum);
  sprintf(tmgfile, "%s.timing", name);

  time_t t = time(0);
  fprintf(stderr, "====> Reading %s at %s", ckpfile, ctime(&t));

  errno = 0;
  FILE *F = fopen(tmgfile, "a");
  if (errno == 0) {
    fprintf(F, "====> Reading %s at %s", ckpfile, ctime(&t));
    fclose(F);
  }

  errno = 0;
  F = fopen(ckpfile, "r");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for reading checkpoint: %s\n", ckpfile, strerror(errno)), exit(1);

  ScaffoldGraph = (ScaffoldGraphT *)safe_calloc(1, sizeof(ScaffoldGraphT));
  int    i, status;

  status = AS_UTL_safeRead(F, ScaffoldGraph->name, "LoadScaffoldGraphFromCheckpoint", sizeof(char), 256);
  assert(status == 256);

  ScaffoldGraph->iidToFragIndex = CreateFromFileVA_InfoByIID(F);
  ScaffoldGraph->CIFrags        = CreateFromFileVA_CIFragT(F);
  ScaffoldGraph->Dists          = CreateFromFileVA_DistT(F);

  ScaffoldGraph->CIGraph       = LoadGraphCGWFromStream(F);
  ScaffoldGraph->ContigGraph   = LoadGraphCGWFromStream(F);
  ScaffoldGraph->ScaffoldGraph = LoadGraphCGWFromStream(F);

  CheckGraph(ScaffoldGraph->CIGraph);
  CheckGraph(ScaffoldGraph->ContigGraph);
  CheckGraph(ScaffoldGraph->ScaffoldGraph);

  //  Erase the bogus pointer in dists
  //
  for (i=0; i<GetNumDistTs(ScaffoldGraph->Dists); i++)
    GetDistT(ScaffoldGraph->Dists, i)->histogram = NULL;

  // Temporary
  ScaffoldGraph->ChunkInstances = ScaffoldGraph->CIGraph->nodes;
  ScaffoldGraph->Contigs        = ScaffoldGraph->ContigGraph->nodes;
  ScaffoldGraph->CIScaffolds    = ScaffoldGraph->ScaffoldGraph->nodes;
  ScaffoldGraph->CIEdges        = ScaffoldGraph->CIGraph->edges;
  ScaffoldGraph->ContigEdges    = ScaffoldGraph->ContigGraph->edges;
  ScaffoldGraph->SEdges         = ScaffoldGraph->ScaffoldGraph->edges;
  ScaffoldGraph->overlapper     = ScaffoldGraph->CIGraph->overlapper;

  status  = AS_UTL_safeRead(F, &ScaffoldGraph->doRezOnContigs,            "LoadScaffoldGraphFromCheckpoint", sizeof(int32), 1);
  status += AS_UTL_safeRead(F, &ScaffoldGraph->checkPointIteration,       "LoadScaffoldGraphFromCheckpoint", sizeof(int32), 1);
  status += AS_UTL_safeRead(F, &ScaffoldGraph->numContigs,                "LoadScaffoldGraphFromCheckpoint", sizeof(int32), 1);
  status += AS_UTL_safeRead(F, &ScaffoldGraph->numDiscriminatorUniqueCIs, "LoadScaffoldGraphFromCheckpoint", sizeof(int32), 1);
  status += AS_UTL_safeRead(F, &ScaffoldGraph->numOriginalCIs,            "LoadScaffoldGraphFromCheckpoint", sizeof(int32), 1);
  status += AS_UTL_safeRead(F, &ScaffoldGraph->numLiveCIs,                "LoadScaffoldGraphFromCheckpoint", sizeof(int32), 1);
  status += AS_UTL_safeRead(F, &ScaffoldGraph->numLiveScaffolds,          "LoadScaffoldGraphFromCheckpoint", sizeof(int32), 1);
  assert(status == 7);

  if(ScaffoldGraph->doRezOnContigs){
    ScaffoldGraph->RezGraph = ScaffoldGraph->ContigGraph;
  }else{
    ScaffoldGraph->RezGraph = ScaffoldGraph->CIGraph;
  }

  ReportMemorySize(ScaffoldGraph,stderr);

  //  Check that the iteration we were told to load is what we loaded.
  //  Generally not a good thing if these disagree (SeqStore will
  //  probably be wrong).
  //
  if (ScaffoldGraph->checkPointIteration != (checkPointNum + 1)){
    fprintf(stderr,"ERROR:  Checkpoint claims to be at iteration %d, but we wanted to load iteration %d!\n",
	    ScaffoldGraph->checkPointIteration, checkPointNum);
    exit(1);
  }

  //  Open the seqStore
  sprintf(ckpfile, "%s.SeqStore", name);
  ScaffoldGraph->sequenceDB = SequenceDB = openSequenceDB(ckpfile, readWrite, checkPointNum);

  //  Open the gkpStore
  ScaffoldGraph->gkpStore = new gkStore(GlobalData->gkpStoreName, FALSE, FALSE);

  //  Check (and cleanup?) scaffolds
  SetCIScaffoldTLengths(ScaffoldGraph, TRUE);
  CheckCIScaffoldTs(ScaffoldGraph);

  fclose(F);
}



void
CheckpointScaffoldGraph(const char *logicalname, const char *location) {
  char ckpfile[FILENAME_MAX];
  char tmgfile[FILENAME_MAX];

  sprintf(ckpfile, "%s.ckp.%d", GlobalData->outputPrefix, ScaffoldGraph->checkPointIteration++);
  sprintf(tmgfile, "%s.timing", GlobalData->outputPrefix);

  errno = 0;
  FILE *F = fopen(ckpfile, "w");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for writing checkpoint: %s\n", ckpfile, strerror(errno)), exit(1);

  AS_UTL_safeWrite(F, ScaffoldGraph->name, "CheckpointScaffoldGraph", sizeof(char), 256);

  CopyToFileVA_InfoByIID(ScaffoldGraph->iidToFragIndex, F);
  CopyToFileVA_CIFragT(ScaffoldGraph->CIFrags, F);
  CopyToFileVA_DistT(ScaffoldGraph->Dists, F);

  SaveGraphCGWToStream(ScaffoldGraph->CIGraph,F);
  SaveGraphCGWToStream(ScaffoldGraph->ContigGraph,F);
  SaveGraphCGWToStream(ScaffoldGraph->ScaffoldGraph,F);

  saveSequenceDB(ScaffoldGraph->sequenceDB);

  AS_UTL_safeWrite(F, &ScaffoldGraph->doRezOnContigs,            "CheckpointScaffoldGraph", sizeof(int32), 1);
  AS_UTL_safeWrite(F, &ScaffoldGraph->checkPointIteration,       "CheckpointScaffoldGraph", sizeof(int32), 1);
  AS_UTL_safeWrite(F, &ScaffoldGraph->numContigs,                "CheckpointScaffoldGraph", sizeof(int32), 1);
  AS_UTL_safeWrite(F, &ScaffoldGraph->numDiscriminatorUniqueCIs, "CheckpointScaffoldGraph", sizeof(int32), 1);
  AS_UTL_safeWrite(F, &ScaffoldGraph->numOriginalCIs,            "CheckpointScaffoldGraph", sizeof(int32), 1);
  AS_UTL_safeWrite(F, &ScaffoldGraph->numLiveCIs,                "CheckpointScaffoldGraph", sizeof(int32), 1);
  AS_UTL_safeWrite(F, &ScaffoldGraph->numLiveScaffolds,          "CheckpointScaffoldGraph", sizeof(int32), 1);

  fclose(F);

  time_t t = time(0);
  fprintf(stderr, "====> Writing %s (logical %s) %s at %s", ckpfile, logicalname, location, ctime(&t));

  errno = 0;
  F = fopen(tmgfile, "a");
  if (errno == 0) {
    fprintf(F, "====> Writing %s (logical %s) %s at %s", ckpfile, logicalname, location, ctime(&t));
    fclose(F);
  }
}





/***************************************************************************/
// TEMPORARY HACK
//
void InsertRepeatCIsInScaffolds(ScaffoldGraphT *sgraph){
  CDS_CID_t cid;
  int scaffolded = 0;
  LengthT NullLength = {0,0.0};
  CIScaffoldT CIScaffold;
  GraphNodeIterator nodes = {0};
  NodeCGW_T *CI;

  memset(&CIScaffold, 0, sizeof(CIScaffoldT));

  CIScaffold.info.Scaffold.AEndCI = NULLINDEX;
  CIScaffold.info.Scaffold.BEndCI = NULLINDEX;
  CIScaffold.info.Scaffold.numElements = 0;
  CIScaffold.edgeHead = NULLINDEX;
  CIScaffold.bpLength = NullLength;
  CIScaffold.id = NULLINDEX;
  InitializeScaffold(&CIScaffold, OUTPUT_SCAFFOLD);
  CIScaffold.flags.bits.isDead = FALSE;
  CIScaffold.numEssentialA = CIScaffold.numEssentialB = 0;
  CIScaffold.essentialEdgeA = CIScaffold.essentialEdgeB = NULLINDEX;

  InitGraphNodeIterator(&nodes, sgraph->RezGraph, GRAPH_NODE_DEFAULT);
  while(NULL != (CI = NextGraphNodeIterator(&nodes))){
    ChunkInstanceType type = CI->type;
    cid = CI->id;

    // Skip scaffolded CIs
    if(CI->scaffoldID != NULLINDEX)
      continue;

    scaffolded++;
    CI->scaffoldID = CIScaffold.id = GetNumGraphNodes(sgraph->ScaffoldGraph);

    AppendGraphNode(sgraph->ScaffoldGraph, &CIScaffold);

    InsertCIInScaffold(sgraph, cid, CIScaffold.id, NullLength,
                       CI->bpLength,  FALSE, FALSE);
    CI->flags.bits.isUnique = FALSE;
    CI->type = type; // restore type for graphics
  }

  fprintf(stderr,
          "* Inserting %d Unscaffolded CIs into newly created scaffolds\n",
	  scaffolded);

}

/****************************************************************************/

ScaffoldGraphT *CreateScaffoldGraph(int rezOnContigs, char *name) {
  char buffer[2000];
  ScaffoldGraphT *sgraph = (ScaffoldGraphT *)safe_calloc(1,sizeof(ScaffoldGraphT));

  strcpy(sgraph->name, name);

  sprintf(buffer,"%s.SeqStore", name);
  sgraph->sequenceDB    = createSequenceDB(buffer);

  sgraph->CIGraph       = CreateGraphCGW(CI_GRAPH, 16 * 1024, 16 * 1024);
  sgraph->ContigGraph   = CreateGraphCGW(CONTIG_GRAPH, 1, 1);
  sgraph->ScaffoldGraph = CreateGraphCGW(SCAFFOLD_GRAPH, NULLINDEX, NULLINDEX);

  sgraph->gkpStore = new gkStore(GlobalData->gkpStoreName, FALSE, FALSE);

  int numFrags = sgraph->gkpStore->gkStore_getNumFragments();
  int numDists = sgraph->gkpStore->gkStore_getNumLibraries();

  sgraph->iidToFragIndex = CreateVA_InfoByIID(numFrags);
  sgraph->Dists          = CreateVA_DistT(numDists);
  sgraph->CIFrags        = CreateVA_CIFragT(numFrags);

  if(rezOnContigs){
    // switch this after doing CI overlaps
    sgraph->RezGraph = sgraph->CIGraph;
    //    sgraph->overlapper = sgraph->ContigGraph->overlapper;
    sgraph->overlapper = sgraph->CIGraph->overlapper;
  }else{
    sgraph->RezGraph = sgraph->CIGraph;
    sgraph->overlapper = sgraph->CIGraph->overlapper;
  }

  sgraph->doRezOnContigs = rezOnContigs;

  // Temporary
  sgraph->ChunkInstances = sgraph->CIGraph->nodes;
  sgraph->CIEdges        = sgraph->CIGraph->edges;

  sgraph->Contigs        = sgraph->ContigGraph->nodes;
  sgraph->ContigEdges    = sgraph->ContigGraph->edges;

  sgraph->CIScaffolds    = sgraph->ScaffoldGraph->nodes;
  sgraph->SEdges         = sgraph->ScaffoldGraph->edges;

  return sgraph;
}

/* Destructor */
void DestroyScaffoldGraph(ScaffoldGraphT *sgraph){
  int i;

  deleteSequenceDB(sgraph->sequenceDB);

  DeleteGraphCGW(sgraph->CIGraph);
  DeleteGraphCGW(sgraph->ContigGraph);
  DeleteGraphCGW(sgraph->ScaffoldGraph);

  delete sgraph->gkpStore;

  DeleteVA_CIFragT(sgraph->CIFrags);

  DeleteVA_DistT(sgraph->Dists);
  DeleteVA_InfoByIID(sgraph->iidToFragIndex);

  safe_free(sgraph);
}


void ReportMemorySize(ScaffoldGraphT *graph, FILE *stream){

  size_t totalMemorySize = 0;

  fprintf(stream,"* ScaffoldGraph Memory Report:\n");
  fprintf(stream,"*\t     \tnumElements\tAllocatedElements\tAllocatedSize\n");
  totalMemorySize += ReportMemorySize_VA(graph->Dists, "Dists", stream);
  totalMemorySize += ReportMemorySize_VA(graph->CIFrags, "CIFrags", stream);
  totalMemorySize += ReportMemorySize_VA(graph->iidToFragIndex, "iidToFrag", stream);
  totalMemorySize += ReportMemorySizeGraphCGW(graph->CIGraph, stream);
  totalMemorySize += ReportMemorySizeGraphCGW(graph->ContigGraph, stream);
  totalMemorySize += ReportMemorySizeGraphCGW(graph->ScaffoldGraph, stream);
  fprintf(stream,"*\tTotalMemorySize = " F_SIZE_T "\n", totalMemorySize);
}



void ScaffoldSanity(CIScaffoldT *scaffold, ScaffoldGraphT *graph){
  int numElements = 0;
  CIScaffoldTIterator CIs;
  ChunkInstanceT *CI;
  double scaffoldMinPos = (double)INT32_MAX;
  double scaffoldMaxPos = (double)INT32_MIN;
  double scratch;
  int    scaffoldInsane = 0;

  assert(scaffold->flags.bits.isScaffold);
  if(scaffold->type != REAL_SCAFFOLD)
    return;

  InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &CIs);
  while(NULL != (CI = NextCIScaffoldTIterator(&CIs))){
    scratch = MIN(CI->offsetAEnd.mean, CI->offsetBEnd.mean);
    scaffoldMinPos = MIN(scaffoldMinPos, scratch);
    scratch = MAX(CI->offsetAEnd.mean, CI->offsetBEnd.mean);
    scaffoldMaxPos = MAX(scaffoldMaxPos, scratch);
  }

  if(scaffold->bpLength.mean < 0 ||scaffold->bpLength.variance < 0 ){
    fprintf(stderr,
            "*!!! Sanity  scaffold " F_CID " length (%g,%g) screwed up!\n",
            scaffold->id,
            scaffold->bpLength.mean, scaffold->bpLength.variance);
  }
  if(scaffold->info.Scaffold.AEndCI != NULLINDEX &&
     scaffold->bpLength.mean > 0 &&
     abs(scaffold->bpLength.mean - (scaffoldMaxPos - scaffoldMinPos)) > 100.0){
    fprintf(stderr,
            "*!!! Sanity  scaffold " F_CID " length %g not equal to (max - min) %g\n",
            scaffold->id,
            scaffold->bpLength.mean, (scaffoldMaxPos - scaffoldMinPos));
#ifdef STRICT_SCAFFOLD_CHECKING
    scaffoldInsane = 1;
#endif
    scaffold->bpLength.mean = scaffoldMaxPos - scaffoldMinPos;
  }
  InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &CIs);
  while(NULL != (CI = NextCIScaffoldTIterator(&CIs))){

    if(CI->scaffoldID == scaffold->id &&
       CI->flags.bits.isUnique &&
       !CI->flags.bits.isDead){
      numElements++;
    }else{
      fprintf(stderr,
              "* Scaffold Sanity: node " F_CID " is screwed: %s%s%sin scaffold " F_CID "\n",
              CI->id,
              (CI->scaffoldID == scaffold->id?"":"scaffoldID is wrong "),
              (!CI->flags.bits.isDead?"":"Dead "),
              (CI->flags.bits.isUnique?"":"Not Unique "),
              scaffold->id);
#ifdef STRICT_SCAFFOLD_CHECKING
      scaffoldInsane = 1;
#endif
    }
  }

  if(numElements != scaffold->info.Scaffold.numElements){
    fprintf(stderr,
            "* numElements = %d scaffold says it has %d elements\n",
            numElements, scaffold->info.Scaffold.numElements);
    scaffoldInsane = 1;
  }

  if (scaffoldInsane) {
    DumpCIScaffold(stderr,graph, scaffold, FALSE);
    assert(!scaffoldInsane);
  }

  // CheckScaffoldOrder(scaffold, graph);
}




// CheckScaffoldOrder checks whether all the ahangs in a scaffold are
// positive (ie, that the contigs are ordered by their distance from the
// A end of the scaffold)
void CheckScaffoldOrder(CIScaffoldT *scaffold, ScaffoldGraphT *graph)
{
  double currentMinPos = (double)INT32_MAX;
  CIScaffoldTIterator CIs;
  ChunkInstanceT *CI, *prevCI = NULL;

  assert(scaffold->flags.bits.isScaffold);
  if(scaffold->type != REAL_SCAFFOLD)
    return;

  CI = GetGraphNode(graph->RezGraph, scaffold->info.Scaffold.AEndCI);

  currentMinPos = MIN( CI->offsetAEnd.mean, CI->offsetBEnd.mean);
  InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &CIs);
  while(NULL != (CI = NextCIScaffoldTIterator(&CIs)))
    {
      if( MIN( CI->offsetAEnd.mean, CI->offsetBEnd.mean) < currentMinPos)
        {
          fprintf( stderr,
                   "CIs " F_CID " and " F_CID "are out of order\n",
                   CI->id, prevCI->id);
          fprintf( stderr,
                   "CI " F_CID ": AEndOffset.mean: %f, AEndOffset.mean: %f\n",
                   CI->id, CI->offsetAEnd.mean, CI->offsetBEnd.mean);
          fprintf( stderr,
                   "CI " F_CID ": AEndOffset.mean: %f, AEndOffset.mean: %f\n",
                   prevCI->id, prevCI->offsetAEnd.mean, prevCI->offsetBEnd.mean);
          DumpCIScaffold(stderr, graph, scaffold, FALSE);

          // allow for a base of rounding error, but fix it
          if( MIN( CI->offsetAEnd.mean, CI->offsetBEnd.mean) - currentMinPos < 1.0)
            {
              CI->offsetAEnd.mean += 1.0;
              CI->offsetBEnd.mean += 1.0;
              fprintf( stderr,
                       "shifted pos of CI " F_CID " to (%f, %f)\n",
                       CI->id, CI->offsetAEnd.mean, CI->offsetBEnd.mean);
            }
          else
            assert(0);
        }
      currentMinPos = MIN( CI->offsetAEnd.mean, CI->offsetBEnd.mean);
      prevCI = CI;
    }
}



#if 0
void DumpScaffoldGraph(ScaffoldGraphT *graph){
  static int version = 0;
  char fileName[FILENAME_MAX];
  FILE *dumpFile = NULL;

  sprintf(fileName,"%s.sg%c.%d",
          graph->name, (graph->doRezOnContigs?'C':'c'),version++);

  fprintf(stderr,"* Dumping graph file %s\n", fileName);

  dumpFile = fopen(fileName,"w");

  //  DumpGraph(ScaffoldGraph->CIGraph, dumpFile);
  DumpGraph(ScaffoldGraph->RezGraph, dumpFile);
  DumpGraph(ScaffoldGraph->ScaffoldGraph, dumpFile);

  fclose(dumpFile);
}
#endif


int GetCoverageStat(ChunkInstanceT *CI){
  // If this is a ChunkInstance, return its coverage stat
  // If this is an unscaffolded (singleton) contig, return its
  // lone ChunkInstance's coverage stat else, assert.

  if(CI->flags.bits.isCI)
    return CI->info.CI.coverageStat;

  if(CI->flags.bits.isContig){
    ChunkInstanceT *ci;

    if(CI->info.Contig.numCI == 1){
      ci = GetGraphNode(ScaffoldGraph->CIGraph, CI->info.Contig.AEndCI);
      AssertPtr(ci);

      return ci->info.CI.coverageStat;
    }else{ // if they have stuff in 'em, they're unique-equivalent
      return GlobalData->cgbUniqueCutoff;
    }
  }
  assert(0);
  return(0);
}

/* *********************************************************************** */
/* Add a fixed amount to the offsetAEnd and offsetBEnd starting from a given
   CI to the end of the Scaffold                                */
/* *********************************************************************** */
void AddDeltaToScaffoldOffsets(ScaffoldGraphT *graph,
                CDS_CID_t scaffoldIndex,
                CDS_CID_t indexOfCI,
                int aEndToBEnd,
                int verbose,
                LengthT delta){
AddDeltaToScaffoldOffsets(graph, scaffoldIndex, indexOfCI, aEndToBEnd, verbose, delta, 0);
}

void AddDeltaToScaffoldOffsets(ScaffoldGraphT *graph,
			       CDS_CID_t scaffoldIndex,
			       CDS_CID_t indexOfCI,
			       int aEndToBEnd,
			       int verbose,
			       LengthT delta,
                uint32 mark){
  CIScaffoldT *scaffold;
  CIScaffoldTIterator Nodes;
  NodeCGW_T *thisNode;

#ifdef DEBUG_DETAILED
  fprintf(stderr,
          "##### Adding delta (%g,%g) to scaffold " F_CID " at CI " F_CID " #######\n",
	  delta.mean, delta.variance, scaffoldIndex, indexOfCI);
#endif

  scaffold = GetGraphNode(graph->ScaffoldGraph, scaffoldIndex);

  InitCIScaffoldTIteratorFromCI(graph, scaffold, indexOfCI, aEndToBEnd,
				verbose, &Nodes);
  while((thisNode = NextCIScaffoldTIterator(&Nodes)) != NULL){
    thisNode->offsetAEnd.mean += delta.mean;
    thisNode->offsetAEnd.variance += delta.variance;
    thisNode->offsetBEnd.mean += delta.mean;
    thisNode->offsetBEnd.variance += delta.variance;

#ifdef TRY_UNDO_JIGGLE_POSITIONS    
    thisNode->flags.bits.isJiggled = mark;
    thisNode->offsetDelta.mean = delta.mean;
    thisNode->offsetDelta.variance = delta.variance;
#endif
  }

  scaffold->bpLength.mean += delta.mean;
  scaffold->bpLength.variance += delta.variance;
  return;
}




void
CheckCITypes(ScaffoldGraphT *sgraph){
  NodeCGW_T *CI;
  GraphNodeIterator nodes;

  //  An alternate iterator
  //  for(i = 0; i < GetNumGraphNodes(sgraph->RezGraph); i++)
  //    ChunkInstanceT *CI = GetGraphNode(sgraph->RezGraph,i);

  InitGraphNodeIterator(&nodes, sgraph->ContigGraph, GRAPH_NODE_DEFAULT);
  while((CI = NextGraphNodeIterator(&nodes)) != NULL){
    if(CI->flags.bits.isUnique ){
      if(CI->scaffoldID == NULLINDEX){
        CIEdgeT *nA = NULL, *nB = NULL;
        CDS_CID_t nidA = NULLINDEX, nidB = NULLINDEX;
        ChunkInstanceT *CIa = NULL, *CIb = NULL;

        fprintf(stderr,"* Dump SuspiciousCI " F_CID " of type %d**\n", CI->id, CI->type);
        DumpChunkInstance(stderr, ScaffoldGraph, CI, FALSE, FALSE, FALSE, FALSE);

        fprintf(stderr,"* numEssentialA:%d essentialA:" F_CID "  numEssentialB:%d essentialB:" F_CID "\n",
                CI->numEssentialA, CI->essentialEdgeA,
                CI->numEssentialB, CI->essentialEdgeB);

        fprintf(stderr,"* Essential Edges *\n");
        if(CI->essentialEdgeA != NULLINDEX){
          nA = GetCIEdgeT(ScaffoldGraph->CIEdges, CI->essentialEdgeA);
          nidA = (nA->idA == CI->id? nA->idB: nA->idA);
          CIa = GetChunkInstanceT (ScaffoldGraph->ChunkInstances, nidA);
          PrintCIEdgeT(stderr, ScaffoldGraph, " ", nA , nidA);
        }

        if(CI->essentialEdgeB != NULLINDEX){
          nB = GetCIEdgeT(ScaffoldGraph->CIEdges, CI->essentialEdgeB);
          nidB = (nB->idA == CI->id? nB->idB: nB->idA);
          CIb = GetChunkInstanceT (ScaffoldGraph->ChunkInstances, nidB);
          PrintCIEdgeT(stderr, ScaffoldGraph, " ", nB, nidB);
        }

        fprintf(stderr,"* Essential Neighbors *\n");
        if(CIa){
          fprintf(stderr,"* Chunk " F_CID " in Scaffold " F_CID " of type %d\n", CIa->id, CIa->scaffoldID, CIa->type);
          DumpChunkInstance(stderr,ScaffoldGraph, CIa, FALSE, FALSE, FALSE, FALSE);
        }

        if(CIb){
          fprintf(stderr,"* Chunk " F_CID " in Scaffold " F_CID " type %d\n", CIb->id, CIb->scaffoldID, CIb->type);
          DumpChunkInstance(stderr, ScaffoldGraph, CIb, FALSE, FALSE, FALSE, FALSE);
        }

	//assert(0);
      }
    }else{
      assert((CI->type != DISCRIMINATORUNIQUECHUNK_CGW) && (CI->type != UNIQUECHUNK_CGW));
      assert(CI->scaffoldID == NULLINDEX);
    }
  }
}



#define  MAX_REZ_ITERATIONS     3
#define  FILL_GAPS_THRESHHOLD   10
#define  DO_CONTAINED_ROCKS     1
#define  DO_CONTAINED_STONES    0

/* Repeat resolution */
int RepeatRez(int repeatRezLevel, char *name){
  int didSomething = FALSE;

  if(repeatRezLevel > 0) {
    int  iter = 0;
    int  normal_inserts, contained_inserts, contained_stones;

    //
    // repeat Fill_Gaps until we are not able to insert anymore
    //
    CheckEdgesAgainstOverlapper(ScaffoldGraph->RezGraph);
    // commented out next call to allow mouse_20010307 run to succeed
    CheckCITypes(ScaffoldGraph);
    CheckCIScaffoldTs(ScaffoldGraph);

    do
      {
        normal_inserts = Fill_Gaps (name, repeatRezLevel, iter);
        if  (normal_inserts > 0)
          {
            didSomething = TRUE;
            CheckCIScaffoldTs(ScaffoldGraph);

            TidyUpScaffolds (ScaffoldGraph);
            CheckEdgesAgainstOverlapper(ScaffoldGraph->RezGraph);

            CheckCIScaffoldTs(ScaffoldGraph);

            //GeneratePlacedContigGraphStats("rocks", iter);
            //GenerateLinkStats(ScaffoldGraph->ContigGraph,"rocks",iter);
            //GenerateScaffoldGraphStats("rocks",iter);
          }

        if  (DO_CONTAINED_ROCKS && iter == 0)
          {
            contained_inserts = Hurl_Contained_Rocks (name, repeatRezLevel, iter);
            if  (contained_inserts > 0)
              {
                didSomething = TRUE;
                TidyUpScaffolds (ScaffoldGraph);
                CheckEdgesAgainstOverlapper(ScaffoldGraph->RezGraph);
                //GeneratePlacedContigGraphStats("controcks", iter);
                //GenerateLinkStats(ScaffoldGraph->ContigGraph,"controcks",iter);
                //GenerateScaffoldGraphStats("controcks",iter);
              }
          }
        else
          contained_inserts = 0;

        if  (DO_CONTAINED_STONES && iter == 0)
          {
            contained_stones = Toss_Contained_Stones (name, repeatRezLevel, iter);
            if  (contained_stones > 0)
              {
                didSomething = TRUE;
                TidyUpScaffolds (ScaffoldGraph);
                CheckEdgesAgainstOverlapper(ScaffoldGraph->RezGraph);
                //GeneratePlacedContigGraphStats("contstones", iter);
                //GenerateLinkStats(ScaffoldGraph->ContigGraph,"contstones",iter);
                //GenerateScaffoldGraphStats("contstones",iter);
              }
          }
        else
          contained_stones = 0;

        iter++;
        if (iter >= MAX_REZ_ITERATIONS) {
          fprintf (stderr, "Maximum number of REZ iterations reached\n");
          break;
        }
      }  while  (normal_inserts + contained_inserts
                 + contained_stones > FILL_GAPS_THRESHHOLD);

    fprintf(stderr,"**** AFTER repeat rez ****\n");
  }
  return didSomething;
}





/***************************************************************************/
void RebuildScaffolds(ScaffoldGraphT *ScaffoldGraph,
                      int markShakyBifurcations){

#ifdef DEBUG_BUCIS
  BuildUniqueCIScaffolds(ScaffoldGraph, markShakyBifurcations,TRUE);
  fprintf(stderr,"** After BuildUniqueCIScaffolds **\n");
  DumpChunkInstances(stderr, ScaffoldGraph,
                     FALSE, TRUE, TRUE, FALSE);
  DumpCIScaffolds(stderr,ScaffoldGraph, TRUE);
#else
  BuildUniqueCIScaffolds(ScaffoldGraph, markShakyBifurcations, FALSE);
#endif
  CheckEdgesAgainstOverlapper(ScaffoldGraph->RezGraph);

#ifdef DEBUG_BUCIS
  LeastSquaresGapEstimates(ScaffoldGraph, TRUE, FALSE, TRUE, CHECK_CONNECTIVITY, TRUE);
#else
  LeastSquaresGapEstimates(ScaffoldGraph, TRUE, FALSE, TRUE, CHECK_CONNECTIVITY, FALSE);
#endif

  TidyUpScaffolds (ScaffoldGraph);

  return;
}


/***************************************************************************/
//  Used to be in middle of  RebuildScaffolds .
void  TidyUpScaffolds(ScaffoldGraphT *ScaffoldGraph)
{
  // Contig now!
  CleanupScaffolds(ScaffoldGraph, FALSE, NULLINDEX, FALSE);

  // We want tor recompute the contig coordinates
  LeastSquaresGapEstimates(ScaffoldGraph, TRUE, FALSE, TRUE,
			   CHECK_CONNECTIVITY, FALSE);

  if(GlobalData->debugLevel > 0)
    CheckAllContigFragments();

  CheckAllTrustedEdges(ScaffoldGraph);

  BuildSEdges(ScaffoldGraph, TRUE, FALSE);
  MergeAllGraphEdges(ScaffoldGraph->ScaffoldGraph, TRUE, FALSE);

  //clearCacheSequenceDB(ScaffoldGraph->sequenceDB);
}


