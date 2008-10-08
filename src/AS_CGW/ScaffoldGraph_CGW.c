
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
static char *rcsid = "$Id: ScaffoldGraph_CGW.c,v 1.34 2008-10-08 22:02:55 brianwalenz Exp $";

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
#include "Checkpoints_CGW.h"

ScaffoldGraphT *ScaffoldGraph = NULL;
tSequenceDB *SequenceDB = NULL;

void ClearChunkInstance(ChunkInstanceT *ci){
  memset(ci, 0, sizeof(ChunkInstanceT));
  ci->id = ci->scaffoldID = NULLINDEX;
  ci->info.CI.contigID = NULLINDEX;
  ci->flags.all = 0;
}

ScaffoldGraphT *LoadScaffoldGraphFromCheckpoint( char *name,
                                                 int32 checkPointNum,
                                                 int readWrite){
  char buffer[1024];
  FILE *inStream;
  ScaffoldGraphT *graph;

  sprintf(buffer,"%s.SeqStore", name);
  SequenceDB = openSequenceDB(buffer, readWrite, checkPointNum);

  sprintf(buffer,"%s.ckp.%d",name,checkPointNum);
  inStream = File_Open(buffer,"r",TRUE);
  graph = LoadScaffoldGraphFromStream(inStream);
  fclose(inStream);

  if(graph->checkPointIteration != (checkPointNum + 1)){
    fprintf(stderr,"**** Loaded checkpoint has checkpoint iteration %d, but we tried to load checkpoint %d...fishy!\n",
	    graph->checkPointIteration, checkPointNum);
    graph->checkPointIteration = checkPointNum + 1;
  }

  AssertPtr(graph);
  return graph;
}



void CheckpointScaffoldGraph(ScaffoldGraphT *graph, int logicalCheckpoint){
  char buffer[1024];
  FILE *outStream;
  char *name = GlobalData->File_Name_Prefix;
  time_t t;

  sprintf(buffer,"%s.ckp.%d",name,graph->checkPointIteration++);

  t = time(0);
  fprintf(GlobalData->timefp,"\n");
  fprintf(GlobalData->timefp, "====> Saving %s at %s", buffer, ctime(&t));

  outStream = File_Open(buffer,"w",TRUE);

  SaveScaffoldGraphToStream(graph, outStream);

  fclose(outStream);

  t = time(0);
  fprintf(GlobalData->timefp, "====> Done with checkpoint %d (logical %d) at %s", graph->checkPointIteration - 1, logicalCheckpoint, ctime(&t));
  fprintf(GlobalData->timefp, "\n");
}


void SaveScaffoldGraphToStream(ScaffoldGraphT *sgraph, FILE *stream){
  int status;

  AS_UTL_safeWrite(stream, sgraph->name, "SaveScaffoldGraphToStream", sizeof(char), 256);


  CopyToFileVA_InfoByIID(sgraph->iidToFragIndex, stream);
  CopyToFileVA_CIFragT(sgraph->CIFrags, stream);
#ifdef AS_ENABLE_SOURCE
  CopyToFileVA_char(sgraph->SourceFields, stream);
#endif
  CopyToFileVA_DistT(sgraph->Dists, stream);

  SaveGraphCGWToStream(sgraph->CIGraph,stream);
  SaveGraphCGWToStream(sgraph->ContigGraph,stream);
  SaveGraphCGWToStream(sgraph->ScaffoldGraph,stream);

  saveSequenceDB(sgraph->sequenceDB);

  AS_UTL_safeWrite(stream, &sgraph->doRezOnContigs, "SaveScaffoldGraphToStream", sizeof(int32), 1);
  AS_UTL_safeWrite(stream, &sgraph->checkPointIteration, "SaveScaffoldGraphToStream", sizeof(int32), 1);
  AS_UTL_safeWrite(stream, &sgraph->numContigs, "SaveScaffoldGraphToStream", sizeof(int32), 1);
  AS_UTL_safeWrite(stream, &sgraph->numDiscriminatorUniqueCIs, "SaveScaffoldGraphToStream", sizeof(int32), 1);
  AS_UTL_safeWrite(stream, &sgraph->numOriginalCIs, "SaveScaffoldGraphToStream", sizeof(int32), 1);
  AS_UTL_safeWrite(stream, &sgraph->numLiveCIs, "SaveScaffoldGraphToStream", sizeof(int32), 1);
  AS_UTL_safeWrite(stream, &sgraph->numLiveScaffolds, "SaveScaffoldGraphToStream", sizeof(int32), 1);
}

ScaffoldGraphT * LoadScaffoldGraphFromStream(FILE *stream){
  ScaffoldGraphT *sgraph =
    (ScaffoldGraphT *)safe_calloc(1, sizeof(ScaffoldGraphT));
  int    i, status;

  ScaffoldGraph = sgraph;
  sgraph->sequenceDB = SequenceDB;

  sgraph->gkpStore = openGateKeeperStore(GlobalData->Gatekeeper_Store_Name, FALSE);
  if(sgraph->gkpStore == NULL){
    fprintf(stderr,"**** Failure to open gatekeeper store %s ...exiting\n", GlobalData->Gatekeeper_Store_Name);
    exit(1);
  }

  status = AS_UTL_safeRead(stream, sgraph->name, "LoadScaffoldGraphFromStream", sizeof(char), 256);
  assert(status == 256);
  fprintf(GlobalData->stderrc,"* Reading graph %s *\n", sgraph->name);


  sgraph->iidToFragIndex = CreateFromFileVA_InfoByIID(stream);
  sgraph->CIFrags        = CreateFromFileVA_CIFragT(stream);
#ifdef AS_ENABLE_SOURCE
  sgraph->SourceFields   = CreateFromFileVA_char(stream);
#endif
  sgraph->Dists          = CreateFromFileVA_DistT(stream);

  sgraph->CIGraph       = LoadGraphCGWFromStream(stream);
  sgraph->ContigGraph   = LoadGraphCGWFromStream(stream);
  sgraph->ScaffoldGraph = LoadGraphCGWFromStream(stream);

  CheckGraph(sgraph->CIGraph);
  CheckGraph(sgraph->ContigGraph);
  CheckGraph(sgraph->ScaffoldGraph);

  //  Erase the bogus pointer in dists
  //
  for (i=0; i<GetNumDistTs(ScaffoldGraph->Dists); i++)
    GetDistT(ScaffoldGraph->Dists, i)->histogram = NULL;

  // Temporary
  sgraph->ChunkInstances = sgraph->CIGraph->nodes;
  sgraph->Contigs        = sgraph->ContigGraph->nodes;
  sgraph->CIScaffolds    = sgraph->ScaffoldGraph->nodes;
  sgraph->CIEdges        = sgraph->CIGraph->edges;
  sgraph->ContigEdges    = sgraph->ContigGraph->edges;
  sgraph->SEdges         = sgraph->ScaffoldGraph->edges;
  sgraph->overlapper     = sgraph->CIGraph->overlapper;


  status  = AS_UTL_safeRead(stream, &sgraph->doRezOnContigs, "LoadScaffoldGraphFromStream", sizeof(int32), 1);
  status += AS_UTL_safeRead(stream, &sgraph->checkPointIteration, "LoadScaffoldGraphFromStream", sizeof(int32), 1);
  status += AS_UTL_safeRead(stream, &sgraph->numContigs,"LoadScaffoldGraphFromStream",  sizeof(int32), 1);
  status += AS_UTL_safeRead(stream, &sgraph->numDiscriminatorUniqueCIs, "LoadScaffoldGraphFromStream", sizeof(int32), 1);
  status += AS_UTL_safeRead(stream, &sgraph->numOriginalCIs, "LoadScaffoldGraphFromStream", sizeof(int32), 1);
  status += AS_UTL_safeRead(stream, &sgraph->numLiveCIs, "LoadScaffoldGraphFromStream", sizeof(int32), 1);
  status += AS_UTL_safeRead(stream, &sgraph->numLiveScaffolds, "LoadScaffoldGraphFromStream", sizeof(int32), 1);
  assert(status == 7);

  if(sgraph->doRezOnContigs){
    sgraph->RezGraph = sgraph->ContigGraph;
  }else{
    sgraph->RezGraph = sgraph->CIGraph;
  }

  ReportMemorySize(sgraph,GlobalData->stderrc);

  SetCIScaffoldTLengths(sgraph, TRUE);
  CheckCIScaffoldTs(sgraph);

  return sgraph;
}





/***************************************************************************/
// TEMPORARY HACK
//
void InsertRepeatCIsInScaffolds(ScaffoldGraphT *sgraph){
  CDS_CID_t cid;
  int scaffolded = 0;
  LengthT NullLength = {0,0.0};
  CIScaffoldT CIScaffold = {0};
  GraphNodeIterator nodes = {0};
  NodeCGW_T *CI;

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
  CIScaffold.aEndCoord = CIScaffold.bEndCoord = -1;

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

  fprintf(GlobalData->stderrc,
          "* Inserting %d Unscaffolded CIs into newly created scaffolds\n",
	  scaffolded);

}

/****************************************************************************/

ScaffoldGraphT *CreateScaffoldGraph(int rezOnContigs, char *name,
                                    int32 numNodes, int32 numEdges){
  char buffer[2000];
  ScaffoldGraphT *sgraph = (ScaffoldGraphT *)safe_calloc(1,sizeof(ScaffoldGraphT));
  int numFrags = 0;
  int numDists = 0;
  sgraph->name[0] = '\0';
  strcpy(sgraph->name, name);
  fprintf(GlobalData->stderrc,"* Created scaffold graph %s\n", sgraph->name);

  sprintf(buffer,"%s.SeqStore", name);
  sgraph->sequenceDB = createSequenceDB(buffer);

#ifdef AS_ENABLE_SOURCE
  sgraph->SourceFields = CreateVA_char(1024);
#endif

  sgraph->CIGraph = CreateGraphCGW(CI_GRAPH, numNodes, numEdges);

  sgraph->ContigGraph = CreateGraphCGW(CONTIG_GRAPH, 1, 1);
  sgraph->ScaffoldGraph = CreateGraphCGW(SCAFFOLD_GRAPH, NULLINDEX, NULLINDEX);

  sgraph->gkpStore = openGateKeeperStore(GlobalData->Gatekeeper_Store_Name, FALSE);
  if(sgraph->gkpStore == NULL){
    fprintf(stderr,"**** Failure to open gatekeeper store %s ...exiting\n",buffer);
    exit(1);
  }

  numFrags = getNumGateKeeperFragments(sgraph->gkpStore);
  numDists = getNumGateKeeperLibraries(sgraph->gkpStore);

  sgraph->iidToFragIndex = CreateVA_InfoByIID(numFrags);
  sgraph->Dists = CreateVA_DistT(numDists);
  sgraph->CIFrags = CreateVA_CIFragT(numFrags);

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


  sgraph->checkPointIteration = 0;
  sgraph->numContigs = 0;
  sgraph->numDiscriminatorUniqueCIs = 0;
  sgraph->numLiveScaffolds = 0;
  return sgraph;
}

/* Destructor */
void DestroyScaffoldGraph(ScaffoldGraphT *sgraph){
  int i;

  deleteSequenceDB(sgraph->sequenceDB);

  DeleteGraphCGW(sgraph->CIGraph);
  DeleteGraphCGW(sgraph->ContigGraph);
  DeleteGraphCGW(sgraph->ScaffoldGraph);

  closeGateKeeperStore(sgraph->gkpStore);

  DeleteVA_CIFragT(sgraph->CIFrags);
#ifdef AS_ENABLE_SOURCE
  DeleteVA_char(sgraph->SourceFields);
#endif

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
#ifdef AS_ENABLE_SOURCE
  totalMemorySize += ReportMemorySize_VA(graph->SourceFields, "SourceFields", stream);
#endif
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
  double scaffoldMinPos = (double)CDS_COORD_MAX;
  double scaffoldMaxPos = (double)CDS_COORD_MIN;
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
    fprintf(GlobalData->stderrc,
            "*!!! Sanity  scaffold " F_CID " length (%g,%g) screwed up!\n",
            scaffold->id,
            scaffold->bpLength.mean, scaffold->bpLength.variance);
  }
  if(scaffold->info.Scaffold.AEndCI != NULLINDEX &&
     scaffold->bpLength.mean > 0 &&
     abs(scaffold->bpLength.mean - (scaffoldMaxPos - scaffoldMinPos)) > 100.0){
    fprintf(GlobalData->stderrc,
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
      fprintf(GlobalData->stderrc,
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
    fprintf(GlobalData->stderrc,
            "* numElements = %d scaffold says it has %d elements\n",
            numElements, scaffold->info.Scaffold.numElements);
    scaffoldInsane = 1;
  }

  if (scaffoldInsane) {
    DumpCIScaffold(GlobalData->stderrc,graph, scaffold, FALSE);
    assert(!scaffoldInsane);
  }

  // CheckScaffoldOrder(scaffold, graph);
}




// CheckScaffoldOrder checks whether all the ahangs in a scaffold are
// positive (ie, that the contigs are ordered by their distance from the
// A end of the scaffold)
void CheckScaffoldOrder(CIScaffoldT *scaffold, ScaffoldGraphT *graph)
{
  double currentMinPos = (double)CDS_COORD_MAX;
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
          fprintf( GlobalData->stderrc,
                   "CIs " F_CID " and " F_CID "are out of order\n",
                   CI->id, prevCI->id);
          fprintf( GlobalData->stderrc,
                   "CI " F_CID ": AEndOffset.mean: %f, AEndOffset.mean: %f\n",
                   CI->id, CI->offsetAEnd.mean, CI->offsetBEnd.mean);
          fprintf( GlobalData->stderrc,
                   "CI " F_CID ": AEndOffset.mean: %f, AEndOffset.mean: %f\n",
                   prevCI->id, prevCI->offsetAEnd.mean, prevCI->offsetBEnd.mean);
          DumpCIScaffold(GlobalData->stderrc, graph, scaffold, FALSE);

          // allow for a base of rounding error, but fix it
          if( MIN( CI->offsetAEnd.mean, CI->offsetBEnd.mean) - currentMinPos < 1.0)
            {
              CI->offsetAEnd.mean += 1.0;
              CI->offsetBEnd.mean += 1.0;
              fprintf( GlobalData->stderrc,
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



void DumpScaffoldGraph(ScaffoldGraphT *graph){
  static int version = 0;
  char fileName[FILENAME_MAX];
  FILE *dumpFile = NULL;

  sprintf(fileName,"%s.sg%c.%d",
          graph->name, (graph->doRezOnContigs?'C':'c'),version++);

  fprintf(GlobalData->stderrc,"* Dumping graph file %s\n", fileName);

  dumpFile = fopen(fileName,"w");

  //  DumpGraph(ScaffoldGraph->CIGraph, dumpFile);
  DumpGraph(ScaffoldGraph->RezGraph, dumpFile);
  DumpGraph(ScaffoldGraph->ScaffoldGraph, dumpFile);

  fclose(dumpFile);
}


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
  CIScaffoldT *scaffold;
  CIScaffoldTIterator Nodes;
  NodeCGW_T *thisNode;

#ifdef DEBUG_DETAILED
  fprintf(GlobalData->stderrc,
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
  }

  scaffold->bpLength.mean += delta.mean;
  scaffold->bpLength.variance += delta.variance;
  return;
}




static
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
    ReportMemorySize(ScaffoldGraph,GlobalData->stderrc);
    fflush(GlobalData->stderrc);
    fprintf(GlobalData->stderrc,"**** BEFORE repeat rez ****\n");
    //
    // repeat Fill_Gaps until we are not able to insert anymore
    //
    CheckEdgesAgainstOverlapper(ScaffoldGraph->RezGraph);
    // commented out next call to allow mouse_20010307 run to succeed
    CheckCITypes(ScaffoldGraph);


    fprintf(stderr,"* Calling CheckCIScaffoldTs\n");
    fflush(stderr);
    CheckCIScaffoldTs(ScaffoldGraph);
    fprintf(stderr,"* Done with CheckCIScaffoldTs\n");
    fflush(stderr);

    do
      {
        normal_inserts = Fill_Gaps (GlobalData, name, repeatRezLevel, iter);
        if  (normal_inserts > 0)
          {
            didSomething = TRUE;
            fprintf(stderr,"* Calling CheckCIScaffoldTs\n");
            fflush(stderr);
            CheckCIScaffoldTs(ScaffoldGraph);
            fprintf(stderr,"* Done with CheckCIScaffoldTs\n");
            fflush(stderr);

            TidyUpScaffolds (ScaffoldGraph);
            CheckEdgesAgainstOverlapper(ScaffoldGraph->RezGraph);

            fprintf(stderr,"* Calling CheckCIScaffoldTs\n");
            fflush(stderr);
            CheckCIScaffoldTs(ScaffoldGraph);
            fprintf(stderr,"* Done with CheckCIScaffoldTs\n");
            fflush(stderr);

            GeneratePlacedContigGraphStats("rocks", iter);
            GenerateLinkStats(ScaffoldGraph->ContigGraph,"rocks",iter);
            GenerateScaffoldGraphStats("rocks",iter);
          }

        if  (DO_CONTAINED_ROCKS && iter == 0)
          {
            contained_inserts = Hurl_Contained_Rocks (name, repeatRezLevel, iter);
            if  (contained_inserts > 0)
              {
                didSomething = TRUE;
                TidyUpScaffolds (ScaffoldGraph);
                CheckEdgesAgainstOverlapper(ScaffoldGraph->RezGraph);
                GeneratePlacedContigGraphStats("controcks", iter);
                GenerateLinkStats(ScaffoldGraph->ContigGraph,"controcks",iter);
                GenerateScaffoldGraphStats("controcks",iter);
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
                GeneratePlacedContigGraphStats("contstones", iter);
                GenerateLinkStats(ScaffoldGraph->ContigGraph,"contstones",iter);
                GenerateScaffoldGraphStats("contstones",iter);
              }
          }
        else
          contained_stones = 0;

        iter++;
        if (iter >= MAX_REZ_ITERATIONS) {
          fprintf (GlobalData->stderrc,
                   "Maximum number of REZ iterations reached\n");
          break;
        }
      }  while  (normal_inserts + contained_inserts
                 + contained_stones > FILL_GAPS_THRESHHOLD);

    fprintf(GlobalData->stderrc,"**** AFTER repeat rez ****\n");
  }
  return didSomething;
}






void DumpScaffoldSnapshot(char *name){
  char temp[200];
  FILE *dumpfp;
  MarkMisplacedContigs();
  CelamyAssembly(name);
  CelamyCIScaffolds(name, ScaffoldGraph);
  sprintf(temp,"%s.dump", name);
  dumpfp = fopen(temp,"w");
  DumpContigs(dumpfp,ScaffoldGraph, FALSE);
  DumpCIScaffolds(dumpfp, ScaffoldGraph, FALSE);
  fclose(dumpfp);
}


/***************************************************************************/
void RebuildScaffolds(ScaffoldGraphT *ScaffoldGraph,
                      int markShakyBifurcations){

#ifdef DEBUG_BUCIS
  BuildUniqueCIScaffolds(ScaffoldGraph, markShakyBifurcations,TRUE);
  fprintf(GlobalData->stderrc,"** After BuildUniqueCIScaffolds **\n");
  DumpChunkInstances(GlobalData->stderrc, ScaffoldGraph,
                     FALSE, TRUE, TRUE, FALSE);
  DumpCIScaffolds(GlobalData->stderrc,ScaffoldGraph, TRUE);
#else
  BuildUniqueCIScaffolds(ScaffoldGraph, markShakyBifurcations, FALSE);
#endif
  CheckEdgesAgainstOverlapper(ScaffoldGraph->RezGraph);
  fprintf(GlobalData->stderrc,"* Report Memory Size in RebuildScaffolds\n");
  ReportMemorySize(ScaffoldGraph, GlobalData->stderrc);
  fflush(NULL);

#ifdef DEBUG_BUCIS
  LeastSquaresGapEstimates(ScaffoldGraph, TRUE, FALSE, TRUE,
			   CHECK_CONNECTIVITY, TRUE);
#else
  LeastSquaresGapEstimates(ScaffoldGraph, TRUE, FALSE, TRUE,
			   CHECK_CONNECTIVITY, FALSE);
#endif

#ifdef DEBUG_BUCIS
  if(markShakyBifurcations){
    DumpScaffoldSnapshot("InitialScaffolds");
  }
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

  BuildSEdges(ScaffoldGraph, FALSE);
  MergeAllGraphEdges(ScaffoldGraph->ScaffoldGraph, TRUE);

  //clearCacheSequenceDB(ScaffoldGraph->sequenceDB);
}



/***************************************************************************/
void BuildScaffoldsFromFirstPriniciples(ScaffoldGraphT *ScaffoldGraph,
                                        int skipInitialScaffolding){
  int changedByRepeatRez;

  GenerateScaffoldGraphStats("initial",1);

  if (skipInitialScaffolding) {

    //  This is a continuation of a previous partial result, or a
    //  continuation of rocks, in either case, we have a checkpoint,
    //  so we skip the creation of a new one after TidyUpScaffolds(),
    //  and just get on with our work.

    TidyUpScaffolds (ScaffoldGraph);

  } else {
    // This includes CleanupScaffolds
    RebuildScaffolds(ScaffoldGraph, TRUE); // Transitive reduction of RezGraph followed by construction of SEdges

    //  Hooray!  We have scaffolds!  Checkpoint!

    fprintf(GlobalData->timefp,"* Checkpoint %d, After Building Initial Unique CI Scaffolds and before Tidying up there are %d scaffolds\n",
            ScaffoldGraph->checkPointIteration,
            (int) GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph));
    if(GlobalData->dumpScaffoldSnapshots){
      char temp[2000];
      sprintf(temp,"Rebuild%d",ScaffoldGraph->checkPointIteration);
      DumpScaffoldSnapshot(temp);
    }
    CheckpointScaffoldGraph(ScaffoldGraph, CHECKPOINT_AFTER_BUILDING_SCAFFOLDS);
  }

  GenerateScaffoldGraphStats("initial",1);
  GeneratePlacedContigGraphStats("initial",1);
  GenerateLinkStats(ScaffoldGraph->ContigGraph,"initial",1);

#define  MAX_OUTER_REZ_ITERATIONS   10

  if(GlobalData->repeatRezLevel > 0){
    int iter = 0;
    int ctme = time(0);

    fprintf(GlobalData->stderrc,"** Running Level 1 Repeat Rez **\n");
    do{
      changedByRepeatRez = FALSE;
      CheckEdgesAgainstOverlapper(ScaffoldGraph->RezGraph);

      // commented out next call to allow mouse_20010307 run to succeed
      CheckCITypes(ScaffoldGraph);

      changedByRepeatRez = RepeatRez(GlobalData->repeatRezLevel,
                                     GlobalData->File_Name_Prefix);
      if(changedByRepeatRez){
	CheckCIScaffoldTs(ScaffoldGraph);
        // merge in stuff placed by rocks, assuming its position is correct!
	CleanupScaffolds(ScaffoldGraph, FALSE, NULLINDEX, FALSE);

        // Transitive reduction of RezGraph followed by construction of SEdges
	RebuildScaffolds(ScaffoldGraph, FALSE);

        //  This checkpoint used to be included in RebuildScaffolds,
        //  but checkpointing after every iteration generates far too
        //  many checkpoints.  On a large mammal, on Opteron 2.2GHz,
        //  an iteration takes 15 - 20 minutes.  Microbes take seconds.
        //
        //  So, if we've been running for 2 hours, AND we've not just completed
        //  the last iteration, checkpoint.
        //
        if ((time(0) - ctme > 120 * 60) && (iter+1 < MAX_OUTER_REZ_ITERATIONS)) {
          ctme = time(0);
          fprintf(GlobalData->timefp, "* After RebuildScaffolds Rocks %d there are %d scaffolds\n",
                  iter, (int)GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph));
          if(GlobalData->dumpScaffoldSnapshots){
            char temp[2000];
            sprintf(temp,"Rebuild%d",ScaffoldGraph->checkPointIteration);
            DumpScaffoldSnapshot(temp);
          }
          CheckpointScaffoldGraph(ScaffoldGraph, CHECKPOINT_AFTER_BUILDING_SCAFFOLDS);
        }

        iter ++;
      }

#ifdef DEBUG_BUCIS
      fprintf(GlobalData->stderrc,"** After Gap Filling **\n");
      DumpChunkInstances(GlobalData->stderrc, ScaffoldGraph, FALSE, TRUE, TRUE, FALSE);
      DumpCIScaffolds(GlobalData->stderrc,ScaffoldGraph, FALSE);
#endif
    }while(changedByRepeatRez && iter < MAX_OUTER_REZ_ITERATIONS);


    //  Hooray!  We have cleaned scaffolds!

    fprintf(GlobalData->timefp,"* Checkpoint %d, After Gap Filling\n",
            ScaffoldGraph->checkPointIteration);
    if(GlobalData->dumpScaffoldSnapshots){
      char temp[2000];
      sprintf(temp,"Rebuild%d",ScaffoldGraph->checkPointIteration);
      DumpScaffoldSnapshot(temp);
    }
    CheckpointScaffoldGraph(ScaffoldGraph, CHECKPOINT_AFTER_BUILDING_AND_CLEANING_SCAFFOLDS);
  }


#ifdef DEBUG_BUCIS
  fprintf(GlobalData->stderrc,"* After building scaffolds \n");
  DumpCIScaffolds(GlobalData->stderrc,ScaffoldGraph,FALSE);
#endif
}


