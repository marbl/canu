
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
static char *rcsid = "$Id: ScaffoldGraph_CGW.c,v 1.66 2012-08-08 19:25:48 brianwalenz Exp $";

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

void ClearChunkInstance(ChunkInstanceT *ci){
  memset(ci, 0, sizeof(ChunkInstanceT));
  ci->id = ci->scaffoldID = NULLINDEX;
  ci->info.CI.contigID = NULLINDEX;
  ci->flags.all = 0;
}

//  This is a ginormous hack.  These two variables are private globals in AS_CNS, but we need to set
//  them before anything in consensus works.  Hopefully, we can build a 'consensus' object and set
//  them at construction time.  Hopefully, we can also get some thread safety at the same
//  time....dreaming....
//
extern gkStore               *gkpStore;
extern MultiAlignStore       *tigStore;

void
LoadScaffoldGraphFromCheckpoint(char   *name,
                                int32   checkPointNum,
                                int     writable){
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

  int    status;

  status = AS_UTL_safeRead(F, ScaffoldGraph->name, "LoadScaffoldGraphFromCheckpoint", sizeof(char), 256);
  assert(status == 256);

  ScaffoldGraph->CIFrags        = CreateFromFileVA_CIFragT(F);
  ScaffoldGraph->Dists          = CreateFromFileVA_DistT(F);

  ScaffoldGraph->CIGraph       = LoadGraphCGWFromStream(F);
  ScaffoldGraph->ContigGraph   = LoadGraphCGWFromStream(F);
  ScaffoldGraph->ScaffoldGraph = LoadGraphCGWFromStream(F);

  ScaffoldGraph->ChunkOverlaps = LoadChunkOverlapperFromStream(F);

  CheckGraph(ScaffoldGraph->CIGraph);
  CheckGraph(ScaffoldGraph->ContigGraph);
  CheckGraph(ScaffoldGraph->ScaffoldGraph);

  //  Load distance estimate histograms
  //
  for (int32 i=0; i<GetNumDistTs(ScaffoldGraph->Dists); i++) {
    DistT *dptr = GetDistT(ScaffoldGraph->Dists, i);

    dptr->histogram = (int32 *)safe_malloc(sizeof(int32) * dptr->bnum);

    status = AS_UTL_safeRead(F, dptr->histogram, "LoadScaffoldGraphFromCheckpoint", sizeof(int32), dptr->bnum);
    assert(status == dptr->bnum);
  }

  // Temporary
  ScaffoldGraph->ChunkInstances = ScaffoldGraph->CIGraph->nodes;
  ScaffoldGraph->Contigs        = ScaffoldGraph->ContigGraph->nodes;
  ScaffoldGraph->CIScaffolds    = ScaffoldGraph->ScaffoldGraph->nodes;
  ScaffoldGraph->CIEdges        = ScaffoldGraph->CIGraph->edges;
  ScaffoldGraph->ContigEdges    = ScaffoldGraph->ContigGraph->edges;
  ScaffoldGraph->SEdges         = ScaffoldGraph->ScaffoldGraph->edges;

  status  = AS_UTL_safeRead(F, &ScaffoldGraph->checkPointIteration,       "LoadScaffoldGraphFromCheckpoint", sizeof(int32), 1);
  status += AS_UTL_safeRead(F, &ScaffoldGraph->numContigs,                "LoadScaffoldGraphFromCheckpoint", sizeof(int32), 1);
  status += AS_UTL_safeRead(F, &ScaffoldGraph->numDiscriminatorUniqueCIs, "LoadScaffoldGraphFromCheckpoint", sizeof(int32), 1);
  status += AS_UTL_safeRead(F, &ScaffoldGraph->numOriginalCIs,            "LoadScaffoldGraphFromCheckpoint", sizeof(int32), 1);
  status += AS_UTL_safeRead(F, &ScaffoldGraph->numLiveCIs,                "LoadScaffoldGraphFromCheckpoint", sizeof(int32), 1);
  status += AS_UTL_safeRead(F, &ScaffoldGraph->numLiveScaffolds,          "LoadScaffoldGraphFromCheckpoint", sizeof(int32), 1);
  assert(status == 6);

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
  ScaffoldGraph->tigStore = tigStore = new MultiAlignStore(GlobalData->tigStoreName, checkPointNum, 0, 0, writable, FALSE);

  //  Open the gkpStore
  ScaffoldGraph->gkpStore = gkpStore = new gkStore(GlobalData->gkpStoreName, FALSE, writable);

  //  Do NOT check and cleanup scaffolds on load.  Do that BEFORE we save!

  fclose(F);
}



void
CheckpointScaffoldGraph(const char *logicalname, const char *location) {
  char ckpfile[FILENAME_MAX];
  char tmgfile[FILENAME_MAX];

  //  Check (and cleanup?) scaffolds.  This is done BEFORE we save, so that we (you know) SAVE the
  //  corrections made.  This solves a (rare?) problem in terminator where it would want to merge
  //  two contigs, but fails because everything (tigStore) is opened read-only.
  //
  //  The last checkpoint doesn't have a tigStore available as we already partitioned it and closed
  //  it, and so we cannot do the checks.  We don't expect any changes to have been made since the
  //  last checkpoint, so no real loss.  See AS_CGW_main.c, not that it'll help the causal observer
  //  any.
  //
  if (ScaffoldGraph->tigStore) {
    SetCIScaffoldTLengths(ScaffoldGraph);
    ScaffoldSanity(ScaffoldGraph);
  }

  sprintf(ckpfile, "%s.ckp.%d", GlobalData->outputPrefix, ScaffoldGraph->checkPointIteration++);
  sprintf(tmgfile, "%s.timing", GlobalData->outputPrefix);

  errno = 0;
  FILE *F = fopen(ckpfile, "w");
  if (errno)
    fprintf(stderr, "Failed to open '%s' for writing checkpoint: %s\n", ckpfile, strerror(errno)), exit(1);

  AS_UTL_safeWrite(F, ScaffoldGraph->name, "CheckpointScaffoldGraph", sizeof(char), 256);

  CopyToFileVA_CIFragT(ScaffoldGraph->CIFrags, F);
  CopyToFileVA_DistT(ScaffoldGraph->Dists, F);

  SaveGraphCGWToStream(ScaffoldGraph->CIGraph,F);
  SaveGraphCGWToStream(ScaffoldGraph->ContigGraph,F);
  SaveGraphCGWToStream(ScaffoldGraph->ScaffoldGraph,F);

  SaveChunkOverlapperToStream(ScaffoldGraph->ChunkOverlaps, F);

  //  Save the distance estimate histograms -- terminator needs these to output
  //
  for (int32 i=0; i<GetNumDistTs(ScaffoldGraph->Dists); i++) {
    DistT *dptr = GetDistT(ScaffoldGraph->Dists, i);

    AS_UTL_safeWrite(F, dptr->histogram, "LoadScaffoldGraphFromCheckpoint", sizeof(int32), dptr->bnum);
  }

  AS_UTL_safeWrite(F, &ScaffoldGraph->checkPointIteration,       "CheckpointScaffoldGraph", sizeof(int32), 1);
  AS_UTL_safeWrite(F, &ScaffoldGraph->numContigs,                "CheckpointScaffoldGraph", sizeof(int32), 1);
  AS_UTL_safeWrite(F, &ScaffoldGraph->numDiscriminatorUniqueCIs, "CheckpointScaffoldGraph", sizeof(int32), 1);
  AS_UTL_safeWrite(F, &ScaffoldGraph->numOriginalCIs,            "CheckpointScaffoldGraph", sizeof(int32), 1);
  AS_UTL_safeWrite(F, &ScaffoldGraph->numLiveCIs,                "CheckpointScaffoldGraph", sizeof(int32), 1);
  AS_UTL_safeWrite(F, &ScaffoldGraph->numLiveScaffolds,          "CheckpointScaffoldGraph", sizeof(int32), 1);

  fclose(F);

  if (ScaffoldGraph->tigStore)
    ScaffoldGraph->tigStore->nextVersion();

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
  GraphNodeIterator nodes;
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

  InitGraphNodeIterator(&nodes, sgraph->ContigGraph, GRAPH_NODE_DEFAULT);
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

ScaffoldGraphT *CreateScaffoldGraph(char *name) {
  char buffer[2000];
  ScaffoldGraphT *sgraph = (ScaffoldGraphT *)safe_calloc(1,sizeof(ScaffoldGraphT));

  strcpy(sgraph->name, name);

  int32 tigStoreInputVersion = 5;

  sgraph->checkPointIteration = tigStoreInputVersion + 1;
  sgraph->tigStore      = tigStore = new MultiAlignStore(GlobalData->tigStoreName,
                                                         tigStoreInputVersion,
                                                         0,
                                                         0,
                                                         TRUE, FALSE);

  sgraph->CIGraph       = CreateGraphCGW(CI_GRAPH, 16 * 1024, 16 * 1024);
  sgraph->ContigGraph   = CreateGraphCGW(CONTIG_GRAPH, 1, 1);
  sgraph->ScaffoldGraph = CreateGraphCGW(SCAFFOLD_GRAPH, NULLINDEX, NULLINDEX);

  sgraph->ChunkOverlaps = CreateChunkOverlapper();

  sgraph->gkpStore      = gkpStore = new gkStore(GlobalData->gkpStoreName, FALSE, TRUE);

  int numFrags = sgraph->gkpStore->gkStore_getNumFragments();
  int numDists = sgraph->gkpStore->gkStore_getNumLibraries();

  sgraph->Dists          = CreateVA_DistT(numDists);
  sgraph->CIFrags        = CreateVA_CIFragT(numFrags);

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

  delete sgraph->tigStore;

  DeleteGraphCGW(sgraph->CIGraph);
  DeleteGraphCGW(sgraph->ContigGraph);
  DeleteGraphCGW(sgraph->ScaffoldGraph);

  DestroyChunkOverlapper(sgraph->ChunkOverlaps);

  for (int32 i=0; i<GetNumDistTs(ScaffoldGraph->Dists); i++) {
    DistT *dptr = GetDistT(ScaffoldGraph->Dists, i);
    safe_free(dptr->histogram);
  }

  delete sgraph->gkpStore;

  DeleteVA_CIFragT(sgraph->CIFrags);

  DeleteVA_DistT(sgraph->Dists);

  safe_free(sgraph);
}


void ReportMemorySize(ScaffoldGraphT *graph, FILE *stream){

  size_t totalMemorySize = 0;

  fprintf(stream,"* ScaffoldGraph Memory Report:\n");
  fprintf(stream,"*\t     \tnumElements\tAllocatedElements\tAllocatedSize\n");
  totalMemorySize += ReportMemorySize_VA(graph->Dists, "Dists", stream);
  totalMemorySize += ReportMemorySize_VA(graph->CIFrags, "CIFrags", stream);
  totalMemorySize += ReportMemorySizeGraphCGW(graph->CIGraph, stream);
  totalMemorySize += ReportMemorySizeGraphCGW(graph->ContigGraph, stream);
  totalMemorySize += ReportMemorySizeGraphCGW(graph->ScaffoldGraph, stream);
  fprintf(stream,"*\tTotalMemorySize = " F_SIZE_T "\n", totalMemorySize);
}




#if 0
void DumpScaffoldGraph(ScaffoldGraphT *graph){
  static int version = 0;
  char fileName[FILENAME_MAX];
  FILE *dumpFile = NULL;

  sprintf(fileName,"%s.sg.%d",
          graph->name, version++);

  fprintf(stderr,"* Dumping graph file %s\n", fileName);

  dumpFile = fopen(fileName,"w");

  //  DumpGraph(ScaffoldGraph->CIGraph, dumpFile);
  DumpGraph(ScaffoldGraph->ContigGraph, dumpFile);
  DumpGraph(ScaffoldGraph->ScaffoldGraph, dumpFile);

  fclose(dumpFile);
}
#endif


//  If this is a ChunkInstance, return its coverage stat.  If this is an singleton contig, return
//  its lone ChunkInstance's coverage stat else, assert.
//
int
GetCoverageStat(ChunkInstanceT *CI) {

  if (CI->flags.bits.isCI)
    return ScaffoldGraph->tigStore->getUnitigCoverageStat(CI->id);

  if (CI->flags.bits.isContig) {
    if (CI->info.Contig.numCI == 1)
      return ScaffoldGraph->tigStore->getUnitigCoverageStat(CI->info.Contig.AEndCI);
    else
      //  Multi-unitig contig, unique-equivalent
      return GlobalData->cgbUniqueCutoff;
  }
  assert(0);
  return(0);
}


//  If this is a ChunkInstance, return the number of surrogate instances.  If this is an singleton
//  contig, return its lone ChunkInstance's number of surrogate instances, else, return 0.
//
int
GetNumInstances(ChunkInstanceT *CI) {

  if (CI->flags.bits.isCI)
    return CI->info.CI.numInstances;

  if (CI->flags.bits.isContig) {
    if (CI->info.Contig.numCI == 1) {
      ChunkInstanceT *ci = GetGraphNode(ScaffoldGraph->CIGraph, CI->info.Contig.AEndCI);
      return ci->info.CI.numInstances;
    } else {
      //  Multi-unitig contig, not a surrogate.
      return 0;
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
                               LengthT delta,
                               uint32 mark){
  CIScaffoldT *scaffold;
  CIScaffoldTIterator Nodes;
  NodeCGW_T *thisNode;

  scaffold = GetGraphNode(graph->ScaffoldGraph, scaffoldIndex);

  InitCIScaffoldTIteratorFromCI(graph, scaffold, indexOfCI, aEndToBEnd, FALSE, &Nodes);
  while((thisNode = NextCIScaffoldTIterator(&Nodes)) != NULL){
    thisNode->offsetAEnd.mean += delta.mean;
    thisNode->offsetAEnd.variance += delta.variance;
    thisNode->offsetBEnd.mean += delta.mean;
    thisNode->offsetBEnd.variance += delta.variance;

    thisNode->flags.bits.isJiggled = mark;
    thisNode->offsetDelta.mean = delta.mean;
    thisNode->offsetDelta.variance = delta.variance;
  }

  scaffold->bpLength.mean += delta.mean;
  scaffold->bpLength.variance += delta.variance;
}




void
CheckCITypes(ScaffoldGraphT *sgraph){
  NodeCGW_T *CI;
  GraphNodeIterator nodes;

  //  An alternate iterator
  //  for(i = 0; i < GetNumGraphNodes(sgraph->ContigGraph); i++)
  //    ChunkInstanceT *CI = GetGraphNode(sgraph->ContigGraph,i);

  InitGraphNodeIterator(&nodes, sgraph->ContigGraph, GRAPH_NODE_DEFAULT);
  while((CI = NextGraphNodeIterator(&nodes)) != NULL){
    if(CI->flags.bits.isUnique ){
      if(CI->scaffoldID == NULLINDEX){
        // This code used to print several warning messages when it found a unique Contig not in a Scaffold
        // However, the code being used to print was not up to date with all CGW changes and was incorrect
        // Now this code asserts
	     assert(0);
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
    CheckEdgesAgainstOverlapper(ScaffoldGraph->ContigGraph);
    // commented out next call to allow mouse_20010307 run to succeed
    CheckCITypes(ScaffoldGraph);
    ScaffoldSanity(ScaffoldGraph);

    do
      {
        normal_inserts = Fill_Gaps (name, repeatRezLevel, iter);
        if  (normal_inserts > 0)
          {
            didSomething = TRUE;

            ScaffoldSanity(ScaffoldGraph);

            TidyUpScaffolds (ScaffoldGraph);
            CheckEdgesAgainstOverlapper(ScaffoldGraph->ContigGraph);

            ScaffoldSanity(ScaffoldGraph);

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
                CheckEdgesAgainstOverlapper(ScaffoldGraph->ContigGraph);
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
                CheckEdgesAgainstOverlapper(ScaffoldGraph->ContigGraph);
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




static
void
CheckAllTrustedEdges(ScaffoldGraphT * sgraph) {
  GraphNodeIterator nodes;
  ChunkInstanceT   *contig;

  InitGraphNodeIterator(&nodes, sgraph->ContigGraph, GRAPH_NODE_DEFAULT);
  while ((contig = NextGraphNodeIterator(&nodes)) != NULL) {
    if (contig->scaffoldID == NULLINDEX)
      continue;

    ChunkInstanceT *thisC = GetGraphNode(sgraph->ContigGraph, contig->id);

    assert(contig == thisC);

    GraphEdgeIterator  edges(sgraph->ContigGraph, contig->id, ALL_END, ALL_TRUSTED_EDGES);
    CIEdgeT           *edge;

    while((edge = edges.nextMerged()) != NULL){
      CDS_CID_t next = (contig->id == edge->idA) ? edge->idB : edge->idA;

      ChunkInstanceT *thatC = GetGraphNode(ScaffoldGraph->ContigGraph, next);

      if (thatC->scaffoldID != thisC->scaffoldID)
        fprintf(stderr,"-=> BAD edge id:"F_SIZE_T" "F_CID"("F_CID")->"F_CID"("F_CID") (weight %d, status %d)\n",
                GetVAIndex_CIEdgeT(sgraph->ContigGraph->edges, edge),
                contig->id, thisC->scaffoldID,
                thatC->id, thatC->scaffoldID,
                edge->edgesContributing, edge->flags.bits.edgeStatus);
    }
  }
}

void
TidyUpScaffolds(ScaffoldGraphT *ScaffoldGraph) {

  CleanupScaffolds(ScaffoldGraph, FALSE, NULLINDEX, FALSE);

  for (int32 sID=0; sID < GetNumCIScaffoldTs(ScaffoldGraph->CIScaffolds); sID++)
    LeastSquaresGapEstimates(ScaffoldGraph, GetCIScaffoldT(ScaffoldGraph->CIScaffolds, sID));

  if(GlobalData->debugLevel > 0)
    CheckAllContigFragments();

  CheckAllTrustedEdges(ScaffoldGraph);

  BuildSEdges(ScaffoldGraph, TRUE, FALSE);
  MergeAllGraphEdges(ScaffoldGraph->ScaffoldGraph, TRUE, FALSE);
}
