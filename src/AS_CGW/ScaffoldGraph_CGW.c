
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
static char CM_ID[] = "$Id: ScaffoldGraph_CGW.c,v 1.1.1.1 2004-04-14 13:51:01 catmandew Exp $";

//#define DEBUG 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <sys/types.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "UtilsREZ.h"
#include "AS_PER_SafeIO.h"
#include "AS_CGW_dataTypes.h"
#include "Globals_CGW.h"
#include "DiagnosticsCGW.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"
#include "RepeatRez.h"
#include "CommonREZ.h"
#include "GreedyOverlapREZ.h"
#include "Stats_CGW.h"

ScaffoldGraphT *ScaffoldGraph = NULL;
tSequenceDB *SequenceDB = NULL;

void ClearChunkInstance(ChunkInstanceT *ci){
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
  SequenceDB = OpenSequenceDB(buffer, readWrite, checkPointNum);

  sprintf(buffer,"%s.ckp.%d",name,checkPointNum);
  inStream = File_Open(buffer,"r",TRUE);
  graph = LoadScaffoldGraphFromStream(inStream);
  fclose(inStream);

#if 0
  sprintf(buffer,"%s.ckpl",name);
  fprintf(GlobalData->stderrc,
          "* Trying to open gatekeeper link store %s\n", buffer);
  graph->gkplStore = loadStore(buffer);
#endif
  if(graph->checkPointIteration != (checkPointNum + 1)){
    fprintf(stderr,"**** Loaded checkpoint has checkpoint iteration %d, but we tried to load checkpoint %d...fishy!\n",
	    graph->checkPointIteration, checkPointNum);
            graph->checkPointIteration = checkPointNum + 1;
  }
  
  AssertPtr(graph);
  return graph;
}


void CheckpointOnDemand(int whatToDoAfter)
{
  FILE * fp;
  if((fp = fopen(CHECKPOINT_DEMAND_FILE, "r")) != NULL)
  {
    char command[1024];
    
    fclose(fp);
    
    CleanupScaffolds(ScaffoldGraph, FALSE, NULLINDEX, FALSE);
    fprintf( GlobalData->stderrc, "Checkpoint %d written during MergeScaffoldsAggressive demanded by user\n", ScaffoldGraph->checkPointIteration);
    fprintf( GlobalData->timefp,"Checkpoint %d written during MergeScaffoldsAggressive demanded by user\n", ScaffoldGraph->checkPointIteration);
    CheckpointScaffoldGraph(ScaffoldGraph);
    
    sprintf(command, "rm -f %s", CHECKPOINT_DEMAND_FILE);
    if(system(command) != 0)
    {
      fprintf(GlobalData->stderrc,
              "ALERT!!!!! Failed to remove checkpoint-on-demand file, %s\n",
              CHECKPOINT_DEMAND_FILE);
      fprintf(GlobalData->stderrc, "Please remove it immediately!!!!\n");
    }

    switch(whatToDoAfter)
    {
      case EXIT_AFTER_CHECKPOINTING:
        fprintf(GlobalData->stderrc,
                "FYI: Exiting after checkpointing on user demand!!\n");
        exit(0);
        break;
      case RETURN_AFTER_CHECKPOINTING:
        return;
        break;
      default:
        fprintf(GlobalData->stderrc,
                "WARNING: Checkpoint written on user demand with invalid instruction.\n");
        fprintf(GlobalData->stderrc,
                "\t\tOptions are to exit or continue running.\n");
        fprintf(GlobalData->stderrc,
                "Continuing with run.\n");
        break;
    }
  }
}


void CheckpointScaffoldGraph(ScaffoldGraphT *graph){
  char buffer[1024];
  FILE *outStream;
  char *name = GlobalData->File_Name_Prefix;

#if 0
  sprintf(buffer,"%s.ckpl",name);

  storeStore(graph->gkplStore, buffer);
#endif

  sprintf(buffer,"%s.ckp.%d",name,graph->checkPointIteration++);
  
  {
    time_t t;
    t = time(0);
    fprintf(GlobalData->timefp,">>>>*************************************************************************<<<<\n");    
    fprintf(GlobalData->timefp, "====> Saving %s at %s\n", buffer, ctime(&t));
  }
  {
    long cycles;
    fprintf(GlobalData->timefp,"******* Dumping Checkpoint %d \n", graph->checkPointIteration - 1);
    fprintf(GlobalData->timefp,"* Time in Chunk Selection %g seconds (%ld calls)\n",
	    TotalTimerT(&GlobalData->ChooseChunksTimer, &cycles), cycles);
    fprintf(GlobalData->timefp,"* Time in Consistency Check %g seconds (%ld calls)\n",
	    TotalTimerT(&GlobalData->ConsistencyCheckTimer, &cycles), cycles);
    fprintf(GlobalData->timefp,"* Time in Update %g seconds (%ld calls)\n",
	    TotalTimerT(&GlobalData->UpdateTimer, &cycles), cycles);
    fprintf(GlobalData->timefp,"* Time in WalkUpdate %g seconds (%ld calls)\n",
	    TotalTimerT(&GlobalData->WalkUpdateTimer, &cycles), cycles);
    fprintf(GlobalData->timefp,"* Time in RecomputeOffsets %g seconds (%ld calls)\n",
	    TotalTimerT(&GlobalData->RecomputeOffsetsTimer, &cycles), cycles);
    fprintf(GlobalData->timefp,"* Time in MergeScaffolds %g seconds (%ld calls)\n",
	    TotalTimerT(&GlobalData->MergeScaffoldsTimer, &cycles), cycles);
    fprintf(GlobalData->timefp,"* Time in Gap Fill %g seconds (%ld calls)\n",
	    TotalTimerT(&GlobalData->GapFillTimer, &cycles), cycles);
    fprintf(GlobalData->timefp,"* Time in Gap Walking %g seconds (%ld calls)\n",
	    TotalTimerT(&GlobalData->GapWalkerTimer, &cycles), cycles);
    fprintf(GlobalData->timefp,"* Time in Stone Throwing %g seconds (%ld calls)\n",
	    TotalTimerT(&GlobalData->StoneThrowingTimer, &cycles), cycles);
    fprintf(GlobalData->timefp,"* Time in Consensus %g seconds (%ld calls)\n",
	    TotalTimerT(&GlobalData->ConsensusTimer, &cycles), cycles);
  }
  fflush(NULL);

  outStream = File_Open(buffer,"w",TRUE);
  
  SaveScaffoldGraphToStream(graph, outStream);

  {
    time_t t;
    t = time(0);
    fprintf(GlobalData->timefp, "====> Done with checkpoint at %s\n", ctime(&t));
    fprintf(GlobalData->timefp,">>>>*************************************************************************<<<<\n\n\n");    
    fflush(NULL);
  }
  fclose(outStream);
}


void SaveScaffoldGraphToStream(ScaffoldGraphT *sgraph, FILE *stream){
  int status;

  //  fprintf(GlobalData->stderrc,"* Saving graph %s\n", sgraph->name);
  status = safeWrite(stream, sgraph->name, 256 * sizeof(char));
  assert(status == FALSE);


  CopyToFileVA_InfoByIID(sgraph->iidToFragIndex, stream);
  CopyToFileVA_CIFragT(sgraph->CIFrags, stream);
  CopyToFileVA_char(sgraph->SourceFields, stream);
  CopyToFileVA_DistT(sgraph->Dists, stream);
    
  //  fprintf(GlobalData->stderrc,"* Saving CIGraph\n");
  SaveGraphCGWToStream(sgraph->CIGraph,stream);
  //  fprintf(GlobalData->stderrc,"* Saving ContigGraph\n");
  SaveGraphCGWToStream(sgraph->ContigGraph,stream);
  SaveGraphCGWToStream(sgraph->ScaffoldGraph,stream);

  SaveSequenceDB(sgraph->sequenceDB);

  status = safeWrite(stream, &sgraph->doRezOnContigs, sizeof(int32));
  assert(status == FALSE);
  status = safeWrite(stream, &sgraph->checkPointIteration, sizeof(int32));
  assert(status == FALSE);
  status = safeWrite(stream, &sgraph->numContigs, sizeof(int32));
  assert(status == FALSE);
  status = safeWrite(stream, &sgraph->numDiscriminatorUniqueCIs, sizeof(int32));
  assert(status == FALSE);
  status = safeWrite(stream, &sgraph->numOriginalCIs, sizeof(int32));
  assert(status == FALSE);
  status = safeWrite(stream, &sgraph->numLiveCIs, sizeof(int32));
  assert(status == FALSE);
  status = safeWrite(stream, &sgraph->numLiveScaffolds, sizeof(int32));
  assert(status == FALSE);

}

ScaffoldGraphT * LoadScaffoldGraphFromStream(FILE *stream){
  ScaffoldGraphT *sgraph =
    (ScaffoldGraphT *)safe_calloc(1, sizeof(ScaffoldGraphT));
  int status;
  TimerT timer;

  InitTimerT(&timer);
  StartTimerT(&timer);
  
  ScaffoldGraph = sgraph;
  sgraph->sequenceDB = SequenceDB;

  if(strlen(GlobalData->Frag_Store_Name) > 1)
  {
    char buffer[256];
    sprintf(buffer,"%s", GlobalData->Frag_Store_Name);
    fprintf(GlobalData->stderrc,"* Trying to open %s\n", buffer);
#if 0
    sgraph->fragStore = loadFragStore( buffer);
#else
    sgraph->fragStore = openFragStore( buffer, "rw+");
#endif

    if(sgraph->fragStore == NULLSTOREHANDLE){
      fprintf(stderr,"**** Failure to open frag store %s ...exiting\n",buffer);
      exit(1);
    }else{
      fprintf(stderr,"*** Succeeded to open frag Store.\n");
    }

    sprintf(buffer,"%s", GlobalData->Gatekeeper_Store_Name);
    InitGateKeeperStore(&sgraph->gkpStore, buffer);
    if(OpenReadOnlyGateKeeperStore(&sgraph->gkpStore)){
      fprintf(stderr,"*** Failure to open Gatekeeper Store...exiting\n");
      exit(1);
    }else{
      fprintf(stderr,"*** Succeeded to open Gatekeeper Store.\n");
    }
  }
  else
    sgraph->fragStore = NULLSTOREHANDLE;

  status = safeRead(stream, sgraph->name, 256 * sizeof(char));
  assert(status == FALSE);
  fprintf(GlobalData->stderrc,"* Reading graph %s *\n", sgraph->name);


  sgraph->iidToFragIndex = CreateFromFileVA_InfoByIID( stream,0);
  sgraph->CIFrags =        CreateFromFileVA_CIFragT(stream,0);
  sgraph->SourceFields =  CreateFromFileVA_char(stream,0);
  sgraph->Dists =          CreateFromFileVA_DistT( stream,0);
  {
    int i;
    for( i = 0; i < GetNumDistTs(sgraph->Dists); i++){
      DistT *dist = GetDistT(sgraph->Dists,i);
      dist->samples = NULL;
    }
  }
  fflush(stderr);
  //  fprintf(GlobalData->stderrc,"* Loading CIGraph \n");
  sgraph->CIGraph = LoadGraphCGWFromStream(stream);
  //  fprintf(GlobalData->stderrc,"* Loading ContigGraph \n");
  sgraph->ContigGraph = LoadGraphCGWFromStream(stream);
  sgraph->ScaffoldGraph = LoadGraphCGWFromStream(stream);

  CheckGraph(sgraph->CIGraph);
  CheckGraph(sgraph->ContigGraph);
  CheckGraph(sgraph->ScaffoldGraph);

  // Temporary
  sgraph->ChunkInstances = sgraph->CIGraph->nodes;
  sgraph->Contigs        = sgraph->ContigGraph->nodes;
  sgraph->CIScaffolds    = sgraph->ScaffoldGraph->nodes;
  sgraph->CIEdges =        sgraph->CIGraph->edges;
  sgraph->ContigEdges    = sgraph->ContigGraph->edges;
  sgraph->SEdges         = sgraph->ScaffoldGraph->edges;
  sgraph->overlapper = sgraph->CIGraph->overlapper;


  if(strlen(GlobalData->Frag_Store_Name) > 1)
  {
    char buffer[256];
    sprintf(buffer,"%s", GlobalData->Frag_Store_Name);
    fprintf(GlobalData->stderrc,"* Trying to open %s\n", buffer);
#if 0
    sgraph->fragStore = loadFragStore( buffer);
#else
    sgraph->fragStore = openFragStore( buffer, "rw+");
#endif

    if(sgraph->fragStore == NULLSTOREHANDLE){
      fprintf(stderr,"**** Failure to open frag store %s ...exiting\n",buffer);
      exit(1);
    }

    sprintf(buffer,"%s", GlobalData->Gatekeeper_Store_Name);
    InitGateKeeperStore(&sgraph->gkpStore, buffer);
    if(OpenReadOnlyGateKeeperStore(&sgraph->gkpStore)){
      fprintf(stderr,"*** Failure to open Gatekeeper Store...exiting\n");
      exit(1);
    }
  }
  else
    sgraph->fragStore = NULLSTOREHANDLE;

  status = safeRead(stream, &sgraph->doRezOnContigs, sizeof(int32));
  assert(status == FALSE);
  status = safeRead(stream, &sgraph->checkPointIteration, sizeof(int32));
  assert(status == FALSE);
  status = safeRead(stream, &sgraph->numContigs, sizeof(int32));
  assert(status == FALSE);
  status = safeRead(stream, &sgraph->numDiscriminatorUniqueCIs, sizeof(int32));
  assert(status == FALSE);
  status = safeRead(stream, &sgraph->numOriginalCIs, sizeof(int32));
  assert(status == FALSE);
  status = safeRead(stream, &sgraph->numLiveCIs, sizeof(int32));
  assert(status == FALSE);
  status = safeRead(stream, &sgraph->numLiveScaffolds, sizeof(int32));
  assert(status == FALSE);

  if(sgraph->doRezOnContigs){
    sgraph->RezGraph = sgraph->ContigGraph;
  }else{
    sgraph->RezGraph = sgraph->CIGraph;
  }

  StopTimerT(&timer);
  fprintf(GlobalData->stderrc,
          "*** Loading Scaffold Graph %s took %g seconds\n",
	  sgraph->name, TotalTimerT(&timer, NULL));

  ReportMemorySize(sgraph,GlobalData->stderrc);

  fprintf(stderr,"* Calling SetCIScaffoldTLengths\n");
  SetCIScaffoldTLengths(sgraph, TRUE);
  fprintf(stderr,"* Calling CheckCIScaffoldTs\n");
  fflush(stderr);
  CheckCIScaffoldTs(sgraph);
  fprintf(stderr,"* Done with CheckCIScaffoldTs\n");
  fflush(stderr);

  return sgraph;

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
  sgraph->sequenceDB = CreateSequenceDB(buffer, numNodes, TRUE);

#if 0
  sgraph->gkplStore = createGateKeeperLinkStore(NULL, "lnk", 1); // in memory
#endif
  sgraph->SourceFields = CreateVA_char(1024);

  sgraph->CIGraph = CreateGraphCGW(CI_GRAPH, numNodes, numEdges);
  // Get a little space, we'll expand it later
  sgraph->ContigGraph = CreateGraphCGW(CONTIG_GRAPH, 1, 1);
  sgraph->ScaffoldGraph = CreateGraphCGW(SCAFFOLD_GRAPH, NULLINDEX, NULLINDEX);
  {
    char buffer[256];
    sprintf(buffer,"%s", GlobalData->Frag_Store_Name);
#if  0
    sgraph->fragStore = loadFragStore( buffer);
#else
    sgraph->fragStore = openFragStore( buffer, "rw+");
#endif
    if(sgraph->fragStore == NULLSTOREHANDLE){
      fprintf(stderr,"**** Failure to open frag store %s ...exiting\n",buffer);
      exit(1);
    }else{
      fprintf(stderr,"*** Succeeded to open frag Store.\n");
    }

    sprintf(buffer,"%s", GlobalData->Gatekeeper_Store_Name);
    InitGateKeeperStore(&sgraph->gkpStore, buffer);
    OpenReadOnlyGateKeeperStore(&sgraph->gkpStore);
    
    numFrags = getNumGateKeeperFragments(sgraph->gkpStore.frgStore);
    numDists = getNumGateKeeperDistances(sgraph->gkpStore.dstStore);
    sgraph->iidToFragIndex = CreateVA_InfoByIID(numFrags);
    sgraph->Dists = CreateVA_DistT(numDists);
    sgraph->CIFrags = CreateVA_CIFragT(numFrags);

  }

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
  sgraph->Contigs = sgraph->ContigGraph->nodes;
  sgraph->CIScaffolds = sgraph->ScaffoldGraph->nodes;
  sgraph->CIEdges = sgraph->CIGraph->edges;
  sgraph->ContigEdges = sgraph->ContigGraph->edges;
  sgraph->SEdges = sgraph->ScaffoldGraph->edges;

  sgraph->checkPointIteration = 0;
  sgraph->numContigs = 0;
  sgraph->numDiscriminatorUniqueCIs = 0;
  sgraph->numLiveScaffolds = 0;
  return sgraph;
}

/* Destructor */
void DestroyScaffoldGraph(ScaffoldGraphT *sgraph){
#if 0
  closeGateKeeperLinkStore(sgraph->gkplStore);
#endif
  DeleteSequenceDB(sgraph->sequenceDB);

  DeleteGraphCGW(sgraph->CIGraph);
  DeleteGraphCGW(sgraph->ContigGraph);
  DeleteGraphCGW(sgraph->ScaffoldGraph);

  if(sgraph->fragStore != NULLSTOREHANDLE)
  {
    CloseGateKeeperStore(&(sgraph->gkpStore));
    closeFragStore(sgraph->fragStore);
  }

  DeleteVA_CIFragT(sgraph->CIFrags);
  DeleteVA_char(sgraph->SourceFields);

  {
    int i;
    for( i = 0; i < GetNumDistTs(sgraph->Dists); i++){
      DistT *dist = GetDistT(sgraph->Dists,i);
      if(dist->samples)
	DeleteVA_CDS_COORD_t(dist->samples);
    }
  }
  DeleteVA_DistT(sgraph->Dists);
  DeleteVA_InfoByIID(sgraph->iidToFragIndex);

  free(sgraph);
}


void ReportMemorySize(ScaffoldGraphT *graph, FILE *stream){

  size_t totalMemorySize = 0;

  fprintf(stream,"* ScaffoldGraph Memory Report:\n");
  fprintf(stream,"*\t     \tnumElements\tAllocatedElements\tAllocatedSize\n");
  totalMemorySize += ReportMemorySize_VA(graph->Dists, "Dists", stream);
  totalMemorySize += ReportMemorySize_VA(graph->CIFrags, "CIFrags", stream);
  totalMemorySize += ReportMemorySize_VA(graph->SourceFields, "SourceFields", stream);
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
  assert(scaffold->flags.bits.isScaffold);
  if(scaffold->type != REAL_SCAFFOLD)
    return;
  
  InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &CIs);
  while(NULL != (CI = NextCIScaffoldTIterator(&CIs))){
    scratch = min(CI->offsetAEnd.mean, CI->offsetBEnd.mean);
    scaffoldMinPos = min(scaffoldMinPos, scratch);
    scratch = max(CI->offsetAEnd.mean, CI->offsetBEnd.mean);
    scaffoldMaxPos = max(scaffoldMaxPos, scratch);
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
    DumpCIScaffold(GlobalData->stderrc,ScaffoldGraph, scaffold, FALSE);
#ifdef STRICT_SCAFFOLD_CHECKING	   
    assert(0);
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
              "* Scaffold Sanity: node " F_CID " is screwed %s %s %s in scaffold " F_CID "\n",
              CI->id,
              (CI->scaffoldID == scaffold->id?"":"scaffoldID is wrong "),
              (!CI->flags.bits.isDead?"":"Dead "),
              (CI->flags.bits.isUnique?"":"Not Unique "),
              scaffold->id);
      DumpCIScaffold(GlobalData->stderrc,graph, scaffold, FALSE);
#ifdef STRICT_SCAFFOLD_CHECKING
      assert(0 /* Scaffold Sanity */);	  
#endif
    }
  }
  if(numElements != scaffold->info.Scaffold.numElements){
    fprintf(GlobalData->stderrc,
            "* numElements = %d scaffold says it has %d elements\n",
            numElements, scaffold->info.Scaffold.numElements);
    DumpCIScaffold(GlobalData->stderrc,graph, scaffold, FALSE);
    assert(0);
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

  currentMinPos = min( CI->offsetAEnd.mean, CI->offsetBEnd.mean);
  InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &CIs);
  while(NULL != (CI = NextCIScaffoldTIterator(&CIs)))
  {
    if( min( CI->offsetAEnd.mean, CI->offsetBEnd.mean) < currentMinPos)
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
      if( min( CI->offsetAEnd.mean, CI->offsetBEnd.mean) - currentMinPos < 1.0)
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
    currentMinPos = min( CI->offsetAEnd.mean, CI->offsetBEnd.mean);
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

#if 0


void CheckCITypes(ScaffoldGraphT *sgraph){
  int numCIs = GetNumChunkInstanceTs(sgraph->ChunkInstances);
  int i;
  for(i = 0; i < numCIs; i++){
    ChunkInstanceT *CI = GetChunkInstanceT(sgraph->ChunkInstances, i);
    if(CI->flags.bits.isUnique ){
      assert( (CI->type == DISCRIMINATORUNIQUECHUNK_CGW ||
               CI->type == UNIQUECHUNK_CGW));
      if(CI->scaffoldID == NULLINDEX){
	DumpSuspiciousCI(CI);
	//assert(0);
      }
    }else{
      assert( !(CI->type == DISCRIMINATORUNIQUECHUNK_CGW ||
                CI->type == UNIQUECHUNK_CGW));
      assert(CI->scaffoldID == NULLINDEX);
    }
  }
}


void CheckAllowedCITypes(ScaffoldGraphT *sgraph){
  int numCIs = GetNumChunkInstanceTs(sgraph->ChunkInstances);
  int i;
  for(i = 0; i < numCIs; i++){
    ChunkInstanceT *CI = GetChunkInstanceT(sgraph->ChunkInstances, i);
    if(CI->flags.bits.isUnique )
      assert( CI->type == DISCRIMINATORUNIQUECHUNK_CGW );
    else
      assert( CI->type == UNRESOLVEDCHUNK_CGW);
  }
}


void DumpSuspiciousCI(ChunkInstanceT *CI){
  CIEdgeT *nA = NULL, *nB = NULL;
  CDS_CID_t nidA = NULLINDEX, nidB = NULLINDEX;
  ChunkInstanceT *CIa = NULL, *CIb = NULL;
  fprintf(GlobalData->stderrc,
          "* Dump SuspiciousCI " F_CID " of type %d**\n", CI->id, CI->type);
  DumpChunkInstance(GlobalData->stderrc, ScaffoldGraph, CI,
                    FALSE, FALSE, FALSE, FALSE);
  fprintf(GlobalData->stderrc,
          "* numEssentialA:%d essentialA: " F_CID
          " numEssentialB:%d essentialB:" F_CID "\n",
	  CI->numEssentialA, CI->essentialEdgeA,
	  CI->numEssentialB, CI->essentialEdgeB);
  fprintf(GlobalData->stderrc,"* Essential Edges *\n");
  if(CI->essentialEdgeA != NULLINDEX){
    nA = GetCIEdgeT(ScaffoldGraph->CIEdges, CI->essentialEdgeA);
    nidA = (nA->idA == CI->id? nA->idB: nA->idA);
    CIa = GetChunkInstanceT (ScaffoldGraph->ChunkInstances, nidA);
    PrintCIEdgeT(GlobalData->stderrc, ScaffoldGraph, " ", nA , nidA);
  }
  if(CI->essentialEdgeB != NULLINDEX){
    nB = GetCIEdgeT(ScaffoldGraph->CIEdges, CI->essentialEdgeB);
    nidB = (nB->idA == CI->id? nB->idB: nB->idA);
    CIb = GetChunkInstanceT (ScaffoldGraph->ChunkInstances, nidB);
    PrintCIEdgeT(GlobalData->stderrc, ScaffoldGraph, " ", nB, nidB);
  }
  fprintf(GlobalData->stderrc,"* Essential Neighbors *\n");
  if(CIa){
    fprintf(GlobalData->stderrc,
            "* Chunk " F_CID " in Scaffold " F_CID " of type %d\n",
            CIa->id, CIa->scaffoldID, CIa->type);
    DumpChunkInstance(GlobalData->stderrc,ScaffoldGraph, CIa,
                      FALSE, FALSE, FALSE, FALSE);
  }    
  if(CIb){
    fprintf(GlobalData->stderrc,
            "* Chunk " F_CID " in Scaffold " F_CID " type %d\n",
            CIb->id, CIb->scaffoldID, CIb->type);
    DumpChunkInstance(GlobalData->stderrc, ScaffoldGraph, CIb,
                      FALSE, FALSE, FALSE, FALSE);
  }    
  fflush(GlobalData->stderrc);
}
#endif




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
    fprintf(GlobalData->logfp,"**** BEFORE repeat rez ****\n");
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
    
    fprintf(GlobalData->logfp,"**** AFTER repeat rez ****\n");
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
  fprintf(GlobalData->logfp,"** After BuildUniqueCIScaffolds **\n");
  DumpChunkInstances(GlobalData->logfp, ScaffoldGraph,
                     FALSE, TRUE, TRUE, FALSE);
  DumpCIScaffolds(GlobalData->logfp,ScaffoldGraph, TRUE);
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

  fprintf(GlobalData->stderrc,"* RebuildScaffolds save:%d  markShaky:%d\n",
	  GlobalData->saveCheckPoints, markShakyBifurcations);
  fflush(GlobalData->stderrc);

  if( markShakyBifurcations){
    char temp[2000];
    fprintf(GlobalData->timefp,"* Checkpoint %d, After Building Initial Unique CI Scaffolds and before Tidying up there are %d scaffolds\n",
	    ScaffoldGraph->checkPointIteration,
            (int) GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph));
    if(GlobalData->dumpScaffoldSnapshots){
      sprintf(temp,"Rebuild%d",ScaffoldGraph->checkPointIteration);
      DumpScaffoldSnapshot(temp);
    }
    CheckpointScaffoldGraph(ScaffoldGraph);
  }

  TidyUpScaffolds (ScaffoldGraph);

  if(GlobalData->saveCheckPoints){
    char temp[2000];

    fprintf(GlobalData->timefp,"* After RebuildScaffolds #%d there are %d scaffolds\n",
	    ScaffoldGraph->checkPointIteration,
            (int) GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph));
    if(GlobalData->dumpScaffoldSnapshots){
      sprintf(temp,"Rebuild%d",ScaffoldGraph->checkPointIteration);
      DumpScaffoldSnapshot(temp);
    }
    CheckpointScaffoldGraph(ScaffoldGraph);
  }

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
  
  if(GlobalData->debugLevel > 0){
    CheckAllContigFragments();
  }
  CheckAllTrustedEdges(ScaffoldGraph);
  
  /*** Build the scaffold edges from the raw edges in the graph of CIEdges **/
  BuildSEdges(ScaffoldGraph, FALSE);
    
#ifdef DEBUG_CGW
  DumpCIScaffolds(GlobalData->stderrc,ScaffoldGraph, TRUE);
#endif
  /* Merge the SEdges */
  MergeAllGraphEdges(ScaffoldGraph->ScaffoldGraph, TRUE);// Merge 'em

#ifdef DEBUG_CGW
  fprintf(GlobalData->logfp,"**** AFTER MERGESEdges ****\n");
  fflush(GlobalData->stderrc);
  DumpCIScaffolds(GlobalData->logfp,ScaffoldGraph, FALSE);
#endif

  ClearCacheSequenceDB(ScaffoldGraph->sequenceDB, FALSE);
  }



/***************************************************************************/
void BuildScaffoldsFromFirstPriniciples(ScaffoldGraphT *ScaffoldGraph,
                                        int skipInitialScaffolding){
  int changedByRepeatRez;

  GenerateScaffoldGraphStats("initial",1);

  if(skipInitialScaffolding){
    // This is a continuation of a previous partial result, or
    // a continuation of rocks
    fprintf(stderr,"* Tidying Up...\n");
    TidyUpScaffolds (ScaffoldGraph);
    fprintf(stderr,"* Done Tidying Up...\n");

    if(GlobalData->saveCheckPoints){
      
      fprintf(GlobalData->timefp,
              "* Checkpoint %d after tidying up %d scaffolds\n",
	      ScaffoldGraph->checkPointIteration,
              (int) GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph));
      if(GlobalData->dumpScaffoldSnapshots){
	char temp[2000];
	sprintf(temp,"Rebuild%d",ScaffoldGraph->checkPointIteration);
	DumpScaffoldSnapshot(temp);
      }
      CheckpointScaffoldGraph(ScaffoldGraph);
    }
  }else{
    // This includes CleanupScaffolds
    RebuildScaffolds(ScaffoldGraph, TRUE); // Transitive reduction of RezGraph followed by construction of SEdges
  }
  GenerateScaffoldGraphStats("initial",1);
  GeneratePlacedContigGraphStats("initial",1);
  GenerateLinkStats(ScaffoldGraph->ContigGraph,"initial",1);

#define  MAX_OUTER_REZ_ITERATIONS   10

  if(GlobalData->repeatRezLevel > 0){
    int iter = 0;

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
#if 0
	// We want tor recompute the contig coordinates
	// SAK **** Make sure our gap estimates are good before
        // SAK **** rebuilding scaffolds
	// SAK **** (This proved ineffectual)
	LeastSquaresGapEstimates(ScaffoldGraph, TRUE, FALSE, TRUE,
				 CHECK_CONNECTIVITY, FALSE);
#endif
	// Build Scaffolds of Discriminator Uniques
	fprintf (GlobalData -> timefp,
                 "Checkpoint %d written after Rocks %d\n",
                 ScaffoldGraph -> checkPointIteration, iter);
        
        // Transitive reduction of RezGraph followed by construction of SEdges
	RebuildScaffolds(ScaffoldGraph, FALSE); 

	iter ++;
      }
    
#ifdef DEBUG_BUCIS
      fprintf(GlobalData->logfp,"** After Gap Filling **\n");
      DumpChunkInstances(GlobalData->logfp, ScaffoldGraph,
                         FALSE, TRUE, TRUE, FALSE);
      DumpCIScaffolds(GlobalData->logfp,ScaffoldGraph, FALSE);
#endif   
    }while(changedByRepeatRez && iter < MAX_OUTER_REZ_ITERATIONS);
      
    fprintf(GlobalData->timefp,
            "\n\nCheckpoint %d written after Gap Filling\n",
            ScaffoldGraph->checkPointIteration);
    CheckpointScaffoldGraph(ScaffoldGraph);
    
  }

  
    

  // ++++++++++++++++++++++++++++++++
  // + CHANGED ++++++++++++++++++++++
  // ++++++++++++++++++++++++++++++++
  //     checkEdgeQuality(ScaffoldGraph->RezGraph,100,10,filePrefix);
  //  place_multiply_contained(ScaffoldGraph->ContigGraph,100,10,filePrefix;)

  

  fprintf(GlobalData->stderrc,"* After building scaffolds \n");
  fprintf(GlobalData->logfp,"* After building scaffolds \n");
  DumpCIScaffolds(GlobalData->logfp,ScaffoldGraph,FALSE);

}


