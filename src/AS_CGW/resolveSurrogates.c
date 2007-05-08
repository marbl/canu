
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
static char CM_ID[] = "$Id: resolveSurrogates.c,v 1.16 2007-05-08 15:17:21 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_UTL_timer.h"
#include "AS_CGW_dataTypes.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"
#include "Globals_CGW.h"
#include "ScaffoldGraph_CGW.h"
#include "Output_CGW.h"
#include "GreedyOverlapREZ.h"
#include "CommonREZ.h"
#include "RepeatRez.h"
#include "FbacREZ.h"
#include "PublicAPI_CNS.h"
#include "AS_ALN_aligners.h"
#include "AS_ALN_forcns.h"
#include "fragmentPlacement.h"
#include "AS_UTL_Hash.h"


int my_getChunkInstanceID(ChunkInstanceT *chunk, int index)
{
  if (chunk->info.CI.numInstances == 0)  // chunk is not a surrogate
    {
      if (index == 0)  // just return chunk's id
        return(chunk->id);
      else
        return(-1);
    }
  else  // chunk is a surrogate
    {
      if (chunk->info.CI.numInstances == 1  && index == 0)
        return( chunk->info.CI.instances.in_line.instance1 );
      else if (chunk->info.CI.numInstances == 2 && (index == 0 || index == 1))
	{
	  if (index == 0)
            return( chunk->info.CI.instances.in_line.instance1);
	  else if (index == 1)
            return( chunk->info.CI.instances.in_line.instance2);
	}
      else if (index < chunk->info.CI.numInstances)
        return( * (int32 *) Getint32(chunk->info.CI.instances.va, index));
      else
        return(-1);
    }

  assert(0); //we should never get here
  return -1;
}





int main( int argc, char *argv[])
{
  int placeAllFragsInSinglePlacedSurros = 0;  /* if 1, aggressively place fragments in surrogates that are
						 only used once in the assembly; "aggressively" means place
						 all the fragments in the unitig, regardless of mate status,
					         alignment quality etc */
  float cutoffToInferSingleCopyStatus=1.0;
  Global_CGW *data;
  char *outputPath = NULL;
  int setGatekeeperStore = FALSE;
  int setPrefixName = FALSE;
  int setSingleSid = FALSE;
  CDS_CID_t singleSid = NULLINDEX;
  int ckptNum = NULLINDEX;
  ChunkInstanceT * parentChunk;
  LengthT gap;
  GraphNodeIterator CIGraphIterator;
  VA_TYPE(IntMultiPos) **impLists = NULL;
  int allocedImpLists = 100;
  int i;
  int numReallyPlaced=0;
  int totalNumParentFrags=0;

  if(impLists==NULL){
    impLists = (VA_TYPE(IntMultiPos)**) safe_malloc(allocedImpLists*sizeof(VA_TYPE(IntMultiPos)*));
    for(i=0;i<allocedImpLists;i++)
      impLists[i] = CreateVA_IntMultiPos(20);
  }

  GlobalData  = data = CreateGlobal_CGW();
  data->stderrc = stderr;
  data->timefp = stderr;

  setbuf(stdout,NULL);

  { /* Parse the argument list using "man 3 getopt". */ 
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg && ((ch = getopt(argc, argv,
				    "c:f:g:n:S:1")) != EOF)){
      switch(ch) {
        case 'c':
          {
            strcpy( data->File_Name_Prefix, argv[optind - 1]);
            setPrefixName = TRUE;		  
          }
          break;
        case 'g':
          {
            strcpy( data->Gatekeeper_Store_Name, argv[optind - 1]);
            setGatekeeperStore = TRUE;
          }
          break;	  
        case 'n':
          ckptNum = atoi(argv[optind - 1]);
          break;
        case 'S':
          cutoffToInferSingleCopyStatus=atof(argv[optind-1]);
          assert( cutoffToInferSingleCopyStatus >= 0.0 && cutoffToInferSingleCopyStatus <= 1.0 );
          break;
        case '1':
          placeAllFragsInSinglePlacedSurros = 1;
          break;
        default :
          fprintf(stderr,"Unrecognized option -%c",optopt);
          errflg++;
      }
    }

    if((setPrefixName == FALSE) || (setGatekeeperStore == 0))
      {
	fprintf(stderr, "usage: %s -g <gkp> -c <ckp> -n <num> opts\n",argv[0]);
        fprintf(stderr, "  -S x   place all frags in singly-placed surrogates if\n");
        fprintf(stderr, "         at least fraction x can be placed.\n");
	fprintf(stderr, "  -1     place all frags in singly-placed surrogates\n");
        fprintf(stderr, "         aggressively; equivalent to -S 0.0\n");
	exit(1);
      }
  }

  ScaffoldGraph = LoadScaffoldGraphFromCheckpoint( data->File_Name_Prefix, ckptNum, TRUE);
  GlobalData->aligner=Local_Overlap_AS_forCNS;

  // scan all the chunks

  InitGraphNodeIterator(&CIGraphIterator, ScaffoldGraph->CIGraph, GRAPH_NODE_DEFAULT);
  while ((parentChunk = NextGraphNodeIterator(&CIGraphIterator))!=NULL){

    int numFrgsToPlace=0;
    UIDHashTable_AS *fHash;
    int numInstances = parentChunk->info.CI.numInstances;
    int i,index, numFragmentsInParent;
    MultiAlignT *maParent;

    assert(numInstances <= allocedImpLists);
    // get the MultiAlign for this parentChunk

    // if numInstances >= 1 then it has a surrogate
    if(numInstances==0)continue;


    maParent = LoadMultiAlignTFromSequenceDB( ScaffoldGraph->sequenceDB, parentChunk->id, TRUE);

    AssertPtr (maParent);

    // count fragments and positions
    numFragmentsInParent = GetNumIntMultiPoss(maParent->f_list);

    totalNumParentFrags += numFragmentsInParent;

    for(i=0;i<numInstances;i++){
      ChunkInstanceT *candidateChunk;
      CDS_CID_t sid,ctgiid;

      index = my_getChunkInstanceID(parentChunk, i);
      candidateChunk = GetGraphNode(ScaffoldGraph->CIGraph, index);
      AssertPtr (candidateChunk);

      //  These were historically problems that were not asserts, but
      //  would just skip this instance.
      assert(parentChunk->type == UNRESOLVEDCHUNK_CGW);
      assert(candidateChunk->type == RESOLVEDREPEATCHUNK_CGW);
      assert(parentChunk != candidateChunk);

      if(candidateChunk->info.CI.baseID != parentChunk->id){
	if(candidateChunk==parentChunk){
	  fprintf(stderr,"resolveSurrogates: instance == parent for " F_CID " instance %d\n",
		  parentChunk->id,i);
	  continue;
	} else {
	  assert(candidateChunk->info.CI.baseID == parentChunk->id);
	}
      }

      sid=candidateChunk->scaffoldID;
      ctgiid = candidateChunk->info.CI.contigID;

      if(sid==NULLINDEX)continue;

      // now loop over fragments looking for those with mates in the scaffold
      {
	CGWFragIterator frags;
	CIFragT *nextfrg;
	InitCIFragTInChunkIterator(&frags,parentChunk,FALSE);
	while(NextCIFragTInChunkIterator(&frags, &nextfrg)){
	  CIFragT *mate;
	  ChunkInstanceT *mateChunk;
	  int fragIsGood = 0;

	  if( placeAllFragsInSinglePlacedSurros && numInstances == 1 ) {
	    fragIsGood=1;
	  } else {
	    if(matePlacedOnlyIn(nextfrg,sid,&mate,&mateChunk)){
	      assert(nextfrg->flags.bits.innieMate);
	      if(FragAndMateAreCompatible(nextfrg,candidateChunk,mate,mateChunk,AS_INNIE)){
		fragIsGood= 1;
	      }
	    }
	  }
	  if(fragIsGood){
	    // we're hot to trot ... now do something!
	    IntMultiPos imp;
	    imp.type = nextfrg->type;
	    imp.ident = nextfrg->iid;  
	    imp.position.bgn = nextfrg->offset5p.mean;
	    imp.position.end = nextfrg->offset3p.mean;
	    imp.contained = 0; /* this might be wrong! */
	    imp.delta_length=0;
	    imp.delta=NULL;
	    AppendVA_IntMultiPos(impLists[i],&imp);
	      
	    numFrgsToPlace++;
	  }
	}
	CleanupCIFragTInChunkIterator(&frags);
      }

    }  //  Over all instances


    if(numFrgsToPlace==0)continue;

    fHash = CreateUIDHashTable_AS(numFrgsToPlace);

    for(i=0;i<numInstances;i++){
      int j, numToPlace = GetNumIntMultiPoss(impLists[i]);
      for(j=0;j<numToPlace;j++){
	CDS_CID_t iid = GetIntMultiPos(impLists[i],j)->ident;
	int32 *count_so_far = LookupInUID2IIDHashTable_AS(fHash,(uint64)iid);
	if(count_so_far==NULL){
	  InsertInUID2IIDHashTable_AS(fHash,(uint64)iid,1);
	} else {
	  (*count_so_far)++;
	}
      }
    }

    for(i=0;i<numInstances;i++){

      int j, numToPlace = GetNumIntMultiPoss(impLists[i]);
      VA_TYPE(CDS_CID_t) *toplace;
      if(numToPlace==0)
        continue;

      toplace = CreateVA_CDS_CID_t(numToPlace);

      //  Build the list of fragments to place
      for(j=0;j<numToPlace;j++){
	CDS_CID_t iid = GetIntMultiPos(impLists[i],j)->ident;
	int32 *count_so_far = LookupInUID2IIDHashTable_AS(fHash,(uint64)iid);
	AssertPtr(count_so_far);
	assert(*count_so_far>0);
	if(*count_so_far>1){
	  continue;
	}
	AppendVA_CDS_CID_t(toplace,&iid);
      }


      // now, second-guess ourselves: if a sufficient fraction of reads can be placed with mates, then
      // place all fragments
      { 
	int reallyPlaced;
	float fractionPlaced;

	reallyPlaced = GetNumCDS_CID_ts(toplace);
	fractionPlaced = ((float) reallyPlaced)/ ((float) numFragmentsInParent);

	if( numInstances==1 && fractionPlaced>cutoffToInferSingleCopyStatus ){

	  CGWFragIterator frags;
	  CIFragT *nextfrg;
	  InitCIFragTInChunkIterator(&frags,parentChunk,FALSE);

	  ResetVA_CDS_CID_t(toplace);

	  while(NextCIFragTInChunkIterator(&frags, &nextfrg)){
	    AppendVA_CDS_CID_t(toplace,&(nextfrg->iid));  
	  }

	  CleanupCIFragTInChunkIterator(&frags);
	}

      }


      // now really do the placement
      ReallyAssignFragsToResolvedCI(ScaffoldGraph->CIGraph,
                                    parentChunk->id, 
                                    my_getChunkInstanceID(parentChunk,i),
                                    toplace);

      numReallyPlaced+=GetNumCDS_CID_ts(toplace);

      ResetVA_CDS_CID_t(toplace);

      ResetVA_IntMultiPos(impLists[i]);
    }

    DeleteUIDHashTable_AS(fHash);

  }

  for(i=0;i<allocedImpLists;i++){
    DeleteVA_IntMultiPos(impLists[i]);
  }

  fprintf(data->timefp,"Checkpoint %d written after resolveSurrogates\n",
	  ScaffoldGraph->checkPointIteration);
  CheckpointScaffoldGraph(ScaffoldGraph, -1);

  fprintf(data->stderrc,"Placed %d surrogate fragments out of %d surrogate fragments\n",
          numReallyPlaced, totalNumParentFrags);

  exit(0);
}
