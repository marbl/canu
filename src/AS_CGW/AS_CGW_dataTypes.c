
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
static char CM_ID[] = "$Id: AS_CGW_dataTypes.c,v 1.4 2005-03-22 19:48:33 jason_miller Exp $";

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

#include "AS_global.h"
#include "AS_MSG_pmesg.h"
#include "AS_UTL_Var.h"
#include "UtilsREZ.h"
#include "AS_CGW_dataTypes.h"
#ifdef NEVER
#include "AS_CGW_histo.h"
#endif
#include "ScaffoldGraph_CGW.h"
#include "Globals_CGW.h"

#define INITIAL_NUM_DISTS 10

Global_CGW *CreateGlobal_CGW(void){
  Global_CGW *g;

  fprintf(stderr,"*CreateGlobal_CGW\n");
  g = (Global_CGW *)safe_calloc(1, sizeof(Global_CGW));
  AssertPtr(g);
  g->verbose = FALSE;
  g->outfp = NULL;
  g->logfp = NULL;
  g->gwlogfp = NULL;

#ifdef NEVER
  g->scaffold_repeat = create_extended_histogram(10000,500,0, TRUE, sizeof(ChunkAggregate), 
						 aggregateChunks,
						 printChunks,
						 printChunkAggregate);
  g->scaffold_unique = create_extended_histogram(10000,500,0, TRUE, sizeof(ChunkAggregate),
						 aggregateChunks,
						 printChunks,
						 printChunkAggregate);
#endif
  
  InitTimerT(&g->RecomputeOffsetsTimer);
  InitTimerT(&g->MergeScaffoldsTimer);
  InitTimerT(&g->BuildSEdgesTimer);
  InitTimerT(&g->InputTimer);
  InitTimerT(&g->OutputTimer);
  InitTimerT(&g->OverlapTimer);
  InitTimerT(&g->ConsistencyCheckTimer);
  InitTimerT(&g->WalkUpdateTimer);
  InitTimerT(&g->UpdateTimer);
  InitTimerT(&g->GapWalkerTimer);
  InitTimerT(&g->GapFillTimer);
  InitTimerT(&g->StoneThrowingTimer);
  InitTimerT(&g->BccTimer);
  InitTimerT(&g->ConsensusTimer);

  return(g);
}

#ifdef NEVER
void ResetHistograms_CGW(Global_CGW *g){

  free_histogram(g->scaffold_repeat);
  free_histogram(g->scaffold_unique);

  g->scaffold_repeat = create_extended_histogram(10000,500,0, TRUE, sizeof(ChunkAggregate), 
						 aggregateChunks,
						 printChunks,
						 printChunkAggregate);
  g->scaffold_unique = create_extended_histogram(10000,500,0, TRUE, sizeof(ChunkAggregate),
						 aggregateChunks,
						 printChunks,
						 printChunkAggregate);

}
#endif

void DeleteGlobal_CGW(Global_CGW *g){

  AssertPtr(g);

#ifdef NEVER
  free_histogram(g->scaffold_unique);
  free_histogram(g->scaffold_repeat);
#endif
  
#ifdef CREATE_CHUNK_GRAPH
  DeleteVA_DistInfoT(g->Dists);
  DeleteVA_ScaffoldPO_T(GlobalData->ScaffoldPOs);
  DeleteVA_ScaffoldPOChunkT(GlobalData->ScaffoldPOChunks);
  DeleteVA_ScaffoldPOLinkT(GlobalData->ScaffoldPOLinks);
  DeleteVA_ScaffoldPOOverlapT(GlobalData->ScaffoldPOOverlaps);
#endif
  free(g);
}





void ComputeIntervalLength(LengthT *result, 
			   LengthT *aEndA, LengthT *bEndA,
			   LengthT *aEndB, LengthT *bEndB){

  LengthT *leftMost, *rightMost;

  leftMost = aEndB;
  if( aEndA->mean < leftMost->mean)
    leftMost = aEndA;

  if( bEndA->mean < leftMost->mean)
    leftMost = bEndA;
  
  if( bEndB->mean < leftMost->mean)
    leftMost = bEndB;


  rightMost = aEndB;
  if( aEndA->mean > rightMost->mean)
    rightMost = aEndA;

  if( bEndA->mean > rightMost->mean)
    rightMost = bEndA;
  
  if( bEndB->mean > rightMost->mean)
    rightMost = bEndB;


  result->mean = rightMost->mean - leftMost->mean;
  result->variance = rightMost->variance - leftMost->variance;
  assert(result->variance > 0 && result->mean > 0);

}


// Do the arithmetic and stats on a pair of LengthTs
// If resulting variance is negative assert
void ComputeLength(LengthT *result, 
		   LengthT *length1, LengthT *length2){


  result->mean = length2->mean - length1->mean;
  result->variance = length2->variance - length1->variance;
  assert(result->variance >= 0 && result->mean > 0);  // temp change???

}

