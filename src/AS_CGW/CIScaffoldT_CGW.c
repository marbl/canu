
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
static char CM_ID[] = "$Id: CIScaffoldT_CGW.c,v 1.13 2006-11-06 22:53:42 brianwalenz Exp $";

#undef DEBUG
#undef DEBUG_INSERT
#undef DEBUG_DIAG
#undef DEBUG_SPLIT

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <fcntl.h>
#include <sys/types.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_UTL_interval.h"
#include "AS_CGW_dataTypes.h"
#include "Globals_CGW.h"
#include "DiagnosticsCGW.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"
#include "UnionFind_AS.h"
#include "UtilsREZ.h"
#include "ChiSquareTest_CGW.h"
#include "MultiAlignment_CNS.h"
#include "DataTypesREZ.h"
#include "CommonREZ.h"
#include "Stats_CGW.h"   // for collecting scaffold merging stats



void PrintCINodeFields(FILE * stream, NodeCGW_T * node)
{
  fprintf(stream, "\t\tcontigID:" F_CID "\n", node->info.CI.contigID);
  fprintf(stream, "\t\theadOfFragments:" F_CID "\n", node->info.CI.headOfFragments);
  fprintf(stream, "\t\tnumFragments:%d\n", node->info.CI.numFragments);
  fprintf(stream, "\t\tbranchPointA:" F_COORD "\n", node->info.CI.branchPointA);
  fprintf(stream, "\t\tbranchPointB:" F_COORD "\n", node->info.CI.branchPointB);
  fprintf(stream, "\t\tcoverageStat:%d\n", node->info.CI.coverageStat);
  fprintf(stream, "\t\tbaseID:" F_CID "\n", node->info.CI.baseID);
  fprintf(stream, "\t\tnumInstances:%d\n", node->info.CI.numInstances);
}

void PrintContigNodeFields(FILE * stream, NodeCGW_T * node)
{
  fprintf(stream, "\t\tAEndCI:" F_CID ", BEndCI:" F_CID ", numCI:%d\n",
          node->info.Contig.AEndCI, node->info.Contig.BEndCI,
          node->info.Contig.numCI);
  fprintf(stream, "\t\tbranchPointA:" F_COORD ", branchPointB:" F_COORD "\n",
          node->info.Contig.branchPointA, node->info.Contig.branchPointB);
}

void PrintScaffoldNodeFields(FILE * stream, NodeCGW_T * node)
{
  fprintf(stream, "\t\tAEndCI:" F_CID ", BEndCI:" F_CID ", numElements:%d\n",
          node->info.Scaffold.AEndCI, node->info.Scaffold.BEndCI,
          node->info.Scaffold.numElements);
  fprintf(stream, "\t\tleastSquareError:%.1f, numLeastSquareClones:%d\n",
          node->info.Scaffold.leastSquareError,
          node->info.Scaffold.numLeastSquareClones);
  fprintf(stream, "\t\tinternalEdges:%d, confirmedInternalEdges:%d\n",
          node->info.Scaffold.internalEdges,
          node->info.Scaffold.confirmedInternalEdges);
}

void PrintNodeFlagBits(FILE * stream, NodeCGW_T * node)
{
  fprintf(stream, "\t\tisUnique:%d, smoothSeenAlready:%d\n",
          node->flags.bits.isUnique, node->flags.bits.smoothSeenAlready);
  fprintf(stream, "\t\tcgbType:%d, isDead:%d, isFree:%d, containsCIs:%d\n",
          node->flags.bits.cgbType, node->flags.bits.isDead,
          node->flags.bits.isFree, node->flags.bits.containsCIs);
  fprintf(stream,
          "\t\ttandemOverlaps:%d, isCI:%d, isContig:%d, isScaffold:%d\n",
          node->flags.bits.tandemOverlaps, node->flags.bits.isCI,
          node->flags.bits.isContig, node->flags.bits.isScaffold);
  fprintf(stream, "\t\tisSurrogate:%d, beingContigged:%d, includesFinishedBacFragments:%d\n",
          node->flags.bits.isSurrogate, node->flags.bits.beingContigged,
          node->flags.bits.includesFinishedBacFragments);
  fprintf(stream, "\t\twalkedAlready:%d, walkedTooShort:%d, walkedTooLong:%d\n",
          node->flags.bits.walkedAlready, node->flags.bits.walkedTooShort,
          node->flags.bits.walkedTooLong);
  fprintf(stream, "\t\twalkMaxedOut:%d, walkedTrivial:%d\n",
          node->flags.bits.walkMaxedOut, node->flags.bits.walkedTrivial);
  fprintf(stream, "\t\tisStoneSurrogate:%d, isWalkSurrogate:%d\n",
          node->flags.bits.isStoneSurrogate, node->flags.bits.isWalkSurrogate);
  fprintf(stream, "\t\tfailedToContig:%d, isChaff:%d, isMisplaced:%d, isStone:%d\n",
          node->flags.bits.failedToContig, node->flags.bits.isChaff,
          node->flags.bits.isMisplaced, node->flags.bits.isStone);
  fprintf(stream, "\t\tisWalk:%d, isRock:%d, isPotentialRock:%d, isPotentialStone:%d\n",
          node->flags.bits.isWalk, node->flags.bits.isRock,
          node->flags.bits.isPotentialRock,
          node->flags.bits.isPotentialStone);
  fprintf(stream, "\tall:%d\n", node->flags.all);
}

void PrintNodeFields(FILE * stream, NodeCGW_T * node)
{
  fprintf(stream,"\ttype:%c, outputID:" F_CID ", scaffoldID:" F_CID ", prevScaffoldID:" F_CID "\n",
          node->type, node->outputID, node->scaffoldID, node->prevScaffoldID);
  fprintf(stream,"\tindexInScaffold:%d, smoothExpectedCID:" F_CID "\n",
          node->indexInScaffold, node->smoothExpectedCID);
  fprintf(stream, "\tnumEssentialA:%d, numEssentialB:%d\n",
          node->numEssentialA, node->numEssentialB);
  fprintf(stream, "\tessentialEdgeA:" F_CID ", essentialEdgeB:" F_CID "\n",
          node->essentialEdgeA, node->essentialEdgeB);
  fprintf(stream,"\tAEndNext:" F_CID ", BEndNext:" F_CID "\n",
          node->AEndNext, node->BEndNext);
  fprintf(stream, "\tbpLength:(%.1f,%.1f), offsetAEnd:(%.1f,%.1f), offsetBEnd:(%.1f,%.1f)\n",
          node->bpLength.mean, node->bpLength.variance,
          node->offsetAEnd.mean, node->offsetAEnd.variance,
          node->offsetBEnd.mean, node->offsetBEnd.variance);
  fprintf(stream, "\taEndCoord:" F_COORD ", bEndCoord:" F_COORD "\n",
          node->aEndCoord, node->bEndCoord);
  switch(node->type)
    {
      case DISCRIMINATORUNIQUECHUNK_CGW:
      case UNRESOLVEDCHUNK_CGW:
      case UNIQUECHUNK_CGW:
      case RESOLVEDREPEATCHUNK_CGW:
        //PrintUnitigNodeFields(stream, node);
        break;
      case CONTIG_CGW:
      case UNIQUECONTIG_CGW:
      case RESOLVEDCONTIG_CGW:
      case UNRESOLVEDCONTIG_CGW:
        //PrintContigNodeFields(stream, node);
        break;
      case REAL_SCAFFOLD:
      case OUTPUT_SCAFFOLD:
      case SCRATCH_SCAFFOLD:
        //PrintScaffoldNodeFields(stream, node);
        break;
    }
  PrintNodeFlagBits(stream, node);
  fprintf(stream, "\tedgeHead:" F_CID ", microhetScore:%.1f, setID:" F_CID "\n",
          node->edgeHead, node->microhetScore, node->setID);
}


void PrintCIScaffoldHeader(FILE *stream, ScaffoldGraphT *graph, CIScaffoldT *scaffold){
  fprintf(stream,"\n* CIScaffold " F_CID " numCI:%d (a:" F_CID " b:" F_CID ")  length: %d [" F_COORD "," F_COORD "]\n",
	  scaffold->id, scaffold->info.Scaffold.numElements, scaffold->info.Scaffold.AEndCI, scaffold->info.Scaffold.BEndCI,
	  (int)scaffold->bpLength.mean,
	  scaffold->aEndCoord, scaffold->bEndCoord);
  // PrintNodeFields(stream, scaffold);
}



void DumpCIScaffold(FILE *stream, ScaffoldGraphT *graph, CIScaffoldT *scaffold, int raw){
  ChunkInstanceT *CI;
  CIScaffoldTIterator CIs;
  SEdgeTIterator SEdges;
  SEdgeT *SEdge;

  PrintCIScaffoldHeader(stream, graph, scaffold);
  fprintf(stream, "> Includes CIs\n");
  InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &CIs);

  while((CI = NextCIScaffoldTIterator(&CIs)) != NULL){
    //      double ratio = 0.0;
    //      if(CI->offsetAEnd.mean > 0){
    //	ratio = ComputeFudgeVariance(CI->offsetAEnd.mean)/ CI->offsetAEnd.variance;
    //
    //      }
    fprintf(stream," \t %5d: CI %6" F_CIDP " sid:" F_CID " len:%6d,%6.2g ends:%6d,%6d var:%6.2g,%6.2g orient:%c [%8" F_COORDP ",%8" F_COORDP "]\n",
            CI->indexInScaffold,
            CI->id,
            CI->scaffoldID,
            (int)CI->bpLength.mean, CI->bpLength.variance,
            (int)CI->offsetAEnd.mean, (int)CI->offsetBEnd.mean,
            CI->offsetAEnd.variance, CI->offsetBEnd.variance, 
            //      ratio,
            (GetNodeOrient(CI) == A_B? 'A':'B'),
            CI->aEndCoord, CI->bEndCoord);
  }
  fprintf(stream, "> %s Edges A \n",  (raw?" R ":" M "));
  InitSEdgeTIterator(graph, scaffold->id, raw, FALSE, A_END, FALSE,  &SEdges);
  while((SEdge = NextSEdgeTIterator(&SEdges)) != NULL){
    PrintSEdgeT(stream, graph, " ", SEdge, scaffold->id);
  }


  fprintf(stream, "> %s Edges B\n",
	  (raw?" R ":" M "));
  InitSEdgeTIterator(graph, scaffold->id, raw, FALSE, B_END, FALSE,  &SEdges);
  while((SEdge = NextSEdgeTIterator(&SEdges)) != NULL){
    PrintSEdgeT(stream, graph, " ", SEdge, scaffold->id);
  }
}
/*****************************************************************/
void DumpACIScaffold(FILE *stream, ScaffoldGraphT *graph, CIScaffoldT *scaffold, int raw){
  ChunkInstanceT *CI;
  CIScaffoldTIterator CIs;

  PrintCIScaffoldHeader(stream, graph, scaffold);
  fprintf(stream, "> Includes CIs\n");
  InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &CIs);

  while((CI = NextCIScaffoldTIterator(&CIs)) != NULL){
    //      double ratio = 0.0;
    //      if(CI->offsetAEnd.mean > 0){
    //	ratio = ComputeFudgeVariance(CI->offsetAEnd.mean)/ CI->offsetAEnd.variance;
    //
    //      }
    DumpContig(stream, graph, CI,  raw);
  }
}
/***************************************************************************/
void DumpACIScaffoldNew(FILE *stream, ScaffoldGraphT *graph,
                        CIScaffoldT *scaffold, int raw){
  ChunkInstanceT *CI;
  CIScaffoldTIterator CIs;

  PrintCIScaffoldHeader(stream, graph, scaffold);
  fprintf(stream, "> Includes CIs\n");
  InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &CIs);

  while((CI = NextCIScaffoldTIterator(&CIs)) != NULL){
    //      double ratio = 0.0;
    //      if(CI->offsetAEnd.mean > 0){
    //	ratio = ComputeFudgeVariance(CI->offsetAEnd.mean)/ CI->offsetAEnd.variance;
    //
    //      }
    DumpContigInScfContext(stream, graph, CI,  raw);
  }
}
/***************************************************************************/
void DumpCIScaffolds(FILE *stream, ScaffoldGraphT *graph, int raw){
  CDS_CID_t sid;
  /* For each chunk */
  for(sid = 0; sid < GetNumCIScaffoldTs(graph->CIScaffolds); sid++){
    CIScaffoldT *scaffold = GetGraphNode(graph->ScaffoldGraph,sid);
    if(isDeadCIScaffoldT(scaffold))
      continue;
    DumpCIScaffold(stream, graph, scaffold, raw);
  }
}




void  MarkCIElementsForScaffoldMembership(ChunkInstanceT *chunkInstance,
                                          CDS_CID_t scaffoldID){
  if(chunkInstance->flags.bits.isContig){
    SetContigScaffoldIds(chunkInstance, scaffoldID);

  }else if(chunkInstance->flags.bits.isCI){
    SetCIScaffoldIds(chunkInstance, scaffoldID);
  }else assert(0);

}



/****************************************************************************/
// InsertCIInScaffold
// Insert chunk instance ci int scaffold sid at offset with orientation orient.
// offsetFromAEnd = offset of the end of the CI that is closest to the A end of the scaffold
// orient
int InsertCIInScaffold(ScaffoldGraphT *sgraph,
                       CDS_CID_t ci, CDS_CID_t sid, 
                       LengthT aEndOffset, LengthT bEndOffset,
                       int AEndToBend, int contigNow){

  CIScaffoldT *ciScaffold = GetGraphNode(sgraph->ScaffoldGraph, sid);
  ChunkInstanceT *chunkInstance = GetGraphNode(sgraph->RezGraph, ci);    
  int32 reversed;
  LengthT *maxOffset, *minOffset;

  assert(!chunkInstance->flags.bits.isDead);
  assert(!ciScaffold->flags.bits.isDead);

  // check for bad variances on ends
  if ( fabs (fabs(aEndOffset.variance - bEndOffset.variance) - chunkInstance->bpLength.variance) > 1.0 )
    {
      fprintf( stderr, "*** Variance Fixup Alert ***\n");

      fprintf( stderr, "chunkInstance->id: " F_CID "\n", chunkInstance->id);
      fprintf( stderr, "aEndOffset.mean: %f, bEndOffset.mean: %f\n", aEndOffset.mean, bEndOffset.mean);
      fprintf( stderr, "aEndOffset.variance: %f, bEndOffset.variance: %f\n", aEndOffset.variance, bEndOffset.variance);
      fprintf( stderr, "chunkInstance->bpLength.mean: %f\n", chunkInstance->bpLength.mean);
      fprintf( stderr, "chunkInstance->bpLength.variance: %f\n\n", chunkInstance->bpLength.variance);
	
      fprintf( stderr, "chunkInstance->offsetAEnd.mean: %f, offsetAEnd.variance: %f\n", 
               chunkInstance->offsetAEnd.mean, chunkInstance->offsetAEnd.variance);
      fprintf( stderr, "chunkInstance->offsetBEnd.mean: %f, offsetBEnd.variance: %f\n", 
               chunkInstance->offsetBEnd.mean, chunkInstance->offsetBEnd.variance);
      fprintf( stderr, "chunkInstance->bpLength.mean: %f, bpLength.variance: %f\n", 
               chunkInstance->bpLength.mean, chunkInstance->bpLength.variance);

      // now if we're too far off assert
      // assert ( fabs(aEndOffset.variance - bEndOffset.variance) >= 0.5 * chunkInstance->bpLength.variance);

      // we're in the ballpark, so just fixup the ends
      if (aEndOffset.mean < bEndOffset.mean)
        bEndOffset.variance = aEndOffset.variance + chunkInstance->bpLength.variance;
      else
        aEndOffset.variance = bEndOffset.variance + chunkInstance->bpLength.variance;
    }

#if defined( DEBUG) || defined(DEBUG_INSERT)
  fprintf(stderr,"* Insert ci:" F_CID " sid:" F_CID " contigNow:0x%x [%g,%g]\n",
    	  ci,sid,contigNow, aEndOffset.mean, bEndOffset.mean);
#endif

  {
    double diff = abs(aEndOffset.mean - bEndOffset.mean);
    if( diff > 1.2 * chunkInstance->bpLength.mean){
      fprintf(stderr,"***** SERIOUS: InsertCIInScaffold " F_CID " [%g,%g] = %g, length only %g\n",
	      ci, aEndOffset.mean, bEndOffset.mean, diff, chunkInstance->bpLength.mean);
    }
  }

  // This loop insures that we do not try to insert the same chunk
  // twice into a scaffold.
  {
    CIScaffoldTIterator CIs;
    ChunkInstanceT *CI;
    
    InitCIScaffoldTIterator(sgraph, ciScaffold, TRUE, FALSE, &CIs);
    while((CI = NextCIScaffoldTIterator(&CIs)) != NULL){
      if (CI->id == ci)  // chunk is already in scaffold due to being a member of a previous path
	return(0);
    }
  }

  chunkInstance->offsetAEnd = aEndOffset;
  chunkInstance->offsetBEnd = bEndOffset;

#ifndef DEBUG_INSERT
  if(GlobalData->debugLevel > 0){
#endif
    fprintf(stderr,"*&&& BEFORE inserting ci " F_CID " in scaffold " F_CID " at (%g,%g) (%g,%g)\n",
	    ci,sid, 
	    aEndOffset.mean, aEndOffset.variance,
	    bEndOffset.mean, bEndOffset.variance);
    DumpCIScaffold(stderr, sgraph, ciScaffold, FALSE);
#ifndef DEBUG_INSERT
  }
#endif

#if 0
  if(contigNow & CONTAINMENT_CONTIGGING){
    if( CheckForContainmentContigs(sgraph, ci, sid, aEndOffset, bEndOffset)){

#ifndef DEBUG_INSERT
      if(GlobalData->debugLevel > 0){
#endif
        fprintf(stderr,"* InsertCI " F_CID " InScaffold  " F_CID " returned after CheckForContigs\n", ci, sid);
        fprintf(stderr,"* AFTER\n");
        DumpCIScaffold(stderr,sgraph, ciScaffold, FALSE);
#ifndef DEBUG_INSERT
      }
#endif
      return 0; // we did it already
    }else{
#ifdef DEBUG_INSERT
      fprintf(stderr,"* InsertCI " F_CID " InScaffold  " F_CID " *Old Fashioned way *\n", ci, sid);
#endif
    }
  }
#endif
    
  if(contigNow & DOVETAIL_CONTIGGING){
    if(CheckForContigs(sgraph, ci, sid, aEndOffset, bEndOffset)){
#ifdef DEBUG_INSERT
      fprintf(stderr,"* InsertCI " F_CID " InScaffold  " F_CID " returned after CheckForContigs\n", ci, sid);
      fprintf(stderr,"* AFTER\n");
      DumpCIScaffold(stderr,sgraph, ciScaffold, FALSE);
#endif
      return 0; // we did it already
    }else{
#ifdef DEBUG_INSERT
      fprintf(stderr,"* InsertCI " F_CID " InScaffold  " F_CID " *Old Fashioned way *\n", ci, sid);
#endif
    }
  }

  assert(!chunkInstance->flags.bits.isDead);
  MarkCIElementsForScaffoldMembership(chunkInstance, sid);

  if((aEndOffset.mean > bEndOffset.mean)){ 
    reversed = TRUE;
    minOffset = &bEndOffset;
    maxOffset = &aEndOffset;
  }else{
    reversed = FALSE;
    minOffset = &aEndOffset;
    maxOffset = &bEndOffset;
  }    

  AssertPtr(ciScaffold);
  AssertPtr(chunkInstance);
  
  chunkInstance->scaffoldID = sid;
  chunkInstance->flags.bits.isChaff = FALSE;  // we need this one in the output



  ciScaffold->info.Scaffold.numElements++;

  /* Make sure that this chunk has simulation coordinates, even if the CI on the
     AEnd is a repeat */
  if(chunkInstance->flags.bits.cgbType == UU_CGBTYPE){
    if(ciScaffold->aEndCoord == -1){
      ciScaffold->aEndCoord = (reversed?chunkInstance->bEndCoord:chunkInstance->aEndCoord);
    }
    if(ciScaffold->bEndCoord == -1){
      ciScaffold->bEndCoord = (reversed?chunkInstance->aEndCoord:chunkInstance->bEndCoord);
    }
  }
#ifdef DEBUG_INSERT  
  fprintf(GlobalData->stderrc,"* Inserting cid " F_CID " into scaffold " F_CID " at offset %d, %d\n",
          ci, sid, (int) aEndOffset.mean, (int) bEndOffset.mean);
#endif

  if(ciScaffold->info.Scaffold.AEndCI == NULLINDEX){ /* Empty scaffold */

    ciScaffold->info.Scaffold.AEndCI = ciScaffold->info.Scaffold.BEndCI = ci;

    if(chunkInstance->flags.bits.cgbType == UU_CGBTYPE){
      ciScaffold->aEndCoord = (reversed?chunkInstance->bEndCoord:chunkInstance->aEndCoord);
      ciScaffold->bEndCoord = (reversed?chunkInstance->aEndCoord:chunkInstance->bEndCoord);
    }else{
      ciScaffold->aEndCoord = ciScaffold->bEndCoord = -1;
    }
    chunkInstance->AEndNext = chunkInstance->BEndNext = NULLINDEX;
    ciScaffold->bpLength = *maxOffset;

#ifdef DEBUG_CISCAFFOLD
    fprintf(GlobalData->stderrc,"* Inserted into empty scaffold aEndCoord=" F_COORD " bEndCoord=" F_COORD "\n",
	    ciScaffold->aEndCoord, ciScaffold->bEndCoord);
#endif
#ifdef DEBUG_INSERT
    fprintf(stderr,"* Inserted into empty scaffold aEndCoord=%d bEndCoord=%d\n",
	    ciScaffold->aEndCoord, ciScaffold->bEndCoord);
#endif
  }else{
    CIScaffoldTIterator CIs;
    ChunkInstanceT *CI;
    CDS_COORD_t chunkInstanceMin = (CDS_COORD_t) min(chunkInstance->offsetAEnd.mean,
                                                     chunkInstance->offsetBEnd.mean);
    InitCIScaffoldTIterator(sgraph, ciScaffold, AEndToBend, FALSE, &CIs);
    while((CI = NextCIScaffoldTIterator(&CIs)) != NULL){
      CDS_COORD_t CImin = (CDS_COORD_t) min(CI->offsetAEnd.mean,
                                            CI->offsetBEnd.mean);
      if(CImin > chunkInstanceMin) {

	/*          ------- chunkInstance
                    ------- CI           

                    (including containment of CI or positive gap) */

	ChunkInstanceT *prevCI;

	// WARNING: this condition is ok ONLY if AEndToBend == TRUE
	// When we traverse the list of CIs from the Bend, the condition 
	// (CImin > chunkInstanceMin) will be satisfied immediately and therefore
	// the chunk will end up in the wrong position (SteLo)
	assert(AEndToBend);

#ifdef DEBUG_INSERT
	fprintf(GlobalData->stderrc,"* Inserting before " F_CID "\n", CI->id);
#endif
	chunkInstance->BEndNext = CI->id;
	chunkInstance->AEndNext = CI->AEndNext;
	if (CI->AEndNext != NULLINDEX) {
	  prevCI = GetGraphNode(sgraph->RezGraph, CI->AEndNext);
	  AssertPtr(prevCI);
	  prevCI->BEndNext = ci;
	}
 
	CI->AEndNext = ci;

	if(CI->id == ciScaffold->info.Scaffold.AEndCI){
	  ciScaffold->info.Scaffold.AEndCI = ci;
	  if(chunkInstance->flags.bits.cgbType == UU_CGBTYPE){
	    ciScaffold->aEndCoord = (reversed?chunkInstance->bEndCoord:chunkInstance->aEndCoord);
	  }
#ifdef DEBUG_INSERT
	  fprintf(GlobalData->stderrc,"* Inserted at head\n");
#endif
	}else{
#ifdef DEBUG_INSERT
	  fprintf(GlobalData->stderrc,"* Inserted in middle\n");
#endif
	}


	break;
      }else if(CI->id == ciScaffold->info.Scaffold.BEndCI &&
	       CImin <= chunkInstanceMin ){ // append

	/*   -------- CI (last chunk currently in scaffold)
             -------- chunkInstance

             (includes containment or positive gap) */

	ciScaffold->info.Scaffold.BEndCI = ci;
	CI->BEndNext = ci;
	chunkInstance->AEndNext = CI->id;
	chunkInstance->BEndNext = NULLINDEX;
        if(chunkInstance->flags.bits.cgbType == UU_CGBTYPE){
          ciScaffold->bEndCoord = (reversed?chunkInstance->aEndCoord:chunkInstance->bEndCoord);
        }

	// Due to containments, the CI with the maximal mean does not
	// have the maximal variance.
	if(ciScaffold->bpLength.mean < maxOffset->mean)
	  ciScaffold->bpLength.mean = maxOffset->mean;
	if(ciScaffold->bpLength.variance < maxOffset->variance)
	  ciScaffold->bpLength.variance = maxOffset->variance;


#ifdef DEBUG_INSERT
	fprintf(GlobalData->stderrc,"* Inserted at end aEndCoord:" F_COORD " bEndCoord:" F_COORD " bpLength = %g\n",
		ciScaffold->aEndCoord, ciScaffold->bEndCoord, ciScaffold->bpLength.mean);
#endif

	break;
      }else if(CI->id == ciScaffold->info.Scaffold.AEndCI &&
	       CImin > chunkInstanceMin ){ // prepend

	// ALH: it does not appear we ever get here -- doesn't the first case cover this?
	// Let us find out ...  (9/21/04)
	assert(0);

	ciScaffold->info.Scaffold.AEndCI = ci;
	CI->AEndNext = ci;

	chunkInstance->BEndNext = CI->id;
	chunkInstance->AEndNext = NULLINDEX;
#ifdef DEBUG_INSERT
	fprintf(GlobalData->stderrc,"* CImin " F_COORD ", chunkInstanceMin " F_COORD "\n", CImin, chunkInstanceMin);
#endif
        if(chunkInstance->flags.bits.cgbType == UU_CGBTYPE){
          ciScaffold->aEndCoord = (reversed?chunkInstance->bEndCoord:chunkInstance->aEndCoord);
        }
	fprintf(GlobalData->stderrc,"* Inserted at head aEndCoord:" F_COORD " bEndCoord:" F_COORD "\n",
		ciScaffold->aEndCoord, ciScaffold->bEndCoord);
	break;

      } else { /* fall through to next existing chunk */
      }
    }


#ifndef DEBUG_INSERT
    fprintf(GlobalData->stderrc,"* Inserted CI " F_CID " in scaffold " F_CID "  Anext = " F_CID "   Bnext = " F_CID "\n",
	    chunkInstance->id, sid, chunkInstance->AEndNext, chunkInstance->BEndNext);
#endif

    // Due to containments, the CI with the maximal mean does not
    // have the maximal variance.
    if(ciScaffold->bpLength.mean < maxOffset->mean)
      ciScaffold->bpLength.mean = maxOffset->mean;
    if(ciScaffold->bpLength.variance < maxOffset->variance)
      ciScaffold->bpLength.variance = maxOffset->variance;

    if(GlobalData->debugLevel > 0){
      NodeCGW_T *previous = GetGraphNode(ScaffoldGraph->RezGraph, chunkInstance->AEndNext);

      if(previous && 
	 max(previous->offsetAEnd.variance, previous->offsetBEnd.variance) >
	 min(chunkInstance->offsetAEnd.variance, chunkInstance->offsetBEnd.variance)){
	fprintf(stderr,"**** VARIANCES ARE SCREWED UP ****\n");
      }

      fprintf(stderr,"* AFTER\n");
      DumpCIScaffold(stderr,sgraph, ciScaffold, FALSE);
    }

    return 0;
  }
  return 0;
}


/****************************************************************************/
// RemoveCIFromScaffold
// 
int RemoveCIFromScaffold(ScaffoldGraphT *sgraph, CIScaffoldT *ciScaffold,
                         ChunkInstanceT *CI, int adjustPositions){
  CDS_CID_t cid = CI->id;
  ChunkInstanceT *bnext = NULL, *anext = NULL;
  int middle = TRUE;
  int aend = FALSE;
  int bend = FALSE;
  LengthT base;
  LengthT maxoffset = (CI->offsetAEnd.mean < CI->offsetBEnd.mean? CI->offsetBEnd:CI->offsetAEnd);

#if 0
  if(ciScaffold->info.Scaffold.numElements < 3){;
  fprintf(GlobalData->stderrc,"* Removing CI " F_CID " from scaffold " F_CID ", elements left %d\nBEFORE:\n",
	  CI->id, ciScaffold->id, ciScaffold->info.Scaffold.numElements);

  DumpCIScaffold(stderr,sgraph, ciScaffold, FALSE);
  }
#endif

  assert(ciScaffold->info.Scaffold.AEndCI != NULLINDEX && ciScaffold->info.Scaffold.BEndCI != NULLINDEX);
  assert(ciScaffold && (CI->scaffoldID == ciScaffold->id));

  ciScaffold = GetGraphNode(sgraph->ScaffoldGraph, CI->scaffoldID);

  assert(ciScaffold && !isDeadCIScaffoldT(ciScaffold));

  if(CI->AEndNext != NULLINDEX)
    anext = GetGraphNode(sgraph->RezGraph, CI->AEndNext);

  if(CI->BEndNext != NULLINDEX)
    bnext = GetGraphNode(sgraph->RezGraph, CI->BEndNext);

#if 0
  fprintf(stderr,"* Predecessor is " F_CID " Successor is " F_CID "\n",
	  CI->AEndNext, CI->BEndNext);
#endif

  if(cid == ciScaffold->info.Scaffold.AEndCI){ // We're removing the Contig at the A-end of the scaffold
    ciScaffold->info.Scaffold.AEndCI = CI->BEndNext;
    if(bnext){
      aend = TRUE;
      if(bnext->offsetAEnd.mean < bnext->offsetBEnd.mean){
        base = bnext->offsetAEnd;
      }else{
        base = bnext->offsetBEnd;
      }
      assert(bnext->AEndNext == CI->id);
      bnext->AEndNext = NULLINDEX;
#if 0
      fprintf(stderr,"* bneighbor " F_CID " has AEndNext " F_CID " and BEndNext " F_CID "\n",
              bnext->id, bnext->AEndNext, bnext->BEndNext);
#endif
    }
    middle = FALSE;
  }

  if(cid == ciScaffold->info.Scaffold.BEndCI){
    ciScaffold->info.Scaffold.BEndCI = CI->AEndNext;
    if(anext){
      // NodeCGW_T *prevCI = GetGraphNode(ScaffoldGraph->ContigGraph, CI->AEndNext);
	  
      bend = TRUE;
      assert(anext->BEndNext == CI->id);
      anext->BEndNext = NULLINDEX;
      //  You can't simply look to your scaffold predecessor, 
      //  due to containments the scaffold length may not be determined by
      //  the previous CI in the scaffold...do a scan to determine max offset in scaffold
      //  and use this for scaffold length.

      SetCIScaffoldTLength(ScaffoldGraph, ciScaffold, TRUE);

      //  ciScaffold->bpLength.mean = max( prevCI->offsetAEnd.mean, prevCI->offsetBEnd.mean);
      //  ciScaffold->bpLength.variance = max( prevCI->offsetAEnd.variance, prevCI->offsetBEnd.variance);

	  
#if 0
      fprintf(stderr,"* aneighbor " F_CID " has AEndNext " F_CID " and BEndNext " F_CID "\n",
	      anext->id, anext->AEndNext, anext->BEndNext);
#endif
    }
    middle = FALSE;
  }

  if(middle){

    bnext->AEndNext = CI->AEndNext;
    anext->BEndNext = CI->BEndNext;

    if(abs(maxoffset.mean - ciScaffold->bpLength.mean) < 100){
      SetCIScaffoldTLength(ScaffoldGraph, ciScaffold, TRUE);
    }

#if 0
    fprintf(stderr,"* bneighbor " F_CID " has AEndNext " F_CID " and BEndNext " F_CID "\n",
	    bnext->id, bnext->AEndNext, bnext->BEndNext);
    fprintf(stderr,"* aneighbor " F_CID " has AEndNext " F_CID " and BEndNext " F_CID "\n",
	    anext->id, anext->AEndNext, anext->BEndNext);
#endif
  }  

  CI->scaffoldID = NULLINDEX;
  CI->AEndNext = CI->BEndNext = NULLINDEX;

  // If we are deleting from the a-end, renormalize all scaffold coordinates to start from 0
  if(aend && adjustPositions){
    base.mean = - base.mean;
    base.variance = - base.variance;
    // AddDelta adjusts contig positions and bpLength
    AddDeltaToScaffoldOffsets(ScaffoldGraph, ciScaffold->id, bnext->id, TRUE, FALSE, base);

    fprintf(stderr,"* After RemoveCIFromScaffold ci:" F_CID " scaffold:" F_CID " removing from a-end base = (%g,%g), bpLength = (%g,%g)\n",
	    CI->id,
	    ciScaffold->id,
	    base.mean, base.variance,
	    ciScaffold->bpLength.mean, ciScaffold->bpLength.variance);
    //    CheckCIScaffoldTLength(sgraph, ciScaffold);
  }
  ciScaffold->info.Scaffold.numElements--;
#if 0
  fprintf(stderr,"* Removing CI " F_CID " froms scaffold " F_CID "\nAFTER:\n",
	  CI->id, ciScaffold->id);
  DumpCIScaffold(stderr,sgraph, ciScaffold, FALSE);
#endif
  if(ciScaffold->info.Scaffold.numElements == 0){
    ciScaffold->bpLength.mean = 0.0;
    ciScaffold->bpLength.variance = 0.0;
  }
  return FALSE;
}




/****************************************************************************/
void FindScaffoldComponents(ScaffoldGraphT *graph, int findPaths){
  int32 numScaffolds = GetNumCIScaffoldTs(graph->CIScaffolds);
  UFDataT *UFData;
  int scaffoldMap[numScaffolds];
  CDS_CID_t i;
  int set = 0;
  int setA, setB;
  int numComponents;

  UFData = UFCreateSets(graph->numLiveScaffolds);
  //  scaffoldMap = (int *)malloc(sizeof(int) * numScaffolds);

#ifdef DEBUG_MERGE
  fprintf(stderr,"* ScaffoldGraph has %d scaffolds\n", graph->numLiveScaffolds);
#endif

  for(i = 0; i < numScaffolds; i++){
    CIScaffoldT *scaffold = GetGraphNode(graph->ScaffoldGraph, i);
    UFSetT *scaffoldSet = UFGetSet(UFData, set);
    SEdgeTIterator SEdges;
    SEdgeT *sedge, *sedgeA, *sedgeB;
    int32 aEndEdges, bEndEdges;
    CDS_CID_t essentialEdgeA = NULLINDEX, essentialEdgeB = NULLINDEX;

    if(isDeadCIScaffoldT(scaffold) ||
       scaffold->type != REAL_SCAFFOLD){
      scaffoldMap[i] = NULLINDEX;
      continue;
    }
    AssertPtr(scaffoldSet);
    scaffoldMap[i] = set++;
    scaffoldSet->data = (void *)scaffold;

    // TEMPORARY-- WE ARE FILTERING EDGES WITH MEAN < 0
    // 

    aEndEdges = 0;
    InitSEdgeTIterator(graph, i, FALSE, TRUE, A_END, FALSE, &SEdges);

    // This merging code is simpleminded.  We want to prevent merges
    // that imply an overlap between the two scaffolds, since we don't
    // deal with them in an elegant manner.

#define MIN_MERGE_GAP_SIZE 0

    while((sedge = NextSEdgeTIterator(&SEdges)) != NULL){
      //  At some point in the past, someone tried subtracting 3.0 *
      //  sqrt(sedge->distance.variance) from the mean, and then
      //  ifdef'd it out.
      if(sedge->distance.mean > MIN_MERGE_GAP_SIZE){
        aEndEdges++;
        essentialEdgeA = GetVAIndex_SEdgeT(graph->SEdges, sedge);
        sedgeA = sedge;
      }
    }

    bEndEdges = 0;
    InitSEdgeTIterator(graph, i, FALSE, TRUE, B_END, FALSE, &SEdges);
    while((sedge = NextSEdgeTIterator(&SEdges)) != NULL){
      //  At some point in the past, someone tried subtracting 3.0 *
      //  sqrt(sedge->distance.variance) from the mean, and then
      //  ifdef'd it out.
      if(sedge->distance.mean > MIN_MERGE_GAP_SIZE){
        bEndEdges++;
        essentialEdgeB = GetVAIndex_SEdgeT(graph->SEdges, sedge);
        sedgeB = sedge;
      }
    }

    scaffold->numEssentialA = aEndEdges;
    scaffold->numEssentialB = bEndEdges;
    
    scaffold->essentialEdgeA = (aEndEdges == 1?essentialEdgeA:NULLINDEX);
    scaffold->essentialEdgeB = (bEndEdges == 1?essentialEdgeB:NULLINDEX);

#ifdef DEBUG_MERGE
    fprintf(stderr,"* Scaffold " F_CID " has (%d,%d) edges\n", i,aEndEdges, bEndEdges);
    if(aEndEdges)PrintSEdgeT(stderr,graph, " ", sedgeA, i);
    if(bEndEdges)PrintSEdgeT(stderr,graph, " ", sedgeB, i);
#endif
    
  }

  /* Now connect two scaffolds if they each have only a single confirmed edge on the
     side where their connecting edge attaches */

  for(i = 0; i < numScaffolds; i++){
    CIScaffoldT *scaffold = GetGraphNode(graph->ScaffoldGraph, i);
    SEdgeT *sedge;
    SEdgeTIterator SEdges;

    if(isDeadCIScaffoldT(scaffold) ||
       scaffold->type != REAL_SCAFFOLD)
      continue;

    if(findPaths){
      //    fprintf(stderr,"* Scaffold " F_CID " has %d a edges and %d b edges\n",
      //	    i,scaffold->numEssentialA, scaffold->numEssentialB);

      if(scaffold->numEssentialA == 1){
        CIScaffoldT *otherScaffold;
        int end;
        int edgeCount;

        sedge = GetGraphEdge(graph->ScaffoldGraph, scaffold->essentialEdgeA);

        if(sedge->idA != i)
          continue;

        if(sedge->distance.mean < -50000)
          continue;

#ifndef DONT_FIX_BUG
        if(sedge->flags.bits.inAssembly)
          continue;
#endif

        otherScaffold = GetGraphNode(graph->ScaffoldGraph,sedge->idB);
        if(sedge->orient == BA_AB){
          end = A_END;
          edgeCount = otherScaffold->numEssentialA;
        }else{
          end = B_END;
          edgeCount = otherScaffold->numEssentialB;
        }
#ifdef DEBUG_MERGE
        fprintf(stderr,"* Considering edge (" F_CID "," F_CID ") between sets (%d, %d)\n",
                sedge->idA, sedge->idB,
                scaffoldMap[sedge->idA],
                scaffoldMap[sedge->idB]);
#endif	      
        if(edgeCount == 1){
#ifndef DONT_FIX_BUG
          sedge->flags.bits.inAssembly = TRUE;
#endif
          setA = scaffoldMap[sedge->idA];
          setB = scaffoldMap[sedge->idB];
          //	fprintf(stderr,"* Union (%d,%d)\n", setA, setB);
          UFUnion(UFData, setA, setB);
        }else{
          /*
            fprintf(stderr,"* Scaffold " F_CID " has edgeCount = %d on %s end\n",
            sedge->idB,edgeCount,(end == A_END?"A":"B"));
          */
        }
      }else{
        /* We've ruled this edge out...so remove reference to it */
        scaffold->essentialEdgeA = NULLINDEX;
        scaffold->numEssentialA = 0;
      }

      if(scaffold->numEssentialB == 1){
        CIScaffoldT *otherScaffold;
        int end;
        int edgeCount = 0;

        sedge = GetGraphEdge(graph->ScaffoldGraph, scaffold->essentialEdgeB);

        if(sedge->idA != i)
          continue;

        if(sedge->distance.mean < -50000)
          continue;

#ifdef DONT_FIX_BUG
        if(sedge->flags.bits.inAssembly)
          continue;
#endif
#ifdef DEBUG_MERGE
        fprintf(stderr,"* Considering edge (" F_CID "," F_CID ") between sets (%d, %d)\n",
                sedge->idA, sedge->idB,
                scaffoldMap[sedge->idA],
                scaffoldMap[sedge->idB]);
#endif
        otherScaffold = GetGraphNode(graph->ScaffoldGraph,sedge->idB);
        if(sedge->orient == AB_AB){
          end = A_END;
          edgeCount = otherScaffold->numEssentialA;
        }else{
          end = B_END;
          edgeCount = otherScaffold->numEssentialB;
        }

        if(edgeCount == 1){
#ifdef DONT_FIX_BUG
          sedge->flags.bits.inAssembly = TRUE;
#endif
          setA = scaffoldMap[sedge->idA];
          setB = scaffoldMap[sedge->idB];
          //	fprintf(stderr,"* Union (%d,%d)\n", setA, setB);
          UFUnion(UFData, setA, setB);
        }else{
          /*
            fprintf(stderr,"* Scaffold " F_CID " has edgeCount = %d on %s end\n",
            sedge->idB, edgeCount,(end == A_END?"A":"B"));
          */
        }
      }else{
        /* We've ruled this edge out...so remove reference to it */
        scaffold->essentialEdgeB = NULLINDEX;
        scaffold->numEssentialB = 0;
      }


    }else{

      InitSEdgeTIterator(graph, i, FALSE, TRUE, ALL_END, FALSE, &SEdges);
      while((sedge = NextSEdgeTIterator(&SEdges)) != NULL){
        if(sedge->idA != i)
          continue;
        setA = scaffoldMap[sedge->idA];
        setB = scaffoldMap[sedge->idB];
        //	fprintf(stderr,"* Union (%d,%d)\n", setA, setB);
        UFUnion(UFData, setA, setB);
      }
    }
  }
  numComponents = UFRenumberSets(UFData);
  fprintf(stderr," Scaffold Graph has %d subPaths\n", numComponents);

  for(set = 0; set < UFData->numSets; set++){
    UFSetT *scaffoldSet = UFGetSet(UFData, set);
    CIScaffoldT *scaffold = (CIScaffoldT *)scaffoldSet->data;

    scaffold->setID = scaffoldSet->component;
#ifdef DEBUG_DIAG
    fprintf(GlobalData->stderrc,"* Scaffold " F_CID " has component " F_CID "\n",
            scaffold->id, scaffold->setID);
#endif
  }

  UFFreeSets(UFData);
  //  free(scaffoldMap);

}



/***************************************************************************/
int IsScaffoldInternallyConnected(ScaffoldGraphT *sgraph,
                                  CIScaffoldT *scaffold, int32 edgeTypes) {
  //
  // returns the number of connected components of the <scaffold>
  // NOTE: it considers ONLY trusted edges
  // Will modify the setId field of the NodeCGW_T structure to reflect
  // which component a node belongs to.
  //
  UFDataT
    * UFData = UFCreateSets(scaffold->info.Scaffold.numElements);
  CIEdgeT
    * edge;
  ChunkInstanceT
    * chunk;
  GraphEdgeIterator   edges;
  CIScaffoldTIterator CIs;
  int set = 0;
  int numComponents;

  assert(UFData != NULL);
  assert(scaffold != NULL);
  assert(sgraph != NULL);

  //
  // make a set for each vertex
  //
  InitCIScaffoldTIterator(sgraph, scaffold, TRUE,  FALSE, &CIs);
  while ((chunk = NextCIScaffoldTIterator(&CIs)) != NULL) {
    //
    // create a set
    //
    UFSetT
      * chunkSet = UFGetSet(UFData, set);
    //
    // map the set to a chunk
    //
    chunkSet->data = (void *)chunk;
    //
    // map the chunkId to setId
    //
    chunk->setID = set++;
  }

  //
  // now do the unions: iterate over all trusted/raw edges
  //
  InitCIScaffoldTIterator(sgraph, scaffold, TRUE,
                          FALSE, &CIs);
  while ((chunk = NextCIScaffoldTIterator(&CIs)) != NULL) {
    assert(chunk->setID >= 0);
    InitGraphEdgeIterator(sgraph->RezGraph, chunk->id, 
                          ALL_END, edgeTypes, // ALL_TRUSTED_EDGES, 
                          GRAPH_EDGE_DEFAULT, //GRAPH_EDGE_CONFIRMED_ONLY,
                          &edges);
    while ((edge = NextGraphEdgeIterator(&edges)) != NULL) {
      //
      // get the other end
      //
      ChunkInstanceT
        * otherChunk = GetGraphNode(sgraph->RezGraph,
                                    (chunk->id == edge->idA) ?
                                    edge->idB : edge->idA);
      int32 weight = edge->edgesContributing - (isOverlapEdge(edge));
      assert(otherChunk != NULL);
      
      // See each edge only once
      if(chunk->id != edge->idA)
        continue;

      if(edge->flags.bits.isBridge){
        fprintf(stderr,"* WARNING: chunk " F_CID " weight = %d bridge edge\n",
                chunk->id, weight);
        PrintGraphEdge(stderr, ScaffoldGraph->ContigGraph,
                       "Bridge ", edge, chunk->id);
#ifdef DEBUG
        EdgeCGW_T *e;
        GraphEdgeIterator Edges;
        InitGraphEdgeIterator(sgraph->ContigGraph,chunk->id,ALL_END,
                              ALL_TRUSTED_EDGES,GRAPH_EDGE_DEFAULT,&Edges);
        fprintf(stderr,"Edges out from " F_CID ":\n",chunk->id); 
        while(NULL!= (e = NextGraphEdgeIterator(&Edges)))
          PrintGraphEdge(stderr, ScaffoldGraph->ContigGraph,
                         "DEBUG Bridge ",e, chunk->id);	
#endif
      }
      if(isSingletonOverlapEdge(edge) ||
         (weight == 1 && edge->flags.bits.isBridge))
        continue;

      //
      // if the other end is not in this scaffold
      // ignore it
      //
      if (chunk->scaffoldID != otherChunk->scaffoldID)
        continue;

      //
      // do union
      //
      UFUnion(UFData, chunk->setID, otherChunk->setID);
    }
  }
  
  //
  // clean it up and return the # of components
  //
  numComponents = UFRenumberSets(UFData);

  //
  // renumber the NodeCGW_T setID to reflect component membership
  //
  set = 0;
  InitCIScaffoldTIterator(sgraph, scaffold, TRUE,  FALSE, &CIs);
  while ((chunk = NextCIScaffoldTIterator(&CIs)) != NULL) {
    //
    // create a set
    //
    UFSetT
      * chunkSet = UFGetSet(UFData, set);
    assert(chunkSet->data == (void *)chunk);
    //
    // map the chunkId to setId
    //
    chunk->setID = chunkSet->component;
    set++;
  }
  UFFreeSets(UFData);

  //fprintf(stderr, "IsScaffoldInternallyConnected() sid:"F_CID" %d\n", scaffold->id, numComponents);

  return numComponents;
}

/****************************************************************************/
int IsScaffoldInternallyConnectedCheck(ScaffoldGraphT *sgraph,
                                       CIScaffoldT *scaffold,
                                       int32 edgeTypes,
                                       CDS_CID_t ignoredChunkID)
{
  //
  // returns the number of connected components of the <scaffold>
  // NOTE: it considers ONLY trusted edges
  // Will modify the setId field of the NodeCGW_T structure to reflect
  // which component a node belongs to.
  //
  UFDataT  * UFData = UFCreateSets(scaffold->info.Scaffold.numElements);
  CIEdgeT  * edge;
  ChunkInstanceT  * chunk;
  GraphEdgeIterator   edges;
  CIScaffoldTIterator CIs;
  int32 numComponents;
  int set = 0;

  assert(UFData != NULL);
  assert(scaffold != NULL);
  assert(sgraph != NULL);

  //
  // make a set for each vertex
  //
  InitCIScaffoldTIterator(sgraph, scaffold, TRUE,  FALSE, &CIs);
  while ((chunk = NextCIScaffoldTIterator(&CIs)) != NULL) {
    //
    // create a set
    //
    UFSetT
      * chunkSet = UFGetSet(UFData, set);

    if (chunk->id == ignoredChunkID)
      continue;
	
    //
    // map the set to a chunk
    //
    chunkSet->data = (void *)chunk;
    //
    // map the chunkId to setId
    //
    chunk->setID = set++;
  }

  //
  // now do the unions: iterate over all trusted/raw edges
  //
  InitCIScaffoldTIterator(sgraph, scaffold, TRUE,
                          FALSE, &CIs);
  while ((chunk = NextCIScaffoldTIterator(&CIs)) != NULL) {
    if (chunk->id == ignoredChunkID)
      continue;
    assert(chunk->setID >= 0);
    InitGraphEdgeIterator(sgraph->RezGraph, chunk->id, 
                          ALL_END, edgeTypes, // ALL_TRUSTED_EDGES, 
                          GRAPH_EDGE_DEFAULT, //GRAPH_EDGE_CONFIRMED_ONLY,
                          &edges);
    while ((edge = NextGraphEdgeIterator(&edges)) != NULL) {
      //
      // get the other end
      //
      ChunkInstanceT
        * otherChunk = GetGraphNode(sgraph->RezGraph,
                                    (chunk->id == edge->idA) ?
                                    edge->idB : edge->idA);
      int32 weight = edge->edgesContributing - (isOverlapEdge(edge));
      assert(otherChunk != NULL);
      
      // See each edge only once
      if(chunk->id != edge->idA)
        continue;

      if(edge->flags.bits.isBridge){
        fprintf(stderr,"* WARNING: chunk " F_CID " weight = %d bridge edge\n", chunk->id, weight);
        PrintGraphEdge(stderr, ScaffoldGraph->ContigGraph, "Bridge ",edge, chunk->id);	
      }
      if(isSingletonOverlapEdge(edge) ||
         (weight == 1 && edge->flags.bits.isBridge))
        continue;

      //
      // if the other end is not in this scaffold
      // ignore it
      //
      if (chunk->scaffoldID != otherChunk->scaffoldID)
        continue;

      //
      // do union
      //
      UFUnion(UFData, chunk->setID, otherChunk->setID);
    }
  }

  //
  // clean it up and return the # of components
  //
  numComponents = UFRenumberSets(UFData);

  //
  // renumber the NodeCGW_T setID to reflect component membership
  //
  set = 0;
  InitCIScaffoldTIterator(sgraph, scaffold, TRUE,  FALSE, &CIs);
  while ((chunk = NextCIScaffoldTIterator(&CIs)) != NULL) {
    //
    // create a set
    //
    UFSetT
      * chunkSet = UFGetSet(UFData, set);
    assert(chunkSet->data == (void *)chunk);
    //
    // map the chunkId to setId
    //
    chunk->setID = chunkSet->component;
    set++;
  }
  UFFreeSets(UFData);
  return numComponents;
}


void
killScaffoldIfOnlySurrogate(CDS_CID_t scaffoldID) {
  CIScaffoldT     *scaffold  = GetGraphNode(ScaffoldGraph->ScaffoldGraph, scaffoldID);
  ContigT         *contig    = NULL;
  ChunkInstanceT  *chunk     = NULL;
  ChunkInstanceT  *basechunk = NULL;

  if (scaffold->flags.bits.isDead)
    return;

  if (scaffold->info.Scaffold.numElements > 1)
    return;

  contig = GetGraphNode(ScaffoldGraph->ContigGraph, scaffold->info.Scaffold.AEndCI);

  if (contig->info.Contig.numCI == 1) {
    chunk = GetGraphNode(ScaffoldGraph->CIGraph, contig->info.Contig.AEndCI);

    if (chunk->flags.bits.isStoneSurrogate) {
      basechunk = GetGraphNode(ScaffoldGraph->CIGraph, chunk->info.CI.baseID);

      fprintf(stderr, "WARNING!  scaffold %d contians just a surrogate (contig=%d chunk=%d base=%d)!\n",
              scaffold->id,
              contig->id,
              chunk->id,
              basechunk->id);
          
      //  See ya!

      //  Clean up the base unitig
      //
      if (basechunk->info.CI.numInstances == 1) {
        basechunk->flags.bits.isChaff = FALSE;
        if (basechunk->info.CI.numFragments == 1)
          basechunk->flags.bits.isChaff = TRUE;
        basechunk->info.CI.instances.in_line.instance1 = -1;
        basechunk->info.CI.instances.in_line.instance2 = -1;
        basechunk->info.CI.numInstances = 0;
      } else if (basechunk->info.CI.numInstances == 2) {
        if (basechunk->info.CI.instances.in_line.instance1 == chunk->id)
          basechunk->info.CI.instances.in_line.instance1 = basechunk->info.CI.instances.in_line.instance2;
        basechunk->info.CI.instances.in_line.instance2 = -1;
        basechunk->info.CI.numInstances = 1;
      } else if (basechunk->info.CI.numInstances == 3) {
        CDS_CID_t  a = *GetCDS_CID_t(basechunk->info.CI.instances.va, 0);
        CDS_CID_t  b = *GetCDS_CID_t(basechunk->info.CI.instances.va, 1);
        CDS_CID_t  c = *GetCDS_CID_t(basechunk->info.CI.instances.va, 2);

        if (a == chunk->id)
          a = c;
        if (b == chunk->id)
          b = c;

        DeleteVA_CDS_CID_t(basechunk->info.CI.instances.va);
        basechunk->info.CI.instances.in_line.instance1 = a;
        basechunk->info.CI.instances.in_line.instance2 = b;
        basechunk->info.CI.numInstances = 2;
      } else {
        //  Find which one is this chunk, move the last one over it.
        int  index = 0;

        for (index=0; index<basechunk->info.CI.numInstances; index++)
          if (*GetCDS_CID_t(basechunk->info.CI.instances.va, index) == chunk->id)
            SetCDS_CID_t(basechunk->info.CI.instances.va,
                         index,
                         GetCDS_CID_t(basechunk->info.CI.instances.va, basechunk->info.CI.numInstances-1));
        basechunk->info.CI.numInstances--;
      }

      //  Kill the unitig, contig, scaffold edges and scaffold.
      DeleteGraphNode(ScaffoldGraph->CIGraph, chunk);
      DeleteGraphNode(ScaffoldGraph->ContigGraph, contig);

      DeleteScaffoldEdgesForScaffold(ScaffoldGraph, scaffold);

      scaffold->flags.bits.isDead         = TRUE; 
      scaffold->info.Scaffold.AEndCI      = NULLINDEX;
      scaffold->info.Scaffold.BEndCI      = NULLINDEX;
      scaffold->info.Scaffold.numElements = 0;
      scaffold->bpLength.mean             = 0.0;
      scaffold->bpLength.variance         = 0.0;
    }
  }
}


/***********************************************************************/
int32 CheckScaffoldConnectivityAndSplit(ScaffoldGraphT *graph, CDS_CID_t scaffoldID, int32 edgeTypes, int verbose){
  CIScaffoldT  *scaffold      = GetCIScaffoldT(graph->CIScaffolds, scaffoldID);
  int           numComponents = IsScaffoldInternallyConnected(graph, scaffold, edgeTypes);
  int32         numNodes      = scaffold->info.Scaffold.numElements;

  // Expected case, Scaffold is connected
  if(numComponents > 1){
    CDS_CID_t nodes[numNodes];
    int inode;

    // IsScaffoldInternalyConnected does a connected component analysis, marking the contigs with their component number
    // the following code leverages this marking to break up the scaffold.

    int component;
    int  *nodesEnd;
    NodeCGW_T *thisNode;
    CIScaffoldTIterator scaffoldNodes;

    fprintf(stderr, "WARNING! Scaffold " F_CID " is not connected has %d components\n", scaffoldID, numComponents);
    fprintf(stderr, "Splitting into scaffolds: (search for \"Splitting "F_CID" into scaffold\" to get the new scaffolds)\n", scaffoldID);

    if(verbose)
      DumpACIScaffold(stderr,graph, scaffold, FALSE);

#ifdef DEBUG_SPLIT
    fprintf(stderr,"Prior to split ...");
    DumpACIScaffoldNew(stderr,graph, scaffold, FALSE);
#endif 

    nodesEnd = nodes + numNodes;
    InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &scaffoldNodes);
	
    inode = 0;
    while((thisNode = NextCIScaffoldTIterator(&scaffoldNodes)) != NULL){
      assert(inode < numNodes);
      nodes[inode++] = thisNode->id;
    }

    // For each component, create a scaffold and insert the relevant
    // contigs
    for(component = 0; component < numComponents; component++){
      LengthT NullLength = {0.0, 0.0};
      LengthT firstOffset;
      int seenFirstOffset;
      CIScaffoldT CIScaffold;
      CDS_CID_t newScaffoldID;

      InitializeScaffold(&CIScaffold, REAL_SCAFFOLD);
      CIScaffold.info.Scaffold.AEndCI = NULLINDEX;
      CIScaffold.info.Scaffold.BEndCI = NULLINDEX;
      CIScaffold.info.Scaffold.numElements = 0;
      CIScaffold.edgeHead = NULLINDEX;
      CIScaffold.bpLength = NullLength;
      newScaffoldID = CIScaffold.id = GetNumGraphNodes(graph->ScaffoldGraph);
      CIScaffold.flags.bits.isDead = FALSE;
      CIScaffold.aEndCoord = CIScaffold.bEndCoord = -1;
      CIScaffold.numEssentialA = CIScaffold.numEssentialB = 0;
      CIScaffold.essentialEdgeB = CIScaffold.essentialEdgeA = NULLINDEX;
      AppendGraphNode(graph->ScaffoldGraph, &CIScaffold);

      scaffold = GetCIScaffoldT(graph->CIScaffolds, scaffoldID);

      for(inode = 0, seenFirstOffset = FALSE; inode < numNodes; inode++){
        NodeCGW_T *thisNode = GetGraphNode(graph->RezGraph, nodes[inode]);
        if(thisNode->setID == component){
          LengthT offsetAEnd, offsetBEnd;
          if(!seenFirstOffset){
            if(GetNodeOrient(thisNode) == A_B){
              firstOffset = thisNode->offsetAEnd;
            }else{
              firstOffset = thisNode->offsetBEnd;
            }
            seenFirstOffset = TRUE;
          }
          offsetAEnd.mean     = thisNode->offsetAEnd.mean     - firstOffset.mean;
          offsetAEnd.variance = thisNode->offsetAEnd.variance - firstOffset.variance;
          offsetBEnd.mean     = thisNode->offsetBEnd.mean     - firstOffset.mean;
          offsetBEnd.variance = thisNode->offsetBEnd.variance - firstOffset.variance;

          RemoveCIFromScaffold(graph, scaffold, thisNode, FALSE);
          InsertCIInScaffold(graph, thisNode->id, newScaffoldID, offsetAEnd, offsetBEnd, TRUE, FALSE);
        }
      }
      assert((GetGraphNode(graph->ScaffoldGraph, newScaffoldID))->info.Scaffold.numElements > 0);

      fprintf(stderr, "Splitting "F_CID" into scaffold "F_CID"\n", scaffoldID, newScaffoldID);

      //  Make sure that our new scaffold contains more than just a single surrogate.
      //
      killScaffoldIfOnlySurrogate(newScaffoldID);

#ifdef DEBUG_SPLIT
      fprintf(stderr,"... post split ...");
      DumpACIScaffoldNew(stderr,graph,
                         GetGraphNode(graph->ScaffoldGraph,newScaffoldID), 
                         TRUE);
#endif 

    }

    // Delete any remaining edges
    DeleteScaffoldEdgesForScaffold(graph, scaffold);

    // Mark the old scaffold dead
    scaffold->flags.bits.isDead         = TRUE;
    scaffold->info.Scaffold.AEndCI      = NULLINDEX;
    scaffold->info.Scaffold.BEndCI      = NULLINDEX;
    scaffold->info.Scaffold.numElements = 0;
    scaffold->bpLength.mean             = 0.0;
    scaffold->bpLength.variance         = 0.0;
  }
  return numComponents;
}
    

/*****************************************************************************/

void CheckTrustedEdges(ScaffoldGraphT * sgraph,  CDS_CID_t cid) {
  //
  // iterates over all trusted edges of CI cid and check if
  // the other end is in the same scaffold
  //
  GraphEdgeIterator  edges;
  CIEdgeT * edge;
  CDS_CID_t next;
  ChunkInstanceT * next_chunk;
  ChunkInstanceT * this_chunk = GetGraphNode(sgraph->RezGraph, cid);
  CDS_CID_t sid = this_chunk->scaffoldID;

  InitGraphEdgeIterator(sgraph->RezGraph, cid, 
                        ALL_END, ALL_TRUSTED_EDGES, 
                        GRAPH_EDGE_DEFAULT,
                        &edges);
  while((edge = NextGraphEdgeIterator(&edges)) != NULL){
    assert(edge != NULL);
    //
    // get the other end
    //
    if (cid == edge->idA)
      next = edge->idB;
    else
      //      continue;  // avoid double checking of (i,j) and (j,i)
      next = edge->idA;

    next_chunk = GetGraphNode(ScaffoldGraph->RezGraph, next);
    assert(next_chunk != NULL);

    if (next_chunk->scaffoldID != sid)
#if 1
      fprintf(stderr,"-=> BAD edge id:" F_CID " " F_CID "(" F_CID ")->" F_CID "(" F_CID ") (weight %d, status %d)\n",
              (CDS_CID_t) GetVAIndex_CIEdgeT(sgraph->RezGraph->edges, edge),
              cid,
              sid,
              next,
              next_chunk->scaffoldID,
              edge->edgesContributing,
              edge->flags.bits.edgeStatus);
#endif
  }
}


/*****************************************************************************/

void CheckAllTrustedEdges(ScaffoldGraphT * sgraph){
  GraphNodeIterator nodes;
  ChunkInstanceT *contig;
  
  InitGraphNodeIterator(&nodes, sgraph->RezGraph, GRAPH_NODE_DEFAULT);
  while((contig = NextGraphNodeIterator(&nodes)) != NULL){

    if(contig->scaffoldID == NULLINDEX)
      continue;

    CheckTrustedEdges(sgraph, contig->id);

  }
}

/*****************************************************************************/

int CheckAllEdges(ScaffoldGraphT * sgraph,  CDS_CID_t sid, CDS_CID_t cid) {
  //
  // iterates over all edges of CI cid and check if
  // the other end is in the same scaffold
  // returns the number of edges (of any type except UNTRUSTED_EDGE_STATUS)
  // to chunks in another scaffold
  //
  GraphEdgeIterator edges;
  EdgeCGW_T * edge;
  CDS_CID_t next;
  NodeCGW_T * next_chunk;
  int out_of_sid_links = 0;

  InitGraphEdgeIterator(sgraph->RezGraph, cid, ALL_END, ALL_EDGES, GRAPH_EDGE_DEFAULT, &edges);
  while((edge = NextGraphEdgeIterator(&edges)) != NULL){
    assert(edge != NULL);

    if (!(edge->flags.bits.edgeStatus & (TRUSTED_EDGE_STATUS  || TENTATIVE_TRUSTED_EDGE_STATUS)) )
      continue;
    //
    // get the other end
    //
    if (cid == edge->idA)
      next = edge->idB;
    else
      next = edge->idA;

    next_chunk = GetGraphNode(ScaffoldGraph->RezGraph, next);
    assert(next_chunk != NULL);

    if ((next_chunk->scaffoldID != sid) && (next_chunk->scaffoldID != -1)) {
      out_of_sid_links++;
      /*** mjf ***/
      fprintf(stderr, "in CheckAllEdges -=> BAD edge id:" F_CID " " F_CID "(" F_CID ")->" F_CID "(" F_CID ") (weight %d, status %d)\n",
              (CDS_CID_t) GetVAIndex_EdgeCGW_T(sgraph->RezGraph->edges, edge),
              cid,
              sid,
              next,
              next_chunk->scaffoldID,
              edge->edgesContributing,
              edge->flags.bits.edgeStatus);
    }
  }
  return out_of_sid_links;
}


/*****************************************************************************/


void CheckCIScaffoldTs(ScaffoldGraphT *sgraph){
  GraphNodeIterator scaffolds;
  CIScaffoldT *scaffold;

  InitGraphNodeIterator(&scaffolds, sgraph->ScaffoldGraph, GRAPH_NODE_DEFAULT);
  while((scaffold = NextGraphNodeIterator(&scaffolds)) != NULL){
    if(scaffold->type != REAL_SCAFFOLD)
      continue;
    assert(!scaffold->flags.bits.isDead);
    CheckCIScaffoldT(sgraph, scaffold);
  }
}


void CheckCIScaffoldTLengths(ScaffoldGraphT *sgraph){
  GraphNodeIterator scaffolds;
  CIScaffoldT *scaffold;

  InitGraphNodeIterator(&scaffolds, sgraph->ScaffoldGraph, GRAPH_NODE_DEFAULT);
  while((scaffold = NextGraphNodeIterator(&scaffolds)) != NULL){
    if(scaffold->type != REAL_SCAFFOLD)
      continue;
    assert(!scaffold->flags.bits.isDead);
    CheckCIScaffoldTLength(sgraph, scaffold);
  }
}


void  SetCIScaffoldTLength(ScaffoldGraphT *sgraph, CIScaffoldT *scaffold, int32 verbose){
  CIScaffoldTIterator CIs;
  ChunkInstanceT *chunk;
  LengthT maxOffset = {0.0,0.0};

  InitCIScaffoldTIterator(sgraph, scaffold, TRUE, FALSE, &CIs);
  while((chunk = NextCIScaffoldTIterator(&CIs)) != NULL){

    if(chunk->offsetAEnd.mean > maxOffset.mean){
      maxOffset.mean = chunk->offsetAEnd.mean;
    }
    if(chunk->offsetBEnd.mean > maxOffset.mean){
      maxOffset.mean = chunk->offsetBEnd.mean;
    }
    if(chunk->offsetAEnd.variance > maxOffset.variance){
      maxOffset.variance = chunk->offsetAEnd.variance;
    }
    if(chunk->offsetBEnd.variance > maxOffset.variance){
      maxOffset.variance = chunk->offsetBEnd.variance;
    }
  }

  if(verbose && (scaffold->bpLength.mean != maxOffset.mean ||
                 scaffold->bpLength.variance != maxOffset.variance)){
    fprintf(GlobalData->stderrc, "SetCIScaffoldTLength adjusted length of scaffold " F_CID " from (%g,%g) to (%g,%g)\n",
            scaffold->id, scaffold->bpLength.mean, scaffold->bpLength.variance,
            maxOffset.mean, maxOffset.variance);
    if (verbose > 2)
      DumpCIScaffold(GlobalData->stderrc, sgraph, scaffold, FALSE);
  }

  scaffold->bpLength = maxOffset;
}

void SetCIScaffoldTLengths(ScaffoldGraphT *sgraph, int verbose){
  GraphNodeIterator scaffolds;
  CIScaffoldT *scaffold;

  InitGraphNodeIterator(&scaffolds, sgraph->ScaffoldGraph, GRAPH_NODE_DEFAULT);
  while((scaffold = NextGraphNodeIterator(&scaffolds)) != NULL){
    if(scaffold->type != REAL_SCAFFOLD)
      continue;
    assert(!scaffold->flags.bits.isDead);
    SetCIScaffoldTLength(sgraph, scaffold, verbose);
  }

}


//  bpw -- for each chunk in a scaffold, checks that the chunk
//  position is sane, printing warnings if not.  Possibly extends the
//  scaffold size to fit all chunks.
//
void CheckCIScaffoldTLength(ScaffoldGraphT *sgraph, CIScaffoldT *scaffold){
  CIScaffoldTIterator CIs;
  ChunkInstanceT *chunk, *prevChunk;
  CDS_CID_t sid;
  // double mean, variance;
  int cnt = 0;

  sid = scaffold->id;

  InitCIScaffoldTIterator(sgraph, scaffold, TRUE,
                          FALSE, &CIs);
  prevChunk = NULL;
  while((chunk = NextCIScaffoldTIterator(&CIs)) != NULL){
    assert(chunk->scaffoldID == sid);
       
    if(cnt++ == 0){
      if((chunk->offsetAEnd.mean > 0.1 &&
          chunk->offsetBEnd.mean > 0.1) ||
         (chunk->offsetAEnd.variance <0.0 ||
          chunk->offsetBEnd.variance < 0.0) ){
        fprintf(stderr,"*# First Chunk " F_CID " in scaffold " F_CID " is not positioned at scaffold start, but at (%g,%g)...fixing\n",
                chunk->id, scaffold->id, chunk->offsetAEnd.mean, chunk->offsetBEnd.mean);
#ifdef STRICT_SCAFFOLD_CHECKING
        assert(0);
#endif
        chunk->offsetAEnd.mean = chunk->offsetAEnd.variance = 0.0;
        chunk->offsetBEnd = chunk->bpLength;
      }
    }
    /*	fprintf(stderr,"* Chunk " F_CID " [%g,%g] mean = %g\n",
        chunk->id, 
        chunk->offsetAEnd.mean,
        chunk->offsetBEnd.mean,
        mean); 
    */
    if(chunk->offsetAEnd.mean > scaffold->bpLength.mean ||
       chunk->offsetBEnd.mean > scaffold->bpLength.mean){

      fprintf(stderr,"*# Chunk " F_CID " has end point [%g,%g] outside of scaffold length %g\n",
              chunk->id,
              chunk->offsetAEnd.mean,
              chunk->offsetBEnd.mean,
              scaffold->bpLength.mean);

    }

    if(chunk->offsetAEnd.mean < 0 || chunk->offsetBEnd.mean < 0 ||
       chunk->offsetAEnd.variance < 0 || chunk->offsetBEnd.variance < 0){
      fprintf(stderr,"*# Chunk " F_CID " in scaffold " F_CID " at odd position: (%g,%g) (%g,%g)\n",
              chunk->id, scaffold->id,
              chunk->offsetAEnd.mean, chunk->offsetAEnd.variance,
              chunk->offsetBEnd.mean, chunk->offsetBEnd.variance);
#ifdef STRICT_SCAFFOLD_CHECKING	   
      assert(0);
#endif
      chunk->offsetAEnd.mean = chunk->offsetAEnd.variance = 0.0;
      chunk->offsetBEnd = chunk->bpLength;
    }

    if(chunk->offsetAEnd.mean > scaffold->bpLength.mean || chunk->offsetAEnd.variance > scaffold->bpLength.variance){
      fprintf(stderr,"*# A Scaffold " F_CID " bpLength TOO Short (%g < %g: %g)...variances (%g, %g: %g)corrected\n",
              scaffold->id, 
              scaffold->bpLength.mean, chunk->offsetAEnd.mean, 
              fabs(scaffold->bpLength.mean- chunk->offsetAEnd.mean), 
              scaffold->bpLength.variance, chunk->offsetAEnd.variance,
              fabs(scaffold->bpLength.variance- chunk->offsetAEnd.variance));
#ifdef STRICT_SCAFFOLD_CHECKING
      assert(0);
#endif
      scaffold->bpLength.mean = chunk->offsetAEnd.mean;
      scaffold->bpLength.variance = chunk->offsetAEnd.variance;
    }
    if(chunk->offsetBEnd.mean > scaffold->bpLength.mean || chunk->offsetBEnd.variance > scaffold->bpLength.variance){
      fprintf(stderr,"*# B Scaffold " F_CID " bpLength TOO Short (%g < %g:%g)...variances (%g, %g:%g) corrected\n",
              scaffold->id, 
              scaffold->bpLength.mean, chunk->offsetBEnd.mean, 
              fabs(scaffold->bpLength.mean- chunk->offsetAEnd.mean), 
              scaffold->bpLength.variance, chunk->offsetBEnd.variance,
              fabs(scaffold->bpLength.variance- chunk->offsetBEnd.variance));
#ifdef STRICT_SCAFFOLD_CHECKING
      assert(0);
#endif
      scaffold->bpLength.mean = chunk->offsetBEnd.mean;
      scaffold->bpLength.variance = chunk->offsetBEnd.variance;
    }

  }
}




/* Check means and variance in a single */
void CheckCIScaffoldT(ScaffoldGraphT *sgraph, CIScaffoldT *scaffold){
  CIScaffoldTIterator CIs;
  ChunkInstanceT *chunk, *prevChunk;
  CDS_CID_t sid;
  double mean, variance;
  int cgwError;
  int iterationError;
  int iteration = 0;
  double LSE;
  double improvement = 1.0;

  static  VA_TYPE(PtrT) *chunksToBeRemoved = NULL;

  if(chunksToBeRemoved == NULL){
    chunksToBeRemoved = CreateVA_PtrT(2048);
  }else{
    ResetVA_PtrT(chunksToBeRemoved);
  }


  sid = scaffold->id;
  LSE = scaffold->info.Scaffold.leastSquareError;


  ScaffoldSanity(scaffold, sgraph);

  if(scaffold->info.Scaffold.numElements == 1){
#ifdef DEBUG
    fprintf(GlobalData->stderrc, " Early end to CheckCIScaffoldT() for scaffold " F_CID ": numElements ==1\n",scaffold->id);
#endif
    return;
  }
   
  cgwError = TRUE;
  iteration = 0;
  while(cgwError && improvement > 0.005 && iteration++ < 20){
    int status;
    int numChecked = 0;
    cgwError = FALSE;
    mean = - DBL_MAX;
    ResetVA_PtrT(chunksToBeRemoved);
    InitCIScaffoldTIterator(sgraph, scaffold, TRUE,
                            FALSE, &CIs);
    prevChunk = NULL;
    while((chunk = NextCIScaffoldTIterator(&CIs)) != NULL){
      assert(chunk->scaffoldID == sid);

      /*	fprintf(stderr,"* Chunk " F_CID " [%g,%g] mean = %g\n",
		chunk->id, 
		chunk->offsetAEnd.mean,
		chunk->offsetBEnd.mean,
		mean); 
      */
      if(chunk->offsetAEnd.mean > scaffold->bpLength.mean ||
         chunk->offsetBEnd.mean > scaffold->bpLength.mean){

        fprintf(stderr,"* Chunk " F_CID " has end point [%g,%g] outside of scaffold length %g\n",
                chunk->id,
                chunk->offsetAEnd.mean,
                chunk->offsetBEnd.mean,
                scaffold->bpLength.mean);

      }

      /* Check for Chimeric Scaffolds */
      if(prevChunk && 
         prevChunk->aEndCoord >= 0 &&
         chunk->aEndCoord >=0){
        CDS_COORD_t calcDiff = (CDS_COORD_t) fabs(chunk->offsetAEnd.mean - prevChunk->offsetAEnd.mean );
        CDS_COORD_t realDiff = abs(chunk->aEndCoord - prevChunk->aEndCoord);

        if(realDiff > 250000 || ((realDiff > calcDiff) && (realDiff - calcDiff > 100000))){
          fprintf(stderr,"*** Scaffold " F_CID " is CHIMERIC at point between CIs " F_CID " and " F_CID "\n",
                  sid, prevChunk->id, chunk->id);
          fprintf(stderr,"*** CI " F_CID " [" F_COORD "," F_COORD "]  CI " F_CID " [" F_COORD "," F_COORD "] scaffold gap is " F_COORD "\n",
                  prevChunk->id, 
                  prevChunk->aEndCoord, prevChunk->bEndCoord,
                  chunk->id,
                  chunk->aEndCoord, chunk->bEndCoord,
                  calcDiff);
        }
      }
      prevChunk = chunk;

      if(chunk->offsetAEnd.mean < 0 || chunk->offsetBEnd.mean < 0 ||
         chunk->offsetAEnd.variance < 0 || chunk->offsetBEnd.variance < 0){
        fprintf(stderr,"* Chunk " F_CID " in scaffold " F_CID " at odd position: (%g,%g) (%g,%g) fixing\n",
                chunk->id, scaffold->id,
                chunk->offsetAEnd.mean, chunk->offsetAEnd.variance,
                chunk->offsetBEnd.mean, chunk->offsetBEnd.variance);
        chunk->offsetAEnd.mean = chunk->offsetAEnd.variance = 0.0;
        chunk->offsetBEnd = chunk->bpLength;

      }

      if(chunk->offsetAEnd.mean > scaffold->bpLength.mean || chunk->offsetAEnd.variance > scaffold->bpLength.variance){
        fprintf(stderr,"* A Scaffold " F_CID " bpLength TOO Short (%g < %g)...variances (%g, %g)corrected\n",
                scaffold->id, scaffold->bpLength.mean, chunk->offsetAEnd.mean, scaffold->bpLength.variance, chunk->offsetAEnd.variance);
        scaffold->bpLength.mean = chunk->offsetAEnd.mean;
        scaffold->bpLength.variance = chunk->offsetAEnd.variance;
      }
      if(chunk->offsetBEnd.mean > scaffold->bpLength.mean || chunk->offsetBEnd.variance > scaffold->bpLength.variance){
        fprintf(stderr,"* B Scaffold " F_CID " bpLength TOO Short (%g < %g)...variances (%g, %g) corrected\n",
                scaffold->id, scaffold->bpLength.mean, chunk->offsetBEnd.mean, scaffold->bpLength.variance, chunk->offsetBEnd.variance);
        scaffold->bpLength.mean = chunk->offsetBEnd.mean;
        scaffold->bpLength.variance = chunk->offsetBEnd.variance;
      }

      if(chunk->offsetAEnd.mean < mean ||
         chunk->offsetBEnd.mean < mean){
        fprintf(GlobalData->stderrc,"* Screwed up scaffold " F_CID ": Chunk " F_CID " has bad mean\n",
                sid, chunk->id);

        AppendPtrT(chunksToBeRemoved, (const void *) &chunk);

      }else{
        // Only update the mean if we are moving along as expected
        mean = min(chunk->offsetAEnd.mean, chunk->offsetBEnd.mean);
      }
      numChecked++;
    }
    cgwError = GetNumPtrTs(chunksToBeRemoved);
    if(cgwError){
      int i;
      fprintf(stderr,"* Scaffold " F_CID " iteration %d  %d bad CIs out of %d LSE:%g improvement:%g\n",
              sid, iteration, cgwError, numChecked, LSE, improvement);
      for(i = 0; i < GetNumPtrTs(chunksToBeRemoved); i++){

        chunk = *(ChunkInstanceT **)GetPtrT(chunksToBeRemoved, i);
        fprintf(stderr,"* Reinserting chunk " F_CID " in scaffold " F_CID "\n",
                chunk->id,sid);
        RemoveCIFromScaffold(sgraph, scaffold, chunk, FALSE);
        InsertCIInScaffold(sgraph, chunk->id, sid,
                           chunk->offsetAEnd, chunk->offsetBEnd, TRUE, FALSE);
      }

#if 0
      // special one-time hack added for mouse_20010307 run
      MarkInternalEdgeStatus(sgraph, scaffold, PAIRWISECHI2THRESHOLD_CGW,
                             SLOPPY_EDGE_VARIANCE_THRESHHOLD, TRUE, TRUE, 0, TRUE);	   
#endif
  

      status = RecomputeOffsetsInScaffold(sgraph, scaffold, TRUE, TRUE /* was FALSE*/,FALSE);
      if (status != RECOMPUTE_OK) {
        fprintf(stderr, "RecomputeOffsetsInScaffold failed (%d) for scaffold " F_CID " in CheckScaffolds\n",
                status, sid);
        break;
      }
      // This is how much the LSE improved
      improvement = 1.0; // first time through
      if(LSE > 0.0){
        improvement = (LSE - scaffold->info.Scaffold.leastSquareError)/LSE;
      }
      fprintf(stderr,"* improvement = %g LSE = %g\n", improvement, LSE);
      LSE = scaffold->info.Scaffold.leastSquareError;
    }
  }

  if(iteration > 20){
    fprintf(stderr,"* Took %d iterations to patch up scaffold " F_CID "\n",
            iteration, sid);
    cgwError = TRUE;
  }
  variance = -1.0;
  mean = -1.0;
  InitCIScaffoldTIterator(sgraph, scaffold, TRUE,
                          FALSE, &CIs);
  while((chunk = NextCIScaffoldTIterator(&CIs)) != NULL){
    
    if((chunk->offsetAEnd.mean > mean && chunk->offsetAEnd.variance < variance) ||
       (chunk->offsetBEnd.mean > mean && chunk->offsetBEnd.variance < variance)){
      cgwError = TRUE;
      fprintf(GlobalData->stderrc,"* Screwed up scaffold " F_CID ": Chunk " F_CID " has bad variance\n",
              sid, chunk->id);
    }

    if( mean < chunk->offsetAEnd.mean || mean < chunk->offsetBEnd.mean)
      {
        mean = max(chunk->offsetAEnd.mean, chunk->offsetBEnd.mean);
        variance = max(chunk->offsetAEnd.variance, chunk->offsetBEnd.variance);
      }
  }
  if(cgwError){
    fprintf(GlobalData->stderrc,"* Screwed up scaffold " F_CID "\n", sid);
    DumpCIScaffold(stderr,sgraph, scaffold, FALSE);
    iterationError = TRUE;
  }
  ScaffoldSanity(scaffold, sgraph);


}

void FixupLengthsScaffoldTs(ScaffoldGraphT *sgraph){
  GraphNodeIterator scaffolds;
  CIScaffoldT *scaffold;

  InitGraphNodeIterator(&scaffolds, sgraph->ScaffoldGraph, GRAPH_NODE_DEFAULT);
  while((scaffold = NextGraphNodeIterator(&scaffolds)) != NULL){
    if(scaffold->type != REAL_SCAFFOLD)
      continue;
    FixupLengthScaffoldT(sgraph, scaffold);
  }

}

/* Check means and variance in a single */
void FixupLengthScaffoldT(ScaffoldGraphT *sgraph, CIScaffoldT *scaffold){
  CIScaffoldTIterator CIs;
  ChunkInstanceT *CI;
  LengthT minOffset, maxOffset;
  LengthT computedLength;
  CDS_COORD_t aEndCoord = -1 , bEndCoord = -1;
  int aEndFound = FALSE;
  int capture;

  minOffset.variance = minOffset.mean = maxOffset.variance = DBL_MAX;
  maxOffset.mean = -DBL_MAX;

  InitCIScaffoldTIterator(sgraph, scaffold, TRUE,   FALSE, &CIs);
  while((CI = NextCIScaffoldTIterator(&CIs)) != NULL){

    assert(CI->scaffoldID == scaffold->id);

    /*	fprintf(stderr,"* Chunk " F_CID " [%g,%g] mean = %g\n",
        chunk->id, 
        chunk->offsetAEnd.mean,
        chunk->offsetBEnd.mean,
        mean); 
    */
    if(CI->flags.bits.cgbType == UU_CGBTYPE){
      capture = TRUE;
      assert(CI->aEndCoord >= 0 && CI->bEndCoord >= 0);
    }else{
      capture = FALSE;
    }

    if(CI->offsetAEnd.mean > CI->offsetBEnd.mean){
      if(!aEndFound && capture){
        aEndCoord = CI->bEndCoord;
        aEndFound = TRUE;
      }
      if(CI->offsetAEnd.mean > maxOffset.mean){
        maxOffset = CI->offsetAEnd;
        if(capture){
          bEndCoord = CI->aEndCoord;
        }
      }
      if(CI->offsetBEnd.mean < minOffset.mean){
        minOffset = CI->offsetBEnd;
      }
    }else{
      if(!aEndFound && capture){
        aEndFound = TRUE;
        aEndCoord = CI->aEndCoord;
      }
      if(CI->offsetBEnd.mean > maxOffset.mean){
        maxOffset = CI->offsetBEnd;
        if(capture){
          bEndCoord = CI->bEndCoord;
        }
      }
      if(CI->offsetAEnd.mean < minOffset.mean){
        minOffset = CI->offsetAEnd;
      }
    }
  }
  ComputeLength(&computedLength, &minOffset, &maxOffset); 

  if(computedLength.mean > scaffold->bpLength.mean){
    fprintf(stderr,"* Adjusting scaffold " F_CID " bplength from %g to %g [%g,%g]\n",
            scaffold->id, scaffold->bpLength.mean, computedLength.mean,
            minOffset.mean, maxOffset.mean);

    scaffold->bpLength = computedLength;

  }

  if(aEndCoord < 0 ||
     bEndCoord < 0 ){
    fprintf(stderr,">>>>* Scaffold " F_CID " has aEndCoord :" F_COORD " bEndCoord :" F_COORD "  minOffset = %g maxOffset = %g\n",
            scaffold->id, aEndCoord, bEndCoord, minOffset.mean, maxOffset.mean);
  }
  scaffold->bEndCoord = bEndCoord;
  scaffold->aEndCoord = aEndCoord;
}



//  DemoteSmallSingletonScaffolds
//
//  We want to demote the contigs/unitigs in small singleton scaffolds
//  so that they can be candidates for stone/rock throwing.
//
void DemoteSmallSingletonScaffolds(void) {
  GraphNodeIterator   scaffolds;
  CIScaffoldT        *scaffold;
  int                 numScaffolds = 0;
  int                 numSingletonScaffolds = 0;
  int                 numDemoted= 0;

  InitGraphNodeIterator(&scaffolds, ScaffoldGraph->ScaffoldGraph, GRAPH_NODE_DEFAULT);
  while ((scaffold = NextGraphNodeIterator(&scaffolds)) != NULL) {
    ContigT        *contig;
    ChunkInstanceT *CI;

    if (scaffold->info.Scaffold.numElements > 1)
      continue;

    contig = GetGraphNode(ScaffoldGraph->ContigGraph, scaffold->info.Scaffold.AEndCI);

    numScaffolds++;
    if (contig->info.Contig.numCI > 1)
      continue;

    CI = GetGraphNode(ScaffoldGraph->CIGraph, contig->info.Contig.AEndCI);

    numSingletonScaffolds++;
    if ((CI->info.CI.coverageStat > GlobalData->cgbDefinitelyUniqueCutoff) ||
        (CI->bpLength.mean > 2000.0))
      continue;

    // We've found a victim!!!

    numDemoted++;

    fprintf(GlobalData->stderrc,
            "** Demoting Contig/Unitig " F_CID "/" F_CID " with coverage stat %d length %g scaffold " F_CID "\n",
            contig->id, CI->id, CI->info.CI.coverageStat, scaffold->bpLength.mean, scaffold->id);

    // Remove the Contig from the Scaffold.  We don't need to use the
    // RemoveCIFromScaffold machinery, since we are dealing with a
    // pathological case
    //
    contig->flags.bits.isUnique = 0;
    contig->scaffoldID          = NULLINDEX;
    contig->AEndNext            = NULLINDEX;
    contig->BEndNext            = NULLINDEX;
    
    // Delete any remaining edges
    DeleteScaffoldEdgesForScaffold(ScaffoldGraph, scaffold);

    // Mark the scaffold dead
    scaffold->flags.bits.isDead         = TRUE;
    scaffold->info.Scaffold.AEndCI      = NULLINDEX;
    scaffold->info.Scaffold.BEndCI      = NULLINDEX;
    scaffold->info.Scaffold.numElements = 0;
    scaffold->bpLength.mean             = 0.0;
    scaffold->bpLength.variance         = 0.0;

    // Mark the Underlying Unitig as un-scaffolded, and not-unique
    SetNodeType(CI, UNRESOLVEDCHUNK_CGW);
  }

  //  If we removed any scaffolds, rebuild all the edges.
  //
  if (numDemoted > 0)
    BuildNewScaffoldEdges(ScaffoldGraph, 0);

  fprintf(GlobalData->stderrc,
          "# Considered %d scaffolds of which %d were single and %d (%g%%) were demoted\n",
          numScaffolds, numSingletonScaffolds, numDemoted, 
          (numSingletonScaffolds > 0? ((double)(numDemoted)/(double)(numSingletonScaffolds)): 0.0));
}
