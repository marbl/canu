
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
static char *rcsid = "$Id: CIScaffoldT_CGW.c,v 1.72 2012-08-22 01:24:26 jasonmiller9704 Exp $";

#undef DEBUG_INSERT
#undef DEBUG_DIAG
#undef DEBUG_SPLIT

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_UTL_interval.h"
#include "AS_CGW_dataTypes.h"
#include "Globals_CGW.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"
#include "UnionFind_AS.h"
#include "UtilsREZ.h"
#include "ChiSquareTest_CGW.h"
#include "MultiAlignment_CNS.h"
#include "DataTypesREZ.h"
#include "CommonREZ.h"
#include "Stats_CGW.h"   // for collecting scaffold merging stats

VA_DEF(PtrT)


void PrintCINodeFields(FILE * stream, NodeCGW_T * node)
{
  fprintf(stream, "\t\tcontigID:"F_CID"\n", node->info.CI.contigID);
  fprintf(stream, "\t\tnumFragments:%d\n", ScaffoldGraph->tigStore->getNumFrags(node->id, TRUE));
  fprintf(stream, "\t\tcoverageStat:%d\n", ScaffoldGraph->tigStore->getUnitigCoverageStat(node->id));
  fprintf(stream, "\t\tbaseID:"F_CID"\n", node->info.CI.baseID);
  fprintf(stream, "\t\tnumInstances:%d\n", node->info.CI.numInstances);
}

void PrintContigNodeFields(FILE * stream, NodeCGW_T * node)
{
  fprintf(stream, "\t\tAEndCI:"F_CID", BEndCI:"F_CID", numCI:%d\n",
          node->info.Contig.AEndCI, node->info.Contig.BEndCI,
          node->info.Contig.numCI);
}

void PrintScaffoldNodeFields(FILE * stream, NodeCGW_T * node)
{
  fprintf(stream, "\t\tAEndCI:"F_CID", BEndCI:"F_CID", numElements:%d\n",
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
  fprintf(stream, "\t\tisDead:%d, isFree:%d, containsCIs:%d\n",
          node->flags.bits.isDead,
          node->flags.bits.isFree, node->flags.bits.containsCIs);
  fprintf(stream,
          "\t\tisCI:%d, isContig:%d, isScaffold:%d\n",
          node->flags.bits.isCI,
          node->flags.bits.isContig, node->flags.bits.isScaffold);
  fprintf(stream, "\t\tisSurrogate:%d, beingContigged:%d\n",
          node->flags.bits.isSurrogate, node->flags.bits.beingContigged);
  fprintf(stream, "\t\twalkedAlready:%d, walkedTooShort:%d, walkedTooLong:%d\n",
          node->flags.bits.walkedAlready, node->flags.bits.walkedTooShort,
          node->flags.bits.walkedTooLong);
  fprintf(stream, "\t\twalkMaxedOut:%d, walkedTrivial:%d\n",
          node->flags.bits.walkMaxedOut, node->flags.bits.walkedTrivial);
  fprintf(stream, "\t\tisStoneSurrogate:%d, isWalkSurrogate:%d\n",
          node->flags.bits.isStoneSurrogate, node->flags.bits.isWalkSurrogate);
  fprintf(stream, "\t\tfailedToContig:%d, isChaff:%d, isStone:%d\n",
          node->flags.bits.failedToContig, node->flags.bits.isChaff,
          node->flags.bits.isStone);
  fprintf(stream, "\t\tisWalk:%d, isRock:%d, isPotentialRock:%d, isPotentialStone:%d\n",
          node->flags.bits.isWalk, node->flags.bits.isRock,
          node->flags.bits.isPotentialRock,
          node->flags.bits.isPotentialStone);
  fprintf(stream, "\tall:%d\n", node->flags.all);
}

void PrintNodeFields(FILE * stream, NodeCGW_T * node)
{
  fprintf(stream,"\ttype:%c, scaffoldID:"F_CID", prevScaffoldID:"F_CID"\n",
          node->type, node->scaffoldID, node->prevScaffoldID);
  fprintf(stream,"\tindexInScaffold:%d, smoothExpectedCID:"F_CID"\n",
          node->indexInScaffold, node->smoothExpectedCID);
  fprintf(stream, "\tnumEssentialA:%d, numEssentialB:%d\n",
          node->numEssentialA, node->numEssentialB);
  fprintf(stream, "\tessentialEdgeA:"F_CID", essentialEdgeB:"F_CID"\n",
          node->essentialEdgeA, node->essentialEdgeB);
  fprintf(stream,"\tAEndNext:"F_CID", BEndNext:"F_CID"\n",
          node->AEndNext, node->BEndNext);
  fprintf(stream, "\tbpLength:(%.1f,%.1f), offsetAEnd:(%.1f,%.1f), offsetBEnd:(%.1f,%.1f)\n",
          node->bpLength.mean, node->bpLength.variance,
          node->offsetAEnd.mean, node->offsetAEnd.variance,
          node->offsetBEnd.mean, node->offsetBEnd.variance);
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
  fprintf(stream, "\tedgeHead:"F_CID", setID:"F_CID"\n",
          node->edgeHead, node->setID);
}


void PrintCIScaffoldHeader(FILE *stream, ScaffoldGraphT *graph, CIScaffoldT *scaffold){
  fprintf(stream,"\n* CIScaffold "F_CID" numCI:%d (a:"F_CID" b:"F_CID")  length: %d\n",
	  scaffold->id, scaffold->info.Scaffold.numElements, scaffold->info.Scaffold.AEndCI, scaffold->info.Scaffold.BEndCI,
	  (int)scaffold->bpLength.mean);
  // PrintNodeFields(stream, scaffold);
}



void DumpCIScaffold(FILE *stream, ScaffoldGraphT *graph, CIScaffoldT *scaffold, int raw){

  PrintCIScaffoldHeader(stream, graph, scaffold);
  fprintf(stream, "> Includes CIs\n");

  CIScaffoldTIterator CIs;
  ChunkInstanceT     *CI;

  InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &CIs);

  while((CI = NextCIScaffoldTIterator(&CIs)) != NULL){
    fprintf(stream," \t %5d: CI %8d sid %8d len %7d,%9.1f ends %7d,%7d var %9.1f,%9.1f orient %c\n",
            CI->indexInScaffold,
            CI->id,
            CI->scaffoldID,
            (int)CI->bpLength.mean, CI->bpLength.variance,
            (int)CI->offsetAEnd.mean, (int)CI->offsetBEnd.mean,
            CI->offsetAEnd.variance, CI->offsetBEnd.variance,
            GetNodeOrient(CI).toLetter());
  }
  fprintf(stream, "> %s Edges A \n",  (raw?" R ":" M "));

  {
    GraphEdgeIterator  SEdges(graph->ScaffoldGraph, scaffold->id, A_END, ALL_EDGES);
    SEdgeT            *SEdge;

    while((SEdge = (raw) ? SEdges.nextRaw() : SEdges.nextMerged()) != NULL){
      PrintSEdgeT(stream, graph, " ", SEdge, scaffold->id);
    }
  }

  fprintf(stream, "> %s Edges B\n", (raw?" R ":" M "));

  {
    GraphEdgeIterator  SEdges(graph->ScaffoldGraph, scaffold->id, B_END, ALL_EDGES);
    SEdgeT            *SEdge;

    while((SEdge = (raw) ? SEdges.nextRaw() : SEdges.nextMerged()) != NULL){
      PrintSEdgeT(stream, graph, " ", SEdge, scaffold->id);
    }
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
  for (uint32 sid = 0; sid < GetNumCIScaffoldTs(graph->CIScaffolds); sid++){
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




//  Insert chunk instance ci int scaffold sid at offset with
//  orientation orient.
//
//  offsetFromAEnd = offset of the end of the CI that is closest to the A end
//  of the scaffold.
//
//  If the CI has edges that are marked isContigConfirming, it will be
//  merged into a contig with the indicated CIs.  The edges must
//  either be non-tandem overlap singleton overlap edges, or must have
//  a distance variance of less than N base pairs.
//
void
InsertCIInScaffold(ScaffoldGraphT *sgraph,
                   CDS_CID_t ci,
                   CDS_CID_t sid,
                   LengthT aEndOffset,
                   LengthT bEndOffset,
                   int AEndToBend,
                   int contigNow) {

  CIScaffoldT    *ciScaffold = GetGraphNode(sgraph->ScaffoldGraph, sid);
  ChunkInstanceT *chunkInstance = GetGraphNode(sgraph->ContigGraph, ci);

  assert(chunkInstance->flags.bits.isDead == FALSE);
  assert(ciScaffold->flags.bits.isDead    == FALSE);

  //  We used to check that:
  //  --  contig end point variance was within what the actual contig said it was (at most +1 from truth)
  //  --  contig end point distance was within what the actual contig said it was (at most 1.2x from truth)
  //
  //  For variance, we just reset the right end.
  //  For distance, we just ignored it.
  //
  //  Both spitting out much logging.
  //
  //  We now trust the positions.  Any checking would have been relaxed to where they would have been useless.
  //
  //  Motivated by this -- rounding error? -- example:
  //    supplied placement has larger variance (100.00000) than 1.0 + chunk (99.80000).


  //  Check that the CI isn't already in the scaffold.
  {
    CIScaffoldTIterator CIs;
    ChunkInstanceT     *CI;

    InitCIScaffoldTIterator(sgraph, ciScaffold, TRUE, FALSE, &CIs);
    while ((CI = NextCIScaffoldTIterator(&CIs)) != NULL) {
      if (CI->id == ci)
	return;
    }
  }

  chunkInstance->offsetAEnd = aEndOffset;
  chunkInstance->offsetBEnd = bEndOffset;

#ifdef DEBUG_INSERT
  fprintf(stderr,"InsertCIInScaffold()--  Before Inserting:\n");
  DumpCIScaffold(stderr, sgraph, ciScaffold, FALSE);
#endif


  if ((contigNow == TRUE) &&
      (CheckForContigs(sgraph, ci, sid, aEndOffset, bEndOffset))) {
#ifdef DEBUG_INSERT
    fprintf(stderr,"InsertCIInScaffold()--  After Inserting:\n");
    DumpCIScaffold(stderr,sgraph, ciScaffold, FALSE);
#endif
    return;
  }

  //
  //  New contig doesn't intersect any other contig, insert as is.
  //

  //  Reget the pointers.
  ciScaffold    = GetGraphNode(sgraph->ScaffoldGraph, sid);
  chunkInstance = GetGraphNode(sgraph->ContigGraph, ci);

  assert(!chunkInstance->flags.bits.isDead);
  MarkCIElementsForScaffoldMembership(chunkInstance, sid);

  LengthT *minOffset = (aEndOffset.mean > bEndOffset.mean) ? &bEndOffset : &aEndOffset;
  LengthT *maxOffset = (aEndOffset.mean > bEndOffset.mean) ? &aEndOffset : &bEndOffset;

  chunkInstance->scaffoldID         = sid;
  chunkInstance->flags.bits.isChaff = FALSE;  // we need this one in the output

  ciScaffold->info.Scaffold.numElements++;

  if (ciScaffold->info.Scaffold.AEndCI == NULLINDEX) {
    //  Inserting into an empty scaffold
    chunkInstance->AEndNext          = NULLINDEX;
    chunkInstance->BEndNext          = NULLINDEX;
    ciScaffold->info.Scaffold.AEndCI = ci;
    ciScaffold->info.Scaffold.BEndCI = ci;
    ciScaffold->bpLength             = *maxOffset;
    return;
  }

  CIScaffoldTIterator CIs;
  ChunkInstanceT     *CI;

  int32 chunkInstanceMin = (int32)MIN(chunkInstance->offsetAEnd.mean, chunkInstance->offsetBEnd.mean);

  InitCIScaffoldTIterator(sgraph, ciScaffold, AEndToBend, FALSE, &CIs);
  while((CI = NextCIScaffoldTIterator(&CIs)) != NULL){
    int32 CImin = (int32) MIN(CI->offsetAEnd.mean, CI->offsetBEnd.mean);

    if (CImin > chunkInstanceMin) {
      // WARNING: this condition is ok ONLY if AEndToBend == TRUE
      // When we traverse the list of CIs from the Bend, the condition
      // (CImin > chunkInstanceMin) will be satisfied immediately and therefore
      // the chunk will end up in the wrong position (SteLo)
      //
      assert(AEndToBend);

      chunkInstance->BEndNext = CI->id;
      chunkInstance->AEndNext = CI->AEndNext;

      //  Set the previous CI's next pointer to us
      if (CI->AEndNext != NULLINDEX)
        GetGraphNode(sgraph->ContigGraph, CI->AEndNext)->BEndNext = ci;

      CI->AEndNext = ci;

      if (CI->id == ciScaffold->info.Scaffold.AEndCI)
        ciScaffold->info.Scaffold.AEndCI = ci;

      break;

    } else if ((CI->id == ciScaffold->info.Scaffold.BEndCI) &&
               (CImin <= chunkInstanceMin)) {
      // append
      ciScaffold->info.Scaffold.BEndCI = ci;
      CI->BEndNext                     = ci;
      chunkInstance->AEndNext          = CI->id;
      chunkInstance->BEndNext          = NULLINDEX;

      // Due to containments, the CI with the maximal mean does not
      // have the maximal variance.

      if (ciScaffold->bpLength.mean < maxOffset->mean)
        ciScaffold->bpLength.mean = maxOffset->mean;

      if (ciScaffold->bpLength.variance < maxOffset->variance)
        ciScaffold->bpLength.variance = maxOffset->variance;

      break;
    }
  }  //  end of while()


  // Due to containments, the CI with the maximal mean does not
  // have the maximal variance.

  if(ciScaffold->bpLength.mean < maxOffset->mean)
    ciScaffold->bpLength.mean = maxOffset->mean;

  if(ciScaffold->bpLength.variance < maxOffset->variance)
    ciScaffold->bpLength.variance = maxOffset->variance;
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
  LengthT base = {0,0};
  LengthT maxoffset = (CI->offsetAEnd.mean < CI->offsetBEnd.mean? CI->offsetBEnd:CI->offsetAEnd);

#if 0
  if(ciScaffold->info.Scaffold.numElements < 3){;
  fprintf(stderr,"* Removing CI "F_CID" from scaffold "F_CID", elements left %d\nBEFORE:\n",
	  CI->id, ciScaffold->id, ciScaffold->info.Scaffold.numElements);

  DumpCIScaffold(stderr,sgraph, ciScaffold, FALSE);
  }
#endif

  assert(ciScaffold->info.Scaffold.AEndCI != NULLINDEX && ciScaffold->info.Scaffold.BEndCI != NULLINDEX);
  assert(ciScaffold && (CI->scaffoldID == ciScaffold->id));

  ciScaffold = GetGraphNode(sgraph->ScaffoldGraph, CI->scaffoldID);

  assert(ciScaffold && !isDeadCIScaffoldT(ciScaffold));

  if(CI->AEndNext != NULLINDEX)
    anext = GetGraphNode(sgraph->ContigGraph, CI->AEndNext);

  if(CI->BEndNext != NULLINDEX)
    bnext = GetGraphNode(sgraph->ContigGraph, CI->BEndNext);

#if 0
  fprintf(stderr,"* Predecessor is "F_CID" Successor is "F_CID"\n",
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
      fprintf(stderr,"* bneighbor "F_CID" has AEndNext "F_CID" and BEndNext "F_CID"\n",
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

      SetCIScaffoldTLength(ScaffoldGraph, ciScaffold);

      //  ciScaffold->bpLength.mean = MAX( prevCI->offsetAEnd.mean, prevCI->offsetBEnd.mean);
      //  ciScaffold->bpLength.variance = MAX( prevCI->offsetAEnd.variance, prevCI->offsetBEnd.variance);


#if 0
      fprintf(stderr,"* aneighbor "F_CID" has AEndNext "F_CID" and BEndNext "F_CID"\n",
	      anext->id, anext->AEndNext, anext->BEndNext);
#endif
    }
    middle = FALSE;
  }

  if(middle){

    bnext->AEndNext = CI->AEndNext;
    anext->BEndNext = CI->BEndNext;

    if(fabs(maxoffset.mean - ciScaffold->bpLength.mean) < 100){
      SetCIScaffoldTLength(ScaffoldGraph, ciScaffold);
    }

#if 0
    fprintf(stderr,"* bneighbor "F_CID" has AEndNext "F_CID" and BEndNext "F_CID"\n",
	    bnext->id, bnext->AEndNext, bnext->BEndNext);
    fprintf(stderr,"* aneighbor "F_CID" has AEndNext "F_CID" and BEndNext "F_CID"\n",
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
    AddDeltaToScaffoldOffsets(ScaffoldGraph, ciScaffold->id, bnext->id, TRUE, base);

    fprintf(stderr,"* After RemoveCIFromScaffold ci:"F_CID" scaffold:"F_CID" removing from a-end base = (%g,%g), bpLength = (%g,%g)\n",
            CI->id,
            ciScaffold->id,
            base.mean, base.variance,
            ciScaffold->bpLength.mean, ciScaffold->bpLength.variance);
  }
  ciScaffold->info.Scaffold.numElements--;
#if 0
  fprintf(stderr,"* Removing CI "F_CID" froms scaffold "F_CID"\nAFTER:\n",
	  CI->id, ciScaffold->id);
  DumpCIScaffold(stderr,sgraph, ciScaffold, FALSE);
#endif
  if(ciScaffold->info.Scaffold.numElements == 0){
    ciScaffold->bpLength.mean = 0.0;
    ciScaffold->bpLength.variance = 0.0;
  }
  return FALSE;
}




//  Determines whether the scaffold is connected by edges marked TRUSTED and
//  TENTATIVELY_TRUSTED. This is a necessary condition for boths sanity and successful recomputation
//  of positions of Scaffold CI positions.  Also interesting to evaluate this after
//  MarkInternalCIEdgeStatus. edgeTypes defines the set of edges used.  LeastSquares uses
//  ALL_TRUSTED_EDGES, other manipulations use ALL_EDGES.
//
//  returns the number of connected components of the <scaffold>
//
//  Will modify the setId field of the NodeCGW_T structure to reflect which component a node belongs
//  to.
//
int
IsScaffoldInternallyConnected(ScaffoldGraphT *sgraph,
                              CIScaffoldT    *scaffold,
                              bool            useMerged,
                              bool            useTrusted) {

  UFDataT * UFData = UFCreateSets(scaffold->info.Scaffold.numElements);
  ChunkInstanceT * chunk;
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
    UFSetT   *chunkSet = UFGetSet(UFData, set);

    chunkSet->data = chunk;
    chunk->setID   = set++;
  }

  //
  // now do the unions: iterate over all trusted/raw edges
  //
  InitCIScaffoldTIterator(sgraph, scaffold, TRUE,
                          FALSE, &CIs);
  while ((chunk = NextCIScaffoldTIterator(&CIs)) != NULL) {
    GraphEdgeIterator   edges(sgraph->ContigGraph, chunk->id, ALL_END, (useTrusted) ? ALL_TRUSTED_EDGES : ALL_EDGES);
    CIEdgeT            *edge;

    assert(chunk->setID >= 0);

    while ((edge = (useMerged) ? edges.nextMerged() : edges.nextRaw()) != NULL) {
      ChunkInstanceT *otherChunk = GetGraphNode(sgraph->ContigGraph, (chunk->id == edge->idA) ? edge->idB : edge->idA);
      int32           weight = edge->edgesContributing - (isOverlapEdge(edge));

      // See each edge only once
      if(chunk->id != edge->idA)
        continue;

      if(isSingletonOverlapEdge(edge) ||
         (weight == 1 && edge->flags.bits.isBridge))
        continue;

      if (chunk->scaffoldID != otherChunk->scaffoldID)
        continue;

      UFUnion(UFData, chunk->setID, otherChunk->setID);
    }
        
    // merge unions based on closure reads as well (i.e. consider them edges)
    if (chunk->flags.bits.isClosure) {
       MultiAlignT *ma = ScaffoldGraph->tigStore->loadMultiAlign(chunk->id, chunk->flags.bits.isCI);
       int i = 0;
       assert(ma != NULL);
       
       for(i = 0; i < GetNumIntMultiPoss(ma->f_list); i++) {      
          IntMultiPos *mp = GetIntMultiPos(ma->f_list, i);
          
          gkPlacement *gkpl = ScaffoldGraph->gkpStore->gkStore_getReadPlacement(mp->ident);
          if (gkpl == NULL) {
            continue;
          }
          assert(gkpl->bound1);
          assert(gkpl->bound2);
   
          // get the reads indicated by the input line
          CIFragT *leftMate = GetCIFragT(ScaffoldGraph->CIFrags, gkpl->bound1); 
          CIFragT *rightMate = GetCIFragT(ScaffoldGraph->CIFrags, gkpl->bound2);
          if (leftMate->contigID == NULLINDEX || rightMate->contigID == NULLINDEX) {
            continue;
          }
          ChunkInstanceT * begin_chunk = GetGraphNode(ScaffoldGraph->ContigGraph, leftMate->contigID);
          ChunkInstanceT * end_chunk   = GetGraphNode(ScaffoldGraph->ContigGraph, rightMate->contigID);
          
          if (chunk->scaffoldID != begin_chunk->scaffoldID) {
            continue;
          }
          
          if (begin_chunk->scaffoldID != end_chunk->scaffoldID) {
            continue;
          }       
          UFUnion(UFData, chunk->setID, begin_chunk->setID);
          UFUnion(UFData, chunk->setID, end_chunk->setID);
       }
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

static
void
DeleteScaffoldEdgesForScaffold(ScaffoldGraphT * graph, CIScaffoldT * scaffold) {
  while(scaffold->edgeHead != NULLINDEX)
      DeleteGraphEdge(graph->ScaffoldGraph,
                      GetGraphEdge(graph->ScaffoldGraph,
                                   scaffold->edgeHead));
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
        if (ScaffoldGraph->tigStore->getNumFrags(basechunk->id, TRUE) == 1)
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

        assert(basechunk->info.CI.numInstances == GetNumCDS_CID_ts(basechunk->info.CI.instances.va));

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

        assert(basechunk->info.CI.numInstances == GetNumCDS_CID_ts(basechunk->info.CI.instances.va));

        for (index=0; index<basechunk->info.CI.numInstances; index++)
          if (*GetCDS_CID_t(basechunk->info.CI.instances.va, index) == chunk->id)
            SetCDS_CID_t(basechunk->info.CI.instances.va,
                         index,
                         GetCDS_CID_t(basechunk->info.CI.instances.va, basechunk->info.CI.numInstances-1));

        basechunk->info.CI.numInstances--;

        ResetToRangeVA_CDS_CID_t(basechunk->info.CI.instances.va, basechunk->info.CI.numInstances);
        assert(basechunk->info.CI.numInstances == GetNumCDS_CID_ts(basechunk->info.CI.instances.va));
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



int32
CheckScaffoldConnectivityAndSplit(ScaffoldGraphT *graph,
                                  CDS_CID_t       scaffoldID,
                                  int32           edgeTypes,
                                  int             verbose) {
  CIScaffoldT  *scaffold      = GetCIScaffoldT(graph->CIScaffolds, scaffoldID);

  //  Test connectedness, mark bridge edges
  if (IsScaffold2EdgeConnected(ScaffoldGraph, scaffold) == true)
    return(1);

  // Connected component analysis, mark contigs with their component number
  int32  numComponents = IsScaffoldInternallyConnected(graph, scaffold, true, (edgeTypes == ALL_TRUSTED_EDGES));

  assert(numComponents > 0);

  if (numComponents == 1)
    //  Expected case, all connected.
    return(numComponents);

  fprintf(stderr, "Scaffold "F_CID" is not connected, has %d components.\n", scaffoldID, numComponents);
  fprintf(stderr, "Splitting into scaffolds: (search for \"Splitting "F_CID" into scaffold\" to get the new scaffolds)\n", scaffoldID);

  //DumpACIScaffold(stderr,graph, scaffold, FALSE);

  //fprintf(stderr,"Prior to split ...");
  //DumpACIScaffoldNew(stderr,graph, scaffold, FALSE);

  int32               contigMax = 0;
  CDS_CID_t          *contigID  = (CDS_CID_t *)safe_malloc(sizeof(CDS_CID_t) * scaffold->info.Scaffold.numElements);

  CIScaffoldTIterator contigs;
  NodeCGW_T          *contig;

  InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &contigs);

  while ((contig = NextCIScaffoldTIterator(&contigs)) != NULL)
    contigID[contigMax++] = contig->id;

  assert(contigMax == scaffold->info.Scaffold.numElements);

  // For each component, create a scaffold and insert the relevant contigs

  for (int32 component=0; component<numComponents; component++) {
    CIScaffoldT newScaffold;

    InitializeScaffold(&newScaffold, REAL_SCAFFOLD);

    newScaffold.info.Scaffold.AEndCI      = NULLINDEX;
    newScaffold.info.Scaffold.BEndCI      = NULLINDEX;
    newScaffold.info.Scaffold.numElements = 0;
    newScaffold.edgeHead                  = NULLINDEX;
    newScaffold.bpLength.mean             = 0.0;
    newScaffold.bpLength.variance         = 0.0;
    newScaffold.id                        = GetNumGraphNodes(graph->ScaffoldGraph);
    newScaffold.flags.bits.isDead         = FALSE;
    newScaffold.numEssentialA             = 0;
    newScaffold.numEssentialB             = 0;
    newScaffold.essentialEdgeB            = NULLINDEX;
    newScaffold.essentialEdgeA            = NULLINDEX;

    //CDS_CID_t newScaffoldID = newScaffold.id;

    fprintf(stderr, "Splitting "F_CID" into scaffold "F_CID"\n", scaffoldID, newScaffold.id);

    AppendGraphNode(graph->ScaffoldGraph, &newScaffold);

    scaffold = GetCIScaffoldT(graph->CIScaffolds, scaffoldID);

    LengthT  firstOffset = { DBL_MAX, DBL_MAX };

    for (int32 ic = 0; ic < contigMax; ic++) {
      contig = GetGraphNode(graph->ContigGraph, contigID[ic]);

      if (contig->setID != component)
        continue;

      RemoveCIFromScaffold(graph, scaffold, contig, FALSE);

      if (firstOffset.mean == DBL_MAX)
        firstOffset = (GetNodeOrient(contig).isForward()) ? contig->offsetAEnd : contig->offsetBEnd;

      contig->offsetAEnd.mean     -= firstOffset.mean;
      contig->offsetAEnd.variance -= firstOffset.variance;
      contig->offsetBEnd.mean     -= firstOffset.mean;
      contig->offsetBEnd.variance -= firstOffset.variance;

      InsertCIInScaffold(graph, contig->id, newScaffold.id, contig->offsetAEnd, contig->offsetBEnd, TRUE, FALSE);
    }

    assert((GetGraphNode(graph->ScaffoldGraph, newScaffold.id))->info.Scaffold.numElements > 0);

    //  Make sure that our new scaffold contains more than just a single surrogate.
    //
    killScaffoldIfOnlySurrogate(newScaffold.id);
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

  safe_free(contigID);

  return(numComponents);
}



/*****************************************************************************/


void
SetCIScaffoldTLength(ScaffoldGraphT *sgraph, CIScaffoldT *scaffold) {
  CIScaffoldTIterator CIs;
  ChunkInstanceT     *CI;
  LengthT             maxOffset = {0.0,0.0};

  if (scaffold->type != REAL_SCAFFOLD)
    return;

  if (scaffold->flags.bits.isDead)
    //  Should assert, but somebody calls us on dead scaffolds.
    return;

  assert(scaffold->flags.bits.isDead == FALSE);

  InitCIScaffoldTIterator(sgraph, scaffold, TRUE, FALSE, &CIs);
  while((CI = NextCIScaffoldTIterator(&CIs)) != NULL){
    if (CI->offsetAEnd.mean     > maxOffset.mean)     maxOffset.mean     = CI->offsetAEnd.mean;
    if (CI->offsetBEnd.mean     > maxOffset.mean)     maxOffset.mean     = CI->offsetBEnd.mean;
    if (CI->offsetAEnd.variance > maxOffset.variance) maxOffset.variance = CI->offsetAEnd.variance;
    if (CI->offsetBEnd.variance > maxOffset.variance) maxOffset.variance = CI->offsetBEnd.variance;
  }

  if ((scaffold->bpLength.mean     != maxOffset.mean) ||
      (scaffold->bpLength.variance != maxOffset.variance))
    fprintf(stderr, "SetCIScaffoldTLength()-- adjusted scaffold "F_CID" from (%.0f +- %.0f) to (%.0f +- %.0f)\n",
            scaffold->id,
            scaffold->bpLength.mean, scaffold->bpLength.variance,
            maxOffset.mean,          maxOffset.variance);

  scaffold->bpLength = maxOffset;
}


void
SetCIScaffoldTLengths(ScaffoldGraphT *sgraph) {
  GraphNodeIterator scaffolds;
  CIScaffoldT      *scaffold;

  InitGraphNodeIterator(&scaffolds, sgraph->ScaffoldGraph, GRAPH_NODE_DEFAULT);
  while ((scaffold = NextGraphNodeIterator(&scaffolds)) != NULL)
    SetCIScaffoldTLength(sgraph, scaffold);
}




bool
AdjustNegativePositions(ScaffoldGraphT *graph, CIScaffoldT *scaffold);
bool
AdjustNegativeVariances(ScaffoldGraphT *graph, CIScaffoldT *scaffold);



void
ScaffoldSanity(ScaffoldGraphT *graph, CIScaffoldT *scaffold) {

  //sid = scaffold->id;
  //LSE = scaffold->info.Scaffold.leastSquareError;

  //fprintf(stderr, "ScaffoldSanity()-- Checking scaffold %d\n", scaffold->id);

  assert(scaffold->flags.bits.isScaffold == true);

  if (scaffold->flags.bits.isDead == true) {
    fprintf(stderr, "ScaffoldSanity()-- WARNING: Called on dead scaffold %d\n", scaffold->id);
    return;
  }
  assert(scaffold->flags.bits.isDead     == false);

  if (scaffold->type != REAL_SCAFFOLD)
    return;

  if (scaffold->info.Scaffold.numElements == 1)
    return;

  CIScaffoldTIterator  CIs;
  ChunkInstanceT      *CI;

  int32   numCtg  = 0;

  LengthT lastMin = {0, 0};
  LengthT lastMax = {0, 0};

  uint32  hasProblems = 0;

  InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &CIs);

  while (NULL != (CI = NextCIScaffoldTIterator(&CIs))) {
    LengthT thisMin = (CI->offsetAEnd.mean < CI->offsetBEnd.mean) ? CI->offsetAEnd : CI->offsetBEnd;
    LengthT thisMax = (CI->offsetAEnd.mean > CI->offsetBEnd.mean) ? CI->offsetAEnd : CI->offsetBEnd;

    //fprintf(stderr, "CI %d %d %.0f +- %.0f %.0f +- %.0f -- surrogate %d stone %d rock %d\n",
    //        CI->id, CI->type,
    //        thisMin.mean, thisMin.variance,
    //        thisMax.mean, thisMax.variance,
    //        CI->flags.bits.isSurrogate,
    //        CI->flags.bits.isStone,
    //        CI->flags.bits.isRock);

    numCtg++;

    if (CI->scaffoldID != scaffold->id)
      fprintf(stderr, "ScaffoldSanity()--  contig %d in scaffold %d -- claims scaffold %d\n", CI->id, scaffold->id, CI->scaffoldID), hasProblems++;

    if (CI->flags.bits.isUnique != true)
      fprintf(stderr, "ScaffoldSanity()--  contig %d in scaffold %d -- marked as not unique\n", CI->id, scaffold->id), hasProblems++;
    if (CI->flags.bits.isDead   != false)
      fprintf(stderr, "ScaffoldSanity()--  contig %d in scaffold %d -- marked as dead\n", CI->id, scaffold->id), hasProblems++;

    if (0.0 > CI->offsetAEnd.mean)
      fprintf(stderr, "ScaffoldSanity()--  contig %d in scaffold %d -- offsetAEnd.mean %f negative\n", CI->id, scaffold->id, CI->offsetAEnd.mean), hasProblems++;
    if (0.0 > CI->offsetBEnd.mean)
      fprintf(stderr, "ScaffoldSanity()--  contig %d in scaffold %d -- offsetAEnd.mean %f negative\n", CI->id, scaffold->id, CI->offsetBEnd.mean), hasProblems++;

    if (CI->offsetAEnd.mean > scaffold->bpLength.mean)
      fprintf(stderr, "ScaffoldSanity()--  contig %d in scaffold %d -- offsetAEnd.mean %f beyond scaffold end %f\n", CI->id, scaffold->id, CI->offsetAEnd.mean, scaffold->bpLength.mean), hasProblems++;
    if (CI->offsetBEnd.mean > scaffold->bpLength.mean)
      fprintf(stderr, "ScaffoldSanity()--  contig %d in scaffold %d -- offsetBEnd.mean %f beyond scaffold end %f\n", CI->id, scaffold->id, CI->offsetBEnd.mean, scaffold->bpLength.mean), hasProblems++;

    if (0.0 > CI->offsetAEnd.variance)
      fprintf(stderr, "ScaffoldSanity()--  contig %d in scaffold %d -- offsetAEnd.variance %f negative\n", CI->id, scaffold->id, CI->offsetBEnd.variance), hasProblems++;
    if (0.0 > CI->offsetBEnd.variance)
      fprintf(stderr, "ScaffoldSanity()--  contig %d in scaffold %d -- offsetBEnd.variance %f negative\n", CI->id, scaffold->id, CI->offsetBEnd.variance), hasProblems++;

    //  +1 for rounding errors
    if (CI->offsetAEnd.variance > scaffold->bpLength.variance + 1.0)
      fprintf(stderr, "ScaffoldSanity()--  contig %d in scaffold %d -- offsetAEnd.variance %f beyond scaffold end %f\n", CI->id, scaffold->id, CI->offsetAEnd.variance, scaffold->bpLength.variance), hasProblems++;
    if (CI->offsetBEnd.variance > scaffold->bpLength.variance + 1.0)
      fprintf(stderr, "ScaffoldSanity()--  contig %d in scaffold %d -- offsetAEnd.variance %f beyond scaffold end %f\n", CI->id, scaffold->id, CI->offsetBEnd.variance, scaffold->bpLength.variance), hasProblems++;

    if (numCtg == 0)
      if ((CI->offsetAEnd.mean != 0.0) && (CI->offsetBEnd.mean != 0.0))
        fprintf(stderr, "ScaffoldSanity()--  contig %d in scaffold %d -- first contig at %f,%f not zero\n", CI->id, scaffold->id, CI->offsetAEnd.mean, CI->offsetBEnd.mean), hasProblems++;

    //  Variances of rocks/stones/surrogates are not correct (these might be caught by the containment check below).

    if (CI->flags.bits.isRock == true)
      continue;
    if (CI->flags.bits.isStone == true)
      continue;
    if (CI->flags.bits.isSurrogate == true)
      continue;

    //  The +1 is guarding against rounding problems, in particular, this:
    //    (gdb) p lastMin  $5 = {mean = 2515.7648383446817, variance = 41392.301712181623}
    //    (gdb) p thisMin  $6 = {mean = 2515.7648383446813, variance = 41712.021995136231}
    //                                     0.0000000000004
    //
    if (lastMin.mean > thisMin.mean + 1.0)
      fprintf(stderr, "ScaffoldSanity()--  contig %d in scaffold %d -- begins at %f before previous begin at %f\n", CI->id, scaffold->id, thisMin.mean, lastMin.mean), hasProblems++;
    if (lastMin.mean > thisMax.mean + 1.0)
      fprintf(stderr, "ScaffoldSanity()--  contig %d in scaffold %d -- ends at %f before previous begin at %f\n", CI->id, scaffold->id, thisMax.mean, lastMin.mean), hasProblems++;

    //  If this contig is contained in the last, all bets are off on variance.
    //  Thse occur when placing rocks during initial scaffolding (we shouldn't be calling this yet, though).
    //
    if (lastMax.mean < thisMin.mean)
      //  An actual gap -- not a -20 gap or overlapping contigs.
      if (lastMax.variance > thisMin.variance)
        fprintf(stderr, "ScaffoldSanity()--  contig %d in scaffold %d -- negative gap variance %f on positive gap size %f\n", CI->id, scaffold->id, thisMin.variance - lastMax.variance, thisMin.mean - lastMax.mean), hasProblems++;

    lastMin = thisMin;
    lastMax = thisMax;
  }

  // ACTIVATE THE DEFINE IF YOU WOULD RATHER HAVE A BOGUS SCAFFOLD THAN NO ASSEMBLY AT ALL.
  // WHEN THIS IS SET, REMEMBER TO SEARCH LOGS FOR "papered over".
  //#define PAPER_OVER_PROBLEM_SCAFFOLD 
#ifdef PAPER_OVER_PROBLEM_SCAFFOLD
  if (hasProblems > 0) {
    // If this approach becomes necessary, add flag to scaffold struct to mark it bad.
    fprintf(stderr, "ScaffoldSanity()-- scaffold %d has "F_U32" problems.\n", scaffold->id, hasProblems);
    InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &CIs);
    CI = NextCIScaffoldTIterator(&CIs);
    CDS_CID_t  prevCIid     = CI->id;
    LengthT    thisLeftEnd  = { 0, 0 };
    LengthT    thisRightEnd = { 0, 0 };
    bool A_before_B = (CI->offsetAEnd.mean < CI->offsetBEnd.mean);
    LengthT    prevLeftEnd  = (A_before_B) ? CI->offsetAEnd : CI->offsetBEnd;
    LengthT    prevRightEnd = (A_before_B) ? CI->offsetBEnd : CI->offsetAEnd;
    int papered_offset = 0;
    int papered_variance = 0;
    int contig_length = 0;
    int papered_gap_length = 100;  // As we paper over bad scaffold, every gap gets a uniform arbitrary size
    while (NULL != (CI = NextCIScaffoldTIterator(&CIs))) {
      A_before_B = (CI->offsetAEnd.mean < CI->offsetBEnd.mean);
      thisLeftEnd  = (A_before_B) ? CI->offsetAEnd : CI->offsetBEnd;
      thisRightEnd = (A_before_B) ? CI->offsetBEnd : CI->offsetAEnd;
      fprintf(stderr, "ScaffoldSanity()-- Contig %8d at %9.0f +- %-11.0f to %9.0f +- %-11.0f  ctg len %9.0f  gap to next %9.0f +- %-11.0f\n",
              prevCIid,
              prevLeftEnd.mean,  prevLeftEnd.variance,
              prevRightEnd.mean, prevRightEnd.variance,
              prevRightEnd.mean - prevLeftEnd.mean,
              thisLeftEnd.mean - prevRightEnd.mean,
              thisLeftEnd.variance - prevRightEnd.variance);
      // Paper over the problem by placing the contig downstream of the previous contig
      contig_length = abs(thisRightEnd.mean - thisLeftEnd.mean);
      if (contig_length < 1) contig_length = 1;  // enforce a minimum
      if (A_before_B) {
	CI->offsetAEnd.mean = papered_offset;
	papered_offset += contig_length;  
	CI->offsetBEnd.mean = papered_offset;
	papered_offset += papered_gap_length; 
      } else {
	CI->offsetBEnd.mean = papered_offset;
	papered_offset += contig_length;  
	CI->offsetAEnd.mean = papered_offset;
	papered_offset += papered_gap_length; 
      } 
      // Paper over the variance too. Arbitrarily set variance same as mean
      CI->offsetAEnd.variance = CI->offsetAEnd.mean;
      CI->offsetBEnd.variance = CI->offsetBEnd.mean;
      // will we need to adjust scaffold->bpLength.mean ?
      // On to the next: this contig becomes previous contig
      prevCIid      = CI->id;
      prevLeftEnd   = thisLeftEnd;
      prevRightEnd  = thisRightEnd;
    }
    // End case: print last contig without gap
    fprintf(stderr, "ScaffoldSanity()-- Contig %8d at %9.0f +- %-11.0f to %9.0f +- %-11.0f  ctg len %9.0f\n",
            prevCIid,
            prevLeftEnd.mean,  prevLeftEnd.variance,
            prevRightEnd.mean, prevRightEnd.variance,
            prevRightEnd.mean - prevLeftEnd.mean);
    fprintf(stderr, "ScaffoldSanity()-- papered over problems in scaffold %d.\n", scaffold->id);
  }
#else
  assert (hasProblems == 0);
#endif

  assert(numCtg == scaffold->info.Scaffold.numElements);
  assert(scaffold->bpLength.mean     >= 0);
  assert(scaffold->bpLength.variance >= 0);

  // CheckScaffoldOrder checks whether all the ahangs in a scaffold are
  // positive (ie, that the contigs are ordered by their distance from the
  // A end of the scaffold)
  //
  // CheckScaffoldOrder(scaffold, graph);
}


void
ScaffoldSanity(ScaffoldGraphT *graph) {
  GraphNodeIterator   scaffolds;
  CIScaffoldT        *scaffold;

  InitGraphNodeIterator(&scaffolds, graph->ScaffoldGraph, GRAPH_NODE_DEFAULT);
  while ((scaffold = NextGraphNodeIterator(&scaffolds)) != NULL)
    ScaffoldSanity(graph, scaffold);
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

    if (contig->info.Contig.numCI > 1)
       continue;

    CI = GetGraphNode(ScaffoldGraph->CIGraph, contig->info.Contig.AEndCI);

    numSingletonScaffolds++;

    // if we are forced marked unique and we are not allowed to be demoted, continue
    if (ScaffoldGraph->tigStore->getUnitigFUR(CI->id) == AS_FORCED_UNIQUE &&
        GlobalData->allowDemoteMarkedUnitigs == FALSE) {
       continue;
    }

    if ((ScaffoldGraph->tigStore->getUnitigFUR(CI->id) != AS_FORCED_REPEAT &&
         ScaffoldGraph->tigStore->getUnitigCoverageStat(CI->id) > GlobalData->cgbDefinitelyUniqueCutoff) ||
         (CI->bpLength.mean > 2000.0))
       continue;

    // We've found a victim!!!

    numDemoted++;

    fprintf(stderr,
             "** Demoting Contig/Unitig "F_CID"/"F_CID" with coverage stat %d length %.0f scaffold "F_CID"\n",
             contig->id, CI->id, ScaffoldGraph->tigStore->getUnitigCoverageStat(CI->id), scaffold->bpLength.mean, scaffold->id);
    // Mark the Underlying Unitig as un-scaffolded, and not-unique
    SetNodeType(CI, UNRESOLVEDCHUNK_CGW);
  
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
  }

  //  If we removed any scaffolds, rebuild all the edges.
  //
  if (numDemoted > 0) {
    BuildSEdges(ScaffoldGraph, FALSE, TRUE);
    MergeAllGraphEdges(ScaffoldGraph->ScaffoldGraph, TRUE, TRUE);
  }

  fprintf(stderr,
          "# Considered %d scaffolds of which %d were single and %d (%.0f%%) were demoted\n",
          numScaffolds, numSingletonScaffolds, numDemoted,
          (numSingletonScaffolds > 0? ((double)(numDemoted)/(double)(numSingletonScaffolds)): 0.0));
}
