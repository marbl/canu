
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
static char *rcsid = "$Id: CIScaffoldT_CGW.c,v 1.59 2012-08-02 21:56:32 brianwalenz Exp $";

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
  ChunkInstanceT *CI;
  CIScaffoldTIterator CIs;
  SEdgeTIterator SEdges;
  SEdgeT *SEdge;

  PrintCIScaffoldHeader(stream, graph, scaffold);
  fprintf(stream, "> Includes CIs\n");
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

  //fprintf(stderr, "InsertCIInScaffold()--  Insert CI %d into scaffold %d at offset %.0f,%.0f orient %d\n",
  //        ci, sid, aEndOffset.mean, bEndOffset.mean, AEndToBend);

  //  Check for bad variances on ends -- hit in gos bunches of times
  //
  //*** Variance Fixup Alert ***
  //chunkInstance->id: 9653113
  //aEndOffset.mean: 7348.000000, bEndOffset.mean: 8305.995965
  //aEndOffset.variance: 2230432.164665, bEndOffset.variance: 2230432.164665
  //chunkInstance->bpLength.mean: 958.000000
  //chunkInstance->bpLength.variance: 24.908000
  //
  //chunkInstance->offsetAEnd.mean: 0.000000, offsetAEnd.variance: 0.000000
  //chunkInstance->offsetBEnd.mean: 958.000000, offsetBEnd.variance: 24.908000
  //chunkInstance->bpLength.mean: 958.000000, bpLength.variance: 24.908000
  //
  //
  if (fabs(aEndOffset.variance - bEndOffset.variance) > 1.0 + chunkInstance->bpLength.variance) {
    fprintf(stderr, "InsertCIInScaffold()-- VARIANCE FIXUP ALERT for CI %d length %.0f variance %.0f\n",
            chunkInstance->id,chunkInstance->bpLength.mean, chunkInstance->bpLength.variance);
    fprintf(stderr, "                       supplied placement has larger variance (%.0f) than 1.0 + chunk (%.0f).\n",
            fabs(aEndOffset.variance - bEndOffset.variance), 1.0 + chunkInstance->bpLength.variance);
    fprintf(stderr, "                       aEndOffset mean %12.0f variance %12.0f\n", aEndOffset.mean, aEndOffset.variance);
    fprintf(stderr, "                       bEndOffset mean %12.0f variance %12.0f\n", bEndOffset.mean, bEndOffset.variance);

    //  If we're too far off assert
    //assert(fabs(aEndOffset.variance - bEndOffset.variance) >= 0.5 * chunkInstance->bpLength.variance);

    //  Fix the ends.
    if (aEndOffset.mean < bEndOffset.mean)
      bEndOffset.variance = aEndOffset.variance + chunkInstance->bpLength.variance;
    else
      aEndOffset.variance = bEndOffset.variance + chunkInstance->bpLength.variance;
  }

  //  Check that the placement is more or less the same as the length.
  //
  if (fabs(aEndOffset.mean - bEndOffset.mean) >= (1.2 * chunkInstance->bpLength.mean)) {
    fprintf(stderr, "InsertCIInScaffold()-- LENGTH FIXUP ALERT for CI %d length %.0f variance %.0f\n",
            chunkInstance->id,chunkInstance->bpLength.mean, chunkInstance->bpLength.variance);
    fprintf(stderr, "                       supplied placement has larger length (%.0f) than 1.2 * chunk (%.0f).\n",
            fabs(aEndOffset.mean - bEndOffset.mean), 1.2 * chunkInstance->bpLength.mean);
    fprintf(stderr, "                       aEndOffset mean %12.0f variance %12.0f\n", aEndOffset.mean, aEndOffset.variance);
    fprintf(stderr, "                       bEndOffset mean %12.0f variance %12.0f\n", bEndOffset.mean, bEndOffset.variance);

    //  This is probably a bad thing to do...
    //
    //  The one time BPW has seen this was from AS_REZ/Update_Scaffold_Graph() (the second call) on
    //  nasonia vitripennis, 6x sanger, ~300Mbp genome, reads from TraceDB, with little QC on the
    //  reads.
    //
    //SplitUnresolvedCI()--  Cloned surrogate ma of CI 161408 (length 2226) has length (2226) 2226
    //SplitUnresolvedCI()--  Split base CI=161408 (now with 1 instances) into new CI=365714
    //SplitUnresolvedContig()-- Split base CI=161408 contig=161408 into new CI=365714 contig=386443
    //InsertCIInScaffold()-- VARIANCE FIXUP ALERT for CI 386443 length 2226 variance 58
    //                       supplied placement has larger variance (640) than 1.0 + chunk (59).
    //                       aEndOffset mean       275483 variance     35086224
    //                       bEndOffset mean       300090 variance     35086864
    //InsertCIInScaffold()-- LENGTH FIXUP ALERT for CI 386443 length 2226 variance 58
    //                       supplied placement has larger length (24607) than 1.2 * chunk (2671).
    //                       aEndOffset mean       275483 variance     35086224
    //                       bEndOffset mean       300090 variance     35086282
    //
    //  And once in GOSIII
    //
    //YIKES!!  old_len = 342  new_len = 415 for chunk 9959948
    //InsertCIInScaffold()-- LENGTH FIXUP ALERT for CI 9959948 length 342 variance 9
    //                       supplied placement has larger length (415) than 1.2 * chunk (410).
    //                       aEndOffset mean            0 variance           36
    //                       bEndOffset mean          415 variance           45
    //
#if 0
    if (aEndOffset.mean < bEndOffset.mean)
      bEndOffset.mean = aEndOffset.mean + chunkInstance->bpLength.mean;
    else
      aEndOffset.mean = bEndOffset.mean + chunkInstance->bpLength.mean;
#endif
  }
  //assert(fabs(aEndOffset.mean - bEndOffset.mean) < (1.2 * chunkInstance->bpLength.mean));

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

      SetCIScaffoldTLength(ScaffoldGraph, ciScaffold, TRUE);

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
      SetCIScaffoldTLength(ScaffoldGraph, ciScaffold, TRUE);
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
    //    CheckCIScaffoldTLength(sgraph, ciScaffold);
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




/***************************************************************************/
int IsScaffoldInternallyConnected(ScaffoldGraphT *sgraph,
                                  CIScaffoldT *scaffold, int32 edgeTypes) {
  //
  // returns the number of connected components of the <scaffold>
  // NOTE: it considers ONLY trusted edges
  // Will modify the setId field of the NodeCGW_T structure to reflect
  // which component a node belongs to.
  //
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
    GraphEdgeIterator   edges(sgraph->ContigGraph, chunk->id, ALL_END, edgeTypes);
    CIEdgeT            *edge;

    assert(chunk->setID >= 0);

    while ((edge = edges.nextMerged()) != NULL) {
      //
      // get the other end
      //
      ChunkInstanceT
        * otherChunk = GetGraphNode(sgraph->ContigGraph,
                                    (chunk->id == edge->idA) ?
                                    edge->idB : edge->idA);
      int32 weight = edge->edgesContributing - (isOverlapEdge(edge));
      assert(otherChunk != NULL);

      // See each edge only once
      if(chunk->id != edge->idA)
        continue;

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


/***********************************************************************/
int32 CheckScaffoldConnectivityAndSplit(ScaffoldGraphT *graph, CDS_CID_t scaffoldID, int32 edgeTypes, int verbose){
  CIScaffoldT  *scaffold      = GetCIScaffoldT(graph->CIScaffolds, scaffoldID);
  int           numComponents = IsScaffoldInternallyConnected(graph, scaffold, edgeTypes);
  int32         numNodes      = scaffold->info.Scaffold.numElements;

  // Expected case, Scaffold is connected
  if(numComponents > 1){
    CDS_CID_t *nodes = (CDS_CID_t *)safe_malloc(sizeof(CDS_CID_t) * numNodes);
    int inode;

    // IsScaffoldInternalyConnected does a connected component analysis, marking the contigs with their component number
    // the following code leverages this marking to break up the scaffold.

    int component;
    int  *nodesEnd;
    NodeCGW_T *thisNode;
    CIScaffoldTIterator scaffoldNodes;

    fprintf(stderr, "WARNING! Scaffold "F_CID" is not connected has %d components\n", scaffoldID, numComponents);
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
      LengthT firstOffset = {0.0, 0.0};
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
#if 0
      CIScaffold.aEndCoord = CIScaffold.bEndCoord = -1;
#endif
      CIScaffold.numEssentialA = CIScaffold.numEssentialB = 0;
      CIScaffold.essentialEdgeB = CIScaffold.essentialEdgeA = NULLINDEX;
      AppendGraphNode(graph->ScaffoldGraph, &CIScaffold);

      scaffold = GetCIScaffoldT(graph->CIScaffolds, scaffoldID);

      for(inode = 0, seenFirstOffset = FALSE; inode < numNodes; inode++){
        NodeCGW_T *thisNode = GetGraphNode(graph->ContigGraph, nodes[inode]);
        if(thisNode->setID == component){
          LengthT offsetAEnd, offsetBEnd;
          if(!seenFirstOffset){
            if(GetNodeOrient(thisNode).isForward()) {
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

    safe_free(nodes);
  }
  return numComponents;
}


/*****************************************************************************/

void CheckTrustedEdges(ScaffoldGraphT * sgraph,  CDS_CID_t cid) {
  //
  // iterates over all trusted edges of CI cid and check if
  // the other end is in the same scaffold
  //
  GraphEdgeIterator  edges(sgraph->ContigGraph, cid, ALL_END, ALL_TRUSTED_EDGES);
  CIEdgeT * edge;
  CDS_CID_t next;
  ChunkInstanceT * next_chunk;
  ChunkInstanceT * this_chunk = GetGraphNode(sgraph->ContigGraph, cid);
  CDS_CID_t sid = this_chunk->scaffoldID;

  while((edge = edges.nextMerged()) != NULL){
    assert(edge != NULL);
    //
    // get the other end
    //
    if (cid == edge->idA)
      next = edge->idB;
    else
      //      continue;  // avoid double checking of (i,j) and (j,i)
      next = edge->idA;

    next_chunk = GetGraphNode(ScaffoldGraph->ContigGraph, next);
    assert(next_chunk != NULL);

    if (next_chunk->scaffoldID != sid)
#if 1
      fprintf(stderr,"-=> BAD edge id:"F_CID" "F_CID"("F_CID")->"F_CID"("F_CID") (weight %d, status %d)\n",
              (CDS_CID_t) GetVAIndex_CIEdgeT(sgraph->ContigGraph->edges, edge),
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

  InitGraphNodeIterator(&nodes, sgraph->ContigGraph, GRAPH_NODE_DEFAULT);
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
  GraphEdgeIterator edges(sgraph->ContigGraph, cid, ALL_END, ALL_EDGES);
  EdgeCGW_T * edge;
  CDS_CID_t next;
  NodeCGW_T * next_chunk;
  int out_of_sid_links = 0;

  while((edge = edges.nextMerged()) != NULL){
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

    next_chunk = GetGraphNode(ScaffoldGraph->ContigGraph, next);
    assert(next_chunk != NULL);

    if ((next_chunk->scaffoldID != sid) && (next_chunk->scaffoldID != -1)) {
      out_of_sid_links++;
      /*** mjf ***/
      fprintf(stderr, "in CheckAllEdges -=> BAD edge id:"F_CID" "F_CID"("F_CID")->"F_CID"("F_CID") (weight %d, status %d)\n",
              (CDS_CID_t) GetVAIndex_EdgeCGW_T(sgraph->ContigGraph->edges, edge),
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


void
SetCIScaffoldTLength(ScaffoldGraphT *sgraph, CIScaffoldT *scaffold, int32 verbose) {
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
    fprintf(stderr, "SetCIScaffoldTLength()-- adjusted scaffold "F_CID" from (%g,%g) to (%g,%g)\n",
            scaffold->id,
            scaffold->bpLength.mean, scaffold->bpLength.variance,
            maxOffset.mean,          maxOffset.variance);

  scaffold->bpLength = maxOffset;
}


void
SetCIScaffoldTLengths(ScaffoldGraphT *sgraph, int verbose){
  GraphNodeIterator scaffolds;
  CIScaffoldT      *scaffold;

  InitGraphNodeIterator(&scaffolds, sgraph->ScaffoldGraph, GRAPH_NODE_DEFAULT);
  while ((scaffold = NextGraphNodeIterator(&scaffolds)) != NULL)
    SetCIScaffoldTLength(sgraph, scaffold, verbose);
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
        fprintf(stderr,"*# First Chunk "F_CID" in scaffold "F_CID" is not positioned at scaffold start, but at (%g,%g)...fixing\n",
                chunk->id, scaffold->id, chunk->offsetAEnd.mean, chunk->offsetBEnd.mean);
#ifdef STRICT_SCAFFOLD_CHECKING
        assert(0);
#endif
        chunk->offsetAEnd.mean = chunk->offsetAEnd.variance = 0.0;
        chunk->offsetBEnd = chunk->bpLength;
      }
    }
    /*	fprintf(stderr,"* Chunk "F_CID" [%g,%g] mean = %g\n",
        chunk->id,
        chunk->offsetAEnd.mean,
        chunk->offsetBEnd.mean,
        mean);
    */
    if(chunk->offsetAEnd.mean > scaffold->bpLength.mean ||
       chunk->offsetBEnd.mean > scaffold->bpLength.mean){

      fprintf(stderr,"*# Chunk "F_CID" has end point [%g,%g] outside of scaffold length %g\n",
              chunk->id,
              chunk->offsetAEnd.mean,
              chunk->offsetBEnd.mean,
              scaffold->bpLength.mean);

    }

    if(chunk->offsetAEnd.mean < 0 || chunk->offsetBEnd.mean < 0 ||
       chunk->offsetAEnd.variance < 0 || chunk->offsetBEnd.variance < 0){
      fprintf(stderr,"*# Chunk "F_CID" in scaffold "F_CID" at odd position: (%g,%g) (%g,%g)\n",
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
      fprintf(stderr,"*# A Scaffold "F_CID" bpLength TOO Short (%g < %g: %g)...variances (%g, %g: %g)corrected\n",
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
      fprintf(stderr,"*# B Scaffold "F_CID" bpLength TOO Short (%g < %g:%g)...variances (%g, %g:%g) corrected\n",
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
    //fprintf(stderr, " Early end to CheckCIScaffoldT() for scaffold "F_CID": numElements ==1\n",scaffold->id);
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

      /*	fprintf(stderr,"* Chunk "F_CID" [%g,%g] mean = %g\n",
		chunk->id,
		chunk->offsetAEnd.mean,
		chunk->offsetBEnd.mean,
		mean);
      */
      if(chunk->offsetAEnd.mean > scaffold->bpLength.mean ||
         chunk->offsetBEnd.mean > scaffold->bpLength.mean){

        fprintf(stderr,"* Chunk "F_CID" has end point [%g,%g] outside of scaffold length %g\n",
                chunk->id,
                chunk->offsetAEnd.mean,
                chunk->offsetBEnd.mean,
                scaffold->bpLength.mean);

      }

      /* Check for Chimeric Scaffolds */
#if 0
      if(prevChunk &&
         prevChunk->aEndCoord >= 0 &&
         chunk->aEndCoord >=0){
        int32 calcDiff = (int32) fabs(chunk->offsetAEnd.mean - prevChunk->offsetAEnd.mean );
        int32 realDiff = abs(chunk->aEndCoord - prevChunk->aEndCoord);

        if(realDiff > 250000 || ((realDiff > calcDiff) && (realDiff - calcDiff > 100000))){
          fprintf(stderr,"*** Scaffold "F_CID" is CHIMERIC at point between CIs "F_CID" and "F_CID"\n",
                  sid, prevChunk->id, chunk->id);
          fprintf(stderr,"*** CI "F_CID" ["F_S32","F_S32"]  CI "F_CID" ["F_S32","F_S32"] scaffold gap is "F_S32"\n",
                  prevChunk->id,
                  prevChunk->aEndCoord, prevChunk->bEndCoord,
                  chunk->id,
                  chunk->aEndCoord, chunk->bEndCoord,
                  calcDiff);
        }
      }
#endif
      prevChunk = chunk;

      if(chunk->offsetAEnd.mean < 0 || chunk->offsetBEnd.mean < 0 ||
         chunk->offsetAEnd.variance < 0 || chunk->offsetBEnd.variance < 0){
        fprintf(stderr,"* Chunk "F_CID" in scaffold "F_CID" at odd position: (%g,%g) (%g,%g) fixing\n",
                chunk->id, scaffold->id,
                chunk->offsetAEnd.mean, chunk->offsetAEnd.variance,
                chunk->offsetBEnd.mean, chunk->offsetBEnd.variance);
        chunk->offsetAEnd.mean = chunk->offsetAEnd.variance = 0.0;
        chunk->offsetBEnd = chunk->bpLength;

      }

      if(chunk->offsetAEnd.mean > scaffold->bpLength.mean || chunk->offsetAEnd.variance > scaffold->bpLength.variance){
        fprintf(stderr,"* A Scaffold "F_CID" bpLength TOO Short (%g < %g)...variances (%g, %g)corrected\n",
                scaffold->id, scaffold->bpLength.mean, chunk->offsetAEnd.mean, scaffold->bpLength.variance, chunk->offsetAEnd.variance);
        scaffold->bpLength.mean = chunk->offsetAEnd.mean;
        scaffold->bpLength.variance = chunk->offsetAEnd.variance;
      }
      if(chunk->offsetBEnd.mean > scaffold->bpLength.mean || chunk->offsetBEnd.variance > scaffold->bpLength.variance){
        fprintf(stderr,"* B Scaffold "F_CID" bpLength TOO Short (%g < %g)...variances (%g, %g) corrected\n",
                scaffold->id, scaffold->bpLength.mean, chunk->offsetBEnd.mean, scaffold->bpLength.variance, chunk->offsetBEnd.variance);
        scaffold->bpLength.mean = chunk->offsetBEnd.mean;
        scaffold->bpLength.variance = chunk->offsetBEnd.variance;
      }

      if(chunk->offsetAEnd.mean < mean ||
         chunk->offsetBEnd.mean < mean){
        fprintf(stderr,"* Screwed up scaffold "F_CID": Chunk "F_CID" has bad mean\n",
                sid, chunk->id);

        AppendPtrT(chunksToBeRemoved, (void **) &chunk);

      }else{
        // Only update the mean if we are moving along as expected
        mean = MIN(chunk->offsetAEnd.mean, chunk->offsetBEnd.mean);
      }
      numChecked++;
    }
    cgwError = GetNumPtrTs(chunksToBeRemoved);
    if(cgwError){
      int i;
      fprintf(stderr,"* Scaffold "F_CID" iteration %d  %d bad CIs out of %d LSE:%g improvement:%g\n",
              sid, iteration, cgwError, numChecked, LSE, improvement);
      for(i = 0; i < GetNumPtrTs(chunksToBeRemoved); i++){

        chunk = *(ChunkInstanceT **)GetPtrT(chunksToBeRemoved, i);
        fprintf(stderr,"* Reinserting chunk "F_CID" in scaffold "F_CID"\n",
                chunk->id,sid);
        RemoveCIFromScaffold(sgraph, scaffold, chunk, FALSE);
        InsertCIInScaffold(sgraph, chunk->id, sid,
                           chunk->offsetAEnd, chunk->offsetBEnd, TRUE, FALSE);
      }

      if (LeastSquaresGapEstimates(sgraph, scaffold) == false) {
        fprintf(stderr, "LeastSquaresGapEstimates failed for scaffold "F_CID" in CheckScaffolds.\n",
                scaffold->id);
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
    fprintf(stderr,"* Took %d iterations to patch up scaffold "F_CID"\n",
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
      fprintf(stderr,"* Screwed up scaffold "F_CID": Chunk "F_CID" has bad variance\n",
              sid, chunk->id);
    }

    if( mean < chunk->offsetAEnd.mean || mean < chunk->offsetBEnd.mean)
      {
        mean = MAX(chunk->offsetAEnd.mean, chunk->offsetBEnd.mean);
        variance = MAX(chunk->offsetAEnd.variance, chunk->offsetBEnd.variance);
      }
  }
  if(cgwError){
    fprintf(stderr,"* Screwed up scaffold "F_CID"\n", sid);
    DumpCIScaffold(stderr,sgraph, scaffold, FALSE);
    iterationError = TRUE;
  }
  ScaffoldSanity(scaffold, sgraph);


}

void
FixupLengthsScaffoldTs(ScaffoldGraphT *sgraph) {
  GraphNodeIterator scaffolds;
  CIScaffoldT *scaffold;

  InitGraphNodeIterator(&scaffolds, sgraph->ScaffoldGraph, GRAPH_NODE_DEFAULT);
  while ((scaffold = NextGraphNodeIterator(&scaffolds)) != NULL){
    if (scaffold->type != REAL_SCAFFOLD)
      continue;
    FixupLengthScaffoldT(sgraph, scaffold);
  }
}

/* Check means and variance in a single */
void
FixupLengthScaffoldT(ScaffoldGraphT *sgraph, CIScaffoldT *scaffold) {
  CIScaffoldTIterator CIs;
  ChunkInstanceT *CI;
  LengthT minOffset;
  LengthT maxOffset;
  LengthT computedLength;

  minOffset.mean     =  DBL_MAX;
  minOffset.variance =  DBL_MAX;

  maxOffset.mean     = -DBL_MAX;
  maxOffset.variance =  DBL_MAX;

  InitCIScaffoldTIterator(sgraph, scaffold, TRUE,   FALSE, &CIs);
  while ((CI = NextCIScaffoldTIterator(&CIs)) != NULL){

    assert(CI->scaffoldID == scaffold->id);

    if (CI->offsetAEnd.mean > CI->offsetBEnd.mean) {
      if (CI->offsetAEnd.mean > maxOffset.mean)  maxOffset = CI->offsetAEnd;
      if (CI->offsetBEnd.mean < minOffset.mean)  minOffset = CI->offsetBEnd;
    } else {
      if (CI->offsetBEnd.mean > maxOffset.mean)  maxOffset = CI->offsetBEnd;
      if (CI->offsetAEnd.mean < minOffset.mean)  minOffset = CI->offsetAEnd;
    }
  }

  computedLength.mean     = maxOffset.mean     - minOffset.mean;
  computedLength.variance = maxOffset.variance - minOffset.variance;

  assert(computedLength.variance >= 0);
  assert(computedLength.mean     >  0);

  if (computedLength.mean > scaffold->bpLength.mean) {
    fprintf(stderr,"* Adjusting scaffold "F_CID" bplength from %g to %g [%g,%g]\n",
            scaffold->id,
            scaffold->bpLength.mean,
            computedLength.mean,
            minOffset.mean, maxOffset.mean);

    scaffold->bpLength = computedLength;
  }
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
             "** Demoting Contig/Unitig "F_CID"/"F_CID" with coverage stat %d length %g scaffold "F_CID"\n",
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
          "# Considered %d scaffolds of which %d were single and %d (%g%%) were demoted\n",
          numScaffolds, numSingletonScaffolds, numDemoted,
          (numSingletonScaffolds > 0? ((double)(numDemoted)/(double)(numSingletonScaffolds)): 0.0));
}
