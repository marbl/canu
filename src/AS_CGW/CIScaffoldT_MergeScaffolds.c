
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
static char *rcsid = "$Id: CIScaffoldT_MergeScaffolds.c,v 1.1 2012-05-01 04:18:19 brianwalenz Exp $";

#include "AS_global.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"
#include "ChiSquareTest_CGW.h"
#include "InterleavedMerging.h"


//  Define this to check (and assert) if the graph is not internally
//  connected before recomputing offsets.  It's expensive, and if you
//  already know it's OK (debugging, maybe??) you can skip it.
#define CHECKCONNECTED



static
double
GetVarianceOffset(CIScaffoldT * scaffold, double meanOffset, int isAB) {
  double lengthDelta = (isAB) ? 0.0 : scaffold->bpLength.mean;
  double varDelta    = (isAB) ? 0.0 : scaffold->bpLength.variance;

  CIScaffoldTIterator  CIs;
  ChunkInstanceT      *CI;

  LengthT lowerOffset  = {-1, 0};
  LengthT higherOffset = { 0, 0};

  InitCIScaffoldTIterator(ScaffoldGraph, scaffold, isAB, FALSE, &CIs);

  while ((CI = NextCIScaffoldTIterator(&CIs)) != NULL) {
    double aEnd = fabs(lengthDelta - CI->offsetAEnd.mean);
    double bEnd = fabs(lengthDelta - CI->offsetBEnd.mean);

    if (aEnd > meanOffset || bEnd > meanOffset) {
      if (aEnd < meanOffset) {
        lowerOffset  = CI->offsetAEnd;
        higherOffset = CI->offsetBEnd;
        break;
      } else if (bEnd < meanOffset) {
        lowerOffset  = CI->offsetBEnd;
        higherOffset = CI->offsetAEnd;
        break;
      } else {
        higherOffset = (aEnd < bEnd) ? CI->offsetAEnd : CI->offsetBEnd;
        break;
      }
    }

    lowerOffset = (aEnd < bEnd) ? CI->offsetBEnd : CI->offsetAEnd;
  }

  if (higherOffset.mean == 0)
    return(scaffold->bpLength.variance * meanOffset / scaffold->bpLength.mean);

  double vl = fabs(varDelta - lowerOffset.variance);
  double vh = fabs(varDelta - higherOffset.variance);

  double ll = fabs(lengthDelta - lowerOffset.mean);
  double lh = fabs(lengthDelta - higherOffset.mean);

  return(vl + (vh - vl) * (meanOffset - ll) / (lh - ll));
}




static
void
InsertScaffoldContentsIntoScaffold(ScaffoldGraphT *sgraph,
                                   CDS_CID_t newScaffoldID,
                                   CDS_CID_t oldScaffoldID,
                                   SequenceOrient orient,
                                   LengthT * offset,
                                   int contigNow) {
  CIScaffoldT *oldScaffold = GetCIScaffoldT(sgraph->CIScaffolds, oldScaffoldID);
  CIScaffoldT *newScaffold = GetCIScaffoldT(sgraph->CIScaffolds, newScaffoldID);

  // CheckCIScaffoldTLength(sgraph, oldScaffold);
  // CheckCIScaffoldTLength(sgraph, newScaffold);

  fprintf(stderr,"InsertScaffoldContentsIntoScaffold()-- Insert scaffold "F_CID" (%.0fbp) into scaffold "F_CID" (%.0fbp) at offset %.3f +/- %.3f orient %c\n",
          oldScaffoldID, oldScaffold->bpLength.mean,
          newScaffoldID, newScaffold->bpLength.mean,
          offset->mean, sqrt(offset->variance),
          orient.toLetter());

  assert(offset->mean     >= 0.0);
  assert(offset->variance >= 0.0);

  //  Shouldn't be necessary, occasionally (in GOSIII) the length was wrong.
  SetCIScaffoldTLength(sgraph, newScaffold, FALSE);
  SetCIScaffoldTLength(sgraph, oldScaffold, FALSE);

  CIScaffoldTIterator CIs;
  ChunkInstanceT     *CI;

  InitCIScaffoldTIterator(sgraph, oldScaffold, orient.isForward(), FALSE, &CIs);
  while ((CI = NextCIScaffoldTIterator(&CIs)) != NULL) {
    LengthT offsetAEnd;
    LengthT offsetBEnd;

    //  These should be guaranteed by the SetCIScaffoldTLength() above.
    assert(CI->offsetAEnd.mean <= oldScaffold->bpLength.mean);
    assert(CI->offsetBEnd.mean <= oldScaffold->bpLength.mean);

    //  These guaranteed too.  Previous versions would detect the case and simply reset the variance
    //  of the placed contig (that's offsetBEnd) to the variance of the 'offset'.  That is:
    //    offsetBEnd.variance = offset->variance
    //    offsetAEnd.variance = offset->variance
    assert(CI->offsetAEnd.variance <= oldScaffold->bpLength.variance);
    assert(CI->offsetBEnd.variance <= oldScaffold->bpLength.variance);

    if (orient.isForward()) {
      offsetAEnd.mean     = offset->mean     + CI->offsetAEnd.mean;
      offsetAEnd.variance = offset->variance + CI->offsetAEnd.variance;
      offsetBEnd.mean     = offset->mean     + CI->offsetBEnd.mean;
      offsetBEnd.variance = offset->variance + CI->offsetBEnd.variance;
    } else {
      offsetAEnd.mean     = offset->mean     + (oldScaffold->bpLength.mean     - CI->offsetAEnd.mean);
      offsetAEnd.variance = offset->variance + (oldScaffold->bpLength.variance - CI->offsetAEnd.variance);
      offsetBEnd.mean     = offset->mean     + (oldScaffold->bpLength.mean     - CI->offsetBEnd.mean);
      offsetBEnd.variance = offset->variance + (oldScaffold->bpLength.variance - CI->offsetBEnd.variance);
    }

    //  Unfortunately, if oldScaffold has screwed up variances, 
    //assert(offsetAEnd.variance >= 0.0);
    //assert(offsetBEnd.variance >= 0.0);

    fprintf(stderr,"InsertScaffoldContentsIntoScaffold()-- Insert CI "F_CID" (%.0fbp) at offset (%.0f,%.0f); was at (%.0f,%.0f)\n",
            CI->id,
            CI->bpLength.mean,
            offsetAEnd.mean,     offsetBEnd.mean,
            CI->offsetAEnd.mean, CI->offsetBEnd.mean);

    InsertCIInScaffold(sgraph, CI->id, newScaffoldID, offsetAEnd, offsetBEnd, TRUE, contigNow);
  }

  //  Shouldn't be necessary, occasionally (in GOSIII) the length was wrong.
  SetCIScaffoldTLength(sgraph, newScaffold, FALSE);

  //  Check all the contig coordinatess in the new scaffold.
#if 1
  LengthT lastMin   = {0, 0};
  LengthT lastMax   = {0, 0};
  int     thisIdx   = 0;
  int     thisFails = 0;

  InitCIScaffoldTIterator(ScaffoldGraph, newScaffold, TRUE, FALSE, &CIs);
  while ((CI = NextCIScaffoldTIterator(&CIs)) != NULL) {

    if ((CI->offsetAEnd.mean < 0.0) || (CI->offsetAEnd.variance < 0.0) ||
        (CI->offsetBEnd.mean < 0.0) || (CI->offsetBEnd.variance < 0.0)) {
      fprintf(stderr, "InsertScaffoldContentsIntoScaffold()-- Negative contig mean or variance\n");
      thisFails++;
    }

    if ((thisIdx++ != 0) && ((CI->offsetAEnd.mean == 0.0) ||
                             (CI->offsetBEnd.mean == 0.0))) {
      fprintf(stderr, "InsertScaffoldContentsIntoScaffold()-- Zero offset of internal contig\n");
      thisFails++;
    }

    LengthT thisMin = (CI->offsetAEnd.mean < CI->offsetBEnd.mean) ? CI->offsetAEnd : CI->offsetBEnd;
    LengthT thisMax = (CI->offsetAEnd.mean > CI->offsetBEnd.mean) ? CI->offsetAEnd : CI->offsetBEnd;

    if (thisMax.mean < lastMin.mean) {
      fprintf(stderr, "InsertScaffoldContentsIntoScaffold()-- Seriously out of order contigs\n");
      thisFails++;
    }

    if (thisMin.mean < lastMin.mean) {
      fprintf(stderr, "InsertScaffoldContentsIntoScaffold()-- Out of order contigs\n");
      thisFails++;
    }

    lastMin = thisMin;
    lastMax = thisMax;
  }

  if (thisFails)
    DumpCIScaffold(stderr, ScaffoldGraph, newScaffold, FALSE);
  //assert(thisFails == 0);
#endif
}





static
void
MarkUnderlyingRawCIEdgeTrusted(ScaffoldGraphT * sgraph, EdgeCGW_T * raw) {
  CIEdgeT *ciedge      = GetGraphEdge(sgraph->ContigGraph, raw->referenceEdge);
  CIEdgeT *topCIEdge   = GetGraphEdge(sgraph->ContigGraph, ciedge->topLevelEdge);

  assert(topCIEdge->idA != NULLINDEX);
  assert(topCIEdge->idB != NULLINDEX);

  ContigT *contigA = GetGraphNode(sgraph->ContigGraph, topCIEdge->idA);
  ContigT *contigB = GetGraphNode(sgraph->ContigGraph, topCIEdge->idB);

  assert(contigA->scaffoldID == contigB->scaffoldID);

  assert(contigA->flags.bits.isUnique);
  assert(contigB->flags.bits.isUnique);

  SetGraphEdgeStatus(sgraph->ContigGraph, topCIEdge, TRUSTED_EDGE_STATUS);

  if (GlobalData->debugLevel > 0)
    fprintf(stderr,"* Marked contig edge " F_CID " (" F_CID "," F_CID ")%c as trusted(inside scaf " F_CID ")\n",
            topCIEdge->topLevelEdge, topCIEdge->idA, topCIEdge->idB, topCIEdge->orient.toLetter(),
            contigA->scaffoldID);
}





int
MergeScaffolds(InterleavingSpec * iSpec, int32 verbose) {
  int mergedSomething = 0;
  GraphNodeIterator scaffolds;
  CIScaffoldT *thisScaffold;
  CDS_CID_t thisScaffoldID; /* The index of the thisScaffold. We have to be careful about thisScaffold since
                               it is a pointer to a an element of the scaffolds array, which may be reallocated
                               during the loop. */
  CDS_CID_t currentSetID = 0;
  int cntScaffold = 1;
  int cnt;

  InitGraphNodeIterator(&scaffolds, ScaffoldGraph->ScaffoldGraph, GRAPH_NODE_DEFAULT);

  while ((thisScaffold = NextGraphNodeIterator(&scaffolds)) != NULL) {
    const LengthT  nullLength = {0.0, 0.0};
    CIScaffoldT   *AendScaffold = NULL;
    CIScaffoldT   *BendScaffold = NULL;
    SEdgeT        *edgeA = NULL;
    SEdgeT        *edgeB = NULL;
    SEdgeT        *edge = NULL;
    CIScaffoldT   *neighbor = NULL;
    CDS_CID_t      neighborID = -1;
    LengthT        currentOffset = nullLength;
    CIScaffoldT    CIScaffold;
    SequenceOrient orientCI;
    CDS_CID_t      newScaffoldID = -1;
    int            numMerged = 0;

    thisScaffoldID = thisScaffold->id;

    //fprintf(stderr,"* Examining scaffold %d " F_CID "\n",
    //        cntScaffold, thisScaffoldID);

    if (thisScaffold->type != REAL_SCAFFOLD) {
      //fprintf(stderr, "Not a REAL_SCAFFOLD\n");
      continue;
    }

    if (thisScaffold->setID != NULLINDEX) {
      // This Scaffold has already been placed in a Scaffold.
      //fprintf(stderr, "Already placed in a scaffold.\n");
      continue;
    }

    if (thisScaffold->essentialEdgeA != NULLINDEX) {
      edgeA        = GetGraphEdge(ScaffoldGraph->ScaffoldGraph, thisScaffold->essentialEdgeA);
      AendScaffold = GetGraphNode(ScaffoldGraph->ScaffoldGraph,
                                  (edgeA->idA == thisScaffold->id) ? edgeA->idB : edgeA->idA);
    }

    if (thisScaffold->essentialEdgeB != NULLINDEX) {
      edgeB        = GetGraphEdge(ScaffoldGraph->ScaffoldGraph, thisScaffold->essentialEdgeB);
      BendScaffold = GetGraphNode(ScaffoldGraph->ScaffoldGraph,
                                  (edgeB->idA == thisScaffold->id) ? edgeB->idB : edgeB->idA);
    }

    if ((AendScaffold != NULL) &&
        (BendScaffold != NULL)) {
      // This CI is not a starting point for a Scaffold.
      //fprintf(stderr, "CI not a starting point.\n");
      continue;
    }

    if        (BendScaffold != NULL) {
      orientCI.setIsForward();
      edge       = edgeB;
      neighbor   = BendScaffold;
      neighborID = neighbor->id;

    } else if (AendScaffold != NULL) {
      orientCI.setIsReverse();
      edge       = edgeA;
      neighbor   = AendScaffold;
      neighborID = neighbor->id;

    } else {
      // Singleton Scaffold
      //fprintf(stderr, "singleton scaffold.\n");
      continue;
    }

    //  This is guarding against a failure in GOS III; the very first edge in a scaffold merge
    //  had negative variance, which would trigger an assert in the while() loop below.  In general,
    //  we should be checking all edges used in the merge, not just the first edge.
    //
    if (edge->distance.variance <= 0.0)
      continue;

    //  This is our last chance to abort before we create a new scaffold!

    InitializeScaffold(&CIScaffold, REAL_SCAFFOLD);

    CIScaffold.info.Scaffold.AEndCI = NULLINDEX;
    CIScaffold.info.Scaffold.BEndCI = NULLINDEX;
    CIScaffold.info.Scaffold.numElements = 0;
    CIScaffold.edgeHead = NULLINDEX;
    CIScaffold.bpLength = nullLength;
    thisScaffold->setID = currentSetID;
    newScaffoldID = CIScaffold.id = GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph);

    CIScaffold.flags.bits.isDead = FALSE;
    CIScaffold.numEssentialA = CIScaffold.numEssentialB = 0;
    CIScaffold.essentialEdgeB = CIScaffold.essentialEdgeA = NULLINDEX;
    CIScaffold.setID = NULLINDEX;
    thisScaffold->flags.bits.isDead = TRUE;  // Mark the old scaffold dead

    AppendGraphNode(ScaffoldGraph->ScaffoldGraph, &CIScaffold);  /* Potential realloc of ScaffoldGraph->ScaffoldGraph->nodes */

    thisScaffold = GetGraphNode(ScaffoldGraph->ScaffoldGraph, thisScaffoldID);
    neighbor     = GetGraphNode(ScaffoldGraph->ScaffoldGraph, neighborID);

    if (verbose) {
      fprintf(stderr,"* START: Inserting scaffold " F_CID " into scaffold " F_CID "\n", thisScaffold->id, newScaffoldID);
    }
    assert(thisScaffold->bpLength.variance >= 0);
    assert(neighbor->bpLength.variance >= 0);

    cnt = 0;

    //AppendVA_CDS_CID_t(deadScaffoldIDs, &(thisScaffold->id));

    /* Potential realloc of ScaffoldGraph->ScaffoldGraph->nodes */

    //  BPW - being the first merge, all we're doing is copying the
    //  original scaffold into the new scaffold.
    //
    InsertScaffoldContentsIntoScaffold(ScaffoldGraph,
                                       newScaffoldID, thisScaffold->id,
                                       orientCI, &currentOffset,
                                       iSpec->contigNow);

    thisScaffold = GetGraphNode(ScaffoldGraph->ScaffoldGraph, thisScaffoldID);
    neighbor     = GetGraphNode(ScaffoldGraph->ScaffoldGraph, neighborID);

    currentOffset = thisScaffold->bpLength;
    numMerged = 1;
    mergedSomething++;
    while (neighbor != (CIScaffoldT *)NULL) {
      PairOrient edgeOrient = GetEdgeOrientationWRT(edge, thisScaffold->id);

      cnt++;

      assert(edge->distance.variance >= 0);

      assert(currentOffset.mean >= 0);
      assert(currentOffset.variance >= 0);

      assert(thisScaffold->bpLength.variance >= 0);
      assert(neighbor->bpLength.variance >= 0);

      neighborID = neighbor->id;

      assert(orientCI.isUnknown() == false);

      if (orientCI.isForward()) {
        if (edgeOrient.isAB_AB()) {
          orientCI.setIsForward();
        }else if (edgeOrient.isAB_BA()) {
          orientCI.setIsReverse();
        } else {
          assert(0);
        }
      } else {
        if (edgeOrient.isBA_AB()) {
          orientCI.setIsForward();
        }else if (edgeOrient.isBA_BA()) {
          orientCI.setIsReverse();
        } else {
          assert(0);
        }
      }
      thisScaffold = neighbor;
      thisScaffoldID = thisScaffold->id;

      assert(thisScaffold->bpLength.variance >= 0);
      assert(edge->distance.variance >= 0);
      assert(currentOffset.mean >= 0);
      assert(currentOffset.variance >= 0);
      assert(orientCI.isUnknown() == false);

      currentOffset.mean += edge->distance.mean;
      currentOffset.variance += edge->distance.variance;

      if (edge->distance.mean < -CGW_MISSED_OVERLAP) {
        CIScaffoldT * newScaffold =
          GetCIScaffoldT(ScaffoldGraph->CIScaffolds, newScaffoldID);
        if (currentOffset.mean < 0.0) {
          /*
            0
            newScaffold:                --------------->
            thisScaffold:      ------------------
          */
          CIScaffoldTIterator CIs;
          ChunkInstanceT * CI;
          double variance = GetVarianceOffset(thisScaffold, -currentOffset.mean, orientCI.isForward());

          fprintf(stderr, "Adjusting newScaffold offsets by mean - %f variance + %f\n", currentOffset.mean, variance);

          InitCIScaffoldTIterator(ScaffoldGraph, newScaffold,
                                  TRUE, FALSE, &CIs);
          while ((CI = NextCIScaffoldTIterator(&CIs)) != NULL) {
            CI->offsetAEnd.mean -= currentOffset.mean;
            CI->offsetAEnd.variance += variance;
            CI->offsetBEnd.mean -= currentOffset.mean;
            CI->offsetBEnd.variance += variance;
          }
          currentOffset.mean = currentOffset.variance = 0.0;
        } else {
          /*
            0
            newScaffold:     --------------->
            thisScaffold:              ------------------
          */
          currentOffset.variance = GetVarianceOffset(newScaffold, currentOffset.mean, TRUE);
        }
      }
      assert(currentOffset.variance >= 0);

      if (verbose) {
        fprintf(stderr,"* Adding to scaffold " F_CID " scaffold " F_CID " at orient %c offset %g\n",
                newScaffoldID, thisScaffoldID, orientCI.toLetter(), currentOffset.mean);
      }
#ifdef CHECK_CONTIG_ORDERS
      AddScaffoldToContigOrientChecker(ScaffoldGraph, thisScaffold, coc);
#endif

      //AppendVA_CDS_CID_t(deadScaffoldIDs, &(thisScaffold->id));

      InsertScaffoldContentsIntoScaffold(ScaffoldGraph,
                                         newScaffoldID, thisScaffold->id,
                                         orientCI, &currentOffset,
                                         iSpec->contigNow);

      if (iSpec->contigNow == TRUE) {
        // remove references to edges of dead contigs from scaffold edge

        if (edge->flags.bits.isRaw == 0) {
          SEdgeT * raw     = edge;
          SEdgeT * lastRaw = edge;

          while ((raw = GetGraphEdge(ScaffoldGraph->ScaffoldGraph, raw->nextRawEdge)) != NULL) {
            // referenceEdge references inducing contig edge
            // topLevelEdge references top contig edge
            CIEdgeT * ciEdge = GetGraphEdge(ScaffoldGraph->ContigGraph, raw->referenceEdge);
            if (ciEdge->idA == NULLINDEX || ciEdge->idB == NULLINDEX) {
              lastRaw->nextRawEdge = raw->nextRawEdge;
              raw = lastRaw;
            }
          }
        } else {
          // do what?
        }


      } else {
        // need to set edge statuses so scaffold is clone-connected
        // all the function does is check asserts

        if (GlobalData->debugLevel > 0)
          fprintf(stderr,"* MarkUnderlyingCIEdgesTrusted on SEdge (" F_CID "," F_CID ")%c nextRaw = " F_CID "\n",
                  edge->idA, edge->idB, edge->orient.toLetter(), edge->nextRawEdge);
        
        if (edge->flags.bits.isRaw) {
          MarkUnderlyingRawCIEdgeTrusted(ScaffoldGraph, edge);

        } else {
          SEdgeT *raw = edge;
          while ((raw = GetGraphEdge(ScaffoldGraph->ScaffoldGraph, raw->nextRawEdge)) != NULL)
            MarkUnderlyingRawCIEdgeTrusted(ScaffoldGraph, raw);
        }

      }



      neighbor = GetGraphNode(ScaffoldGraph->ScaffoldGraph, neighborID);

      assert(orientCI.isUnknown() == false);
      assert(edge->distance.variance >= 0);
      assert(thisScaffold->bpLength.variance >= 0);
      assert(currentOffset.variance >= 0);

      thisScaffold->setID = currentSetID;
      thisScaffold->flags.bits.isDead = TRUE; // Mark the old scaffold dead
      currentOffset.mean += thisScaffold->bpLength.mean;
      currentOffset.variance += thisScaffold->bpLength.variance;

      assert(thisScaffold->bpLength.variance >= 0);
      assert(currentOffset.mean >= 0);
      assert(currentOffset.variance >= 0);

      numMerged++;
      if (orientCI.isForward()) {
        if (thisScaffold->essentialEdgeB != NULLINDEX) {
          edge = GetGraphEdge(ScaffoldGraph->ScaffoldGraph, thisScaffold->essentialEdgeB);
          neighbor = GetGraphNode(ScaffoldGraph->ScaffoldGraph,
                                  (edge->idA == thisScaffold->id) ? edge->idB : edge->idA);
          neighborID = neighbor->id;
        } else {// End of Scaffold
          edge = (SEdgeT *)NULL;
          neighbor = (CIScaffoldT *)NULL;
          neighborID = NULLINDEX;
        }
      } else {// orientCI == B_A
        if (thisScaffold->essentialEdgeA != NULLINDEX) {
          edge = GetGraphEdge(ScaffoldGraph->ScaffoldGraph, thisScaffold->essentialEdgeA);
          neighbor = GetGraphNode(ScaffoldGraph->ScaffoldGraph,
                                  (edge->idA == thisScaffold->id) ? edge->idB : edge->idA);
          neighborID = neighbor->id;
        } else {// End of Scaffold
          edge = (SEdgeT *)NULL;
          neighbor = (CIScaffoldT *)NULL;
          neighborID = NULLINDEX;
        }
      }
    }  //  while (neighbor != (CIScaffoldT *)NULL)

    // New scaffold fully popuated now.

    if (iSpec->contigNow == TRUE &&
        GetGraphNode(ScaffoldGraph->ScaffoldGraph,
                     newScaffoldID)->info.Scaffold.numElements > 1) {

      int32  numScaffoldsBefore = GetNumCIScaffoldTs(ScaffoldGraph->CIScaffolds);
      int32  status             = RECOMPUTE_SINGULAR;

      int32 recomputeIteration = 0;

      while (recomputeIteration < 3 &&
             (status == RECOMPUTE_SINGULAR ||
              status == RECOMPUTE_CONTIGGED_CONTAINMENTS)) {
        // need to make sure scaffold is connected with trusted raw edges

        MarkInternalEdgeStatus(ScaffoldGraph,
                               GetGraphNode(ScaffoldGraph->ScaffoldGraph,
                                            newScaffoldID),
                               PAIRWISECHI2THRESHOLD_CGW,
                               1000.0 * SLOPPY_EDGE_VARIANCE_THRESHHOLD,
                               TRUE, TRUE, 0, TRUE);

#ifdef CHECKCONNECTED
        assert(IsScaffoldInternallyConnected(ScaffoldGraph,
                                             GetGraphNode(ScaffoldGraph->ScaffoldGraph,
                                                          newScaffoldID),
                                             ALL_EDGES));
#endif

        status = RecomputeOffsetsInScaffold(ScaffoldGraph,
                                            newScaffoldID,
                                            TRUE, TRUE, FALSE);
        recomputeIteration++;
      }

      if (status != RECOMPUTE_OK)
        fprintf(stderr, "ReomputeOffsetsInScaffold failed (%d) for scaffold " F_CID " in MergeScaffolds\n", status, newScaffoldID);

      int32  numScaffoldsAfter = GetNumCIScaffoldTs(ScaffoldGraph->CIScaffolds);

      //  If numScaffoldsBefore is not the same as numScaffoldsAfter, then RecomputeOffsetsInScaffold() decided
      //  to split a scaffold we just merged.  This (once) led to an infinite loop in scaffold merging where
      //  a two contig and a one contig scaffold were merged, then ROIS() would unmerge them into the original
      //  layout.

      if (status == RECOMPUTE_FAILED_CONTIG_DELETED)
        fprintf(stderr, "MergeScaffolds()-- WARNING:  ROIS failed, contig deleted.\n");

      if         (numScaffoldsAfter - numScaffoldsBefore == numMerged - 1) {
        fprintf(stderr, "MergeScaffolds()-- WARNING:  merged %d scaffolds into 1, then unmerged into %d new scaffolds.  NO MERGING DONE?!\n",
                numMerged, numScaffoldsAfter - numScaffoldsBefore);
        mergedSomething--;

      } else  if (numScaffoldsAfter != numScaffoldsBefore) {
        fprintf(stderr, "MergeScaffolds()-- WARNING:  merged %d scaffolds into 1, then unmerged into %d new scaffolds.  Partial merging done?\n",
                numMerged, numScaffoldsAfter - numScaffoldsBefore);
      }
    }  //  if (iSpec->contigNow == TRUE && .....)

#warning numLiveScaffolds is probably broken, and is nearly unused
    ScaffoldGraph->numLiveScaffolds += (1 - numMerged);

    currentSetID++;
  }

  return mergedSomething;
}
