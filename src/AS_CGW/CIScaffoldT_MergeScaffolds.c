
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
static char *rcsid = "$Id: CIScaffoldT_MergeScaffolds.c,v 1.14 2012-09-10 10:55:44 brianwalenz Exp $";

#include "CIScaffoldT_MergeScaffolds.h"


#include <vector>
#include <set>

using namespace std;


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
InsertScaffoldContentsIntoScaffold(ScaffoldGraphT           *sgraph,
                                   CDS_CID_t                 newScaffoldID,
                                   CDS_CID_t                 oldScaffoldID,
                                   SequenceOrient            orient,
                                   LengthT                  *offset,
                                   int                       currentSetID,
                                   int                       contigNow) {
  CIScaffoldT *oldScaffold = GetCIScaffoldT(sgraph->CIScaffolds, oldScaffoldID);
  CIScaffoldT *newScaffold = GetCIScaffoldT(sgraph->CIScaffolds, newScaffoldID);

  CIScaffoldTIterator CIs;
  ChunkInstanceT     *CI;

  fprintf(stderr,"InsertScaffoldContentsIntoScaffold()-- Insert scaffold "F_CID" (%.0fbp) into scaffold "F_CID" (%.0fbp) at offset %.3f +/- %.3f orient %c\n",
          oldScaffoldID, oldScaffold->bpLength.mean,
          newScaffoldID, newScaffold->bpLength.mean,
          offset->mean, sqrt(offset->variance),
          orient.toLetter());

  assert(offset->mean     >= 0.0);
  assert(offset->variance >= 0.0);

  //  Shouldn't be necessary, occasionally (in GOSIII) the length was wrong.
  SetCIScaffoldTLength(sgraph, newScaffold);
  SetCIScaffoldTLength(sgraph, oldScaffold);

  //  Mark the old one as dead

  oldScaffold->setID             = currentSetID;
  oldScaffold->flags.bits.isDead = TRUE;

  //  Move contigs.

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

#if 1
    if (offsetAEnd.mean < offsetBEnd.mean)
      fprintf(stderr,"InsertScaffoldContentsIntoScaffold()-- Insert CI %8"F_CIDP" %6.0fbp fwd %9.0f +- %-11.0f %9.0f +- %-11.0f was %9.0f +- %-11.0f %9.0f +- %-11.0f\n",
              CI->id,
              CI->bpLength.mean,
              offsetAEnd.mean, offsetAEnd.variance,
              offsetBEnd.mean, offsetBEnd.variance,
              CI->offsetAEnd.mean, CI->offsetAEnd.variance,
              CI->offsetBEnd.mean, CI->offsetBEnd.variance);
    else
      fprintf(stderr,"InsertScaffoldContentsIntoScaffold()-- Insert CI %8"F_CIDP" %6.0fbp rev %9.0f +- %-11.0f %9.0f +- %-11.0f was %9.0f +- %-11.0f %9.0f +- %-11.0f\n",
              CI->id,
              CI->bpLength.mean,
              offsetBEnd.mean, offsetBEnd.variance,
              offsetAEnd.mean, offsetAEnd.variance,
              CI->offsetBEnd.mean, CI->offsetBEnd.variance,
              CI->offsetAEnd.mean, CI->offsetAEnd.variance);
#endif


    //  Contig now?  If we do, and something contigs, we cannot back out this merge.
#warning NOT CONTIGGING NOW
    InsertCIInScaffold(sgraph, CI->id, newScaffoldID, offsetAEnd, offsetBEnd, TRUE, contigNow);
  }

  //  Blindly adjust variances.
#if 1
  {
  LengthT lastMax   = {0.0, 0.0};

  vector<double>  varOld;
  uint32          gap = 0;

  //  Save the old gap variances.

  InitCIScaffoldTIterator(ScaffoldGraph, newScaffold, TRUE, FALSE, &CIs);
  while ((CI = NextCIScaffoldTIterator(&CIs)) != NULL) {
    LengthT thisMin = (CI->offsetAEnd.mean < CI->offsetBEnd.mean) ? CI->offsetAEnd : CI->offsetBEnd;
    LengthT thisMax = (CI->offsetAEnd.mean > CI->offsetBEnd.mean) ? CI->offsetAEnd : CI->offsetBEnd;

    varOld.push_back(thisMin.variance - lastMax.variance);

    lastMax = thisMax;
  }

  //  If any of those old ones are negative, reset them to something sloppy.

  bool  wasAdjusted = false;

  lastMax.mean     = 0.0;
  lastMax.variance = 0.0;

  InitCIScaffoldTIterator(ScaffoldGraph, newScaffold, TRUE, FALSE, &CIs);
  while ((CI = NextCIScaffoldTIterator(&CIs)) != NULL) {
    LengthT *thisMin = (CI->offsetAEnd.mean < CI->offsetBEnd.mean) ? &CI->offsetAEnd : &CI->offsetBEnd;
    LengthT *thisMax = (CI->offsetAEnd.mean > CI->offsetBEnd.mean) ? &CI->offsetAEnd : &CI->offsetBEnd;

    //double ctgVar = thisMax->variance - thisMin->variance;
    double ctgMea = CI->bpLength.mean;
    double ctgVar = CI->bpLength.variance;

    assert(ctgVar > 0);

    if (varOld[gap] < 0) {
      fprintf(stderr, "RESET neg var for gap %u from %.0f to %.0f\n",
              gap, varOld[gap], SLOPPY_EDGE_VARIANCE_THRESHHOLD);

      wasAdjusted = true;

      thisMin->variance = lastMax.variance + SLOPPY_EDGE_VARIANCE_THRESHHOLD;
      thisMax->variance = lastMax.variance + SLOPPY_EDGE_VARIANCE_THRESHHOLD + ctgVar;

    } else {
      thisMin->variance = lastMax.variance + varOld[gap];
      thisMax->variance = lastMax.variance + varOld[gap] + ctgVar;
    }

    gap++;

    lastMax = *thisMax;
  }

  if (wasAdjusted) {
    lastMax.mean     = 0.0;
    lastMax.variance = 0.0;

    InitCIScaffoldTIterator(ScaffoldGraph, newScaffold, TRUE, FALSE, &CIs);
    while ((CI = NextCIScaffoldTIterator(&CIs)) != NULL) {
      LengthT *thisMin = (CI->offsetAEnd.mean < CI->offsetBEnd.mean) ? &CI->offsetAEnd : &CI->offsetBEnd;
      LengthT *thisMax = (CI->offsetAEnd.mean > CI->offsetBEnd.mean) ? &CI->offsetAEnd : &CI->offsetBEnd;

#if 0
      fprintf(stderr, "VARADJ()-- Contig %8d at %9.0f +- %-11.0f to %9.0f +- %-11.0f  ctg len %9.0f  gap to prev %9.0f +- %-11.0f\n",
              CI->id,
              thisMin->mean, thisMin->variance,
              thisMax->mean, thisMax->variance,
              thisMax->mean - thisMin->mean,
              thisMin->mean     - lastMax.mean,
              thisMin->variance - lastMax.variance);
#endif

      lastMax = *thisMax;
    }
  }

  }
#endif

  //  If we adjust variances above, we need to reset scaffold length.  Just be nice and do it for everyone.
  SetCIScaffoldTLength(sgraph, newScaffold);

  ScaffoldSanity(sgraph, newScaffold);
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
    fprintf(stderr,"* Marked contig edge "F_CID" ("F_CID","F_CID")%c as trusted(inside scaf "F_CID")\n",
            topCIEdge->topLevelEdge, topCIEdge->idA, topCIEdge->idB, topCIEdge->orient.toLetter(),
            contigA->scaffoldID);
}





int
MergeScaffolds(InterleavingSpec               *iSpec,
               set<EdgeCGWLabel_T>            &bEdges) {
  int               mergedSomething = 0;
  GraphNodeIterator scaffolds;
  CIScaffoldT      *thisScaffold;
  CDS_CID_t         thisScaffoldID;
  CDS_CID_t         currentSetID = 0;

  InitGraphNodeIterator(&scaffolds, ScaffoldGraph->ScaffoldGraph, GRAPH_NODE_DEFAULT);

  vector<CDS_CID_t>  oldScaffolds;  //  These scaffolds were merged, and no longer exist
  vector<CDS_CID_t>  newScaffolds;  //  These scaffolds were created

  while ((thisScaffold = NextGraphNodeIterator(&scaffolds)) != NULL) {
    const LengthT  nullLength = {0.0, 0.0};
    CIScaffoldT   *AendScaffold = NULL;
    CIScaffoldT   *BendScaffold = NULL;
    SEdgeT        *edgeA = NULL;
    SEdgeT        *edgeB = NULL;
    SEdgeT        *edge = NULL;
    EdgeCGWLabel_T edgeLabel;
    CIScaffoldT   *neighbor = NULL;
    CDS_CID_t      neighborID = -1;
    LengthT        currentOffset = nullLength;
    SequenceOrient orientCI;
    CDS_CID_t      newScaffoldID = -1;
    int            numMerged = 0;

    thisScaffoldID = thisScaffold->id;

    if (thisScaffold->type != REAL_SCAFFOLD)
      continue;

    if (thisScaffold->setID != NULLINDEX)
      // This Scaffold has already been placed in a Scaffold.
      continue;

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
        (BendScaffold != NULL))
      // This CI is not a starting point for a Scaffold.
      continue;

   if        (BendScaffold != NULL) {
      orientCI.setIsForward();
      edge               = edgeB;
      edgeA              = NULL;
      edgeLabel.idA      = edgeB->idA;
      edgeLabel.idB      = edgeB->idB;
      edgeLabel.orient   = edgeB->orient;
      edgeLabel.distance = edgeB->distance;
      neighbor           = BendScaffold;
      neighborID         = neighbor->id;

    } else if (AendScaffold != NULL) {
      orientCI.setIsReverse();
      edge               = edgeA;
      edgeB              = NULL;
      edgeLabel.idA      = edgeA->idA;
      edgeLabel.idB      = edgeA->idB;
      edgeLabel.orient   = edgeA->orient;
      edgeLabel.distance = edgeA->distance;
      neighbor           = AendScaffold;
      neighborID         = neighbor->id;

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

#if 0
    //fprintf(stderr, "TESTING EDGE: %d <-> %d orient %c distance %.1f +- %.1f\n",
    //        edgeLabel.idA, edgeLabel.idB, edgeLabel.orient.toLetter(), edgeLabel.distance.mean, sqrt(edgeLabel.distance.variance));

    //  This shouldn't happen.  Somehow we evaluated an edge that is bogus.
    //
    if (bEdges.find(edgeLabel) != bEdges.end()) {
      fprintf(stderr, "SKIPPING PREVIOUSLY BAD EDGE: %d <-> %d orient %c distance %.1f +- %.1f\n",
              edgeLabel.idA, edgeLabel.idB, edgeLabel.orient.toLetter(), edgeLabel.distance.mean, sqrt(edgeLabel.distance.variance));
      continue;
    }
    assert(bEdges.find(edgeLabel) == bEdges.end());
#endif

    //  This is our last chance to abort before we create a new scaffold!
    //  Save enough of the old scaffolds so we can reconstruct if it all goes bad.

#ifdef TRACK_MATE_PAIR_TEST
    uint32   edgeID = GetVAIndex_CIEdgeT(graph->ContigGraph->edges, edge);

    matePairTestResult[edgeID].mergeAccepted = true;
#endif

    //  Build a new scaffold.

    {
      CIScaffoldT    ns;

      InitializeScaffold(&ns, REAL_SCAFFOLD);

      ns.info.Scaffold.AEndCI      = NULLINDEX;
      ns.info.Scaffold.BEndCI      = NULLINDEX;
      ns.info.Scaffold.numElements = 0;
      ns.bpLength                  = nullLength;
      ns.id                        = GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph);

      ns.flags.bits.isDead         = FALSE;
      ns.numEssentialA             = 0;
      ns.numEssentialB             = 0;
      ns.essentialEdgeB            = NULLINDEX;
      ns.essentialEdgeA            = NULLINDEX;
      ns.setID                     = NULLINDEX;

      newScaffoldID                = ns.id;

      AppendGraphNode(ScaffoldGraph->ScaffoldGraph, &ns);

      //  Ensure that there are no edges, and that the edgeList is allocated.
      assert(ScaffoldGraph->ScaffoldGraph->edgeLists[ns.id].empty() == true);

      thisScaffold = GetGraphNode(ScaffoldGraph->ScaffoldGraph, thisScaffoldID);
      neighbor     = GetGraphNode(ScaffoldGraph->ScaffoldGraph, neighborID);

      newScaffolds.push_back(newScaffoldID);
    }

    assert(thisScaffold->bpLength.variance >= 0);
    assert(neighbor->bpLength.variance >= 0);

    //fprintf(stderr, "\n\nDEBUG - merging on this edge:\n");
    //if (edgeA)  PrintGraphEdge(stderr, ScaffoldGraph->ScaffoldGraph, "A", edgeA, -1);
    //if (edgeB)  PrintGraphEdge(stderr, ScaffoldGraph->ScaffoldGraph, "B", edgeB, -1);


    //  Being the first merge, all we're doing is copying the original scaffold into the new
    //  scaffold.
    //
    oldScaffolds.push_back(thisScaffold->id);

    InsertScaffoldContentsIntoScaffold(ScaffoldGraph,
                                       newScaffoldID, thisScaffold->id,
                                       orientCI, &currentOffset,
                                       currentSetID,
                                       iSpec->contigNow);

    thisScaffold = GetGraphNode(ScaffoldGraph->ScaffoldGraph, thisScaffoldID);
    neighbor     = GetGraphNode(ScaffoldGraph->ScaffoldGraph, neighborID);

    currentOffset = thisScaffold->bpLength;
    numMerged = 1;
    mergedSomething++;

    while (neighbor != (CIScaffoldT *)NULL) {
      PairOrient edgeOrient = GetEdgeOrientationWRT(edge, thisScaffold->id);

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
      thisScaffold   = neighbor;
      thisScaffoldID = thisScaffold->id;

      assert(thisScaffold->bpLength.variance >= 0);
      assert(edge->distance.variance >= 0);
      assert(currentOffset.mean >= 0);
      assert(currentOffset.variance >= 0);
      assert(orientCI.isUnknown() == false);

      currentOffset.mean     += edge->distance.mean;
      currentOffset.variance += edge->distance.variance;

      if (edge->distance.mean < -CGW_MISSED_OVERLAP) {
        CIScaffoldT * newScaffold = GetCIScaffoldT(ScaffoldGraph->CIScaffolds, newScaffoldID);

        if (currentOffset.mean < 0.0) {
          //  newScaffold:                --------------->
          //  thisScaffold:      ------------------
          CIScaffoldTIterator CIs;
          ChunkInstanceT     *CI;
          double              variance = GetVarianceOffset(thisScaffold, -currentOffset.mean, orientCI.isForward());

          fprintf(stderr, "Shift existing contigs by %.0f +- %.0f\n", -currentOffset.mean, variance);

          InitCIScaffoldTIterator(ScaffoldGraph, newScaffold, TRUE, FALSE, &CIs);

          while ((CI = NextCIScaffoldTIterator(&CIs)) != NULL) {
            CI->offsetAEnd.mean      += -currentOffset.mean;
            CI->offsetAEnd.variance  += variance;
            CI->offsetBEnd.mean      += -currentOffset.mean;
            CI->offsetBEnd.variance  += variance;
          }

          currentOffset.mean     = 0.0;
          currentOffset.variance = 0.0;
        } else {
          //  newScaffold:     --------------->
          //  thisScaffold:              ------------------
          currentOffset.variance = GetVarianceOffset(newScaffold, currentOffset.mean, TRUE);
        }
      }
      assert(currentOffset.variance >= 0);

      oldScaffolds.push_back(thisScaffold->id);

      InsertScaffoldContentsIntoScaffold(ScaffoldGraph,
                                         newScaffoldID, thisScaffold->id,
                                         orientCI, &currentOffset,
                                         currentSetID,
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

      currentOffset.mean     += thisScaffold->bpLength.mean;
      currentOffset.variance += thisScaffold->bpLength.variance;

      assert(thisScaffold->bpLength.variance >= 0);
      assert(currentOffset.mean >= 0);
      assert(currentOffset.variance >= 0);

      numMerged++;

      //  Assume no more scaffolds to merge in...
      edge       = NULL;
      neighbor   = NULL;
      neighborID = NULLINDEX;

      //  ...unless there are.
      if ((orientCI.isForward() == true) && (thisScaffold->essentialEdgeB != NULLINDEX)) {
        edge       = GetGraphEdge(ScaffoldGraph->ScaffoldGraph, thisScaffold->essentialEdgeB);
        neighbor   = GetGraphNode(ScaffoldGraph->ScaffoldGraph, (edge->idA == thisScaffold->id) ? edge->idB : edge->idA);
        neighborID = neighbor->id;
      }

      if ((orientCI.isForward() == false) && (thisScaffold->essentialEdgeA != NULLINDEX)) {
        edge       = GetGraphEdge(ScaffoldGraph->ScaffoldGraph, thisScaffold->essentialEdgeA);
        neighbor   = GetGraphNode(ScaffoldGraph->ScaffoldGraph, (edge->idA == thisScaffold->id) ? edge->idB : edge->idA);
        neighborID = neighbor->id;
      }
    }  //  while (neighbor != (CIScaffoldT *)NULL)

    // New scaffold fully popuated now.

    if (iSpec->contigNow == TRUE) {
      CIScaffoldT *scaffold           = GetGraphNode(ScaffoldGraph->ScaffoldGraph, newScaffoldID);
      int32        numScaffoldsBefore = GetNumCIScaffoldTs(ScaffoldGraph->CIScaffolds);

      if (LeastSquaresGapEstimates(ScaffoldGraph, scaffold, LeastSquares_Cleanup) == false) {
        //fprintf(stderr, "ADDING EDGE TO IGNORE: %d <-> %d orient %c distance %.1f +- %.1f\n",
        //        edgeLabel.idA,
        //        edgeLabel.idB,
        //        edgeLabel.orient.toLetter(),
        //        edgeLabel.distance.mean, sqrt(edgeLabel.distance.variance));
      }  //  Bad LeastSquaresEstimate

      ScaffoldSanity(ScaffoldGraph, scaffold);
    }  //  if (iSpec->contigNow == TRUE && .....)

#warning numLiveScaffolds is probably broken, and is nearly unused
    ScaffoldGraph->numLiveScaffolds += (1 - numMerged);

    currentSetID++;
  }


  //  For all the old scaffolds, kill the edges.

  RemoveDeadSEdges(oldScaffolds);

  //  For all the new scaffolds, create new edges.

  {
    vector<CDS_CID_t> rawEdges;

    BuildSEdges(newScaffolds, rawEdges, TRUE);
    MergeAllGraphEdges(ScaffoldGraph->ScaffoldGraph, rawEdges, TRUE, FALSE);
  } 


  return mergedSomething;
}
