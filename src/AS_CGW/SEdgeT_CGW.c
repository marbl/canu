
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
static char *rcsid = "$Id: SEdgeT_CGW.c,v 1.15 2008-10-08 22:02:55 brianwalenz Exp $";

//#define DEBUG 1
//#define TRY_IANS_SEDGES

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>

#include "AS_global.h"
#include "AS_UTL_Var.h"
#include "AS_UTL_interval.h"
#include "AS_CGW_dataTypes.h"
#include "Globals_CGW.h"
#include "GraphCGW_T.h"
#include "ScaffoldGraph_CGW.h"
#include "ScaffoldGraphIterator_CGW.h"

#ifdef TRY_IANS_SEDGES
#include "AS_CGW_EdgeDiagnostics.h"
#endif

#define  LOG_LINE_LEN  255

/* ************************************************************/
void PrintSEdgeT(FILE *fp, ScaffoldGraphT *graph, char *label, SEdgeT *edge, CDS_CID_t sid){
  char actualOverlap[LOG_LINE_LEN + 1];
  int actual = 0;
  int delta = 0;
  char* flag = "  ";
  CIScaffoldT *scaffoldA = GetCIScaffoldT(graph->CIScaffolds, edge->idA);
  CIScaffoldT *scaffoldB = GetCIScaffoldT(graph->CIScaffolds, edge->idB);

  assert(scaffoldA && scaffoldB);

  if(edge->flags.bits.isDeleted) return;

  if(edge->flags.bits.isBogus){
    if(edge->flags.bits.isProbablyBogus)
      strcpy(actualOverlap," *Bogus and Prob Bogus*");
    else
      strcpy(actualOverlap," *Bogus*");
  }else if(scaffoldA->aEndCoord > 0 && scaffoldB->aEndCoord > 0){
    actual = -IntervalsOverlap(scaffoldA->aEndCoord, scaffoldA->bEndCoord,
                               scaffoldB->aEndCoord, scaffoldB->bEndCoord,-500000);
    delta = edge->distance.mean - actual;
    sprintf(actualOverlap,"actual = %d(%d) %s", actual,delta, "");
    assert (strlen (actualOverlap) < LOG_LINE_LEN);
  }else
    *actualOverlap = '\0';

  if(edge->flags.bits.hasContributingOverlap){
    if(edge->flags.bits.isPossibleChimera)
      flag = "$C";
    else
      flag = "$O";
  }

  fprintf(fp,"\t  cidA:" F_CID " cidB:" F_CID " weight:%d %s ori:%c con:%d distance:%d stddev:%g %s (" F_CID "," F_CID ")\n",
          edge->idA, edge->idB,
          edge->edgesContributing,
          flag,
          GetEdgeOrientationWRT(edge, sid),
          edge->flags.bits.hasContainmentOverlap,
          (int)edge->distance.mean, sqrt(edge->distance.variance),
          actualOverlap, edge->fragA, edge->fragB);

#ifdef NEVER
  fprintf(fp,
          "\tidA:" F_CID " idB:" F_CID " orient:%s edgesContributing:%d quality:%d all:" F_S64 "\n",
          edge->idA,
          edge->idB,
          edge->orient == AB_AB ? "AB_AB" :
          (edge->orient == AB_BA ? "AB_BA" :
           (edge->orient == BA_AB ? "BA_AB" :
            (edge->orient == BA_BA ? "BA_BA" : "XX_XX"))),
          edge->edgesContributing,
          (int) edge->quality,
          edge->flags.all);
  fprintf(fp, "\t\tisInferred:%d\n", edge->flags.bits.isInferred);
  fprintf(fp, "\t\tisTentative:%d\n", edge->flags.bits.isTentative);
  fprintf(fp, "\t\tisLeastSquares:%d\n", edge->flags.bits.isLeastSquares);
  fprintf(fp, "\t\tisEssential:%d\n", edge->flags.bits.isEssential);
  fprintf(fp, "\t\twasEssential:%d\n", edge->flags.bits.wasEssential);
  fprintf(fp, "\t\tisActive:%d\n", edge->flags.bits.isActive);
  fprintf(fp, "\t\tisConfirmed:%d\n", edge->flags.bits.isConfirmed);
  fprintf(fp, "\t\tisContigConfirming:%d\n", edge->flags.bits.isContigConfirming);
  fprintf(fp, "\t\tisUniquetoUnique:%d\n", edge->flags.bits.isUniquetoUnique);
  fprintf(fp, "\t\tisTransitivelyRemoved:%d\n", edge->flags.bits.isTransitivelyRemoved);
  fprintf(fp, "\t\tisInferredRemoved:%d\n", edge->flags.bits.isInferredRemoved);
  fprintf(fp, "\t\tisRedundantRemoved:%d\n", edge->flags.bits.isRedundantRemoved);
  fprintf(fp, "\t\tisDeleted:%d\n", edge->flags.bits.isDeleted);
  fprintf(fp, "\t\tisPossibleChimera:%d\n", edge->flags.bits.isPossibleChimera);
  fprintf(fp, "\t\tinducedByUnknownOrientation:%d\n", edge->flags.bits.inducedByUnknownOrientation);
  fprintf(fp, "\t\thasContributingOverlap:%d\n", edge->flags.bits.hasContributingOverlap);
  fprintf(fp, "\t\taContainsB:%d\n", edge->flags.bits.aContainsB);
  fprintf(fp, "\t\tbContainsA:%d\n", edge->flags.bits.bContainsA);
  fprintf(fp, "\t\tmustOverlap:%d\n", edge->flags.bits.mustOverlap);
  fprintf(fp, "\t\thasTransChunk:%d\n", edge->flags.bits.hasTransChunk);
  fprintf(fp, "\t\thasContainmentOverlap:%d\n", edge->flags.bits.hasContainmentOverlap);
  fprintf(fp, "\t\tisRaw:%d\n", edge->flags.bits.isRaw);
  fprintf(fp, "\t\thasExtremalAFrag:%d\n", edge->flags.bits.hasExtremalAFrag);
  fprintf(fp, "\t\thasExtremalBFrag:%d\n", edge->flags.bits.hasExtremalBFrag);
  fprintf(fp, "\t\trangeTruncated:%d\n", edge->flags.bits.rangeTruncated);
  fprintf(fp, "\t\tinAssembly:%d\n", edge->flags.bits.inAssembly);
  fprintf(fp, "\t\tisBogus:%d\n", edge->flags.bits.isBogus);
  fprintf(fp, "\t\tisProbablyBogus:%d\n", edge->flags.bits.isProbablyBogus);
  fprintf(fp, "\t\thasConfirmingPath:%d\n", edge->flags.bits.hasConfirmingPath);
  fprintf(fp, "\t\tedgeStatus:%d\n", edge->flags.bits.edgeStatus);
  fprintf(fp, "\t\tisMarkedForDeletion:%d\n", edge->flags.bits.isMarkedForDeletion);
  fprintf(fp, "\t\tMeanChangedByWalking:%d\n", edge->flags.bits.MeanChangedByWalking);
  fprintf(fp, "\t\thighQualityA:%d\n", edge->flags.bits.highQualityA);
  fprintf(fp, "\t\thighQualityB:%d\n", edge->flags.bits.highQualityB);
  fprintf(fp, "\t\tisSloppy:%d\n", edge->flags.bits.isSloppy);
  fprintf(fp, "\t\tisBridge:%d\n", edge->flags.bits.isBridge);
  fprintf(fp, "\tdistance.mean:%.2f\n", edge->distance.mean);
  fprintf(fp, "\tdistance.variance:%.2f\n", edge->distance.variance);
  fprintf(fp, "\tnextALink:" F_CID "\n", edge->nextALink);
  fprintf(fp, "\tnextBLink:" F_CID "\n", edge->nextBLink);
  fprintf(fp, "\tprevALink:" F_CID "\n", edge->prevALink);
  fprintf(fp, "\tprevBLink:" F_CID "\n", edge->prevBLink);
  fprintf(fp, "\tminDistance:%.2f\n", edge->minDistance);
  fprintf(fp, "\tfragA:" F_CID "\n", edge->fragA);
  fprintf(fp, "\tfragB:" F_CID "\n", edge->fragB);
  fprintf(fp, "\tdistIndex:" F_CID "\n", edge->distIndex);
  fprintf(fp, "\tnextRawEdge:" F_CID "\n", edge->nextRawEdge);
  fprintf(fp, "\ttopLevelEdge:" F_CID "\n", edge->topLevelEdge);
  fprintf(fp, "\treferenceEdge:" F_CID "\n",edge->referenceEdge);
#endif
}

/* CorrectEdgeVariance calculates and returns
   the amount of variance to add when flipping an edge by subtracting the
   old offset variance and adding the new offset variance. */

double CorrectEdgeVariance(ScaffoldGraphT *graph, CIEdgeT *edge){
  LengthT oldOffsetA, newOffsetA, oldOffsetB, newOffsetB;
  FragOrient fragOrient;
  CDS_CID_t extremalFrag;
  NodeCGW_T *nodeA = GetGraphNode(graph->RezGraph, edge->idA);
  NodeCGW_T *nodeB = GetGraphNode(graph->RezGraph, edge->idB);
  CIFragT *fragA = GetCIFragT(graph->CIFrags, edge->fragA);
  CIFragT *fragB = GetCIFragT(graph->CIFrags, edge->fragB);

  if(FragOffsetAndOrientation(fragA, nodeA, &oldOffsetA, &fragOrient,
                              &extremalFrag, TRUE) == FALSE ||
     FragOffsetAndOrientation(fragA, nodeA, &newOffsetA, &fragOrient,
                              &extremalFrag, FALSE) == FALSE ||
     FragOffsetAndOrientation(fragB, nodeB, &oldOffsetB, &fragOrient,
                              &extremalFrag, TRUE) == FALSE ||
     FragOffsetAndOrientation(fragB, nodeB, &newOffsetB, &fragOrient,
                              &extremalFrag, FALSE) == FALSE)
    assert(0);
  return((newOffsetA.variance -  oldOffsetA.variance) +
	 (newOffsetB.variance -  oldOffsetB.variance));
}


/* Compute the offset and orientation of a fragment in its chunk and scaffold
   Offset is from 5p end of fragment to the end of the chunk/scaffold end in the
   direction of the 3p end of the fragment.

*/

int CIOffsetAndOrientation(ScaffoldGraphT *graph,
			   CDS_CID_t         cid, // input
			   CIOrient   chunkOrient,    // input orientation of fragment in CI
			   //			   int        *sid,// output:  scaffold in which this fragment resides
			   LengthT    *ciOffset,    // CI's offset within Scaffold
			   LengthT    *ciFlipOffset,    // CI's offset within Scaffold if we flip the edge
			   CIOrient   *ciOrient){    // CI's orientation within Scaffold

  ChunkInstanceT *CI = GetGraphNode(graph->RezGraph, cid);
  CIScaffoldT *CIS;
  CIOrient CIInScaffoldOrient;

  if(CI->scaffoldID == NULLINDEX)
    return FALSE;

  //  *sid = CI->scaffoldID;
  CIS = GetCIScaffoldT(graph->CIScaffolds, CI->scaffoldID);

  CIInScaffoldOrient = GetNodeOrient(CI);

#ifdef DEBUG_SEDGE
  fprintf(stderr,"* CI " F_CID " is oriented %c in scaffold " F_CID " (FragOrient = %c)\n",
	  cid, CIInScaffoldOrient, CI->scaffoldID, chunkOrient);
#endif
  /* Find the offset of the B end of the chunk within its chunkOfScaffolds.  The fragments
     offset is additive with this */

  switch(chunkOrient){
    case A_B:  // Fragment is A_B in CI
      switch(CIInScaffoldOrient){
        case A_B:// Chunk is A_B in COS
          /*    COS    ------------------------------------------------------->
                CI      A_B     ---------------------->
                Frag       A_B            ------>
                Offset                                |==========================|
          */

          *ciOrient = A_B;
          ciOffset->mean =
            (CIS->bpLength.mean -  CI->offsetBEnd.mean);
          ciFlipOffset->mean = CI->offsetBEnd.mean;

          ciOffset->variance =
            (CIS->bpLength.variance - CI->offsetBEnd.variance);
          ciFlipOffset->variance = CI->offsetAEnd.variance;

          if(ciOffset->variance < 0.0){
            fprintf(stderr,"* A_B Negative offset variance %g for position of CI " F_CID " in scaffold " F_CID "==> set to 1\n",
                    ciOffset->variance, CI->id, CIS->id);
            ciOffset->variance = 1.0;
          }
          break;
        case B_A: // Chunk is B_A in COS
          /*    COS    ------------------------------------------------------->
                Chunk     B_A  <----------------------
                Frag      A_B         <------
                Offset |=======|
          */
          *ciOrient = B_A;
          ciOffset->mean = CI->offsetBEnd.mean;
          ciFlipOffset->mean =
            (CIS->bpLength.mean - CI->offsetBEnd.mean);
          ciOffset->variance = CI->offsetBEnd.variance;
          ciFlipOffset->variance =
            (CIS->bpLength.variance - CI->offsetAEnd.variance);

          if(ciFlipOffset->variance < 0.0){
            fprintf(stderr,"* A_B Negative Flip offset variance %g for position of CI " F_CID " in scaffold " F_CID "==> set to 1\n",
                    ciFlipOffset->variance, CI->id, CIS->id);
            ciFlipOffset->variance = 1.0;
          }
          break;
        default:
          assert(0);
      }
      break;

    case B_A:  // Fragment is B_A in Chunk
      switch(CIInScaffoldOrient){
        case A_B:// Chunk is A_B in COS
          /*    COS    ------------------------------------------------------->
                Chunk      A_B     ---------------------->
                Frag       B_A            <------
                Offset |===========|
          */

          *ciOrient = B_A; // in scaffold
          ciOffset->mean = CI->offsetAEnd.mean;
          ciFlipOffset->mean =
            (CIS->bpLength.mean - CI->offsetAEnd.mean);
          ciOffset->variance = CI->offsetAEnd.variance;
          ciFlipOffset->variance =
            (CIS->bpLength.variance - CI->offsetBEnd.variance);

          if(ciFlipOffset->variance < 0.0){
            fprintf(stderr,"* B_A Negative Flip offset variance %g for position of CI " F_CID " in scaffold " F_CID "==> set to 1\n",
                    ciFlipOffset->variance, CI->id, CIS->id);
            ciFlipOffset->variance = 1.0;
          }
          break;
        case B_A: // Chunk is B_A in COS
          /*    COS    ------------------------------------------------------->
                Chunk           <----------------------   B_A
                Frag                   ------>    B_A
                Offset                                |========================|
          */
          *ciOrient = A_B; // in scaffold
          ciOffset->mean =
            (CIS->bpLength.mean - CI->offsetAEnd.mean);
          ciFlipOffset->mean = CI->offsetAEnd.mean;

          ciOffset->variance =
            (CIS->bpLength.variance - CI->offsetAEnd.variance);
          ciFlipOffset->variance = CI->offsetBEnd.variance;

          if(ciOffset->variance < 0.0){
            fprintf(stderr,"* B_A Negative offset variance %g for position of CI " F_CID " in scaffold " F_CID "==> set to 1\n",
                    ciOffset->variance, CI->id, CIS->id);
            ciOffset->variance = 1.0;
          }
          break;

        default:
          assert(0);
      }
      break;
    default:
      assert(0);
  }

  return TRUE;

}

void PopulateReverseEdge(EdgeCGW_T * reverseEdge, EdgeCGW_T * forwardEdge)
{
  *reverseEdge = *forwardEdge;

  reverseEdge->idA = forwardEdge->idB;
  reverseEdge->idB = forwardEdge->idA;

  switch(forwardEdge->orient)
    {
      case AB_BA:
      case BA_AB:
        break;
      case AB_AB:
        reverseEdge->orient = BA_BA;
        break;
      case BA_BA:
        reverseEdge->orient = AB_AB;
        break;
      default:
        assert(0);
        break;
    }

  reverseEdge->flags.bits.aContainsB = forwardEdge->flags.bits.bContainsA;
  reverseEdge->flags.bits.bContainsA = forwardEdge->flags.bits.aContainsB;

  reverseEdge->flags.bits.highQualityA = forwardEdge->flags.bits.highQualityB;
  reverseEdge->flags.bits.highQualityB = forwardEdge->flags.bits.highQualityA;

  reverseEdge->nextALink = forwardEdge->nextBLink;
  reverseEdge->nextBLink = forwardEdge->nextALink;
  reverseEdge->prevALink = forwardEdge->prevBLink;
  reverseEdge->prevBLink = forwardEdge->prevALink;
}


/************************************************************************************/
int BuildSEdgeFromChunkEdge(ScaffoldGraphT * graph,
                            ChunkInstanceT * thisCI,
                            ChunkInstanceT * otherCI,
                            CIEdgeT * edge,
                            int canonicalOnly)
{
  CIOrient ciOrient, mciOrient;
  ChunkOrientationType sedgeOrient;
  LengthT ciOffset, ciFlipOffset, mciOffset, mciFlipOffset;
  LengthT distance, flipDistance;
  SEdgeT sedge = {0};
  ChunkOrientationType edgeOrient;
  FragOrient orient;
  int CIok, mCIok;
  CIScaffoldT * scaffold =
    GetCIScaffoldT(graph->CIScaffolds, thisCI->scaffoldID);
  CIScaffoldT * otherScaffold =
    GetCIScaffoldT(graph->CIScaffolds, otherCI->scaffoldID);

  if(canonicalOnly && otherCI->scaffoldID < thisCI->scaffoldID)
    return TRUE;

  InitGraphEdge(&sedge);

#ifndef TRY_IANS_SEDGES
  edgeOrient = GetEdgeOrientationWRT(edge, otherCI->id);
  orient = ((edgeOrient == AB_BA || edgeOrient == AB_AB)?A_B:B_A);

#ifdef DEBUG_SEDGE
  fprintf(stderr,"* Edge %s (" F_CID "," F_CID ") %c dist: %d in scaffolds (" F_CID "," F_CID ") orient = %c\n",
          (edge->flags.bits.isBogus?"*Bogus*":"     "),
          thisCI->id, otherCI->id, edgeOrient, (int)edge->distance.mean,
          thisCI->scaffoldID, otherCI->scaffoldID, orient);
#endif

  mCIok = CIOffsetAndOrientation(graph,
                                 otherCI->id, // input
                                 orient,    // input orientation of fragment in CI
                                 &mciOffset,    // CI's offset within Scaffold
                                 &mciFlipOffset,    // CI's offset within Scaffold if we flip the edge
                                 &mciOrient);
  if(!mCIok)
    return FALSE;

  edgeOrient = GetEdgeOrientationWRT(edge, thisCI->id);
  orient = ((edgeOrient == AB_BA || edgeOrient == AB_AB)?A_B:B_A);
  CIok = CIOffsetAndOrientation(graph,
                                thisCI->id, // input
                                orient,    // input orientation of fragment in CI
                                &ciOffset,    // CI's offset within Scaffold
                                &ciFlipOffset,    // CI's offset within Scaffold if we flip the edge
                                &ciOrient);
  assert(CIok);

#ifdef DEBUG_SEDGE
  fprintf(stderr,"* mciOffset = %d mciOrient = %c  ciOffset = %d ciOrient = %c\n",
          (int)mciOffset.mean, mciOrient, (int)ciOffset.mean, ciOrient);
#endif

  /* Mate pairs must be oriented in opposite directions.
     So, if they are oriented the same wrt
     their own chunk, the chunks must be oriented opposite one another */
  switch(ciOrient){
    //
    case A_B:

      switch(mciOrient){
        case A_B:
          //           length - 5'             gap            length - 5'
          //      |------------------------||---------------||-----------|
          //  A --------------------------- B               B --------------------------- A
          //    5'----->                                           <------5'
          //      |-------------------------------------------------------|
          //                             mate distance
          //
          sedgeOrient = AB_BA;
          break;
        case B_A:
          //           length - 5'             gap                5'
          //      |------------------------||---------------||-----------|
          //  A --------------------------- B               A --------------------------- B
          //    5'----->                                           <------5'
          //      |-------------------------------------------------------|
          //                             mate distance
          //
          sedgeOrient = AB_AB;
          break;
        default:
          assert(0);
          break;
      }
      break;
    case B_A:

      switch(mciOrient){
        case A_B:
          //                     5'             gap            length - 5'
          //      |------------------------||---------------||-----------|
          //  B --------------------------- A               B --------------------------- A
          //    5'----->                                           <------5'
          //      |-------------------------------------------------------|
          //                             mate distance
          //
          sedgeOrient = BA_BA;
          break;
        case B_A:
          //                     5'             gap                5'
          //      |------------------------||---------------||-----------|
          //  B --------------------------- A               A --------------------------- B
          //    5'----->                                           <------5'
          //      |-------------------------------------------------------|
          //                             mate/guide distance
          //
          sedgeOrient = BA_AB;
          break;
        default:
          assert(0);
          break;
      }
      break;
    default:
      assert(0);
      break;
  }
  distance.mean = edge->distance.mean - ciOffset.mean - mciOffset.mean;
  flipDistance.mean = - (edge->distance.mean + ciFlipOffset.mean +
                         mciFlipOffset.mean);
  if(distance.mean < flipDistance.mean){
    distance.mean = flipDistance.mean;
    // Since the two offsets and the dist are independent we SUM their
    // variances but we need to correct for the variance of the offset
    // already included in the edge->distance.variance
    distance.variance = edge->distance.variance + ciFlipOffset.variance +
      mciFlipOffset.variance + CorrectEdgeVariance(graph, edge);
    sedgeOrient = InvertEdgeOrient(sedgeOrient);
  }else{
    // Since the two offsets and the dist are independent we SUM their variances
    distance.variance = edge->distance.variance + ciOffset.variance +
      mciOffset.variance;
  }
  sedge.orient = sedgeOrient;
  sedge.distance = distance;

#else
  /* Ian's attempt to simplify, goes straight to fragments
     get edge between fragments, innie or outtie, distance
  */
  PopulateChunkEdgeBasics(graph,
                          GetCIFragT(graph->CIFrags,
                                     (edge->idA == thisCI->id ?
                                      edge->fragA: edge->fragB)),
                          scaffold,
                          GetCIFragT(graph->CIFrags,
                                     (edge->idB == otherCI->id ?
                                      edge->fragB: edge->fragA)),
                          otherScaffold,
                          GetDistT(graph->Dists, edge->distIndex),
                          &sedge);
  distance = sedge.distance;
#endif

  sedge.idA = thisCI->scaffoldID;
  sedge.idB = otherCI->scaffoldID;
  sedge.fragA = edge->fragA;
  sedge.fragB = edge->fragB;
  sedge.flags = edge->flags;
  sedge.flags.bits.aContainsB = sedge.flags.bits.bContainsA = FALSE;
  sedge.flags.bits.hasContainmentOverlap = FALSE;

  // BE CAREFUL: hasContainmentOverlap is advisory.  aContainsB and bContainsA say
  // this is an overlap edge
  if(distance.mean < 0.0){
    double absmean = -distance.mean;
    sedge.flags.bits.hasContainmentOverlap = TRUE;
    if(absmean > scaffold->bpLength.mean){
      ;// sedge.flags.bits.bContainsA = TRUE;
    }else if( absmean > otherScaffold->bpLength.mean){
      ; //sedge.flags.bits.aContainsB = TRUE;
    }else{
      sedge.flags.bits.hasContainmentOverlap = FALSE;
    }
  }

  sedge.edgesContributing = edge->edgesContributing;
  sedge.topLevelEdge = NULLINDEX;
  sedge.referenceEdge =
    (CDS_CID_t)GetVAIndex_CIEdgeT(graph->RezGraph->edges, edge);
#if 0
  fprintf(stderr,"*SEdge (" F_CID "," F_CID ") induced by edge#" F_CID " (" F_CID "," F_CID ")\n",
          sedge.idA, sedge.idB,
          sedge.referenceEdge,
          edge->idA, edge->idB);
  {
    CIEdgeT *edge = GetGraphEdge(graph->RezGraph, sedge.referenceEdge);
    CIEdgeT *topEdge = GetGraphEdge(graph->RezGraph, edge->topLevelEdge);
    assert(edge->idA == topEdge->idA);
    assert(edge->idB == topEdge->idB);
  }
#endif
  sedge.nextALink = sedge.nextBLink = sedge.prevALink = sedge.prevBLink = NULLINDEX;

  if(sedge.idA > sedge.idB)
    {
      SEdgeT reverseEdge = {0};

      assert(canonicalOnly == FALSE);
      PopulateReverseEdge(&reverseEdge, &sedge);
      sedge = reverseEdge;
    }


  {
    SEdgeT * newEdge = GetFreeGraphEdge(graph->ScaffoldGraph);
    sedge.topLevelEdge =
      GetVAIndex_EdgeCGW_T(graph->ScaffoldGraph->edges, newEdge);
    *newEdge = sedge;
    InsertGraphEdge(graph->ScaffoldGraph, newEdge->topLevelEdge, FALSE);
  }
  return TRUE;
}


void BuildSEdgesForScaffold(ScaffoldGraphT * graph,
                            CIScaffoldT * scaffold,
                            int canonicalOnly,
                            int includeNegativeEdges)
{
  CIScaffoldTIterator CIs;
  GraphEdgeIterator edges;
  CIEdgeT *edge;
  ChunkInstanceT *thisCI;

  if(isDeadCIScaffoldT(scaffold) ||
     scaffold->type != REAL_SCAFFOLD)
    return;

  InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &CIs);

  while((thisCI = NextCIScaffoldTIterator(&CIs)) != NULL){
    InitGraphEdgeIterator(graph->RezGraph,  thisCI->id,   ALL_END,   ALL_EDGES, GRAPH_EDGE_RAW_ONLY, &edges);// ONLY RAW

    while((edge = NextGraphEdgeIterator(&edges)) != NULL){
      int isA = (edge->idA == thisCI->id);
      ChunkInstanceT *otherCI = GetGraphNode(graph->RezGraph,
                                             (isA? edge->idB: edge->idA));

      // RAW EDGES ONLY
      assert(edge->flags.bits.isRaw);

      // Only non-overlap  edges (i,j) s.t. i < j are interesting
      if(isOverlapEdge(edge) ||
         (!includeNegativeEdges && edge->distance.mean < -10000) ||
         isProbablyBogusEdge(edge) ||                // filter junk
         otherCI->scaffoldID == NULLINDEX)             // not scaffolded
        continue;

      if(canonicalOnly)
        {
          if(otherCI->scaffoldID < thisCI->scaffoldID)       // not canonical
            continue;
        }
      else
        {
          if(thisCI->scaffoldID < otherCI->scaffoldID)       // not canonical
            continue;
        }

      if(otherCI->scaffoldID == thisCI->scaffoldID)  // internal
        continue;

      BuildSEdgeFromChunkEdge(graph, thisCI, otherCI, edge, canonicalOnly);
    }
  }
}


void PrintSEdges(ScaffoldGraphT * graph){
  int32 numEdges = (int32) GetNumGraphEdges(graph->ScaffoldGraph);
  CDS_CID_t i;
  for(i = 0; i < numEdges; i++){
    PrintSEdgeT(GlobalData->stderrc, graph, "\t",
                GetGraphEdge(graph->ScaffoldGraph, i), i);
  }
}


void PrintSEdgesForScaffold(ScaffoldGraphT * graph,
                            CIScaffoldT * scaffold)
{
  CIScaffoldTIterator CIs;
  ChunkInstanceT *thisCI;
  int32 numEdges = 0;
  FILE * fp;
  char filename[1024];

  sprintf(filename, "scf%010" F_CIDP "Edges.txt", scaffold->id);
  fp = fopen(filename, "w");
  assert(fp != NULL);

  fprintf(fp, "********************************\n");
  fprintf(fp, "Printing edges for scaffold " F_CID "\n",
          scaffold->id);
  fprintf(fp, "\nPrinting inter-scaffold contig edges for scaffold " F_CID "\n",
          scaffold->id);

  InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &CIs);
  while((thisCI = NextCIScaffoldTIterator(&CIs)) != NULL)
    {
      GraphEdgeIterator edges;
      CIEdgeT *edge;

      InitGraphEdgeIterator(graph->RezGraph,
                            thisCI->id,
                            ALL_END,
                            ALL_EDGES,
                            GRAPH_EDGE_DEFAULT,
                            &edges);
      while((edge = NextGraphEdgeIterator(&edges)) != NULL)
        {
          if(edge->idA != edge->idB)
            {
              int isA = (edge->idA == thisCI->id);
              ChunkInstanceT *otherCI = GetGraphNode(graph->RezGraph,
                                                     (isA? edge->idB: edge->idA));

              fprintf(fp, "" F_CID " to " F_CID " in scaffold " F_CID ". %s",
                      thisCI->id, otherCI->id, otherCI->scaffoldID,
                      (edge->flags.bits.isRaw) ? "Raw\n" : "\n");
              numEdges++;
            }
        }
    }
  fprintf(fp, "Done printing %d inter-scaffold contig edges for " F_CID "\n",
          numEdges, scaffold->id);


  {
    SEdgeTIterator SEdges = {0};
    SEdgeT * sEdge;

    fprintf(fp,
            "\nPrinting raw inter-scaffold scaffold edges for scaffold " F_CID "\n",
            scaffold->id);

    numEdges = 0;
    InitSEdgeTIterator(ScaffoldGraph, scaffold->id,
                       TRUE, FALSE, ALL_END, FALSE, &SEdges);
    while((sEdge = NextSEdgeTIterator(&SEdges)) != NULL)
      {
        fprintf(fp, "" F_CID " to " F_CID ", weight %d\n",
                sEdge->idA, sEdge->idB, sEdge->edgesContributing);
        numEdges++;
      }

    fprintf(fp,
            "\nPrinting merged inter-scaffold scaffold edges for scaffold " F_CID "\n",
            scaffold->id);

    numEdges = 0;
    InitSEdgeTIterator(ScaffoldGraph, scaffold->id,
                       FALSE, FALSE, ALL_END, FALSE, &SEdges);
    while((sEdge = NextSEdgeTIterator(&SEdges)) != NULL)
      {
        fprintf(fp, "" F_CID " to " F_CID ", weight %d\n",
                sEdge->idA, sEdge->idB, sEdge->edgesContributing);
        numEdges++;
      }

    fprintf(fp, "Done printing %d inter-scaffold scaffold edges for " F_CID "\n",
            numEdges, scaffold->id);
  }
  fprintf(fp, "********************************\n");
  fclose(fp);
}


void BuildSEdges(ScaffoldGraphT *graph, int includeNegativeEdges)
{
  CDS_CID_t sid;

  /* Recycle the SEdge VA */
  ResetEdgeCGW_T(graph->ScaffoldGraph->edges);
  graph->ScaffoldGraph->tobeFreeEdgeHead = NULLINDEX;
  graph->ScaffoldGraph->freeEdgeHead = NULLINDEX;

  /* Reset scaffold edge heads */
  for(sid = 0; sid < GetNumGraphNodes(graph->ScaffoldGraph); sid++){
    CIScaffoldT *scaffold = GetGraphNode(graph->ScaffoldGraph,sid);
    scaffold->edgeHead = NULLINDEX;
    scaffold->essentialEdgeA = scaffold->essentialEdgeB = NULLINDEX;
    scaffold->numEssentialA = scaffold->numEssentialB = 0;
    scaffold->flags.bits.smoothSeenAlready = FALSE;
    scaffold->flags.bits.walkedAlready = FALSE;
  }

  /* Build the scaffold edges */
  for(sid = 0; sid < GetNumGraphNodes(graph->ScaffoldGraph); sid++)
    BuildSEdgesForScaffold(graph,
                           GetGraphNode(graph->ScaffoldGraph, sid),
                           TRUE, includeNegativeEdges);

  fprintf(GlobalData->stderrc,"* BuildSEdges: %d edges on completion\n",
          (int) GetNumGraphEdges(graph->ScaffoldGraph));
}
