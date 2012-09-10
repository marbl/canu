
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
static char *rcsid = "$Id: SEdgeT_CGW.c,v 1.30 2012-09-10 10:55:44 brianwalenz Exp $";

//#define DEBUG 1

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
  }else
    *actualOverlap = '\0';

  if(edge->flags.bits.hasContributingOverlap){
    if(edge->flags.bits.isPossibleChimera)
      flag = "$C";
    else
      flag = "$O";
  }

  fprintf(fp,"\t  cidA:"F_CID" cidB:"F_CID" weight:%d %s ori:%c con:%d distance:%d stddev:%g %s ("F_CID","F_CID")\n",
          edge->idA, edge->idB,
          edge->edgesContributing,
          flag,
          GetEdgeOrientationWRT(edge, sid).toLetter(),
          edge->flags.bits.hasContainmentOverlap,
          (int)edge->distance.mean, sqrt(edge->distance.variance),
          actualOverlap, edge->fragA, edge->fragB);
}

/* CorrectEdgeVariance calculates and returns
   the amount of variance to add when flipping an edge by subtracting the
   old offset variance and adding the new offset variance. */

double CorrectEdgeVariance(ScaffoldGraphT *graph, CIEdgeT *edge){
  LengthT oldOffsetA, newOffsetA, oldOffsetB, newOffsetB;
  SequenceOrient fragOrient;
  CDS_CID_t extremalFrag;
  NodeCGW_T *nodeA = GetGraphNode(graph->ContigGraph, edge->idA);
  NodeCGW_T *nodeB = GetGraphNode(graph->ContigGraph, edge->idB);
  CIFragT *fragA = GetCIFragT(graph->CIFrags, edge->fragA);
  CIFragT *fragB = GetCIFragT(graph->CIFrags, edge->fragB);

  FragOffsetAndOrientation(fragA, nodeA, &oldOffsetA, &fragOrient, &extremalFrag, TRUE);
  FragOffsetAndOrientation(fragA, nodeA, &newOffsetA, &fragOrient, &extremalFrag, FALSE);
  FragOffsetAndOrientation(fragB, nodeB, &oldOffsetB, &fragOrient, &extremalFrag, TRUE);
  FragOffsetAndOrientation(fragB, nodeB, &newOffsetB, &fragOrient, &extremalFrag, FALSE);

  return((newOffsetA.variance -  oldOffsetA.variance) + (newOffsetB.variance -  oldOffsetB.variance));
}


/* Compute the offset and orientation of a fragment in its chunk and scaffold
   Offset is from 5p end of fragment to the end of the chunk/scaffold end in the
   direction of the 3p end of the fragment.

*/

void
CIOffsetAndOrientation(ScaffoldGraphT *graph,
                       CDS_CID_t         cid, // input
                       SequenceOrient   chunkOrient,    // input orientation of fragment in CI
                       LengthT    *ciOffset,    // CI's offset within Scaffold
                       LengthT    *ciFlipOffset,    // CI's offset within Scaffold if we flip the edge
                       SequenceOrient   *ciOrient){    // CI's orientation within Scaffold
  
  ChunkInstanceT *CI = GetGraphNode(graph->ContigGraph, cid);
  CIScaffoldT *CIS;
  SequenceOrient CIInScaffoldOrient;

  assert(CI->scaffoldID != NULLINDEX);

  //  *sid = CI->scaffoldID;
  CIS = GetCIScaffoldT(graph->CIScaffolds, CI->scaffoldID);

  CIInScaffoldOrient = GetNodeOrient(CI);

#ifdef DEBUG_SEDGE
  fprintf(stderr,"* CI "F_CID" is oriented %c in scaffold "F_CID" (FragOrient = %c)\n",
	  cid, CIInScaffoldOrient.toLetter(), CI->scaffoldID, chunkOrient.toLetter());
#endif
  /* Find the offset of the B end of the chunk within its chunkOfScaffolds.  The fragments
     offset is additive with this */

  assert(chunkOrient.isUnknown() == false);
  assert(CIInScaffoldOrient.isUnknown() == false);

  if (chunkOrient.isForward()) {
    if (CIInScaffoldOrient.isForward()) {
      /*    COS    ------------------------------------------------------->
            CI      A_B     ---------------------->
            Frag       A_B            ------>
            Offset                                |==========================|
      */

      ciOrient->setIsForward();
      ciOffset->mean = (CIS->bpLength.mean -  CI->offsetBEnd.mean);
      ciFlipOffset->mean = CI->offsetBEnd.mean;
      
      ciOffset->variance = (CIS->bpLength.variance - CI->offsetBEnd.variance);
      ciFlipOffset->variance = CI->offsetAEnd.variance;

      if(ciOffset->variance < 0.0){
        fprintf(stderr,"* A_B Negative offset variance %g for position of CI "F_CID" in scaffold "F_CID"==> set to 1\n",
                ciOffset->variance, CI->id, CIS->id);
        ciOffset->variance = 1.0;
      }
    } else {
      /*    COS    ------------------------------------------------------->
            Chunk     B_A  <----------------------
            Frag      A_B         <------
            Offset |=======|
      */
      ciOrient->setIsReverse();
      ciOffset->mean = CI->offsetBEnd.mean;
      ciFlipOffset->mean = (CIS->bpLength.mean - CI->offsetBEnd.mean);
      ciOffset->variance = CI->offsetBEnd.variance;
      ciFlipOffset->variance = (CIS->bpLength.variance - CI->offsetAEnd.variance);

      if(ciFlipOffset->variance < 0.0){
        fprintf(stderr,"* A_B Negative Flip offset variance %g for position of CI "F_CID" in scaffold "F_CID"==> set to 1\n",
                ciFlipOffset->variance, CI->id, CIS->id);
        ciFlipOffset->variance = 1.0;
      }
    }

  } else {
    if (CIInScaffoldOrient.isForward()) {
      /*    COS    ------------------------------------------------------->
            Chunk      A_B     ---------------------->
            Frag       B_A            <------
            Offset |===========|
      */
      
      ciOrient->setIsReverse(); // in scaffold
      ciOffset->mean = CI->offsetAEnd.mean;
      ciFlipOffset->mean = (CIS->bpLength.mean - CI->offsetAEnd.mean);
      ciOffset->variance = CI->offsetAEnd.variance;
      ciFlipOffset->variance = (CIS->bpLength.variance - CI->offsetBEnd.variance);

      if(ciFlipOffset->variance < 0.0){
        fprintf(stderr,"* B_A Negative Flip offset variance %g for position of CI "F_CID" in scaffold "F_CID"==> set to 1\n",
                ciFlipOffset->variance, CI->id, CIS->id);
        ciFlipOffset->variance = 1.0;
      }
    } else {
      /*    COS    ------------------------------------------------------->
            Chunk           <----------------------   B_A
            Frag                   ------>    B_A
            Offset                                |========================|
      */
      ciOrient->setIsForward(); // in scaffold
      ciOffset->mean = (CIS->bpLength.mean - CI->offsetAEnd.mean);
      ciFlipOffset->mean = CI->offsetAEnd.mean;

      ciOffset->variance = (CIS->bpLength.variance - CI->offsetAEnd.variance);
      ciFlipOffset->variance = CI->offsetBEnd.variance;

      if(ciOffset->variance < 0.0){
        fprintf(stderr,"* B_A Negative offset variance %g for position of CI "F_CID" in scaffold "F_CID"==> set to 1\n",
                ciOffset->variance, CI->id, CIS->id);
        ciOffset->variance = 1.0;
      }
    }
  }

  assert(ciOrient->isUnknown() == false);
}



/************************************************************************************/
static
CDS_CID_t
BuildSEdgeFromChunkEdge(ChunkInstanceT  *thisCI,
                        ChunkInstanceT  *thatCI,
                        CIEdgeT         *edge) {
  SequenceOrient ciOrient;
  SequenceOrient miOrient;

  LengthT ciOffset;
  LengthT ciFlipOffset;

  LengthT miOffset;
  LengthT miFlipOffset;

  SEdgeT sedge;

  CIScaffoldT  *thisScf  = GetCIScaffoldT(ScaffoldGraph->CIScaffolds, thisCI->scaffoldID);
  CIScaffoldT  *thatScf  = GetCIScaffoldT(ScaffoldGraph->CIScaffolds, thatCI->scaffoldID);

  assert(thisCI->scaffoldID < thatCI->scaffoldID);

  InitGraphEdge(&sedge);

  {
    PairOrient      edgeOrient = GetEdgeOrientationWRT(edge, thatCI->id);
    SequenceOrient  orient;

    orient.setIsForward(edgeOrient.isAB_BA() || edgeOrient.isAB_AB());

    CIOffsetAndOrientation(ScaffoldGraph, thatCI->id, orient, &miOffset, &miFlipOffset, &miOrient);
  }

  {
    PairOrient      edgeOrient = GetEdgeOrientationWRT(edge, thisCI->id);
    SequenceOrient  orient;

    orient.setIsForward(edgeOrient.isAB_BA() || edgeOrient.isAB_AB());

    CIOffsetAndOrientation(ScaffoldGraph, thisCI->id, orient, &ciOffset, &ciFlipOffset, &ciOrient);
  }

  assert(ciOrient.isUnknown() == false);
  assert(miOrient.isUnknown() == false);

#warning ASSUMES INNIE ORIENT

  PairOrient sedgeOrient;

  if ((ciOrient.isForward() == true)  && (miOrient.isForward()) == true)    sedgeOrient.setIsAB_BA();
  if ((ciOrient.isForward() == true)  && (miOrient.isForward()) == false)   sedgeOrient.setIsAB_AB();
  if ((ciOrient.isForward() == false) && (miOrient.isForward()) == true)    sedgeOrient.setIsBA_BA();
  if ((ciOrient.isForward() == false) && (miOrient.isForward()) == false)   sedgeOrient.setIsBA_AB();

  LengthT distance;
  LengthT flipDistance;

  distance.mean     =    edge->distance.mean - ciOffset.mean - miOffset.mean;
  flipDistance.mean = - (edge->distance.mean + ciFlipOffset.mean + miFlipOffset.mean);

  if (distance.mean < flipDistance.mean) {
    distance.mean = flipDistance.mean;

    // Since the two offsets and the dist are independent we SUM their
    // variances but we need to correct for the variance of the offset
    // already included in the edge->distance.variance

    distance.variance = edge->distance.variance + ciFlipOffset.variance + miFlipOffset.variance + CorrectEdgeVariance(ScaffoldGraph, edge);
    sedgeOrient.invert();

  } else {
    // Since the two offsets and the dist are independent we SUM their variances
    distance.variance = edge->distance.variance + ciOffset.variance + miOffset.variance;
  }

  sedge.orient = sedgeOrient;
  sedge.distance = distance;

  sedge.idA                              = thisCI->scaffoldID;
  sedge.idB                              = thatCI->scaffoldID;
  sedge.fragA                            = edge->fragA;
  sedge.fragB                            = edge->fragB;
  sedge.flags                            = edge->flags;
  sedge.flags.bits.aContainsB            = FALSE;
  sedge.flags.bits.bContainsA            = FALSE;
  sedge.flags.bits.hasContainmentOverlap = FALSE;

  // BE CAREFUL: hasContainmentOverlap is advisory.  aContainsB and bContainsA say
  // this is an overlap edge

  if(distance.mean < 0.0){
    double absmean = -distance.mean;

    sedge.flags.bits.hasContainmentOverlap = TRUE;

    if (absmean > thisScf->bpLength.mean)
      ;// sedge.flags.bits.bContainsA = TRUE;

    else if (absmean > thatScf->bpLength.mean)
      ; //sedge.flags.bits.aContainsB = TRUE;

    else
      sedge.flags.bits.hasContainmentOverlap = FALSE;
  }

  sedge.edgesContributing = edge->edgesContributing;
  sedge.topLevelEdge = NULLINDEX;
  sedge.referenceEdge = (CDS_CID_t)GetVAIndex_CIEdgeT(ScaffoldGraph->ContigGraph->edges, edge);

  sedge.nextRawEdge = NULLINDEX;

  SEdgeT * newEdge = GetFreeGraphEdge(ScaffoldGraph->ScaffoldGraph);

  sedge.topLevelEdge = GetVAIndex_EdgeCGW_T(ScaffoldGraph->ScaffoldGraph->edges, newEdge);

  *newEdge = sedge;

  //  DON'T INSERT!  We return these edges to the caller for insertion into a vector.
  //InsertGraphEdge(ScaffoldGraph->ScaffoldGraph, newEdge->topLevelEdge, FALSE);

  return(sedge.topLevelEdge);
}



void
RemoveDeadSEdges(vector<CDS_CID_t> &oldScaffolds) {

  for (uint32 si=0; si<oldScaffolds.size(); si++) {
    CDS_CID_t    sid      = oldScaffolds[si];
    CIScaffoldT *scaffold = GetGraphNode(ScaffoldGraph->ScaffoldGraph, sid);

    assert(isDeadCIScaffoldT(scaffold) == true);

    //  Move edges from the dead scaffold to the unused edge list

    GraphEdgeIterator edges(ScaffoldGraph->ScaffoldGraph, sid, ALL_END, ALL_EDGES);
    CIEdgeT          *edge;

    while ((edge = edges.nextMerged()) != NULL) {
      UnlinkGraphEdge(ScaffoldGraph->ScaffoldGraph, edge);  //  Unlink it from the graph
      FreeGraphEdge(ScaffoldGraph->ScaffoldGraph, edge);    //  And recycle the components
    }

    ScaffoldGraph->ScaffoldGraph->edgeLists[sid].clear();
  }
}



static
void
BuildSEdgesForScaffold(CIScaffoldT        *scaffold,
                       vector<CDS_CID_t>  &rawEdges,
                       bool                includeNegativeEdges,
                       bool                includeAllEdges) {

  if ((isDeadCIScaffoldT(scaffold) == true) ||
      (scaffold->type != REAL_SCAFFOLD))
    return;

  CIScaffoldTIterator CIs;
  ChunkInstanceT     *CI;

  InitCIScaffoldTIterator(ScaffoldGraph, scaffold, TRUE, FALSE, &CIs);

  while((CI = NextCIScaffoldTIterator(&CIs)) != NULL){
    GraphEdgeIterator edges(ScaffoldGraph->ContigGraph, CI->id, ALL_END, ALL_EDGES);
    CIEdgeT          *edge;

    while ((edge = edges.nextRaw()) != NULL) {
      ChunkInstanceT *MI = GetGraphNode(ScaffoldGraph->ContigGraph, ((edge->idA == CI->id) ? edge->idB : edge->idA));

      assert(edge->flags.bits.isRaw);

      if (isOverlapEdge(edge) == true)
        continue;

      if ((includeNegativeEdges == false) && (edge->distance.mean < -10000))
        continue;

      if (isProbablyBogusEdge(edge) == true)
        continue;

      if (MI->scaffoldID == NULLINDEX)
        continue;

      if (MI->scaffoldID == CI->scaffoldID)
        //  Internal edge
        continue;

      if ((includeAllEdges == false) && (MI->scaffoldID < CI->scaffoldID))
        //  Non-canonical edge.
        continue;

      if ((includeAllEdges == true) &&
          (GetGraphNode(ScaffoldGraph->ScaffoldGraph, CI->scaffoldID)->flags.bits.isHavingEdgesBuilt == true) &&
          (GetGraphNode(ScaffoldGraph->ScaffoldGraph, MI->scaffoldID)->flags.bits.isHavingEdgesBuilt == true) &&
          (MI->scaffoldID < CI->scaffoldID))
        //  Non-canonical edge between two new scaffolds.
        continue;

      if (CI->scaffoldID < MI->scaffoldID)
        rawEdges.push_back(BuildSEdgeFromChunkEdge(CI, MI, edge));
      else
        rawEdges.push_back(BuildSEdgeFromChunkEdge(MI, CI, edge));
    }
  }
}



void
BuildSEdges(vector<CDS_CID_t>  &newScaffolds,
            vector<CDS_CID_t>  &rawEdges,
            bool                includeNegativeEdges) {

  for (uint32 si=0; si<newScaffolds.size(); si++) {
    CDS_CID_t    sid      = newScaffolds[si];
    CIScaffoldT *scaffold = GetGraphNode(ScaffoldGraph->ScaffoldGraph, sid);

    assert(isDeadCIScaffoldT(scaffold) == false);
    assert(scaffold->type              == REAL_SCAFFOLD);

    //  Reset scaffold edge heads

    scaffold->essentialEdgeA                = NULLINDEX;
    scaffold->essentialEdgeB                = NULLINDEX;
    scaffold->numEssentialA                 = 0;
    scaffold->numEssentialB                 = 0;
    scaffold->flags.bits.smoothSeenAlready  = FALSE;
    scaffold->flags.bits.walkedAlready      = FALSE;

    //  Make note of the fact that we're building edges for this scaffold.  If we encounter an edge
    //  between two 'new' scaffolds, we'll only use the canonical version.  Otherwise, all edges are
    //  added.

    scaffold->flags.bits.isHavingEdgesBuilt = true;

    ScaffoldGraph->ScaffoldGraph->edgeLists[sid].clear();
  }


  for (uint32 si=0; si<newScaffolds.size(); si++) {
    CDS_CID_t    sid      = newScaffolds[si];
    CIScaffoldT *scaffold = GetGraphNode(ScaffoldGraph->ScaffoldGraph, sid);

    fprintf(stderr, "BuildSEdges()-- build scaffold edges for scaffold %d\n", sid);

    //  Build SEdges involving the scaffold.  We want to include all edges, regardless of canonical,
    //  since we're only building edges for specific scaffolds.

    BuildSEdgesForScaffold(scaffold, rawEdges, includeNegativeEdges, true);
  }


  for (uint32 si=0; si<newScaffolds.size(); si++) {
    CDS_CID_t    sid      = newScaffolds[si];
    CIScaffoldT *scaffold = GetGraphNode(ScaffoldGraph->ScaffoldGraph, sid);

    scaffold->flags.bits.isHavingEdgesBuilt = false;
  }

}



void
BuildSEdges(vector<CDS_CID_t>  &rawEdges,
            bool                includeNegativeEdges) {

  //  Recycle the SEdge VA

  ResetEdgeCGW_T(ScaffoldGraph->ScaffoldGraph->edges);

  ScaffoldGraph->ScaffoldGraph->edgeLists.clear();

  ResizeEdgeList(ScaffoldGraph->ScaffoldGraph);

  ScaffoldGraph->ScaffoldGraph->tobeFreeEdgeHead = NULLINDEX;
  ScaffoldGraph->ScaffoldGraph->freeEdgeHead     = NULLINDEX;

  //  Reset scaffold edge heads

  for (CDS_CID_t sid=0; sid<GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph); sid++) {
    CIScaffoldT *scaffold = GetGraphNode(ScaffoldGraph->ScaffoldGraph, sid);

    scaffold->essentialEdgeA                = NULLINDEX;
    scaffold->essentialEdgeB                = NULLINDEX;
    scaffold->numEssentialA                 = 0;
    scaffold->numEssentialB                 = 0;
    scaffold->flags.bits.smoothSeenAlready  = FALSE;
    scaffold->flags.bits.walkedAlready      = FALSE;

    ScaffoldGraph->ScaffoldGraph->edgeLists[sid].clear();
  }

  //  Build the scaffold edges.  We're building all edges from scratch, and we
  //  want to filter out non-canonical edges.

  for (CDS_CID_t sid=0; sid<GetNumGraphNodes(ScaffoldGraph->ScaffoldGraph); sid++) {
    CIScaffoldT *scaffold = GetGraphNode(ScaffoldGraph->ScaffoldGraph, sid);

    BuildSEdgesForScaffold(scaffold, rawEdges, includeNegativeEdges, false);
  }
}
