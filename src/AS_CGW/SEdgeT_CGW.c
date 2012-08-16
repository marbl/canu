
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
static char *rcsid = "$Id: SEdgeT_CGW.c,v 1.25 2012-08-16 03:39:43 brianwalenz Exp $";

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
			   SequenceOrient   chunkOrient,    // input orientation of fragment in CI
			   //			   int        *sid,// output:  scaffold in which this fragment resides
			   LengthT    *ciOffset,    // CI's offset within Scaffold
			   LengthT    *ciFlipOffset,    // CI's offset within Scaffold if we flip the edge
			   SequenceOrient   *ciOrient){    // CI's orientation within Scaffold

  ChunkInstanceT *CI = GetGraphNode(graph->ContigGraph, cid);
  CIScaffoldT *CIS;
  SequenceOrient CIInScaffoldOrient;

  if(CI->scaffoldID == NULLINDEX)
    return FALSE;

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

  return TRUE;
}



/************************************************************************************/
int BuildSEdgeFromChunkEdge(ScaffoldGraphT * graph,
                            ChunkInstanceT * thisCI,
                            ChunkInstanceT * otherCI,
                            CIEdgeT * edge,
                            int canonicalOnly)
{
  SequenceOrient ciOrient, mciOrient;
  PairOrient sedgeOrient;
  LengthT ciOffset, ciFlipOffset, mciOffset, mciFlipOffset;
  LengthT distance, flipDistance;
  SEdgeT sedge;
  PairOrient edgeOrient;
  SequenceOrient orient;
  int CIok, mCIok;
  CIScaffoldT * scaffold      = GetCIScaffoldT(graph->CIScaffolds, thisCI->scaffoldID);
  CIScaffoldT * otherScaffold = GetCIScaffoldT(graph->CIScaffolds, otherCI->scaffoldID);

  if(canonicalOnly && otherCI->scaffoldID < thisCI->scaffoldID)
    return TRUE;

  InitGraphEdge(&sedge);

#ifndef TRY_IANS_SEDGES
  edgeOrient = GetEdgeOrientationWRT(edge, otherCI->id);
  orient.setIsForward(edgeOrient.isAB_BA() || edgeOrient.isAB_AB());

#ifdef DEBUG_SEDGE
  fprintf(stderr,"* Edge %s ("F_CID","F_CID") %c dist: %d in scaffolds ("F_CID","F_CID") orient = %c\n",
          (edge->flags.bits.isBogus?"*Bogus*":"     "),
          thisCI->id, otherCI->id, edgeOrient.toLetter(), (int)edge->distance.mean,
          thisCI->scaffoldID, otherCI->scaffoldID, orient.toLetter());
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
  orient.setIsForward(edgeOrient.isAB_BA() || edgeOrient.isAB_AB());
  CIok = CIOffsetAndOrientation(graph,
                                thisCI->id, // input
                                orient,    // input orientation of fragment in CI
                                &ciOffset,    // CI's offset within Scaffold
                                &ciFlipOffset,    // CI's offset within Scaffold if we flip the edge
                                &ciOrient);
  assert(CIok);

#ifdef DEBUG_SEDGE
  fprintf(stderr,"* mciOffset = %d mciOrient = %c  ciOffset = %d ciOrient = %c\n",
          (int)mciOffset.mean, mciOrient.toLetter(), (int)ciOffset.mean, ciOrient.toLetter());
#endif

  /* Mate pairs must be oriented in opposite directions.
     So, if they are oriented the same wrt
     their own chunk, the chunks must be oriented opposite one another */

  assert(ciOrient.isUnknown() == false);
  assert(mciOrient.isUnknown() == false);

  if (ciOrient.isForward()) {
    if (mciOrient.isForward()) {
      //           length - 5'             gap            length - 5'
      //      |------------------------||---------------||-----------|
      //  A --------------------------- B               B --------------------------- A
      //    5'----->                                           <------5'
      //      |-------------------------------------------------------|
      //                             mate distance
      //
      sedgeOrient.setIsAB_BA();
    } else {
      //           length - 5'             gap                5'
      //      |------------------------||---------------||-----------|
      //  A --------------------------- B               A --------------------------- B
      //    5'----->                                           <------5'
      //      |-------------------------------------------------------|
      //                             mate distance
      //
      sedgeOrient.setIsAB_AB();
    }
  } else {
    if (mciOrient.isForward()) {
      //                     5'             gap            length - 5'
      //      |------------------------||---------------||-----------|
      //  B --------------------------- A               B --------------------------- A
      //    5'----->                                           <------5'
      //      |-------------------------------------------------------|
      //                             mate distance
      //
      sedgeOrient.setIsBA_BA();
    } else {
      //                     5'             gap                5'
      //      |------------------------||---------------||-----------|
      //  B --------------------------- A               A --------------------------- B
      //    5'----->                                           <------5'
      //      |-------------------------------------------------------|
      //                             mate/guide distance
      //
      sedgeOrient.setIsBA_AB();
    }
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
    sedgeOrient.invert();
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
  sedge.flags.bits.aContainsB = FALSE;
  sedge.flags.bits.bContainsA = FALSE;
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
  sedge.referenceEdge = (CDS_CID_t)GetVAIndex_CIEdgeT(graph->ContigGraph->edges, edge);

#if 0
  fprintf(stderr,"*SEdge ("F_CID","F_CID") induced by edge#"F_CID" ("F_CID","F_CID")\n",
          sedge.idA, sedge.idB,
          sedge.referenceEdge,
          edge->idA, edge->idB);
  {
    CIEdgeT *edge = GetGraphEdge(graph->ContigGraph, sedge.referenceEdge);
    CIEdgeT *topEdge = GetGraphEdge(graph->ContigGraph, edge->topLevelEdge);
    assert(edge->idA == topEdge->idA);
    assert(edge->idB == topEdge->idB);
  }
#endif

  sedge.nextALink = NULLINDEX;
  sedge.nextBLink = NULLINDEX;
  sedge.prevALink = NULLINDEX;
  sedge.prevBLink = NULLINDEX;

  if(sedge.idA > sedge.idB) {
    CDS_CID_t  idA = sedge.idA;
    CDS_CID_t  idB = sedge.idB;
    int32      ab  = sedge.flags.bits.aContainsB;
    int32      ba  = sedge.flags.bits.bContainsA;
    int32      hqa = sedge.flags.bits.highQualityA;
    int32      hqb = sedge.flags.bits.highQualityB;

    assert(canonicalOnly == FALSE);

    sedge.idA = idB;
    sedge.idB = idA;

    sedge.orient.flip();

    sedge.flags.bits.aContainsB = ba;
    sedge.flags.bits.bContainsA = ab;

    sedge.flags.bits.highQualityA = hqa;
    sedge.flags.bits.highQualityB = hqb;
  }


  {
    SEdgeT * newEdge = GetFreeGraphEdge(graph->ScaffoldGraph);
    sedge.topLevelEdge = GetVAIndex_EdgeCGW_T(graph->ScaffoldGraph->edges, newEdge);
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
  ChunkInstanceT *thisCI;

  if(isDeadCIScaffoldT(scaffold) ||
     scaffold->type != REAL_SCAFFOLD)
    return;

  InitCIScaffoldTIterator(graph, scaffold, TRUE, FALSE, &CIs);

  while((thisCI = NextCIScaffoldTIterator(&CIs)) != NULL){
    GraphEdgeIterator edges(graph->ContigGraph, thisCI->id, ALL_END, ALL_EDGES);
    CIEdgeT          *edge;

    while((edge = edges.nextRaw()) != NULL){
      int isA = (edge->idA == thisCI->id);
      ChunkInstanceT *otherCI = GetGraphNode(graph->ContigGraph,
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



void BuildSEdges(ScaffoldGraphT *graph, int canonicalOnly, int includeNegativeEdges)
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
                           canonicalOnly,
                           includeNegativeEdges);
}
