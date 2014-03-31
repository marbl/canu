
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
static char *rcsid = "$Id$";

#undef DEBUG_CONTIG

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>

#include "AS_global.H"
#include "AS_UTL_Var.H"
#include "AS_UTL_interval.H"
#include "AS_UTL_reverseComplement.H"
#include "AS_CGW_dataTypes.H"
#include "Globals_CGW.H"
#include "ScaffoldGraph_CGW.H"
#include "ScaffoldGraphIterator_CGW.H"


void CheckContigs()
{
  GraphNodeIterator CIs;
  NodeCGW_T * contig;

  InitGraphNodeIterator(&CIs, ScaffoldGraph->ContigGraph, GRAPH_NODE_DEFAULT);
  while((contig = NextGraphNodeIterator(&CIs)) != NULL)
    {
      double meanDelta = fabs(contig->offsetAEnd.mean - contig->offsetBEnd.mean);
      double varianceDelta = fabs(contig->offsetAEnd.variance -
                                  contig->offsetBEnd.variance);
      if(meanDelta + .5 < contig->bpLength.mean ||
         meanDelta - .5 > contig->bpLength.mean ||
         varianceDelta + .5 < contig->bpLength.variance ||
         varianceDelta - .5 > contig->bpLength.variance)
        {
          fprintf(stderr,
                  "Contig "F_CID" length (%f,%f) doesn't match offset difference (%f,%f)\n",
                  contig->id, contig->bpLength.mean, contig->bpLength.variance,
                  meanDelta, varianceDelta);
        }
    }

}

void
dumpContigInfo(ChunkInstanceT *contig) {
  int           contigOrientation;
  MultiAlignT  *ma;
  char         *seq1;
  int           len1;

  VA_TYPE(char) *consensus = CreateVA_char(2048);
  VA_TYPE(char) *quality   = CreateVA_char(2048);

  fprintf( stderr, "*********************** contig analysis **************************\n");
  fprintf( stderr, "analyzing contig: %d\n", contig->id);

  if (contig->offsetAEnd.mean < contig->offsetBEnd.mean)
    contigOrientation = 0;
  else
    contigOrientation = 1;

  fprintf(stderr, "contig orientation: %d\t length: %d  contig offsetAEnd: %d\t offsetBEnd: %d\n",
          contigOrientation,
          (int)contig->bpLength.mean,
          (int)contig->offsetAEnd.mean,
          (int)contig->offsetBEnd.mean);

  ma = ScaffoldGraph->tigStore->loadMultiAlign(contig->id, ScaffoldGraph->ContigGraph->type == CI_GRAPH);

  // Get the consensus sequences for the contig from the Store
  GetConsensus(ScaffoldGraph->ContigGraph, contig->id, consensus, quality);

  seq1 = Getchar(consensus, 0);
  len1 = strlen(seq1);

  if (contigOrientation == 1)
    reverseComplementSequence(seq1, len1);

  if (len1 < 5000) {
    fprintf( stderr, ">contig%d consensus seq (flipped to reflect scaff orientation)\n", contig->id);
    fprintf( stderr, "%s\n", seq1);
  } else {
    char tmpchar = seq1[2500];
    seq1[2500] = '\0';

    fprintf( stderr, ">contig%d left end\n", contig->id);
    fprintf( stderr, "%s\n", seq1);

    seq1[2500] = tmpchar;

    fprintf( stderr, ">contig%d right end\n", contig->id);
    fprintf( stderr, "%s\n", seq1 + len1 - 2501);
  }

#if 1
  int numUnitigs = GetNumIntUnitigPoss(ma->u_list);
  fprintf( stderr, "number unitigs: %d\n", numUnitigs);

  int i;
  for (i = 0; i < numUnitigs; i++) {
    IntUnitigPos *upos = GetIntUnitigPos( ma->u_list, i);
    ChunkInstanceT *unitig = GetGraphNode( ScaffoldGraph->CIGraph, upos->ident);
    MultiAlignT *uma = ScaffoldGraph->tigStore->loadMultiAlign(unitig->id, ScaffoldGraph->CIGraph->type == CI_GRAPH);
    IntMultiPos *ump;
    int icntfrag;

    fprintf( stderr, "  unitig: %d\t num frags: %ld surrogate: %d\n", unitig->id, GetNumIntMultiPoss(uma->f_list),
             (unitig->flags.bits.isStoneSurrogate || unitig->flags.bits.isWalkSurrogate));

    if (unitig->flags.bits.isStoneSurrogate ||
        unitig->flags.bits.isWalkSurrogate) {
      fprintf (stderr, "  surrogate unitig offsetAEnd: %f, offsetBEnd: %f\n", unitig->offsetAEnd.mean, unitig->offsetBEnd.mean);

      unitig = GetGraphNode( ScaffoldGraph->CIGraph, unitig->info.CI.baseID);
      fprintf ( stderr, "  using original unitig: %d\n", unitig->id);
      uma = ScaffoldGraph->tigStore->loadMultiAlign(unitig->id,
                                                    ScaffoldGraph->CIGraph->type == CI_GRAPH);
    }

    // now print out info on the frags in the unitig
    for (icntfrag = 0; icntfrag < GetNumIntMultiPoss(uma->f_list); icntfrag++) {
      IntMultiPos *imp = GetIntMultiPos(uma->f_list, icntfrag);
      CIFragT     *frag = GetCIFragT(ScaffoldGraph->CIFrags, imp->ident);

      fprintf(stderr, "    frag: %6d\t contig pos (5p, 3p): %6d, %6d\n",
              imp->ident, (int) frag->contigOffset5p.mean, (int) frag->contigOffset3p.mean);
    }
  }
#endif


#if 1
  CIEdgeT * e;
  GraphEdgeIterator edges(ScaffoldGraph->ContigGraph, contig->id, ALL_END, ALL_EDGES);

  //  FALSE == ITERATOR_VERBOSE

  while((e = edges.nextRaw()) != NULL)
    PrintGraphEdge( stderr, ScaffoldGraph->ContigGraph, "Analyzing edge", e, 0);
#endif

  DeleteVA_char(consensus);
  DeleteVA_char(quality);
}


void
GetContigPositionInScaffold(ChunkInstanceT *contig, int *left_end, int *right_end,
                            int *contigScaffoldOrientation) {
  if (contig->offsetAEnd.mean <= contig->offsetBEnd.mean) {
    *left_end = contig->offsetAEnd.mean;
    *right_end = contig->offsetBEnd.mean;
    *contigScaffoldOrientation = 0;
  } else {
    *left_end = contig->offsetBEnd.mean;
    *right_end = contig->offsetAEnd.mean;
    *contigScaffoldOrientation = 1;
  }
}


static
void
GetFragmentPositionInScaffoldFromContig(CIFragT *frag, int *left_end, int *right_end,
                                        int *fragmentScaffoldOrientation,
                                        int contigLeftEnd, int contigRightEnd, int contigScaffoldOrientation) {
  if (contigScaffoldOrientation == 0) {
    // contig is direct in scaffold
    if (frag->contigOffset5p.mean < frag->contigOffset3p.mean) {
      // frag is direct in contig
      *left_end = contigLeftEnd + frag->contigOffset5p.mean;
      *right_end = contigLeftEnd + frag->contigOffset3p.mean;
      *fragmentScaffoldOrientation = 0;
    } else {
      // frag is reversed in contig
      *left_end = contigLeftEnd + frag->contigOffset3p.mean;
      *right_end = contigLeftEnd + frag->contigOffset5p.mean;
      *fragmentScaffoldOrientation = 1;
    }
  } else {
    // contig is reversed in scaffold
    if (frag->contigOffset5p.mean < frag->contigOffset3p.mean) {
      // frag is direct in contig
      *left_end = contigRightEnd - frag->contigOffset3p.mean;
      *right_end = contigRightEnd - frag->contigOffset5p.mean;
      *fragmentScaffoldOrientation = 1;
    } else {
      // frag is reversed in contig
      *left_end = contigRightEnd - frag->contigOffset5p.mean;
      *right_end = contigRightEnd - frag->contigOffset3p.mean;
      *fragmentScaffoldOrientation = 0;
    }
  }
}


void
GetFragmentPositionInScaffold(CIFragT *frag, int *left_end, int *right_end,
                              int *fragmentScaffoldOrientation) {
  ContigT *containingContig = GetGraphNode(ScaffoldGraph->ContigGraph, frag->contigID);
  int contigLeftEnd, contigRightEnd, contigScaffoldOrientation;

  GetContigPositionInScaffold( containingContig, &contigLeftEnd, &contigRightEnd, &contigScaffoldOrientation);

  GetFragmentPositionInScaffoldFromContig( frag, left_end, right_end, fragmentScaffoldOrientation,
                                           contigLeftEnd, contigRightEnd, contigScaffoldOrientation);
}


void DumpContig(FILE *stream, ScaffoldGraphT *graph, ContigT *contig, int raw){
  int numCI;
  ContigTIterator CIs;
  ChunkInstanceT *CI;

  assert(contig->type == CONTIG_CGW);

  fprintf(stream, "* Contig "F_CID" sc:"F_CID" aoff:%d boff:%d ANext:"F_CID" BNext:"F_CID" len:%d AEnd:"F_CID" BEnd:"F_CID" numCI:%d\nCIs:\n",
          contig->id,
          contig->scaffoldID,
          (int)contig->offsetAEnd.mean,
          (int)contig->offsetBEnd.mean,
          contig->AEndNext, contig->BEndNext,
          (int)contig->bpLength.mean,
          contig->info.Contig.AEndCI,
          contig->info.Contig.BEndCI,
          contig->info.Contig.numCI);

#ifdef DEBUG_CONTIG
  fprintf(stderr, "* Contig "F_CID" ANext:"F_CID" BNext:"F_CID" length:%d AEnd:"F_CID" BEnd:"F_CID" numCI:%d\nCIs:\n",
          contig->id,
          contig->AEndNext, contig->BEndNext, (int)contig->bpLength.mean,
          contig->info.Contig.AEndCI,
          contig->info.Contig.BEndCI,
          contig->info.Contig.numCI);
  fflush(stderr);
#endif
  InitContigTIterator(graph,contig->id, TRUE, FALSE, &CIs);
  numCI = 0;
  while((CI = NextContigTIterator(&CIs)) != NULL){
    int isSurrogate = IsSurrogateNode(CI);
    numCI++;
    fprintf(stream,"\t%c:"F_CID" cov:%d aoff:%d boff:%d anxt:"F_CID" bnxt:"F_CID" ctg:"F_CID" scf:"F_CID"\n",
            (isSurrogate?'s':'u'),
            CI->id,
            ScaffoldGraph->tigStore->getUnitigCoverageStat(CI->id),
            (int)CI->offsetAEnd.mean,
            (int)CI->offsetBEnd.mean,
            CI->AEndNext,
            CI->BEndNext,
            CI->info.CI.contigID,
            CI->scaffoldID);
#ifdef DEBUG_CONTIG
    fprintf(stderr,"\t%c:"F_CID" cov:%d aoff:%d boff:%d anxt:"F_CID" bnxt:"F_CID" ctg:"F_CID" scf:"F_CID"\n",
            (isSurrogate?'s':'u'),
            CI->id,
            ScaffoldGraph->tigStore->getUnitigCoverageStat(CI->id),
            (int)CI->offsetAEnd.mean,
            (int)CI->offsetBEnd.mean,
            CI->AEndNext,
            CI->BEndNext,
            CI->info.CI.contigID,
            CI->scaffoldID);
    fflush(stderr);
#endif
  }

  fprintf(stream,"\t%s Edges from A End:\n", (raw?" Raw ":" Merged "));

  {
    GraphEdgeIterator edges(graph->ContigGraph, contig->id, A_END, ALL_EDGES);
    CIEdgeT *edge;

    while((edge = (raw) ? edges.nextRaw() : edges.nextMerged()) != NULL)
      PrintGraphEdge(stream,graph->ContigGraph,"\t",edge, contig->id);
  }

  fprintf(stream,"\t%s Edges from B End:\n", (raw?" Raw ":" Merged "));

  {
    GraphEdgeIterator edges(graph->ContigGraph, contig->id, B_END, ALL_EDGES);
    CIEdgeT *edge;

    while((edge = (raw) ? edges.nextRaw() : edges.nextMerged()) != NULL)
      PrintGraphEdge(stream,graph->ContigGraph,"\t",edge, contig->id);
  }
}

void DumpContigInScfContext(FILE *stream, ScaffoldGraphT *graph,
                            ContigT *contig, int raw){
  int numCI;
  ContigTIterator CIs;
  ChunkInstanceT *CI;

  assert(contig->type == CONTIG_CGW);

  fprintf(stream, "* Contig "F_CID" sc:"F_CID" aoff:(%d,%e) boff:(%d,%e) ANext:"F_CID" BNext:"F_CID" len:%d AEnd:"F_CID" BEnd:"F_CID" numCI:%d\nCIs:\n",
          contig->id,
          contig->scaffoldID,
          (int32)contig->offsetAEnd.mean,
          contig->offsetAEnd.variance,
          (int32)contig->offsetBEnd.mean,
          contig->offsetBEnd.variance,
          contig->AEndNext, contig->BEndNext,
          (int)contig->bpLength.mean,
          contig->info.Contig.AEndCI,
          contig->info.Contig.BEndCI,
          contig->info.Contig.numCI);

#ifdef DEBUG_CONTIG
  fprintf(stderr, "* Contig "F_CID" ANext:"F_CID" BNext:"F_CID" length:%d AEnd:"F_CID" BEnd:"F_CID" numCI:%d\nCIs:\n",
          contig->id,
          contig->AEndNext, contig->BEndNext, (int)contig->bpLength.mean,
          contig->info.Contig.AEndCI,
          contig->info.Contig.BEndCI,
          contig->info.Contig.numCI);
  fflush(stderr);
#endif
  InitContigTIterator(graph,contig->id, TRUE, FALSE, &CIs);
  numCI = 0;
  while((CI = NextContigTIterator(&CIs))!=NULL){
    int isSurrogate = IsSurrogateNode(CI);
    numCI++;
    fprintf(stream,"\t%c:"F_CID" cov:%d aoff:(%d,%e) boff:(%d,%e) anxt:"F_CID" bnxt:"F_CID" ctg:"F_CID" scf:"F_CID"\n",
            (isSurrogate?'s':'u'),
            CI->id,
            ScaffoldGraph->tigStore->getUnitigCoverageStat(CI->id),
            (int)CI->offsetAEnd.mean,
            CI->offsetAEnd.variance,
            (int)CI->offsetBEnd.mean,
            CI->offsetBEnd.variance,
            CI->AEndNext,
            CI->BEndNext,
            CI->info.CI.contigID,
            CI->scaffoldID);
#ifdef DEBUG_CONTIG
    fprintf(stderr,"\t%c:"F_CID" cov:%d aoff:%d boff:%d anxt:"F_CID" bnxt:"F_CID" ctg:"F_CID" scf:"F_CID"\n",
            (isSurrogate?'s':'u'),
            CI->id,
            ScaffoldGraph->tigStore->getUnitigCoverageStat(CI->id),
            (int)CI->offsetAEnd.mean,
            (int)CI->offsetBEnd.mean,
            CI->AEndNext,
            CI->BEndNext,
            CI->info.CI.contigID,
            CI->scaffoldID);
    fflush(stderr);
#endif
  }
  fprintf(stream,"\t%s Edges from A End:\n", (raw?" Raw ":" Merged "));

  {
    GraphEdgeIterator edges(graph->ContigGraph, contig->id, A_END, ALL_EDGES);
    CIEdgeT *edge;

    while((edge = (raw) ? edges.nextRaw() : edges.nextMerged()) != NULL)
      PrintContigEdgeInScfContext(stream,graph->ContigGraph,"\t",edge, contig->id);
  }

  fprintf(stream,"\t%s Edges from B End:\n", (raw?" Raw ":" Merged "));

  {
    GraphEdgeIterator edges(graph->ContigGraph, contig->id, B_END, ALL_EDGES);
    CIEdgeT *edge;

    while((edge = (raw) ? edges.nextRaw() : edges.nextMerged()) != NULL)
      PrintContigEdgeInScfContext(stream,graph->ContigGraph,"\t",edge, contig->id);
  }
}

void DumpContigs(FILE *stream, ScaffoldGraphT *graph, int raw){

  GraphNodeIterator CIs;
  NodeCGW_T *contig;

  fprintf(stream,"************Dumping Contigs ***********\n");

  InitGraphNodeIterator(&CIs, graph->ContigGraph, GRAPH_NODE_DEFAULT);
  while((contig = NextGraphNodeIterator(&CIs)) != NULL){
    DumpContig(stream, graph, contig, raw);
  }
}




void
BuildInitialContigs(ScaffoldGraphT *graph) {

  //  Resize the ContigGraph to the same size as the CI Graph

  fprintf(stderr,"BuildInitialContigs()-- converting %d unitigs with %d edges to contigs.\n",
          GetNumGraphNodes(graph->CIGraph),
          GetNumGraphEdges(graph->CIGraph));

  DeleteVA_NodeCGW_T(graph->ContigGraph->nodes);
  DeleteVA_EdgeCGW_T(graph->ContigGraph->edges);

  graph->ContigGraph->nodes = CreateVA_NodeCGW_T(GetNumGraphNodes(graph->CIGraph));
  graph->ContigGraph->edges = CreateVA_EdgeCGW_T(GetNumGraphEdges(graph->CIGraph));

  EnableRange_VA(graph->ContigGraph->nodes, GetNumGraphNodes(graph->CIGraph));

  graph->ContigGraph->edgeLists.clear();

  ResizeEdgeList(graph->ContigGraph);

  //  Clear contigs.

  for (int32 cid=0; cid < GetNumGraphNodes(graph->ContigGraph); cid++) {
    NodeCGW_T *ctg = GetGraphNode(graph->ContigGraph, cid);

    ctg->flags.all           = 0;
    ctg->flags.bits.isContig = TRUE;
    ctg->flags.bits.isDead   = TRUE;

    //ctg->edgeHead            = NULLINDEX;

    graph->ContigGraph->edgeLists[cid].clear();
  }

  //  And copy.

  GraphNodeIterator CIs;
  NodeCGW_T        *CI;

  InitGraphNodeIterator(&CIs, graph->CIGraph, GRAPH_NODE_DEFAULT);

  while ((CI = NextGraphNodeIterator(&CIs)) != NULL){
    assert(CI->flags.bits.isDead == 0);

    //  Reset the unitig.

    CI->AEndNext                = NULLINDEX;
    CI->BEndNext                = NULLINDEX;
    CI->info.CI.contigID        = CI->id;

    //  Copy to a new contig

    ContigT contig = *CI;

    contig.type                 = CONTIG_CGW;
    contig.id                   = CI->id;
    contig.scaffoldID           = NULLINDEX;
    contig.smoothExpectedCID    = NULLINDEX;
    contig.numEssentialA        = 0;
    contig.numEssentialB        = 0;
    contig.essentialEdgeA       = NULLINDEX;
    contig.essentialEdgeB       = NULLINDEX;
    contig.info.Contig.AEndCI   = CI->id;
    contig.info.Contig.BEndCI   = CI->id;
    contig.info.Contig.numCI    = 1;
    contig.indexInScaffold      = NULLINDEX;
    contig.flags.bits.isCI      = FALSE;
    contig.flags.bits.isContig  = TRUE;
    contig.flags.bits.isChaff   = CI->flags.bits.isChaff;
    contig.flags.bits.isClosure = CI->flags.bits.isClosure;
    //contig.edgeHead             = NULLINDEX;

    SetNodeCGW_T(graph->ContigGraph->nodes, contig.id, &contig);

    //  Ensure that there are no edges, and that the edgeList is allocated.
    assert(graph->ContigGraph->edgeLists[contig.id].empty() == true);
  }

  graph->numContigs = GetNumGraphNodes(graph->ContigGraph);

  //  Now, work on the edges.

  uint32   nRawSkipped = 0;
  uint32   nMerged     = 0;
  uint32   nTopRaw     = 0;
  uint32   nRaw        = 0;

  for (uint32 i=0; i<GetNumGraphEdges(graph->CIGraph); i++) {
    CIEdgeT  *edge = GetGraphEdge(graph->CIGraph, i);

    if (edge->flags.bits.isDeleted)
      continue;

    //  If this isn't a top-level edge, skip it.
    //  It must also be raw, and therefore already added.

    if (edge->topLevelEdge != GetVAIndex_CIEdgeT(graph->CIGraph->edges, edge)) {
      assert(edge->flags.bits.isRaw == true);
      nRawSkipped++;
      continue;
    }

    //  Is it a top-level raw edge?

    if (edge->flags.bits.isRaw == true) {
      CIEdgeT   newEdge     = *edge;

      newEdge.referenceEdge = i;
      newEdge.topLevelEdge  = GetNumGraphEdges(graph->ContigGraph);

      AppendGraphEdge(graph->ContigGraph, &newEdge);

      InsertGraphEdgeInList(graph->ContigGraph, newEdge.topLevelEdge, newEdge.idA);
      InsertGraphEdgeInList(graph->ContigGraph, newEdge.topLevelEdge, newEdge.idB);

      nTopRaw++;

      continue;
    }

    //  Otherwise, it must be a top-level merged edge

    assert(edge->nextRawEdge != NULLINDEX);

    if (edge->flags.bits.isRaw == FALSE) {
      CIEdgeT   newEdge     = *edge;
      CIEdgeT   rawEdge;

      newEdge.topLevelEdge  = GetNumGraphEdges(graph->ContigGraph);
      newEdge.nextRawEdge   = GetNumGraphEdges(graph->ContigGraph) + 1;  //  Must be raw edges!

      AppendGraphEdge(graph->ContigGraph, &newEdge);

      InsertGraphEdgeInList(graph->ContigGraph, newEdge.topLevelEdge, newEdge.idA);
      InsertGraphEdgeInList(graph->ContigGraph, newEdge.topLevelEdge, newEdge.idB);

      nMerged++;

      //  And copy over all the raw edges that compose this merged edge.

      while (edge->nextRawEdge != NULLINDEX) {
        CIEdgeT  *redge       = GetGraphEdge(graph->CIGraph, edge->nextRawEdge);  //  Grab the raw CI edge

        rawEdge               = *redge;

        //  These used to be assignments, but they should be correct as is.
        assert(rawEdge.idA == newEdge.idA);
        assert(rawEdge.idB == newEdge.idB);

        rawEdge.topLevelEdge  = newEdge.topLevelEdge;   //  The ID of the new contig top level edge
        rawEdge.referenceEdge = edge->nextRawEdge;      //  The ID of the current CI raw edge

        //  rawEdge.nextRawEdge is currently the next CI raw edge.  If that is defined,
        //  reset it to the next edge we'd add to the contig graph.
        if (rawEdge.nextRawEdge != NULLINDEX)
          rawEdge.nextRawEdge = GetNumGraphEdges(graph->ContigGraph) + 1;

        AppendGraphEdge(graph->ContigGraph, &rawEdge);

        nRaw++;

        edge = redge;
      }
    }
  }

  fprintf(stderr,"BuildInitialContigs()-- converted "F_U32" merged edges with "F_U32" raw edges; skipped "F_U32" raw edges in merged edges; converted "F_U32" top level raw edges.\n",
          nMerged, nRaw, nRawSkipped, nTopRaw);

  assert(nRawSkipped == nRaw);
}

int GetConsensus(GraphCGW_T *graph, CDS_CID_t CIindex,
                 VA_TYPE(char) *consensusVA, VA_TYPE(char) *qualityVA){
  // Return value is length of unitig or contig  sequence/quality (-1 if failure)
  ChunkInstanceT *CI = GetGraphNode(graph, CIindex);
  MultiAlignT *MA = NULL;

  ResetVA_char(consensusVA);
  ResetVA_char(qualityVA);
  if(CI->flags.bits.isCI){
    // Get it from the store of Unitig multi alignments
    MA = ScaffoldGraph->tigStore->loadMultiAlign(CIindex, TRUE);
  }else if(CI->flags.bits.isContig){// Get it from the store of Contig multi alignments
    assert(graph->type == CONTIG_GRAPH);
    MA = ScaffoldGraph->tigStore->loadMultiAlign(CIindex, FALSE);
  }else assert(0);

  GetMultiAlignUngappedConsensus(MA, consensusVA, qualityVA);

  return GetNumchars(consensusVA);
}


void SetCIScaffoldIds(ChunkInstanceT *CI, CDS_CID_t scaffoldID){
  // Set the scaffold ID of this CI
  CI->scaffoldID = scaffoldID;
  if(CI->flags.bits.isChaff){ // This can only happen once
    MultiAlignT *ma = ScaffoldGraph->tigStore->loadMultiAlign(CI->id, TRUE);
    CIFragT *frag = GetCIFragT(ScaffoldGraph->CIFrags, GetIntMultiPos(ma->f_list,0)->ident);
    assert(frag->flags.bits.isSingleton);
    frag->flags.bits.isChaff = FALSE;
    CI->flags.bits.isChaff = FALSE;
    if(GlobalData->debugLevel > 0)
      fprintf(stderr,"* SetCIScaffoldIDs ci "F_CID" and frag "F_CID" are NOT chaff\n",
              CI->id, frag->read_iid);
  }
}

void SetCIContigIds(ChunkInstanceT *CI, CDS_CID_t contigID){
  // Set the contig ID of this CI
  CI->info.CI.contigID = contigID;
}

void SetContigScaffoldIds(ContigT *contig, CDS_CID_t scaffoldID){
  ContigTIterator CIs;
  ChunkInstanceT *CI;

  // Set the scaffold ID of this CI
  contig->scaffoldID = scaffoldID;

  // Set the isUnique bit
  contig->flags.bits.isUnique = TRUE;
  contig->flags.bits.isChaff = FALSE;
  if(GlobalData->debugLevel > 0)
    fprintf(stderr,"* SetContigScaffoldIDs contig "F_CID" is NOT chaff\n",
            contig->id);

  // Set the scaffold ID of all of the contained CIs
  InitContigTIterator(ScaffoldGraph, contig->id, TRUE, FALSE, &CIs);
  while((CI = NextContigTIterator(&CIs)) != NULL){
    assert(CI->flags.bits.isCI);
    SetCIScaffoldIds(CI,scaffoldID);
  }
}



void CheckAllContigFragments(void){
  GraphNodeIterator contigs;
  NodeCGW_T *contig;

  InitGraphNodeIterator(&contigs, ScaffoldGraph->ContigGraph, GRAPH_NODE_DEFAULT);
  while((contig = NextGraphNodeIterator(&contigs)) != NULL){
    MultiAlignT *ma  = ScaffoldGraph->tigStore->loadMultiAlign(contig->id, FALSE);
    int i;
    if(!ma){
      fprintf(stderr,"*CheckAllContigFragments -- Contig "F_CID" is missing\n", contig->id);
      continue;
    }
    for(i = 0; i < GetNumIntMultiPoss(ma->f_list); i++){
      IntMultiPos *mp = GetIntMultiPos(ma->f_list,i);
      CIFragT *frag = GetCIFragT(ScaffoldGraph->CIFrags, mp->ident);
      assert(frag->contigID == contig->id);
    }
  }
}

CDS_CID_t GetOriginalContigID(CDS_CID_t contigID){
  NodeCGW_T *contig = GetGraphNode(ScaffoldGraph->ContigGraph, contigID);
  NodeCGW_T *ci;

  if(contig->info.Contig.numCI > 1)
    return contig->id;

  // Get the CI in this contig
  ci = GetGraphNode(ScaffoldGraph->CIGraph, contig->info.Contig.AEndCI);

  // If this is not a surrogate, we're done
  if(ci->type != RESOLVEDREPEATCHUNK_CGW)
    return contig->id;

  // Get the parent of the surrogate
  ci = GetGraphNode(ScaffoldGraph->CIGraph, ci->info.CI.baseID);

  // Return its contigID
  return ci->info.CI.contigID;

}


int IsDefinitelyUniqueContig(ContigT *contig){
  NodeCGW_T *ci;

  /*
    fprintf(stderr,"* IsDefinitelyUniqueContig "F_CID" numCI:%d\n",
    contig->id, contig->info.Contig.numCI);
  */
  // If it has been contigged already, it is unique
  if(contig->info.Contig.numCI > 1)
    return TRUE;

  // Get the CI in this contig
  ci = GetGraphNode(ScaffoldGraph->CIGraph, contig->info.Contig.AEndCI);
  /*
    fprintf(stderr,"* type = %d  coverage = %d  (%d)\n",
    ci->type, ScaffoldGraph->tigStore->getUnitigCoverageStat(CI->id),
    GlobalData->cgbDefinitelyUniqueCutoff);
  */

  // when the flag says we are unique, we always return true
  if (ScaffoldGraph->tigStore->getUnitigForceUnique(ci->id) == true) {
    return TRUE;
  }

  // when the flag says we are repeat, we always return false
  if (ScaffoldGraph->tigStore->getUnitigForceRepeat(ci->id) == true) {
    return FALSE;
  }

  // If this is not a surrogate, we're done
  return( ScaffoldGraph->tigStore->getUnitigCoverageStat(ci->id) > GlobalData->cgbDefinitelyUniqueCutoff);

}


int IsShakyContigAtScaffoldEnd(ContigT *contig){

  if(contig->numEssentialB <= 1 &&
     contig->numEssentialA <= 1){
    return FALSE;
  }

#if 0
  fprintf(stderr,"* IsShakyContigAtScaffoldEnd "F_CID" numA:%d numB:%d %s\n",
          contig->id,
          contig->numEssentialA,
          contig->numEssentialB,
          (IsDefinitelyUniqueContig(contig)?"unique":"shaky"));
#endif

  return(!IsDefinitelyUniqueContig(contig));
}
