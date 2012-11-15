
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

static char *rcsid = "$Id: GraphCGW_T.c,v 1.116 2012-11-15 02:17:45 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>

#include "AS_global.h"
#include "AS_UTL_fileIO.h"
#include "AS_UTL_interval.h"
#include "ScaffoldGraphIterator_CGW.h"
#include "ScaffoldGraph_CGW.h"
#include "GraphCGW_T.h"
#include "Instrument_CGW.h"
#include "MultiAlignment_CNS.h"
#include "UtilsREZ.h"

#include "Input_CGW.h"

#include "omp.h"

#include <vector>
#include <algorithm>

using namespace std;

VA_DEF(PtrT)

//  If enabled, output logging of the resulting merged edges.
#undef EDGELOG

void InitializeChunkInstance(ChunkInstanceT *ci, ChunkInstanceType type){
  ClearChunkInstance(ci);
  ci->flags.bits.isCI = TRUE;
  SetNodeType(ci, type);
}
void InitializeContig(ContigT *contig, ChunkInstanceType type){
  ClearChunkInstance(contig);
  contig->flags.bits.isContig = TRUE;
  SetNodeType(contig, type);
}
void InitializeScaffold(CIScaffoldT *scaffold, ChunkInstanceType type){
  ClearChunkInstance(scaffold);
  scaffold->flags.bits.isScaffold = TRUE;
  SetNodeType(scaffold, type);
}

int isDeadCIScaffoldT(CIScaffoldT *scaffold){
  return scaffold->flags.bits.isDead;
}


/* Make sure there is enough space for the number of edges we expect */
void ReallocGraphEdges(GraphCGW_T *graph, int32 numEdges){
  MakeRoom_VA(graph->edges, (size_t) numEdges);
}

GraphCGW_T *
CreateGraphCGW(GraphType type, int32 numNodes, int32 numEdges) {
  GraphCGW_T *graph = new GraphCGW_T;

  memset(graph, 0, sizeof(GraphCGW_T));

  graph->type  = type;

  graph->nodes = CreateVA_NodeCGW_T((size_t) MAX(1024, numNodes));
  graph->edges = CreateVA_EdgeCGW_T((size_t) MAX(1024, numNodes));

  graph->edgeLists.clear();

  graph->numActiveNodes   = 0;
  graph->numActiveEdges   = 0;

  graph->freeEdgeHead     = NULLINDEX;
  graph->tobeFreeEdgeHead = NULLINDEX;
  graph->freeNodeHead     = NULLINDEX;
  graph->tobeFreeNodeHead = NULLINDEX;
  graph->deadNodeHead     = NULLINDEX;

  return graph;
}








bool
utgGraphEdgeCompare(CDS_CID_t a, CDS_CID_t b) {
  return(edgeCompareForStoring(GetGraphEdge(ScaffoldGraph->CIGraph, a),
                               GetGraphEdge(ScaffoldGraph->CIGraph, b)));
}

bool
ctgGraphEdgeCompare(CDS_CID_t a, CDS_CID_t b) {
  return(edgeCompareForStoring(GetGraphEdge(ScaffoldGraph->ContigGraph, a),
                               GetGraphEdge(ScaffoldGraph->ContigGraph, b)));
}

bool
scfGraphEdgeCompare(CDS_CID_t a, CDS_CID_t b) {
  return(edgeCompareForStoring(GetGraphEdge(ScaffoldGraph->ScaffoldGraph, a),
                               GetGraphEdge(ScaffoldGraph->ScaffoldGraph, b)));
}


void ResizeEdgeList(GraphCGW_T *graph) {
  uint32   newSize = 0;

  if (graph->edgeLists.size() <= GetNumNodeCGW_Ts(graph->nodes))
    newSize = 2 * GetNumNodeCGW_Ts(graph->nodes) + 16;

  if (newSize == 0)
    return;

  fprintf(stderr, "ResizeEdgeList()-- sizing list for '%c' to "F_U32" (currently "F_SIZE_T")\n",
          graph->type, newSize, graph->edgeLists.size());

  if (graph->type == CI_GRAPH) {
    //set<CDS_CID_t,utgGraphEdgeCompare>  c;
    set<CDS_CID_t,bool(*)(CDS_CID_t,CDS_CID_t)> c(utgGraphEdgeCompare);
    graph->edgeLists.resize(newSize, c);
  }

  if (graph->type == CONTIG_GRAPH) {
    //set<CDS_CID_t,ctgGraphEdgeCompare>  c;
    set<CDS_CID_t,bool(*)(CDS_CID_t,CDS_CID_t)> c(ctgGraphEdgeCompare);
    graph->edgeLists.resize(newSize, c);
  }

  if (graph->type == SCAFFOLD_GRAPH) {
    //set<CDS_CID_t,scfGraphEdgeCompare>  c;
    set<CDS_CID_t,bool(*)(CDS_CID_t,CDS_CID_t)> c(scfGraphEdgeCompare);
    graph->edgeLists.resize(newSize, c);
  }
}








void SaveGraphCGWToStream(GraphCGW_T *graph, FILE *stream){
  int32 i;
  CopyToFileVA_NodeCGW_T(graph->nodes, stream);

  // Save lists of indicies, if they exist
  if(graph->type == CI_GRAPH)
    for(i = 0; i < GetNumGraphNodes(graph); i++){
      NodeCGW_T *node = GetGraphNode(graph,i);

      if (!node->flags.bits.isCI)
        continue;
      if (node->info.CI.numInstances < 3)
        continue;

      assert(node->info.CI.numInstances == GetNumCDS_CID_ts(node->info.CI.instances.va));

      CopyToFileVA_CDS_CID_t(node->info.CI.instances.va, stream);
    }

  CopyToFileVA_EdgeCGW_T(graph->edges, stream);

  AS_UTL_safeWrite(stream, &graph->type,             "SaveGraphCGWToStream", sizeof(int32),     1);
  AS_UTL_safeWrite(stream, &graph->numActiveNodes,   "SaveGraphCGWToStream", sizeof(int32),     1);
  AS_UTL_safeWrite(stream, &graph->numActiveEdges,   "SaveGraphCGWToStream", sizeof(int32),     1);
  AS_UTL_safeWrite(stream, &graph->freeEdgeHead,     "SaveGraphCGWToStream", sizeof(CDS_CID_t), 1);
  AS_UTL_safeWrite(stream, &graph->tobeFreeEdgeHead, "SaveGraphCGWToStream", sizeof(CDS_CID_t), 1);
  AS_UTL_safeWrite(stream, &graph->freeNodeHead,     "SaveGraphCGWToStream", sizeof(CDS_CID_t), 1);
  AS_UTL_safeWrite(stream, &graph->tobeFreeNodeHead, "SaveGraphCGWToStream", sizeof(CDS_CID_t), 1);
  AS_UTL_safeWrite(stream, &graph->deadNodeHead,     "SaveGraphCGWToStream", sizeof(CDS_CID_t), 1);  //  UNUSED
}


GraphCGW_T *LoadGraphCGWFromStream(FILE *stream){
  CDS_CID_t i;
  int       status = 0;
  GraphCGW_T *graph = new GraphCGW_T;

  memset(graph, 0, sizeof(GraphCGW_T));

  graph->nodes = CreateFromFileVA_NodeCGW_T(stream);

  // Load lists of indicies, if they exist
  for(i = 0; i < GetNumGraphNodes(graph); i++){
    NodeCGW_T *node = GetGraphNode(graph,i);

    if (!node->flags.bits.isCI)
      continue;
    if (node->info.CI.numInstances < 3)
      continue;

    node->info.CI.instances.va = CreateFromFileVA_CDS_CID_t(stream);

    assert(node->info.CI.numInstances == GetNumCDS_CID_ts(node->info.CI.instances.va));
  }

  graph->edges = CreateFromFileVA_EdgeCGW_T( stream);

  status  = AS_UTL_safeRead(stream, &graph->type,             "LoadGraphCGWFromStream", sizeof(int32),     1);
  status += AS_UTL_safeRead(stream, &graph->numActiveNodes,   "LoadGraphCGWFromStream", sizeof(int32),     1);
  status += AS_UTL_safeRead(stream, &graph->numActiveEdges,   "LoadGraphCGWFromStream", sizeof(int32),     1);
  status += AS_UTL_safeRead(stream, &graph->freeEdgeHead,     "LoadGraphCGWFromStream", sizeof(CDS_CID_t), 1);
  status += AS_UTL_safeRead(stream, &graph->tobeFreeEdgeHead, "LoadGraphCGWFromStream", sizeof(CDS_CID_t), 1);
  status += AS_UTL_safeRead(stream, &graph->freeNodeHead,     "LoadGraphCGWFromStream", sizeof(CDS_CID_t), 1);
  status += AS_UTL_safeRead(stream, &graph->tobeFreeNodeHead, "LoadGraphCGWFromStream", sizeof(CDS_CID_t), 1);
  status += AS_UTL_safeRead(stream, &graph->deadNodeHead,     "LoadGraphCGWFromStream", sizeof(CDS_CID_t), 1);  //  UNUSED
  assert(status == 8);

  return graph;
}


void
RebuildGraphEdges(GraphCGW_T *graph) {

  ResizeEdgeList(graph);

#if 0
  //  BPW couldn't get scaffold edges to load properly.  They triggered a ton of issues
  //  with null scaffolds and edges to deleted scaffolds.
  //
  //  We never restart in the middle of using scaffold edges; they're always rebuilt.
  //
  if (graph->type == SCAFFOLD_GRAPH) {
    fprintf(stderr, "RebuildGraphEdges()-- NOT rebuilding "F_SIZE_T" edges in "F_SIZE_T" %s.\n",
            GetNumEdgeCGW_Ts(graph->edges),
            GetNumEdgeCGW_Ts(graph->nodes),
            (graph->type == 'c') ? "unitigs" : ((graph->type == 'C') ? "contigs" : "scaffolds"));

    ResetEdgeCGW_T(graph->edges);

    return;
  }
#endif

  fprintf(stderr, "RebuildGraphEdges()-- Rebuilding "F_SIZE_T" edges in "F_SIZE_T" %s.\n",
          GetNumEdgeCGW_Ts(graph->edges),
          GetNumEdgeCGW_Ts(graph->nodes),
          (graph->type == 'c') ? "unitigs" : ((graph->type == 'C') ? "contigs" : "scaffolds"));

  for (uint32 ei=0; ei<GetNumEdgeCGW_Ts(graph->edges); ei++) {
    EdgeCGW_T  *edge = GetGraphEdge(graph, ei);

    if (edge->topLevelEdge != ei)
      continue;

    if (edge->flags.bits.isDeleted == true)
      continue;

    NodeCGW_T  *A = GetGraphNode(graph, edge->idA);
    NodeCGW_T  *B = GetGraphNode(graph, edge->idB);

    if ((A == NULL) || (B == NULL)) {
      //  Who left this edge in here?
      fprintf(stderr, "WARNING: edge eid %d cid %d %c ptr %p cid %d %c ptr %p\n",
              ei,
              edge->idA, graph->type, A,
              edge->idB, graph->type, B);
    }
    assert(A != NULL);
    assert(B != NULL);

    if ((A->flags.bits.isDead == true) ||
        (B->flags.bits.isDead == true)) {
      //  Who left this edge in here?
      fprintf(stderr, "WARNING: edge eid %d cid %d %c deleted %d cid %d %c deleted %d\n",
              ei,
              edge->idA, graph->type, A->flags.bits.isDead,
              edge->idB, graph->type, B->flags.bits.isDead);
    }
    assert(A->flags.bits.isDead == false);
    assert(B->flags.bits.isDead == false);


    InsertGraphEdgeInList(graph, ei, edge->idA);
    InsertGraphEdgeInList(graph, ei, edge->idB);
  }
}


void DeleteGraphCGW(GraphCGW_T *graph){
  int32 i;

  // Delete lists of indicies, if they exist
  if(graph->type == CI_GRAPH)
    for(i = 0; i < GetNumGraphNodes(graph); i++){
      NodeCGW_T *node = GetGraphNode(graph,i);
      if (!node->flags.bits.isCI)
        continue;
      if (node->info.CI.numInstances < 3)
        continue;
      DeleteVA_CDS_CID_t(node->info.CI.instances.va);
    }

  DeleteVA_NodeCGW_T(graph->nodes);
  DeleteVA_EdgeCGW_T(graph->edges);

  delete graph;
}

/* Diagnostic */
size_t ReportMemorySizeGraphCGW(GraphCGW_T *graph, FILE *stream){
  char *nodeName = NULL;
  char *edgeName = NULL;
  size_t totalMemorySize = 0;

  switch(graph->type){
    case CI_GRAPH:
      nodeName = "CIs    ";
      edgeName = "CIEdges";
      break;
    case CONTIG_GRAPH:
      nodeName = "Contigs    ";
      edgeName = "ConEdges";
      break;
    case SCAFFOLD_GRAPH:
      nodeName = "Scaffolds";
      edgeName = "SEdges";
      break;
    default:
      assert(0);
  }
  totalMemorySize += ReportMemorySize_VA(graph->nodes, nodeName, stream);
  totalMemorySize += ReportMemorySize_VA(graph->edges, edgeName, stream);
  return totalMemorySize;
}




void
DumpGraphEdges(GraphCGW_T *graph, char *outname) {
  FILE *outfile = fopen(outname, "w");

  fprintf(stderr, "DumpGraphEdges()-- to '%s'\n", outname);

  for (uint32 i=0; i<GetNumGraphEdges(graph); i++) {
    EdgeCGW_T  *edge = GetGraphEdge(graph, i);

    fprintf(outfile, "eid "F_U32" deleted %d triplet "F_U32" "F_U32" %c contributing %d %.0f +- %.8f frag "F_U32":"F_U32" nextRaw "F_U32" topLevel "F_U32"\n",
            i,
            edge->flags.bits.isDeleted,
            edge->idA, edge->idB, edge->orient.toLetter(),
            edge->edgesContributing,
            edge->distance.mean, edge->distance.variance,
            edge->fragA, edge->fragB,
            edge->nextRawEdge,
            edge->topLevelEdge);
  }

  for (uint32 i=0; i<GetNumGraphNodes(graph); i++) {
    GraphEdgeIterator  edges(graph, i, A_END, ALL_EDGES);
    EdgeCGW_T         *edge;

    while ((edge = edges.nextMerged()) != NULL) {
    }
    edges.reset();
    while ((edge = edges.nextRaw()) != NULL) {
    }
  }

  for (uint32 i=0; i<GetNumGraphNodes(graph); i++) {
    GraphEdgeIterator  edges(graph, i, B_END, ALL_EDGES);
    EdgeCGW_T         *edge;

    while ((edge = edges.nextMerged()) != NULL) {
    }
    edges.reset();
    while ((edge = edges.nextRaw()) != NULL) {
    }
  }

  for (uint32 i=0; i<GetNumGraphNodes(graph); i++) {
    NodeCGW_T  *node = GetGraphNode(graph, i);

    GraphEdgeIterator  edges(graph, i, A_END, ALL_EDGES);
    EdgeCGW_T         *edge;

    while ((edge = edges.nextMerged()) != NULL) {
      fprintf(outfile, "mergedEdge node "F_U32" eid "F_SIZE_T" triplet "F_U32" "F_U32" %c\n",
              i,
              GetVAIndex_EdgeCGW_T(graph->edges, edge),
              edge->idA, edge->idB, edge->orient.toLetter());
    }
  }

  for (uint32 i=0; i<GetNumGraphNodes(graph); i++) {
    NodeCGW_T  *node = GetGraphNode(graph, i);

    GraphEdgeIterator  edges(graph, i, B_END, ALL_EDGES);
    EdgeCGW_T         *edge;

    while ((edge = edges.nextRaw()) != NULL) {
      fprintf(outfile, "rawEdge node "F_U32" eid "F_SIZE_T" triplet "F_U32" "F_U32" %c\n",
              i,
              GetVAIndex_EdgeCGW_T(graph->edges, edge),
              edge->idA, edge->idB, edge->orient.toLetter());
    }
  }

  fclose(outfile);
}


/* Check that edge with index eid is properly wired in the graph:
   - we find it when looking for it in both lists which it is supposed to be a member
   - we can walk backwards from the edge to the heads of the two lists
*/
void
GraphEdgeSanity(GraphCGW_T *graph, CDS_CID_t eid) {
  CIEdgeT  *edge      = GetGraphEdge(graph, eid);
  CIEdgeT  *edgeFromA = NULL;
  CIEdgeT  *edgeFromB = NULL;
  bool      fails     = false;

  {
    GraphEdgeIterator edges(graph, edge->idA, ALL_END, ALL_EDGES);
    while (NULL != (edgeFromA = edges.nextMerged())){
      //fprintf(stderr, "GraphEdgeSanity()--  testA Found edge %16p idx %d (%d,%d)\n",
      //        edgeFromA, GetVAIndex_EdgeCGW_T(graph->edges, edgeFromA), edgeFromA->idA, edgeFromA->idB);
      if (edgeFromA == edge)
        break;
    }
    if (edgeFromA == NULL) {
      fprintf(stderr,"GraphEdgeSanity()-- Couldn't find edge %16p "F_CID" ("F_CID","F_CID") looking from A "F_CID"\n",
              edge, eid, edge->idA, edge->idB, edge->idA);
      edges.reset();
      while (NULL != (edgeFromA = edges.nextMerged())){
        fprintf(stderr, "GraphEdgeSanity()--  Found edge %16p idx "F_SIZE_T" (%d,%d)\n",
                edgeFromA, GetVAIndex_EdgeCGW_T(graph->edges, edgeFromA), edgeFromA->idA, edgeFromA->idB);
      }
      fails = true;
    }
  }

  {
    GraphEdgeIterator edges(graph, edge->idB, ALL_END, ALL_EDGES);
    while (NULL != (edgeFromB = edges.nextMerged())){
      //fprintf(stderr, "GraphEdgeSanity()--  testB Found edge %16p idx %d (%d,%d)\n",
      //        edgeFromB, GetVAIndex_EdgeCGW_T(graph->edges, edgeFromB), edgeFromB->idA, edgeFromB->idB);
      if (edgeFromB == edge)
        break;
    }
    if(edgeFromB == NULL) {
      fprintf(stderr,"GraphEdgeSanity()-- Couldn't find edge %16p "F_CID" ("F_CID","F_CID") looking from B "F_CID"\n",
              edge, eid, edge->idA, edge->idB, edge->idB);
      edges.reset();
      while (NULL != (edgeFromB = edges.nextMerged())){
        fprintf(stderr, "GraphEdgeSanity()--  Found edge %16p idx "F_SIZE_T" (%d,%d)\n",
                edgeFromB, GetVAIndex_EdgeCGW_T(graph->edges, edgeFromB), edgeFromB->idA, edgeFromB->idB);
      }
      fails = true;
    }
  }

  assert(fails == false);

  //PrintGraphEdge(stderr,graph," ", edge, idA);
}



void
InsertGraphEdgeInList(GraphCGW_T *graph,
                      CDS_CID_t   edgeID,
                      CDS_CID_t   ciID){

  CIEdgeT *newEdge = GetGraphEdge(graph, edgeID);

  assert(newEdge->idA == ciID || newEdge->idB == ciID);
  assert(newEdge->flags.bits.isDeleted == false);

  ChunkInstanceT *ci = GetGraphNode(graph, ciID);

  if (ci->flags.bits.isDead) {
    fprintf(stderr, "InsertGraphEdgeInList()--  CI="F_CID" isDead!\n", ciID);
    assert(!ci->flags.bits.isDead);
    return;
  }

  //  Shouldn't be needed here.  If the node exists, the edgeList must exist.
  //ResizeEdgeList(graph);

  //  Might need to keep edges sorted using edgeCompare, after grabbing the
  //  CIEdgeT from the graph.

  //fprintf(stderr, "InsertGraphEdgeInList()-- edge eid %d cid %d %c cid %d %c orient %c sloppy %d overlap %d distance %f ciID %d %c\n",
  //        edgeID, newEdge->idA, graph->type, newEdge->idB, graph->type, newEdge->orient.toLetter(),
  //        isSloppyEdge(newEdge), isOverlapEdge(newEdge), newEdge->distance.mean,
  //        ciID, graph->type);

  if (graph->edgeLists[ciID].count(edgeID) > 0) {
    fprintf(stderr, "InsertGraphEdgeInList()--  WARNING:  Edge %d cid %d %c cid %d %c orient %c ciID %d %c already exists.\n",
            edgeID, newEdge->idA, graph->type, newEdge->idB, graph->type, newEdge->orient.toLetter(), ciID, graph->type);
    for (set<CDS_CID_t>::iterator  it = graph->edgeLists[ciID].begin(); it != graph->edgeLists[ciID].end(); it++) {
      EdgeCGW_T *e = GetGraphEdge(graph, *it);
      fprintf(stderr, "eid "F_U32" %d-%d orient %c sloppy %d overlap %d distance %f\n",
              *it, e->idA, e->idB, e->orient.toLetter(), isSloppyEdge(e), isOverlapEdge(e), e->distance.mean);
    }
  }
  assert(graph->edgeLists[ciID].count(edgeID) == 0);

  graph->edgeLists[ciID].insert(edgeID);
}


// Initialize the status flags for the given edge.
void InitGraphEdgeFlags(GraphCGW_T *graph, EdgeCGW_T *edge){
  ChunkInstanceT *CIA = GetGraphNode(graph, edge->idA);
  ChunkInstanceT *CIB = GetGraphNode(graph, edge->idB);

  assert(edge && (edge->idA < edge->idB) && CIA && CIB);

  edge->flags.bits.isActive = FALSE;
  edge->flags.bits.isUniquetoUnique = FALSE;
  edge->flags.bits.isConfirmed = FALSE;

  // Set the isUniquetoUnique flag for edges where both CIA and CIB are unique.
  if(CIA->flags.bits.isUnique && CIB->flags.bits.isUnique){
    edge->flags.bits.isUniquetoUnique = TRUE;

    // Set the isActive flag for edges which are not overlap only
    // Don't include high variance edges
    // **SAK** These edges to contained contigs were confusing the transitive reduction code
    if((isOverlapEdge(edge) && (edge->edgesContributing == 1)) ||
       edge->flags.bits.isProbablyBogus ||
       isSloppyEdge(edge) ||
       edge->distance.mean < -1000
       ){
      edge->flags.bits.isActive = FALSE;
      edge->flags.bits.isConfirmed = FALSE;
    }else{
      edge->flags.bits.isActive = TRUE;
      // Set the isConfirmed flag for edges with more than one pair of mates
      if((edge->edgesContributing > MIN_EDGES) ||
         (!isOverlapEdge(edge) && (edge->edgesContributing == MIN_EDGES))){
        edge->flags.bits.isConfirmed = TRUE;
      }else{
        edge->flags.bits.isConfirmed = FALSE;
      }
    }
  }

  edge->flags.bits.isEssential = FALSE;
  edge->flags.bits.wasEssential = FALSE;
  edge->flags.bits.isInferred = FALSE;
  edge->flags.bits.isInferredRemoved = FALSE;
  edge->flags.bits.isRedundantRemoved = FALSE;
  edge->flags.bits.isTransitivelyRemoved = FALSE;
  edge->flags.bits.isDeleted = FALSE;
}

// Insert a copy of the tobeInserted Mate Edge
// Returns index of CIEdge if successful
// Edge must be 'canonical' (idA < idB)
//
void
InsertGraphEdge(GraphCGW_T *graph, CDS_CID_t edgeID) {
  EdgeCGW_T *edge = GetGraphEdge(graph, edgeID);

  assert(GetEdgeStatus(edge) != INVALID_EDGE_STATUS);
  assert(edge->idA < edge->idB);
  assert(edge->topLevelEdge == edgeID);

  InitGraphEdgeFlags(graph, edge);

  InsertGraphEdgeInList(graph, edgeID, edge->idA);
  InsertGraphEdgeInList(graph, edgeID, edge->idB);
}


void
FreeGraphEdge(GraphCGW_T *graph,  EdgeCGW_T *edge) {
  CDS_CID_t eid = GetVAIndex_EdgeCGW_T(graph->edges, edge);

  assert(edge->flags.bits.isDeleted == false);

  //fprintf(stderr, "FreeGraphEdgeByEID()-- edge eid %d cid %d %c cid %d %c orient %c sloppy %d overlap %d distance %f\n",
  //        eid, edge->idA, graph->type, edge->idB, graph->type, edge->orient.toLetter(),
  //        isSloppyEdge(edge), isOverlapEdge(edge), edge->distance.mean);

  if (edge->nextRawEdge != NULLINDEX) {
    //assert(edge->flags.bits.isRaw == false);
    EdgeCGW_T *rawEdge = GetGraphEdge(graph, edge->nextRawEdge);

    while (rawEdge) {
      CDS_CID_t  rid     = GetVAIndex_EdgeCGW_T(graph->edges, rawEdge);
      CDS_CID_t  nxt     = rawEdge->nextRawEdge;

      assert(rawEdge->idA    == edge->idA);
      assert(rawEdge->idB    == edge->idB);
      assert(rawEdge->orient == edge->orient);

      assert(graph->edgeLists[edge->idA].count(rid) == 0);
      assert(graph->edgeLists[edge->idB].count(rid) == 0);

      assert(eid != graph->tobeFreeEdgeHead);

      rawEdge->nextRawEdge            = NULLINDEX;
      rawEdge->flags.bits.isDeleted   = TRUE;

      rawEdge->referenceEdge          = graph->tobeFreeEdgeHead;
      rawEdge->nextRawEdge            = NULLINDEX;
      rawEdge->flags.bits.isRaw       = TRUE;
      rawEdge->idA                    = NULLINDEX;
      rawEdge->idB                    = NULLINDEX;

      graph->tobeFreeEdgeHead         = rid;

      rawEdge = GetGraphEdge(graph, nxt);
    }
  }

  //  This is done by Unlink.
  //graph->edgeList[edge->idA].erase(eid);
  //graph->edgeList[edge->idB].erase(eid);

  assert(graph->edgeLists[edge->idA].count(eid) == 0);
  assert(graph->edgeLists[edge->idB].count(eid) == 0);

  assert(eid != graph->tobeFreeEdgeHead);

  edge->flags.bits.isDeleted  = TRUE;
  edge->referenceEdge         = graph->tobeFreeEdgeHead;
  edge->nextRawEdge           = NULLINDEX;
  edge->flags.bits.isRaw      = TRUE;
  edge->idA                   = NULLINDEX;
  edge->idB                   = NULLINDEX;

  graph->tobeFreeEdgeHead     = eid;
}


// Initialize a Graph Edge
void InitGraphEdge(EdgeCGW_T *edge){
  memset(edge, 0, sizeof(EdgeCGW_T));

  edge->idA           = NULLINDEX;
  edge->idB           = NULLINDEX;
  edge->flags.all     = 0;
  edge->quality       = 1.0;
  edge->nextRawEdge   = NULLINDEX;
  edge->referenceEdge = NULLINDEX;
  edge->topLevelEdge  = NULLINDEX;
}

// Get the edge from the free list
EdgeCGW_T *GetFreeGraphEdge(GraphCGW_T *graph){
  EdgeCGW_T *freeEdge = GetGraphEdge(graph, graph->freeEdgeHead);

  //fprintf(stderr, "GetFreeGraphEdge()-- %d %d\n", graph->freeEdgeHead, GetNumGraphEdges(graph));

  if (freeEdge == NULL) {
    EdgeCGW_T edge;
    memset(&edge, 0, sizeof(EdgeCGW_T));
    AppendEdgeCGW_T(graph->edges, &edge);
    freeEdge = GetGraphEdge(graph, GetNumGraphEdges(graph) - 1);
  } else {
    graph->freeEdgeHead = freeEdge->referenceEdge;
  }

  InitGraphEdge(freeEdge);

  freeEdge->flags.all         = 0;
  freeEdge->flags.bits.isRaw  = TRUE;
  freeEdge->nextRawEdge       = NULLINDEX;
  freeEdge->edgesContributing = 1;
  freeEdge->referenceEdge     = NULLINDEX;
  freeEdge->topLevelEdge      = GetVAIndex_EdgeCGW_T(graph->edges, freeEdge);

  SetGraphEdgeStatus(graph, freeEdge, UNKNOWN_EDGE_STATUS);

  return freeEdge;
}




// Delete a graph node and all edges incident on this node
void DeleteGraphNode(GraphCGW_T *graph, NodeCGW_T *node){
  GraphEdgeIterator edges(graph, node->id, ALL_END, ALL_EDGES);
  EdgeCGW_T        *edge;

  //  Remove edges from the graph
  while (NULL != (edge = edges.nextMerged()))
    DeleteGraphEdge(graph, edge);

  //  Remove edges from the scaffold
  graph->edgeLists[node->id].clear();

  //  Mark the node dead and
  node->flags.bits.isDead = TRUE;

  //  Unreference the consensus for this contig
  if(graph->type != SCAFFOLD_GRAPH)
    ScaffoldGraph->tigStore->deleteMultiAlign(node->id, graph->type == CI_GRAPH);

}


void UnlinkGraphEdge(GraphCGW_T *graph, EdgeCGW_T *edge){

  //  If this is a top level edge, deletion is straight forward, just mark
  //  the edge as deleted and remove from both sets.

  CDS_CID_t eid = GetVAIndex_EdgeCGW_T(graph->edges, edge);

  //fprintf(stderr, "UnlinkGraphEdge()-- edge eid %d cid %d %c cid %d %c orient %c sloppy %d overlap %d distance %f\n",
  //        eid, edge->idA, graph->type, edge->idB, graph->type, edge->orient.toLetter(), isSloppyEdge(edge), isOverlapEdge(edge), edge->distance.mean);

  assert(edge->flags.bits.isDeleted == false);

  if (edge->topLevelEdge != eid)
    fprintf(stderr, "ERROR:  Attempt to delete a raw edge from a merged edge.\n");
  assert (edge->topLevelEdge == eid);

  uint32 Aok = graph->edgeLists[edge->idA].erase(eid);
  uint32 Bok = graph->edgeLists[edge->idB].erase(eid);
    
  if ((Aok != 1) || (Bok != 1)) {
    fprintf(stderr, "UnlinkGraphEdge()-- WARNING: already unlinked eid %d cid %d %c cid %d %c orient %c sloppy %d overlap %d distance %f (%d %d) %p\n",
            eid, edge->idA, graph->type, edge->idB, graph->type, edge->orient.toLetter(), isSloppyEdge(edge), isOverlapEdge(edge), edge->distance.mean, Aok, Bok, edge);
    fprintf(stderr, "A Edges\n");
    for (set<CDS_CID_t>::iterator  it = graph->edgeLists[edge->idA].begin(); it != graph->edgeLists[edge->idA].end(); it++) {
      EdgeCGW_T *e = GetGraphEdge(graph, *it);
      fprintf(stderr, "eid "F_U32" %d-%d orient %c sloppy %d overlap %d distance %f %p compare %d %d\n",
              *it, e->idA, e->idB, e->orient.toLetter(), isSloppyEdge(e), isOverlapEdge(e), e->distance.mean, e,
              edgeCompareForStoring(edge, e), edgeCompareForStoring(e, edge));
    }
    fprintf(stderr, "B Edges\n");
    for (set<CDS_CID_t>::iterator  it = graph->edgeLists[edge->idB].begin(); it != graph->edgeLists[edge->idB].end(); it++) {
      EdgeCGW_T *e = GetGraphEdge(graph, *it);
      fprintf(stderr, "eid "F_U32" %d-%d orient %c sloppy %d overlap %d distance %f %p compare %d %d\n",
              *it, e->idA, e->idB, e->orient.toLetter(), isSloppyEdge(e), isOverlapEdge(e), e->distance.mean, e,
              edgeCompareForStoring(edge, e), edgeCompareForStoring(e, edge));
    }
  }
  assert(Aok == 1);
  assert(Bok == 1);
}



// Delete a top level CIEdgeT by unlinking it from the graph - remember
// this will unlink any dangling raw edges as well - and mark it as deleted.
void  DeleteGraphEdge(GraphCGW_T *graph,  EdgeCGW_T *edge){

  assert(edge->flags.bits.isDeleted == false);

  UnlinkGraphEdge(graph, edge);
  FreeGraphEdge(graph,edge);
}


void PrintGraphEdge(FILE *fp, GraphCGW_T *graph,
                    const char *label, EdgeCGW_T *edge, CDS_CID_t cid){
  char actualOverlap[256];
  int32 actual = 0;
  int32 delta = 0;
  char *flagTrans = "  ";
  char flagbuf[32];
  ChunkInstanceT *ChunkInstanceA = GetGraphNode(graph, edge->idA);
  ChunkInstanceT *ChunkInstanceB = GetGraphNode(graph, edge->idB);
  CDS_CID_t eid = GetVAIndex_EdgeCGW_T(graph->edges, edge);

  if(edge->flags.bits.isDeleted){
    fprintf(fp,"***FOLLOWING EDGE IS DELETED!!!!***\n");
  }
  assert(ChunkInstanceA && ChunkInstanceB);

  *actualOverlap = '\0';

  if(edge->flags.bits.isBogus){
    if(edge->flags.bits.isProbablyBogus)
      strcpy(actualOverlap," *Bogus and Prob Bogus*");
    else
      strcpy(actualOverlap," *Bogus*");
  }
  strcpy(flagbuf,"");
  if(edge->flags.bits.hasContributingOverlap){
    if(edge->flags.bits.isPossibleChimera)
      strcat(flagbuf,"$?");
    else
      strcat(flagbuf,"$0");
  }
  if(edge->flags.bits.hasTransChunk){
    strcat(flagbuf,"$t");
  }
  if(edge->flags.bits.aContainsB){
    strcat(flagbuf,"$C");
  }else if(edge->flags.bits.bContainsA){
    strcat(flagbuf,"$I");
  }
  if(edge->flags.bits.MeanChangedByWalking){
    strcat(flagbuf,"@M");
  }
  if(edge->flags.bits.isEssential){
    if(edge->flags.bits.isInferred){
      flagTrans = "&I";
    }else{
      flagTrans = "&E";
    }
  }else if(edge->flags.bits.isTransitivelyRemoved){
    flagTrans = "&T";
  }else if(edge->flags.bits.isRedundantRemoved){
    flagTrans = "&R";
  }else if(edge->flags.bits.isInferredRemoved){
    flagTrans = "&t";
  }else if(edge->flags.bits.isActive){
    if(edge->flags.bits.isConfirmed){
      flagTrans = "&C";
    }else{
      flagTrans = "&A";
    }
  }else if(edge->flags.bits.isProbablyBogus){
    flagTrans = "&B";
  }

  fprintf(fp,"%s eid:"F_CID" A:"F_CID" B:"F_CID" wgt:%d %s %s %s %s%s ori:%c qua:%3.2g trstd:%d con:%d dst:%d std:%g %s ",
          label,
          eid,
          edge->idA, edge->idB, edge->edgesContributing, flagbuf,
          (edge->flags.bits.hasExtremalAFrag?"$A":"  "),
          (edge->flags.bits.hasExtremalBFrag?"$B":"  "),
          flagTrans, (edge->flags.bits.isLeastSquares?"L":" "),
          GetEdgeOrientationWRT(edge, edge->idA).toLetter(),
          edge->quality, edge->flags.bits.edgeStatus & TRUSTED_EDGE_STATUS,
          edge->flags.bits.hasContainmentOverlap,
          (int)edge->distance.mean,
          (edge->distance.variance > .001?sqrt(edge->distance.variance):0.0),
          actualOverlap);

  if(edge->flags.bits.isRaw && !isOverlapEdge(edge))
    fprintf(fp,"("F_CID","F_CID")", edge->fragA, edge->fragB);
  fprintf(fp,"\n");
}


void PrintContigEdgeInScfContext(FILE *fp, GraphCGW_T *graph,
                                 char *label, EdgeCGW_T *edge,
                                 CDS_CID_t cid){
  char actualOverlap[256];
  int32 actual = 0;
  int32 delta = 0;
  char *flagTrans = "  ";
  char flagbuf[32];
  ChunkInstanceT *ChunkInstanceA = GetGraphNode(graph, edge->idA);
  ChunkInstanceT *ChunkInstanceB = GetGraphNode(graph, edge->idB);
  CDS_CID_t eid = GetVAIndex_EdgeCGW_T(graph->edges, edge);

  if(edge->flags.bits.isDeleted){
    fprintf(fp,"***FOLLOWING EDGE IS DELETED!!!!***\n");
  }
  assert(ChunkInstanceA && ChunkInstanceB);

  *actualOverlap = '\0';

  if(edge->flags.bits.isBogus){
    if(edge->flags.bits.isProbablyBogus)
      strcpy(actualOverlap," *Bogus and Prob Bogus*");
    else
      strcpy(actualOverlap," *Bogus*");
  }
  strcpy(flagbuf,"");
  if(edge->flags.bits.hasContributingOverlap){
    if(edge->flags.bits.isPossibleChimera)
      strcat(flagbuf,"$?");
    else
      strcat(flagbuf,"$0");
  }
  if(edge->flags.bits.hasTransChunk){
    strcat(flagbuf,"$t");
  }
  if(edge->flags.bits.aContainsB){
    strcat(flagbuf,"$C");
  }else if(edge->flags.bits.bContainsA){
    strcat(flagbuf,"$I");
  }
  if(edge->flags.bits.MeanChangedByWalking){
    strcat(flagbuf,"@M");
  }
  if(edge->flags.bits.isEssential){
    if(edge->flags.bits.isInferred){
      flagTrans = "&I";
    }else{
      flagTrans = "&E";
    }
  }else if(edge->flags.bits.isTransitivelyRemoved){
    flagTrans = "&T";
  }else if(edge->flags.bits.isRedundantRemoved){
    flagTrans = "&R";
  }else if(edge->flags.bits.isInferredRemoved){
    flagTrans = "&t";
  }else if(edge->flags.bits.isActive){
    if(edge->flags.bits.isConfirmed){
      flagTrans = "&C";
    }else{
      flagTrans = "&A";
    }
  }else if(edge->flags.bits.isProbablyBogus){
    flagTrans = "&B";
  }

  fprintf(fp,"%s eid:"F_CID" A:"F_CID" B:"F_CID" wgt:%d %s %s %s %s%s ori:%c qua:%3.2g trstd:%d con:%d dst:%d std:%g %s sameScf: %c Apos [%d,%d] Bpos [%d,%d]",
          label,
          eid,
          edge->idA, edge->idB, edge->edgesContributing, flagbuf,
          (edge->flags.bits.hasExtremalAFrag?"$A":"  "),
          (edge->flags.bits.hasExtremalBFrag?"$B":"  "),
          flagTrans, (edge->flags.bits.isLeastSquares?"L":" "),
          GetEdgeOrientationWRT(edge, edge->idA).toLetter(),
          edge->quality, edge->flags.bits.edgeStatus & TRUSTED_EDGE_STATUS,
          edge->flags.bits.hasContainmentOverlap,
          (int)edge->distance.mean,
          (edge->distance.variance > .001?sqrt(edge->distance.variance):0.0),
          actualOverlap,
          ( ChunkInstanceA->scaffoldID == ChunkInstanceB->scaffoldID ? 'T' : 'F'),
          (int) ChunkInstanceA->offsetAEnd.mean,
          (int) ChunkInstanceA->offsetBEnd.mean,
          (int) ChunkInstanceB->offsetAEnd.mean,
          (int) ChunkInstanceB->offsetBEnd.mean	);

  if(edge->flags.bits.isRaw && !isOverlapEdge(edge))
    fprintf(fp,"("F_CID","F_CID")", edge->fragA, edge->fragB);
  fprintf(fp,"\n");
}



CDS_CID_t
AddGraphEdge(GraphCGW_T *graph,
             CDS_CID_t cidA, CDS_CID_t cidB,
             CDS_CID_t fragidA, CDS_CID_t fragidB,
             CDS_CID_t dist,
             LengthT distance,
             double   quality,
             int32 fudgeDistance,
             PairOrient orientation,
             int isInducedByUnknownOrientation,
             int isOverlap,
             int isAContainsB,
             int isBContainsA,
             int isTransChunk,
             int isExtremalA,
             int isExtremalB,
             EdgeStatus status,
             int collectOverlap,
             int insert){
  CIEdgeT *ciedge = GetFreeGraphEdge(graph);
  CDS_CID_t ciedgeIndex = GetVAIndex_EdgeCGW_T(graph->edges, ciedge);
  CDS_CID_t idA, idB;
  CDS_CID_t fidA, fidB;
  PairOrient orient;
  ChunkInstanceT *CIa, *CIb;
  int aContainsB, bContainsA;
  int extremalA, extremalB;

  assert (dist != 0);

  /* First make this into a canonical format edge */
  idA = cidA;
  idB = cidB;

  fidA = fragidA;
  fidB = fragidB;
  orient = orientation;
  extremalA = isExtremalA;
  extremalB = isExtremalB;
  bContainsA = isBContainsA;
  aContainsB = isAContainsB;
  if(idA > idB){
    idA = cidB;
    idB = cidA;
    fidA = fragidB;
    fidB = fragidA;
    extremalA = isExtremalB;
    extremalB = isExtremalA;
    aContainsB = isBContainsA;
    bContainsA = isAContainsB;
    orient.flip();
  }

  CIa = GetGraphNode(graph, idA);
  CIb = GetGraphNode(graph, idB);

  assert(distance.variance >= 0.0);

  ciedge->idA = idA;
  ciedge->idB = idB;
  ciedge->fragA = fidA;
  ciedge->fragB = fidB;
  ciedge->distIndex = dist;
  ciedge->edgesContributing = 1;
  ciedge->distance = distance;
  ciedge->quality = quality;
  ciedge->minDistance = distance.mean - 3 * sqrt(distance.variance);
  ciedge->orient = orient;
  ciedge->flags.all = 0;
  ciedge->flags.bits.isRaw = TRUE;
  ciedge->nextRawEdge = NULLINDEX;
  ciedge->referenceEdge = NULLINDEX;
  ciedge->topLevelEdge = ciedgeIndex;
  ciedge->flags.bits.inducedByUnknownOrientation = isInducedByUnknownOrientation;
  ciedge->flags.bits.hasContributingOverlap = isOverlap;
  ciedge->flags.bits.mustOverlap = isOverlap ; /// TEMPORARY
  ciedge->flags.bits.aContainsB = aContainsB;
  ciedge->flags.bits.bContainsA = bContainsA;
  ciedge->flags.bits.hasTransChunk = isTransChunk;
  ciedge->flags.bits.hasExtremalAFrag = extremalA;
  ciedge->flags.bits.hasExtremalBFrag = extremalB;
  ciedge->flags.bits.hasContainmentOverlap = aContainsB | bContainsA;
  ciedge->flags.bits.isSloppy = (distance.variance > SLOPPY_EDGE_VARIANCE_THRESHHOLD);

  SetGraphEdgeStatus(graph, ciedge, status);

  ciedge->flags.bits.isProbablyBogus = FALSE;

  //  Determine whether this is potentally bogus edge
  //
  if(ciedge->fragA != NULLINDEX && ciedge->fragB != NULLINDEX){
    // Not an overlap fragment

    // Determine whether this LOOKS like a bogus edge.  Our criteria are:
    //  1) Implied overlap is > 500 base pairs + 5 sigma
    //  2) Implied distance is outside mean of distribution + 5 sigma

    DistT *distRecord = GetDistT(ScaffoldGraph->Dists, dist);
    int32 minDistance = -500 - distRecord->sigma * 5;
    int32 maxDistance = distRecord->mu + distRecord->sigma * 5;

    //  The isBogus flag is from the now-dead simulator, but it's
    //  overloaded in scaffold merging.

    ciedge->flags.bits.isBogus         = FALSE;
    ciedge->flags.bits.isProbablyBogus = (ciedge->distance.mean < minDistance ||
                                                                  ciedge->distance.mean > maxDistance);

#ifdef DEBUG_CIEDGES
    if (ciedge->flags.bits.isProbablyBogus) {
      CIFragT        *fragA = GetCIFragT(ScaffoldGraph->CIFrags, ciedge->fragA);
      CIFragT        *fragB = GetCIFragT(ScaffoldGraph->CIFrags, ciedge->fragB);
      ChunkInstanceT *CIA   = GetGraphNode(graph, ciedge->idA);
      ChunkInstanceT *CIB   = GetGraphNode(graph, ciedge->idB);

      fprintf(stderr,"* Edge ("F_CID","F_CID") %c distance:%f ["F_S32","F_S32"] marked probably bogus --> IS REAL?\nrawDistance "F_S32" ["F_S32","F_S32"] ["F_S32","F_S32"] simChunks ["F_S32","F_S32"] ["F_S32","F_S32"]!\n",
              ciedge->idA, ciedge->idB, ciedge->orient,
              ciedge->distance.mean, minDistance, maxDistance,
              distance,
              fragA->aEndCoord, fragA->bEndCoord,
              fragB->aEndCoord, fragB->bEndCoord,
              CIA->aEndCoord, CIA->bEndCoord,
              CIB->aEndCoord, CIB->bEndCoord );
    }
#endif
  }

  if(insert)
    InsertGraphEdge(graph, ciedgeIndex);

  //  isPotentialRock is set in ComputeMatePairStatisticsRestricted()

  if(collectOverlap &&
     CIa->flags.bits.isPotentialRock &&
     CIb->flags.bits.isPotentialRock){

    if(isOverlapEdge(ciedge)){ // Add an entry to the overlapper's table
      CollectChunkOverlap(graph,
                          ciedge->idA,
                          ciedge->idB,
                          ciedge->orient,
                          -ciedge->distance.mean,
                          (double) fudgeDistance,
                          ciedge->quality,
                          TRUE,
                          TRUE,
                          FALSE);
    }else{
      CollectChunkOverlap(graph,
                          ciedge->idA,
                          ciedge->idB,
                          ciedge->orient,
                          -ciedge->distance.mean,
                          sqrt(distance.variance),
                          ciedge->quality, // bad quality
                          FALSE, // not bayesian
                          FALSE,
                          FALSE);
    }
  }

  return ciedgeIndex;
}




void DumpGraph(GraphCGW_T *graph, FILE *stream){

  char *graphType = "";
  GraphNodeIterator nodes;
  NodeCGW_T *node = NULL;

  if(graph->type == CI_GRAPH)
    graphType = "CI";
  else if (graph->type == CONTIG_GRAPH)
    graphType = "Contig";
  else if (graph->type == SCAFFOLD_GRAPH)
    graphType = "Scaffold";

  fprintf(stream,"* Dumping %s Graph nodes:%d edges:%d\n",
          graphType, (int) GetNumGraphNodes(graph),
          (int) GetNumGraphEdges(graph));

  InitGraphNodeIterator(&nodes, graph, GRAPH_NODE_DEFAULT);
  while(NULL != (node = NextGraphNodeIterator(&nodes))){
    GraphEdgeIterator edges(graph, node->id, ALL_END, ALL_EDGES);
    EdgeCGW_T        *edge = NULL;

    fprintf(stream,"Node "F_CID"\n", node->id);

    while(NULL != (edge = edges.nextMerged())){
      fprintf(stream,"\t("F_CID","F_CID") %c wt:%d %s\n",
              edge->idA, edge->idB, edge->orient.toLetter(), edge->edgesContributing,
              (isOverlapEdge(edge)?"*O*":" "));
    }
  }
}




/* Find an overlap edge between A,B */
EdgeCGW_T *FindGraphOverlapEdge(GraphCGW_T *graph,
                                CDS_CID_t idA, CDS_CID_t idB,
                                PairOrient orient){

  int end;
  int found = FALSE;

  assert(orient.isUnknown() == false);

  if (orient.isAB_AB() || orient.isAB_BA())
    end = B_END;
  else
    end = A_END;

  GraphEdgeIterator edges(graph, idA, end, ALL_EDGES);
  EdgeCGW_T        *edge = NULL;
  
  while(NULL != (edge = edges.nextMerged())) {
    PairOrient corient = GetEdgeOrientationWRT(edge, idA);
    if(!isOverlapEdge(edge)){
      continue;
    }
    if(corient != orient){
      continue;
    }
    if(edge->idA != idB &&
       edge->idB != idB){
      continue;
    }
    found = TRUE;
    break;
  }

  if(found)
    return edge;

  return NULL;
}

/* Find an overlap edge between A,B */
EdgeCGW_T *FindGraphEdge(GraphCGW_T *graph,
                         CDS_CID_t idA, CDS_CID_t idB,
                         PairOrient orient){

  int end = ALL_END;
  int found = FALSE;

  //  ALL_END if orient is Unknown

  if (orient.isAB_AB() || orient.isAB_BA())
    end = B_END;
  if (orient.isBA_AB() || orient.isBA_BA())
    end = A_END;

  GraphEdgeIterator edges(graph, idA, end, ALL_EDGES);
  EdgeCGW_T        *edge;

  while(NULL != (edge = edges.nextMerged())){
    PairOrient corient = GetEdgeOrientationWRT(edge, idA);
    if(corient != orient)
      continue;
    if(edge->idA != idB &&
       edge->idB != idB)
      continue;
    found = TRUE;
    break;
  }

  if(found)
    return edge;

  return NULL;

}

// Find an overlap edge (assumed to exist) between the two CIs.
// Unlink it from the graph
void  DeleteGraphOverlapEdge(GraphCGW_T *graph,
                             CDS_CID_t idA, CDS_CID_t idB,
                             PairOrient orient){
  DeleteGraphEdge(graph, FindGraphOverlapEdge(graph, idA,idB, orient));
}




bool
MergeAllGraphEdges_EdgeCompare(EdgeCGW_T *A, EdgeCGW_T *B) {

#if 1
  if (A->idA < B->idA)   return(true);
  if (A->idA > B->idA)   return(false);

  if (A->idB < B->idB)   return(true);
  if (A->idB > B->idB)   return(false);

  if (A->orient.toLetter() < B->orient.toLetter())  return(true);
  if (A->orient.toLetter() > B->orient.toLetter())  return(false);


#else
  if (((A->idA <  B->idA)) ||
      ((A->idA == B->idA) && (A->idB <  B->idB)) ||
      ((A->idA == B->idA) && (A->idB == B->idB) && (A->orient.toLetter() < B->orient.toLetter())))
    return(true);
#endif

  return(false);
}

bool
MergeAllGraphEdges_UTG_EdgeIDCompare(CDS_CID_t const idA, CDS_CID_t const idB) {
  return(edgeCompareForMerging(GetGraphEdge(ScaffoldGraph->CIGraph, idA),
                               GetGraphEdge(ScaffoldGraph->CIGraph, idB)));
}

bool
MergeAllGraphEdges_CTG_EdgeIDCompare(CDS_CID_t const idA, CDS_CID_t const idB) {
  return(edgeCompareForMerging(GetGraphEdge(ScaffoldGraph->ContigGraph, idA),
                               GetGraphEdge(ScaffoldGraph->ContigGraph, idB)));
}

bool
MergeAllGraphEdges_SCF_EdgeIDCompare(CDS_CID_t const idA, CDS_CID_t const idB) {
  return(edgeCompareForMerging(GetGraphEdge(ScaffoldGraph->ScaffoldGraph, idA),
                               GetGraphEdge(ScaffoldGraph->ScaffoldGraph, idB)));
}

void
MergeAllGraphEdges(GraphCGW_T         *graph,
                   vector<CDS_CID_t>  &rawEdges,
                   bool                includeGuides,
                   bool                mergeAll) {
  vector<CDS_CID_t>  newEdges;
  vector<CDS_CID_t>  chain;

  uint32             minReportSize = 25000;
  uint32             maxTestSize   = 46340;  //  Otherwise we overflow an int

  newEdges.reserve(rawEdges.size());  //  Slight overkill.
  chain.reserve(rawEdges.size());     //  Massive overkill.

#ifdef EDGELOG
  static int  EdgeLog_FileNumber             = 0;
  char        EdgeLog_FileName[FILENAME_MAX] = {0};
#endif

  //  Crud, we have a list of edge IDs as input, but we need to sort using the
  //  edge itself.  As we're building rawEdges, the edge storage array is probably
  //  reallocated, and so we cannot store pointers in rawEdges.  We either need
  //  to use three different sort functions (one per graph) or rebuild rawEdges
  //  using pointers.
  //
  if      (graph == ScaffoldGraph->CIGraph) {
    if (rawEdges.size() > minReportSize)
      fprintf(stderr, "MergeAllGraphEdges()--  Working on unitig edges; includeGuides=%c mergeAll=%c.\n",
              (includeGuides) ? 'T' : 'F', (mergeAll) ? 'T' : 'F');
    sort(rawEdges.begin(), rawEdges.end(), MergeAllGraphEdges_UTG_EdgeIDCompare);

#ifdef EDGELOG
    sprintf(EdgeLog_FileName, "%s.unitigEdges.%d", ScaffoldGraph->name, ++EdgeLog_FileNumber);
#endif

  } else if (graph == ScaffoldGraph->ContigGraph) {
    if (rawEdges.size() > minReportSize)
      fprintf(stderr, "MergeAllGraphEdges()--  Working on contig edges; includeGuides=%c mergeAll=%c.\n",
              (includeGuides) ? 'T' : 'F', (mergeAll) ? 'T' : 'F');
    sort(rawEdges.begin(), rawEdges.end(), MergeAllGraphEdges_CTG_EdgeIDCompare);

  } else if (graph == ScaffoldGraph->ScaffoldGraph) {
    if (rawEdges.size() > minReportSize)
      fprintf(stderr, "MergeAllGraphEdges()--  Working on scaffold edges; includeGuides=%c mergeAll=%c.\n",
              (includeGuides) ? 'T' : 'F', (mergeAll) ? 'T' : 'F');
    sort(rawEdges.begin(), rawEdges.end(), MergeAllGraphEdges_SCF_EdgeIDCompare);

#ifdef EDGELOG
    sprintf(EdgeLog_FileName, "%s.scaffoldEdges.%d", ScaffoldGraph->name, ++EdgeLog_FileNumber);
#endif

  } else {
    fprintf(stderr, "MergeAllGraphEdges()--  Invalid edges.\n");
    assert(0);
  }


  //  Walk down the now sorted list of edges, merging anything with the same
  //  (idA,idB,orient) triplet.

  for (uint32 rid=0; rid < rawEdges.size(); ) {
    EdgeCGW_T  *A = GetGraphEdge(graph, rawEdges[rid]);

    chain.clear();
    chain.push_back(rawEdges[rid++]);

    EdgeCGW_T  *B = (rid < rawEdges.size()) ? GetGraphEdge(graph, rawEdges[rid]) : NULL;

    //  Keep adding to the merge set if the triplet agrees, or we can include guides or the edge is sloppy.

    while ((B != NULL) &&
           (A->idA == B->idA) &&
           (A->idB == B->idB) &&
           (A->orient.toLetter() == B->orient.toLetter()) &&
           ((includeGuides == TRUE) || (isSloppyEdge(B) == FALSE))) {
      chain.push_back(rawEdges[rid++]);

      //(!includeGuides && isSloppyEdge(edge))  /* don't merge guides */ ){

      B = (rid < rawEdges.size()) ? GetGraphEdge(graph, rawEdges[rid]) : NULL;
    }

    if (chain.size() > minReportSize)
      fprintf(stderr, "MergeAllGraphEdges()--  processing chain of size "F_SIZE_T" for triplet "F_U32" "F_U32" %c.\n",
              chain.size(), A->idA, A->idB, A->orient.toLetter());

    if (chain.size() < maxTestSize)
      MergeGraphEdges(graph, chain);
    else
      fprintf(stderr, "MergeAllGraphEdges()--  chain of size "F_SIZE_T" too large to merge.\n",
              chain.size());

    //newEdges.insert(newEdges.begin(), chain.begin(), chain.end());
    for (uint32 cc=0; cc<chain.size(); cc++)
      newEdges.push_back(chain[cc]);
  }

  //  Now add the new edges to the graph.

  if (rawEdges.size() > minReportSize)
    fprintf(stderr, "MergeAllGraphEdges()--  processed "F_SIZE_T" raw edges into "F_SIZE_T" merged edges.\n",
            rawEdges.size(), newEdges.size());

  for (uint32 ee=0; ee<newEdges.size(); ee++)
    InsertGraphEdge(graph, newEdges[ee]);

#ifdef EDGELOG
  if (EdgeLog_FileName[0] != 0) {
    FILE *fp = fopen(EdgeLog_FileName, "w");

    for (uint32 ee=0; ee<newEdges.size(); ee++) {
      EdgeCGW_T *edge = GetGraphEdge(graph, newEdges[ee]);

      fprintf(fp, "EDGE "F_CID" "F_CID" %c frag "F_CID" "F_CID" weight %d overlap %d dist %.0f +- %.0f lib "F_CID"\n",
              edge->idA,
              edge->idB,
              edge->orient.toLetter(),
              edge->fragA,
              edge->fragB,
              edge->edgesContributing,
              isOverlapEdge(edge),
              edge->distance.mean, sqrt(edge->distance.variance),
              edge->distIndex);

      while (edge->nextRawEdge != NULLINDEX) {
        edge = GetGraphEdge(graph, edge->nextRawEdge);

        fprintf(fp, "     "F_CID" "F_CID" %c frag "F_CID" "F_CID" weight %d overlap %d dist %.0f +- %.0f lib "F_CID" RAW\n",
                edge->idA,
                edge->idB,
                edge->orient.toLetter(),
                edge->fragA,
                edge->fragB,
                edge->edgesContributing,
                isOverlapEdge(edge),
                edge->distance.mean, sqrt(edge->distance.variance),
                edge->distIndex);
      }
    }

    fclose(fp);
  }
#endif
}





/* CheckEdgesAgainstOverlapper
   ALL CIEdges that have a contributing overlap should have an overlap entry in the
   ScaffoldGraph overlap hashtable, and the length of the overlap should be > 0.

   This routine checks for this invariant, and asserts if it is not satisfied.

*/
void CheckEdgesAgainstOverlapper(GraphCGW_T *graph){
  EdgeCGW_T *edge;
  GraphNodeIterator nodes;
  NodeCGW_T *node;
  ChunkOverlapCheckT olap;
  int count    = 0;
  int failures = 0;

  fprintf(stderr,"**** Calling CheckEdgesAgainstOverlapper ****\n");
  InitGraphNodeIterator(&nodes, graph, GRAPH_NODE_DEFAULT);

  while(NULL != (node = NextGraphNodeIterator(&nodes))){
    GraphEdgeIterator edgeMates(graph, node->id, ALL_END, ALL_EDGES);

    while(NULL != (edge = edgeMates.nextMerged())){
      if(isOverlapEdge(edge)){
        count++;
        if (LookupOverlap(graph,
                          edge->idA, edge->idB,
                          edge->orient, &olap) == FALSE) {
          fprintf(stderr,"* Failure in CheckEdgeAgainstOverlapper for edge #%d:\n", count);
          PrintGraphEdge(stderr, graph, "F", edge, node->id);
          failures++;
        }
      }
    }
  }

  //  assert(failures == 0);

  fprintf(stderr,"**** Survived CheckEdgesAgainstOverlapper with %d failures****\n", failures);
}


/* Update all the Unitigs belonging to the multiAlignment for contig
   contigID so that their membership and offsets are recorded properly. */

static VA_TYPE(int32) *UngappedOffsets = NULL;

void UpdateNodeUnitigs(MultiAlignT *ma, ContigT *contig){
  int32 i;
  int32 numCIs = GetNumIntUnitigPoss(ma->u_list);
  NodeCGW_T *previous = NULL;
  int32 *offsets;

  contig->flags.bits.isChaff = FALSE;  // We need this one in the output

  if(!UngappedOffsets){
    UngappedOffsets = CreateVA_int32(1000);
  }

  GetMultiAlignUngappedOffsets(ma, UngappedOffsets);
  offsets = Getint32(UngappedOffsets,0);

  assert(numCIs > 0);

  for (int32 i=0; i<numCIs ; i++){
    IntUnitigPos *pos = GetIntUnitigPos(ma->u_list, i);
    NodeCGW_T *node = GetGraphNode(ScaffoldGraph->CIGraph, pos->ident);
    int flip = (pos->position.end < pos->position.bgn);
    int32 bgn, end;

    // mp->position is an interval.  We used to think we needed to subtract one from
    // the upper end of the interval
    if(flip){
      bgn = pos->position.bgn /* - 1*/;
      end = pos->position.end;
    }else{
      bgn = pos->position.bgn;
      end = pos->position.end /* - 1 */;
    }
    // Set the contigID
    node->info.CI.contigID = contig->id;

    // Set the scaffoldID
    node->scaffoldID = contig->scaffoldID;
    node->flags.bits.isChaff = FALSE;  // We need this one in the output

    // Set the intra-contig position
    node->offsetAEnd.mean = offsets[bgn];
    node->offsetAEnd.variance = ComputeFudgeVariance(node->offsetAEnd.mean);

    node->offsetBEnd.mean = offsets[end];
    node->offsetBEnd.variance = ComputeFudgeVariance(node->offsetBEnd.mean);


    // Link to neighbors
    if(previous != NULL){
      previous->BEndNext = pos->ident;
      node->AEndNext = previous->id;
    }else{
      contig->info.Contig.AEndCI = pos->ident;
      node->AEndNext = NULLINDEX;
    }
    previous = node;
  }
  previous->BEndNext = NULLINDEX;
  contig->info.Contig.BEndCI = previous->id;

}


/* Update all the fragments belonging to the multiAlignment for Chunk cid
   so that their membership and offsets in their CIFragT record are recorded
   properly.
   Call this on either a CI or a Contig after it is created */
void UpdateNodeFragments(GraphCGW_T *graph, CDS_CID_t cid,
                         int markFragmentsPlaced, int markUnitigAndContig){
  MultiAlignT *ma;
  static VA_TYPE(int32) *ungappedOffsets = NULL;
  NodeCGW_T *node = GetGraphNode(graph, cid);
  int isPlaced = markFragmentsPlaced || (node->scaffoldID != NULLINDEX);
  CDS_CID_t i;
  CDS_CID_t extremalA = NULLINDEX;
  CDS_CID_t extremalB = NULLINDEX;
  int32 minOffset = INT32_MAX;
  int32 maxOffset = INT32_MIN;

  if(ungappedOffsets == NULL)
    ungappedOffsets = CreateVA_int32(10);

  ma = ScaffoldGraph->tigStore->loadMultiAlign(cid, graph->type == CI_GRAPH);
  assert(ma != NULL);

  GetMultiAlignUngappedOffsets(ma, ungappedOffsets);

  /* Determine extremal fragments so we can label the fragments */
  for(i = 0; i < GetNumIntMultiPoss(ma->f_list); i++){
    IntMultiPos *mp = GetIntMultiPos(ma->f_list, i);
    int32 end = MAX( mp->position.end, mp->position.bgn);
    int32 beg = MIN( mp->position.end, mp->position.bgn);

    if(minOffset > beg){
      minOffset = beg;
      extremalA = i;
    }
    if(maxOffset < end){
      maxOffset = end;
      extremalB = i;
    }

  }

  for(i = 0; i < GetNumIntMultiPoss(ma->f_list); i++){
    IntMultiPos *mp = GetIntMultiPos(ma->f_list, i);
    CIFragT *frag = GetCIFragT(ScaffoldGraph->CIFrags, mp->ident);
    LengthT offset3p, offset5p;
    int flip = (mp->position.end < mp->position.bgn);

    //  This baloney is to adjust unitig fragment positions.  They should be asserts, but we'll just
    //  fix up positions so they are within bounds.

    int32 utgLen = GetNumint32s(ungappedOffsets) - 1;
    int32 bgn    = mp->position.bgn;
    int32 end    = mp->position.end;

    if (mp->position.bgn < mp->position.end) {
      int32  origLen = mp->position.end - mp->position.bgn;

      if (mp->position.bgn < 0) {
        bgn = 0;
        end = origLen;
      }
      if (utgLen < end) {
        bgn = utgLen - origLen;
        end = utgLen;
      }
      if (bgn < 0)
        bgn = 0;

    } else {
      int32  origLen = mp->position.bgn - mp->position.end;

      if (mp->position.end < 0) {
        end = 0;
        bgn = origLen;
      }
      if (utgLen < bgn) {
        end = utgLen - origLen;
        bgn = utgLen;
      }
      if (end < 0)
        end = 0;
    }

    assert(0 <= bgn);
    assert(bgn < ungappedOffsets->numElements);

    assert(0 <= end);
    assert(end < ungappedOffsets->numElements);

    int32 ubgn = *Getint32(ungappedOffsets, bgn);
    int32 uend = *Getint32(ungappedOffsets, end);

    //  if we fall entirely within a gap, give it at least a positive size.

    if (ubgn == uend) {
      if (mp->position.bgn < mp->position.end) {
        if (ubgn == 0)
          uend++;
        else
          ubgn--;
      } else {
        if (uend == 0)
          ubgn++;
        else
          uend--;
      }
    }

    offset5p.mean     = ubgn;
    offset5p.variance = ComputeFudgeVariance(offset5p.mean);

    offset3p.mean     = uend;
    offset3p.variance = ComputeFudgeVariance(offset3p.mean);

    frag->flags.bits.isPlaced = isPlaced;

    // On initial call, each unitig corresponds to a contig,
    // so we mark them together.
    if(node->flags.bits.isContig || markUnitigAndContig){
      frag->contigID = node->id;
      frag->contigOffset3p = offset3p;
      frag->contigOffset5p = offset5p;
    }
    if(node->flags.bits.isCI || markUnitigAndContig){
      frag->CIid = frag->cid = node->id;
      frag->offset3p = offset3p;
      frag->offset5p = offset5p;
    }

    if(i == extremalA){
      frag->flags.bits.chunkLabel = AS_INTERCHUNK_A;
    }else if(i == extremalB){
      frag->flags.bits.chunkLabel = AS_INTERCHUNK_B;  /*  A->B? */
    }else{
      frag->flags.bits.chunkLabel = AS_INTRACHUNK;
    }

    // for placement constraints, mark fragments as constrained
    if (ScaffoldGraph->gkpStore->gkStore_getFRGtoPLC(frag->read_iid) != 0) {
      if (node->scaffoldID == -1) {
        node->flags.bits.isClosure = TRUE;
      }
    }
  }
}





//  Compute the offset and orientation of a fragment in its chunk/contig orientIsOpposite == TRUE
//  Offset is from 5p end of fragment to the end of the chunk in the direction of the 3p end of the fragment.
//
//  orientIsOpposite == FALSE
//  Offset is from 5p end of fragment to the end of the chunk in the direction of the 5p end of the
//  fragment.
//
//  orientIsOpposite -- TRUE if offset should be calculated from 5' towards end of chunk closest to
//  3' end of fragment.  FALSE if offset should be calculated from 5' towards end of chunk closest
//  to 5' end.  See comments below */


void
FragOffsetAndOrientation(CIFragT          *frag,
                         ChunkInstanceT   *chunk,
                         LengthT          *chunkOffset,    // output
                         SequenceOrient   *chunkOrient,    // output
                         int32            *extremal,       // output
                         int32             orientIsOpposite) {
  LengthT offset3p={0,0}, offset5p={0,0};
  assert(frag && chunk && chunkOffset && chunkOrient && extremal);

  if(chunk->flags.bits.isCI){
    offset3p = frag->offset3p;
    offset5p = frag->offset5p;
    *chunkOrient = getCIFragOrient(frag);

  }else if(chunk->flags.bits.isContig){
    offset3p = frag->contigOffset3p;
    offset5p = frag->contigOffset5p;
    *chunkOrient = GetContigFragOrient(frag);
  }else assert(0);

  *extremal = FALSE;

  assert(chunkOrient->isUnknown() == false);

  if(orientIsOpposite){   /* This is the case for ALL Mates and Bac Ends */
    if (chunkOrient->isForward()) {
      /*
        Chunk                        ----------------------------------->
        Frag                           --------->
        ChunkOffset                    |--------------------------------|
        ChunkOrient                         A_B
      */
      chunkOffset->mean = chunk->bpLength.mean - offset5p.mean;
      if(frag->flags.bits.chunkLabel == AS_INTERCHUNK_B || frag->flags.bits.chunkLabel == AS_SINGLETON)
        *extremal = TRUE;
    } else {
      /*
        Chunk                        ----------------------------------->
        Frag                           <---------
        ChunkOffset                  |----------|
        ChunkOrient                     B_A
      */
      chunkOffset->mean = offset5p.mean;
      if(frag->flags.bits.chunkLabel == AS_INTERCHUNK_A || frag->flags.bits.chunkLabel == AS_SINGLETON)
        *extremal = TRUE;
    }
  }else{
    if (chunkOrient->isForward()) {
      /*
        Chunk                        ----------------------------------->
        Frag                           --------->
        ChunkOffset                  |-|
        ChunkOrient                         A_B
      */
      chunkOffset->mean = offset5p.mean;
      if(frag->flags.bits.chunkLabel == AS_INTERCHUNK_B || frag->flags.bits.chunkLabel == AS_SINGLETON)
        *extremal = TRUE;
    } else {
      /*
        Chunk                        ----------------------------------->
        Frag                           <---------
        ChunkOffset                              |----------------------|
        ChunkOrient                        B_A
      */
      chunkOffset->mean = chunk->bpLength.mean - offset5p.mean;
      if(frag->flags.bits.chunkLabel == AS_INTERCHUNK_A || frag->flags.bits.chunkLabel == AS_SINGLETON)
        *extremal = TRUE;
    }
  }

  chunkOffset->variance = ComputeFudgeVariance(chunkOffset->mean);
}


PairOrient
ciEdgeOrientFromFragment(int          orient,
                         SequenceOrient  ciOrient,
                         SequenceOrient  mciOrient) {
  PairOrient  ciEdgeOrient;

  // The following triply nested switch/if statement captures all of the cases that arise from
  // different relative alignments of the fragments in the LKG relationship, and their alignment
  // with their respective chunks.
  //
  //  The original version had 'pictures' which contributed a lot of useless information (like the
  //  distance between mates), but didn't help at clarifying what is going on here...and probably
  //  hurt just because of the volume of noise.  See version 1.83 if you want the pictures.
  //
  //  One picture, for AS_READ_ORIENT_INNIE, ciOrient.isForward(), mciOrient.isForward() (aka, the
  //  first case):
  //
  //           length - 5'             gap            length - 5'
  //      |------------------------||---------------||-----------|
  //  A --------------------------- B               B --------------------------- A
  //    5'----->                                           <------5'
  //      |-------------------------------------------------------|
  //                             mate distance
  //
  //
  //  BPW's best guess here is that ciOrient and mciOrient are the orientation of the fragment
  //  relative to the containing ci.  In this case, we know that the frags are innies, and that they
  //  are both oriented in the same direction as the ci's.  Thus, we know that the ci's are AB_BA.
  //
  //  In the second case, the mate fragment is reversed from the mci, which will reverse the mci
  //  itself, hence AB_AB.

  assert(ciOrient.isUnknown() == false);
  assert(mciOrient.isUnknown() == false);

  switch(orient){
    case AS_READ_ORIENT_INNIE:
      if (ciOrient.isForward()) {
        if (mciOrient.isForward()) {ciEdgeOrient.setIsAB_BA();} else {ciEdgeOrient.setIsAB_AB();}
      } else {
        if (mciOrient.isForward()) {ciEdgeOrient.setIsBA_BA();} else {ciEdgeOrient.setIsBA_AB();}
      }
      break;

    case AS_READ_ORIENT_NORMAL:
      if (ciOrient.isForward()) {
        if (mciOrient.isReverse()) {ciEdgeOrient.setIsAB_BA();} else {ciEdgeOrient.setIsAB_AB();}
      } else {
        if (mciOrient.isForward()) {ciEdgeOrient.setIsBA_AB();} else {ciEdgeOrient.setIsBA_BA();}
      }
      break;

    case AS_READ_ORIENT_ANTINORMAL:
      if (ciOrient.isForward()) {
        if (mciOrient.isReverse()) {ciEdgeOrient.setIsBA_AB();} else {ciEdgeOrient.setIsBA_BA();}
      } else{
        if (mciOrient.isForward()) {ciEdgeOrient.setIsAB_BA();} else {ciEdgeOrient.setIsBA_BA();}
      }
      break;

    case AS_READ_ORIENT_OUTTIE:
      if (ciOrient.isForward()) {
        if (mciOrient.isReverse()) {ciEdgeOrient.setIsBA_BA();} else {ciEdgeOrient.setIsBA_AB();}
      } else {
        if (mciOrient.isReverse()) {ciEdgeOrient.setIsAB_BA();} else {ciEdgeOrient.setIsAB_AB();}
      }
      break;

    default:
      assert(0);
  }

  return(ciEdgeOrient);
}




// Create All raw link-based graph edges
void
BuildGraphEdgesDirectly(GraphCGW_T         *graph,
                        vector<CDS_CID_t>  &rawEdges) {

  GraphEdgeStatT    stat;
  GraphNodeIterator Nodes;
  NodeCGW_T        *node;
  uint32            tf = 0;
  uint32            nf = 0;

  fprintf(stderr,"BuildGraphEdgesDirectly()--\n");

  InitGraphEdgeStatT(&stat);

  InitGraphNodeIterator(&Nodes, graph, GRAPH_NODE_DEFAULT);
  while(NULL != (node = NextGraphNodeIterator(&Nodes))){
    if (node->flags.bits.isChaff && GlobalData->ignoreChaffUnitigs)
      continue;

    assert(node->flags.bits.isDead == 0);

    MultiAlignT  *ma = ScaffoldGraph->tigStore->loadMultiAlign(node->id, graph->type == CI_GRAPH);

    nf += GetNumIntMultiPoss(ma->f_list);
    tf += GetNumIntMultiPoss(ma->f_list);

    if (nf > 10000000) {
      fprintf(stderr, "BuildGraphEdgesDirectly()-- at unitig %d with "F_U32" total fragments.\n",
              node->id, tf);
      nf = 0;
    }

    BuildGraphEdgesFromMultiAlign(graph, node, ma, &stat, FALSE, &rawEdges);
  }

  double f = (stat.totalMatePairs > 0) ? (100.0 * stat.totalExternalMatePairs / stat.totalMatePairs) : (0.0);

  fprintf(stderr,"BuildGraphEdgesDirectly()-- Found %d fragments\n", stat.totalFragments);
  fprintf(stderr,"BuildGraphEdgesDirectly()-- Found %d/%d BacEnd pairs Unique-Unique\n", stat.totalUUBacPairs, stat.totalBacPairs);
  fprintf(stderr,"BuildGraphEdgesDirectly()-- Found %d/%d (%.2f%%) mate pairs node external\n", stat.totalExternalMatePairs, stat.totalMatePairs, f);

  //  Make sure that there are no edges associated with the objects yet.
  InitGraphNodeIterator(&Nodes, graph, GRAPH_NODE_DEFAULT);
  while (NULL != (node = NextGraphNodeIterator(&Nodes)))
    assert(graph->edgeLists[node->id].empty() == true);
}



// Create the raw link-based edges
void
BuildGraphEdgesFromMultiAlign(GraphCGW_T         *graph,
                              NodeCGW_T          *node,
                              MultiAlignT        *ma,
                              GraphEdgeStatT     *stat,
                              int                 buildAll,
                              vector<CDS_CID_t>  *rawEdges) {

  if (stat)
    stat->totalFragments += GetNumIntMultiPoss(ma->f_list);

  for (uint32 i=0; i<GetNumIntMultiPoss(ma->f_list); i++) {
    IntMultiPos *mp     = GetIntMultiPos(ma->f_list, i);
    CDS_CID_t    frgID  = mp->ident;
    CIFragT     *frg    = GetCIFragT(ScaffoldGraph->CIFrags, mp->ident);
    NodeCGW_T   *frgCtg = GetGraphNode(graph, (graph->type == CI_GRAPH) ? frg->cid : frg->contigID);

    if ((frg->flags.bits.hasMate == 0) ||
        (frg->mate_iid == 0))
      //  Not mated.
      continue;

    //  Skip it if the frag was discovered, previously, to have its mate in the same contig or unitig

    if ((frgCtg->flags.bits.isContig) && (frg->flags.bits.hasInternalOnlyContigLinks))
      continue;

    if ((frgCtg->flags.bits.isCI) && (frg->flags.bits.hasInternalOnlyCILinks))
      continue;

    CDS_CID_t    mrgID  = frg->mate_iid;
    CIFragT     *mrg    = GetCIFragT(ScaffoldGraph->CIFrags, mrgID);
    NodeCGW_T   *mrgCtg = GetGraphNode(graph, (graph->type == CI_GRAPH) ? mrg->cid : mrg->contigID);

    assert(frgID != mrgID);

    //  Guard against bogus data.  These should never occur.  If they do, your input unitigs are
    //  missing fragments....and we should have already caught that in Input_CGW.c.  Possibly, we can
    //  simply return FALSE instead of failing, but your input IS messed up, and that's bad.

    if (frg->flags.bits.isDeleted)
      fprintf(stderr, "BuildGraphEdgesFromMultiAlign()-- WARNING: frag %d (mate %d) is DELETED!\n",
              frg->read_iid, frg->mate_iid);

    if (mrg->flags.bits.isDeleted)
      fprintf(stderr, "BuildGraphEdgesFromMultiAlign()-- WARNING: frag %d (mate %d) is DELETED!\n",
              mrg->read_iid, mrg->mate_iid);

    if (frgCtg == NULL)
      fprintf(stderr, "BuildGraphEdgesFromMultiAlign()-- WARNING: node is NULL!  isCI=%d  node %d/%d/%d mnode %d/%d/%d\n",
              (graph->type == CI_GRAPH), frg->read_iid, frg->cid, frg->contigID, mrg->read_iid, mrg->cid, mrg->contigID);

    if (mrgCtg == NULL)
      fprintf(stderr, "BuildGraphEdgesFromMultiAlign()-- WARNING: mnode is NULL!  isCI=%d  node %d/%d/%d mnode %d/%d/%d\n",
              (graph->type == CI_GRAPH), frg->read_iid, frg->cid, frg->contigID, mrg->read_iid, mrg->cid, mrg->contigID);

    assert(frg->flags.bits.isDeleted == 0);
    assert(mrg->flags.bits.isDeleted == 0);
    assert(frgCtg != NULL);
    assert(mrgCtg != NULL);

    if (stat)
      stat->totalMatePairs++;

    assert((graph->type == CI_GRAPH) || (graph->type == CONTIG_GRAPH));
    assert(graph->type != SCAFFOLD_GRAPH);


    //  Don't add edges to chaff, or mates in the same object
    if ((GlobalData->ignoreChaffUnitigs && (frgCtg->flags.bits.isChaff || mrgCtg->flags.bits.isChaff)) ||
        (frgCtg->id == mrgCtg->id)) {

      if (node->flags.bits.isContig) {
        frg->flags.bits.hasInternalOnlyContigLinks = TRUE;
      }
      
      if (node->flags.bits.isCI) {
        frg->flags.bits.hasInternalOnlyCILinks     = TRUE;
        frg->flags.bits.hasInternalOnlyContigLinks = TRUE;
      }

      continue;
    }

    //  Only build specific orientations.

    if ((buildAll == TRUE) || (frgCtg->id < mrgCtg->id)) {
      LengthT        ciOffset;
      LengthT        miOffset;
      int32          extremalA;
      int32          extremalB;
      SequenceOrient ciOrient;
      SequenceOrient miOrient;

      assert(frg->flags.bits.innieMate == mrg->flags.bits.innieMate);

      FragOffsetAndOrientation(frg, frgCtg, &ciOffset, &ciOrient, &extremalA, frg->flags.bits.innieMate);
      FragOffsetAndOrientation(mrg, mrgCtg, &miOffset, &miOrient, &extremalB, frg->flags.bits.innieMate);

      PairOrient ciEdgeOrient = ciEdgeOrientFromFragment(frg->flags.bits.innieMate, ciOrient, miOrient);

      // Since the two offsets and the dist are independent we SUM their variances

      DistT      *dist = GetDistT(ScaffoldGraph->Dists, frg->dist);
      LengthT     distance;

      distance.mean     = dist->mu - ciOffset.mean - miOffset.mean;
      distance.variance = dist->sigma * dist->sigma + ciOffset.variance + miOffset.variance;

      EdgeStatus status = AS_CGW_SafeConvert_uintToEdgeStatus(frg->flags.bits.edgeStatus);

      // Insert a chunk Overlap in preparation for ComputeOverlaps when we are building the extended
      // Unitig Graph.  With contigs, overlaps will be discovered as needed, and no call to
      // ComputeOverlaps is performed.

      //  Do not insert if rawEdges is supplied.  Add them to the vector of raw edges instead.

      CDS_CID_t cid = AddGraphEdge(graph,
                                   frgCtg->id,
                                   mrgCtg->id,
                                   frgID,
                                   mrgID,
                                   frg->dist,
                                   distance,
                                   1.0,
                                   sqrt(ciOffset.variance + miOffset.variance), // This is used by collectOverlap as the fudge distance
                                   ciEdgeOrient,
                                   FALSE,
                                   FALSE,    // isOverlap
                                   FALSE,    // isAContainsB
                                   FALSE,    // isBContainsA
                                   FALSE,    // isTransChunk
                                   extremalA,
                                   extremalB,
                                   status,
                                   graph->type == CI_GRAPH,
                                   (rawEdges == NULL));

      if (rawEdges)
        rawEdges->push_back(cid);
    }

    //  Even if we didn't build, mark the reads as having external mates.

    if (node->flags.bits.isContig) {
      frg->flags.bits.hasInternalOnlyContigLinks = FALSE;
      mrg->flags.bits.hasInternalOnlyContigLinks = FALSE;
    }

    if (node->flags.bits.isCI) {
      mrg->flags.bits.hasInternalOnlyCILinks     = FALSE;
      frg->flags.bits.hasInternalOnlyCILinks     = FALSE;

      frg->flags.bits.hasInternalOnlyContigLinks = FALSE;
      mrg->flags.bits.hasInternalOnlyContigLinks = FALSE;
    }

    if (stat)
      stat->totalExternalMatePairs++;
  }
}  




/* PropagateRawEdgeStatusToFrag */
void PropagateRawEdgeStatusToFrag(EdgeCGW_T *edge){
  CIFragT *frag1, *frag2;
  EdgeStatus status = GetEdgeStatus(edge);
  assert(edge->flags.bits.isRaw);

  if(isOverlapEdge(edge))
    return;


  frag1 = GetCIFragT(ScaffoldGraph->CIFrags, edge->fragA);
  frag2 = GetCIFragT(ScaffoldGraph->CIFrags, edge->fragB);

  frag1->flags.bits.edgeStatus = status;
  frag2->flags.bits.edgeStatus = status;

#if 0
  fprintf(stderr,"*PropagateRawEdgeStatus ("F_CID","F_CID",%c) to fragments "F_CID" and "F_CID" status:%d (from %d)\n",
          edge->idA, edge->idB, edge->orient,
          edge->fragA, edge->fragB, frag1->flags.bits.edgeStatus, status);
#endif
}

/*  PropagateEdgeStatusToFrag
    Iterate through all raw edges constituting this edge, and mark
    their fragments mate status
*/
void PropagateEdgeStatusToFrag(GraphCGW_T *graph, EdgeCGW_T *edge){
  if(edge->flags.bits.isRaw){
    if(edge->fragA == NULLINDEX||edge->fragB == NULLINDEX){
      PrintGraphEdge(stderr,graph," PESTF1 : ",edge,edge->idA);
      fprintf(stderr," * .... isInferred %d\n",
              edge->flags.bits.isInferred);
      assert(0);
    }

    PropagateRawEdgeStatusToFrag(edge);
    return;
  }
  // Propagate to raw edges that are attached
  while(NULL != (edge = GetGraphEdge(graph, edge->nextRawEdge))){

    if(edge->fragA == NULLINDEX||edge->fragB == NULLINDEX){
      PrintGraphEdge(stderr,graph," PESTF2 : ",edge,edge->idB);
      fprintf(stderr," * .... isInferred %d\n",
              edge->flags.bits.isInferred);
      assert(0);
    }

    PropagateRawEdgeStatusToFrag(edge);
  }
  return;
}




static VA_TYPE(IntMultiPos) *f_list_CI = NULL;
static VA_TYPE(IntMultiPos) *f_list_Contig = NULL;


/*
  Assign a subset of the fragments in an unresolved CI to one of its surrogates.
  The fragments listed are marked for membership (via their CIid field) int he new element
*/
void AssignFragsToResolvedCI(GraphCGW_T *graph,
                             CDS_CID_t fromID, CDS_CID_t toID,
                             VA_TYPE(CDS_CID_t) *fragments){
  int i;
  int32 numFrags = GetNumCDS_CID_ts(fragments);
  NodeCGW_T *fromCI = GetGraphNode(graph, fromID);
  NodeCGW_T *toCI = GetGraphNode(graph, toID);
  ContigT *toContig = GetGraphNode(ScaffoldGraph->ContigGraph, toCI->info.CI.contigID);
  int32 surrogateAOffset = toCI->offsetAEnd.mean;
  int32 surrogateBOffset = toCI->offsetBEnd.mean;
  int flipped = (surrogateAOffset > surrogateBOffset);
  IntMultiPos fragPos;
  assert(fromCI->type == UNRESOLVEDCHUNK_CGW);
  assert(toCI->type == RESOLVEDREPEATCHUNK_CGW);
  assert(toCI->scaffoldID != NULLINDEX);

  if(f_list_CI){
    ResetVA_IntMultiPos(f_list_CI);
    ResetVA_IntMultiPos(f_list_Contig);
  }else{
    f_list_CI = CreateVA_IntMultiPos(GetNumCDS_CID_ts(fragments));
    f_list_Contig = CreateVA_IntMultiPos(GetNumCDS_CID_ts(fragments));
  }

  fragPos.delta = NULL;
  fragPos.delta_length = 0;


  /* Check that the fragments CIid == node->id ! */
  /* Then, assign it to the new node and create a degenerate MultiAlignT */
  for(i = 0; i < numFrags; i++){
    CDS_CID_t fragID = *GetCDS_CID_t(fragments,i);
    CIFragT  *frag   = GetCIFragT(ScaffoldGraph->CIFrags, fragID);

    assert(frag->CIid == fromID);
    assert(frag->cid  == fromID);

    frag->CIid     = toID;          // Assign the fragment to the surrogate
    frag->contigID = toContig->id;  // Assign the fragment to the contig

    fragPos.type         = AS_READ;
    fragPos.position.bgn = frag->offset5p.mean;
    fragPos.position.end = frag->offset3p.mean;

    AppendIntMultiPos(f_list_CI, &fragPos);

    // Now figure out fragment position in target contig
    if(flipped){
      fragPos.position.bgn = surrogateAOffset - frag->offset5p.mean;
      fragPos.position.end = surrogateAOffset - frag->offset5p.mean;
    }else{
      fragPos.position.bgn = surrogateAOffset + frag->offset5p.mean;
      fragPos.position.end = surrogateAOffset + frag->offset5p.mean;
    }
    // We shouldn't need this!
    frag->contigOffset5p.variance = 1.0;
    frag->contigOffset5p.mean = fragPos.position.bgn;
    frag->contigOffset3p.variance = 1.0;
    frag->contigOffset3p.mean = fragPos.position.end;

    AppendIntMultiPos(f_list_Contig, &fragPos);
  }

  /* Do not Rebuild the Mate Edges of the target CI tor reflect the changes in fragment membership */

  /* Do NOT Rebuild the Mate Edges of the target CI's Contig to reflect the changes in fragment membership.
     We will rebuild all mate edges when all fragments have been placed */
}


/* Split an unresolved CI, moving a subset of its fragments to the new node.
   Returns the index of the new node, or NULLINDEX if failure.
   The new CI gets copies of all of the overlap edges of its parent, and inherits some of the
   mate edges, as a function of fragments selected.
   The fragments listed are marked for membership (via their CIid field) int he new element
   An empty fragments array is handled the same as fragments == NULL.
*/
CDS_CID_t SplitUnresolvedCI(GraphCGW_T *graph,
                            CDS_CID_t oldNodeID,
                            VA_TYPE(CDS_CID_t) *fragments){
  NodeCGW_T    *newNode   = CreateNewGraphNode(graph);
  CDS_CID_t     newNodeID = newNode->id;

  NodeCGW_T    *oldNode   = GetGraphNode(graph, oldNodeID);

  MultiAlignT  *oldMA     = ScaffoldGraph->tigStore->loadMultiAlign(oldNodeID, TRUE);
  MultiAlignT  *newMA     = NULL;

  assert(graph->type == CI_GRAPH);

  assert(oldNode->flags.bits.isCI);

  newNode->offsetAEnd              = oldNode->offsetAEnd;
  newNode->offsetBEnd              = oldNode->offsetBEnd;
  newNode->flags                   = oldNode->flags;
  newNode->flags.bits.isSurrogate  = TRUE;
  newNode->info                    = oldNode->info;
  newNode->info.CI.numInstances    = 0;
  //newNode->info.CI.numFragments    = (fragments == NULL) ? 0 : GetNumCDS_CID_ts(fragments);
  newNode->info.CI.baseID          = oldNodeID;

  SetNodeType(newNode, RESOLVEDREPEATCHUNK_CGW);

  assert(oldNode->info.CI.numInstances >= 0);
  assert(oldNode->type == UNRESOLVEDCHUNK_CGW);

  /* Bookkeeping on instances derived from the base */

  if (oldNode->info.CI.numInstances > 2)
    assert(oldNode->info.CI.numInstances == GetNumCDS_CID_ts(oldNode->info.CI.instances.va));

  if(oldNode->info.CI.numInstances == 0){
    if( oldNode->flags.bits.isChaff == TRUE){
      CIFragT *frag = GetCIFragT(ScaffoldGraph->CIFrags,  GetIntMultiPos(oldMA->f_list,0)->ident);
      assert(frag->flags.bits.isSingleton);
      newNode->flags.bits.isChaff = FALSE;  // Making a surrogate causes the parent & surrogate to become !chaff
      oldNode->flags.bits.isChaff = FALSE;
      frag->flags.bits.isChaff = FALSE;
    }
    oldNode->info.CI.instances.in_line.instance1 = newNodeID;

  }else if(oldNode->info.CI.numInstances == 1){
    oldNode->info.CI.instances.in_line.instance2 = newNodeID;

  }else if(oldNode->info.CI.numInstances == 2){
    VA_TYPE(CDS_CID_t) *instances = CreateVA_CDS_CID_t(16);
    AppendCDS_CID_t(instances, &oldNode->info.CI.instances.in_line.instance1);
    AppendCDS_CID_t(instances, &oldNode->info.CI.instances.in_line.instance2);
    AppendCDS_CID_t(instances, &(newNodeID));
    oldNode->info.CI.instances.va = instances;

  }else if(oldNode->info.CI.numInstances > 2){
    AppendCDS_CID_t(oldNode->info.CI.instances.va, &(newNodeID));
  }

  oldNode->info.CI.numInstances++;

  if (oldNode->info.CI.numInstances > 2)
    assert(oldNode->info.CI.numInstances == GetNumCDS_CID_ts(oldNode->info.CI.instances.va));

  //  Create Surrogate MultiAlignT -- immediately after the clone, we
  //  own the multialign, but then we hand it over to the store, and
  //  we no longer own it.
  //
  newMA = CloneSurrogateOfMultiAlignT(oldMA, newNodeID);
  assert(newMA->maID == newNodeID);
  ScaffoldGraph->tigStore->insertMultiAlign(newMA, TRUE, TRUE);

  newNode->bpLength.mean     = GetMultiAlignUngappedLength(newMA);
  newNode->bpLength.mean     = newNode->offsetBEnd.mean - newNode->offsetAEnd.mean;

  newNode->bpLength.variance = ComputeFudgeVariance(newNode->bpLength.mean);
  newNode->bpLength.variance = newNode->offsetBEnd.variance - newNode->offsetAEnd.variance;

  assert(GetMultiAlignLength(newMA) == newNode->bpLength.mean);

  /* Copy all of the overlap edges from the original node to the surrogate */

#if 0
  // There is no need to copy the edges on the surrogate unitig
  // they are not putput and not used for anything
  GraphEdgeIterator edges(graph, oldNodeID, ALL_END, ALL_EDGES);
  EdgeCGW_T        *edge;

  while(edge = edges.nextRaw()){
    if(isOverlapEdge(edge)){
      EdgeCGW_T *newEdge = GetFreeGraphEdge(graph);
      CDS_CID_t eid = GetVAIndex_EdgeCGW_T(graph->edges, newEdge);

      CDS_CID_t otherCID = (edge->idA == oldNodeID? edge->idB:edge->idA);
      *newEdge = *edge;
      newEdge->topLevelEdge = eid;

      // Make it canonical WRT otherCID and oldNodeID

      if(oldNodeID < otherCID){
        newEdge->idA = oldNodeID;
        newEdge->idB = otherCID;
      }else{
        newEdge->idB = oldNodeID;
        newEdge->idA = otherCID;
      }
      InsertGraphEdge(graph, eid);
    }
  }
#endif

  /* Assign Fragments to the surrogate */

  if ((fragments) && (GetNumCDS_CID_ts(fragments) > 0))
    AssignFragsToResolvedCI(graph, oldNodeID, newNodeID, fragments);

  //fprintf(stderr,"SplitUnresolvedCI()--  Cloned surrogate ma of CI "F_CID" (length %d) has length (%d) %g\n",
  //        oldNodeID, (int) GetMultiAlignUngappedLength(oldMA),
  //        (int) GetMultiAlignLength(newMA), newNode->bpLength.mean);

  //fprintf(stderr, "SplitUnresolvedCI()--  Split base CI="F_CID" (now with %d instances) into new CI="F_CID"\n",
  //        oldNodeID, oldNode->info.CI.numInstances, newNodeID);

  return newNodeID;
}





/* Split an unresolved Contig.
   This involves splitting its underling CI, and creating a new Contig containing
   the split CI.
   The new CI gets copies of all of the overlap edges of its parent, and inherits some of the
   mate edges, as a function of fragments selected.
   Returns the index of the new node, or NULLINDEX if failure.
   The fragments listed are marked for membership (via their contigID field) int he new element
   An empty fragments array is handled the same as fragments == NULL.
*/
CDS_CID_t SplitUnresolvedContig(GraphCGW_T         *graph,
                                CDS_CID_t           oldID,
                                VA_TYPE(CDS_CID_t) *fragments,
                                int32               copyAllOverlaps){
  NodeCGW_T     *newNode   = CreateNewGraphNode(graph);
  CDS_CID_t      newID     = newNode->id;

  NodeCGW_T     *oldNode   = GetGraphNode(graph, oldID);
  MultiAlignT   *oldMA     = ScaffoldGraph->tigStore->loadMultiAlign(oldID, FALSE);

  CDS_CID_t      baseCIid  = GetIntUnitigPos(oldMA->u_list,0)->ident;
  NodeCGW_T     *baseCI    = GetGraphNode(ScaffoldGraph->CIGraph, baseCIid);

  CDS_CID_t      newCIid   = 0;
  NodeCGW_T     *newCI     = NULL;

  assert(graph->type == CONTIG_GRAPH);
  assert(oldNode->flags.bits.isContig);
  assert(oldNode->scaffoldID == NULLINDEX);
  assert(GetNumIntUnitigPoss(oldMA->u_list) == 1);

  // Split the base CI, creating a surrogate
  newCIid = SplitUnresolvedCI(ScaffoldGraph->CIGraph, baseCIid, fragments);

  baseCI  = GetGraphNode(ScaffoldGraph->CIGraph, baseCIid);
  newCI   = GetGraphNode(ScaffoldGraph->CIGraph, newCIid);

  assert(newCI->type == RESOLVEDREPEATCHUNK_CGW);

  assert(newNode == GetGraphNode(graph, newID));
  assert(oldNode == GetGraphNode(graph, oldID));

  newCI->flags.bits.isSurrogate        = TRUE;
  newCI->flags.bits.isStoneSurrogate   = copyAllOverlaps;
  newCI->flags.bits.isWalkSurrogate    = !copyAllOverlaps;

  newCI->info.CI.contigID              = newID;

  //  What this does: We load unitig 'newCIid', make a copy of it, and
  //  insert it as contig newID.
  //
  //  The unitig instance is already copied in SplitUnresolvedCI.
  //
  MultiAlignT *newMA = CopyMultiAlignT(NULL, ScaffoldGraph->tigStore->loadMultiAlign(newCIid, TRUE));

  newMA->maID = newID;

  ScaffoldGraph->tigStore->insertMultiAlign(newMA, FALSE, TRUE);

  // Create the new surrogate contig

  newNode->bpLength               = newCI->bpLength;
  newNode->flags                  = oldNode->flags;
  newNode->flags.bits.isSurrogate = TRUE;
  newNode->info.Contig.AEndCI     = newCIid;
  newNode->info.Contig.BEndCI     = newCIid;
  newNode->info.Contig.numCI      = 1;

  SetNodeType(newNode, CONTIG_CGW);

  // length of new CI may differ from length of original due to gapped vs ungapped
  {
    int32 newLength = GetMultiAlignUngappedLength(ScaffoldGraph->tigStore->loadMultiAlign(newCIid, TRUE));

    assert(newNode->bpLength.mean == newLength);
    assert(newCI->bpLength.mean   == newLength);
  }


  // Set the fragments offset and contig membership
  {
    int numFrags = (fragments == NULL?0:GetNumCDS_CID_ts(fragments));
    int i;

    for(i = 0; i < numFrags; i++){
      CDS_CID_t fragID = *GetCDS_CID_t(fragments, i);
      CIFragT *frag = GetCIFragT(ScaffoldGraph->CIFrags, fragID);
      frag->contigOffset5p = frag->offset5p;
      frag->contigOffset3p = frag->offset3p;
      frag->contigID = newID;
    }
  }

  /* Copy all of the overlap edges from the original contig to the surrogate */
  if(copyAllOverlaps) {
    EdgeCGW_T *edge;

    // copyOfEdge - Bug fix by Jason Miller, 5/26/06.
    // The call to GetFreeGraphEdge can realloc the graph.
    // In that case, the pointer named edge became invalid.
    // Keeping a copy of the contents solved the problem.

    //    fprintf(stderr,"* Copying edges\n");

    GraphEdgeIterator edges(graph, oldID, ALL_END, ALL_EDGES);
    while(NULL != (edge = edges.nextRaw())){
      if(isOverlapEdge(edge)){
        EdgeCGW_T *newEdge;
        CDS_CID_t eid;
        CDS_CID_t otherCID = (edge->idA == oldID ? edge->idB : edge->idA);

        EdgeCGW_T copyOfEdge = *edge;
        newEdge = GetFreeGraphEdge(graph); // has side effects!
        eid = GetVAIndex_EdgeCGW_T(graph->edges, newEdge);

        *newEdge = copyOfEdge;
        newEdge->topLevelEdge = eid;
        newEdge->nextRawEdge  = NULLINDEX;

        // Make it canonical WRT otherCID and node->id

        if(newID < otherCID){
          newEdge->idA = newID;
          newEdge->idB = otherCID;
          if(oldID > otherCID){
            // If the order of the nodes is different than it used to be, flip
            // the edge orientation.  This shouldn't happen, since  the
            // surrogate ids are higher than any other node
            newEdge->orient.flip();
          }
        }else{
          newEdge->idA = otherCID;
          newEdge->idB = newID;
          if(oldID < otherCID){
            // If the order of the nodes is different than it used to be, flip
            // the edge orientation.  This should occur with probability 0.5
            newEdge->orient.flip();
          }
        }
        // insert the graph edge in the graph
        InsertGraphEdge(graph, eid);

        CreateChunkOverlapFromEdge(graph, newEdge); // add a hashtable entry
      }
    }
  }

  fprintf(stderr,"SplitUnresolvedContig()-- Split base CI="F_CID" (now with %d refs) contig="F_CID" into new CI="F_CID" contig="F_CID"\n",
          baseCIid, baseCI->info.CI.numInstances, oldID, newCIid, newID);

  //DumpContig(stderr,ScaffoldGraph, newNode, FALSE);

  return newID;
}


/***** existsContainmentRelationship *****/

UnitigOverlapType existsContainmentRelationship(NodeCGW_T *ci,
                                                NodeCGW_T *otherCI){
  int32 overlap = IntervalsOverlap(otherCI->offsetAEnd.mean,
                                   otherCI->offsetBEnd.mean,
                                   ci->offsetAEnd.mean,
                                   ci->offsetBEnd.mean,
                                   -15000);

  if(overlap <= 0) {
    return AS_NO_OVERLAP; }
  else if(overlap >= ci->bpLength.mean) {
    return AS_1_CONTAINS_2_OVERLAP; }
  else if(overlap >= otherCI->bpLength.mean) {
    return AS_2_CONTAINS_1_OVERLAP; }
  //else
  return AS_OVERLAP;

}

int compDists( const void *s1, const void *s2)
{
  const MateInfoT * t1 = (const MateInfoT *) s1;
  const MateInfoT * t2 = (const MateInfoT *) s2;
  assert( t1 == s1 );
  assert( t2 == s2 );

  if (t1->samples < t2->samples)
    return -1;
  else if (t1->samples > t2->samples)
    return 1;
  else
    return 0;
}




/* NOTE: This routine should make a collection of distances for each
   dist so that the bucketizing can run on the precomputed collection,
   rather than recomputing.

   Also, we need to iterate to deal with the case that the original
   estimate is bad.
   Low priority.
*/

//#include "obsolete/computematepairstatistics"



int compareInt (const void * a, const void * b) {
  return ( *(int*)a - *(int*)b );
}



void ComputeMatePairDetailedStatus(void) {

  GraphCGW_T *graph = ScaffoldGraph->CIGraph;
  GraphNodeIterator nodes;
  NodeCGW_T *node;
  DistT *dptr;

  int i;

  int numTotalFrags = 0;
  int numReverse    = 0;
  int numGood       = 0;
  int numShort      = 0;
  int numLong       = 0;
  int numSame       = 0;
  int numOuttie     = 0;
  int numNoMate     = 0;
  int numBothChaff  = 0;
  int numChaff      = 0;
  int numBothDegen  = 0;
  int numDegen      = 0;
  int numBothSurr   = 0;
  int numSurrogate  = 0;
  int numDiffScaf   = 0;
  int numInUnresolv = 0;
  int numInResRep   = 0;
  int numInNullScaf = 0;
  int numZero       = 0;
  int numOne        = 0;
  int numMore       = 0;
  int numSkipChaff  = 0;
  int numSkipSurr   = 0;
  int numSkipDegen  = 0;
  int numDeadCtg    = 0;

  int duCI = 0;
  int urCI = 0;
  int uCI  = 0;
  int rrCI = 0;
  int ctig  = 0;
  int uctig = 0;
  int rctig = 0;
  int urctig= 0;
  int rScaf = 0;
  int oScaf = 0;
  int sScaf = 0;

  int fduCI = 0;
  int furCI = 0;
  int fuCI  = 0;
  int frrCI = 0;
  int fctig  = 0;
  int fuctig = 0;
  int frctig = 0;
  int furctig= 0;
  int frScaf = 0;
  int foScaf = 0;
  int fsScaf = 0;

  HashTable_AS *surrHash = CreateScalarHashTable_AS();

  InitGraphNodeIterator(&nodes, graph, GRAPH_NODE_DEFAULT);

  while(NULL != (node = NextGraphNodeIterator(&nodes))) {
    switch (node->type) {
      case DISCRIMINATORUNIQUECHUNK_CGW: duCI++; break;
      case UNRESOLVEDCHUNK_CGW:          urCI++; break;
      case UNIQUECHUNK_CGW:               uCI++; break;
      case RESOLVEDREPEATCHUNK_CGW:      rrCI++; break;
      case CONTIG_CGW:                   ctig++; break;
      case UNIQUECONTIG_CGW:            uctig++; break;
      case RESOLVEDCONTIG_CGW:          rctig++; break;
      case UNRESOLVEDCONTIG_CGW:       urctig++; break;
      case REAL_SCAFFOLD:               rScaf++; break;
      case OUTPUT_SCAFFOLD:             oScaf++; break;
      case SCRATCH_SCAFFOLD:            sScaf++; break;
      default:
        assert(0);
    }

    if (node->flags.bits.isDead)
      continue;

    MultiAlignT *ma = ScaffoldGraph->tigStore->loadMultiAlign(node->id, graph->type == CI_GRAPH);

    int numFrags  = GetNumIntMultiPoss(ma->f_list);

    numTotalFrags += numFrags;

    for( i = 0; i < numFrags; i++) {
      IntMultiPos *mp = GetIntMultiPos(ma->f_list, i);
      CIFragT *frag, *mate;
      int32 dist;

      frag = GetCIFragT(ScaffoldGraph->CIFrags, mp->ident);
      assert(frag->read_iid == mp->ident);
      if (frag->flags.bits.hasMate == 0) {
        numNoMate++;
        frag->flags.bits.mateDetail = NO_MATE;
        continue;
      }
      mate = GetCIFragT(ScaffoldGraph->CIFrags,frag->mate_iid);
      if (mate == NULL) {
        numNoMate++;
        frag->flags.bits.mateDetail = NO_MATE;
        continue;
      }
      if (mate->mate_iid != mp->ident) {
        fprintf (stderr, "ERROR: The gkpStore appears corrupt!\n");
        fprintf (stderr, "INFO: Expect (mate(mate(read))==read.\n");
        fprintf (stderr, "INFO: Alignment %d of %d references read IID %d\n", i, numFrags, mp->ident);
        fprintf (stderr, "INFO: Retrieved read IID %d\n", frag->read_iid);
        fprintf (stderr, "INFO: That read references mate %d\n", frag->mate_iid);
        fprintf (stderr, "INFO: Retrieved read IID %d\n", mate->read_iid);
        fprintf (stderr, "INFO: That read references mate  %d\n", mate->mate_iid);
        assert(mate->mate_iid == mp->ident);
      }
      if(frag->flags.bits.isChaff) {
        if (mate->flags.bits.isChaff) {
          if (mate->flags.bits.mateDetail != BOTH_CHAFF_MATE) {
            mate->flags.bits.mateDetail  = BOTH_CHAFF_MATE;
            frag->flags.bits.mateDetail  = BOTH_CHAFF_MATE;
            numBothChaff+=2;
          }
        } else {
          numChaff+=2;
          mate->flags.bits.mateDetail = CHAFF_MATE;
          frag->flags.bits.mateDetail = CHAFF_MATE;
        }
        continue;
      }
      if(mate->flags.bits.isChaff || frag->flags.bits.mateDetail == CHAFF_MATE) {
        numSkipChaff++;
        continue;
      }

      NodeCGW_T *fragContig, *mateContig;
      fragContig = GetGraphNode( ScaffoldGraph->ContigGraph, frag->contigID);
      AssertPtr(fragContig);

      if (fragContig->flags.bits.isDead) {
        numDeadCtg++;
      }

      switch (fragContig->type) {
        case DISCRIMINATORUNIQUECHUNK_CGW: fduCI++; break;
        case UNRESOLVEDCHUNK_CGW:          furCI++; break;
        case UNIQUECHUNK_CGW:               fuCI++; break;
        case RESOLVEDREPEATCHUNK_CGW:      frrCI++; break;
        case CONTIG_CGW:                   fctig++; break;
        case UNIQUECONTIG_CGW:            fuctig++; break;
        case RESOLVEDCONTIG_CGW:          frctig++; break;
        case UNRESOLVEDCONTIG_CGW:       furctig++; break;
        case REAL_SCAFFOLD:               frScaf++; break;
        case OUTPUT_SCAFFOLD:             foScaf++; break;
        case SCRATCH_SCAFFOLD:            fsScaf++; break;
        default:
          assert(0);
      }

      mateContig = GetGraphNode( ScaffoldGraph->ContigGraph, mate->contigID);
      AssertPtr(mateContig);

      NodeCGW_T *mchunk = GetGraphNode( ScaffoldGraph->CIGraph, mate->CIid);
      if ( node->type == RESOLVEDREPEATCHUNK_CGW)
        numInResRep++;

      if ( node->type != DISCRIMINATORUNIQUECHUNK_CGW ) {
        if (node->info.CI.numInstances == 0) {
          numZero++;
        } else if (node->info.CI.numInstances == 1) {
          numOne++;
        } else {
          numMore++;
        }

        if ( (node->info.CI.numInstances > 1 || node->flags.bits.isStoneSurrogate) && !frag->flags.bits.isPlaced) {
          int fragExists  = ExistsInHashTable_AS(surrHash, frag->read_iid, 0);
          int mateExists  = ExistsInHashTable_AS(surrHash, mate->read_iid, 0);

          if (fragExists && mateExists) {
            if (mate->flags.bits.mateDetail != BOTH_SURR_MATE) {
              numSurrogate-=1;
              numBothSurr+=2;
              mate->flags.bits.mateDetail = BOTH_SURR_MATE;
              frag->flags.bits.mateDetail = BOTH_SURR_MATE;
              //fprintf(stderr, "* Both surr %d,%d from chunks %d.\n",
              //        frag->read_iid, mate->read_iid, node->id, mchunk->id);
            }
            continue;
          } else {
            numSurrogate++;
            mate->flags.bits.mateDetail = SURR_MATE;
            frag->flags.bits.mateDetail = SURR_MATE;
            if (!fragExists) {
              //fprintf(stderr, "* Add frag %d from repeat chunk %d to surr hash, num: %d\n",
              //        frag->read_iid, node->id, node->info.CI.numInstances );
              InsertInHashTable_AS(surrHash, frag->read_iid, 0, 0, 0);
            }
          }
      
          if( (mchunk->info.CI.numInstances > 1 || mchunk->flags.bits.isStoneSurrogate) && !mate->flags.bits.isPlaced ) {
            int fragExists = ExistsInHashTable_AS(surrHash, frag->read_iid, 0);
      
            if (fragExists) { // already seen
              if (frag->flags.bits.mateDetail != BOTH_SURR_MATE) {
                numSurrogate-=1;
                numBothSurr+=2;
                mate->flags.bits.mateDetail = BOTH_SURR_MATE;
                frag->flags.bits.mateDetail = BOTH_SURR_MATE;
                InsertInHashTable_AS(surrHash, mate->read_iid, 0, 0, 0);
              }
              continue;
            }
            int mateExists = ExistsInHashTable_AS(surrHash, mate->read_iid, 0);
            if (!mateExists) {
              //fprintf(stderr, "* Add frag %d from chunk %d to surr hash, num: %d\n",
              //        mate->read_iid, mchunk->id, mchunk->info.CI.numInstances );
              InsertInHashTable_AS(surrHash, mate->read_iid, 0, 0, 0);
              numSurrogate++;
              mate->flags.bits.mateDetail = SURR_MATE;
              frag->flags.bits.mateDetail = SURR_MATE;
            }
          }
          continue;
        }
      }
      
      if ( node->type == UNRESOLVEDCHUNK_CGW)
        numInUnresolv++;
      if ( fragContig->scaffoldID == NULLINDEX )
        numInNullScaf++;

      if( fragContig->scaffoldID == NULLINDEX ) {
        assert( node->type != DISCRIMINATORUNIQUECHUNK_CGW );
        if (mateContig->scaffoldID == NULLINDEX && mchunk->info.CI.numInstances <= 1) {
          if (mate->flags.bits.mateDetail != BOTH_DEGEN_MATE) {
            mate->flags.bits.mateDetail  = BOTH_DEGEN_MATE;
            frag->flags.bits.mateDetail  = BOTH_DEGEN_MATE;
            numBothDegen+=2;
          }
        } else {
          numDegen+=2;
          mate->flags.bits.mateDetail = DEGEN_MATE;
          frag->flags.bits.mateDetail = DEGEN_MATE;
        }
        continue;
      }
      if(mateContig->scaffoldID == NULLINDEX) {
        numSkipDegen++;
        continue;
      }
      if ( fragContig->scaffoldID != mateContig->scaffoldID ) {
        // we want them to be in the same scaffold
        numDiffScaf++;
        mate->flags.bits.mateDetail = DIFF_SCAFF_MATE;
        frag->flags.bits.mateDetail = DIFF_SCAFF_MATE;
        continue;
      }

      int32 fragLeftEnd, fragRightEnd;
      int32 mateLeftEnd, mateRightEnd;
      int mateScaffoldOrientation, fragScaffoldOrientation;

      GetFragmentPositionInScaffold( frag, &fragLeftEnd, &fragRightEnd, &fragScaffoldOrientation);
      GetFragmentPositionInScaffold( mate, &mateLeftEnd, &mateRightEnd, &mateScaffoldOrientation);

      if (fragScaffoldOrientation == mateScaffoldOrientation) {
        numSame++;
        frag->flags.bits.mateDetail = SAME_ORIENT_MATE;
        mate->flags.bits.mateDetail = SAME_ORIENT_MATE;
        continue;
      }

      if (fragScaffoldOrientation == 0) {
        // frag ---->  <---- mate
        if ( fragLeftEnd < mateRightEnd ) {
          //innie
          dist = mateRightEnd - fragLeftEnd;
        } else {
          //outtie if they're psat each other
          numOuttie++;
          frag->flags.bits.mateDetail = OUTTIE_ORIENT_MATE;
          mate->flags.bits.mateDetail = OUTTIE_ORIENT_MATE;
          continue; // real outtie mates not currently supported
        }
      } else {
        // mate ---->  <---- frag
        numReverse++;
        if ( mateLeftEnd < fragRightEnd ) {
          //innie
          dist = fragRightEnd - mateLeftEnd;
        } else {
          //outtie
          numOuttie++;
          frag->flags.bits.mateDetail = OUTTIE_ORIENT_MATE;
          mate->flags.bits.mateDetail = OUTTIE_ORIENT_MATE;
          continue; // real outtie mates not currently supported
        }
      }

      dptr = GetDistT(ScaffoldGraph->Dists, frag->dist);
      if (dist < dptr->lower ) {
        numShort++;
        frag->flags.bits.mateDetail = BAD_SHORT_MATE;
        mate->flags.bits.mateDetail = BAD_SHORT_MATE;
      } else if ( dist > dptr->upper) {
        numLong++;
        frag->flags.bits.mateDetail = BAD_LONG_MATE;
        mate->flags.bits.mateDetail = BAD_LONG_MATE;
      } else {
        numGood++;
        frag->flags.bits.mateDetail = GOOD_MATE;
        mate->flags.bits.mateDetail = GOOD_MATE;
      }
    }
  }

#if 0
  fprintf(stderr,"\n");
  fprintf(stderr,"* Mate counts from ComputeMatePairDetailedStatus()\n");
  fprintf(stderr,"* num Frags                          %d\n",numTotalFrags);
  fprintf(stderr,"* num reverse frags                  %d\n",numReverse);
  fprintf(stderr,"* num frags in dead ctgs             %d\n",numDeadCtg);
  fprintf(stderr,"* num frags in unresolved chunks     %d\n",numInUnresolv);
  fprintf(stderr,"* num frags in repeat chunks         %d\n",numInResRep);
  fprintf(stderr,"* num frags in zero instance chunks  %d\n",numZero);
  fprintf(stderr,"* num frags in one instance chunks   %d\n",numOne);
  fprintf(stderr,"* num frags in >=two instance chunks %d\n",numMore);
  fprintf(stderr,"* num frags in NULL scafs            %d\n",numInNullScaf);
  fprintf(stderr,"* num no mates                       %d\n",numNoMate);
  fprintf(stderr,"* num good mates                     %d\n",numGood);
  fprintf(stderr,"* num bad short mates                %d\n",numShort);
  fprintf(stderr,"* num bad long mates                 %d\n",numLong);
  fprintf(stderr,"* num same orientation mates         %d\n",numSame);
  fprintf(stderr,"* num outtie mates                   %d\n",numOuttie);
  fprintf(stderr,"* num both chaff mates               %d\n",numBothChaff);
  fprintf(stderr,"* num chaff mates                    %d\n",numChaff);
  fprintf(stderr,"* num skiped chaff                   %d\n",numSkipChaff);
  fprintf(stderr,"* num both degen mates               %d\n",numBothDegen);
  fprintf(stderr,"* num degen mates                    %d\n",numDegen);
  fprintf(stderr,"* num skiped degen                   %d\n",numSkipDegen);
  fprintf(stderr,"* num both surrogate mates           %d\n",numBothSurr);
  fprintf(stderr,"* num surrogate mates                %d\n",numSurrogate);
  fprintf(stderr,"* num other scaffold                 %d\n",numDiffScaf);

  int sum = (numNoMate + numGood + numShort + numLong + numSame + numOuttie +
             numBothChaff + numChaff + numBothSurr + numSurrogate + numBothDegen +
             numDegen + numDiffScaf);
  fprintf(stderr,"* sum of frag mate status            %d\n\n",sum);

  //  assert( sum == numTotalFrags );

  fprintf(stderr,"* Counts of top level node type\n");
  fprintf(stderr,"* num DISCRIMINATORUNIQUECHUNK_CGW   %d\n", duCI);
  fprintf(stderr,"* num UNRESOLVEDCHUNK_CGW            %d\n", urCI);
  fprintf(stderr,"* num UNIQUECHUNK_CGW                %d\n", uCI);
  fprintf(stderr,"* num RESOLVEDREPEATCHUNK_CGW        %d\n", rrCI);
  fprintf(stderr,"* num CONTIG_CGW                     %d\n", ctig);
  fprintf(stderr,"* num UNIQUECONTIG_CGW               %d\n", uctig);
  fprintf(stderr,"* num RESOLVEDCONTIG_CGW             %d\n", rctig);
  fprintf(stderr,"* num UNRESOLVEDCONTIG_CGW           %d\n", urctig);
  fprintf(stderr,"* num REAL_SCAFFOLD                  %d\n", rScaf);
  fprintf(stderr,"* num OUTPUT_SCAFFOLD                %d\n", oScaf);
  fprintf(stderr,"* num SCRATCH_SCAFFOLD               %d\n", sScaf);
  fprintf(stderr,"\n");
  fprintf(stderr,"* Counts of contig level node type\n");
  fprintf(stderr,"* num DISCRIMINATORUNIQUECHUNK_CGW   %d\n", fduCI);
  fprintf(stderr,"* num UNRESOLVEDCHUNK_CGW            %d\n", furCI);
  fprintf(stderr,"* num UNIQUECHUNK_CGW                %d\n", fuCI);
  fprintf(stderr,"* num RESOLVEDREPEATCHUNK_CGW        %d\n", frrCI);
  fprintf(stderr,"* num CONTIG_CGW                     %d\n", fctig);
  fprintf(stderr,"* num UNIQUECONTIG_CGW               %d\n", fuctig);
  fprintf(stderr,"* num RESOLVEDCONTIG_CGW             %d\n", frctig);
  fprintf(stderr,"* num UNRESOLVEDCONTIG_CGW           %d\n", furctig);
  fprintf(stderr,"* num REAL_SCAFFOLD                  %d\n", frScaf);
  fprintf(stderr,"* num OUTPUT_SCAFFOLD                %d\n", foScaf);
  fprintf(stderr,"* num SCRATCH_SCAFFOLD               %d\n", fsScaf);
  fprintf(stderr,"\n");
#endif

  DeleteHashTable_AS(surrHash);
}



void ComputeMatePairStatisticsRestricted(int operateOnNodes,
                                         int32 minSamplesForOverride,
                                         char *instance_label) {
  GraphCGW_T *graph = NULL;
  GraphNodeIterator nodes;
  NodeCGW_T *node;

  int numPotentialRocks = 0;
  int numPotentialStones = 0;

  int NN = GetNumDistTs(ScaffoldGraph->Dists);

  DistT                 *dwork        = new DistT [NN];
  VA_TYPE(int32)       **dworkSamples = new VA_TYPE(int32) *     [NN];
  VA_TYPE(CDS_CID_t)   **dworkFrags   = new VA_TYPE(CDS_CID_t) * [NN];
  VA_TYPE(CDS_CID_t)   **dworkMates   = new VA_TYPE(CDS_CID_t) * [NN];

  int i, j;

  AS_UTL_mkdir("stat");

  if (operateOnNodes == UNITIG_OPERATIONS)
    graph = ScaffoldGraph->CIGraph;
  if (operateOnNodes == CONTIG_OPERATIONS)
    graph = ScaffoldGraph->ContigGraph;
  if (operateOnNodes == SCAFFOLD_OPERATIONS)
    graph = ScaffoldGraph->CIGraph;

  //  Copy all distances to our work structure.  Any that have valid
  //  changes will be copied back at the end.  This is needed because
  //  we call ComputeMatePairStatisticsRestricted() twice at the end
  //  of a scaffolding run -- once on scaffolds, then on contigs.  The
  //  contig computation destroys lots of what we just computed for
  //  scaffolds, and if a contig doesn't have enough samples for an
  //  update, we just screwed up the scaffold estimate.

  for(i = 1; i < GetNumDistTs(ScaffoldGraph->Dists); i++) {
    DistT *dorig = GetDistT(ScaffoldGraph->Dists, i);

    dwork[i].mu             = 0.0;
    dwork[i].sigma          = 0.0;
    dwork[i].numSamples     = 0;
    dwork[i].min            = INT32_MAX;
    dwork[i].max            = INT32_MIN;
    dwork[i].bnum           = 0;
    dwork[i].bsize          = 0;
    dwork[i].histogram      = NULL;
    dwork[i].lower          = dorig->mu - CGW_CUTOFF * dorig->sigma;
    dwork[i].upper          = dorig->mu + CGW_CUTOFF * dorig->sigma;
    dwork[i].numBad         = 0;
    dwork[i].allowUpdate    = dorig->allowUpdate;

    dworkSamples[i]         = CreateVA_int32(1024);
    dworkFrags[i]           = CreateVA_CDS_CID_t(1024);
    dworkMates[i]           = CreateVA_CDS_CID_t(1024);
  }

  int numChaff = 0;
  int numSingle = 0;
  int numNolink = 0;
  int numCtgNotInternal = 0;
  int numReverse = 0;
  int numOtherUtg = 0;
  int numSameOrient = 0;
  int numNot5stddev = 0;
  int numDiffScaf = 0;
  int numTotalFrags = 0;
  int numNullMate = 0;
  int numMateNotSource = 0;

  InitGraphNodeIterator(&nodes, graph, GRAPH_NODE_DEFAULT);
  while (NULL != (node = NextGraphNodeIterator(&nodes))) {
    CDS_CID_t i;
    int       numFrags;
    int       numExternalLinks = 0;

    // Don't waste time loading singletons for this
    if(node->flags.bits.isChaff) {
      numChaff++;
      continue;
    }

    MultiAlignT *ma = ScaffoldGraph->tigStore->loadMultiAlign(node->id, graph->type == CI_GRAPH);

    numFrags  = GetNumIntMultiPoss(ma->f_list);
    numTotalFrags += numFrags;

    if (numFrags < 2) {
      numSingle++;
      continue;
    }

    for (i = 0; i < numFrags; i++) {
      IntMultiPos *mp    = NULL;
      CIFragT     *frag  = NULL;
      CIFragT     *mate  = NULL;
      int32        dist  = 0;
      DistT       *dorig = NULL;
      DistT       *dfrg  = NULL;

      mp   = GetIntMultiPos(ma->f_list, i);
      frag = GetCIFragT(ScaffoldGraph->CIFrags, mp->ident);

      assert(frag->read_iid == mp->ident);

      if (frag->flags.bits.hasMate == 0) {
        numNolink++;
        continue;
      }

      if (frag->flags.bits.hasMate == 1 &&  // the typical case
          (operateOnNodes == CONTIG_OPERATIONS && !frag->flags.bits.hasInternalOnlyContigLinks)) {
        //     ||    (!operateOnContigs && !frag->flags.bits.hasInternalOnlyCILinks))
        numCtgNotInternal++;
        continue;
      }

      mate = GetCIFragT(ScaffoldGraph->CIFrags,frag->mate_iid);

      if ((operateOnNodes == UNITIG_OPERATIONS) &&
          (mate != NULL) &&
          (mate->cid != frag->cid))
        numExternalLinks++;

      if (mate == NULL) {
        numNullMate++;
        continue;
      }
      if (mate->mate_iid != mp->ident) {
        numMateNotSource++;
        continue;
      }

      //  dorig -- the original distance, used to get the existing
      //  value of mu and sigma, used below to check that a mate is
      //  the correct size.
      //
      //  dfrg  -- the temporary distance we are computing in
      //
      dorig = GetDistT(ScaffoldGraph->Dists, frag->dist);
      dfrg  = &dwork[frag->dist];

      if (operateOnNodes == UNITIG_OPERATIONS) {
        NodeCGW_T *unitig = GetGraphNode( ScaffoldGraph->CIGraph, frag->cid);

        if (frag->cid != mate->cid) {
          numOtherUtg++;
          continue;
        }

        if (getCIFragOrient(mate) == getCIFragOrient(frag)) {
          frag->flags.bits.edgeStatus = UNTRUSTED_EDGE_STATUS;
          mate->flags.bits.edgeStatus = UNTRUSTED_EDGE_STATUS;
          dfrg->numBad++;
          numSameOrient++;
          continue;
        }

        // now make sure the 5p end is less than the 3p end
        if (frag->offset5p.mean > frag->offset3p.mean) {
          numReverse++;
          continue;
        }

        // try to sample fairly by only doing mates where ones of any length could live
        if (frag->offset5p.mean + dorig->mu + CGW_CUTOFF * dorig->sigma > unitig->bpLength.mean) {
          numNot5stddev++;
          continue;
        }

        dist = mate->offset5p.mean - frag->offset5p.mean;

        if((frag->flags.bits.innieMate && getCIFragOrient(frag).isReverse()) ||
           (!frag->flags.bits.innieMate && getCIFragOrient(frag).isForward()) )
          dist = -dist;
      } else if (operateOnNodes == CONTIG_OPERATIONS) {
        ContigT *contig = GetGraphNode( ScaffoldGraph->ContigGraph, frag->contigID);

        assert(frag->contigID == mate->contigID);
        if(GetContigFragOrient(mate) == GetContigFragOrient(frag)) {
          //  fprintf(stderr,"* ("F_CID","F_CID") is bad due to orientation problems\n",      frag->read_iid, mate->read_iid);
          frag->flags.bits.edgeStatus = UNTRUSTED_EDGE_STATUS;
          mate->flags.bits.edgeStatus = UNTRUSTED_EDGE_STATUS;
          dfrg->numBad++;
          numSameOrient++;
          continue;
        }

        // now make sure the 5p end is less than the 3p end
        if ( frag->contigOffset5p.mean > frag->contigOffset3p.mean) {
          numReverse++;
          continue;
        }

        // try to sample fairly by only doing mates where ones of any length could live
        if ( frag->contigOffset5p.mean + dorig->mu + CGW_CUTOFF * dorig->sigma > contig->bpLength.mean) {
          numNot5stddev++;
          continue;
        }

        dist =  mate->contigOffset5p.mean - frag->contigOffset5p.mean;
        //   -------------------->          <----------------------
        //     mate                 innie            frag

        //   <-------------------           ---------------------->
        //     mate                 outie            frag

        if((frag->flags.bits.innieMate && GetContigFragOrient(frag).isReverse()) ||
           (!frag->flags.bits.innieMate && GetContigFragOrient(frag).isForward()) )
          dist =  -dist;

      } else if (operateOnNodes == SCAFFOLD_OPERATIONS) {
        NodeCGW_T *fragContig, *mateContig;
        int32 fragLeftEnd, fragRightEnd;
        int32 mateLeftEnd, mateRightEnd;
        int mateScaffoldOrientation, fragScaffoldOrientation;

        fragContig = GetGraphNode( ScaffoldGraph->ContigGraph, frag->contigID);
        AssertPtr(fragContig);

        mateContig = GetGraphNode( ScaffoldGraph->ContigGraph, mate->contigID);
        AssertPtr(mateContig);

        // we want them to be in the same scaffold
        if ( fragContig->scaffoldID != mateContig->scaffoldID || fragContig->scaffoldID == -1) {
          numDiffScaf++;
          continue;
        }

        GetFragmentPositionInScaffold( frag, &fragLeftEnd, &fragRightEnd, &fragScaffoldOrientation);
        GetFragmentPositionInScaffold( mate, &mateLeftEnd, &mateRightEnd, &mateScaffoldOrientation);

        if (fragScaffoldOrientation == mateScaffoldOrientation) {
          // fprintf(stderr,"* ("F_CID","F_CID") is bad due to orientation problems\n",      frag->read_iid, mate->read_iid);
          frag->flags.bits.edgeStatus = UNTRUSTED_EDGE_STATUS;
          mate->flags.bits.edgeStatus = UNTRUSTED_EDGE_STATUS;
          dfrg->numBad++;
          numSameOrient++;
          continue;
        }

        // now make sure the 5p end is less than the 3p end
        if ( fragScaffoldOrientation == 1 ) {
          numReverse++;
          continue;
        }

        // try to sample fairly by only doing mates where ones of any length could live
        {
          NodeCGW_T *scaff, *extremeContig;
          int32 contigLeftEnd, contigRightEnd;
          int contigScaffoldOrientation;

          // grab scaffold
          scaff = GetGraphNode( ScaffoldGraph->ScaffoldGraph, fragContig->scaffoldID);
          extremeContig = GetGraphNode( ScaffoldGraph->ContigGraph, scaff->info.Scaffold.BEndCI);
          GetContigPositionInScaffold ( extremeContig, &contigLeftEnd, &contigRightEnd, &contigScaffoldOrientation);

          if ( fragLeftEnd + dorig->mu + CGW_CUTOFF * dorig->sigma > contigRightEnd) {
            numNot5stddev++;
            continue;
          }
        }

        if (frag->flags.bits.innieMate) {
          if (fragScaffoldOrientation == 0) // frag ---->  <---- mate
            dist = mateRightEnd - fragLeftEnd;
          else                              // mate ---->  <---- frag
            dist = fragRightEnd - mateLeftEnd;
        } else { // outtie pair
          if (fragScaffoldOrientation == 0) // mate <----  ----> frag
            dist = fragRightEnd - mateLeftEnd;
          else                              // frag <----  ----> mate
            dist = mateRightEnd - fragLeftEnd;
        }

        //if (dist < 0)
        //  fprintf( stderr, "frag, mate: "F_CID", "F_CID" have negative dist: "F_S32"\n",
        //           frag->read_iid, mate->read_iid, dist);
      }  //  end of operateOnNodes if..elseif..elseif block

      if (dist > 0 && dist < dfrg->min)
        dfrg->min = dist;
      if (dist > dfrg->max)
        dfrg->max = dist;

      Appendint32(dworkSamples[frag->dist], &dist);
      AppendCDS_CID_t(dworkFrags[frag->dist], &frag->read_iid);
      AppendCDS_CID_t(dworkMates[frag->dist], &mate->read_iid);

      // See if the mate distance implied is outside of a 5-sigma range

      if ((dist < dfrg->lower) || (dist > dfrg->upper)) {
        frag->flags.bits.edgeStatus = UNTRUSTED_EDGE_STATUS;
        mate->flags.bits.edgeStatus = UNTRUSTED_EDGE_STATUS;
        dfrg->numBad++;
      } else {
        frag->flags.bits.edgeStatus = TRUSTED_EDGE_STATUS;
        mate->flags.bits.edgeStatus = TRUSTED_EDGE_STATUS;

        //if (frag->dist == 12)
        //  fprintf(stderr, "lib %d sample %d dist "F_S32"\n", frag->dist, dfrg->numSamples, dist);

        dfrg->numSamples++;
        dfrg->mu    += dist;
        dfrg->sigma += (double)dist * (double)dist;
      }
    }  //  over all frags

    ma = NULL;

    // Mark unitigs as potential Rocks and Stones
    //
    //  As of oct 4 2012 this is no longer needed; it's done during unitig loading.  Left in for the
    //  asserts.

    if (operateOnNodes == UNITIG_OPERATIONS) {
      int rock = FALSE;
      int stone = FALSE;
      switch(numExternalLinks){
        case 0:
          stone = TRUE;
          rock = FALSE;
          break;
        case 1:
          stone = TRUE;
          rock = FALSE;
          numPotentialStones++;
          break;
        case 2:
        default:
          stone = rock = TRUE;
          numPotentialRocks++;
          break;
      }

      assert(node->flags.bits.isPotentialRock == rock);
      assert(node->flags.bits.isPotentialStone == stone);
    }
  }  //  over all graph nodes

#if 0
  fprintf(stderr, "* ComputeMatePairStats some mate data:\n");
  fprintf(stderr, "* num total - chaff                 %d\n",numTotalFrags);
  fprintf(stderr, "* num chaff nodes                   %d\n",numChaff);
  fprintf(stderr, "* num singleton                     %d\n",numSingle);
  fprintf(stderr, "* num no link                       %d\n",numNolink);
  fprintf(stderr, "* num ctg not internal              %d\n",numCtgNotInternal);
  fprintf(stderr, "* num 5' > 3' so skip               %d\n",numReverse);
  fprintf(stderr, "* num different unitig              %d\n",numOtherUtg);
  fprintf(stderr, "* num same orientation              %d\n",numSameOrient);
  fprintf(stderr, "* num greater then 5 stddev distant %d\n",numNot5stddev);
  fprintf(stderr, "* num different scaffold            %d\n",numDiffScaf);
  fprintf(stderr, "* num mate != sourceInt             %d\n",numMateNotSource);
  fprintf(stderr, "* num NULL mate                     %d\n",numNullMate);
#endif

#if 0
  if (operateOnNodes == UNITIG_OPERATIONS)
    fprintf(stderr,
            "* ComputeMatePairStats has marked %d/%d unitigs as potential rocks +  %d/%d as potential stones\n",
            numPotentialRocks, (int) GetNumGraphNodes(graph),
            numPotentialStones, (int) GetNumGraphNodes(graph));
#endif

  // now sort the samples, mates, and frags arrays, based on samples
  for (i = 1; i < GetNumDistTs(ScaffoldGraph->Dists); i++) {
    MateInfoT   *matePairs  = NULL;
    int          icnt       = 0;
    int32        newLower   = 0;
    int32        newUpper   = 0;
    int32        median     = 0;
    int32        lowerSigma = 0;
    int32        upperSigma = 0;
    DistT       *dfrg       = &dwork[i];
    int          numSamples = GetNumint32s(dworkSamples[i]);

    if (dfrg->numSamples == 0 || dfrg->numSamples == 1)
      continue;

    matePairs = (MateInfoT *)safe_malloc(sizeof(MateInfoT) * numSamples);

    for ( icnt = 0; icnt<numSamples; icnt++) {
      matePairs[icnt].samples = *Getint32(dworkSamples[i], icnt);
      matePairs[icnt].frags   = *GetCDS_CID_t(dworkFrags[i], icnt);
      matePairs[icnt].mates   = *GetCDS_CID_t(dworkMates[i], icnt);
    }

    qsort(matePairs, numSamples, sizeof(MateInfoT), compDists);

    median     =  matePairs[numSamples / 2].samples;
    lowerSigma =  median - matePairs[ (int) ((0.5 - 0.34) * numSamples)].samples;
    upperSigma = -median + matePairs[ (int) ((0.5 + 0.34) * numSamples)].samples;

    newLower = median - CGW_CUTOFF * MAX(lowerSigma, upperSigma);
    newUpper = median + CGW_CUTOFF * MAX(lowerSigma, upperSigma);
    if (newLower < 0)
      newLower = 0;

#if 0
    {
      DistT       *dorig      = GetDistT(ScaffoldGraph->Dists, i);  //  only for an output
      double       mu         = dfrg->mu / dfrg->numSamples;
      double       sigma      = sqrt((dfrg->sigma - dfrg->numSamples * mu * mu) / (dfrg->numSamples - 1));

      fprintf( stderr, "lib "F_CID", numSamples: %d, orig mean, sig: ( %.2f, %.2f), calc mean, sig: (%.2f, %.2f) median: "F_S32"\n",
               i,
               dfrg->numSamples,
               dfrg->mu,
               dfrg->sigma,
               mu,
               sigma,
               median);
      fprintf( stderr, "dfrg->lower: "F_S32"  dfrg->upper: "F_S32"\n", dfrg->lower, dfrg->upper);
      fprintf( stderr, "dfrg->min  : "F_S32"  dfrg->max  : "F_S32"\n", dfrg->min, dfrg->max);
      fprintf( stderr, "lowerSigma : "F_S32"  upperSigma : "F_S32"\n", lowerSigma, upperSigma);
      fprintf( stderr, "newLower   : "F_S32"  newUpper   : "F_S32"\n", newLower, newUpper);
    }
#endif

    // now reset the trusted flag if necessary
    // first see if there are edges marked untrusted that are now considered trusted
    // lower set
    if (dfrg->lower > newLower) {
      for ( icnt = 0; icnt<numSamples; icnt++) {
        if (matePairs[icnt].samples < dfrg->lower && matePairs[icnt].samples > newLower) {
          CIFragT *frag, *mate;

          frag = GetCIFragT( ScaffoldGraph->CIFrags, matePairs[icnt].frags);
          mate = GetCIFragT( ScaffoldGraph->CIFrags, matePairs[icnt].mates);

          //fprintf( stderr, "1 reclassifying samples[%d] ("F_CID", "F_CID") from UNTRUSTED to TRUSTED\n",
          //         icnt, frag->read_iid, mate->read_iid);

          frag->flags.bits.edgeStatus = TRUSTED_EDGE_STATUS;
          mate->flags.bits.edgeStatus = TRUSTED_EDGE_STATUS;
          dfrg->numBad--;
          dfrg->numSamples++;
          dfrg->mu    += matePairs[icnt].samples;
          dfrg->sigma += (double)matePairs[icnt].samples * (double)matePairs[icnt].samples;
        }
      }
    }

    // upper set
    if (dfrg->upper < newUpper) {
      for ( icnt = 0; icnt<numSamples; icnt++) {
        if (matePairs[icnt].samples > dfrg->upper && matePairs[icnt].samples < newUpper) {
          CIFragT *frag, *mate;

          frag = GetCIFragT( ScaffoldGraph->CIFrags, matePairs[icnt].frags);
          mate = GetCIFragT( ScaffoldGraph->CIFrags, matePairs[icnt].mates);

          //fprintf( stderr, "2 reclassifying samples[%d] ("F_CID", "F_CID") from UNTRUSTED to TRUSTED\n",
          //         icnt, frag->read_iid, mate->read_iid);

          frag->flags.bits.edgeStatus = TRUSTED_EDGE_STATUS;
          mate->flags.bits.edgeStatus = TRUSTED_EDGE_STATUS;
          dfrg->numBad--;
          dfrg->numSamples++;
          dfrg->mu    += matePairs[icnt].samples;
          dfrg->sigma += (double)matePairs[icnt].samples * (double)matePairs[icnt].samples;
        }
      }
    }


    // now see if there are edges marked trusted that are now considered untrusted
    // lower set
    if (dfrg->lower < newLower) {
      for ( icnt = 0; icnt<numSamples; icnt++) {
        if (matePairs[icnt].samples > dfrg->lower && matePairs[icnt].samples < newLower) {
          CIFragT *frag, *mate;

          frag = GetCIFragT( ScaffoldGraph->CIFrags, matePairs[icnt].frags);
          mate = GetCIFragT( ScaffoldGraph->CIFrags, matePairs[icnt].mates);

          //fprintf( stderr, "3 reclassifying samples[%d] ("F_CID", "F_CID") from TRUSTED to UNTRUSTED\n",
          //         icnt, frag->read_iid, mate->read_iid);

          frag->flags.bits.edgeStatus = UNTRUSTED_EDGE_STATUS;
          mate->flags.bits.edgeStatus = UNTRUSTED_EDGE_STATUS;
          dfrg->numBad++;
          dfrg->numSamples--;
          dfrg->mu    -= matePairs[icnt].samples;
          dfrg->sigma -= (double)matePairs[icnt].samples * (double)matePairs[icnt].samples;
        }
      }
    }
    // upper set
    if (dfrg->upper > newUpper) {
      for ( icnt = 0; icnt<numSamples; icnt++) {
        if (matePairs[icnt].samples < dfrg->upper && matePairs[icnt].samples > newUpper) {
          CIFragT *frag, *mate;

          frag = GetCIFragT( ScaffoldGraph->CIFrags, matePairs[icnt].frags);
          mate = GetCIFragT( ScaffoldGraph->CIFrags, matePairs[icnt].mates);

          //fprintf( stderr, "4 reclassifying samples[%d] ("F_CID", "F_CID") from TRUSTED to UNTRUSTED\n",
          //         icnt, frag->read_iid, mate->read_iid);

          frag->flags.bits.edgeStatus = UNTRUSTED_EDGE_STATUS;
          mate->flags.bits.edgeStatus = UNTRUSTED_EDGE_STATUS;
          dfrg->numBad++;
          dfrg->numSamples--;
          dfrg->mu    -= matePairs[icnt].samples;
          dfrg->sigma -= (double)matePairs[icnt].samples * (double)matePairs[icnt].samples;
        }
      }
    }
    safe_free(matePairs);
  }  //  over all distances

  //
  // now set mean, stddev, size & number of buckets, and make an update for gatekeeper
  //

  char   updName[FILENAME_MAX];

  sprintf(updName, "stat/%s.distupdate", instance_label);

  FILE *updFile = fopen(updName, "w");

  for (i=1; i<GetNumDistTs(ScaffoldGraph->Dists); i++) {
    DistT        *dfrg = &dwork[i];

    int           numSamples             = 0;
    int32        *samples                = NULL;

    //  If we don't have enough samples for an update, don't do anything to dfrg

    if (dfrg->numSamples <= minSamplesForOverride) {
      fprintf(stderr, "ComputeMatePairStatisticsRestricted()-- LIB IID %d distance %.2f %.2f -- only %u samples, can't reestimate\n", i, dfrg->mu, dfrg->sigma, dfrg->numSamples);

    } else {
      DistT  *dupd       = GetDistT(ScaffoldGraph->Dists, i);

      double  newmu      = dfrg->mu / dfrg->numSamples;
      double  newsigma   = sqrt((dfrg->sigma - dfrg->mu * dfrg->mu / dfrg->numSamples) / (dfrg->numSamples - 1));
      
      if (dupd->allowUpdate) {
        dupd->mu           = newmu;
        dupd->sigma        = newsigma;
      }
      dupd->numSamples     = dfrg->numSamples;
      dupd->min            = dfrg->min;
      dupd->max            = dfrg->max;
      dupd->bnum           = 0;
      dupd->bsize          = 0;
      dupd->lower          = dupd->mu - CGW_CUTOFF * dupd->sigma;
      dupd->upper          = dupd->mu + CGW_CUTOFF * dupd->sigma;
      dupd->numBad         = dfrg->numBad;
      dupd->allowUpdate    = dfrg->allowUpdate;

      dupd->bnum           = 1;
      dupd->bsize          = dfrg->max - dfrg->min;

      if (dupd->sigma > CGW_NUM_BUCKETS)
        dupd->bnum = dupd->bsize * CGW_NUM_BUCKETS / dupd->sigma + 1;

      dupd->bsize /= dupd->bnum;

      //  Remove any existing histogram, then reallocate (and clear)
      //  one big enough.

      safe_free(dupd->histogram);
      dupd->histogram = (int32 *)safe_calloc(dupd->bnum, sizeof(int32));

      // output a histogram file for each library

      numSamples = GetNumint32s(dworkSamples[i]);
      samples    = Getint32(dworkSamples[i],0);

      for (j=0; j<numSamples ; j++) {
        int32 binNum = (samples[j] - dupd->min) / (double)dupd->bsize;

        binNum = MIN(binNum, dupd->bnum - 1);
        binNum = MAX(binNum,0);

        dupd->histogram[binNum]++;
      }

      numSamples = GetNumint32s(dworkSamples[i]);
      samples    = Getint32(dworkSamples[i],0);

      qsort(samples, numSamples, sizeof(int32), &compareInt);

      sprintf(updName, "stat/%s.distlib_%d.cgm", instance_label, i);
      FILE *cgmFile = fopen(updName, "w");

      if (cgmFile) {
        fprintf(cgmFile, "lib %d mu %g sigma %g%s\n", i, newmu, newsigma, dupd->allowUpdate ? "" : " (NOT UPDATED)");
        for (j=0; j<numSamples; j++)
          fprintf(cgmFile, "%d\n", samples[j]);
        fclose(cgmFile);
      }

      if (updFile) {
        if (dupd->allowUpdate)
          fprintf(updFile, "lib iid %d distance %.2f %.2f\n", i, dupd->mu, dupd->sigma);
        else
          fprintf(updFile, "#lib iid %d distance %.2f %.2f # would have updated to %.2f %.2f but not allowed\n", i, dupd->mu, dupd->sigma, newmu, newsigma);
      }

      if (dupd->allowUpdate)
        fprintf(stderr, "ComputeMatePairStatisticsRestricted()-- LIB IID %d distance %.2f %.2f\n", i, dupd->mu, dupd->sigma);
      else
        fprintf(stderr, "ComputeMatePairStatisticsRestricted()-- LIB IID %d distance %.2f %.2f -- would have updated to %.2f %.2f but not allowed\n", i, dupd->mu, dupd->sigma, newmu, newsigma);
    }  //  end of update

    DeleteVA_CDS_CID_t(dworkSamples[i]);
    DeleteVA_CDS_CID_t(dworkFrags[i]);
    DeleteVA_CDS_CID_t(dworkMates[i]);
  }  //  over all distances

  fclose(updFile);

  delete [] dwork;
  delete [] dworkSamples;
  delete [] dworkFrags;
  delete [] dworkMates;
}





// Clean hash table
// Move deleted nodes and edges to their respective free lists
void RecycleDeletedGraphElements(GraphCGW_T *graph){
  EdgeCGW_T *end = NULL;

  // Find the end of the tobefree list (this should be MUCH shorter
  // than looking for teh end of the free list )
  CDS_CID_t nextIndex = graph->tobeFreeEdgeHead;

  while(nextIndex != NULLINDEX){
    end = GetGraphEdge(graph, nextIndex);
    assert(end->flags.bits.isDeleted);
    assert(nextIndex != end->referenceEdge); // We are in a cycle!
    nextIndex = end->referenceEdge;
  }

  // Tack the tobe freed edges onto the head of the free list
  if(end){
    end->referenceEdge = graph->freeEdgeHead;
    graph->freeEdgeHead = graph->tobeFreeEdgeHead;
  }else{
    ; // nothing to do, no
  }

  graph->tobeFreeEdgeHead = NULLINDEX;

  // Flush Deleted nodes from the hashtable
  // Free tobe freed nodes
}




#define SEGLEN 50
static void dumpFastaRecord(FILE *stream, char *header, char *sequence){
  int i;
  int32 FragLen = (int32) strlen(sequence);

  fprintf(stream,"> %s\n", header);
  for (i = 0; i < FragLen; i += SEGLEN)
    if (FragLen-i < SEGLEN)
      fprintf(stream,"%.*s\n",FragLen-i,sequence+i);
    else
      fprintf(stream,"%.*s\n",SEGLEN,sequence+i);
}






/* Restore the Edge Means we squirreled away during walking */
void RestoreEdgeMeans(GraphCGW_T *graph){
  int i;
  EdgeCGW_T *edge;
  for(i = 0; i < GetNumGraphEdges(graph); i++){
    edge = GetGraphEdge(graph,i);
    if(edge->flags.bits.isDeleted ||
       !edge->flags.bits.MeanChangedByWalking)
      continue;

    edge->distance.mean = edge->minDistance;
    edge->flags.bits.MeanChangedByWalking = FALSE;
  }
}



int32 EdgeDegree(GraphCGW_T *graph, EdgeCGW_T *edge) {
  EdgeCGW_T *e = edge;
  int32 tdegree = 0;

  if(edge->flags.bits.isRaw){
    tdegree = 1;
  }else{
    for(e = GetGraphEdge(graph, edge->nextRawEdge);
        e != NULL;
        e = GetGraphEdge(graph, e->nextRawEdge))
      tdegree++;
  }
  return(tdegree);
}



/* ---------------------------------------------
 * Convert values of type: unsigned int
 * to the enumerated type: EdgeStatus
 *
 * Returns INVALID_EDGE_STATUS on bad input.
 ---------------------------------------------- */
EdgeStatus AS_CGW_SafeConvert_uintToEdgeStatus(unsigned int input) {
  EdgeStatus output =            INVALID_EDGE_STATUS;
  if (input<=4) {
    if (input==1) output =       UNKNOWN_EDGE_STATUS;
    else if (input==2) output =  UNTRUSTED_EDGE_STATUS;
    else if (input==4) output =  TENTATIVE_UNTRUSTED_EDGE_STATUS;
  } else {
    if (input==8) output =       TENTATIVE_TRUSTED_EDGE_STATUS;
    else if (input==16) output = TRUSTED_EDGE_STATUS;
    else if (input==32) output = LARGE_VARIANCE_EDGE_STATUS;
    else if (input==64) output = INTER_SCAFFOLD_EDGE_STATUS;
  }
  return output;
}



// some checks on unitigs and their frags
void  CheckUnitigs(CDS_CID_t startUnitigID) {
  GraphNodeIterator Nodes;
  NodeCGW_T *node;
  int32 i;
  int totalUnitigs = (int) GetNumGraphNodes( ScaffoldGraph->CIGraph );

  InitGraphNodeIterator( &Nodes, ScaffoldGraph->CIGraph, GRAPH_NODE_DEFAULT);
  while(NULL != ( node = NextGraphNodeIterator(&Nodes))) {
    if (node->id < startUnitigID)
      continue;

    if (node->flags.bits.isChaff && GlobalData->ignoreChaffUnitigs)
      continue;

    MultiAlignT *ma = ScaffoldGraph->tigStore->loadMultiAlign(node->id, TRUE);

    for(i = 0; i < GetNumIntMultiPoss(ma->f_list); i++) {
      IntMultiPos *mp = GetIntMultiPos(ma->f_list, i);
      CIFragT *frag = GetCIFragT(ScaffoldGraph->CIFrags, mp->ident);

      if (frag->cid != node->id) {
        fprintf(stderr, "* CheckUnitigs: frag "F_CID" with iid "F_CID" of cid "F_CID" has ci = "F_CID"!!!\n",
                mp->ident, frag->read_iid, node->id, frag->cid);
        //assert(0);
      }
    }
  }
}


// Check on the sanity of the surrogate unitigs
void  CheckSurrogateUnitigs() {
  GraphNodeIterator Nodes;
  NodeCGW_T *curChunk;

  InitGraphNodeIterator( &Nodes, ScaffoldGraph->CIGraph, GRAPH_NODE_DEFAULT);

  while (NULL != (curChunk = NextGraphNodeIterator(&Nodes))) {
    if (curChunk->type != UNRESOLVEDCHUNK_CGW)
      continue;
    if (curChunk->info.CI.numInstances < 3)
      continue;

    if (curChunk->info.CI.numInstances != GetNumCDS_CID_ts(curChunk->info.CI.instances.va))
      fprintf(stderr, "CheckSurrogateUnitigs()--  CI.numInstances=%d != instances.va="F_SIZE_T"\n",
              curChunk->info.CI.numInstances,
              GetNumCDS_CID_ts(curChunk->info.CI.instances.va));
    assert(curChunk->info.CI.numInstances != GetNumCDS_CID_ts(curChunk->info.CI.instances.va));
  }
}

