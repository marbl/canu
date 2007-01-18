
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
/*
  cc -g -pg -qfullpath   -qstrict -qbitfields=signed -qchars=signed -qlanglvl=ext -qignerrno -qupconv -qcpluscmt -qmaxmem=-1 -D_ILS_MACROS -D_C_LOCALE_ONLY -c
  -o /work/assembly/rbolanos/IBM_PORT_CDS/ibm_migration_work_dir/cds/AS/obj/GraphCGW_T.o GraphCGW_T.c
*/

static char CM_ID[] = "$Id: GraphCGW_T.c,v 1.25 2007-01-18 20:30:59 eliv Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <fcntl.h>
#include <sys/types.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

#include "ScaffoldGraphIterator_CGW.h"
#include "ScaffoldGraph_CGW.h"
#include "GraphCGW_T.h"
#include "ChunkOverlap_CGW.h"
#include "AS_PER_SafeIO.h"
#include "AS_UTL_interval.h"
#include "MultiAlignment_CNS.h"
#include "AS_MSG_Utility.h"
#include "FbacREZ.h"
#include "UtilsREZ.h"

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
  fprintf(stderr," Reallocating %d edges, currently %d\n",
          numEdges, (int) GetNumGraphEdges(graph));
  MakeRoom_VA(graph->edges, (size_t) numEdges, FALSE);

}

GraphCGW_T *CreateGraphCGW(GraphType type,
                           int32 numNodes, int32 numEdges){
  GraphCGW_T *graph = (GraphCGW_T *)safe_malloc(sizeof(GraphCGW_T));

  InitializeGraph(graph);
  graph->type = type;
  graph->nodes = CreateVA_NodeCGW_T((size_t) MAX(1024, numNodes));
  graph->edges = CreateVA_EdgeCGW_T((size_t) MAX(1024, numNodes));
#ifdef USE_UNITIG_FILE
  graph->chunkFile =   CreateUnitigFile();
  assert(graph->chunkFile);
#endif
  if(graph->type == SCAFFOLD_GRAPH){
    //    graph->maStore = NULL;
    graph->overlapper = NULL;
  }else{
    //    graph->maStore =   CreateMultiAlignStoreT();
    //    AssertPtr(graph->maStore);
    graph->overlapper = CreateChunkOverlapper();
  }

  return graph;

}


void SaveGraphCGWToStream(GraphCGW_T *graph, FILE *stream){
  int status;
  int32 i;
  CopyToFileVA_NodeCGW_T(graph->nodes, stream);

  // Save lists of indicies, if they exist
  if(graph->type == CI_GRAPH)
    for(i = 0; i < GetNumGraphNodes(graph); i++){
      NodeCGW_T *node = GetGraphNode(graph,i);

      if(!node->flags.bits.isCI ||
	 node->info.CI.numInstances < 3)
	continue;
      CopyToFileVA_CDS_CID_t(node->info.CI.instances.va, stream);
    }

  CopyToFileVA_EdgeCGW_T(graph->edges, stream);
  
  status = safeWrite(stream, &graph->type, sizeof(int32));
  assert(status == FALSE);
  status = safeWrite(stream, &graph->numActiveNodes, sizeof(int32));
  assert(status == FALSE);
  status = safeWrite(stream, &graph->numActiveEdges, sizeof(int32));
  assert(status == FALSE);
  status = safeWrite(stream, &graph->freeEdgeHead, sizeof(CDS_CID_t));
  assert(status == FALSE);
  status = safeWrite(stream, &graph->tobeFreeEdgeHead, sizeof(CDS_CID_t));
  assert(status == FALSE);
  status = safeWrite(stream, &graph->freeNodeHead, sizeof(CDS_CID_t));
  assert(status == FALSE);
  status = safeWrite(stream, &graph->tobeFreeNodeHead, sizeof(CDS_CID_t));
  assert(status == FALSE);
  status = safeWrite(stream, &graph->deadNodeHead, sizeof(CDS_CID_t));
  assert(status == FALSE);

  // Save the multiAlignStore
  if(graph->type != SCAFFOLD_GRAPH){
    //    SaveMultiAlignStoreTToStream(graph->maStore, stream, (graph->type == CI_GRAPH));
    SaveChunkOverlapperToStream(graph->overlapper,stream);
  }

}


GraphCGW_T *LoadGraphCGWFromStream(FILE *stream){
  int status;
  CDS_CID_t i;
  GraphCGW_T *graph = (GraphCGW_T *)safe_malloc(sizeof(GraphCGW_T));

  graph->nodes =        CreateFromFileVA_NodeCGW_T(stream,0);

  // Save lists of indicies, if they exist
  for(i = 0; i < GetNumGraphNodes(graph); i++){
    NodeCGW_T *node = GetGraphNode(graph,i);
    if(!node->flags.bits.isCI )
      continue;
    if(node->info.CI.numInstances < 3)
      continue;
    fprintf(GlobalData->stderrc,"* Loading instance list for CI " F_CID " of length %d from VA\n",
	    i,node->info.CI.numInstances);
    node->info.CI.instances.va = CreateFromFileVA_CDS_CID_t(stream,0);
  }



  graph->edges =          CreateFromFileVA_EdgeCGW_T( stream,0);

  status = safeRead(stream, &graph->type, sizeof(int32));
  assert(status == FALSE);
  status = safeRead(stream, &graph->numActiveNodes, sizeof(int32));
  assert(status == FALSE);
  status = safeRead(stream, &graph->numActiveEdges, sizeof(int32));
  assert(status == FALSE);
  status = safeRead(stream, &graph->freeEdgeHead, sizeof(CDS_CID_t));
  assert(status == FALSE);
  status = safeRead(stream, &graph->tobeFreeEdgeHead, sizeof(CDS_CID_t));
  assert(status == FALSE);
  status = safeRead(stream, &graph->freeNodeHead, sizeof(CDS_CID_t));
  assert(status == FALSE);
  status = safeRead(stream, &graph->tobeFreeNodeHead, sizeof(CDS_CID_t));
  assert(status == FALSE);
  status = safeRead(stream, &graph->deadNodeHead, sizeof(CDS_CID_t));
  assert(status == FALSE);

  // Load the multiAlignStore
  if(graph->type == CI_GRAPH){
    //    graph->maStore = LoadMultiAlignStoreTFromStream(stream);
  }else if(graph->type == CONTIG_GRAPH){
    //    graph->maStore = LoadMultiAlignStoreTFromStreamWithReferences(stream, ScaffoldGraph->CIGraph->maStore);
  }// else graph->maStore = NULL;

  graph->overlapper = NULL;
  if(graph->type != SCAFFOLD_GRAPH)
    graph->overlapper = LoadChunkOverlapperFromStream(stream);

  return graph;

}


void DeleteGraphCGW(GraphCGW_T *graph){
  int32 i;

  // Delete lists of indicies, if they exist
  if(graph->type == CI_GRAPH)
    for(i = 0; i < GetNumGraphNodes(graph); i++){
      NodeCGW_T *node = GetGraphNode(graph,i);
      if(!node->flags.bits.isCI ||
	 node->info.CI.numInstances < 3)
	continue;
      DeleteVA_CDS_CID_t(node->info.CI.instances.va);
    }

  DeleteVA_NodeCGW_T(graph->nodes);
  DeleteVA_EdgeCGW_T(graph->edges);
  if(graph->type != SCAFFOLD_GRAPH){
    fprintf(GlobalData->stderrc,"* Destroying Overlapper\n");
    DestroyChunkOverlapper(graph->overlapper);
    //    fprintf(GlobalData->stderrc,"* Destroying Multi Align Store\n");
    //    DeleteMultiAlignStoreT(graph->maStore);
  }
  free(graph);
}

/* Diagnostic */
size_t ReportMemorySizeGraphCGW(GraphCGW_T *graph, FILE *stream){
  char *nodeName;
  char *edgeName;
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
#ifdef USE_UNITIG_FILE
  totalMemorySize += StatsUnitigFile(graph->chunkFile, stream);
#endif
  if(graph->type != SCAFFOLD_GRAPH){
    //    totalMemorySize += StatsMultiAlignStore(graph->maStore, stream, /* CI_GRAPH always owns its contents */(graph->type == CI_GRAPH));
    totalMemorySize += ReportMemorySize_CO(graph->overlapper, stream);
  }
  return totalMemorySize;
}


int FindGraphEdgeChain(GraphCGW_T *graph, CDS_CID_t eid,
                       VA_TYPE(CDS_CID_t) *chain,
                       int extractSingletons, int includeGuides){
  CIEdgeT *edge, *head, *tail;
  int32 count = 0;
  CDS_CID_t headID, edgeID;
  edge = head = GetGraphEdge(graph, eid);
  AssertPtr(edge);
  headID = eid;

  edgeID = eid;

  assert(head->idA <= head->idB); // Only interesting on canonical chains

  /*  fprintf(GlobalData->stderrc,"* FindGraphEdgeCHain eid:" F_CID " (" F_CID "," F_CID ",%c) prev:(" F_CID "," F_CID ") next:(" F_CID "," F_CID ")\n",
      eid, edge->idA, edge->idB,edge->orient,
      edge->prevALink, edge->prevBLink,
      edge->nextALink, edge->nextBLink);
  */
  // Make sure this edge is correctly hooked into the graph
  assert(GraphEdgeSanity(graph, eid));

  /* A chain is characterized by both the A and B links pointing to the next
     element in the chain */
  while(1) {
    count++;
    AppendCDS_CID_t(chain, &edgeID);
    //	fprintf(GlobalData->stderrc,
    //		"* Appended " F_CID " (" F_CID "," F_CID ",%d)to chain\n", edgeID, edge->idA, edge->idB, edge->orient);

    tail = edge;

    if(edge->nextALink == NULLINDEX ||
       edge->nextBLink == NULLINDEX) 
      break;

    if( edge->nextALink != edge->nextBLink)
      break;

    // Guides are always at the end of these chains
    if(!includeGuides && isSloppyEdge(edge))
      break;

    edgeID = edge->nextALink;
    edge = GetGraphEdge(graph, edgeID);

    if( edge->orient != head->orient)
      break;

    // Guides are always at the end of these chains
    if(!includeGuides && isSloppyEdge(edge))
      break;

  }
  
  if(count > 1 || extractSingletons){  // Unlink the entire chain of suckers
    ChunkInstanceT *chunkA, *chunkB;
    CIEdgeT *prevA, *prevB, *nextA, *nextB;
    // head->prevALink should be hooked to tail->nextALink
    // head->prevBLink should be hooked to tail->nextBlink
    // if head->prev?Link is NULLINDEX, we need to update the head of the list
    chunkA = GetGraphNode(graph, edge->idA);
    chunkB = GetGraphNode(graph, edge->idB);
    assert(chunkA && chunkB);
    
    /* fprintf(GlobalData->stderrc,"* Unlinking %d edges... \n\t* ", count);
       {
       CDS_CID_t *edges = GetCDS_CID_t(chain,0);
       int i;
       for(i = 0; i < count; i++){
       fprintf(GlobalData->stderrc,F_CID " ", edges[i]);
       }
       fprintf(GlobalData->stderrc,"\n");
       }
    */
    nextA = NULL;
    if(tail->nextALink != NULLINDEX)
      nextA = GetGraphEdge(graph, tail->nextALink);

    nextB = NULL;
    if(tail->nextBLink != NULLINDEX)
      nextB = GetGraphEdge(graph, tail->nextBLink);

    prevA = NULL;
    if(head->prevALink != NULLINDEX)
      prevA = GetGraphEdge(graph, head->prevALink);

    prevB = NULL;
    if(head->prevBLink != NULLINDEX)
      prevB = GetGraphEdge(graph, head->prevBLink);

    // Unlink A Chain
    if(!prevA){
      chunkA->edgeHead = tail->nextALink;
    }else{
      if(prevA->idA == head->idA){
	prevA->nextALink = tail->nextALink;
      }else{
	prevA->nextBLink = tail->nextALink;
      }
    }
    if(nextA){
      if(nextA->idA == head->idA){
	nextA->prevALink = head->prevALink;
      }else{
	nextA->prevBLink = head->prevALink;
      }
    }
    // Unlink B Chain
    if(!prevB){
      chunkB->edgeHead = tail->nextBLink;
    }else{
      if(prevB->idB == head->idB){
	prevB->nextBLink = tail->nextBLink;
      }else{
	prevB->nextALink = tail->nextBLink;
      }
    }

    if(nextB){
      if(nextB->idB == head->idB){
	nextB->prevBLink = head->prevBLink;
      }else{
	nextB->prevALink = head->prevBLink;
      }
    }


  }


  return count;

}


/* Check that edge with index eid is properly wired in the graph:
   - we find it when looking for it in both lists which it is supposed to be a member
   - we can walk backwards from the edge to the heads of the two lists
*/
int GraphEdgeSanity(GraphCGW_T *graph, CDS_CID_t eid){
  CIEdgeT *edge = GetGraphEdge(graph,eid);
  CIEdgeT *edgeFromA, *edgeFromB;
  GraphEdgeIterator edges;
  CDS_CID_t idA = edge->idA;
  CDS_CID_t idB = edge->idB;
  
  InitGraphEdgeIterator(graph, idA, ALL_END, ALL_EDGES,
                        GRAPH_EDGE_DEFAULT , &edges);
  while(NULL != (edgeFromA = NextGraphEdgeIterator(&edges))){
    if(edgeFromA == edge)
      break;
  }

  if(edgeFromA == NULL){
    fprintf(GlobalData->stderrc,"* Couldn't find edge " F_CID " (" F_CID "," F_CID ") looking from " F_CID "\n",
	    eid, idA, idB, idA);
  }

  InitGraphEdgeIterator(graph, idB, ALL_END, ALL_EDGES,
                        GRAPH_EDGE_DEFAULT , &edges);
  while(NULL != (edgeFromB = NextGraphEdgeIterator(&edges))){
    if(edgeFromB == edge)
      break;
  }

  if(edgeFromB == NULL){
    fprintf(GlobalData->stderrc,"* Couldn't find edge " F_CID " (" F_CID "," F_CID ") looking from " F_CID "\n",
	    eid, idA, idB, idB);
  }

  if(edgeFromB == NULL ||
     edgeFromA == NULL){
    PrintGraphEdge(GlobalData->stderrc,graph," ", edge, idA);
    return FALSE;
  }

  return TRUE;
}


void SetEdgeCGWLinks(CIEdgeT *newEdge, int isA,
                     CDS_CID_t prev, CDS_CID_t next){
  assert(prev >= 0 || prev == NULLINDEX);
  assert(next >= 0 || next == NULLINDEX);
  assert(!newEdge->flags.bits.isDeleted);

  if(isA){
    newEdge->nextALink = next;
    newEdge->prevALink = prev;
  }else{
    newEdge->nextBLink = next;
    newEdge->prevBLink = prev;
  }
}

int EdgeCGWCompare(EdgeCGW_T *edge1, EdgeCGW_T *edge2){
  int diff;

  assert(edge1->idA <= edge1->idB &&
	 edge2->idA <= edge2->idB &&
	 !edge1->flags.bits.isDeleted &&
	 !edge2->flags.bits.isDeleted);


  diff = edge1->idA - edge2->idA;
  if(diff)
    return diff;

  diff = edge1->idB - edge2->idB;
  if(diff)
    return diff;

  diff = edge1->orient - edge2->orient;
  if(diff)
    return diff;

  // We want guide edges AFTER all other edges
  diff = isSloppyEdge(edge1) - isSloppyEdge(edge2);
  if(diff)
    return diff;
  
  // We want overlap edges AFTER non-overlap edges...
  diff = isOverlapEdge(edge1) - isOverlapEdge(edge2);
  if(diff)
    return diff;

  diff = edge1->distance.mean - edge2->distance.mean;
  return diff;
}


void InsertGraphEdgeInList(GraphCGW_T *graph,
                           CDS_CID_t CIedgeID, CDS_CID_t cid, int verbose){
  ChunkInstanceT *ci;
  GraphEdgeIterator edges;
  CIEdgeT *edge, *prevEdge, *newEdge;
  int isA;
  int flags;
  newEdge = GetGraphEdge(graph, CIedgeID);

  assert(newEdge->idA == cid || newEdge->idB == cid);
  assert(!newEdge->flags.bits.isDeleted);

  isA = (newEdge->idA == cid);


  if(isA){
    newEdge->prevALink = NULLINDEX;
    newEdge->nextALink = NULLINDEX;
  }else{
    newEdge->prevBLink = NULLINDEX;
    newEdge->nextBLink = NULLINDEX;
  }

  // First insert into list for idA
  // If idA's list is empty, this is easy
  
  ci = GetGraphNode(graph, cid);
  AssertPtr(ci);
  assert(!ci->flags.bits.isDead);

  
  /* Should we insert at head? */
  flags = GRAPH_EDGE_DEFAULT;

  InitGraphEdgeIterator(graph, cid, ALL_END, ALL_EDGES, flags, &edges);
  edge = NextGraphEdgeIterator(&edges);

  /* Should we insert at head? */
  if(!edge ||
     EdgeCGWCompare(newEdge, edge) <= 0){
    SetEdgeCGWLinks(newEdge, isA, NULLINDEX, ci->edgeHead);
    ci->edgeHead = CIedgeID;
    if(edge){
      if(edge->idA == cid){
	edge->prevALink = CIedgeID;
      }else{
	assert(edge->idB == cid);
	edge->prevBLink = CIedgeID;
      }
    }
    if(verbose)
      fprintf(GlobalData->stderrc,"* cid:" F_CID " edgeID:" F_CID " Insert at head (" F_CID "," F_CID ") -> (" F_CID "," F_CID ")\n",
	      cid, CIedgeID,
	      newEdge->idA, newEdge->idB,
              (edge?edge->idA:-1), (edge?edge->idB:-1));
  }else{
    /* Look down list until we either hit the end, or are > the edge */
    CDS_CID_t eid = GetVAIndex_EdgeCGW_T(graph->edges, edge);

    if(ci->edgeHead != eid){
      EdgeCGW_T *ehead = GetGraphEdge(graph, ci->edgeHead);
      EdgeCGW_T *eprev = GetGraphEdge(graph, edges.prev);
      fprintf(GlobalData->stderrc,"*  Screwup!!!\n");
      PrintGraphEdge(GlobalData->stderrc,graph,"head ", ehead, cid);
      PrintGraphEdge(GlobalData->stderrc,graph,"1st  ", edge, cid);
      if(eprev)
	PrintGraphEdge(GlobalData->stderrc,graph,"1st  ", eprev, cid);

      assert(0);
    }

    assert(!edge->flags.bits.isDeleted);
  
    while(NULL != (edge = NextGraphEdgeIterator(&edges))){
      assert(!edge->flags.bits.isDeleted);
      assert(edge != newEdge);
      assert(edge->idA == cid || edge->idB == cid);
      if( EdgeCGWCompare(newEdge,edge) <= 0)
	break;
    }
    /* Now insert */
    prevEdge = GetGraphEdge(graph, edges.prev);
    assert(prevEdge != newEdge && !prevEdge->flags.bits.isDeleted);
    if(verbose){
      if(edge){
	fprintf(GlobalData->stderrc,"* cid:" F_CID " edge:" F_CID " Inserting in middle " F_CID "(" F_CID "," F_CID ") - >* " F_CID "(" F_CID "," F_CID ") -> * " F_CID "(" F_CID "," F_CID ")\n",
		cid, CIedgeID,
		edges.prev,
		prevEdge->idA, prevEdge->idB,
		CIedgeID,
		newEdge->idA, newEdge->idB,
		edges.curr,
		edge->idA, edge->idB);
      }else{
	fprintf(GlobalData->stderrc,"* cid:" F_CID " edge:" F_CID " Appending after (" F_CID "," F_CID ") - >* (" F_CID "," F_CID ")\n",
		cid, CIedgeID,
		prevEdge->idA, prevEdge->idB,
		newEdge->idA, newEdge->idB);
      }
    }
    AssertPtr(prevEdge);
    if(prevEdge->idA == cid){
      SetEdgeCGWLinks(newEdge, isA, edges.prev, prevEdge->nextALink);
      prevEdge->nextALink = CIedgeID;
    }else{
      assert(prevEdge->idB == cid);
      SetEdgeCGWLinks(newEdge, isA, edges.prev, prevEdge->nextBLink);
      prevEdge->nextBLink = CIedgeID;
    }
    if(edge){
      if(edge->idA == cid){
	edge->prevALink = CIedgeID;
      }else{
	edge->prevBLink = CIedgeID;
      }
    }
    if(prevEdge){
      assert(prevEdge->nextALink == NULLINDEX ||
             (prevEdge->nextALink != prevEdge->prevALink));
      assert(prevEdge->nextBLink == NULLINDEX ||
             (prevEdge->nextBLink != prevEdge->prevBLink));
      assert( (prevEdge->nextALink != edges.prev &&
               prevEdge->nextBLink != edges.prev));
      assert( (prevEdge->prevALink != edges.prev &&
               prevEdge->prevBLink != edges.prev));
    }

#ifdef NEVER
    // previous version of asserts
    assert(newEdge->nextALink != CIedgeID && newEdge->nextBLink != CIedgeID);
    assert(newEdge->prevALink != CIedgeID && newEdge->prevBLink != CIedgeID);
#else
    if(isA)
      assert(newEdge->prevALink != CIedgeID && newEdge->nextALink != CIedgeID);
    else
      assert(newEdge->prevBLink != CIedgeID && newEdge->nextBLink != CIedgeID);
#endif
  }
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
      if((edge->edgesContributing > 2) ||
	 (!isOverlapEdge(edge) && (edge->edgesContributing == 2))){
	edge->flags.bits.isConfirmed = TRUE;
      }else{
	edge->flags.bits.isConfirmed = FALSE;
      }
    }
  }

  // Default state for isEssential is FALSE will be set later.
  edge->flags.bits.isEssential = FALSE;
  // Default state for wasEssential is FALSE will be set later.
  edge->flags.bits.wasEssential = FALSE;
  // Default state for isInferred is FALSE will be set later.
  edge->flags.bits.isInferred = FALSE;
  // Default state for isInferredRemoved is FALSE will be set later.
  edge->flags.bits.isInferredRemoved = FALSE;
  // Default state for isRedundantRemoved is FALSE will be set later.
  edge->flags.bits.isRedundantRemoved = FALSE;
  // Default state for isTransitivelyRemoved is FALSE will be set later.
  edge->flags.bits.isTransitivelyRemoved = FALSE;
  // Default state for isDeleted is FALSE will be set later.
  edge->flags.bits.isDeleted = FALSE;

  return;
}

// Insert a copy of the tobeInserted Mate Edge 
// Returns index of CIEdge if successful
// Edge must be 'canonical' (idA < idB)
//
int32 InsertGraphEdge(GraphCGW_T *graph, CDS_CID_t edgeID, int verbose){
  EdgeCGW_T *edge = GetGraphEdge(graph, edgeID);

  assert(GetEdgeStatus(edge) != INVALID_EDGE_STATUS);
  assert(edge->idA < edge->idB);
  assert(edge->topLevelEdge == edgeID);

  // Initialize the status flags for this edge
  InitGraphEdgeFlags(graph, edge);
  
  // ALD: set  [prev|next][AB]Link 's to NULL since
  // InsertGraphEdgeInList  isn't very smart about keeping the
  // A|B versions separate
  edge -> prevALink = edge -> nextALink
    = edge -> prevBLink = edge -> nextBLink = NULLINDEX;

  // Now insert into list for idA
  InsertGraphEdgeInList(graph,  edgeID, edge->idA, verbose);
  // Now insert into list for idB
  InsertGraphEdgeInList(graph,  edgeID, edge->idB, verbose);

  return TRUE;
}


void FreeGraphEdgeByEID(GraphCGW_T *graph, CDS_CID_t eid){
  int count = 0;
  EdgeCGW_T *edge = GetGraphEdge(graph, eid);
#ifdef DEBUG 
  PrintGraphEdge(GlobalData->stderrc,graph,"Freeing", edge, edge->idA);
#endif
  assert(!edge->flags.bits.isDeleted);

  if(!edge->flags.bits.isRaw)
    {
      EdgeCGW_T *rawEdge = GetGraphEdge(graph, edge->nextRawEdge);
      while(rawEdge){
	EdgeCGW_T *nextEdge = GetGraphEdge(graph, rawEdge->nextRawEdge);
#ifdef DEBUG
	PrintGraphEdge(GlobalData->stderrc,graph,"Freeing raw",
                       rawEdge, rawEdge->idA);
#endif      
	count++;
	//      rawEdge->flags.bits.isDeleted = TRUE;
	//      rawEdge->flags.bits.isRaw = TRUE;

	assert(rawEdge->idA == edge->idA &&
	       rawEdge->idB == edge->idB &&
	       rawEdge->orient == edge->orient );

	rawEdge->nextRawEdge = NULLINDEX;
	FreeGraphEdge(graph,rawEdge);
	rawEdge = nextEdge;
      }
      /*    if(count != edge->edgesContributing){
            fprintf(GlobalData->stderrc,
            "* contributing = %d only deleted %d raw edges!!!\n",
            edge->edgesContributing, count);
	    }
      */
    }
  edge->flags.bits.isDeleted = TRUE;
    
  assert(eid != graph->tobeFreeEdgeHead);

  edge->referenceEdge = graph->tobeFreeEdgeHead;
  edge->nextRawEdge = NULLINDEX;
  edge->flags.bits.isRaw = TRUE;
  edge->idA = edge->idB = NULLINDEX;
  graph->tobeFreeEdgeHead = eid;
}

// Initialize a Graph Edge
void InitGraphEdge(EdgeCGW_T *edge){
  edge->idA = edge->idB = NULLINDEX;
  edge->flags.all = 0;
  edge->nextALink = edge->nextBLink = NULLINDEX;
  edge->prevALink = edge->prevBLink = NULLINDEX;
  edge->quality   = 1.0;
  edge->nextRawEdge = NULLINDEX;
  edge->referenceEdge = NULLINDEX;
  edge->topLevelEdge = NULLINDEX;
}

// Get the edge from the free list
EdgeCGW_T *GetFreeGraphEdge(GraphCGW_T *graph){
  CDS_CID_t eid = graph->freeEdgeHead;
  EdgeCGW_T *freeEdge = GetGraphEdge(graph, eid);
  CDS_CID_t newEid;

  if(!freeEdge){
    freeEdge = CreateNewGraphEdge(graph);
  }else{
    graph->freeEdgeHead = freeEdge->referenceEdge;
  }
  InitGraphEdge(freeEdge);
  freeEdge->flags.all = 0;
  freeEdge->flags.bits.isRaw = TRUE;
  freeEdge->nextRawEdge = NULLINDEX;
  freeEdge->edgesContributing = 1;
  freeEdge->referenceEdge = NULLINDEX;
  newEid = GetVAIndex_EdgeCGW_T(graph->edges, freeEdge);
  freeEdge->topLevelEdge = newEid;
  SetGraphEdgeStatus(graph, freeEdge, UNKNOWN_EDGE_STATUS);
#ifdef DEBUG
  if (graph->freeEdgeHead != NULLINDEX)
    assert (graph->freeEdgeHead != graph->tobeFreeEdgeHead);
  fprintf(GlobalData->stderrc,"* GetFreeGraphEdge returned edge " F_CID "\n", newEid);
#endif
  return freeEdge;
}


// Delete a graph node and all edges incident on this node
void DeleteGraphNode(GraphCGW_T *graph, NodeCGW_T *node){
  GraphEdgeIterator edges;
  EdgeCGW_T *edge;
  int i;
  VA_TYPE(PtrT) *edgeList = CreateVA_PtrT(100);



  InitGraphEdgeIterator(graph, node->id,
                        ALL_END, ALL_EDGES, GRAPH_EDGE_DEFAULT , &edges);
  while(NULL != (edge = NextGraphEdgeIterator(&edges))){
    ///AppendPtrT(edgeList, (PtrT)&edge);
    AppendPtrT(edgeList, (const void *) &edge);
  }

  for(i = 0; i < GetNumPtrTs(edgeList); i++){
    EdgeCGW_T *edge = *(EdgeCGW_T **)GetPtrT(edgeList,i);
    DeleteGraphEdge(graph, edge);
  }

  // Mark the node dead and
  node->flags.bits.isDead = TRUE;
  // Add to dead list
  // When we clean up the hashtable, we'll move it to the free list
  node->essentialEdgeA = graph->deadNodeHead;
  graph->deadNodeHead = node->id;

  // Unreference the consensus for this contig
  if(graph->type != SCAFFOLD_GRAPH){
    DeleteMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB,
                                    node->id, graph->type == CI_GRAPH);
    //    RemoveMultiAlignFromStore(graph->maStore, node->id);
  }
  DeleteVA_PtrT(edgeList);
}


void UnlinkGraphEdge(GraphCGW_T *graph, EdgeCGW_T *edge){

  ChunkInstanceT *chunkA, *chunkB;
  CIEdgeT *prevA, *prevB, *nextA, *nextB;
  // head->prevALink should be hooked to tail->nextALink
  // head->prevBLink should be hooked to tail->nextBlink
  // if head->prev?Link is NULLINDEX, we need to update the head of the list


  // First determine if this is a toplevel edge, or a rawEdge linked
  // to a top level edge
  {
    CDS_CID_t eid = GetVAIndex_EdgeCGW_T(graph->edges, edge);
    // This is a raw edge hanging off of a merged edge
    if(edge->flags.bits.isRaw && edge->topLevelEdge != eid){ 
      EdgeCGW_T *topEdge = GetGraphEdge(graph, edge->topLevelEdge);
      EdgeCGW_T *rawEdge, *prevEdge;
      AssertPtr(topEdge);
      assert(edge->flags.bits.isRaw);
      for(prevEdge = topEdge,
            rawEdge = GetGraphEdge(graph,topEdge->nextRawEdge); 
	  rawEdge != NULL; 
	  prevEdge = rawEdge,
            rawEdge = GetGraphEdge(graph,rawEdge->nextRawEdge)){
	assert(!rawEdge->flags.bits.isDeleted);
	if(rawEdge == edge){
	  prevEdge->nextRawEdge = edge->nextRawEdge;
	  edge->nextRawEdge = 0;
	  edge->topLevelEdge = eid;  // make this a top level edge
	  topEdge->edgesContributing--;
	  return;
	}
      }

      PrintGraphEdge(GlobalData->stderrc,graph,"*looking for* ",
                     edge, edge->idA );
      PrintGraphEdge(GlobalData->stderrc,graph,"*top level  * ",
                     topEdge, topEdge->idA );
      for(rawEdge = GetGraphEdge(graph,topEdge->nextRawEdge);
	  rawEdge != NULL;
	  rawEdge = GetGraphEdge(graph,rawEdge->nextRawEdge)){
	assert(!rawEdge->flags.bits.isDeleted);
	PrintGraphEdge(GlobalData->stderrc,graph,"*raw* ",
                       rawEdge, rawEdge->idA );
      }
      assert(0 /* We didn't find a raw edge we were trying to unlink */);
    }
  }  

  chunkA = GetGraphNode(graph, edge->idA);
  chunkB = GetGraphNode(graph, edge->idB);
  assert(chunkA && chunkB);
  
  /*
    {
    char * prefix = (edge->flags.bits.isRaw?" top level merged ": " top level raw ");
  
    PrintGraphEdge(GlobalData->stderrc,graph,prefix, edge, edge->idA );
    }
  */
  nextA = NULL;
  if(edge->nextALink != NULLINDEX)
    nextA = GetGraphEdge(graph, edge->nextALink);
  
  nextB = NULL;
  if(edge->nextBLink != NULLINDEX)
    nextB = GetGraphEdge(graph, edge->nextBLink);
  
  prevA = NULL;
  if(edge->prevALink != NULLINDEX)
    prevA = GetGraphEdge(graph, edge->prevALink);
  
  prevB = NULL;
  if(edge->prevBLink != NULLINDEX)
    prevB = GetGraphEdge(graph, edge->prevBLink);
  
  assert(!edge->flags.bits.isDeleted);
  assert(!prevA || !prevA->flags.bits.isDeleted);
  assert(!prevB || !prevB->flags.bits.isDeleted);
  assert(!nextA || !nextA->flags.bits.isDeleted);
  assert(!nextB || !nextB->flags.bits.isDeleted);

  // Unlink A Chain
  if(!prevA){
    chunkA->edgeHead = edge->nextALink;
  }else{
    if(prevA->idA == edge->idA){
      prevA->nextALink = edge->nextALink;
    }else{
      assert(prevA->idB == edge->idA);
      prevA->nextBLink = edge->nextALink;
    }
  }
  // Fix up the backpointers from nextA
  // This should handle the case where edge is at the head of the A List
  if(nextA){
    if(nextA->idA == edge->idA){
      nextA->prevALink = edge->prevALink;
    }else{
      assert(nextA->idA == edge->idB);
      nextA->prevBLink = edge->prevALink;
    }
  }

  // Unlink B Chain
  if(!prevB){
    chunkB->edgeHead = edge->nextBLink;
  }else{
    if(prevB->idB == edge->idB){
      prevB->nextBLink = edge->nextBLink;
    }else{
      assert(prevB->idA == edge->idB);
      prevB->nextALink = edge->nextBLink;
    }
  }
  
  if(nextB){
    if(nextB->idB == edge->idB){
      nextB->prevBLink = edge->prevBLink;
    }else{
      assert(nextB->idA == edge->idB);
      nextB->prevALink = edge->prevBLink;
    }
  }

}

// Delete a top level CIEdgeT by unlinking it from the graph - remember
// this will unlink any dangling raw edges as well - and mark it as deleted.
void  DeleteGraphEdge(GraphCGW_T *graph,  EdgeCGW_T *edge){

  assert(!  edge->flags.bits.isDeleted );

  //  PrintGraphEdge(GlobalData->stderrc,graph,"Deleting ", edge, edge->idA);

  // Unlink the edge from the graph
  UnlinkGraphEdge(graph, edge);

  //  edge->flags.bits.isDeleted = TRUE;
  edge->prevALink = edge->prevBLink = edge->nextALink = edge->nextBLink = NULLINDEX;
  FreeGraphEdge(graph,edge);
  return;
}

/* ************************************************************/
void PrintGraphEdge(FILE *fp, GraphCGW_T *graph,
                    char *label, EdgeCGW_T *edge, CDS_CID_t cid){
  char actualOverlap[256];
  CDS_COORD_t actual = 0;
  CDS_COORD_t delta = 0;
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
  }else if(ChunkInstanceA->aEndCoord > 0 && ChunkInstanceB->aEndCoord > 0){
    ChunkInstanceT *containing, *contained;
    int aContainsB = FALSE, bContainsA= FALSE;
    if(edge->flags.bits.bContainsA || (-edge->distance.mean) > ChunkInstanceA->bpLength.mean){
      bContainsA = TRUE;
    }
    if(edge->flags.bits.aContainsB || (-edge->distance.mean) > ChunkInstanceB->bpLength.mean){
      aContainsB = TRUE;
    }

    if(aContainsB || bContainsA){
      if(aContainsB){
	containing = ChunkInstanceA;
	contained = ChunkInstanceB;

	switch(edge->orient){
          case AB_AB:
            //     ------------------------>
            //            =======>----------
            //            |----------------| 
            // distance between A end of contained and B end of containing
            actual = -abs(containing->bEndCoord - contained->aEndCoord);
            break;

          case BA_BA:
            //     ------------------------>
            //     ------=======>
            //     |------------| 
            // distance between A end of containing and B end of contained
            actual = -abs(containing->aEndCoord - contained->bEndCoord);
            break;
          case AB_BA:
            //     ------------------------>
            //            <========--------
            //            |----------------| 
            // distance between B end of containing and B end of contained
            actual = -abs(containing->bEndCoord - contained->bEndCoord);
            break;
          case BA_AB:
            //     <-----------------------
            //             =======>--------
            //            |----------------| 
            // distance between A end of containing and A end of contained
            actual = -abs(containing->aEndCoord - contained->aEndCoord);
            break;
          default:
            assert(0);
            break;
	}

      }else{
	contained = ChunkInstanceA;
	containing = ChunkInstanceB;
	switch(edge->orient){
          case AB_AB:
            //     ------------------------>
            //     ------=======>
            //     |------------| 
            // distance between A end of containing and B end of contained
            actual = -abs(containing->aEndCoord - contained->bEndCoord);
            break;

          case BA_BA:
            //     ------------------------>
            //            =======>----------
            //            |----------------| 
            // distance between A end of contained and B end of containing
            actual = -abs(containing->bEndCoord - contained->aEndCoord);
            break;
          case AB_BA:
            //     ------------------------>
            //            <========--------
            //            |----------------| 
            // distance between B end of containing and B end of contained
            actual = -abs(containing->bEndCoord - contained->bEndCoord);
            break;
          case BA_AB:
            //     <-----------------------
            //             =======>--------
            //            |----------------| 
            // distance between A end of containing and A end of contained
            actual = -abs(containing->aEndCoord - contained->aEndCoord);
            break;
          default:
            assert(0);
            break;
	}
      }
    }else{
      actual = -IntervalsOverlap(ChunkInstanceA->aEndCoord,
                                 ChunkInstanceA->bEndCoord, 
				 ChunkInstanceB->aEndCoord,
                                 ChunkInstanceB->bEndCoord,-500000);
    }
    delta = edge->distance.mean - actual;
    if(actual != 0)
      sprintf(actualOverlap,"actual = " F_COORD "(" F_COORD ") %s",
              actual,delta, (edge->flags.bits.isProbablyBogus?"*PB*":""));
  }

  assert(!(edge->flags.bits.hasContributingOverlap &&
           edge->flags.bits.hasTandemOverlap));

  strcpy(flagbuf,"");
  if(edge->flags.bits.hasContributingOverlap){
    if(edge->flags.bits.isPossibleChimera)
      strcat(flagbuf,"$?");
    else 
      strcat(flagbuf,"$0");
  }else if(edge->flags.bits.hasRepeatOverlap){
    strcat(flagbuf,"$R");

  }else if(edge->flags.bits.hasTandemOverlap){
    strcat(flagbuf,"$T");
  }
  if(edge->flags.bits.hasTransChunk){
    strcat(flagbuf,"$t");
  }
#if 0
  if(edge->flags.bits.mustOverlap){
    strcat(flagbuf,"$m");
  }
#endif
  if(edge->flags.bits.aContainsB){
    strcat(flagbuf,"$C");
  }else if(edge->flags.bits.bContainsA){
    strcat(flagbuf,"$I");
  }
  if(edge->flags.bits.hasGuide){
    strcat(flagbuf,"$G");
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

  fprintf(fp,"%s eid:" F_CID " A:" F_CID " B:" F_CID " wgt:%d %s %s %s %s%s ori:%c qua:%3.2g trstd:%d con:%d dst:%d std:%g %s ",
          label, 
          eid,
          edge->idA, edge->idB, edge->edgesContributing, flagbuf,
          (edge->flags.bits.hasExtremalAFrag?"$A":"  "),
          (edge->flags.bits.hasExtremalBFrag?"$B":"  "),
          flagTrans, (edge->flags.bits.isLeastSquares?"L":" "),
          GetEdgeOrientationWRT(edge, edge->idA),
          edge->quality, edge->flags.bits.edgeStatus & TRUSTED_EDGE_STATUS,
          edge->flags.bits.hasContainmentOverlap,
          (int)edge->distance.mean,
          (edge->distance.variance > .001?sqrt(edge->distance.variance):0.0),
          actualOverlap);

  if(edge->flags.bits.isRaw && !isOverlapEdge(edge))
    fprintf(fp,"(" F_CID "," F_CID ")", edge->fragA, edge->fragB);
  fprintf(fp,"\n");
}


/* ************************************************************/
void PrintContigEdgeInScfContext(FILE *fp, GraphCGW_T *graph,
                                 char *label, EdgeCGW_T *edge,
                                 CDS_CID_t cid){
  char actualOverlap[256];
  CDS_COORD_t actual = 0;
  CDS_COORD_t delta = 0;
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
  }else if(ChunkInstanceA->aEndCoord > 0 && ChunkInstanceB->aEndCoord > 0){
    ChunkInstanceT *containing, *contained;
    int aContainsB = FALSE, bContainsA= FALSE;
    if(edge->flags.bits.bContainsA ||
       (-edge->distance.mean) > ChunkInstanceA->bpLength.mean){
      bContainsA = TRUE;
    }
    if(edge->flags.bits.aContainsB ||
       (-edge->distance.mean) > ChunkInstanceB->bpLength.mean){
      aContainsB = TRUE;
    }

    if(aContainsB || bContainsA){
      if(aContainsB){
	containing = ChunkInstanceA;
	contained = ChunkInstanceB;

	switch(edge->orient){
          case AB_AB:
            //     ------------------------>
            //            =======>----------
            //            |----------------| 
            // distance between A end of contained and B end of containing
            actual = -abs(containing->bEndCoord - contained->aEndCoord);
            break;

          case BA_BA:
            //     ------------------------>
            //     ------=======>
            //     |------------| 
            // distance between A end of containing and B end of contained
            actual = -abs(containing->aEndCoord - contained->bEndCoord);
            break;
          case AB_BA:
            //     ------------------------>
            //            <========--------
            //            |----------------| 
            // distance between B end of containing and B end of contained
            actual = -abs(containing->bEndCoord - contained->bEndCoord);
            break;
          case BA_AB:
            //     <-----------------------
            //             =======>--------
            //            |----------------| 
            // distance between A end of containing and A end of contained
            actual = -abs(containing->aEndCoord - contained->aEndCoord);
            break;
          default:
            assert(0);
            break;
	}

      }else{
	contained = ChunkInstanceA;
	containing = ChunkInstanceB;
	switch(edge->orient){
          case AB_AB:
            //     ------------------------>
            //     ------=======>
            //     |------------| 
            // distance between A end of containing and B end of contained
            actual = -abs(containing->aEndCoord - contained->bEndCoord);
            break;

          case BA_BA:
            //     ------------------------>
            //            =======>----------
            //            |----------------| 
            // distance between A end of contained and B end of containing
            actual = -abs(containing->bEndCoord - contained->aEndCoord);
            break;
          case AB_BA:
            //     ------------------------>
            //            <========--------
            //            |----------------| 
            // distance between B end of containing and B end of contained
            actual = -abs(containing->bEndCoord - contained->bEndCoord);
            break;
          case BA_AB:
            //     <-----------------------
            //             =======>--------
            //            |----------------| 
            // distance between A end of containing and A end of contained
            actual = -abs(containing->aEndCoord - contained->aEndCoord);
            break;
          default:
            assert(0);
            break;
	}
      }
    }else{
      actual = -IntervalsOverlap(ChunkInstanceA->aEndCoord,
                                 ChunkInstanceA->bEndCoord, 
				 ChunkInstanceB->aEndCoord,
                                 ChunkInstanceB->bEndCoord,-500000);
    }
    delta = edge->distance.mean - actual;
    if(actual != 0)
      sprintf(actualOverlap,"actual = " F_COORD "(" F_COORD ") %s",
              actual,delta, (edge->flags.bits.isProbablyBogus?"*PB*":""));
  }

  assert(!(edge->flags.bits.hasContributingOverlap &&
           edge->flags.bits.hasTandemOverlap));

  strcpy(flagbuf,"");
  if(edge->flags.bits.hasContributingOverlap){
    if(edge->flags.bits.isPossibleChimera)
      strcat(flagbuf,"$?");
    else 
      strcat(flagbuf,"$0");
  }else if(edge->flags.bits.hasRepeatOverlap){
    strcat(flagbuf,"$R");

  }else if(edge->flags.bits.hasTandemOverlap){
    strcat(flagbuf,"$T");
  }
  if(edge->flags.bits.hasTransChunk){
    strcat(flagbuf,"$t");
  }
#if 0
  if(edge->flags.bits.mustOverlap){
    strcat(flagbuf,"$m");
  }
#endif
  if(edge->flags.bits.aContainsB){
    strcat(flagbuf,"$C");
  }else if(edge->flags.bits.bContainsA){
    strcat(flagbuf,"$I");
  }
  if(edge->flags.bits.hasGuide){
    strcat(flagbuf,"$G");
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

  fprintf(fp,"%s eid:" F_CID " A:" F_CID " B:" F_CID " wgt:%d %s %s %s %s%s ori:%c qua:%3.2g trstd:%d con:%d dst:%d std:%g %s sameScf: %c Apos [%d,%d] Bpos [%d,%d]",
          label, 
          eid,
          edge->idA, edge->idB, edge->edgesContributing, flagbuf,
          (edge->flags.bits.hasExtremalAFrag?"$A":"  "),
          (edge->flags.bits.hasExtremalBFrag?"$B":"  "),
          flagTrans, (edge->flags.bits.isLeastSquares?"L":" "),
          GetEdgeOrientationWRT(edge, edge->idA),
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
    fprintf(fp,"(" F_CID "," F_CID ")", edge->fragA, edge->fragB);
  fprintf(fp,"\n");
}


/****************************************************************************/
// collectOverlap argument is ignored.
// We NEVER collect overlaps! SAK 10/13/2000
// This means overlaps are never pre-computed, which was needed for dros
// gap walking, but instead, discovered as needed.
CDS_CID_t AddGraphEdge(GraphCGW_T *graph,
                       CDS_CID_t cidA, CDS_CID_t cidB, 
                       CDS_CID_t fragidA, CDS_CID_t fragidB, 
                       CDS_CID_t dist,
                       LengthT distance,
                       float32 quality,
                       CDS_COORD_t fudgeDistance,
                       OrientType orientation,
                       int isInducedByUnknownOrientation,
                       int isGuide,   
                       int isSTSGuide,
                       int isMayJoin,
                       int isMustJoin,
                       int isOverlap,
                       int isRepeatOverlap,
                       int isTandemOverlap,
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
  OrientType orient;
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
    switch(orientation){
      case AB_BA:
      case BA_AB:
	orient = orientation;
	break;
      case AB_AB:
	orient = BA_BA;
	break;
      case BA_BA:
	orient = AB_AB;
	break;
      default:
	assert(0);
    }
  }
  
  CIa = GetGraphNode(graph, idA);
  CIb = GetGraphNode(graph, idB);
  
  assert(! (isTandemOverlap  && isOverlap));
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
#if 0
  ciedge->fudgeDistance = fudgeDistance;
  ciedge->minDistance = distance.mean - fudgeDistance;
  ciedge->maxDistance = distance.mean + fudgeDistance;
#endif
  ciedge->orient = orient;
  ciedge->flags.all = 0;
  ciedge->prevBLink = ciedge->prevALink = NULLINDEX;
  ciedge->nextBLink = ciedge->nextALink = NULLINDEX;
  ciedge->flags.bits.isRaw = TRUE;
  ciedge->nextRawEdge = NULLINDEX;
  ciedge->referenceEdge = NULLINDEX;
  ciedge->topLevelEdge = ciedgeIndex;
  ciedge->flags.bits.inducedByUnknownOrientation = isInducedByUnknownOrientation;
  ciedge->flags.bits.hasGuide = isGuide;
  ciedge->flags.bits.hasSTSGuide = isSTSGuide;
  ciedge->flags.bits.hasMayJoin = isMayJoin;
  ciedge->flags.bits.hasMustJoin = isMustJoin;
  ciedge->flags.bits.hasContributingOverlap = isOverlap;
  ciedge->flags.bits.hasRepeatOverlap = isRepeatOverlap;
  ciedge->flags.bits.hasTandemOverlap = isTandemOverlap;
  ciedge->flags.bits.mustOverlap = isOverlap ; /// TEMPORARY
  ciedge->flags.bits.aContainsB = aContainsB;
  ciedge->flags.bits.bContainsA = bContainsA;
  ciedge->flags.bits.hasTransChunk = isTransChunk;
  ciedge->flags.bits.hasExtremalAFrag = extremalA;
  ciedge->flags.bits.hasExtremalBFrag = extremalB;
  ciedge->flags.bits.hasContainmentOverlap = aContainsB | bContainsA;
  ciedge->flags.bits.isSloppy = isGuide ||
    (distance.variance > SLOPPY_EDGE_VARIANCE_THRESHHOLD);
  
  
  SetGraphEdgeStatus(graph, ciedge, status);
  
#if 0 // DEBUG_CIEDGES
  fprintf(GlobalData->stderrc,"* Add %c Edge (" F_CID "," F_CID ") orient:%c distance:%g +/- %g guide:%d overlap:%d tandem:%d repeat:%d extremeA:%d B:%d\n",
          graph->type,
          idA, idB,
          orient,
          distance.mean,
          sqrt(distance.variance),
          isGuide,
          isOverlap,
          isTandemOverlap,
          isRepeatOverlap,
          isExtremalA,
          isExtremalB);
  fprintf(GlobalData->stderrc,"* ... isIndByUnkOri %d isSTSGuide %d isMayJoin %d isMustJoin %d isAContainsB %d isBContainsA %d isTransChunk %d\n",
          isInducedByUnknownOrientation,
          isSTSGuide,
          isMayJoin,
          isMustJoin,
          isAContainsB,
          isBContainsA,
          isTransChunk);
#endif
  
#ifdef DEBUG_DATA
  // Determine whether this is a bogus edge generated by the simulator
  // or seems bogus otherwise.
  if(ciedge->fragA != NULLINDEX && ciedge->fragB != NULLINDEX){
    // Not an overlap fragment
    CIFragT *fragA = GetCIFragT(ScaffoldGraph->CIFrags, ciedge->fragA);
    CIFragT *fragB = GetCIFragT(ScaffoldGraph->CIFrags, ciedge->fragB);
    int isProbablyBogus = FALSE;
    int bogus = (fragA->flags.bits.hasFalseMate ||  fragB->flags.bits.hasFalseMate);
#ifdef DEBUG_CIEDGES
    CDS_COORD_t distance = MAX(fragA->aEndCoord, fragB->aEndCoord) -
      MIN(fragA->aEndCoord, fragB->aEndCoord);
#endif
    
    ciedge->flags.bits.isBogus = bogus;
    // Determine whether this LOOKS like a bogus edge.  Our criteria are:
    //  1) Implied overlap is > 500 base pairs + 5 sigma
    //  2) Implied distance is outside mean of distribution + 5 sigma
    
    {
      DistT *distRecord = GetDistT(ScaffoldGraph->Dists, dist);
      CDS_COORD_t minDistance = -500 - distRecord->stddev * 5;
      CDS_COORD_t maxDistance = distRecord->mean + distRecord->stddev * 5;
#ifdef DEBUG_CIEDGES
      ChunkInstanceT *CIA = GetGraphNode(graph, ciedge->idA);
      ChunkInstanceT *CIB = GetGraphNode(graph, ciedge->idB);
#endif
      
      isProbablyBogus = 
        ( ciedge->distance.mean < minDistance ||
          ciedge->distance.mean > maxDistance);
      if(isProbablyBogus && !ciedge->flags.bits.isBogus){
#ifdef DEBUG_CIEDGES
        fprintf(GlobalData->stderrc,"* Edge (" F_CID "," F_CID ") %c distance:%f [" F_COORD "," F_COORD "] marked probably bogus --> IS REAL?\nrawDistance " F_COORD " [" F_COORD "," F_COORD "] [" F_COORD "," F_COORD "] simChunks [" F_COORD "," F_COORD "] [" F_COORD "," F_COORD "]!\n",
                ciedge->idA, ciedge->idB, ciedge->orient,
                ciedge->distance.mean, minDistance, maxDistance,
                distance,
                fragA->aEndCoord, fragA->bEndCoord,
                fragB->aEndCoord, fragB->bEndCoord,
                CIA->aEndCoord, CIA->bEndCoord,
                CIB->aEndCoord, CIB->bEndCoord );
#endif
      }
      ciedge->flags.bits.isProbablyBogus = isProbablyBogus;
    }
  }else{
    ciedge->flags.bits.isProbablyBogus = FALSE;
  }
#else
  ciedge->flags.bits.isProbablyBogus = FALSE;
#endif
  
  assert(!(ciedge->flags.bits.hasContributingOverlap &&
           ciedge->flags.bits.hasTandemOverlap));
  /*
    AppendGraphEdge(graph, &ciedge);
    ciedgeIndex = GetNumGraphEdges(graph) - 1;
  */
  if(insert){
    //MarkTandemEdge(graph, &ciedge);
    
    //    fprintf(GlobalData->stderrc,"* Inserting graph edge " F_CID "\n",
    //	    ciedgeIndex);
    InsertGraphEdge(graph, ciedgeIndex, FALSE);
  }
  
#if 0
  fprintf(stderr,"* AddGraphEdge (" F_CID "," F_CID ",%c) %s %s\n",
          ciedge->idA, ciedge->idB, ciedge->orient,
          (CIa->flags.bits.isPotentialRock?"YES":"NO"),
          (CIb->flags.bits.isPotentialRock?"YES":"NO"));
#endif
  if(collectOverlap &&
     CIa->flags.bits.isPotentialRock &&
     CIb->flags.bits.isPotentialRock){ 
    if(!isGuide){
      if(isOverlapEdge(ciedge)){ // Add an entry to the overlapper's table
        CollectChunkOverlap(graph,
                            ciedge->idA, 
                            ciedge->idB,
                            ciedge->orient, 
                            -ciedge->distance.mean, (float64) fudgeDistance,
                            ciedge->quality,
                            ScaffoldGraph->alignOverlaps,
                            TRUE,
                            FALSE);
      }else{ //  if(ciedge->minDistance < CGW_MISSED_OVERLAP){
        //	  fprintf(stderr,"* Collected Overlap (" F_CID "," F_CID ",%c)\n", ciedge->idA, ciedge->idB, ciedge->orient);
        //	  PrintGraphEdge(GlobalData->stderrc, graph, " ", ciedge, ciedge->idA);
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
  }
  
  return ciedgeIndex;
}




void DumpGraph(GraphCGW_T *graph, FILE *stream){

  char *graphType = "";
  GraphNodeIterator nodes;
  GraphEdgeIterator edges;
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
    EdgeCGW_T *edge = NULL;
    
    fprintf(stream,"Node " F_CID "\n", node->id);
    InitGraphEdgeIterator(graph, node->id, ALL_END, ALL_EDGES,
                          GRAPH_EDGE_DEFAULT , &edges);
    while(NULL != (edge = NextGraphEdgeIterator(&edges))){
      fprintf(stream,"\t(" F_CID "," F_CID ") %c wt:%d %s\n",
              edge->idA, edge->idB, edge->orient, edge->edgesContributing,
              (isOverlapEdge(edge)?"*O*":" "));
    }
  }
}




/* Find an overlap edge between A,B */
EdgeCGW_T *FindGraphOverlapEdge(GraphCGW_T *graph,
                                CDS_CID_t idA, CDS_CID_t idB,
                                ChunkOrientationType orient){
  
  GraphEdgeIterator edges;
  EdgeCGW_T *edge;
  int end;
  int found = FALSE;

  switch(orient){
    case AB_AB:
    case AB_BA:
      end = B_END;
      break;
    case BA_AB:
    case BA_BA:
      end = A_END;
      break;
    default:
      assert(0);
  }
  
#if 0
  fprintf(GlobalData->stderrc,
          "* Looking for edges " F_CID "," F_CID " orient:%c\n",
          idA, idB,orient);
  fprintf(GlobalData->stderrc,
          "* Looking for edges of " F_CID " at end %s\n",
          idA, (end == A_END? "A_END":"B_END"));
#endif
  InitGraphEdgeIterator(graph, idA, end, ALL_EDGES,
                        GRAPH_EDGE_DEFAULT, &edges);
  
  while(NULL != (edge = NextGraphEdgeIterator(&edges))){
    ChunkOrientationType corient = GetEdgeOrientationWRT(edge, idA);
    if(!isOverlapEdge(edge)){
#if 0
      fprintf(GlobalData->stderrc,"* Found, but not an overlap edge\n");
#endif
      continue;
    }
    if(corient != orient){
#if 0
      fprintf(GlobalData->stderrc,"* Orientation %c is not %c\n",
              corient, orient);
#endif
      continue;
    }
    if(edge->idA != idB &&
       edge->idB != idB){
#if 0
      fprintf(GlobalData->stderrc,
              "* ida = " F_CID " idb = " F_CID " NOT Equal to " F_CID "\n",
              edge->idA, edge->idB, idB);
#endif
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
                         ChunkOrientationType orient){

  GraphEdgeIterator edges;
  EdgeCGW_T *edge;
  int end;
  int found = FALSE;
  
  switch(orient){
    case AB_AB:
    case AB_BA:
      end = B_END;
      break;
    case BA_AB:
    case BA_BA:
      end = A_END;
      break;
    case XX_XX:
      end = ALL_END;
      break;
    default:
      assert(0);
  }
  
#if 0
  fprintf(GlobalData->stderrc,
          "* Looking for edges " F_CID "," F_CID " orient:%c\n",
          idA, idB,orient);
  fprintf(GlobalData->stderrc,
          "* Looking for edges of " F_CID " at end %s\n",
          idA, (end == A_END? "A_END":"B_END"));
#endif
  InitGraphEdgeIterator(graph, idA, end, ALL_EDGES,
                        GRAPH_EDGE_DEFAULT, &edges);
  
  while(NULL != (edge = NextGraphEdgeIterator(&edges))){
    ChunkOrientationType corient = GetEdgeOrientationWRT(edge, idA);
    if(corient != orient){
#if 0
      fprintf(GlobalData->stderrc,"* Orientation %c is not %c\n",
              corient, orient);
#endif
      continue;
    }
    if(edge->idA != idB &&
       edge->idB != idB){
#if 0
      fprintf(GlobalData->stderrc,
              "* ida = " F_CID " idb = " F_CID " NOT Equal to " F_CID "\n",
              edge->idA, edge->idB, idB);
#endif
      continue;
    }
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
                             ChunkOrientationType orient){
  CIEdgeT *edge = FindGraphOverlapEdge(graph, idA,idB, orient);
  
  AssertPtr(edge); // If we don't find it...assert
  DeleteGraphEdge(graph, edge);
  return;
}


// Scratch space for MergeNodeGraphEdges
//
static  VA_TYPE(CDS_CID_t) *mergeCandidates = NULL;
static   VA_TYPE(CDS_CID_t) *chain = NULL;

/* MergeNodeGraphEdges
 *   Merge the edges incident on a particular node
 *
 */
void MergeNodeGraphEdges(GraphCGW_T *graph, NodeCGW_T *node,
                         int includeGuides, int mergeAll, int debug){
  CDS_CID_t id;
  EdgeCGW_T *head, *edge;
  CDS_CID_t *candidates;
  GraphEdgeIterator edges;
  id = node->id;
  
  if(mergeCandidates == NULL){
    mergeCandidates = CreateVA_CDS_CID_t(1024);
    chain = CreateVA_CDS_CID_t(1024);
  }else{
    ResetVA_CDS_CID_t(mergeCandidates);
    ResetVA_CDS_CID_t(chain);
  }
  
  InitGraphEdgeIterator(graph, id, ALL_END, ALL_EDGES, GRAPH_EDGE_DEFAULT, &edges);
  head = NULL;
  while(NULL != (edge = NextGraphEdgeIterator(&edges))){
    if(!head || 
       head->idA != edge->idA || 
       head->idB != edge->idB ||
       head->orient != edge->orient // ){
       || (!includeGuides && isSloppyEdge(edge))  /* don't merge guides */ ){
      
      if(debug){
        fprintf(GlobalData->stderrc,"* Appending edge " F_CID " since (" F_CID "," F_CID ") is different than (" F_CID "," F_CID ")\n",
		edges.curr, (head?head->idA:-1), (head?head->idB:-1), edge->idA, edge->idB);
      }
      head = edge;
      if(head->idA == id || mergeAll) // avoid double counting
        AppendCDS_CID_t(mergeCandidates, &edges.curr);
    }
  }
  
  
  // Now, iterate over the candidates
  candidates = GetCDS_CID_t(mergeCandidates, 0);
  
  if(debug){
    fprintf(GlobalData->stderrc,"* Found %d merge Candidates from node " F_CID "\n",
            (int) GetNumint32s(mergeCandidates), id);
  }
  
  for(id = 0; id < GetNumCDS_CID_ts(mergeCandidates); id++){
    int length;
    ResetCDS_CID_t(chain);
    if(debug){
      fprintf(GlobalData->stderrc,"* Checking CIEdge Chain from edge " F_CID " " F_CID "\n",
              id,candidates[id]); 
    }
    length = FindGraphEdgeChain(graph, candidates[id], chain, FALSE, includeGuides);
    
    if(length > 1){
      CIEdgeT *edge = GetGraphEdge(graph, candidates[id]);
      int edgesAdded;
      int i;
      
      if(debug){
        fprintf(GlobalData->stderrc,"* Passing the following chain of %d edges to merge\n",
                length);
        for(i = 0; i < GetNumCDS_CID_ts(chain); i ++){
          CDS_CID_t edgeID = *GetCDS_CID_t(chain,i);
          CIEdgeT *edge = GetGraphEdge(graph, edgeID);
          
          PrintGraphEdge(GlobalData->stderrc, graph, " ", edge, edge->idA);
        }
      }
      
      edgesAdded = MergeGraphEdges(graph, chain);
      if(debug && edgesAdded > 0){
        fprintf(GlobalData->stderrc,"* Found a chain of length %d between (" F_CID "," F_CID ") orient:%c ==> Added %d edges\n",
                length, edge->idA, edge->idB, edge->orient, edgesAdded); 
      }
      if(debug && edgesAdded > 0){
        for(i = 0; i < GetNumCDS_CID_ts(chain); i ++){
          CDS_CID_t edgeID = *GetCDS_CID_t(chain,i);
          CIEdgeT *edge = GetGraphEdge(graph, edgeID);
          
          fprintf(GlobalData->stderrc,"* Edge %d returned is edge " F_CID " -- %s\n",
                  i, edgeID, (edge->flags.bits.isRaw?" Raw ":" Merged "));
        }
      }
      
      for(i = 0; i < GetNumCDS_CID_ts(chain); i ++){
        CDS_CID_t edgeID = *GetCDS_CID_t(chain,i);
#if 0 // Not currently supported for ScaffoldGraphT
        TruncateEdgeRange(newEdge, overlapEdgeAll, graph);
#endif
        InsertGraphEdge(graph,edgeID,FALSE);
      }
      
    }
  }
}

// Merge all of the edges
void MergeAllGraphEdges(GraphCGW_T *graph, int includeGuides){
  GraphNodeIterator nodes;
  NodeCGW_T *CI;
  
  InitGraphNodeIterator(&nodes, graph, GRAPH_NODE_DEFAULT);
  while(NULL != (CI = NextGraphNodeIterator(&nodes))){
    
    MergeNodeGraphEdges(graph, CI, includeGuides, FALSE, FALSE);
  }
}




/* isRepeatOverlap
   Check whether an overlap penetrates into the Unique meat of a chunks, or is simply
   a repeat overlap.
*/

#define AS_CGW_MIN_REPEAT_OVERLAP 20

int IsRepeatOverlap(GraphCGW_T *graph,
                    CDS_CID_t cid1, CDS_CID_t cid2,
                    ChunkOrientationType orient, LengthT overlap){
  ChunkInstanceT *chunk1 = GetGraphNode(graph, cid1);
  ChunkInstanceT *chunk2 = GetGraphNode(graph, cid2);
  CDS_COORD_t branchPoint1, branchPoint2;
  int ret;
  assert(chunk1 && chunk1);
  
  branchPoint1 = -1;
  branchPoint2 = -1;
  
  switch(orient){
    case AB_BA:
      branchPoint1 = chunk1->info.CI.branchPointB;
      branchPoint2 = chunk2->info.CI.branchPointB;
      break;
    case AB_AB:
      branchPoint1 = chunk1->info.CI.branchPointB;
      branchPoint2 = chunk2->info.CI.branchPointA;
      break;
    case BA_BA:
      branchPoint1 = chunk1->info.CI.branchPointA;
      branchPoint2 = chunk2->info.CI.branchPointB;
      break;
    case BA_AB:
      branchPoint1 = chunk1->info.CI.branchPointA;
      branchPoint2 = chunk2->info.CI.branchPointA;
      break;
    default:
      assert(0);
  }
  
  ret = FALSE;
  if(-overlap.mean  - branchPoint1 < AS_CGW_MIN_REPEAT_OVERLAP  ||
     -overlap.mean  - branchPoint2 < AS_CGW_MIN_REPEAT_OVERLAP)
    ret = TRUE;
  
#if 0
  fprintf(GlobalData->stderrc,"* isRepeatCIOverlap:%d cid1:" F_CID " cid2:" F_CID " o:%c bp1:" F_COORD " bp2:" F_COORD " overlap:%d\n",
          ret, cid1, cid2, orient, branchPoint1, branchPoint2, (int) -overlap.mean);
#endif
  return ret;
  
}


/* CheckEdgesAgainstOverlapper
   ALL CIEdges that have a contributing overlap should have an overlap entry in the
   ScaffoldGraph overlap hashtable, and the length of the overlap should be > 0.
   
   This routine checks for this invariant, and asserts if it is not satisfied.
   
*/
void CheckEdgesAgainstOverlapper(GraphCGW_T *graph){
  EdgeCGW_T *edge;
  GraphEdgeIterator edgeMates;
  GraphNodeIterator nodes;
  NodeCGW_T *node;
  ChunkOverlapCheckT olap;
  int count    = 0;
  int failures = 0;

  fprintf(GlobalData->stderrc,"**** Calling CheckEdgesAgainstOverlapper ****\n");
  InitGraphNodeIterator(&nodes, graph, GRAPH_NODE_DEFAULT);

  while(NULL != (node = NextGraphNodeIterator(&nodes))){
    InitGraphEdgeIterator(graph, node->id, ALL_END, ALL_EDGES, GRAPH_EDGE_DEFAULT, &edgeMates);

    while(NULL != (edge = NextGraphEdgeIterator(&edgeMates))){
      if(isOverlapEdge(edge)){
        count++;
        if (LookupOverlap(graph, 
                          edge->idA, edge->idB,
                          edge->orient, &olap) == FALSE) {
          fprintf(GlobalData->stderrc,"* Failure in CheckEdgeAgainstOverlapper for edge #%d:\n", count);
          PrintGraphEdge(GlobalData->stderrc, graph, " ", edge, node->id);
          failures++;
        }
      }
    }
  }

  assert(failures == 0);

  fprintf(GlobalData->stderrc,"**** Survived CheckEdgesAgainstOverlapper with %d failures****\n", failures);
}

/****************************************************************************/
/* Update all the Unitigs belonging to the multiAlignment for contig
   contigID so that their membership and offsets are recorded properly. */

static VA_TYPE(CDS_COORD_t) *UngappedOffsets = NULL;

void UpdateNodeUnitigs(MultiAlignT *ma, ContigT *contig){
  int32 i;
  int32 numCIs = GetNumIntUnitigPoss(ma->u_list);
  NodeCGW_T *previous = NULL;
  CDS_CID_t contigID = contig->id;
  CDS_COORD_t *offsets;
  
  contig->flags.bits.isChaff = FALSE;  // We need this one in the output
  
  if(!UngappedOffsets){
    UngappedOffsets = CreateVA_CDS_COORD_t(1000);
  }
  /*
    fprintf(GlobalData->stderrc,
    "* UpdateNodeUnitigs for contig " F_CID "\n", contig->id);
  */
  
  GetMultiAlignUngappedOffsets(ma, UngappedOffsets);
  offsets = GetCDS_COORD_t(UngappedOffsets,0);
  
  for(i = 0; i < numCIs ; i++){
    IntUnitigPos *pos = GetIntUnitigPos(ma->u_list, i);
    NodeCGW_T *node = GetGraphNode(ScaffoldGraph->CIGraph, pos->ident);
    int flip = (pos->position.end < pos->position.bgn);
    CDS_COORD_t bgn, end;
    
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
    node->info.CI.contigID = contigID;
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

/****************************************************************************/
/* Update all the fragments belonging to the multiAlignment for Chunk cid
   so that their membership and offsets in their CIFragT record are recorded
   properly.
   Call this on either a CI or a Contig after it is created */
void UpdateNodeFragments(GraphCGW_T *graph, CDS_CID_t cid,
                         int markFragmentsPlaced, int markUnitigAndContig){
  MultiAlignT *ma;
  static VA_TYPE(CDS_COORD_t) *ungappedOffsets = NULL;
  NodeCGW_T *node = GetGraphNode(graph, cid);
  int isPlaced = markFragmentsPlaced || (node->scaffoldID != NULLINDEX);
  CDS_CID_t i;
  CDS_CID_t extremalA = NULLINDEX;
  CDS_CID_t extremalB = NULLINDEX;
  CDS_COORD_t minOffset = CDS_COORD_MAX;
  CDS_COORD_t maxOffset = CDS_COORD_MIN;
  
  if(ungappedOffsets == NULL)
    ungappedOffsets = CreateVA_CDS_COORD_t(10);

  /*
    fprintf(GlobalData->stderrc,
    "* UpdateNodeFragments on chunk " F_CID "\n", cid);
    fprintf(GlobalData->stderrc,
    "* UpdateNodeFragments ... isPlaced = %d\n", isPlaced);
  */
  
  //  ma = GetMultiAlignInStore(graph->maStore, cid);
  ma = LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, cid,
                                     graph->type == CI_GRAPH);
  GetMultiAlignUngappedOffsets(ma, ungappedOffsets);
  
  /* Determine extremal fragments so we can label the fragments */
  for(i = 0; i < GetNumIntMultiPoss(ma->f_list); i++){
    IntMultiPos *mp = GetIntMultiPos(ma->f_list, i);
    CDS_COORD_t end = MAX( mp->position.end, mp->position.bgn);
    CDS_COORD_t beg = MIN( mp->position.end, mp->position.bgn);

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
    CDS_CID_t fragID = (CDS_CID_t)mp->sourceInt; // GetInfoByIID(ScaffoldGraph->iidToFragIndex, mp->ident)->fragIndex;
    CIFragT *frag = GetCIFragT(ScaffoldGraph->CIFrags, fragID);
    LengthT offset3p, offset5p;
    CDS_COORD_t ubgn, uend;
    int flip = (mp->position.end < mp->position.bgn);
    CDS_COORD_t bgn, end;
    

    // the original code here called for changing from space-based to 0-base-based
    // coordinates ... but this seems incompatible with the behavior of GetMultiAlignUngappedOffsets???
    // So, let's stop adjusting the interval coordinates!
#undef CONVERT_INTERVAL_TO_BASES_WHICH_IS_A_BAD_IDEA
#ifdef CONVERT_INTERVAL_TO_BASES_WHICH_IS_A_BAD_IDEA
#define INTERVAL_TO_BASE_END_CONVERSION 1
#else
#define INTERVAL_TO_BASE_END_CONVERSION 0
#endif

    // mp->position is an interval.  We need to subtract one from
    // the upper end of the interval
    if(flip){
      bgn = mp->position.bgn - INTERVAL_TO_BASE_END_CONVERSION;
      end = mp->position.end;
    }else{
      bgn = mp->position.bgn;
      end = mp->position.end - INTERVAL_TO_BASE_END_CONVERSION;
    }

    if(bgn>=ungappedOffsets->numElements){
      fprintf(stderr,"WARNING: fragment %d falls off end of multialignment %d; fudging ...\n",
	      mp->ident,ma->id);
      bgn=ungappedOffsets->numElements-1;
    }
    ubgn = *GetCDS_COORD_t(ungappedOffsets, bgn);
    if(end>=ungappedOffsets->numElements){
      fprintf(stderr,"WARNING: fragment %d falls off end of multialignment %d; fudging ...\n",
	      mp->ident,ma->id);
      end=ungappedOffsets->numElements-1;
    }
    uend = *GetCDS_COORD_t(ungappedOffsets, end);
    
    if(ubgn == uend){
      fprintf(GlobalData->stderrc,"* Fragment " F_CID " now has ungapped length = %d (" F_COORD "," F_COORD ")...from gapped (" F_COORD "," F_COORD ")...either bad multi-alignment\n"
	      "* or a fragment fully contained within a gap in the consensus due to a bubble\n",
              mp->ident,(int)abs(ubgn-uend), ubgn, uend, bgn,end);
      if(!markUnitigAndContig){
	fprintf(GlobalData->stderrc,"* Details: fragment in CI " F_CID " [" F_COORD "," F_COORD "] CtgID " F_CID " [" F_COORD "," F_COORD "]\n",
		frag->cid,(int)(frag->offset5p.mean),(int)(frag->offset3p.mean),
		frag->contigID,(int)(frag->contigOffset5p.mean),(int)(frag->contigOffset3p.mean));
      }
      fprintf(GlobalData->stderrc,"* More details: graph node " F_CID "\n",cid);

      fprintf(GlobalData->stderrc,"* Setting length to ONE for now...\n");
      if(bgn < end){
	if(ubgn<1){
	  uend++;
	} else {
	  ubgn--;
	}
      } else {
	assert(end<bgn);
	if(uend<1){
	  ubgn++;
	} else {
	  uend--;
	}
      }
    }
    offset5p.mean = ubgn;
    offset5p.variance = ComputeFudgeVariance(offset5p.mean);
    
    offset3p.mean =  uend;
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

#if 0
    fprintf(stderr,"* result: " F_CID " frag->offset to " F_COORD "," F_COORD " and frag->contigOffset " F_COORD "," F_COORD " of " F_CID " (frag->cid " F_CID " CIid " F_CID " contigID " F_CID " )\n",
	    frag->iid,
	    (int) frag->offset5p.mean,
	    (int) frag->offset3p.mean,
	    (int) frag->contigOffset5p.mean,
	    (int) frag->contigOffset3p.mean,
	    node->id,
	    frag->cid,
	    frag->CIid,
	    frag->contigID);
#endif

    if(i == extremalA){
      frag->label = AS_INTERCHUNK_A; 
    }else if(i == extremalB){
      frag->label = AS_INTERCHUNK_B;  /*  A->B? */
    }else{
      frag->label = AS_INTRACHUNK; 
    }
    
  }
  
}





/* Compute the offset and orientation of a fragment in its chunk/contig
   orientIsOpposite == TRUE
   Offset is from 5p end of fragment to the end of the chunk in the
   direction of the 3p end of the fragment.
   orientIsOpposite == FALSE
   Offset is from 5p end of fragment to the end of the chunk in the
   direction of the 5p end of the fragment.
*/

int FragOffsetAndOrientation(CIFragT     *frag,
                             ChunkInstanceT *chunk,
                             LengthT    *chunkOffset, // output
                             FragOrient *chunkOrient, // output
                             int32 *extremal,         // output
                             int32 orientIsOpposite
                             /* TRUE if offset should be calculated from 5' towards end of
                                chunk closest to 3' end of fragment.  
                                FALSE if offset should be calculated from 5' towards end of
                                chunk closest to 5' end.
                                See comments below */
                             
                             )
{
  LengthT offset3p, offset5p;
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
  
  if(orientIsOpposite){   /* This is the case for ALL Mates and Bac Ends */
    switch(*chunkOrient){
      case A_B:
        /*
          Chunk                        ----------------------------------->
          Frag                           --------->
          ChunkOffset                    |--------------------------------|  
          ChunkOrient                         A_B
        */
        chunkOffset->mean = chunk->bpLength.mean - offset5p.mean;
        if(frag->label == AS_INTERCHUNK_B || frag->label == AS_SINGLETON)
          *extremal = TRUE;
        break;
      case B_A:
        /*
          Chunk                        ----------------------------------->
          Frag                           <---------
          ChunkOffset                  |----------|  
          ChunkOrient                     B_A
        */
        chunkOffset->mean = offset5p.mean;
        if(frag->label == AS_INTERCHUNK_A || frag->label == AS_SINGLETON)
          *extremal = TRUE;
        break;
      default:
        assert(0);
    }
  }else{
    switch(*chunkOrient){
      case A_B:
        /*
          Chunk                        ----------------------------------->
          Frag                           --------->
          ChunkOffset                  |-|  
          ChunkOrient                         A_B
        */
        chunkOffset->mean = offset5p.mean;
        if(frag->label == AS_INTERCHUNK_B || frag->label == AS_SINGLETON)
          *extremal = TRUE;
        break;
      case B_A:
        /*
          Chunk                        ----------------------------------->
          Frag                           <---------
          ChunkOffset                              |----------------------|  
          ChunkOrient                        B_A
        */
        chunkOffset->mean = chunk->bpLength.mean - offset5p.mean;
        if(frag->label == AS_INTERCHUNK_A || frag->label == AS_SINGLETON)
          *extremal = TRUE;
        break;
      default:
        assert(0);
    }
  }
  
#ifdef DEBUG_CIEDGES_1
  fprintf(GlobalData->stderrc,"* Frag " F_CID " offset(%g,%g) clength:%g orient:%c [" F_COORD "," F_COORD "]\n",
          frag->iid, frag->offset5p.mean, frag->offset3p.mean,
          chunk->bpLength.mean, *chunkOrient, 
          frag->aEndCoord, frag->bEndCoord);
#endif
  chunkOffset->variance = ComputeFudgeVariance(chunkOffset->mean);
  
  return TRUE;
  
}



int CreateGraphEdge(GraphCGW_T *graph,  
                    CIFragT *frag, 
                    CIFragT *mfrag, 
                    DistT *dist, 
                    int type, 
                    int orient, 
                    int inducedByUnknownOrientation,
                    GraphEdgeStatT *stat,
                    int buildAll){
  int isCI;
  NodeCGW_T *node, *mnode;
  CDS_CID_t fragID = GetVAIndex_CIFragT(ScaffoldGraph->CIFrags, frag);
  CDS_CID_t mfragID = GetVAIndex_CIFragT(ScaffoldGraph->CIFrags, mfrag);
  LengthT ciOffset, mciOffset, distance;
  int32 extremalA, extremalB;
  ChunkOrient ciOrient, mciOrient;
  ChunkOrientationType ciEdgeOrient;
  EdgeStatus status;
  int insert = TRUE ; // (type == AS_MATE); // BOGUS
  
#if 0
  fprintf(GlobalData->stderrc,"* CreateGraphEdge frags (" F_CID "," F_CID ") dist:%g orient %c\n",
          fragID, mfragID, dist->mean, orient);
#endif
  
  switch(graph->type){
    case CI_GRAPH:
      if(frag->cid == mfrag->cid){
        /*
          fprintf(GlobalData->stderrc,"* Frags " F_CID " and " F_CID " are both in same unitig " F_CID "!\n",
          fragID, mfragID, frag->cid);
        */
	return FALSE;
      }
      isCI = TRUE;
      node = GetGraphNode(graph, frag->cid);
      mnode = GetGraphNode(graph, mfrag->cid);
      break;
    case CONTIG_GRAPH:
      if(frag->contigID == mfrag->contigID)
	return FALSE;
      isCI = FALSE;
      node = GetGraphNode(graph, frag->contigID);
      mnode = GetGraphNode(graph, mfrag->contigID);
      break;
    default:
      assert(0);
  }
  
  // Don't add edges to chaff
  if(GlobalData->ignoreChaffUnitigs && (
                                        mnode->flags.bits.isChaff ||
                                        node->flags.bits.isChaff)){
    return FALSE;
  }
  // Don't build double
  if(!buildAll && (node->id >  mnode->id)){
    /*
      fprintf(GlobalData->stderrc,"* not building (" F_CID "," F_CID ")\n",
      node->id, mnode->id);
    */
    return TRUE; // there ARE external edges, we just aren't building them
  }
  if(type == AS_MATE ||
     type == AS_BAC_GUIDE){
    assert(mfragID == frag->mateOf);
    if(getCIFragMateStatus(frag) == MATE_NONE ){
      fprintf(GlobalData->stderrc,
              "* Chunk " F_CID " has Bogus mate " F_CID " of frag " F_CID "\n",
              frag->cid, frag->mateOf, frag->iid);
      assert(0 /* Invalid mateStatus*/);
    }
  }
  
  
  if(GlobalData->verbose)
    fprintf(GlobalData->stderrc,
            "* Found mate of frag " F_CID " (" F_CID ") in chunk:" F_CID " mate:" F_CID " (" F_CID ") chunk:" F_CID " distID:" F_CID "\n",
            fragID, frag->iid,
            frag->cid, 
            mfragID, mfrag->iid,
            mfrag->cid, 
            (CDS_CID_t) GetVAIndex_DistT(ScaffoldGraph->Dists, dist));
  
  
  
  
  if(!FragOffsetAndOrientation(frag,
                               node,
                               &ciOffset,
                               &ciOrient,
                               &extremalA,
                               (orient == AS_GKP_INNIE ||
                                orient == AS_GKP_NORMAL))){
    fprintf(GlobalData->stderrc,
            "* CIFragOffset returned false for frag " F_CID " (" F_CID ")\n",
            fragID, frag->iid);
    return FALSE;
  }
#ifdef DEBUG_CIEDGES
  if(frag->type == AS_EBAC)
    fprintf(GlobalData->stderrc,
            "* Frag " F_CID " (" F_CID ") in chunk" F_CID " at offset %g +- %g with orient %c\n",
            fragID, frag->iid, frag->cid,
            ciOffset.mean, sqrt(ciOffset.variance), ciOrient);
#endif
  if(!FragOffsetAndOrientation(mfrag,
                               mnode,
                               &mciOffset,
                               &mciOrient,
                               &extremalB,
                               (orient == AS_GKP_INNIE ||
                                orient == AS_GKP_ANTINORMAL))){
    fprintf(GlobalData->stderrc,"* CIFragOffset returned false for frag " F_CID " (" F_CID ")\n",
            mfragID, mfrag->iid);
    return FALSE;
  }
  
  
  
#ifdef DEBUG_CIEDGES
  if(mfrag->type == AS_EBAC)
    fprintf(GlobalData->stderrc,
            "* Frag " F_CID " (" F_CID ") in chunk" F_CID " at offset %g +- %g with orient %c\n",
            mfragID, mfrag->iid, mfrag->cid,
            mciOffset.mean, sqrt(mciOffset.variance), mciOrient);
#endif
  
  
  /* The following triply nested case statement captures all of the cases that arise from different
     relative alignments of the fragments in the LKG relationship, and their alignment with their 
     respective chunks.
     There is probably a better way to do this, but I think this is the clearest way to codify the relationships,
     complete with 'drawings'
  */
  
  if(frag->type == AS_EBAC){
    if(node->flags.bits.isUnique && mnode->flags.bits.isUnique)
      stat->totalUUBacPairs++;
    
    stat->totalBacPairs++;
  }
  
  switch(orient){
    
    case AS_GKP_INNIE: /********* AB_BA *******************************/
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
              ciEdgeOrient = AB_BA;  // A
              break;
            case B_A:
              //           length - 5'             gap                5'
              //      |------------------------||---------------||-----------|
              //  A --------------------------- B               A --------------------------- B
              //    5'----->                                           <------5'
              //      |-------------------------------------------------------|
              //                             mate distance
              //
              ciEdgeOrient = AB_AB; // N
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
              ciEdgeOrient = BA_BA; // I
              break;
            case B_A:
              //                     5'             gap                5'
              //      |------------------------||---------------||-----------|
              //  B --------------------------- A               A --------------------------- B
              //    5'----->                                           <------5'
              //      |-------------------------------------------------------|
              //                             mate/guide distance
              //
              ciEdgeOrient = BA_AB; // O
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
      break;
      
    case AS_GKP_NORMAL: /******* AB_AB *******************************/
      switch(ciOrient){
	//
        case A_B:
          
          switch(mciOrient){
            case B_A:
              //           length - 5'             gap              Length - 5'
              //      |------------------------||---------------||-----------|
              //  A --------------------------- B               B --------------------------- A
              //    5'----->                                                  5'------>
              //      |-------------------------------------------------------|
              //                             mate distance
              //
              ciEdgeOrient = AB_BA;  // A
              break;
            case A_B:
              //           length - 5'             gap                5'
              //      |------------------------||---------------||------|
              //  A --------------------------- B               A --------------------------- B
              //    5'----->                                           5'------>
              //      |-------------------------------------------------|
              //                             mate distance
              //
              ciEdgeOrient = AB_AB; // N
              break;
            default:
              assert(0);
              break;
          }
          break;
        case B_A:
          
          switch(mciOrient){
            case A_B:
              //                     5'             gap            5'
              //      |------------------------||---------------||----|
              //  B --------------------------- A               A --------------------------- B
              //    5'----->                                          5'------>
              //      |-----------------------------------------------|
              //                             mate distance
              //
              ciEdgeOrient = BA_AB; // O
              break;
            case B_A:
              //                     5'             gap                Length - 5'
              //      |------------------------||---------------||----|
              //  B --------------------------- A               B --------------------------- A
              //    5'----->                                          5'------>
              //      |-----------------------------------------------|
              //                             mate/guide distance
              //
              ciEdgeOrient = BA_BA; // A
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
      break;
    case AS_GKP_ANTINORMAL: /***** BA_BA *******************************/
      switch(ciOrient){
	//
        case A_B:
          
          switch(mciOrient){
            case B_A:
              //                 5'             gap                     5'
              //          |---------------------||---------------||------------------|
              //  B --------------------------- A               A --------------------------- B
              //    <-----5'                                                  <------5'
              //          |-----------------------------------------------------------|
              //                             mate distance
              //
              ciEdgeOrient = BA_AB; 
              break;
            case A_B:
              //                 5'             gap                     Length - 5'
              //          |---------------------||---------------||------------------|
              //  B --------------------------- A               B --------------------------- A
              //    <-----5'                                                  <------5'
              //          |-----------------------------------------------------------|
              //                             mate distance
              //
              ciEdgeOrient = BA_BA; // N
              break;
            default:
              assert(0);
              break;
          }
          break;
        case B_A:
          
          switch(mciOrient){
            case A_B:
              //                 Length - 5'             gap               Length - 5'
              //          |---------------------||---------------||------------------|
              //  A --------------------------- B               B --------------------------- A
              //    <-----5'                                                  <------5'
              //          |-----------------------------------------------------------|
              //                             mate distance
              //
              
              ciEdgeOrient = AB_BA; 
              break;
            case B_A:
              //                 length - 5'           gap                     5'
              //          |---------------------||---------------||------------------|
              //  A --------------------------- B               A --------------------------- B
              //    <-----5'                                                  <------5'
              //          |-----------------------------------------------------------|
              //                             mate distance
              //
              ciEdgeOrient = BA_BA; // A
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
      break;
    case AS_GKP_OUTTIE: /******** BA_AB *******************************/
      switch(ciOrient){
	//
        case A_B:
          
          switch(mciOrient){
            case B_A:
              //                 5'             gap                     Length - 5'
              //          |---------------------||---------------||------------ |
              //  B --------------------------- A               B --------------------------- A
              //    <-----5'                                                  5'------>
              //          |-----------------------------------------------------|
              //                             mate distance
              //
              ciEdgeOrient = BA_BA; 
              break;
            case A_B:
              //                 5'             gap                      5'
              //          |---------------------||---------------||-----------|
              //  B --------------------------- A               A --------------------------- B
              //    <-----5'                                                  5' ------>
              //          |---------------------------------------------------|
              //                             mate distance
              //
              ciEdgeOrient = BA_AB; // N
              break;
            default:
              assert(0);
              break;
          }
          break;
        case B_A:
          
          switch(mciOrient){
            case B_A:
              //                 Length - 5'          gap                     Length - 5'
              //          |---------------------||---------------||------------ |
              //  A --------------------------- B               B --------------------------- A
              //    <-----5'                                                  5'------>
              //          |-----------------------------------------------------|
              //                             mate distance
              //
              ciEdgeOrient = AB_BA; 
              break;
            case A_B:
              //                 Length - 5'            gap                      5'
              //          |---------------------||---------------||-----------|
              //  A --------------------------- B               A --------------------------- B
              //    <-----5'                                                  5' ------>
              //          |---------------------------------------------------|
              //                             mate distance
              //
              ciEdgeOrient = AB_AB; // N
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
      break;
    case XX_XX:
    default:
      assert(0);
  }
  
  
  
  
#ifdef DEBUG_CIEDGES_1
  if(frag->type == AS_EBAC)
    fprintf(GlobalData->stderrc,"* Distance record mean %d cioffset %g mciOffset %g\n",
            (int)dist->mean, ciOffset.mean, mciOffset.mean);
#endif
  distance.mean = (CDS_COORD_t) dist->mean - ciOffset.mean - mciOffset.mean;
  // Since the two offsets and the dist are independent we SUM their variances
  distance.variance = dist->stddev * dist->stddev +
    ciOffset.variance + mciOffset.variance;
  
#ifdef DEBUG_CIEDGES_1
  fprintf(GlobalData->stderrc,"* Adding edge (" F_CID "," F_CID ",%c) insert:%d induced by fragments (" F_CID "," F_CID ",%c) with distance %g and status %d\n",
          node->id, mnode->id, ciEdgeOrient, insert, frag->iid, mfrag->iid, (frag->flags.bits.innieMate == 1?'I':'O'),
          distance.mean, frag->flags.bits.edgeStatus);
#endif
  
  //status = frag->flags.bits.edgeStatus;
  status = AS_CGW_SafeConvert_uintToEdgeStatus(frag->flags.bits.edgeStatus);
  
  AddGraphEdge(graph, 
               node->id, 
               mnode->id, 
               fragID,
               mfragID,
               GetVAIndex_DistT(ScaffoldGraph->Dists, dist),
               distance,
               1.0,
               (int)sqrt(ciOffset.variance + mciOffset.variance), // This is used by collectOverlap as the fudge distance
               ciEdgeOrient,
               inducedByUnknownOrientation,
               type == AS_BAC_GUIDE, // isGuide
               type == AS_STS_GUIDE, // isSTSGuide
               type == AS_MAY_JOIN,  // isMayJoin
               type == AS_MUST_JOIN,  // isMustJoin
               FALSE,                        // isOverlap
               FALSE,                        // isRepeatOverlap
               FALSE,                        // isTandemOverlap
               FALSE,                        // isAContainsB
               FALSE,                        // isBContainsA
               FALSE,                        // isTransChunk
               extremalA,
               extremalB,
               status,
               graph->type == CI_GRAPH, 
               // Insert a chunk Overlap in preparation for ComputeOverlaps when we are building
               // the extended Unitig Graph.  With contigs, overlaps will be discovered as needed, and no call to
               // ComputeOverlaps is performed.
               insert); // insert only if this is a mate
  
  return TRUE;
}



// Create All raw link-based graph edges
void  BuildGraphEdgesDirectly(GraphCGW_T *graph){
  GraphEdgeStatT stat;
  GraphNodeIterator Nodes;
  NodeCGW_T *node;
  MultiAlignT *ma = CreateEmptyMultiAlignT();
  
  fprintf(stderr,"* BuildGraphEdgesDirectly\n");
  fflush(NULL);
  
  InitGraphEdgeStatT(&stat);
  
  InitGraphNodeIterator(&Nodes, graph, GRAPH_NODE_DEFAULT);
  while(NULL != (node = NextGraphNodeIterator(&Nodes))){
    if(node->flags.bits.isChaff && GlobalData->ignoreChaffUnitigs)
      continue;
    
    if((node->id % 100000) == 0){
      fprintf(stderr,"* Node " F_CID "\n", node->id);
      fflush(NULL);
    }
    ReLoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, ma, node->id, graph->type == CI_GRAPH);    // This will load the ma from the disk file
    BuildGraphEdgesFromMultiAlign(graph, node, ma, &stat, FALSE);
  }
  
  fprintf(GlobalData->stderrc,"*** BuildGraphEdgesDirectly Operated on %d fragments\n",
          stat.totalFragments);
  fprintf(GlobalData->stderrc,"\tfound %d/%d BacEnd pairs are Unique-Unique\n",
          stat.totalUUBacPairs, stat.totalBacPairs);	  
  if(stat.totalMatePairs > 0){
    fprintf(GlobalData->stderrc,"\tfound %d/%d (%g%%)mate pairs are node external\n",
            stat.totalExternalMatePairs, 
            stat.totalMatePairs, 100.* (double)stat.totalExternalMatePairs/(double)(stat.totalMatePairs));
  }
  
  fflush(NULL);
  DeleteMultiAlignT(ma);
  
}


/***************************************************************************/
// Create the raw link-based edges
void  BuildGraphEdgesFromMultiAlign(GraphCGW_T *graph, NodeCGW_T *node,
                                    MultiAlignT *ma, GraphEdgeStatT *stat,
                                    int buildAll){
  int i;
  //  MultiAlignT *ma = GetMultiAlignInStore(graph->maStore, node->id);
  //  MultiAlignT *ma = LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, node->id, graph->type == CI_GRAPH);
  int32 numFrags;
  DistT *dist = NULL;
  
  assert(ma && node && graph && ma);
  
  numFrags = GetNumIntMultiPoss(ma->f_list);
  
#if 0
  fprintf(GlobalData->stderrc,"* BuildGraphEdgesFromMultiAlign on node " F_CID " with %d frags\n",
          node->id, numFrags);
#endif
  stat->totalFragments += numFrags;
  
  for(i = 0; i < numFrags; i++){
    IntMultiPos *mp = GetIntMultiPos(ma->f_list, i);
    CIFragT *frag = GetCIFragT(ScaffoldGraph->CIFrags, (CDS_CID_t)mp->sourceInt);
    int numLinks = frag->numLinks;
    int hasExternalLinks = FALSE;
    CDS_CID_t mfragID;
    CDS_CID_t fragID = (CDS_CID_t)mp->sourceInt;
    CIFragT *mfrag;
    
    /* If this fragment has no constraints... continue */
    if(numLinks == 0){
      assert(frag->mateOf == NULLINDEX);
      //	  fprintf(GlobalData->stderrc,"* skipping frag " F_CID " (fragID " F_CID ") of CI " F_CID " since it has no links\n",
      //		  frag->iid, fragID, node->id);
      continue;
    }
    // If this fragment only has links to fragments within this contig...continue 
    if(node->flags.bits.isContig){
      if(frag->contigID != node->id){
        fprintf(GlobalData->stderrc,"* Frag " F_CID " with iid " F_CID " of Contig " F_CID " has contigOf = " F_CID "!!!\n",
                fragID,
                frag->iid, node->id, frag->contigID);
        assert(0);
      }
#define OPTIMIZE_CGW 1
#if OPTIMIZE_CGW
      if( frag->flags.bits.hasInternalOnlyContigLinks){
        //	fprintf(GlobalData->stderrc,"* Skipping frag " F_CID " of Contig " F_CID " since it has only internal contig links\n",
        //		fragID, node->id);
        continue;
      }
#endif
    }else if(node->flags.bits.isCI){
      // If this fragment only has links to fragments within this CI...continue 
      if(frag->cid != node->id){
        fprintf(GlobalData->stderrc,"* Frag " F_CID " with iid " F_CID " of cid " F_CID " has ci = " F_CID "!!!\n",
                fragID,
                frag->iid, node->id, frag->cid);
        assert(0);
      }
      
#if 1
      if( frag->flags.bits.hasInternalOnlyCILinks){
        //	  fprintf(GlobalData->stderrc,"* skipping frag " F_CID " (fragID " F_CID ") of CI " F_CID " since it has only internal CI links\n",
        //		  frag->iid, fragID, node->id);
        continue;
      }
#endif
    }
    
    if(frag->flags.bits.getLinksFromStore){
      GateKeeperLinkRecordIterator GKPLinks;
      GateKeeperLinkRecord GKPLink;
      /* If we have to, get the links from the GateKeeperLink store.  This should
         be the EXCEPTION to the rule, so we don't care too much how much it costs */
      
      CDS_CID_t linkHead = frag->linkHead ; //GetInfoByIID(ScaffoldGraph->iidToFragIndex, frag->iid)->linkHead;
      assert(linkHead != NULLINDEX);
      assert(0);
      CreateGateKeeperLinkRecordIterator(ScaffoldGraph->gkpStore.lnkStore,
                                         linkHead,
                                         fragID, &GKPLinks);
      while(NextGateKeeperLinkRecordIterator(&GKPLinks, &GKPLink))
        { /* Loop over all constraints on this fragment */
        
        
        
          assert(fragID == GKPLink.frag1 ||
                 fragID == GKPLink.frag2 );
          /* Avoid double counting. This should be the same as cid < mcid, since the fragments
             are enumerated unitig by unitig */
          if(fragID != GKPLink.frag1)
            continue;
        
          mfragID = GKPLink.frag2; // GetInfoByIID(ScaffoldGraph->iidToFragIndex, GKPLink.frag2)->fragIndex;
        
          mfrag = GetCIFragT(ScaffoldGraph->CIFrags, mfragID);
        
          if(mfrag->linkType == AS_MATE)
            stat->totalMatePairs++;
        
        
          AssertPtr(mfrag);
          if(mfrag->cid == NULLINDEX){
            fprintf(GlobalData->stderrc,"* Serious error: mate of frag " F_CID " from chunk " F_CID ", frag " F_CID " has chunk " F_CID "\n",
                    fragID,
                    frag->cid, mfragID, NULLINDEX);
            assert(0);
          }
        
          dist = GetDistT(ScaffoldGraph->Dists, GKPLink.distance);
          AssertPtr(dist);
          if(GKPLink.orientation == AS_GKP_UNKNOWN){
            /* Insert edges with all four orientations of the fragments, since we don't know which one is correct */
            hasExternalLinks |= CreateGraphEdge(graph,  frag, mfrag, dist, GKPLink.type, AS_GKP_INNIE, TRUE, stat, buildAll);
            hasExternalLinks |= CreateGraphEdge(graph,  frag, mfrag, dist, GKPLink.type, AS_GKP_OUTTIE, TRUE, stat, buildAll );
            hasExternalLinks |= CreateGraphEdge(graph,  frag, mfrag, dist, GKPLink.type, AS_GKP_NORMAL, TRUE, stat, buildAll);
            hasExternalLinks |= CreateGraphEdge(graph,  frag, mfrag, dist, GKPLink.type, AS_GKP_ANTINORMAL, TRUE, stat, buildAll);
          }else{
            hasExternalLinks |= CreateGraphEdge(graph,  frag, mfrag, dist, GKPLink.type, GKPLink.orientation, FALSE, stat, buildAll);
          }
        
        }
    }else{
      /* This is the case we've optimized for */
      mfragID = frag->mateOf;
      if(mfragID != NULLINDEX){ // this could happen if there are links, but they are rereads, for example
        mfrag = GetCIFragT(ScaffoldGraph->CIFrags, mfragID);
        if(mfrag->linkType == AS_MATE)
          stat->totalMatePairs++;
        AssertPtr(mfrag);
        dist = GetDistT(ScaffoldGraph->Dists, mfrag->dist);
        hasExternalLinks |= CreateGraphEdge(graph, frag, mfrag, dist, mfrag->linkType, 
                                            (frag->flags.bits.innieMate?AS_GKP_INNIE:AS_GKP_OUTTIE), 
                                            FALSE, stat, buildAll);
      }
      //else{
      //		fprintf(GlobalData->stderrc,"* frag " F_CID " has no mate\n",fragID);
      //      }
    }
    // If we didn't insert any links to other nodes, remember this, since we can save time
    // later
    if(!hasExternalLinks){
      //	fprintf(GlobalData->stderrc,"* Marking frag " F_CID " as internal only\n", fragID);
      if(node->flags.bits.isContig){
        frag->flags.bits.hasInternalOnlyContigLinks = TRUE;
      } else if(node->flags.bits.isCI){
        frag->flags.bits.hasInternalOnlyCILinks = TRUE;
        frag->flags.bits.hasInternalOnlyContigLinks = TRUE;
      } else {
        fprintf(GlobalData->stderrc,"* Marking frag " F_CID " as internal only -- NOT isContig or isCI??\n", fragID);
        PrintFragment(frag, fragID, GlobalData->stderrc);
        assert(0);
      }
    }else{
      if(node->flags.bits.isContig){
        if(frag->flags.bits.hasInternalOnlyContigLinks){
          fprintf(GlobalData->stderrc,"* An internal only fragment (" F_CID ") in node (" F_CID ") induced an edge\n",
                  fragID, node->id);
        }
        mfrag->flags.bits.hasInternalOnlyContigLinks = FALSE;
        frag->flags.bits.hasInternalOnlyContigLinks = FALSE;
      } else if(node->flags.bits.isCI){
        if(frag->flags.bits.hasInternalOnlyCILinks){
          fprintf(GlobalData->stderrc,"* An internal only fragment (" F_CID ") in node (" F_CID ") induced an edge\n",
                  fragID, node->id);
        }
        mfrag->flags.bits.hasInternalOnlyContigLinks = FALSE;
        mfrag->flags.bits.hasInternalOnlyCILinks = FALSE;
        frag->flags.bits.hasInternalOnlyContigLinks = FALSE;
        frag->flags.bits.hasInternalOnlyCILinks = FALSE;
        
      }
      
      stat->totalExternalMatePairs++;
    }
  }
}



#if 0
/*** Insert all of the Guide edges into the graph.  Initially they are added to the CIEdge array,
     but not inserted, so that they do not get merged
***/

void InsertGuideEdges(GraphCGW_T *graph){
  CDS_CID_t edgeID;
  
  for(edgeID = 0; edgeID < GetNumGraphEdges(graph); edgeID++){
    CIEdgeT *edge = GetGraphEdge(graph, edgeID);
    
    if(edge->flags.bits.hasGuide || 
       edge->flags.bits.hasSTSGuide ||
       edge->flags.bits.hasMayJoin || 
       edge->flags.bits.hasMustJoin )
      InsertGraphEdge(graph, edgeID ,FALSE);      
  }
}
#endif



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
  fprintf(GlobalData->stderrc,"*PropagateRawEdgeStatus (" F_CID "," F_CID ",%c) to fragments " F_CID " and " F_CID " status:%d (from %d)\n",
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



/******************************************************************************
 * MarkTandemEdge
 *   Propagate tandem repeat marks between edges<->CIs.
 *   Returns TRUE if something was marked, FALSE otherwise
 *****************************************************************************/
int  MarkTandemEdge(GraphCGW_T *graph, EdgeCGW_T *edge){
  ChunkInstanceT *chunkA = GetGraphNode(graph, edge->idA);
  ChunkInstanceT *chunkB = GetGraphNode(graph, edge->idB);
  int isTandemA = FALSE;
  int isTandemB = FALSE;
  int setOverlaps = FALSE;
  int edgePenetratesBranchPointA = FALSE;
  int edgePenetratesBranchPointB = FALSE;
  
  if(!isOverlapEdge(edge) ||
     isContainmentEdge(edge))  // Dangerous to propagate through containments
    return FALSE;
  
  switch(edge->orient){
    
    case AB_AB:
      isTandemA = chunkA->flags.bits.tandemOverlaps & BEND_TANDEM_OVERLAP;
      isTandemB = chunkB->flags.bits.tandemOverlaps & AEND_TANDEM_OVERLAP;
      if(!edge->flags.bits.hasTandemOverlap && !isContainmentEdge(edge) && chunkA->flags.bits.isUnique && chunkB->flags.bits.isUnique ){
	if(chunkA->info.CI.branchPointB > 0 && ((-edge->distance.mean) > chunkA->info.CI.branchPointB + 10)){
	  edgePenetratesBranchPointA = TRUE;
	}
	if(chunkB->info.CI.branchPointA > 0 && ((-edge->distance.mean) > chunkB->info.CI.branchPointA + 10)){
	  edgePenetratesBranchPointB = TRUE;
	}
      }
      
      break;
    case BA_AB:
      isTandemA = chunkA->flags.bits.tandemOverlaps & AEND_TANDEM_OVERLAP;
      isTandemB = chunkB->flags.bits.tandemOverlaps & AEND_TANDEM_OVERLAP;
      if(!edge->flags.bits.hasTandemOverlap && !isContainmentEdge(edge) && chunkA->flags.bits.isUnique && chunkB->flags.bits.isUnique ){
	if(chunkA->info.CI.branchPointA > 0 && ((-edge->distance.mean) > chunkA->info.CI.branchPointA + 10 )){
	  edgePenetratesBranchPointA = TRUE;
	}
	if(chunkB->info.CI.branchPointA > 0 && ((-edge->distance.mean) > chunkB->info.CI.branchPointA + 10 )){
	  edgePenetratesBranchPointB = TRUE;
	}
      }
      break;
    case AB_BA:
      isTandemA = chunkA->flags.bits.tandemOverlaps & BEND_TANDEM_OVERLAP;
      isTandemB = chunkB->flags.bits.tandemOverlaps & BEND_TANDEM_OVERLAP;
      if(!edge->flags.bits.hasTandemOverlap && !isContainmentEdge(edge) && chunkA->flags.bits.isUnique && chunkB->flags.bits.isUnique ){
	if(chunkA->info.CI.branchPointB > 0 && ((-edge->distance.mean) > chunkA->info.CI.branchPointB + 10)){
	  edgePenetratesBranchPointA = TRUE;
	}
	if(chunkB->info.CI.branchPointB > 0 && ((-edge->distance.mean) > chunkB->info.CI.branchPointB + 10)){
	  edgePenetratesBranchPointB = TRUE;
	}
      }
      break;
    case BA_BA:
      isTandemA = chunkA->flags.bits.tandemOverlaps & AEND_TANDEM_OVERLAP;
      isTandemB = chunkB->flags.bits.tandemOverlaps & BEND_TANDEM_OVERLAP;
      if(!edge->flags.bits.hasTandemOverlap && !isContainmentEdge(edge) && chunkA->flags.bits.isUnique && chunkB->flags.bits.isUnique ){
	if(chunkA->info.CI.branchPointA > 0 && ((-edge->distance.mean) > chunkA->info.CI.branchPointA + 10)){
	  edgePenetratesBranchPointA = TRUE;
	}
	if(chunkB->info.CI.branchPointB > 0 && ((-edge->distance.mean) > chunkB->info.CI.branchPointB + 10 )){
	  edgePenetratesBranchPointB = TRUE;
	}
      }
      break;
    default:
      break;
  }
  
  if(edgePenetratesBranchPointA &&
     edgePenetratesBranchPointB ){
    
    PrintGraphEdge(GlobalData->stderrc, graph, "Pentrates Branch Points ", edge, edge->idA);
  }
  // If the edge doesn't think it's tandem, and it is, mark it, and the chunks
  // If either chunks has a tandem mark on the appropriate end, they both should have it.
  if((edge->flags.bits.hasTandemOverlap == FALSE) && 
     (isTandemA || isTandemB)  &&
     !(edgePenetratesBranchPointA && edgePenetratesBranchPointB)){
    setOverlaps = TRUE;
#ifdef DEBUG_TANDEM1
    fprintf(GlobalData->stderrc,"* Edge (" F_CID "," F_CID ",%c) not marked as tandem, chunks (" F_CID ",%d) (" F_CID ",%d)\n",
            edge->idA, edge->idB, edge->orient, chunkA->id, chunkA->flags.bits.tandemOverlaps,
            chunkB->id, chunkB->flags.bits.tandemOverlaps);
#endif
  }else if ((edge->flags.bits.hasTandemOverlap == TRUE) && !(isTandemA && isTandemB)){
    setOverlaps = TRUE;
#ifdef DEBUG_TANDEM
    fprintf(GlobalData->stderrc,"* Edge (" F_CID "," F_CID ",%c) marked as tandem, chunks NOT (" F_CID ",%d) (" F_CID ",%d)\n",
            edge->idA, edge->idB, edge->orient, chunkA->id, chunkA->flags.bits.tandemOverlaps,
            chunkB->id, chunkB->flags.bits.tandemOverlaps);
#endif
  }
  
  if(setOverlaps){
    edge->distance.variance = TANDEM_OVERLAP_VARIANCE;
    edge->flags.bits.hasTandemOverlap = TRUE;    
    edge->flags.bits.hasContributingOverlap = FALSE;
    edge->flags.bits.hasRepeatOverlap = FALSE;
    edge->flags.bits.mustOverlap = FALSE;
    
    assert(!(edge->flags.bits.hasContributingOverlap && edge->flags.bits.hasTandemOverlap));
    
    switch(edge->orient){
      case AB_AB:
        chunkA->flags.bits.tandemOverlaps |= BEND_TANDEM_OVERLAP;
        chunkB->flags.bits.tandemOverlaps |= AEND_TANDEM_OVERLAP;
        break;
      case BA_AB:
        chunkA->flags.bits.tandemOverlaps |= AEND_TANDEM_OVERLAP;
        chunkB->flags.bits.tandemOverlaps |= AEND_TANDEM_OVERLAP;
        break;
      case AB_BA:
        chunkA->flags.bits.tandemOverlaps |= BEND_TANDEM_OVERLAP;
        chunkB->flags.bits.tandemOverlaps |= BEND_TANDEM_OVERLAP;
        break;
      case BA_BA:
        chunkA->flags.bits.tandemOverlaps |= AEND_TANDEM_OVERLAP;
        chunkB->flags.bits.tandemOverlaps |= BEND_TANDEM_OVERLAP;
        break;
      default:
        assert(0);
    }
    
    /*
      fprintf(GlobalData->stderrc,"* Edge (" F_CID "," F_CID ",%c) marked as tandem, chunks Marked (" F_CID ",%d) (" F_CID ",%d)\n",
      edge->idA, edge->idB, edge->orient, chunkA->id, chunkA->flags.bits.tandemOverlaps,
      chunkB->id, chunkB->flags.bits.tandemOverlaps);
    */
  }
  return(setOverlaps);
}

/****************************************************************************/
/* PropagateTandemMarks
   This should be replaced with a UnionFind algorithm to mark all of the CI ends that
   are tandem repeats, followed by one pass over the edges to mark any edge that touches
   a tandem repeat.
*/
void PropagateTandemMarks(GraphCGW_T *graph){
  int iterations;
  int todo ;
  int i;
  int numTandemEdges = 0;
  
  for(i = 0;i < GetNumGraphEdges(graph); i++){
    CIEdgeT *edge = GetGraphEdge(graph, i);
    if(edge->flags.bits.hasTandemOverlap)
      numTandemEdges++;
    //      MarkTandemHighWaterMark(graph,edge);
  }
  
  fprintf(GlobalData->stderrc,"**** Initially %d/%d tandem edges \n",
          numTandemEdges, (int) GetNumGraphEdges(graph));
  
  for(iterations = 0, todo = 1; iterations < 20 && todo>0; iterations++){
    todo = 0;
    numTandemEdges = 0;
    for(i = 0; i < GetNumGraphEdges(graph); i++){
      CIEdgeT *edge = GetGraphEdge(graph, i);
      
      if(edge->flags.bits.hasTandemOverlap)
	numTandemEdges++;
      
      if(edge->flags.bits.isMarkedForDeletion ||
	 edge->flags.bits.isDeleted)
	continue;
      
      todo += MarkTandemEdge(graph, edge);
    }
    fprintf(GlobalData->stderrc,"* Propagate: iteration %d todo = %d\n",iterations,todo);
  }
  if(todo>0){
    fprintf(GlobalData->stderrc,"* PropagateTandemMarks couldn't converge after %d iterations -- marked %d edges\n",
	    iterations, todo);
  }else{
    fprintf(GlobalData->stderrc,"* PropagateTandemMarks converged after %d iterations\n",
	    iterations);
  }
  fprintf(GlobalData->stderrc,"**** Finally %d/%d tandem edges \n", numTandemEdges, (int) GetNumGraphEdges(graph));
  
}

#if 1
/****************************************************************************/
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
#if 0
  MultiAlignT *toCIMA = LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, toID, graph->type == CI_GRAPH);
  //                        GetMultiAlignInStore(ScaffoldGraph->CIGraph->maStore, toID);
  MultiAlignT *contigMA = LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, toContig->id, graph->type == CI_GRAPH);
#endif
  // GetMultiAlignInStore(ScaffoldGraph->ContigGraph->maStore, toContig->id);
  CDS_COORD_t surrogateAOffset = toCI->offsetAEnd.mean;
  CDS_COORD_t surrogateBOffset = toCI->offsetBEnd.mean;
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
    CIFragT *frag = GetCIFragT(ScaffoldGraph->CIFrags, fragID);
    assert(frag->CIid == fromID);
    assert(frag->cid == fromID);
    frag->CIid = toID;              // Assign the fragment to the surrogate
    frag->contigID = toContig->id;  // Assign the fragment to the contig
    
    //fragPos.type = frag->type;
    fragPos.type = AS_MSG_SafeConvert_charToFragType(frag->type,TRUE);
    
    fragPos.sourceInt = fragID;
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
  
#if 0
  /* Copy IntMultiPos records from the source to destination CI, adjusting consensus sequence */
  PlaceFragmentsInMultiAlignT(toCIMA, f_list_CI);
  UpdateNodeFragments(ScaffoldGraph->CIGraph, toID, FALSE);
  
  /* Copy IntMultiPos records to destination Contig, adjusting consensus sequence */
  PlaceFragmentsInMultiAlignT(contigMA, f_list_Contig);
  UpdateNodeFragments(ScaffoldGraph->ContigGraph, toContig->id,FALSE);
  UpdateNodeUnitigs(contigMA, toContig);
  
#endif
  
  /* Do not Rebuild the Mate Edges of the target CI tor reflect the changes in fragment membership */
  
  /* Do NOT Rebuild the Mate Edges of the target CI's Contig to reflect the changes in fragment membership.
     We will rebuild all mate edges when all fragments have been placed */
  
  
}

/**************************************************************************/
/* Split an unresolved CI, moving a subset of its fragments to the new node.
   Returns the index of the new node, or NULLINDEX if failure.
   The new CI gets copies of all of the overlap edges of its parent, and inherits some of the
   mate edges, as a function of fragments selected.
   The fragments listed are marked for membership (via their CIid field) int he new element
   An empty fragments array is handled the same as fragments == NULL.
*/
CDS_CID_t SplitUnresolvedCI(GraphCGW_T *graph,
                            CDS_CID_t nodeID,
                            VA_TYPE(CDS_CID_t) *fragments){
  NodeCGW_T *newNode = CreateNewGraphNode(graph);
  NodeCGW_T *node = GetGraphNode(graph, nodeID);
  int numFrags = (fragments == NULL?0:GetNumCDS_CID_ts(fragments));
  MultiAlignT *oldMA = LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, node->id, graph->type == CI_GRAPH);
  //                   GetMultiAlignInStore(graph->maStore, node->id);
  MultiAlignT *newMA;
  int verbose = FALSE;
  
  assert(graph->type == CI_GRAPH);
  assert(node->flags.bits.isCI);
  
  /* FOR NOW...we should recalc this from the fragments */
  newNode->aEndCoord = node->aEndCoord;
  newNode->bEndCoord = node->bEndCoord;
  
  newNode->offsetAEnd = node->offsetAEnd;
  newNode->offsetBEnd = node->offsetBEnd;
  newNode->flags = node->flags;
  newNode->flags.bits.isSurrogate = TRUE;
  newNode->info = node->info;
  newNode->info.CI.numInstances = 0;
  newNode->info.CI.numFragments = numFrags;
  newNode->info.CI.baseID = node->id;
  SetNodeType(newNode, RESOLVEDREPEATCHUNK_CGW);
  
  assert(node->info.CI.numInstances >=0);
  assert(node->type == UNRESOLVEDCHUNK_CGW);
  assert(newNode->type == RESOLVEDREPEATCHUNK_CGW);
  
  /* Bookkeeping on instances derived from the base */
  
  if(node->info.CI.numInstances == 0){
    if( node->flags.bits.isChaff == TRUE){
      CIFragT *frag = GetCIFragT(ScaffoldGraph->CIFrags,  (int)GetIntMultiPos(oldMA->f_list,0)->sourceInt);
      assert(frag->flags.bits.isSingleton);
      newNode->flags.bits.isChaff = FALSE;  // Making a surrogate causes the parent & surrogate to become !chaff
      node->flags.bits.isChaff = FALSE;
      frag->flags.bits.isChaff = FALSE;
    }
    node->info.CI.instances.in_line.instance1 = newNode->id;
  }else if(node->info.CI.numInstances == 1){
    node->info.CI.instances.in_line.instance2 = newNode->id;
  }else if(node->info.CI.numInstances == 2){
    VA_TYPE(CDS_CID_t) *instances = CreateVA_CDS_CID_t(16);
    AppendCDS_CID_t(instances, &node->info.CI.instances.in_line.instance1);
    AppendCDS_CID_t(instances, &node->info.CI.instances.in_line.instance2);
    AppendCDS_CID_t(instances, &(node->id));
    node->info.CI.instances.va = instances;
  }else if(node->info.CI.numInstances > 2){
    AppendCDS_CID_t(node->info.CI.instances.va, &(node->id));
  }    
  node->info.CI.numInstances++;
  
#ifdef DEBUG
  fprintf(GlobalData->stderrc,"* Node " F_CID " now has %d instances\n",
	  node->id, node->info.CI.numInstances);
#endif
  
  /* Create Surrogate MultiAlignT */
  newMA = CloneSurrogateOfMultiAlignT(oldMA, newNode->id);
  //  SetMultiAlignInStore(graph->maStore, newNode->id, newMA);
  InsertMultiAlignTInSequenceDB(ScaffoldGraph->sequenceDB, newNode->id, graph->type == CI_GRAPH,newMA, TRUE);
  
  newNode->bpLength.mean = GetMultiAlignUngappedLength(newMA);
  newNode->bpLength.variance = ComputeFudgeVariance(newNode->bpLength.mean);
  
  if(verbose){
    fprintf(GlobalData->stderrc,"* Cloned surrogate ma of CI " F_CID " (length %d) has length (%d) %g\n",
            node->id, (int) GetMultiAlignUngappedLength(oldMA),
            (int) GetMultiAlignLength(newMA), newNode->bpLength.mean);
  }
  assert(GetMultiAlignLength(newMA) == newNode->bpLength.mean);
  
  /* Copy all of the overlap edges from the original node to the surrogate */
  
  
  {
#if 0 
    // There is no need to copy the edges on the surrogate unitig
    // they are not putput and not used for anything
    GraphEdgeIterator edges;
    EdgeCGW_T *edge;
    
    InitGraphEdgeIterator(graph, node->id, ALL_END, ALL_EDGES, GRAPH_EDGE_RAW_ONLY, &edges);
    while(edge = NextGraphEdgeIterator(&edges)){
      if(isOverlapEdge(edge)){
	EdgeCGW_T *newEdge = GetFreeGraphEdge(graph);
	CDS_CID_t eid = GetVAIndex_EdgeCGW_T(graph->edges, newEdge);
        
	CDS_CID_t otherCID = (edge->idA == node->id? edge->idB:edge->idA);
	*newEdge = *edge;
	newEdge->topLevelEdge = eid; 
        
	// Make it canonical WRT otherCID and node->id
        
	if(node->id < otherCID){
	  newEdge->idA = node->id;
	  newEdge->idB = otherCID;
	}else{
	  newEdge->idB = node->id;
	  newEdge->idA = otherCID;
	}
	InsertGraphEdge(graph, eid, FALSE);
      }
    }
#endif
    /* Assign Fragments to the surrogate */
    if(numFrags > 0){
      AssignFragsToResolvedCI(graph, nodeID, newNode->id, fragments);
      
    }
    
    if(verbose){
      fprintf(GlobalData->stderrc,"* Returning new split CI " F_CID "\n", newNode->id);
    }
    return newNode->id;
  }
}




/****************************************************************************/

/* Split an unresolved Contig.
   This involves splitting its underling CI, and creating a new Contig containing
   the split CI.
   The new CI gets copies of all of the overlap edges of its parent, and inherits some of the
   mate edges, as a function of fragments selected.
   Returns the index of the new node, or NULLINDEX if failure.
   The fragments listed are marked for membership (via their contigID field) int he new element
   An empty fragments array is handled the same as fragments == NULL.
*/
CDS_CID_t SplitUnresolvedContig(GraphCGW_T *graph,
                                CDS_CID_t nodeID,
                                VA_TYPE(CDS_CID_t) *fragments,
                                int32 copyAllOverlaps){
  NodeCGW_T *newNode = CreateNewGraphNode(graph);
  NodeCGW_T *node = GetGraphNode(graph, nodeID);
  int numFrags = (fragments == NULL?0:GetNumCDS_CID_ts(fragments));
  MultiAlignT *oldMA = LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, node->id, graph->type == CI_GRAPH); 
  //  GetMultiAlignInStore(graph->maStore, node->id);
  MultiAlignT *newCIMA;
  IntUnitigPos *u = GetIntUnitigPos(oldMA->u_list,0);
  CDS_CID_t newCIid;
  NodeCGW_T *baseCI = GetGraphNode(ScaffoldGraph->CIGraph, u->ident);
  NodeCGW_T *newCI;
  int verbose = FALSE;
  
  if(verbose){
    fprintf(GlobalData->stderrc,
            "* Splitting Contig " F_CID " (%s Copy Overlaps)\n",
            node->id, (copyAllOverlaps?"":"DON'T"));
    DumpContig(GlobalData->stderrc,ScaffoldGraph, node, FALSE);
  }
  assert(graph->type == CONTIG_GRAPH);
  assert(node->flags.bits.isContig);
  assert(node->scaffoldID == NULLINDEX);
  assert(GetNumIntUnitigPoss(oldMA->u_list) == 1);
  
  node = GetGraphNode(graph, nodeID);
  
  // Get the index of the base CI
  u = GetIntUnitigPos(oldMA->u_list,0);
  baseCI = GetGraphNode(ScaffoldGraph->CIGraph, u->ident);
  
  // Split the base CI, creating a surrogate
  newCIid = SplitUnresolvedCI(ScaffoldGraph->CIGraph, baseCI->id, fragments);
  newCI = GetGraphNode(ScaffoldGraph->CIGraph, newCIid);
  newCI->flags.bits.isSurrogate = TRUE;
  newCI->flags.bits.isStoneSurrogate = copyAllOverlaps;
  newCI->flags.bits.isWalkSurrogate = !copyAllOverlaps;

  // regenerate base CI pointer in case the CIGraph gets reallocated (MP)
  baseCI = GetGraphNode(ScaffoldGraph->CIGraph, u->ident);


  assert(newCI->type == RESOLVEDREPEATCHUNK_CGW);
  
  newCI->info.CI.contigID = newNode->id; // Set the contig id
  
  // Set up the MultiAlignT
  //  newCIMA = GetMultiAlignInStore(ScaffoldGraph->CIGraph->maStore, newCIid);
  DuplicateEntryInSequenceDB(ScaffoldGraph->sequenceDB, newCIid, TRUE, newNode->id, FALSE, TRUE);
  newCIMA = LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, newCIid, TRUE);
  
  fprintf(GlobalData->stderrc,"* Inserting split CI " F_CID " (new = " F_CID ") multiAlign for new Contig " F_CID "\n",
	  baseCI->id, newCIid, newNode->id);
  //  SetMultiAlignInStore(ScaffoldGraph->ContigGraph->maStore, newNode->id, newCIMA);
  
  // Create the new surrogate contig
  newNode->bpLength = newCI->bpLength; // length of new CI may differ from length of original due to gapped vs ungapped
  
  {
    int32 newLength = (int32) GetMultiAlignUngappedLength(newCIMA);
    assert(newNode->bpLength.mean == newLength);
    assert(newCI->bpLength.mean == newLength);
  }
  
  
  newNode->aEndCoord = node->aEndCoord;
  newNode->bEndCoord = node->bEndCoord;
  
  newNode->flags = node->flags;
  newNode->flags.bits.isSurrogate = TRUE;
  newNode->info.Contig.AEndCI = newCIid;
  newNode->info.Contig.BEndCI = newCIid;
  newNode->info.Contig.numCI = 1;
  SetNodeType(newNode, CONTIG_CGW);
  
  // Set the fragments offset and contig membership
  {
    int i;
    for(i = 0; i < numFrags; i++){
      CDS_CID_t fragID = *GetCDS_CID_t(fragments, i);
      CIFragT *frag = GetCIFragT(ScaffoldGraph->CIFrags, fragID);
      frag->contigOffset5p = frag->offset5p;
      frag->contigOffset3p = frag->offset3p;
      frag->contigID = newNode->id;
    }
  }
  
  
  
  /* Copy all of the overlap edges from the original contig to the surrogate */
  if(copyAllOverlaps)
    {
      GraphEdgeIterator edges;
      EdgeCGW_T *edge;
      EdgeCGW_T copyOfEdge;
    
      // copyOfEdge - Bug fix by Jason Miller, 5/26/06.
      // The call to GetFreeGraphEdge can realloc the graph.
      // In that case, the pointer named edge became invalid.
      // Keeping a copy of the contents solved the problem.

      //    fprintf(GlobalData->stderrc,"* Copying edges\n");
      InitGraphEdgeIterator(graph, node->id, ALL_END, ALL_EDGES, GRAPH_EDGE_RAW_ONLY, &edges);
      while(NULL != (edge = NextGraphEdgeIterator(&edges))){
        if(isOverlapEdge(edge)){
          EdgeCGW_T *newEdge;
          CDS_CID_t eid;
          CDS_CID_t otherCID = (edge->idA == node->id? edge->idB:edge->idA);
        
          copyOfEdge = *edge;
          newEdge = GetFreeGraphEdge(graph); // has side effects!
          eid = GetVAIndex_EdgeCGW_T(graph->edges, newEdge);
        
          *newEdge = copyOfEdge;
          newEdge->topLevelEdge = eid;
        
          // Make it canonical WRT otherCID and node->id
        
          if(newNode->id < otherCID){
            newEdge->idA = newNode->id;
            newEdge->idB = otherCID;
            if(node->id > otherCID){
              // If the order of the nodes is different than it used to be, flip
              // the edge orientation.  This shouldn't happen, since  the
              // surrogate ids are higher than any other node
              newEdge->orient = FlipEdgeOrient(newEdge->orient);
            }
          }else{
            newEdge->idA = otherCID;
            newEdge->idB = newNode->id;
            if(node->id < otherCID){
              // If the order of the nodes is different than it used to be, flip
              // the edge orientation.  This should occur with probability 0.5
              newEdge->orient = FlipEdgeOrient(newEdge->orient);
            }
          }
          // insert the graph edge in the graph
          InsertGraphEdge(graph, eid, FALSE);
	
          CreateChunkOverlapFromEdge(graph, newEdge, FALSE); // add a hashtable entry
        }
      }
    
    }
  if(verbose){
    fprintf(GlobalData->stderrc,"* >>> NEW SPLIT CONTIG " F_CID "\n",
            newNode->id);
    DumpContig(GlobalData->stderrc,ScaffoldGraph, newNode, FALSE);
  }
  return newNode->id;
}

#endif

/***** existsContainmentRelationship *****/

UnitigOverlapType existsContainmentRelationship(NodeCGW_T *ci,
                                                NodeCGW_T *otherCI){
  CDS_COORD_t overlap = IntervalsOverlap(otherCI->offsetAEnd.mean,
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


/****************************************************************************/


/* NOTE: This routine should make a collection of distances for each
   dist so that the bucketizing can run on the precomputed collection,
   rather than recomputing.
   
   Also, we need to iterate to deal with the case that the original
   estimate is bad.
   Low priority.
*/

//#include "obsolete/computematepairstatistics"




int compareInt (const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}







void ComputeMatePairDetailedStatus() {

  GraphCGW_T *graph = ScaffoldGraph->CIGraph;
  GraphNodeIterator nodes;
  NodeCGW_T *node;
  DistT *dptr;
  MultiAlignT *ma = CreateEmptyMultiAlignT();
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
  int numDeadNode   = 0;
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

  UIDHashTable_AS *surrHash = CreateUIDHashTable_AS( 262144 );

  InitGraphNodeIterator(&nodes, graph, GRAPH_NODE_DEFAULT);
  
  while(NULL != (node = NextGraphNodeIterator(&nodes)))
    {
      switch (node->type)
        {
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

      ReLoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, ma, node->id, graph->type == CI_GRAPH);
      if (node->info.CI.baseID != node->id ) {
        fprintf(GlobalData->stderrc,"* skip surr chunk id %d, baseID %d, type %d, isStone %d\n",node->id,node->info.CI.baseID, node->type, node->flags.bits.isStoneSurrogate);
        //continue;
      }
      int numFrags  = GetNumIntMultiPoss(ma->f_list);
      if (node->flags.bits.isDead) {
        numDeadNode+=numFrags;
        continue;
      }
      numTotalFrags += numFrags;
      int i;
      for( i = 0; i < numFrags; i++)
        {
          IntMultiPos *mp = GetIntMultiPos(ma->f_list, i);
          CIFragT *frag, *mate;
          CDS_COORD_t dist;
      
          frag = GetCIFragT(ScaffoldGraph->CIFrags, (CDS_CID_t)mp->sourceInt);
          assert(frag->iid == mp->ident);
          if (frag->numLinks == 0) {
            numNoMate++;
            frag->flags.bits.mateDetail = NO_MATE;
            continue;
          }
          mate = GetCIFragT(ScaffoldGraph->CIFrags,frag->mateOf);
          if (mate == NULL) {
            numNoMate++;
            frag->flags.bits.mateDetail = NO_MATE;
            continue;
          }
          if ( mate->mateOf != (CDS_CID_t)mp->sourceInt) {
            assert(0);
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

          dptr = GetDistT(ScaffoldGraph->Dists, frag->dist);

          NodeCGW_T *fragContig, *mateContig;
          fragContig = GetGraphNode( ScaffoldGraph->ContigGraph, frag->contigID);
          AssertPtr(fragContig);
          if (fragContig->flags.bits.isDead) {
            numDeadCtg++;
          }
          switch (fragContig->type)
            {
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
            if ( node->info.CI.numInstances > 1 || node->flags.bits.isStoneSurrogate ) {
              int *frgAlready = LookupInUID2IIDHashTable_AS( surrHash, frag->iid );
              int *mateAlready = LookupInUID2IIDHashTable_AS( surrHash, mate->iid );
              if (frgAlready != NULL && mateAlready != NULL) {
                if (mate->flags.bits.mateDetail != BOTH_SURR_MATE) {
                  numSurrogate-=1;
                  numBothSurr+=2;
                  mate->flags.bits.mateDetail = BOTH_SURR_MATE;
                  frag->flags.bits.mateDetail = BOTH_SURR_MATE;
                  /*                      fprintf(GlobalData->stderrc,
                                          "* Both surr %d,%d from chunks %d.\n",
                                          frag->iid, mate->iid, node->id, mchunk->id); 
                  */
                }
                *frgAlready++;
                *mateAlready++;
              } else {
                numSurrogate++;
                mate->flags.bits.mateDetail = SURR_MATE;
                frag->flags.bits.mateDetail = SURR_MATE;
                if (frgAlready == NULL) {
                  /*                      fprintf(GlobalData->stderrc,
                                          "* Add frag %d from repeat chunk %d to surr hash, num: %d\n",
                                          frag->iid, node->id, node->info.CI.numInstances );
                  */
                  InsertInUID2IIDHashTable_AS( surrHash, frag->iid, 1);
                }
              }
              continue;
            } else {
              //              fprintf(GlobalData->stderrc,"* Non surrogate node %d, num %d, frag %d\n",
              //                      node->id, node->info.CI.numInstances, frag->iid);
            }
          }
          if( mchunk->info.CI.numInstances > 1 || mchunk->flags.bits.isStoneSurrogate ) {
            int *frgAlready = LookupInUID2IIDHashTable_AS( surrHash, frag->iid );
            if (frgAlready != NULL) { // already seen
              if (frag->flags.bits.mateDetail != BOTH_SURR_MATE) {
                numSurrogate-=1;
                numBothSurr+=2;
                mate->flags.bits.mateDetail = BOTH_SURR_MATE;
                frag->flags.bits.mateDetail = BOTH_SURR_MATE;
                InsertInUID2IIDHashTable_AS( surrHash, mate->iid, 1);
              }
              *frgAlready++;
              continue;
            } 
            int *mateAlready = LookupInUID2IIDHashTable_AS( surrHash, mate->iid );
            if (mateAlready != NULL) { // already seen
              *mateAlready++;
            } else {
              /*              fprintf(GlobalData->stderrc,
                              "* Add frag %d from chunk %d to surr hash, num: %d\n",
                              mate->iid, mchunk->id, mchunk->info.CI.numInstances );
              */
              InsertInUID2IIDHashTable_AS( surrHash, mate->iid, 1);
              numSurrogate++;
              mate->flags.bits.mateDetail = SURR_MATE;
              frag->flags.bits.mateDetail = SURR_MATE;
            }
            continue;
          }
          if ( node->type == UNRESOLVEDCHUNK_CGW)
            numInUnresolv++;
          if ( fragContig->scaffoldID == NULLINDEX )
            numInNullScaf++;
          else if ( fragContig->scaffoldID != mateContig->scaffoldID )
            // we want them to be in the same scaffold
            {
              numDiffScaf++;
              mate->flags.bits.mateDetail = DIFF_SCAFF_MATE;
              frag->flags.bits.mateDetail = DIFF_SCAFF_MATE;
              continue;
            }

          if( fragContig->scaffoldID == NULLINDEX ) {
            assert( node->type != DISCRIMINATORUNIQUECHUNK_CGW );
            if (mateContig->scaffoldID == NULLINDEX) {
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

          CDS_COORD_t fragLeftEnd, fragRightEnd;
          CDS_COORD_t mateLeftEnd, mateRightEnd;
          int mateScaffoldOrientation, fragScaffoldOrientation;

          GetFragmentPositionInScaffold( frag, &fragLeftEnd, &fragRightEnd, &fragScaffoldOrientation);
          GetFragmentPositionInScaffold( mate, &mateLeftEnd, &mateRightEnd, &mateScaffoldOrientation);

          if (fragScaffoldOrientation == mateScaffoldOrientation) 
            {
              numSame++;
              frag->flags.bits.mateDetail = SAME_ORIENT_MATE;
              mate->flags.bits.mateDetail = SAME_ORIENT_MATE;
              continue;
            }

          if (fragScaffoldOrientation == 0) // frag ---->  <---- mate
            {
              if ( fragLeftEnd < mateRightEnd ) //innie
                dist = mateRightEnd - fragLeftEnd;		  
              else { //outtie if they're psat each other
                numOuttie++;
                frag->flags.bits.mateDetail = OUTTIE_ORIENT_MATE;
                mate->flags.bits.mateDetail = OUTTIE_ORIENT_MATE;
                continue; // real outtie mates not currently supported
              }
            } else {                          // mate ---->  <---- frag
              numReverse++;
              if ( mateLeftEnd < fragRightEnd ) //innie
                dist = fragRightEnd - mateLeftEnd;		  
              else { //outtie
                numOuttie++;
                frag->flags.bits.mateDetail = OUTTIE_ORIENT_MATE;
                mate->flags.bits.mateDetail = OUTTIE_ORIENT_MATE;
                continue; // real outtie mates not currently supported
              }
            }

          if (dist < dptr->lower )
            {
              numShort++;
              frag->flags.bits.mateDetail = BAD_SHORT_MATE;
              mate->flags.bits.mateDetail = BAD_SHORT_MATE;
            }
          else if ( dist > dptr->upper)
            {
              numLong++;
              frag->flags.bits.mateDetail = BAD_LONG_MATE;
              mate->flags.bits.mateDetail = BAD_LONG_MATE;
            } else
              {
                numGood++;
                frag->flags.bits.mateDetail = GOOD_MATE;
                mate->flags.bits.mateDetail = GOOD_MATE;
              }
        }
    }
  fprintf(GlobalData->stderrc,"\n* Mate counts from ComputeMatePairDetailedStatus()\n");
  fprintf(GlobalData->stderrc,"* num Frags %d\n",numTotalFrags);
  fprintf(GlobalData->stderrc,"* num reverse frags %d\n",numReverse);
  fprintf(GlobalData->stderrc,"* num frags in dead chunks %d\n",numDeadNode);
  fprintf(GlobalData->stderrc,"* num frags in dead ctgs %d\n",numDeadCtg);
  fprintf(GlobalData->stderrc,"* num frags in unresolved chunks %d\n",numInUnresolv);
  fprintf(GlobalData->stderrc,"* num frags in repeat chunks %d\n",numInResRep);
  fprintf(GlobalData->stderrc,"* num frags in zero instance chunks %d\n",numZero);
  fprintf(GlobalData->stderrc,"* num frags in one instance chunks %d\n",numOne);
  fprintf(GlobalData->stderrc,"* num frags in >=two instance chunks %d\n",numMore);
  fprintf(GlobalData->stderrc,"* num frags in NULL scafs %d\n",numInNullScaf);
  fprintf(GlobalData->stderrc,"* num no mates %d\n",numNoMate);
  fprintf(GlobalData->stderrc,"* num good mates %d\n",numGood);
  fprintf(GlobalData->stderrc,"* num bad short mates %d\n",numShort);
  fprintf(GlobalData->stderrc,"* num bad long mates %d\n",numLong);
  fprintf(GlobalData->stderrc,"* num same orientation mates %d\n",numSame);
  fprintf(GlobalData->stderrc,"* num outtie mates %d\n",numOuttie);
  fprintf(GlobalData->stderrc,"* num both chaff mates %d\n",numBothChaff);
  fprintf(GlobalData->stderrc,"* num chaff mates %d\n",numChaff);
  fprintf(GlobalData->stderrc,"* num skiped chaff %d\n",numSkipChaff);
  fprintf(GlobalData->stderrc,"* num both degen mates %d\n",numBothDegen);
  fprintf(GlobalData->stderrc,"* num degen mates %d\n",numDegen);
  fprintf(GlobalData->stderrc,"* num skiped degen %d\n",numSkipDegen);
  fprintf(GlobalData->stderrc,"* num both surrogate mates %d\n",numBothSurr);
  fprintf(GlobalData->stderrc,"* num surrogate mates %d\n",numSurrogate);
  fprintf(GlobalData->stderrc,"* num other scaffold %d\n",numDiffScaf);
  int sum = numNoMate + numGood + numShort + numLong + numSame + numOuttie +
    numBothChaff + numChaff + numBothSurr + numSurrogate + numBothDegen +
    numDegen + numDiffScaf;
  fprintf(GlobalData->stderrc,"* sum of frag mate status %d\n\n",sum);
  //  assert( sum == numTotalFrags );
  fprintf(GlobalData->stderrc,"* Counts of top level node type\n");
  fprintf(GlobalData->stderrc,"* num DISCRIMINATORUNIQUECHUNK_CGW %d\n", duCI);
  fprintf(GlobalData->stderrc,"* num UNRESOLVEDCHUNK_CGW          %d\n", urCI);
  fprintf(GlobalData->stderrc,"* num UNIQUECHUNK_CGW              %d\n", uCI);
  fprintf(GlobalData->stderrc,"* num RESOLVEDREPEATCHUNK_CGW      %d\n", rrCI);
  fprintf(GlobalData->stderrc,"* num CONTIG_CGW                   %d\n", ctig);
  fprintf(GlobalData->stderrc,"* num UNIQUECONTIG_CGW             %d\n", uctig);
  fprintf(GlobalData->stderrc,"* num RESOLVEDCONTIG_CGW           %d\n", rctig);
  fprintf(GlobalData->stderrc,"* num UNRESOLVEDCONTIG_CGW         %d\n", urctig);
  fprintf(GlobalData->stderrc,"* num REAL_SCAFFOLD                %d\n", rScaf);
  fprintf(GlobalData->stderrc,"* num OUTPUT_SCAFFOLD              %d\n", oScaf);
  fprintf(GlobalData->stderrc,"* num SCRATCH_SCAFFOLD             %d\n", sScaf);
  fprintf(GlobalData->stderrc,"* Counts of contig level node type\n");
  fprintf(GlobalData->stderrc,"* num DISCRIMINATORUNIQUECHUNK_CGW %d\n", fduCI);
  fprintf(GlobalData->stderrc,"* num UNRESOLVEDCHUNK_CGW          %d\n", furCI);
  fprintf(GlobalData->stderrc,"* num UNIQUECHUNK_CGW              %d\n", fuCI);
  fprintf(GlobalData->stderrc,"* num RESOLVEDREPEATCHUNK_CGW      %d\n", frrCI);
  fprintf(GlobalData->stderrc,"* num CONTIG_CGW                   %d\n", fctig);
  fprintf(GlobalData->stderrc,"* num UNIQUECONTIG_CGW             %d\n", fuctig);
  fprintf(GlobalData->stderrc,"* num RESOLVEDCONTIG_CGW           %d\n", frctig);
  fprintf(GlobalData->stderrc,"* num UNRESOLVEDCONTIG_CGW         %d\n", furctig);
  fprintf(GlobalData->stderrc,"* num REAL_SCAFFOLD                %d\n", frScaf);
  fprintf(GlobalData->stderrc,"* num OUTPUT_SCAFFOLD              %d\n", foScaf);
  fprintf(GlobalData->stderrc,"* num SCRATCH_SCAFFOLD             %d\n\n", fsScaf);

  DeleteUIDHashTable_AS( surrHash );
}

/******************************************************************************
 *****************************************************************************
 *****************************************************************************/


void ComputeMatePairStatisticsRestricted( int operateOnNodes,
                                          int32 minSamplesForOverride,
                                          char *instance_label)
{
  GraphCGW_T *graph = NULL;
  GraphNodeIterator nodes;
  NodeCGW_T *node;
  DistT *dptr;
  int i, ii, maxSamples = 1;
  MultiAlignT *ma = CreateEmptyMultiAlignT();
  int numPotentialRocks = 0;
  int numPotentialStones = 0;
  //VA_TYPE(int32) *dptrFrags[GetNumDistTs(ScaffoldGraph->Dists)];
  //VA_TYPE(int32) *dptrMates[GetNumDistTs(ScaffoldGraph->Dists)];
  int NN = GetNumDistTs(ScaffoldGraph->Dists);
  VA_TYPE(CDS_CID_t) *dptrFrags[NN];
  VA_TYPE(CDS_CID_t) *dptrMates[NN];
  
  // make the stat directory if it doesn't exist already
  {
    char filename[128];
    DIR *dout;
    
    sprintf( filename, "stat");
    fprintf( stderr, "opening directory: %s\n", filename);
    dout = opendir(filename);
    if ( dout == NULL )
      {
        // the original system call has had problems, thus the checking now
        int retval;
        DIR * dirStream;
      
        fprintf( stderr,
                 "could not open directory: %s, creating...\n", filename);
        retval = system( "mkdir stat" );
        assert(retval == 0);

        dirStream = opendir("stat");
        if(dirStream == NULL)
          {
            fprintf(stderr, "Failed to access newly created directory: %s\n",
                    filename);
            assert(0);
          }
        closedir(dirStream);
      }
    else
      closedir(dout);
  }
  
  if (operateOnNodes == UNITIG_OPERATIONS)
    graph = ScaffoldGraph->CIGraph;
  else if (operateOnNodes == CONTIG_OPERATIONS)
    graph = ScaffoldGraph->ContigGraph;
  if (operateOnNodes == SCAFFOLD_OPERATIONS)
    graph = ScaffoldGraph->CIGraph;
  
  // Initialize computed mean/variance and counters, and bounds
  for(i = 1; i < GetNumDistTs(ScaffoldGraph->Dists); i++)
    {
      DistT *dptr = GetDistT(ScaffoldGraph->Dists,i);
      dptr->mu = 0.0;
      dptr->sigma = 0.0;
      dptr->numSamples = 0.0;
      dptr->numReferences = 0;
      dptr->min = CDS_COORD_MAX;
      dptr->max = CDS_COORD_MIN;
      dptr->bsize = 0;
      dptr->numBad = 0;
      dptr->histogram = NULL;
      // Lower and upper are based on nominal mean and stddev
      dptr->lower = dptr->mean - CGW_CUTOFF*dptr->stddev;
      dptr->upper = dptr->mean + CGW_CUTOFF*dptr->stddev;
    
      if(dptr->samples) // starting from scratch
        ResetVA_CDS_COORD_t(dptr->samples);
      else  // if are starting from a checkpoint?
        dptr->samples = CreateVA_CDS_COORD_t(16);
    
#if 1
      dptrFrags[i] = CreateVA_CDS_CID_t(256);
      dptrMates[i] = CreateVA_CDS_CID_t(256);
#endif
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
  
  while(NULL != (node = NextGraphNodeIterator(&nodes)))
    {
      CDS_CID_t i;
      int numFrags;
      int numExternalLinks = 0;
    
      if(node->id % 100000 == 0)
        {
          fprintf(GlobalData->stderrc, "node " F_CID "\n", node->id);
          fflush(GlobalData->stderrc);
        }
      //	fprintf(GlobalData->stderrc,"* Marking mates of %d frags in %s " F_CID "\n",	
      // numFrags, (operateOnNodes == CONTIG_OPERATIONS ? "Contig" : "Unitig"), node->id);
    
      // Don't waste time loading singletons for this
      if(node->flags.bits.isChaff) {
        numChaff++;
        continue;
      }
    
      ReLoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, ma, node->id, graph->type == CI_GRAPH);
      numFrags  = GetNumIntMultiPoss(ma->f_list);
      numTotalFrags += numFrags;
    
      if (numFrags < 2) {
        numSingle++;
        continue;
      }
    
      for( i = 0; i < numFrags; i++)
        {
          IntMultiPos *mp = GetIntMultiPos(ma->f_list, i);
          CIFragT *frag, *mate;
          CDS_COORD_t dist;
      
          frag = GetCIFragT(ScaffoldGraph->CIFrags, (CDS_CID_t)mp->sourceInt);
          assert(frag->iid == mp->ident);

#if 0
          fprintf(GlobalData->stderrc,"* frag " F_CID " (" F_CID ") mate:" F_CID " numLinks:%d extCI:%d extCon:%d\n",
                  (CDS_CID_t)mp->source, frag->iid, frag->mateOf, frag->numLinks, 
                  frag->flags.bits.hasInternalOnlyCILinks,
                  frag->flags.bits.hasInternalOnlyContigLinks);
#endif
          // This is important for keeping our computation as local as possible.
          // We skip fragments that have external links only, or no links
          if (frag->numLinks == 0) {
            numNolink++;
            continue;
          }

          if (frag->numLinks == 1 &&  // the typical case
              (operateOnNodes == CONTIG_OPERATIONS && !frag->flags.bits.hasInternalOnlyContigLinks))
            {
              //     ||    (!operateOnContigs && !frag->flags.bits.hasInternalOnlyCILinks))
              numCtgNotInternal++;
              continue;
            }
      
          mate = GetCIFragT(ScaffoldGraph->CIFrags,frag->mateOf);
      
          if (operateOnNodes == UNITIG_OPERATIONS && mate != NULL && (mate->cid != frag->cid))
            numExternalLinks++;

          if (mate == NULL) {
            numNullMate++;
            continue;
          }
          if ( mate->mateOf != (CDS_CID_t)mp->sourceInt) {
            numMateNotSource++;
            continue;
          }
          dptr = GetDistT(ScaffoldGraph->Dists, frag->dist);
          dptr->numReferences++;
      
      
          if (operateOnNodes == UNITIG_OPERATIONS)
            {
              NodeCGW_T *unitig = GetGraphNode( ScaffoldGraph->CIGraph, frag->cid);

              if (frag->cid != mate->cid) {
                numOtherUtg++;
                continue;
              }
        
              if(	getCIFragOrient(mate) == getCIFragOrient(frag)) 
                {
                  //      fprintf(GlobalData->stderrc,"* (" F_CID "," F_CID ") is bad due to orientation problems\n",      frag->iid, mate->iid);
                  frag->flags.bits.edgeStatus = UNTRUSTED_EDGE_STATUS;
                  mate->flags.bits.edgeStatus = UNTRUSTED_EDGE_STATUS;
                  dptr->numBad++;
                  numSameOrient++;
                  continue;
                }

              // now make sure the 5p end is less than the 3p end
              if ( frag->offset5p.mean > frag->offset3p.mean) {
                numReverse++;
                continue;
              }

              // try to sample fairly by only doing mates where ones of any length could live
              if ( frag->offset5p.mean + dptr->mean + 5 * dptr->stddev > unitig->bpLength.mean)
                {
                  numNot5stddev++;
                  continue;
                }

              dist = mate->offset5p.mean - frag->offset5p.mean;
        
              if((frag->flags.bits.innieMate && getCIFragOrient(frag) == B_A) ||
                 (!frag->flags.bits.innieMate && getCIFragOrient(frag) == A_B) )
                dist = -dist;
            }
          else if (operateOnNodes == CONTIG_OPERATIONS)
            {
              ContigT *contig = GetGraphNode( ScaffoldGraph->ContigGraph, frag->contigID);
        
              assert(frag->contigID == mate->contigID);
              if(GetContigFragOrient(mate) == GetContigFragOrient(frag)) 
                {
                  //  fprintf(GlobalData->stderrc,"* (" F_CID "," F_CID ") is bad due to orientation problems\n",      frag->iid, mate->iid);
                  frag->flags.bits.edgeStatus = UNTRUSTED_EDGE_STATUS;
                  mate->flags.bits.edgeStatus = UNTRUSTED_EDGE_STATUS;
                  dptr->numBad++;
                  numSameOrient++;
                  continue;
                }
        
              // now make sure the 5p end is less than the 3p end
              if ( frag->contigOffset5p.mean > frag->contigOffset3p.mean) {
                numReverse++;
                continue;
              }

              // try to sample fairly by only doing mates where ones of any length could live
              if ( frag->contigOffset5p.mean + dptr->mean + 5 * dptr->stddev > contig->bpLength.mean)
                {
                  numNot5stddev++;
                  continue;
                }
        
              dist =  mate->contigOffset5p.mean - frag->contigOffset5p.mean; 
              //   -------------------->          <----------------------
              //     mate                 innie            frag 
        
              //   <-------------------           ---------------------->
              //     mate                 outie            frag 
        
              if((frag->flags.bits.innieMate && GetContigFragOrient(frag) == B_A) ||
                 (!frag->flags.bits.innieMate && GetContigFragOrient(frag) == A_B) )
                dist =  -dist;
        
            }
          else if (operateOnNodes == SCAFFOLD_OPERATIONS)
            {
              NodeCGW_T *fragContig, *mateContig;
              CDS_COORD_t fragLeftEnd, fragRightEnd;
              CDS_COORD_t mateLeftEnd, mateRightEnd;
              int mateScaffoldOrientation, fragScaffoldOrientation;
        
              fragContig = GetGraphNode( ScaffoldGraph->ContigGraph, frag->contigID);
              AssertPtr(fragContig);
        
              mateContig = GetGraphNode( ScaffoldGraph->ContigGraph, mate->contigID);
              AssertPtr(mateContig);
        
              // we want them to be in the same scaffold
              if ( fragContig->scaffoldID != mateContig->scaffoldID || fragContig->scaffoldID == -1)
                {
                  numDiffScaf++;
                  continue;
                }
        
              GetFragmentPositionInScaffold( frag, &fragLeftEnd, &fragRightEnd, &fragScaffoldOrientation);
              GetFragmentPositionInScaffold( mate, &mateLeftEnd, &mateRightEnd, &mateScaffoldOrientation);

              if (fragScaffoldOrientation == mateScaffoldOrientation) 
                {
                  // fprintf(GlobalData->stderrc,"* (" F_CID "," F_CID ") is bad due to orientation problems\n",      frag->iid, mate->iid);
                  frag->flags.bits.edgeStatus = UNTRUSTED_EDGE_STATUS;
                  mate->flags.bits.edgeStatus = UNTRUSTED_EDGE_STATUS;
                  dptr->numBad++;
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
                CDS_COORD_t contigLeftEnd, contigRightEnd;
                int contigScaffoldOrientation;
          
                // grab scaffold
                scaff = GetGraphNode( ScaffoldGraph->ScaffoldGraph, fragContig->scaffoldID);
                extremeContig = GetGraphNode( ScaffoldGraph->ContigGraph, scaff->info.Scaffold.BEndCI);
                GetContigPositionInScaffold ( extremeContig, &contigLeftEnd, &contigRightEnd, &contigScaffoldOrientation);
          
                if ( fragLeftEnd + dptr->mean + 5 * dptr->stddev > contigRightEnd)
                  {
                    numNot5stddev++;
                    continue;
                  }
              }
        
#if 0
              if (frag->contigID != mate->contigID && fragContig->scaffoldID != -1
                  && mateContig->scaffoldID != -1
                  && mateContig->scaffoldID == fragContig->scaffoldID)
                fprintf( stderr, "frags (" F_CID "," F_CID ") are in contigs (" F_CID "," F_CID ") in scaffold " F_CID " at (" F_COORD ", " F_COORD ", %d) and (" F_COORD ", " F_COORD ", %d)\n", 
                         frag->iid, mate->iid, frag->contigID, mate->contigID,
                         fragContig->scaffoldID,
                         fragLeftEnd, fragRightEnd, fragScaffoldOrientation,
                         mateLeftEnd, mateRightEnd, mateScaffoldOrientation);
#endif
        
              if (frag->flags.bits.innieMate)  
                {
                  if (fragScaffoldOrientation == 0) // frag ---->  <---- mate
                    dist = mateRightEnd - fragLeftEnd;		  
                  else                              // mate ---->  <---- frag
                    dist = fragRightEnd - mateLeftEnd;
                }
              else  // outtie pair
                {
                  if (fragScaffoldOrientation == 0) // mate <----  ----> frag
                    dist = fragRightEnd - mateLeftEnd;		  
                  else                              // frag <----  ----> mate
                    dist = mateRightEnd - fragLeftEnd;
                }
        
              if (dist < 0)
                fprintf( stderr, "frag, mate: " F_CID ", " F_CID " have negative dist: " F_COORD "\n",
                         frag->iid, mate->iid, dist);
            }
      
          if (dist > 0 && dist < dptr->min)
            dptr->min = dist;
          if (dist > dptr->max)
            dptr->max = dist;
      
          AppendCDS_COORD_t( dptr->samples, &dist);
          AppendCDS_CID_t( dptrFrags[frag->dist], &frag->iid);
          AppendCDS_CID_t( dptrMates[frag->dist], &mate->iid);
      
          // fprintf( stderr, "adding (" F_CID ", " F_CID ") to dist " F_COORD " samples\n", frag->iid, mate->iid, dist);
      
          if (1)
            {
              // See if the mate distance implied is outside of a 5-sigma range
        
              if (dist < dptr->lower || dist > dptr->upper) 
                {
                  //if (dist > 5000)
                  // fprintf(GlobalData->stderrc,"* (" F_CID "," F_CID ") lib:" F_CID " is bad due to distance problems " F_COORD " is outside [" F_COORD "," F_COORD "]\n",
                  //	frag->iid, mate->iid, frag->dist, dist, dptr->lower, dptr->upper);
                  frag->flags.bits.edgeStatus = UNTRUSTED_EDGE_STATUS;
                  mate->flags.bits.edgeStatus = UNTRUSTED_EDGE_STATUS;
                  dptr->numBad++;


                }
              else
                {
                  frag->flags.bits.edgeStatus = TRUSTED_EDGE_STATUS;
                  mate->flags.bits.edgeStatus = TRUSTED_EDGE_STATUS;
          
                  //if (dist > 5000)
                  //fprintf(GlobalData->stderrc,"* (" F_CID "," F_CID ") lib:" F_CID " is trusted due to distance " F_COORD " being inside [" F_COORD "," F_COORD "]\n",
                  //	frag->iid, mate->iid, frag->dist, dist, dptr->lower, dptr->upper);
          
                  dptr->numSamples++;
                  dptr->mu += dist;
                  dptr->sigma += (((double)dist)*(double)dist);

                }
            }
        }
    
      //        UnloadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, node->id, graph->type == CI_GRAPH);
    
      // Mark unitigs as potential Rocks and Stones
      if (operateOnNodes == UNITIG_OPERATIONS)
        {
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
          node->flags.bits.isPotentialRock = rock;
          node->flags.bits.isPotentialStone = stone;
      
#if 0
          fprintf(stderr,"* Unitig " F_CID " has %d external links r:%s s:%s\n",
                  node->id, numExternalLinks, 
                  (rock?" YES ": "NO"),
                  (stone?" YES ": "NO"));
#endif
        }
    }
  
  fprintf(GlobalData->stderrc, "* ComputeMatePairStats some mate data:\n");
  fprintf(GlobalData->stderrc, "* num total - chaff                 %d\n",numTotalFrags);
  fprintf(GlobalData->stderrc, "* num chaff nodes                   %d\n",numChaff);
  fprintf(GlobalData->stderrc, "* num singleton                     %d\n",numSingle);
  fprintf(GlobalData->stderrc, "* num no link                       %d\n",numNolink);
  fprintf(GlobalData->stderrc, "* num ctg not internal              %d\n",numCtgNotInternal);
  fprintf(GlobalData->stderrc, "* num 5' > 3' so skip               %d\n",numReverse);
  fprintf(GlobalData->stderrc, "* num different unitig              %d\n",numOtherUtg);
  fprintf(GlobalData->stderrc, "* num same orientation              %d\n",numSameOrient);
  fprintf(GlobalData->stderrc, "* num greater then 5 stddev distant %d\n",numNot5stddev);
  fprintf(GlobalData->stderrc, "* num different scaffold            %d\n",numDiffScaf);
  fprintf(GlobalData->stderrc, "* num mate != sourceInt             %d\n",numMateNotSource);
  fprintf(GlobalData->stderrc, "* num NULL mate                     %d\n",numNullMate);

  if (operateOnNodes == UNITIG_OPERATIONS)
    {
      fprintf(GlobalData->stderrc,
              "* ComputeMatePairStats has marked %d/%d unitigs as potential rocks +  %d/%d as potential stones\n",
              numPotentialRocks, (int) GetNumGraphNodes(graph),
              numPotentialStones, (int) GetNumGraphNodes(graph));
    }
  
#if 1
  // mjf
  
  for (i = 1; i < GetNumDistTs(ScaffoldGraph->Dists); i++)
    {
      fprintf( stderr, "GetNumCDS_COORD_ts( GetDistT(ScaffoldGraph->Dists, " F_CID ")->samples) = %d\n", 
               i, (int) GetNumCDS_COORD_ts( GetDistT(ScaffoldGraph->Dists, i)->samples));
      if ( GetNumCDS_COORD_ts( GetDistT(ScaffoldGraph->Dists, i)->samples) > maxSamples)
        maxSamples = GetNumCDS_COORD_ts( GetDistT(ScaffoldGraph->Dists, i)->samples);
    }
  
  fprintf( stderr, "maxSamples = %d\n", maxSamples);
  
  // now sort the samples, mates, and frags arrays, based on samples
  for (i = 1; i < GetNumDistTs(ScaffoldGraph->Dists); i++)
    {
      // MateInfoT *matePairs;
      MateInfoT matePairs[maxSamples]; // matePairs[maxSamples];
      int icnt;
      CDS_COORD_t newLower, newUpper;
      CDS_COORD_t median = 0, lowerSigma = 0, upperSigma = 0;
    
      dptr = GetDistT(ScaffoldGraph->Dists, i);
      if(dptr->numReferences == 0)
        continue;
      if (dptr->numSamples == 0 || dptr->numSamples == 1)
        continue;
    
      for ( icnt = 0; icnt < GetNumCDS_COORD_ts( dptr->samples ); icnt++)
        {
          matePairs[icnt].samples = *GetCDS_COORD_t( dptr->samples, icnt);
          matePairs[icnt].frags = *GetCDS_CID_t( dptrFrags[i], icnt);
          matePairs[icnt].mates = *GetCDS_CID_t( dptrMates[i], icnt);
        }
    
      qsort( matePairs, GetNumCDS_COORD_ts( dptr->samples ),
             sizeof(MateInfoT), &compDists);
    
      //if (dptr->mean > 5000 && dptr->numSamples > minSamplesForOverride)
      //for ( icnt = 0; icnt < GetNumCDS_COORD_ts( dptr->samples ); icnt++)
      //fprintf( stderr, "dptr: " F_CID ", matePairs[%d].samples: " F_COORD " (" F_CID ", " F_CID ")\n",
      //         i, icnt, matePairs[icnt].samples,
      //	       matePairs[icnt].frags, matePairs[icnt].mates);
    
      // now find median
      median = matePairs[ (int) (GetNumCDS_COORD_ts( dptr->samples ) / 2)].samples;
    
      // find lower and upper "sigmas"
      lowerSigma = median - matePairs[ (int) ((0.5 - 0.34) * GetNumCDS_COORD_ts( dptr->samples ))].samples;
      upperSigma = - median + matePairs[ (int) ((0.5 + 0.34) * GetNumCDS_COORD_ts( dptr->samples ))].samples;	
    
      newLower = median - 5 * MAX (lowerSigma, upperSigma);
      if ( newLower < 0 )
        newLower = 0;
      newUpper = median + 5 * MAX (lowerSigma, upperSigma);
    
      fprintf( stderr, "\nlib " F_CID ", numSamples: %d, orig mean, sig: ( %.2f, %.2f), calc mean, sig: (%.2f, %.2f) median: " F_COORD "\n",
               i, dptr->numSamples, dptr->mean, dptr->stddev, dptr->mu/dptr->numSamples, 
               sqrt((dptr->sigma - (dptr->mu*dptr->mu)/dptr->numSamples) / (dptr->numSamples - 1)),
               median);
      fprintf( stderr, "dptr->lower: " F_COORD ", dptr->upper: " F_COORD "\n", dptr->lower, dptr->upper);
      fprintf( stderr, "  dptr->min: " F_COORD ",   dptr->max: " F_COORD "\n", dptr->min, dptr->max);
      fprintf( stderr, "lowerSigma: " F_COORD ", upperSigma: " F_COORD "\n", lowerSigma, upperSigma);						  
      fprintf( stderr, "newLower: " F_COORD ", newUpper: " F_COORD "\n", newLower, newUpper);
    
      // now reset the trusted flag if necessary
      // first see if there are edges marked untrusted that are now considered trusted
      // lower set
      if (dptr->lower > newLower)
        {
          fprintf( stderr, "GetNumCDS_COORD_ts( dptr->samples ): %d\n",
                   (int) GetNumCDS_COORD_ts( dptr->samples )); 
          for ( icnt = 0; icnt < GetNumCDS_COORD_ts( dptr->samples ); icnt++)
            {
              if (matePairs[icnt].samples < dptr->lower && matePairs[icnt].samples > newLower)
                {
                  CIFragT *frag, *mate;
          
                  frag = GetCIFragT( ScaffoldGraph->CIFrags, 
                                     GetInfoByIID(ScaffoldGraph->iidToFragIndex, matePairs[icnt].frags)->fragIndex);
                  mate = GetCIFragT( ScaffoldGraph->CIFrags,
                                     GetInfoByIID(ScaffoldGraph->iidToFragIndex, matePairs[icnt].mates)->fragIndex);
          
                  fprintf( stderr, "1 reclassifying samples[%d] (" F_CID ", " F_CID ") from UNTRUSTED to TRUSTED\n", 
                           icnt, frag->iid, mate->iid);
          
                  frag->flags.bits.edgeStatus = TRUSTED_EDGE_STATUS;
                  mate->flags.bits.edgeStatus = TRUSTED_EDGE_STATUS;
                  dptr->numBad--;
                  dptr->numSamples++;
                  dptr->mu += matePairs[icnt].samples;
                  dptr->sigma += (((double) matePairs[icnt].samples) * (double) matePairs[icnt].samples);
                }
            }
        }
    
#if 1
      // upper set
      if (dptr->upper < newUpper)
        {
          fprintf( stderr, "GetNumCDS_COORD_ts( dptr->samples ): %d\n",
                   (int) GetNumCDS_COORD_ts( dptr->samples )); 
          for ( icnt = 0; icnt < GetNumCDS_COORD_ts( dptr->samples ); icnt++)
            {	
              if (matePairs[icnt].samples > dptr->upper && matePairs[icnt].samples < newUpper)
                {
                  CIFragT *frag, *mate;
          
                  frag = GetCIFragT( ScaffoldGraph->CIFrags, 
                                     GetInfoByIID(ScaffoldGraph->iidToFragIndex, matePairs[icnt].frags)->fragIndex);
                  mate = GetCIFragT( ScaffoldGraph->CIFrags,
                                     GetInfoByIID(ScaffoldGraph->iidToFragIndex, matePairs[icnt].mates)->fragIndex);
          
                  fprintf( stderr, "2 reclassifying samples[%d] (" F_CID ", " F_CID ") from UNTRUSTED to TRUSTED\n", 
                           icnt, frag->iid, mate->iid);
          
                  frag->flags.bits.edgeStatus = TRUSTED_EDGE_STATUS;
                  mate->flags.bits.edgeStatus = TRUSTED_EDGE_STATUS;
                  dptr->numBad--;
                  dptr->numSamples++;
                  dptr->mu += matePairs[icnt].samples;
                  dptr->sigma += (((double) matePairs[icnt].samples) * (double) matePairs[icnt].samples);
                }
            }
        }
#endif
    
#if 1	
      // now see if there are edges marked trusted that are now considered untrusted
      // lower set
      if (dptr->lower < newLower)
        {
          fprintf( stderr, "GetNumCDS_COORD_ts( dptr->samples ): %d\n",
                   (int) GetNumCDS_COORD_ts( dptr->samples )); 
          for ( icnt = 0; icnt < GetNumCDS_COORD_ts( dptr->samples ); icnt++)
            {
              if (matePairs[icnt].samples > dptr->lower && matePairs[icnt].samples < newLower)
                {
                  CIFragT *frag, *mate;
          
                  frag = GetCIFragT( ScaffoldGraph->CIFrags, 
                                     GetInfoByIID(ScaffoldGraph->iidToFragIndex, matePairs[icnt].frags)->fragIndex);
                  mate = GetCIFragT( ScaffoldGraph->CIFrags,
                                     GetInfoByIID(ScaffoldGraph->iidToFragIndex, matePairs[icnt].mates)->fragIndex);
          
                  fprintf( stderr, "3 reclassifying samples[%d] (" F_CID ", " F_CID ") from TRUSTED to UNTRUSTED\n", 
                           icnt, frag->iid, mate->iid);
          
                  frag->flags.bits.edgeStatus = UNTRUSTED_EDGE_STATUS;
                  mate->flags.bits.edgeStatus = UNTRUSTED_EDGE_STATUS;
                  dptr->numBad++;
                  dptr->numSamples--;
                  dptr->mu -= matePairs[icnt].samples;
                  dptr->sigma -= (((double) matePairs[icnt].samples) * (double) matePairs[icnt].samples);
                }
            }
        }
      // upper set
      if (dptr->upper > newUpper)
        {
          fprintf( stderr, "GetNumCDS_COORD_ts( dptr->samples ): %d\n",
                   (int) GetNumCDS_COORD_ts( dptr->samples )); 
          for ( icnt = 0; icnt < GetNumCDS_COORD_ts( dptr->samples ); icnt++)
            {	
              if (matePairs[icnt].samples < dptr->upper && matePairs[icnt].samples > newUpper)
                {
                  CIFragT *frag, *mate;
          
                  frag = GetCIFragT( ScaffoldGraph->CIFrags, 
                                     GetInfoByIID(ScaffoldGraph->iidToFragIndex, matePairs[icnt].frags)->fragIndex);
                  mate = GetCIFragT( ScaffoldGraph->CIFrags,
                                     GetInfoByIID(ScaffoldGraph->iidToFragIndex, matePairs[icnt].mates)->fragIndex);
          
                  fprintf( stderr, "4 reclassifying samples[%d] (" F_CID ", " F_CID ") from TRUSTED to UNTRUSTED\n", 
                           icnt, frag->iid, mate->iid);
          
                  frag->flags.bits.edgeStatus = UNTRUSTED_EDGE_STATUS;
                  mate->flags.bits.edgeStatus = UNTRUSTED_EDGE_STATUS;
                  dptr->numBad++;
                  dptr->numSamples--;
                  dptr->mu -= matePairs[icnt].samples;
                  dptr->sigma -= (((double) matePairs[icnt].samples) * (double) matePairs[icnt].samples);
                }
            }
        }
#endif
      // free(matePairs);
    }
  
  for (i = 1; i < GetNumDistTs(ScaffoldGraph->Dists); i++)
    {
      DeleteVA_CDS_CID_t( dptrFrags[i] );
      DeleteVA_CDS_CID_t( dptrMates[i] );
    }
  
#endif
  
  /* now set mean, stddev, size & number of buckets */
  for(i = 1; i < GetNumDistTs(ScaffoldGraph->Dists); i++)
    {
      dptr = GetDistT(ScaffoldGraph->Dists, i);
      if(dptr->numReferences == 0)
        {
          //      fprintf(GlobalData->stderrc,"* No references for mate distance record %d\n", i);
          continue;
        }
      if (dptr->numSamples == 0)
        {
          fprintf(GlobalData->stderrc,"* No samples for mate distance record %d... numBad = %d/%d\n", 
                  i, dptr->numBad, dptr->numReferences);
          continue;
        }
      if(dptr->numSamples == 1)
        dptr->sigma = -1.0;
      else
        dptr->sigma = sqrt((dptr->sigma - (dptr->mu*dptr->mu)/dptr->numSamples) /
                           (dptr->numSamples - 1));
    
      dptr->mu = dptr->mu/dptr->numSamples;
    
      if(dptr->numSamples > minSamplesForOverride)
        {
          fprintf(GlobalData->stderrc,
                  "**** >>>>>>>>>  Distance record %d updated from %g +/- %g  ==> %g +/- %g  based on %d samples <<<<<<<< ****\n",
                  i,
                  dptr->mean, dptr->stddev,
                  dptr->mu, dptr->sigma,
                  dptr->numSamples);
      
          dptr->mean = dptr->mu;
          dptr->stddev = dptr->sigma;
          dptr->lower = dptr->mean - CGW_CUTOFF*dptr->stddev;
          dptr->upper = dptr->mean + CGW_CUTOFF*dptr->stddev;
        }
    
      fprintf(GlobalData->stderrc,
              "* Distance %d:  mean:%g(nominal %g, delta %g) stddev:%g (nominal %g, delta %g) numSamples:%d numRef:%d\n",
              i, 
              dptr->mu, dptr->mean, fabs((dptr->mu - dptr->mean)/dptr->mean),
              dptr->sigma, dptr->stddev, fabs((dptr->sigma - dptr->stddev)/dptr->stddev),
              dptr->numSamples, dptr->numReferences);
    
      dptr->bsize = dptr->sigma/CGW_NUM_BUCKETS;
      if (dptr->bsize > 1.) 
        {
          dptr->bnum = (dptr->max-dptr->min)/dptr->bsize + 1;
          dptr->bsize = ((float)(dptr->max-dptr->min))/dptr->bnum;
        } 
      else 
        {
          dptr->bnum = 1;
          dptr->bsize = dptr->max-dptr->min;
        }
#if 0
      fprintf(GlobalData->stderrc,"* id:%d bnum:%d bsize:%g min:" F_COORD " max:" F_COORD "\n",
              i, dptr->bnum, dptr->bsize, dptr->min, dptr->max);
#endif
      free(dptr->histogram);
      dptr->histogram = (int32 *) safe_malloc(sizeof(int32)*dptr->bnum);
      {
        int	j;
        for (j=0; j < dptr->bnum; ++j)
          dptr->histogram[j] = 0;
      }
    }
  
  fprintf(GlobalData->stderrc,"* Bucketizing Data\n");
  fflush(NULL);
  // Now bucketize data
  {
    // int ii;
    // void *ptr = (void *)ScaffoldGraph->sequenceDB->UnitigStore;
    int32 numDists = GetNumDistTs(ScaffoldGraph->Dists);
    DistT *dist = GetDistT(ScaffoldGraph->Dists,1);
    for(ii = 1; ii < numDists ; ii++, dist++){
      int numSamples = GetNumCDS_COORD_ts(dist->samples);
      CDS_COORD_t *samplep = GetCDS_COORD_t(dist->samples,0);
      int j;
      if(!dist->samples || dist->bsize == 0.0)
        continue;
      for(j = 0; j < numSamples ; j++, samplep++)
        {
          CDS_COORD_t sample = *samplep;
          int32 binNum = ((float)sample - (float)dist->min)/dist->bsize;
          binNum = MIN(binNum, dist->bnum -1);
          binNum = MAX(binNum,0); // sample can be negative
          dist->histogram[binNum]++;
        }	  
    }
  }
  {
    // int ii;
    for(ii = 1; ii < GetNumDistTs(ScaffoldGraph->Dists); ii++)
      {
        int j;
        DistT *dist = GetDistT(ScaffoldGraph->Dists, ii);
        if(dist->samples && dist->bsize > 0.0 && dist->numSamples > 0)
          {
            fprintf(GlobalData->stderrc,"* Distance Record %d min:" F_COORD " max:" F_COORD "  mu:%g sigma:%g samples:%d bad:%d references:%d\n",
                    ii, dist->min, dist->max, dist->mu, dist->sigma, dist->numSamples, dist->numBad, dist->numReferences);
            fflush(NULL);
            for(j = 0; j < dist->bnum; j++)
              {
                int32 binVal = dist->histogram[j];
          
                if (binVal > 0 )
                  fprintf(GlobalData->stderrc,"* [%5d,%5d]\t%d\n",
                          (int32)(dist->min + j * dist->bsize), (int32)(dist->min + (j + 1) * dist->bsize), binVal);
              }
          }
      }
  }
  // output a histogram file for each library
  for( ii = 1; ii < GetNumDistTs(ScaffoldGraph->Dists); ii++)
    {
      int j;
      DistT *dist = GetDistT(ScaffoldGraph->Dists,ii);
      if(dist->samples && dist->bsize > 0.0 && dist->numSamples > 0)
        {
          FILE *fout;
          char  filename[1024];

          sprintf( filename, "stat/%s.distlib_%d.cgm", instance_label, ii);
          fprintf( stderr, "writing file %s\n", filename);
          fout = fopen( filename,"w");
          AssertPtr( fout );

          fprintf( fout, "lib %d mu %g sigma %g\n", ii, dist->mu, dist->sigma );
          {
            int numSamples = GetNumCDS_COORD_ts(dist->samples);
            CDS_COORD_t samps[numSamples];
            CDS_COORD_t *samplep = GetCDS_COORD_t(dist->samples,0);
            for( j = 0; j < numSamples; j++, samplep++) {
              samps[ j ] = *samplep;
            }
            qsort( samps, numSamples, sizeof(CDS_COORD_t), &compareInt);
            for( j = 0; j < numSamples; j++) {
              fprintf( fout, "%d\n", samps[j]);
            }
            fclose( fout );
          }
        }
    }
  
  DeleteMultiAlignT(ma);
  fprintf(stderr,"* Done bucketizing data\n");
  fflush(NULL);
  
}



















// Clean hash table
// Move deleted nodes and edges to their respective free lists
void RecycleDeletedGraphElements(GraphCGW_T *graph){
  EdgeCGW_T *end = NULL;
  
#if 0
  CDS_CID_t nextIndex = graph->freeEdgeHead;
  // Find the end of the free list;
#else
  // Find the end of the tobefree list (this should be MUCH shorter than looking for teh end of the free list )
  CDS_CID_t nextIndex = graph->tobeFreeEdgeHead;
#endif
  
  while(nextIndex != NULLINDEX){
    end = GetGraphEdge(graph, nextIndex);
    assert(end->flags.bits.isDeleted);
    assert(nextIndex != end->referenceEdge); // We are in a cycle!
    nextIndex = end->referenceEdge;
  }
  
#if 0
  // Tack the tobe freed edges onto the end of the free list
  if(end){
    end->referenceEdge = graph->tobeFreeEdgeHead;
  }else{
    graph->freeEdgeHead = graph->tobeFreeEdgeHead;
  }
#else
  // Tack the tobe freed edges onto the head of the free list
  if(end){
    end->referenceEdge = graph->freeEdgeHead;
    graph->freeEdgeHead = graph->tobeFreeEdgeHead;
  }else{
    ; // nothing to do, no 
  }
  
#endif
  
  graph->tobeFreeEdgeHead = NULLINDEX;
  
  // Flush Deleted nodes from the hashtable
  
  // Free tobe freed nodes
  
}


int32 GetGappedMultipleCoverageInterval(GraphCGW_T *graph,
                                        CDS_CID_t cid,
                                        SeqInterval *interval,
                                        int end){
  
  MultiAlignT *ma = LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB,
                                                  cid,
                                                  graph->type == CI_GRAPH);
  //  MultiAlignT *ma = GetMultiAlignInStore(graph->maStore, cid);
  CDS_COORD_t length = (CDS_COORD_t) GetNumchars(ma->consensus);
  VA_TYPE(int) *coverageVA = CreateVA_int(GetNumchars(ma->consensus));
  int *coverage;
  SeqInterval range;
  int i;
  int found = FALSE;
  
  if(end == A_END){
    interval->end = range.end = MIN(500,length);
    interval->bgn = range.bgn = 0;
  }else{
    interval->bgn = range.bgn = MAX(0, length - 500);
    interval->end = range.end = length;
  }
  fprintf(GlobalData->stderrc,"Calling GetCoverageInMultiAlignT end %c range [" F_COORD "," F_COORD "] of [0," F_COORD "]\n",
	  (end == A_END?'A':'B'),
	  range.bgn, range.end,
	  length);
  fflush(NULL);
  GetCoverageInMultiAlignT(ma, range, coverageVA, TRUE);
  fprintf(GlobalData->stderrc,"Returned from GetCoverageInMultiAlignT length of coverageVA = %d\n", (int) GetNumints(coverageVA));
  fflush(NULL);
  
  coverage = Getint(coverageVA,0);
  
  
  if(end == A_END){
    fprintf(GlobalData->stderrc,"* A_END %d\n", (int) GetNumints(coverageVA) );
    for(i = 0; i < GetNumints(coverageVA); i++){
      if(coverage[i] > 1){
	interval->bgn = (CDS_COORD_t) i;
	found = TRUE;
	break;
      }
    }
    
  }else{
    found = FALSE;
    fprintf(GlobalData->stderrc,"* B_END %d\n", (int) GetNumints(coverageVA) );
    
    for(i = GetNumints(coverageVA) -1; i > 0; i--){
      if(coverage[i] > 1){
	interval->end = length - (CDS_COORD_t) i;
	fprintf(GlobalData->stderrc,"* BEnd interval end = " F_COORD " length = " F_COORD " i = %d\n", interval->end, length, i);
	found = TRUE;
	break;
      }
    }
  }
  DeleteVA_int(coverageVA);
  
  return found;
  
}

/***********************************************/
#define SEGLEN 50
static void dumpFastaRecord(FILE *stream, char *header, char *sequence){
  int i;
  CDS_COORD_t FragLen = (CDS_COORD_t) strlen(sequence);
  
  fprintf(stream,"> %s\n", header);
  for (i = 0; i < FragLen; i += SEGLEN)
    if (FragLen-i < SEGLEN)
      fprintf(stream,"%.*s\n",FragLen-i,sequence+i);
    else
      fprintf(stream,"%.*s\n",SEGLEN,sequence+i);
}
/***********************************************/

void   DumpNodeEndUngappedToFasta(FILE *fastaFile,
                                  GraphCGW_T *graph,
                                  CDS_CID_t cid,
                                  int nodeEnd,
                                  int reverseComplement){
  MultiAlignT *ma = LoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB,
                                                  cid,
                                                  graph->type == CI_GRAPH);
  //    MultiAlignT *ma = GetMultiAlignInStore(graph->maStore, cid);
  NodeCGW_T *node = GetGraphNode(graph,cid);
  VA_TYPE(char) *consensus = CreateVA_char(1024);
  VA_TYPE(char) *quality = CreateVA_char(1024);
  char buffer[1024];
  char dataBuffer[2048];
  CDS_COORD_t length = (CDS_COORD_t) strlen(Getchar(ma->consensus,0));
  SeqInterval interval;
  if(nodeEnd == A_END){
    interval.bgn = 0;
    interval.end = MIN(500, length);
  }else{
    interval.bgn = MAX(0,length - 500);
    interval.end = length;
  }
  sprintf(buffer,
          " Scaffold " F_CID " %s " F_CID " %s [" F_COORD "," F_COORD "]", 
          node->scaffoldID, 
          (graph->type == CONTIG_GRAPH? "Contig":"Unitig"), 
          node->id, 
          (nodeEnd == A_END? "A_END": "B_END"),
          interval.bgn, interval.end);
  GetMultiAlignUngappedConsensusFromInterval(ma,interval,consensus, quality);
  strcat(dataBuffer,"");
  dumpFastaRecord(fastaFile, buffer, Getchar(consensus,0));
  DeleteVA_char(consensus);
  DeleteVA_char(quality);
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



int isBacOnlyEdge(GraphCGW_T *graph, EdgeCGW_T *edge){
  EdgeCGW_T *e = edge;
  int bacOnly;
  
  if(edge->flags.bits.isRaw){
    if(edge->flags.bits.hasGuide)
      return TRUE;
    else
      return FALSE;
  }
  
  bacOnly = TRUE;
  for(e = GetGraphEdge(graph, edge->nextRawEdge);
      e != NULL;
      e = GetGraphEdge(graph, e->nextRawEdge)){
    if(!edge->flags.bits.hasGuide){
      bacOnly = FALSE;
      break;
    }
  }
  
  return bacOnly;
}


void EdgeDegree(GraphCGW_T *graph, EdgeCGW_T *edge,
                int32 *totalDegree, int32 *noBacDegree){
  EdgeCGW_T *e = edge;
  int32
    tdegree = 0, 
    ndegree = 0;
  
  if(edge->flags.bits.isRaw){
    tdegree = 1;
    if(edge->flags.bits.hasGuide){
      ndegree = 0;
    }else{
      ndegree = 1;
    }
  }else{
    
    for(e = GetGraphEdge(graph, edge->nextRawEdge);
        e != NULL;
        e = GetGraphEdge(graph, e->nextRawEdge)){
      tdegree++;
      if(!edge->flags.bits.hasGuide){
	ndegree++;
      }
    }
  }
  *noBacDegree = ndegree;
  *totalDegree = tdegree;
  
}


void CheckGraph(GraphCGW_T *graph){
  EdgeCGW_T *edge;
  NodeCGW_T *node;
  GraphNodeIterator nodes;
  GraphEdgeIterator edges;
  
  InitGraphNodeIterator(&nodes, graph, GRAPH_NODE_DEFAULT);
  
  while(NULL != (node = NextGraphNodeIterator(&nodes))){
    
    InitGraphEdgeIterator(graph, node->id, ALL_END, ALL_EDGES,
                          GRAPH_EDGE_DEFAULT, &edges);
    while(NULL != (edge = NextGraphEdgeIterator(&edges))){
      continue;
    }
  }
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
void  CheckUnitigs(CDS_CID_t startUnitigID)
{
  // GraphEdgeStatT stat;
  GraphNodeIterator Nodes;
  NodeCGW_T *node;
  MultiAlignT *ma = CreateEmptyMultiAlignT();
  int32 i, numFrags;
  int totalUnitigs = (int) GetNumGraphNodes( ScaffoldGraph->CIGraph );  
  
  fprintf( stderr,"* CheckUnitigs starting at unitig " F_CID "\n",
           startUnitigID);
  fflush( NULL );
  
  InitGraphNodeIterator( &Nodes, ScaffoldGraph->CIGraph, GRAPH_NODE_DEFAULT);
  while(NULL != ( node = NextGraphNodeIterator(&Nodes)))
    {
      if (node->id < startUnitigID)
        continue;
    
      if (node->flags.bits.isChaff && GlobalData->ignoreChaffUnitigs)
        continue;
    
      if((node->id % 100000) == 0)
        {
          fprintf( stderr,"* CheckUnitigs Node " F_CID " out of %d\n",
                   node->id, totalUnitigs);
          fflush( NULL );
        }
      // This will load the ma from the disk file
      ReLoadMultiAlignTFromSequenceDB(ScaffoldGraph->sequenceDB, ma, node->id, TRUE);    
      // BuildGraphEdgesFromMultiAlign(graph, node, ma, &stat, FALSE);
    
      assert(ma && node);
      numFrags = GetNumIntMultiPoss(ma->f_list);
    
#if 0
      fprintf(GlobalData->stderrc,
              "* CheckUnitigs on node " F_CID " with %d frags\n",
              node->id, numFrags);
#endif
    
      for(i = 0; i < numFrags; i++)
        {
          IntMultiPos *mp = GetIntMultiPos(ma->f_list, i);
          CIFragT *frag = GetCIFragT(ScaffoldGraph->CIFrags,
                                     (CDS_CID_t)mp->sourceInt);
          CDS_CID_t fragID = (CDS_CID_t) mp->sourceInt;
      
          if(frag->cid != node->id)
            {
              fprintf(GlobalData->stderrc, "* CheckUnitigs: frag " F_CID " with iid " F_CID " of cid " F_CID " has ci = " F_CID "!!!\n",
                      fragID,	frag->iid, node->id, frag->cid);
              // assert(0);
            }
        }
    }
  
  fflush(NULL);
  DeleteMultiAlignT(ma);  
}
