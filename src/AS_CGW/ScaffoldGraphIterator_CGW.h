
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
/* $Id: ScaffoldGraphIterator_CGW.h,v 1.3 2005-03-22 19:03:59 jason_miller Exp $ */
/*****************************************************************************
 *  ScaffoldGraphIterators
 *  
 *  Saul A. Kravitz 5/99
 *
 *  This is the proposed data structures and operations for the
 *  scaffold graph that will be used for the final phases of the Chunk
 *  Graph walker.
 *
 ******************************************************************************/
#ifndef SCAFFOLD_GRAPH_ITERATOR_H
#define SCAFFOLD_GRAPH_ITERATOR_H

#include "AS_CGW_dataTypes.h"
#include "Globals_CGW.h"
#include "ScaffoldGraph_CGW.h"



// Edge characteristics for edge iterators

#define ALL_TRUSTED_EDGES (TENTATIVE_TRUSTED_EDGE_STATUS | TRUSTED_EDGE_STATUS)
#define ALL_EDGES (TENTATIVE_TRUSTED_EDGE_STATUS | TRUSTED_EDGE_STATUS | \
		   TENTATIVE_UNTRUSTED_EDGE_STATUS | \
		   UNKNOWN_EDGE_STATUS | UNTRUSTED_EDGE_STATUS | \
		   LARGE_VARIANCE_EDGE_STATUS | INTER_SCAFFOLD_EDGE_STATUS)

#define ALL_INTERNAL_EDGES (TENTATIVE_TRUSTED_EDGE_STATUS | TRUSTED_EDGE_STATUS | \
		   TENTATIVE_UNTRUSTED_EDGE_STATUS | \
		   UNKNOWN_EDGE_STATUS | UNTRUSTED_EDGE_STATUS)


/* ************************************************************************ */
/*     Iterate over all CIEdges incident on a ChunkInstance                 */
/* ************************************************************************ */

typedef struct {
  // Indices of MERGED Edges
  CDS_CID_t next;  /* Index into CIEdges */
  CDS_CID_t curr;  /* Index into CIEdges */
  CDS_CID_t prev;  /* Index into CIEdges */
  // Indices of Raw Edges
  CDS_CID_t nextRaw;  /* Index into CIEdges */
  CDS_CID_t currRaw;  /* Index into CIEdges */
  CDS_CID_t prevRaw;  /* Index into CIEdges */
  int confirmedOnly;
  int rawOnly;
  int end; /* A_END, B_END or ALL_END */
  unsigned int edgeStatusSet; /* See CIEdgeStatus and TRUSTED_EDGES and ALL_EDGES */
  int verbose;
  CDS_CID_t cid;
  ScaffoldGraphT *graph;
  VA_TYPE(CIEdgeT) *edges;
}CIEdgeTIterator;


static void InitCIEdgeTIterator(ScaffoldGraphT *graph,
				  CDS_CID_t cid,
			     	  int rawOnly,      // if TRUE, raw Edges, if FALSE merged
				  int confirmedOnly,
				  int end, 
				  int edgeStatusSet,
				  int verbose,
				  CIEdgeTIterator *e){
  ChunkInstanceT *CI;

  assert(graph && e);

  CI = GetChunkInstanceT(graph->ChunkInstances, cid);
  e->curr = e->prev = NULLINDEX ;
  e->next = CI->edgeHead;
  e->currRaw = e->nextRaw = e->prevRaw = NULLINDEX;
  e->graph = graph;
  e->edges = graph->CIEdges;
  assert(end & ALL_END);
  e->confirmedOnly = confirmedOnly;
  e->rawOnly = rawOnly;
  e->verbose = verbose;
  e->end = end;
  e->edgeStatusSet = edgeStatusSet;
  e->cid = cid;
  if(e->verbose)
  fprintf(stderr,"* Iterator for CI " F_CID " end %d head = " F_CID " confirmed:%d raw:%d \n",
	  cid, e->end, e->next,e->confirmedOnly, e->rawOnly);
}


static CIEdgeT *NextCIEdgeTIterator(CIEdgeTIterator *e){
  int isA;
  CIEdgeT *r = (CIEdgeT *)NULL;
  CIEdgeT *retEdge = (CIEdgeT *)NULL;
  assert(e->graph && (e->end & ALL_END));

  if(e->verbose)
    fprintf(stderr,"* NextCIEdgeTIterator nextRaw:" F_CID " prev:" F_CID " curr:" F_CID " next:" F_CID "\n",
	  e->nextRaw, e->prev, e->curr, e->next);

  if(e->nextRaw == NULLINDEX){
    if(e->currRaw != NULLINDEX){ // do this once
      e->prevRaw = e->currRaw;
      e->currRaw = e->nextRaw;
    }
    if(e->next == NULLINDEX){
      if(e->curr != NULLINDEX){ // do this once
	e->prev = e->curr;
	e->curr = e->next;
      }
      if(e->verbose)
	fprintf(stderr,"* Fell off the end \n");
      return retEdge;
    }
  }

  while(retEdge == (CIEdgeT *)NULL &&
	(e->nextRaw != NULLINDEX || e->next != NULLINDEX)){
    ChunkOrientationType orient;
  
    if(e->verbose)
      fprintf(stderr,"* In loop (" F_CID "," F_CID "," F_CID ")\n",
	      e->prev, e->curr, e->next);

    if(e->nextRaw != NULLINDEX){ /* We are iterating through the raw edges of a merged edge */
      // These edges are guaranteed to be the right orientation and end points
      if(e->verbose)
	fprintf(stderr,"* Getting a raw edge\n");

      r = GetCIEdgeT(e->edges, e->nextRaw);
      AssertPtr(r);

      assert(r->idA == e->cid ||
	     r->idB == e->cid);

      e->prevRaw = e->currRaw;
      e->currRaw = e->nextRaw;
      e->nextRaw = r->nextRawEdge;
      retEdge = r;
    }else{              /* We are iterating through the merged and unmerged raw edges */
      r = GetCIEdgeT(e->edges, e->next);
      AssertPtr(r);

      assert(r->idA == e->cid ||
	     r->idB == e->cid);

      // Flip orientation to accomodate canonical graph form
      orient = GetEdgeOrientationWRT(r, e->cid);

      /* Check for correct end and confirmed status */
      switch(orient){
	/* EdgeMate from the A-End */
      case BA_BA:
      case BA_AB:
	if((e->end & A_END) &&
	   (e->confirmedOnly == 0 || isConfirmedEdge(r))){
	  retEdge = r;
	}else{
	  if(e->verbose)
	    fprintf(stderr,"* Skipping edge (" F_CID "," F_CID ") with orient:%c orientWRT %c (e->end&A_END) = %d confirmedOnly = %d\n",
		    r->idA, r->idB, r->orient, orient,e->end&A_END, e->confirmedOnly);
	}
	break;
      case AB_BA:
      case AB_AB:
	if((e->end & B_END) &&
	   (e->confirmedOnly == 0 || isConfirmedEdge(r))){
	  retEdge = r;
	}else{
	  if(e->verbose)
	    fprintf(stderr,"* Skipping edge (" F_CID "," F_CID ") with orient %c orientWRT:%c (e->end&B_END) = %d confirmedOnly = %d\n",
		     r->idA, r->idB, r->orient, orient, e->end&B_END, e->confirmedOnly);
	}
	break;
      default:
	assert(0);
      }
      e->prev = e->curr;
      e->curr = e->next;
      isA = (r->idA == e->cid);
      if(isA){
	e->next = r->nextALink;
      }else{
	e->next = r->nextBLink;
      }

      // Found a top level (not a raw under a merged) that is inappropriate
      // Only top level edges are marked with status values
      if(retEdge && !((uint32)GetEdgeStatus(retEdge) & (uint32)e->edgeStatusSet)){
	  if(e->verbose)
	    fprintf(stderr,"* Looking for status 0x%x only, found an edge with status 0x%x\n",
		    e->edgeStatusSet, GetEdgeStatus(retEdge));
	  retEdge = NULL;    // This stops us from exiting the loop */
      }
      /* If we are iterating over raw edges, dive into teh raw edge list */
      /* assert((!retEdge->flags.bits.isRaw &&
	      (retEdge->nextRawEdge != NULLINDEX)) ||
	     (retEdge->flags.bits.isRaw &&
	     (retEdge->nextRawEdge == NULLINDEX))); */
      if(retEdge && !retEdge->flags.bits.isRaw && e->rawOnly){
	  if(e->verbose)
	    fprintf(stderr,"* Looking for raw only, found a merged edge\n");
	  e->nextRaw = r->nextRawEdge;
	  retEdge = NULL;    // This stops us from exiting the loop */
      }
    }
  }
  if(retEdge && e->verbose)
    fprintf(stderr,"* Found CIEdge (" F_CID "," F_CID ") with orient %c and e->end = %d status:%d (" F_CID "," F_CID "," F_CID ")\n",
	    retEdge->idA, retEdge->idB,
	    retEdge->orient, e->end,
	    GetEdgeStatus(retEdge),
	    e->prev, e->curr, e->next);

  assert(e->curr == NULLINDEX || e->curr != e->next);

  return retEdge;
}

/* *********************************************************************** */
/*       Iterate over all CIEdges incident on a ChunkInstance              */
/* *********************************************************************** */

typedef CIEdgeTIterator SEdgeTIterator;
#define NextSEdgeTIterator NextCIEdgeTIterator


static void InitSEdgeTIterator(ScaffoldGraphT *graph,
				  CDS_CID_t sid,
			     	  int rawOnly,      // if TRUE, raw Edges, if FALSE merged
				  int confirmedOnly,
				  int end, 
				  int verbose,
				  SEdgeTIterator *e){
  CIScaffoldT *CIS;

  assert(graph && e);

  CIS = GetCIScaffoldT(graph->CIScaffolds, sid);
  e->curr = e->prev = NULLINDEX ;
  e->next = CIS->edgeHead;
  e->currRaw = e->nextRaw = e->prevRaw = NULLINDEX;
  e->graph = graph;
  e->edges = graph->SEdges;
  assert(end & ALL_END);
  e->confirmedOnly = confirmedOnly;
  e->rawOnly = rawOnly;
  e->verbose = verbose;
  e->end = end;
  e->edgeStatusSet = ALL_EDGES;
  e->cid = sid;
  if(verbose)
  fprintf(stderr,"* Iterator for Scaffold " F_CID " end %d  head = " F_CID " end = %d  confirmed=%d\n",
	  sid, end, e->next, e->end, e->confirmedOnly);
}

/* *********************************************************************** */
/*       Iterate over all ChunkInstances in a CIScaffoldT                  */
/* *********************************************************************** */

typedef struct {
  CDS_CID_t next;  /* Index into Contigs */
  CDS_CID_t curr;  /* Index into Contigs */
  CDS_CID_t prev;  /* Index into Contigs */
  int aEndToBEnd; /* True if we are scanning from a->b, false if b->a */
  int verbose;
  CDS_CID_t sid;   /* Scaffold ID */
  ScaffoldGraphT *graph;
}CIScaffoldTIterator;

/* Start from the specified CI instead of the end. */

static void InitCIScaffoldTIteratorFromCI(ScaffoldGraphT *graph,
					  CIScaffoldT *scaffold,
					  CDS_CID_t indexOfCI,
					  int aEndToBEnd,
					  int verbose,
					  CIScaffoldTIterator *e){
  ChunkInstanceT *CI = GetGraphNode(graph->RezGraph, indexOfCI);

  assert(graph && e && scaffold);
  assert(CI->scaffoldID == scaffold->id && !CI->flags.bits.isDead);
  e->curr = e->prev = NULLINDEX;
  e->next = indexOfCI;
  e->aEndToBEnd = aEndToBEnd;
  e->graph = graph;
  e->verbose = verbose;
  e->sid = scaffold->id;
  if(verbose) {
    fprintf(stderr,
            "* Iterator for CIScaffold " F_CID " from " F_CID " end = %s  head = " F_CID " scaffold (" F_CID "," F_CID ") \n",
            e->sid, e->next, (aEndToBEnd? "a->b": "b->a"),
            e->next, scaffold->info.Scaffold.AEndCI, scaffold->info.Scaffold.BEndCI);
  }
}


static void InitCIScaffoldTIterator(ScaffoldGraphT *graph,
				    CIScaffoldT *scaffold,
				    int aEndToBEnd,
				    int verbose,
				    CIScaffoldTIterator *e){
  assert(graph && e && scaffold);

  e->curr = e->prev = NULLINDEX;
  e->next = (aEndToBEnd? scaffold->info.Scaffold.AEndCI: scaffold->info.Scaffold.BEndCI);
  e->aEndToBEnd = aEndToBEnd;
  e->graph = graph;
  e->verbose = verbose;
  e->sid = scaffold->id;
  if(verbose)
  fprintf(stderr,"* Iterator for CIScaffold " F_CID " end = %s  head = " F_CID " scaffold (" F_CID "," F_CID ") \n",
	  e->sid, (aEndToBEnd? "a->b": "b->a"), e->next, scaffold->info.Scaffold.AEndCI, scaffold->info.Scaffold.BEndCI);
}


static ChunkInstanceT *NextCIScaffoldTIterator(CIScaffoldTIterator *e){
  ChunkInstanceT *r       = (ChunkInstanceT *)0;
  AssertPtr(e->graph);

  if(e->next <= NULLINDEX){
    if(e->verbose)
      fprintf(stderr,"* NextCIScaffoldTIterator returning NULL\n");
    return r;
  }

    r = GetGraphNode(e->graph->RezGraph, e->next);
    AssertPtr(r);

    e->prev = e->curr;
    e->curr = e->next;
    if(e->aEndToBEnd){
      e->next = r->BEndNext;
    }else{
      e->next = r->AEndNext;
    }

    if(e->verbose){
      if(r)
	  fprintf(stderr,"* Found CI " F_CID " in scaffold " F_CID " next = " F_CID "\n",
		  r->id, e->sid, e->next);
      else
	  fprintf(stderr,"* Found CI NULL in scaffold " F_CID "\n",
		  r->id);
    }

    return r;
}

static ChunkInstanceT *PrevCIScaffoldTIterator(CIScaffoldTIterator *e){
  ChunkInstanceT *r       = (ChunkInstanceT *)0;
  AssertPtr(e->graph);

  if(e->prev <= NULLINDEX){
    if(e->verbose)
      fprintf(stderr,"* PrevCIScaffoldTIterator returning NULL\n");
    return r;
  }

    r = GetGraphNode(e->graph->RezGraph, e->prev);
    AssertPtr(r);

    e->next = e->curr;
    e->curr = e->prev;
    if(e->aEndToBEnd){
      e->prev = r->AEndNext;
    }else{
      e->prev = r->BEndNext;
    }

    if(e->verbose){
      if(r)
	  fprintf(stderr,"* Found CI " F_CID " in scaffold " F_CID " next = " F_CID "\n",
		  r->id, e->sid, e->next);
      else
	  fprintf(stderr,"* Found CI NULL in scaffold " F_CID "\n",
		  r->id);
    }

    return r;
}

static ChunkInstanceT *GetNextFromCIScaffoldT(CIScaffoldTIterator *e){
  ChunkInstanceT *r       = (ChunkInstanceT *)0;
  AssertPtr(e->graph);

  if(e->next <= NULLINDEX){
    if(e->verbose)
      fprintf(stderr,"* GetNextFromCIScaffoldT returning NULL\n");
    return r;
  }

  r = GetGraphNode(e->graph->RezGraph, e->next);
  AssertPtr(r);
  return r;
}

static ChunkInstanceT *GetCurrFromCIScaffoldT(CIScaffoldTIterator *e){
  ChunkInstanceT *r       = (ChunkInstanceT *)0;
  AssertPtr(e->graph);

  if(e->curr <= NULLINDEX){
    if(e->verbose)
      fprintf(stderr,"* GetCurrFromCIScaffoldT returning NULL\n");
    return r;
  }

  r = GetGraphNode(e->graph->RezGraph, e->curr);
  AssertPtr(r);
  return r;
}

static ChunkInstanceT *GetPrevFromCIScaffoldT(CIScaffoldTIterator *e){
  ChunkInstanceT *r       = (ChunkInstanceT *)0;
  AssertPtr(e->graph);

  if(e->prev <= NULLINDEX){
    if(e->verbose)
      fprintf(stderr,"* GetPrevFromCIScaffoldT returning NULL\n");
    return r;
  }

  r = GetGraphNode(e->graph->RezGraph, e->prev);
  AssertPtr(r);
  return r;
}

static ChunkInstanceT *GetNextGivenCIFromCIScaffoldT(CIScaffoldTIterator *e,
						     ChunkInstanceT *given){
  ChunkInstanceT *r       = (ChunkInstanceT *)0;
  CDS_CID_t indexR;
  AssertPtr(e->graph);
  AssertPtr(given);

  if(e->aEndToBEnd){
    indexR = given->BEndNext;
  }else{
    indexR = given->AEndNext;
  }
  if(indexR <= NULLINDEX){
    if(e->verbose)
      fprintf(stderr,"* GetNextGivenCIFromCIScaffoldT returning NULL\n");
    return r;
  }

  r = GetGraphNode(e->graph->RezGraph, indexR);
  AssertPtr(r);
  return r;
}

static ChunkInstanceT *GetPrevGivenCIFromCIScaffoldT(CIScaffoldTIterator *e,
						     ChunkInstanceT *given){
  ChunkInstanceT *r       = (ChunkInstanceT *)0;
  CDS_CID_t indexR;
  AssertPtr(e->graph);
  AssertPtr(given);

  if(e->aEndToBEnd){
    indexR = given->AEndNext;
  }else{
    indexR = given->BEndNext;
  }
  if(indexR <= NULLINDEX){
    if(e->verbose)
      fprintf(stderr,"* GetPrevGivenCIFromCIScaffoldT returning NULL\n");
    return r;
  }

  r = GetGraphNode(e->graph->RezGraph, indexR);
  AssertPtr(r);
  return r;
}


/* *********************************************************************** */
/*       Iterate over all ChunkInstances in a ContigT                      */
/* *********************************************************************** */

typedef struct {
  CDS_CID_t next;  /* Index into ChunkInstances */
  CDS_CID_t curr;  /* Index into ChunkInstances */
  CDS_CID_t prev;  /* Index into ChunkInstances */
  int aEndToBEnd; /* True if we are scanning from a->b, false if b->a */
  int verbose;
  CDS_CID_t cid;   /* Contig ID */
  ScaffoldGraphT *graph;
}ContigTIterator;


static void InitContigTIterator(ScaffoldGraphT *graph,
				    CDS_CID_t cid, 
				    int aEndToBEnd,
				    int verbose,
				    ContigTIterator *e){
  ContigT *contig = GetGraphNode(graph->ContigGraph, cid);

  assert(graph && e && contig);
  assert(contig->type == CONTIG_CGW);
  e->prev = NULLINDEX;
  e->curr = NULLINDEX;
  e->next = (aEndToBEnd? contig->info.Contig.AEndCI: contig->info.Contig.BEndCI);
  e->aEndToBEnd = aEndToBEnd;
  e->graph = graph;
  e->verbose = verbose;
  e->cid = cid;
  if(verbose)
  fprintf(stderr,"* Iterator for Contig " F_CID " end = %s  head = " F_CID " scaffold (" F_CID "," F_CID ") \n",
	  cid, (aEndToBEnd? "a->b": "b->a"), e->next, contig->info.Contig.AEndCI, contig->info.Contig.BEndCI);
}


static ChunkInstanceT *NextContigTIterator(ContigTIterator *e){
  ChunkInstanceT *r       = (ChunkInstanceT *)0;
  AssertPtr(e->graph);

  if(e->next <= NULLINDEX){
    //    fprintf(stderr,"* NextCIScaffoldTIterator returning NULL\n");
    if(e->curr > NULLINDEX){
      e->prev = e->curr;
      e->curr = e->next;
    }
    return r;
  }

    r = GetGraphNode(e->graph->CIGraph, e->next);
    AssertPtr(r);

    e->prev = e->curr;
    e->curr = e->next;
    if(e->aEndToBEnd){
      e->next = r->BEndNext;
    }else{
      e->next = r->AEndNext;
    }

    if(e->verbose){
      if(r)
	  fprintf(stderr,"* Found CI " F_CID " in contig " F_CID " next = " F_CID "\n",
		  r->id, e->cid, e->next);
      else
	  fprintf(stderr,"* Found CI NULL in contig  " F_CID "\n",
		  r->id);
    }

    return r;
}


/* ************************************************************************* */
/* Add a fixed amount to the offsetAEnd and offsetBEnd of all
   CIs in the Scaffold so that the start of the Scaffold is {0.0,0.0}        */
/* ************************************************************************* */

static void NormalizeScaffoldOffsets(ScaffoldGraphT *graph,
				     CIScaffoldT *scaffold, int verbose){
  LengthT curBeginOffset;
  NodeCGW_T *endNodeA = GetGraphNode(graph->RezGraph,
				     scaffold->info.Scaffold.AEndCI);
  if(GetNodeOrient(endNodeA) == A_B){
    curBeginOffset = endNodeA->offsetAEnd;
  }else{
    curBeginOffset = endNodeA->offsetBEnd;
  }
  curBeginOffset.mean = - curBeginOffset.mean;
  curBeginOffset.variance = - curBeginOffset.variance;
  AddDeltaToScaffoldOffsets(graph, scaffold->id, endNodeA->id, TRUE, verbose,
			    curBeginOffset);
  return;
}

/* *********************************************************************** */
/* Call caffoldOffsets for all currently active scaffolds                  */
/* *********************************************************************** */

static void NormalizeAllScaffoldOffsets(ScaffoldGraphT *graph, int verbose){
  GraphNodeIterator scaffolds;
  CIScaffoldT *scaffold;

  if(verbose){
    fprintf(stderr, "Starting NormalizeAllScaffoldOffsets\n");
  }
  InitGraphNodeIterator(&scaffolds, graph->ScaffoldGraph, GRAPH_NODE_DEFAULT);
  while((scaffold = NextGraphNodeIterator(&scaffolds)) != NULL){
    if(isDeadCIScaffoldT(scaffold) || scaffold->type != REAL_SCAFFOLD){
      continue;
    }
    NormalizeScaffoldOffsets(graph, scaffold, verbose);
  }
}

#endif
