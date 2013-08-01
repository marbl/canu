
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

#ifndef SCAFFOLD_GRAPH_ITERATOR_H
#define SCAFFOLD_GRAPH_ITERATOR_H

static const char *rcsid_SCAFFOLD_GRAPH_ITERATOR_H = "$Id: ScaffoldGraphIterator_CGW.h,v 1.12 2012-08-15 15:07:15 brianwalenz Exp $";

#include "AS_CGW_dataTypes.H"
#include "Globals_CGW.H"
#include "ScaffoldGraph_CGW.H"



// Edge characteristics for edge iterators

#define ALL_TRUSTED_EDGES  (TENTATIVE_TRUSTED_EDGE_STATUS | \
                            TRUSTED_EDGE_STATUS)

#define ALL_EDGES          (TENTATIVE_TRUSTED_EDGE_STATUS | \
                            TRUSTED_EDGE_STATUS | \
                            TENTATIVE_UNTRUSTED_EDGE_STATUS | \
                            UNKNOWN_EDGE_STATUS | \
                            UNTRUSTED_EDGE_STATUS | \
                            LARGE_VARIANCE_EDGE_STATUS | \
                            INTER_SCAFFOLD_EDGE_STATUS)

#define ALL_INTERNAL_EDGES (TENTATIVE_TRUSTED_EDGE_STATUS | \
                            TRUSTED_EDGE_STATUS | \
                            TENTATIVE_UNTRUSTED_EDGE_STATUS | \
                            UNKNOWN_EDGE_STATUS | \
                            UNTRUSTED_EDGE_STATUS)



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
  ChunkInstanceT *CI = GetGraphNode(graph->ContigGraph, indexOfCI);

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

  r = GetGraphNode(e->graph->ContigGraph, e->next);
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

  r = GetGraphNode(e->graph->ContigGraph, e->prev);
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

  r = GetGraphNode(e->graph->ContigGraph, e->next);
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

  r = GetGraphNode(e->graph->ContigGraph, e->curr);
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

  r = GetGraphNode(e->graph->ContigGraph, e->prev);
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

  r = GetGraphNode(e->graph->ContigGraph, indexR);
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

  r = GetGraphNode(e->graph->ContigGraph, indexR);
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
  NodeCGW_T *endNodeA = GetGraphNode(graph->ContigGraph, scaffold->info.Scaffold.AEndCI);

  if (GetNodeOrient(endNodeA).isForward())
    curBeginOffset = endNodeA->offsetAEnd;
  else
    curBeginOffset = endNodeA->offsetBEnd;

  curBeginOffset.mean      = -curBeginOffset.mean;
  curBeginOffset.variance  = -curBeginOffset.variance;

  AddDeltaToScaffoldOffsets(graph, scaffold->id, endNodeA->id, TRUE, curBeginOffset);
}

/* *********************************************************************** */
/* Call caffoldOffsets for all currently active scaffolds                  */
/* *********************************************************************** */

static void NormalizeAllScaffoldOffsets(ScaffoldGraphT *graph, int verbose){
  GraphNodeIterator scaffolds;
  CIScaffoldT *scaffold;

  InitGraphNodeIterator(&scaffolds, graph->ScaffoldGraph, GRAPH_NODE_DEFAULT);

  while ((scaffold = NextGraphNodeIterator(&scaffolds)) != NULL) {
    if ((isDeadCIScaffoldT(scaffold)) ||
        (scaffold->type != REAL_SCAFFOLD))
      continue;
    NormalizeScaffoldOffsets(graph, scaffold, verbose);
  }
}

#endif
