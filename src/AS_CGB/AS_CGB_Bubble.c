
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
static char CM_ID[] 
= "$Id: AS_CGB_Bubble.c,v 1.6 2006-09-26 22:21:13 brianwalenz Exp $";

#include <stdlib.h>
#include <stdio.h>
#include "AS_global.h"
#include "AS_CGB_all.h"
#include "AS_CGB_methods.h"
#include "AS_CGB_Bubble_Graph.h"
#include "AS_CGB_Bubble_VertexSet.h"
#include "AS_CGB_Bubble_GraphMethods.h"
#include "AS_CGB_Bubble.h"
#include "AS_CGB_Bubble_Popper.h"
#include "AS_UTL_Hash.h"
#include "AS_CGB_Bubble_log.h"

FILE *BUB_LOG_G = NULL;

typedef struct BVSPair {
  BubVertexSet_t f, r;
} BVSPair;


#define AS_CGB_BUBBLE_ABS(x) (((x) < 0) ? (-(x)) : (x))

int 
_hash_vset_cmp(const void *vset1, const void *vset2)
{
  BVSPair *v1 = (BVSPair *) vset1;
  BVSPair *v2 = (BVSPair *) vset2;
  int comp_f = BVS_compare(v1->f, v2->f);
  
  if (comp_f != 0)
    return comp_f;
  else
    return BVS_compare(v1->r, v2->r);
}


int
_hash_vset_hash(const void *vset, int length)
{
  BVSPair *v = (BVSPair *) vset;

  return (int) (BVS_hash(v->f) + BVS_hash(v->r));
}


int
_is_initiation_node(int in_deg, int out_deg)
{
  return ((in_deg <= 1) && (out_deg > 1) && 
	  (out_deg < AS_CGB_BUBBLE_max_outdegree_G));
}


int
_is_termination_node(int in_deg, int out_deg)
{
  return _is_initiation_node(out_deg, in_deg);
}


AS_CGB_Bubble_List_t
_collect_bubbles(BubGraph_t bg, BubVertexSet *fwd, BubVertexSet *rvs, 
		 IntFragment_ID *top, int num_valid)
{
  IntFragment_ID f, bub_start;
  HashTable_AS *init_nodes = NULL;
  IntFragment_ID *i_node = NULL;
  AS_CGB_Bubble_List result;
  AS_CGB_Bubble_List_t *ins_h = &(result.next);
  BVSPair *bp_ins_keys = NULL, bp_find_key;

  memset(&result,0,sizeof(AS_CGB_Bubble_List));
  init_nodes  = CreateHashTable_AS(num_valid / 2, 
				   _hash_vset_hash,
				   _hash_vset_cmp);
  bp_ins_keys = (BVSPair *)safe_malloc(sizeof(BVSPair) * num_valid );
  result.next = NULL;

  for (f = 0; f < num_valid; ++f) 
    if (_is_initiation_node(BG_inDegree(bg, top[f], AS_CGB_BUBBLE_E_VALID),
			    BG_outDegree(bg, top[f], AS_CGB_BUBBLE_E_VALID)) &&
	!BVS_empty(&(fwd[top[f]])) &&
	!BVS_empty(&(rvs[top[f]]))) {
#if AS_CGB_BUBBLE_VERY_VERBOSE
      fprintf(BUB_LOG_G, "Inserting " F_IID " (" F_IID ") into the table.\n", top[f],
	      get_iid_fragment(BG_vertices(bg), top[f]));
#endif
      bp_ins_keys[f].f = &(fwd[top[f]]);
      bp_ins_keys[f].r = &(rvs[top[f]]);
      InsertInHashTable_AS(init_nodes, (void *) (&(bp_ins_keys[f])), 
			   sizeof(BVSPair), &(top[f]));
    }

  for (f = 0; f < num_valid; ++f) 
    if (_is_termination_node(BG_inDegree(bg, top[f], AS_CGB_BUBBLE_E_VALID),
			     BG_outDegree(bg, top[f], AS_CGB_BUBBLE_E_VALID))&&
	!BVS_empty(&(fwd[top[f]])) &&
	!BVS_empty(&(rvs[top[f]]))) { 
#if AS_CGB_BUBBLE_VERY_VERBOSE
      fprintf(BUB_LOG_G, "Looking for matches for " F_IID " (" F_IID ") in the table.  ", 
	      top[f], get_iid_fragment(BG_vertices(bg), top[f]));
#endif
      bp_find_key.f = &(fwd[top[f]]);
      bp_find_key.r = &(rvs[top[f]]);
      i_node = (IntFragment_ID *)
	LookupInHashTable_AS(init_nodes, (void *) &(bp_find_key), 
				    sizeof(BVSPair));
#if AS_CGB_BUBBLE_VERY_VERBOSE
      if (!i_node)
	fprintf(BUB_LOG_G, "None found.\n");
      else
	fprintf(BUB_LOG_G, "Found init node = " F_IID " (" F_IID ").\n", *i_node,
		get_iid_fragment(BG_vertices(bg), *i_node));
#endif

      if (i_node) {
	AS_CGB_Bubble_List_t new_bub = NULL;
	new_bub = (AS_CGB_Bubble_List *)safe_malloc(sizeof(AS_CGB_Bubble_List));
	bub_start = *i_node;
	new_bub->start = bub_start;
	new_bub->end = top[f];
	*ins_h = new_bub;
	ins_h = &(new_bub->next);
	*ins_h = NULL;
      }
    }

  DeleteHashTable_AS(init_nodes);
  safe_free(bp_ins_keys);
  return result.next;
}
 
 

void
_process_vertex(BubGraph_t bg, IntFragment_ID f, BubVertexSet *bvs,
		BG_E_Iter_t in, int in_deg, BG_E_Iter_t out, int out_deg)
{
  IntFragment_ID opp_f;
  IntEdge_ID e;
    
  if (in_deg > 0) {
    e = BGEI_cur(in);
    opp_f = BG_getOppositeVertex(bg, e, f);
    BVS_copy(&(bvs[opp_f]), &(bvs[f]));
    BVS_age(&(bvs[f]), 1);

    for (e = BGEI_next(bg, in, AS_CGB_BUBBLE_E_VALID); !BGEI_end(in);
	 e = BGEI_next(bg, in, AS_CGB_BUBBLE_E_VALID)) {
      opp_f = BG_getOppositeVertex(bg, e, f);
      BVS_intersect(&(bvs[f]), &(bvs[opp_f]));
    }
  }

  if (_is_initiation_node(in_deg, out_deg))
    BVS_insert(&(bvs[f]), f);
}


void
_forward_collect_sets(BubGraph_t bg, BubVertexSet *fwd, IntFragment_ID *top,
		      int num_valid)
{
  IntFragment_ID f;
  BG_E_Iter in, out;
  int in_deg, out_deg;

  for (f = 0; f < num_valid; ++f) {
    in_deg = BG_inDegree(bg, top[f], AS_CGB_BUBBLE_E_VALID);
    out_deg = BG_outDegree(bg, top[f], AS_CGB_BUBBLE_E_VALID);
    BGEI_bgn(bg, &in, top[f], bgeiIn, AS_CGB_BUBBLE_E_VALID);
    BGEI_bgn(bg, &out, top[f], bgeiOut, AS_CGB_BUBBLE_E_VALID);
    _process_vertex(bg, top[f], fwd, &in, in_deg, &out, out_deg);
  }
}


void
_reverse_collect_sets(BubGraph_t bg, BubVertexSet *fwd, IntFragment_ID *top,
		      int num_valid)
{
  IntFragment_ID f;
  BG_E_Iter in, out;
  int in_deg, out_deg;

  if (num_valid == 0)
    return;

  f = num_valid - 1;
  while (1) {
    in_deg = BG_outDegree(bg, top[f], AS_CGB_BUBBLE_E_VALID);
    out_deg = BG_inDegree(bg, top[f], AS_CGB_BUBBLE_E_VALID);
    BGEI_bgn(bg, &in, top[f], bgeiOut, AS_CGB_BUBBLE_E_VALID);
    BGEI_bgn(bg, &out, top[f], bgeiIn, AS_CGB_BUBBLE_E_VALID);
    _process_vertex(bg, top[f], fwd, &in, in_deg, &out, out_deg);
    if (f == 0)
      break;
    else
      --f;
  }
}


AS_CGB_Bubble_List_t
AS_CGB_Bubble_find_bubbles_with_graph(BubGraph_t bg, int sz, int age,
			   int max_outdegree)
{
  IntFragment_ID *top_order = NULL;
  IntFragment_ID num_frags, num_valid, f;
  BubVertexSet *fwd = NULL;
  BubVertexSet *rvs = NULL;
  AS_CGB_Bubble_List_t result = NULL;
  
  if (sz > 0)
    AS_CGB_BUBBLE_set_size_G = sz;
  if (age > 0)
    AS_CGB_BUBBLE_max_age_G = age;
  if (max_outdegree > 0)
    AS_CGB_BUBBLE_max_outdegree_G = max_outdegree;

  /* Mark the "interesting" fragments, assign them all relative coordinates,
     and make the graph into a DAG. */
  fprintf(BUB_LOG_G, "  * Step 1: Assign coordinates and mark fragments\n");
  AS_CGB_Bubble_dfs(bg);

  num_frags = GetNumFragments(BG_vertices(bg));
  top_order = (IntFragment_ID *)safe_calloc(sizeof(IntFragment_ID), num_frags);

  /* Get a topological ordering of the valid fragments. */
  fprintf(BUB_LOG_G, "  * Step 2: Topological sort of fragment graph\n");
  num_valid = AS_CGB_Bubble_topo_sort(bg, top_order);

#if AS_CGB_BUBBLE_VERY_VERBOSE
  for (f = 0; f < num_valid; ++f) 
    fprintf(BUB_LOG_G, "" F_IID " (" F_IID ")\n", top_order[f], 
	    get_iid_fragment(BG_vertices(bg), top_order[f]));
#endif

  BVS_sysInit();
  fwd = (BubVertexSet *)safe_malloc(sizeof(BubVertexSet) * num_frags);
  rvs = (BubVertexSet *)safe_malloc(sizeof(BubVertexSet) * num_frags);
  for (f = 0; f < num_frags; ++f) {
    BVS_initialize(&(fwd[f]));
    BVS_initialize(&(rvs[f]));
  }

  fprintf(BUB_LOG_G, "  * Step 3: Calculating fragment labels\n");
  fprintf(BUB_LOG_G, "  * Step 3: num_valid = " F_IID "\n", num_valid);
  if( num_valid != 0 ) {
    _forward_collect_sets(bg, fwd, top_order, num_valid);
    _reverse_collect_sets(bg, rvs, top_order, num_valid);
  }

#if AS_CGB_BUBBLE_VERY_VERBOSE
  for (f = 0; f < num_valid; ++f) {
    fprintf(BUB_LOG_G, "" F_IID " (" F_IID "):\t", top_order[f],
	    get_iid_fragment(BG_vertices(bg), top_order[f]));
    BVS_print(&(fwd[top_order[f]]), stderr);
    fprintf(BUB_LOG_G, "  |  ");
    BVS_print(&(rvs[top_order[f]]), stderr);
    fprintf(BUB_LOG_G, "\n");
  }
#endif
  
  fprintf(BUB_LOG_G, "  * Step 4: Finding matching labels\n");
  fprintf(BUB_LOG_G, "  * Step 4: num_valid = " F_IID "\n", num_valid);
  
  if( num_valid != 0 ) { 
    result = _collect_bubbles(bg, fwd, rvs, top_order, num_valid);
  }
  
  {
    AS_CGB_Bubble_List_t ptr = NULL;
    for (ptr = result; NULL != ptr; ptr = ptr->next) {
      assert(247319000 != ptr->start);
      assert(247319000 != ptr->start);
      ptr->start_sx = BG_vertexForward(bg, ptr->start);
      ptr->end_sx = !BG_vertexForward(bg, ptr->end);
    }
  }

  /* Deallocate storage. */
  for (f = 0; f < num_valid; ++f) {
    BVS_destroy(&(fwd[top_order[f]]));
    BVS_destroy(&(rvs[top_order[f]]));
  }

  BVS_sysDone();

  return result;
}


AS_CGB_Bubble_List_t
AS_CGB_Bubble_find_bubbles(Tfragment *frags, Tedge *edges, int sz, int age,
			   int max_outdegree)
{
  BubGraph bg = {0};
  AS_CGB_Bubble_List_t result = NULL;

  BUB_LOG_G = stderr;

  BG_initialize(&bg, frags, edges);
  result = AS_CGB_Bubble_find_bubbles_with_graph(&bg, sz, age, max_outdegree);
  BG_destroy(&bg);

  return result;
}


void
AS_CGB_Bubble_find_and_remove_bubbles
(FragStoreHandle TheFragStore,
 Tfragment *frags, Tedge *edges, 
 TChunkMesg *chunks, TChunkFrag *cfrgs,
 float gar,
 FILE *olap_file, FILE *log_file,
 const char * fileprefix)
{
  BubGraph bg;
  BubblePopper bp;
  AS_CGB_Bubble_List_t bubs = NULL;
  if (log_file)
    BUB_LOG_G = log_file;

  BG_initialize(&bg, frags, edges);
  bubs = AS_CGB_Bubble_find_bubbles_with_graph(&bg, 0, 0, 0);

#ifdef AS_CGB_BUBBLE_VERBOSE
  {
    int num_bubs = 0;
    AS_CGB_Bubble_List_t bptr = NULL;
    for (bptr = bubs; bptr; bptr = bptr->next)
      num_bubs++;
    fprintf(BUB_LOG_G, "  * SPECIAL PREVIEW: Found %d potential bubbles.\n",
	    num_bubs);
  }
#endif

  fprintf(BUB_LOG_G, "  * Processing bubbles.\n");
  BP_init(&bp, &bg, chunks, cfrgs, gar, TheFragStore, fileprefix);
  while (NULL != bubs) {
    AS_CGB_Bubble_List_t bptr = NULL;
    int num_ovl = 0;
    
    OverlapMesg *ovl = AS_CGB_Bubble_pop_bubble
      (&bp, bubs->start, bubs->start_sx,
       bubs->end, bubs->end_sx, &num_ovl);
    if (ovl) {
      int o;
      GenericMesg m;
      m.t = MESG_OVL;
      m.s = sizeof(OverlapMesg);
      for (o = 0; o < num_ovl; ++o) {
	m.m = (void *) &(ovl[o]);
	WriteProtoMesg_AS(olap_file, &m);
      }
    }
    assert(NULL != olap_file);
    fflush(olap_file); // CMM BUG WORK-AROUND

    bptr = bubs;
    bubs = bubs->next;
    safe_free(bptr);
    // release memory as we go.
  }
  // All the memory for bubs is released.

  fprintf(BUB_LOG_G, "  * ====================================================\n");
  fprintf(BUB_LOG_G, "  *                 BUBBLE POPPER STATS\n\n");
  fprintf(BUB_LOG_G, "  * Num Bubbles Processed:   \t\t%d\n", 
	  bp.numBubblesProcessed);
  if (bp.numBubblesProcessed > 0) {
    fprintf(BUB_LOG_G, "  * Num Bubbles Collapsed:   \t\t%d\n", 
	    bp.numBubblesCollapsed);
    fprintf(BUB_LOG_G, "  * Num Rejected by A-stat:  \t\t%d\n", 
	    bp.numRejectedByDiscriminator);
    fprintf(BUB_LOG_G, "  * Num Overlaps Attempted:  \t\t%d\n", 
	    bp.numOlapsComputed);
    fprintf(BUB_LOG_G, "  * Num Overlaps Found:      \t\t%d\n", 
	    bp.numOlapsSuccessful);
    fprintf(BUB_LOG_G, "  * Num Overlaps Output:     \t\t%d\n", 
	    bp.numOlapsRetained);
    fprintf(BUB_LOG_G, "  * Num Frags In Bubbles:    \t\t%d\n", 
	    bp.numFragsInBubbles);
    fprintf(BUB_LOG_G, "  * Num Frags In Collapsed:  \t\t%d\n", 
	    bp.numFragsInCollapsedBubbles);
    fprintf(BUB_LOG_G, "  * Average Bubble Length:   \t\t%d\n", 
	    bp.totalDistSpannedByBubbles / bp.numBubblesProcessed);
    fprintf(BUB_LOG_G, "  * Average Collapsed Length:\t\t%d\n", 
	    bp.totalDistSpannedByBubbles / bp.numBubblesProcessed);
  }
  fprintf(BUB_LOG_G, "  * ====================================================\n");
}
