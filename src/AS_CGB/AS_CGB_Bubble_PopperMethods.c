
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
= "$Id: AS_CGB_Bubble_PopperMethods.c,v 1.1.1.1 2004-04-14 13:49:51 catmandew Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "AS_global.h"
#include "AS_CGB_all.h"
#include "AS_CGB_methods.h"
#include "AS_CGB_Bubble_Graph.h"
#include "AS_CGB_Bubble.h"
#include "AS_CGB_Bubble_log.h"
#include "AS_CGB_Bubble_Popper.h"
#include "AS_CGB_Bubble_PopperMethods.h"


int
BP__insertFragment(BubblePopper_t bp, IntFragment_ID v)
{
  int p;
  int64 v_pos = BG_V_getDistance(bp->bg, v);

  if (bp->curBubSize >= POPPER_MAX_BUBBLE_SIZE)
    return FALSE;

  p = bp->curBubSize;
  while ((p > 0) && BG_V_getDistance(bp->bg, bp->bubFrags[p - 1]) > v_pos) {
    bp->bubFrags[p] = bp->bubFrags[p - 1];
    bp->vidToBid[bp->bubFrags[p]] = p;
    p--;
  }

  bp->bubFrags[p] = v;
  bp->vidToBid[v] = p;
  bp->curBubSize++;

  return TRUE;
}


int
BP_find_bubble_dfs(BubblePopper_t bp, IntFragment_ID start, 
		   IntFragment_ID end)
{
  int32 s_top = -1;
  IntFragment_ID cur_v, dst_v;
  IntEdge_ID cur_e;
  BG_E_Iter_t cur_v_it;
  uint16 e_flags = AS_CGB_BUBBLE_E_DOVETAIL | AS_CGB_BUBBLE_E_UNUSED | 
    AS_CGB_BUBBLE_E_VALID;

  /* Put start vertex onto the stack, and here we go ...*/
  cur_v = start;
  cur_v_it = &(bp->dfsStack[++s_top]);
  BGEI_bgn(bp->bg, cur_v_it, start, bgeiOut, e_flags);
  BG_V_setFlag(bp->bg, start, AS_CGB_BUBBLE_V_STACKED);
  BP__insertFragment(bp, start);

  while (s_top > -1) {
    if (BGEI_end(cur_v_it) || (cur_v == end)) {
#if AS_CGB_BUBBLE_VERY_VERBOSE
      fprintf(BUB_LOG_G, "Done with " F_IID " (" F_IID ").  Backtracking.\n", cur_v,
	      get_iid_fragment(BG_vertices(bp->bg), cur_v));
#endif
      BG_V_clearFlag(bp->bg, cur_v, AS_CGB_BUBBLE_V_STACKED);
      BG_V_setFlag(bp->bg, cur_v, AS_CGB_BUBBLE_V_IN_BUBBLE);
      cur_v_it = &(bp->dfsStack[--s_top]);
      if (s_top != -1) {
	cur_v = cur_v_it->v;
	BGEI_next(bp->bg, cur_v_it, e_flags);
      }
    } 

    else {
      cur_e = BGEI_cur(cur_v_it);
      dst_v = BG_getOppositeVertex(bp->bg, cur_e, cur_v);
#if AS_CGB_BUBBLE_VERY_VERBOSE
      fprintf(BUB_LOG_G, "Processing edge from " F_IID " (" F_IID ") to " F_IID " (" F_IID ").  ", cur_v, 
	      get_iid_fragment(BG_vertices(bp->bg), cur_v), dst_v,
	      get_iid_fragment(BG_vertices(bp->bg), dst_v));
#endif
      BG_E_clearFlagSymmetric(bp->bg, cur_e, AS_CGB_BUBBLE_E_UNUSED); 
      BG_E_setFlagSymmetric(bp->bg, cur_e, AS_CGB_BUBBLE_E_IN_BUBBLE);
      
      if (BG_V_isSetFlag(bp->bg, dst_v, AS_CGB_BUBBLE_V_STACKED) ||
	  BG_V_isSetFlag(bp->bg, dst_v, AS_CGB_BUBBLE_V_IN_BUBBLE)) {
#if AS_CGB_BUBBLE_VERY_VERBOSE
	fprintf(BUB_LOG_G, "Back Edge (or edge to end fragment).\n");
#endif
	BGEI_next(bp->bg, cur_v_it, e_flags);
      }
      else {
	/* Place new vertex on the stack top and make it the current vertex. */
	cur_v = dst_v;
	if (!BP__insertFragment(bp, cur_v)) {
#if AS_CGB_BUBBLE_VERBOSE
	  fprintf(BUB_LOG_G, "Bubble too big (maybe it's not closed?).  Aborting.\n");
#endif
	  return FALSE;
	}

	cur_v_it = &(bp->dfsStack[++s_top]);
	BGEI_bgn(bp->bg, cur_v_it, cur_v, bgeiOut, e_flags);
	BG_V_setFlag(bp->bg, cur_v, AS_CGB_BUBBLE_V_STACKED);
#if AS_CGB_BUBBLE_VERY_VERBOSE
	if (!BG_V_isSetFlag(bp->bg, cur_v, AS_CGB_BUBBLE_V_CONTAINED))
	  fprintf(BUB_LOG_G, "\nGoing to " F_IID "\t ( iid " F_IID ", dist " F_S64 ", forward = %d )\n", 
		  cur_v, 
		  get_iid_fragment(BG_vertices(bp->bg), cur_v),
		  BG_V_getDistance(bp->bg, cur_v), 
		  BG_vertexForward(bp->bg, cur_v));
	else
	  fprintf(BUB_LOG_G, "\nGoing to " F_IID "(C)\t ( iid " F_IID ", dist " F_S64 ", forward = %d )\n", 
		  cur_v, 
		  get_iid_fragment(BG_vertices(bp->bg), cur_v),
		  BG_V_getDistance(bp->bg, cur_v), 
		  BG_vertexForward(bp->bg, cur_v));
#endif
      }
    }
  }

  return TRUE;
}


void
BP_transitive_closure(BubblePopper_t bp)
{
  int src_c, dst_c, src_r, dst_r;

  for (dst_c = 1; dst_c < BP_numFrags(bp); ++dst_c) {
    for (dst_r = 0; dst_r < dst_c; ++dst_r)
      if (BP_getAdj(bp, dst_r, dst_c)) {
	src_c = dst_r;
	for (src_r = 0; src_r <= src_c; ++src_r)
	  if (BP_getAdj(bp, src_r, src_c) && !(BP_getAdj(bp, src_r, dst_c)))
	    BP_setAdj(bp, src_r, dst_c, 2);
      }
  }
}


int
BP_DAG_longest_path(BubblePopper_t bp)
{
  int q_start = 0, q_end = 0, adj;
  int r, c, cur_v;
  int dist;
  int num_frags = BP_numFrags(bp);

  /* Warning: Hack.  Use the vertex field in the dfsStack space to hold
     the topological sort queue.  No need to allocate more space that way. 
     Also, use the distance array to hold the indegrees of the fragments. */

  memset(bp->topDistArray, 0, POPPER_MAX_BUBBLE_SIZE * sizeof(int));
  for (c = 0; c < num_frags; ++c)
    for (r = 0; r < num_frags; ++r)
      if ((r != c) && (BP_getAdj(bp, r, c) != 0))
	bp->topDistArray[c]++;

  if (bp->topDistArray[0] > 0) {
#if AS_CGB_BUBBLE_VERBOSE
    fprintf(BUB_LOG_G, "Start vertex has non-zero in degree!  Aborting.\n");
#endif
    return FALSE;
  }

#if AS_CGB_BUBBLE_VERY_VERBOSE
  fprintf(BUB_LOG_G, "Adding vertex 0 (" F_IID ") as start.\n",
	  get_iid_fragment(BG_vertices(bp->bg), BP_getFrag(bp, 0)));
#endif
  bp->dfsStack[q_end++].v = 0;	/* Assumes bubble start is one and only
				   fragment with indegree 0 (should be). */
  
  while (q_start < q_end) {
#if AS_CGB_BUBBLE_VERY_VERBOSE
    fprintf(BUB_LOG_G, "PROCESSING vertex " F_IID " (" F_IID ").\n", bp->dfsStack[q_start].v,
	    get_iid_fragment(BG_vertices(bp->bg), bp->dfsStack[q_start].v));
#endif

    cur_v = bp->dfsStack[q_start].v;
    for (c = 0; c < num_frags; ++c) {
      if ((cur_v != c) && ((adj = BP_getAdj(bp, cur_v, c)) > 0)) 
	if ((--(bp->topDistArray[c])) == 0) {
	  dist = -1;
	  for (r = 0; r < num_frags; ++r)
	    if ((r != c) && (BP_getAdj(bp, r, c) > 0) && 
		(bp->topDistArray[r] > dist)) 
	      dist = bp->topDistArray[r];
	
	  bp->topDistArray[c] = dist + 1;
	  bp->dfsStack[q_end++].v = c;
	
#if AS_CGB_BUBBLE_VERY_VERBOSE	
	  fprintf(BUB_LOG_G, "Adding vertex " F_IID " (" F_IID ") at distance %d.\n", c,
		  get_iid_fragment(BG_vertices(bp->bg), BP_getFrag(bp, c)),
		  bp->topDistArray[c]);
#endif
	}
    }
    
    q_start++;
  }
  
  if (q_end < BP_numFrags(bp)) {
    fprintf(BUB_LOG_G, "WARNING: Only processed " F_IID " of " F_IID " vertices!  Cyclic graph!\n", q_end, BP_numFrags(bp));
    return 0;
  }

  dist = -1;

  for (r = 0; r < BP_numFrags(bp); ++r) {
    if (bp->topDistArray[r] > dist)
      dist = bp->topDistArray[r];
  }

  return dist;
}


float
BP_discriminator(BubblePopper_t bp)
{
  IntFragment_ID start_bid, end_bid, i;
  IntChunk_ID start_c, end_c;
  BPTYPE total_len;
  int num_rand_frags;
  FragType type;

  start_bid = BP_getFrag(bp, 0);  
  end_bid = BP_getFrag(bp, BP_numFrags(bp) - 1);
  start_c = get_cid_fragment(BG_vertices(bp->bg), start_bid);
  end_c = get_cid_fragment(BG_vertices(bp->bg), end_bid);
  
  total_len = BP_getChunk(bp, start_c)->rho + BP_getChunk(bp, end_c)->rho +
    (BG_V_getDistance(bp->bg, end_bid) - BG_V_getDistance(bp->bg, start_bid));

  num_rand_frags = count_the_randomly_sampled_fragments_in_a_chunk(
                     BG_vertices(bp->bg),
		     BP_chunkFrags(bp),
		     BP_chunks(bp),
		     start_c);

  num_rand_frags += count_the_randomly_sampled_fragments_in_a_chunk(
                      BG_vertices(bp->bg),
		      BP_chunkFrags(bp),
		      BP_chunks(bp),
		      end_c);

  for (i = 0; i < BP_numFrags(bp); ++i) {
    type = get_typ_fragment(BG_vertices(bp->bg), BP_getFrag(bp, i));
    if ((type == AS_READ) || (type == AS_EBAC) || (type == AS_EXTR))
      num_rand_frags++;
  }

#if AS_CGB_BUBBLE_VERY_VERBOSE
  fprintf(BUB_LOG_G, "Computing discriminator with length = " BPFORMAT " and num frags = %d.\n", total_len, num_rand_frags);
#endif

  return compute_coverage_statistic(total_len, num_rand_frags,
				    bp->globalArrivalRate);
}

