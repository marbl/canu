
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

static char *rcsid = "$Id: AS_CGB_Bubble_dfs.c,v 1.9 2008-10-08 22:02:54 brianwalenz Exp $";

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "AS_CGB_all.h"
#include "AS_CGB_Bubble.h"
#include "AS_CGB_Bubble_Graph.h"


static
int64 _get_dist_from_edge(BubGraph_t bg, IntFragment_ID src, IntEdge_ID e)
{
  IntEdge_ID dst;
  int64 d;
  int s_sx, d_sx, s_fwd, d_fwd;
  int walk_is_forward;

  /* Assumes dovetail edges only. */
  assert(BG_E_isSetFlag(bg, e, AS_CGB_BUBBLE_E_DOVETAIL));

  /* Assumes that the source fragment is always the "a" fragment */
  assert(get_avx_edge(BG_edges(bg), e) == src);

  dst = get_bvx_edge(BG_edges(bg), e);

  s_sx = get_asx_edge(BG_edges(bg), e);
  d_sx = get_bsx_edge(BG_edges(bg), e);

  s_fwd = BG_vertexForward(bg, src);
  d_fwd = BG_vertexForward(bg, dst);

  walk_is_forward = (s_fwd ^ !s_sx);

  if (walk_is_forward)
    d = BG_V_getDistance(bg, src) + get_ahg_edge(BG_edges(bg), e);
  else
    d = BG_V_getDistance(bg, src) - get_bhg_edge(BG_edges(bg), e);

  return d;
}


static
int
_dst_edge_is_reversed(int cur_rev, int cur_sx, int dst_sx)
{
	/*
	  cur_rev	cur_sx		dst_sx		dst_rev
	  0		0		0		1
	  0		0 		1		0
	  0		1		0		0
	  0		1		1		1
	  1		0		0		0
	  1		0		1		1
	  1		1		0		1
	  1		1		1		0
	*/

  return !((!cur_rev ^ !cur_sx) ^ dst_sx);
}


static
void
_set_dist_from_edge(BubGraph_t bg, IntFragment_ID src, IntEdge_ID e)
{
  IntFragment_ID dst = BG_getOppositeVertex(bg, e, src);

  BG_V_setDistance(bg, dst, _get_dist_from_edge(bg, src, e));
}


static
int
_back_edge_consistent(BubGraph_t bg, IntEdge_ID cur_e, IntFragment_ID cur_v,
		      IntFragment_ID dst_v)
{
  int cur_sx, dst_sx, cur_fwd;
  int dst_rev_by_edge, dst_rev_by_mark;
  int dir_right;
  int64 cur_pos, dst_pos;

  if (get_avx_edge(BG_edges(bg), cur_e) == cur_v) {
    cur_sx = get_asx_edge(BG_edges(bg), cur_e);
    dst_sx = get_bsx_edge(BG_edges(bg), cur_e);
  }
  else {
    cur_sx = get_bsx_edge(BG_edges(bg), cur_e);
    dst_sx = get_asx_edge(BG_edges(bg), cur_e);
  }

  /* First, check that orientations are consistent. */

  dst_rev_by_edge = _dst_edge_is_reversed(!BG_vertexForward(bg,cur_v),
					  cur_sx, dst_sx);
  dst_rev_by_mark = !BG_vertexForward(bg, dst_v);

  if (dst_rev_by_edge != dst_rev_by_mark)
    return FALSE;

  /* Now check that relative position is correct. */
  cur_fwd = BG_vertexForward(bg, cur_v);
  dir_right = ((cur_fwd && cur_sx) || (!cur_fwd && !cur_sx));
  cur_pos = BG_V_getDistance(bg, cur_v);
  dst_pos = BG_V_getDistance(bg, dst_v);

  if ((dir_right && (cur_pos > dst_pos)) ||
      (!dir_right && (cur_pos < dst_pos)))
    return FALSE;

  /* NOTE: Could insert distance consistency check here.  But
     what should the criteria be? */

  return TRUE;
}


static
void
_do_dfs(BubGraph_t bg, IntFragment_ID start_v, BG_E_Iter stack[])
{
  int32 s_top = -1;
  IntFragment_ID cur_v, dst_v;
  IntEdge_ID cur_e;
  int cur_sx, dst_sx;
  BG_E_Iter_t cur_v_it;
  uint16 e_flags = AS_CGB_BUBBLE_E_DOVETAIL | AS_CGB_BUBBLE_E_UNUSED |
    AS_CGB_BUBBLE_E_VALID;

  /* Put start vertex onto the stack, and here we go ...*/
  cur_v = start_v;
  cur_v_it = &(stack[++s_top]);
  BGEI_bgn(bg, cur_v_it, start_v, bgeiBoth, e_flags);
  BG_V_setFlag(bg, start_v, AS_CGB_BUBBLE_V_STACKED);
  BG_V_setDistance(bg, start_v, 0);

  while (s_top > -1) {
    if (BGEI_end(cur_v_it)) {
#if AS_CGB_BUBBLE_VERY_VERBOSE
      fprintf(BUB_LOG_G, "Done with " F_IID " (" F_IID ").  Backtracking.\n", cur_v,
	      get_iid_fragment(BG_vertices(bg), cur_v));
#endif
      BG_V_clearFlag(bg, cur_v, AS_CGB_BUBBLE_V_STACKED);
      BG_V_setFlag(bg, cur_v, AS_CGB_BUBBLE_V_DONE);
      cur_v_it = &(stack[--s_top]);
      if (s_top != -1) {
	cur_v = cur_v_it->v;
	BGEI_next(bg, cur_v_it, e_flags);
      }
    }
    else {
      cur_e = BGEI_cur(cur_v_it);
      dst_v = BG_getOppositeVertex(bg, cur_e, cur_v);
#if AS_CGB_BUBBLE_VERY_VERBOSE
      fprintf(BUB_LOG_G, "Processing edge from " F_IID " (" F_IID ") to " F_IID " (" F_IID ").  ", cur_v,
	      get_iid_fragment(BG_vertices(bg), cur_v), dst_v,
	      get_iid_fragment(BG_vertices(bg), dst_v));
#endif
      BG_E_clearFlagSymmetric(bg, cur_e, AS_CGB_BUBBLE_E_UNUSED);

      if (BG_V_isSetFlag(bg, dst_v, AS_CGB_BUBBLE_V_STACKED)) {
	/* See if the edge is truly inconsistent, or if it can
	   be oriented in a consistent manner.  */
	if (!_back_edge_consistent(bg, cur_e, cur_v, dst_v)) {
	  /* Bad back edge. Mark edge as invalid. */
	  BG_E_clearFlagSymmetric(bg, cur_e, AS_CGB_BUBBLE_E_VALID);
#if AS_CGB_BUBBLE_VERY_VERBOSE
	  fprintf(BUB_LOG_G, "INCONSISTENT BACK EDGE.\n");
#endif
	}
#if AS_CGB_BUBBLE_VERY_VERBOSE
	else
	  fprintf(BUB_LOG_G, "Consistent Back Edge.\n");
#endif

	BGEI_next(bg, cur_v_it, e_flags);
      }
      else { /* tree edge - no forward or cross edges in undirected search */

	/* Figure out orientation of new vertex. */
	cur_sx = get_asx_edge(BG_edges(bg), cur_e);
	dst_sx = get_bsx_edge(BG_edges(bg), cur_e);

	if (_dst_edge_is_reversed(!BG_vertexForward(bg,cur_v), cur_sx, dst_sx))
	  BG_V_setFlag(bg, dst_v, AS_CGB_BUBBLE_V_REVERSED);

	/* Set position of new vertex. */
	_set_dist_from_edge(bg, cur_v, cur_e);

	/* Place new vertex on the stack top and make it the current vertex. */
	cur_v = dst_v;
	cur_v_it = &(stack[++s_top]);
	BGEI_bgn(bg, cur_v_it, cur_v, bgeiBoth, e_flags);
	BG_V_setFlag(bg, cur_v, AS_CGB_BUBBLE_V_STACKED);
#if AS_CGB_BUBBLE_VERY_VERBOSE
	if (!BG_V_isSetFlag(bg, cur_v, AS_CGB_BUBBLE_V_CONTAINED))
	  fprintf(BUB_LOG_G, "\nGoing to " F_IID "\t ( iid " F_IID ", dist " F_S64 ", forward = %d, ssx = %d, dsx = %d )\n",
		  cur_v,
		  get_iid_fragment(BG_vertices(bg), cur_v),
		  BG_V_getDistance(bg, cur_v),
		  BG_vertexForward(bg, cur_v),
		  cur_sx, dst_sx);
	else
	  fprintf(BUB_LOG_G, "\nGoing to " F_IID "(C)\t ( iid " F_IID ", dist " F_S64 ", forward = %d, ssx = %d, dsx = %d )\n",
		  cur_v,
		  get_iid_fragment(BG_vertices(bg), cur_v),
		  BG_V_getDistance(bg, cur_v),
		  BG_vertexForward(bg, cur_v),
		  cur_sx, dst_sx);
#endif
      }
    }
  }
}


void
AS_CGB_Bubble_dfs(BubGraph_t bg)
{
  BG_E_Iter *stack = NULL;
  IntFragment_ID num_v = GetNumFragments(BG_vertices(bg));
  IntFragment_ID v;
  IntEdge_ID num_e = GetNumEdges(BG_edges(bg));
  IntEdge_ID e;

  /* Mark all fragments as unvisited and edges as unused. */
  for (v = 0; v < num_v; v++)
    BG_V_setFlag(bg, v, AS_CGB_BUBBLE_V_NEW);

  for (e = 0; e < num_e; e++)
    BG_E_setFlag(bg, e, AS_CGB_BUBBLE_E_UNUSED);

  /* Allocate stack space */
  stack = safe_malloc(sizeof(BG_E_Iter) * num_v);

  for (v = 0; v < num_v; v++)
    if ((!BG_V_isSetFlag(bg, v, AS_CGB_BUBBLE_V_DONE)) &&
	(!BG_V_isSetFlag(bg, v, AS_CGB_BUBBLE_V_CONTAINED))) {
#if AS_CGB_BUBBLE_VERY_VERBOSE
      fprintf(BUB_LOG_G, "START: " F_IID " (" F_IID ")\n", v,
	      get_iid_fragment(BG_vertices(bg), v));
#endif
      _do_dfs(bg, v, stack);
    }

  safe_free(stack);
}
