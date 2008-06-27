
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

static char CM_ID[] = "$Id: AS_CGB_Bubble_Graph.c,v 1.7 2008-06-27 06:29:13 brianwalenz Exp $";

#include <stdlib.h>
#include "AS_global.h"
#include "AS_CGB_all.h"
#include "AS_CGB_methods.h"
#include "AS_CGB_Bubble_Graph.h"


//  AS_CGB_edgemate.c
IntEdge_ID find_overlap_edge_mate(const Tfragment frags[],
 const Tedge edges[],
 const IntEdge_ID ie0
);


void
BG_initialize(BubGraph_t bg, Tfragment *frags, Tedge *edges)
{
  IntFragment_ID num_v = GetNumFragments(frags);
  IntEdge_ID e, num_e = GetNumEdges(edges);

  memset(bg,0,sizeof(BubGraph));
  bg->e = edges;
  bg->v = frags;

  bg->vFlags = (uint16 *)safe_calloc(sizeof(uint16), num_v);
  bg->eFlags = (uint16 *)safe_calloc(sizeof(uint16), num_e);
  bg->vPos   = (int64 *)safe_calloc(sizeof(int64),  num_v);

  for (e = 0; e < num_e; ++e)
    switch (get_nes_edge(edges, e)) {

    case AS_CGB_INTERCHUNK_EDGE:
    case AS_CGB_INTRACHUNK_EDGE:
    case AS_CGB_TOUCHES_CONTAINED_EDGE:
      BG_E_setFlag(bg, e, AS_CGB_BUBBLE_E_DOVETAIL);
      BG_E_setFlagSymmetric(bg, e, AS_CGB_BUBBLE_E_VALID);
      break;

    case AS_CGB_DOVETAIL_EDGE:
    case AS_CGB_THICKEST_EDGE:
    case AS_CGB_BETWEEN_CONTAINED_EDGE:
      BG_E_setFlag(bg, e, AS_CGB_BUBBLE_E_DOVETAIL);
      // but otherwise not involved with the chunk graph
      break;

    case AS_CGB_CONTAINED_EDGE:
      BG_E_setFlag(bg, e, AS_CGB_BUBBLE_E_CONTAIN);
      BG_V_setFlag(bg, get_bvx_edge(edges, e), AS_CGB_BUBBLE_V_CONTAINED);
      break;

    case AS_CGB_TOUCHES_CRAPPY_DVT:
    case AS_CGB_BETWEEN_CRAPPY_DVT:
    case AS_CGB_TOUCHES_CRAPPY_CON:
    case AS_CGB_BETWEEN_CRAPPY_CON:
      /* no-op */
      break;

    default:
      fprintf(stderr,"Unrecognized edge label = %d\n", get_nes_edge(edges, e));
      assert(FALSE);
    }
}


void
BG_destroy(BubGraph_t bg)
{
  safe_free(bg->vFlags);
  safe_free(bg->eFlags);
  safe_free(bg->vPos);
}


Tfragment *
BG_vertices(BubGraph_t bg)
{
  return bg->v;
}


Tedge *
BG_edges(BubGraph_t bg)
{
  return bg->e;
}


IntFragment_ID
BG_getOppositeVertex(BubGraph_t bg, IntEdge_ID e, IntFragment_ID v)
{
  IntFragment_ID a = get_avx_edge(bg->e, e);
  if (a == v)
    return get_bvx_edge(bg->e, e);
  else
    return a;
}


int
BG_vertexForward(BubGraph_t bg, IntFragment_ID v)
{
  return !BG_V_isSetFlag(bg, v, AS_CGB_BUBBLE_V_REVERSED);
}


int
BG_Degree(BubGraph_t bg, IntFragment_ID v, uint16 flags)
{
  BG_E_Iter it;
  IntEdge_ID e;
  int result = 0;

  for (e = BGEI_bgn(bg, &it, v, bgeiBoth, flags);
       !BGEI_end(&it);
       e = BGEI_next(bg, &it, flags))
    result++;

  return result;
}


int
BG_outDegree(BubGraph_t bg, IntFragment_ID v, uint16 flags)
{
  BG_E_Iter it;
  IntEdge_ID e;
  int result = 0;

  for (e = BGEI_bgn(bg, &it, v, bgeiOut, flags);
       !BGEI_end(&it);
       e = BGEI_next(bg, &it, flags))
    result++;

  return result;
}


int
BG_inDegree(BubGraph_t bg, IntFragment_ID v, uint16 flags)
{
  BG_E_Iter it;
  IntEdge_ID e;
  int result = 0;

  for (e = BGEI_bgn(bg, &it, v, bgeiIn, flags);
       !BGEI_end(&it);
       e = BGEI_next(bg, &it, flags))
    result++;

  return result;
}


void
BG_V_setFlag(BubGraph_t bg, IntFragment_ID v, uint16 flag)
{
  bg->vFlags[v] |= flag;
}


void
BG_E_setFlag(BubGraph_t bg, IntEdge_ID e, uint16 flag)
{
  bg->eFlags[e] |= flag;
}


void
BG_E_setFlagSymmetric(BubGraph_t bg, IntEdge_ID e, uint16 flag)
{
  IntEdge_ID shadow;

  bg->eFlags[e] |= flag;
  shadow = find_overlap_edge_mate(bg->v, bg->e, e);
  if (AS_CGB_EDGE_NOT_FOUND != shadow) {
    bg->eFlags[shadow] |= flag;
  }
}


void
BG_V_clearFlag(BubGraph_t bg, IntFragment_ID v, uint16 flag)
{
  bg->vFlags[v] &= ~flag;
}


void
BG_E_clearFlag(BubGraph_t bg, IntEdge_ID e, uint16 flag)
{
  bg->eFlags[e] &= ~flag;
}


void
BG_E_clearFlagSymmetric(BubGraph_t bg, IntEdge_ID e, uint16 flag)
{
  IntEdge_ID shadow;

  bg->eFlags[e] &= ~flag;
  shadow = find_overlap_edge_mate(bg->v, bg->e, e);
  if(AS_CGB_EDGE_NOT_FOUND != shadow) {
    bg->eFlags[shadow] &= ~flag;
  }
}


uint16
BG_V_isSetFlag(BubGraph_t bg, IntFragment_ID v, uint16 flag)
{
  assert(NULL != bg);
  return bg->vFlags[v] & flag;
}


uint16
BG_E_isSetFlag(BubGraph_t bg, IntEdge_ID e, uint16 flag)
{
  return bg->eFlags[e] & flag;
}


uint16
BG_V_getFlags(BubGraph_t bg, IntFragment_ID v)
{
  return bg->vFlags[v];
}


uint16
BG_E_getFlags(BubGraph_t bg, IntEdge_ID e)
{
  return bg->eFlags[e];
}


void
BG_V_setDistance(BubGraph_t bg, IntFragment_ID v, int64 p)
{
  bg->vPos[v] = p;
}


int64
BG_V_getDistance(BubGraph_t bg, IntFragment_ID v)
{
  return bg->vPos[v];
}


/*
 * Edge iterator methods
 */

/* Advances to the first iterator position satisfying the flag condition.
   Will not advance if the current position is valid. */
IntEdge_ID
BGEI__advance(BubGraph_t bg, BG_E_Iter_t it, uint16 flags)
{
  int forward;

  if (flags != AS_CGB_BUBBLE_E_ALL)
    while ((it->cur != it->end) &&
	   ((BG_E_getFlags(bg, it->cur) & flags) != flags))
      (it->cur)++;
  else
    (it->cur)++;

  if ((it->cur == it->end) && (it->type == bgeiBoth)) {
    forward = BG_vertexForward(bg, it->v);
    it->type = bgeiIn;
    if (forward) {
      it->cur = get_segstart_vertex(BG_vertices(bg), it->v, 0);
      it->end = it->cur + get_seglen_vertex(BG_vertices(bg), it->v, 0);
    }
    else {
      it->cur = get_segstart_vertex(BG_vertices(bg), it->v, 1);
      it->end = it->cur + get_seglen_vertex(BG_vertices(bg), it->v, 1);
    }
    return BGEI__advance(bg, it, flags);
  }

  return it->cur;
}


IntEdge_ID
BGEI_bgn(BubGraph_t bg, BG_E_Iter_t it, IntFragment_ID v, BG_E_IterType t,
	 uint16 flags)
{
  int forward = BG_vertexForward(bg, v);

  it->type = t;
  it->v = v;

  if (((t == bgeiIn) && forward) ||
      (((t == bgeiOut) || (t == bgeiBoth)) && !forward)) {
    it->cur = get_segstart_vertex (BG_vertices(bg), v, 0);
    it->end = it->cur + get_seglen_vertex(BG_vertices(bg), v, 0);
  }
  else {
    it->cur = get_segstart_vertex(BG_vertices(bg), v, 1);
    it->end = it->cur + get_seglen_vertex(BG_vertices(bg), v, 1);
  }

  return BGEI__advance(bg, it, flags);
}


IntEdge_ID
BGEI_next(BubGraph_t bg, BG_E_Iter_t it, uint16 flags)
{
  if (it->cur == it->end)
    return it->cur;
  else {
    (it->cur)++;
    return BGEI__advance(bg, it, flags);
  }
}


int
BGEI_end(BG_E_Iter_t it)
{
  return ((it->type != bgeiBoth) && (it->cur == it->end));
}


IntEdge_ID
BGEI_cur(BG_E_Iter_t it)
{
  return it->cur;
}
