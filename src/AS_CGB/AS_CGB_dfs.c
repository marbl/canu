
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
= "$Id: AS_CGB_dfs.c,v 1.4 2005-03-22 19:48:27 jason_miller Exp $";
/*********************************************************************
 *
 * Module: AS_CGB_dfs.c
 *
 * Description: Depth First Search Graph traversal for the Chunk Graph
 * Builder.  This module takes the current marked graph and returns a
 * best ranking of the verticies in fragment_rank[].
 *
 * Assumptions:
 *
 * 1. Overview
 * 
 * 2. Memory Usage
 * 
 * 3. Interface
 * 
 * 4. Design
 * 
 * 5. Limitations
 * 
 * 6. Status
 * 
 * Author: Clark Mobarry
 ***********************************************************************/

/*********************************************************************/
/* System include files */
#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/times.h>

/*************************************************************************/
/* Local include files */
#include "AS_global.h"
#include "AS_MSG_pmesg.h"
#include "AS_CGB_all.h"
//#include "AS_CGB_traversal.h"

#undef INTRACHUNK_FIRST
#undef PROCESS_CONTAINED_FRAGS_ALSO             
#define FRAGMENT_SUBSET

/*************************************************************************/
/* Constant definitions; Macro definitions; type definitions */
#define AS_DFS_NOT_RANKED     (CDS_UINT32_MAX)
#define AS_DFS_ROOT           (CDS_UINT32_MAX-1)
#define COORDINATE_SLOP 40           

/*************************************************************************/
/* Static Globals */

/*************************************************************************/
/* Function prototypes for internal static functions */

/*************************************************************************/


typedef struct {
  IntFragment_ID iid;
} DFSstackrec;

typedef struct {
  IntRank next;
  IntRank size;
  DFSstackrec *memory;
} DFSstackobj;

static DFSstackobj * dfs_stack_create(IntRank nfrag) {
  DFSstackobj *  self = NULL;
  SAFE_MALLOC( self, DFSstackobj, 1);
  self->next = 0;
  self->size = nfrag;
  self->memory = NULL;
  SAFE_MALLOC(self->memory, DFSstackrec, nfrag);
  return(self);
}

static void dfs_stack_destroy(DFSstackobj *self) {
  self->next = 0;
  self->size = 0;
  SAFE_FREE(self->memory);
  SAFE_FREE(self);
}

static void dfs_stack_push(DFSstackobj *self, DFSstackrec *x) {
  // assert((self->next) >= 0);
  assert((self->next) < (self->size) - 1 );
  (self->memory)[(self->next)++] = *x;
}

static void dfs_stack_pop(DFSstackobj *self, DFSstackrec *x) {
  assert((self->next) > 0);
  assert((self->next) < (self->size) );
  *x = (self->memory)[--(self->next)];
}

#if 0
static void dfs_stack_print_to(DFSstackobj *self, IntRank ir, FILE *fout) {
  IntFragment_ID ix;
  DFSstackrec rx;
  //#if 0
  assert((self->next) > 0);
  assert((self->next) < (self->size) );
  //#endif  
  for(ix = (self->next); ix > 0 ; ix--) {
    rx = (self->memory)[--(self->next)];
    if( rx.iid == ir ) {
      break;
    }
  }
}
#endif

typedef enum {
  AS_DFS_WHITE,
  AS_DFS_GRAY,
  AS_DFS_BLACK
} DFScolor;


static void found_a_fragment
(
 const IntFragment_ID vid,
 IntRank    *nfound,
 IntRank    *timestamp,
 IntRank    fragment_rank[],
 IntRank    dfs_discovered[],
 IntRank    dfs_finished[],
 DFScolor   dfs_color[]
)
{
  // Convert a WHITE vertex to GRAY. 
  assert(AS_DFS_WHITE == dfs_color[vid]);
  dfs_color[vid] = AS_DFS_GRAY;
  assert(fragment_rank[vid] == AS_DFS_NOT_RANKED);
  fragment_rank[vid] = (*nfound)++;
  assert(dfs_discovered[vid] == AS_DFS_NOT_RANKED);
  dfs_discovered[vid] = (*timestamp)++;
  assert(dfs_finished[vid] == AS_DFS_NOT_RANKED);

}

static void finished_a_fragment(
 const IntFragment_ID vid,
 IntRank    *nfound,
 IntRank    *timestamp,
 IntRank    fragment_rank[],
 IntRank    dfs_discovered[],
 IntRank    dfs_finished[],
 DFScolor   dfs_color[]
)
{
  // Convert a GRAY vertex to BLACK.
  assert(AS_DFS_GRAY == dfs_color[vid]);
  dfs_color[vid] = AS_DFS_BLACK;
  assert(fragment_rank[vid] != AS_DFS_NOT_RANKED);
  assert(dfs_discovered[vid] != AS_DFS_NOT_RANKED);
  assert(dfs_finished[vid] == AS_DFS_NOT_RANKED);
  dfs_finished[vid] = (*timestamp)++;
}


static void as_dfs_recursive
(
 FILE *    fout,
 IntFragment_ID iv0,
 int       is0,
 IntRank   *nfound,
 IntRank   *timestamp,
 IntFragment_ID nfrag, 
 Tfragment  frags[],
 IntEdge_ID nedge, 
 Tedge      edges[], 
 IntRank    dfs_predecessor[],
 IntRank    fragment_rank[],
 IntRank    dfs_discovered[],
 IntRank    dfs_finished[],
 DFScolor   dfs_color[],
 int fragment_direction_is_forward_in_component,
 int sum_of_ahg_in_component,
 int sum_of_bhg_in_component,
 DFSstackobj * dfs_stack
);

static void process_edge
(
 FILE *    fout,
 const IntFragment_ID iv0,
 const int            is0,
 const IntEdge_ID     ie1,
 IntRank   *nfound,
 IntRank   *timestamp,
 IntFragment_ID nfrag, 
 Tfragment  frags[],
 IntEdge_ID nedge, 
 Tedge      edges[], 
 IntRank    dfs_predecessor[],
 IntRank    fragment_rank[],
 IntRank    dfs_discovered[],
 IntRank    dfs_finished[],
 DFScolor   dfs_color[],
 int fragment_direction_is_forward_in_component,
 int sum_of_ahg_in_component,
 int sum_of_bhg_in_component,
 DFSstackobj * dfs_stack
)
{
  const IntFragment_ID iv1 = get_bvx_edge(edges,ie1);
  const int is1 = get_bsx_edge(edges,ie1);

  const int avx = get_avx_edge(edges,ie1);
  const int bvx = get_bvx_edge(edges,ie1);
  const int asx = get_asx_edge(edges,ie1);
  const int bsx = get_bsx_edge(edges,ie1);
  const int ahg = get_ahg_edge(edges,ie1);
  const int bhg = get_bhg_edge(edges,ie1);

  const Tlab lab = get_lab_fragment(frags,iv1);
  
  
  // When following an edge, the new direction flags are:

  const int walk_direction_is_forward_in_component 
    = ! (fragment_direction_is_forward_in_component ^ asx);

  fragment_direction_is_forward_in_component
  = (walk_direction_is_forward_in_component ^ bsx);

  sum_of_ahg_in_component += ( walk_direction_is_forward_in_component
                               ? ahg : - bhg);
  sum_of_bhg_in_component += ( walk_direction_is_forward_in_component
                               ? bhg : - ahg);
  
  assert( iv0 == avx );
  assert( iv1 == bvx );
  assert( is0 == asx );
  assert( is1 == bsx );

  if(
     ((AS_CGB_INTERCHUNK_FRAG == lab) ||
      (AS_CGB_INTRACHUNK_FRAG == lab) ||
      (AS_CGB_THRU_FRAG == lab) ||
      (AS_CGB_HANGING_CHUNK_FRAG == lab) ||
      (AS_CGB_BRANCHMULTICONT_FRAG == lab))
     ) {

  // DFScolor the_color = dfs_color[iv1];
  switch(dfs_color[iv1]) {
  case AS_DFS_WHITE: // Tree edge
    { /* iv1 is unexplored */
      
      // Note that the only necessary stack variables are iv0,
      // is0, in0.
      dfs_predecessor[iv1] = iv0;
      set_o5p_fragment(frags,bvx,sum_of_ahg_in_component);
      set_o3p_fragment(frags,bvx,sum_of_bhg_in_component);
      as_dfs_recursive
        (fout,
         iv1, is1,
         nfound,
         timestamp,
         nfrag, frags, nedge, edges, 
         dfs_predecessor, fragment_rank,
         dfs_discovered, dfs_finished, dfs_color,
         fragment_direction_is_forward_in_component,
         sum_of_ahg_in_component,
         sum_of_bhg_in_component,
         dfs_stack
         );
    }
    break;
  case AS_DFS_GRAY: // Back edge
    {
      const Tnes nes = get_nes_edge(edges,ie1);
      const int sum_of_ahg_in_component1 = get_o5p_fragment(frags,iv1);
      const int sum_of_bhg_in_component1 = get_o3p_fragment(frags,iv1);
      
      if( iv1 == dfs_predecessor[iv0] ) break;
      // The back edge to the predecessor is irrelevant.
      
#if 0
      fprintf(fout,"BACK EDGE: (" F_IID ",%d)-(" F_IID ",%d) %d %d %d %d %d\n",
              get_iid_fragment(frags,iv0),is0,
              get_iid_fragment(frags,iv1),is1,
              nes,
              sum_of_ahg_in_component,
              sum_of_bhg_in_component,
              sum_of_ahg_in_component1,
              sum_of_bhg_in_component1
              );
#endif      
      {
        //  the back-edge is consistent with the component-coordinates 
        if(
           (abs(sum_of_ahg_in_component - sum_of_ahg_in_component1) < COORDINATE_SLOP) &&
           (abs(sum_of_bhg_in_component - sum_of_bhg_in_component1) < COORDINATE_SLOP)
           ) {
          fprintf(fout,"CONSISTENT BACK EDGE: (" F_IID ",%d)-(" F_IID ",%d) %d %d %d %d %d\n",
              get_iid_fragment(frags,iv0),is0,
              get_iid_fragment(frags,iv1),is1,
              nes,
              sum_of_ahg_in_component,
              sum_of_bhg_in_component,
              sum_of_ahg_in_component1,
              sum_of_bhg_in_component1
              );
          // We want to print out the cycle of fragments here.
#if 0
          set_nes_edge(edges,ie1,AS_CGB_MARKED_BY_BRANCH_DVT);
          fix_overlap_edge_mate(frags,edges,ie1);
#endif          
        }
      }
#if 0
      set_nes_edge(edges,ie1,AS_CGB_MARKED_BY_BRANCH_DVT);
      fix_overlap_edge_mate(frags,edges,ie1);
#endif          
    }
    break;
  case AS_DFS_BLACK: // Forward or cross edge.
    {
      fprintf(fout,"FORWARD OR CROSS EDGE: (" F_IID ",%d)-(" F_IID ",%d)\n",
              get_iid_fragment(frags,iv0),is0,
              get_iid_fragment(frags,iv1),is1);
    }
    break;
  default:
    fprintf(stderr,"Unsupported DFS coloring\n");
  }
  }
}

static void as_dfs_recursive
(
 FILE *    fout,
 IntFragment_ID iv0,
 int       is0,
 IntRank   *nfound,
 IntRank   *timestamp,
 IntFragment_ID nfrag, 
 Tfragment  frags[],
 IntEdge_ID nedge, 
 Tedge      edges[], 
 IntRank    dfs_predecessor[],
 IntRank    fragment_rank[],
 IntRank    dfs_discovered[],
 IntRank    dfs_finished[],
 DFScolor   dfs_color[],
 int fragment_direction_is_forward_in_component,
 int sum_of_ahg_in_component,
 int sum_of_bhg_in_component,
 DFSstackobj * dfs_stack
)
{

  {
    IntFragment_ID iid = get_iid_fragment(frags,iv0);
    
    fprintf(fout,
            "TREE FRAG: " F_IID " %d %d %d %d\n",
            iid, is0, fragment_direction_is_forward_in_component,
            sum_of_ahg_in_component, sum_of_bhg_in_component );
  }
  
  // A recursive DFS for the AS_CGB data structures.

  found_a_fragment
    ( iv0, nfound, timestamp,
      fragment_rank, dfs_discovered, dfs_finished, dfs_color);

  {
    /* For each vertex adjacent to the current vertex, explore the
       overlap edges. */
    int isx;
    for(isx=1;isx >= 0; isx--) {
      const int isa = ( is0 ^ isx );
      // Implicitly follow the "fragment edge" first.

      /* Search all vertices "iv1" adjacent from "iv0" */
      const IntEdge_ID ie0 = get_segstart_vertex(frags,iv0,isa);
      const int ne0 = get_seglen_vertex(frags,iv0,isa);
#if 0
      fprintf(fout,"as_dfs: fragment " F_IID " %d " F_IID " %d\n",
              get_iid_fragment(frags,iv0),isa,ie0,ne0);
#endif
      {
        // Process the dovetail edges.
        int in0;
        for(in0=0; in0 < ne0; in0++) {
          IntEdge_ID ie1 = ie0 + in0;
          Tnes nes = get_nes_edge(edges,ie1);
          
          if(
             /* Only follow edges in the sub-graph that we are
                exploring. */
             
             // Ian: This is where you customise the search to whichever
             // phase of the Unitigger's versions of the fragment
             // overlap graph that you are looking for bubbles.
             (AS_CGB_INTERCHUNK_EDGE == nes) ||
             (AS_CGB_INTRACHUNK_EDGE == nes) ||
             (AS_CGB_TOUCHES_CONTAINED_EDGE == nes) ||
             (AS_CGB_BETWEEN_CONTAINED_EDGE == nes)
#ifdef PROCESS_CONTAINED_FRAGS_ALSO             
             || (AS_CGB_TO_CONTAINED_EDGE == nes)
#endif // PROCESS_CONTAINED_FRAGS_ALSO             
             
             ) {
            process_edge(
                         fout,
                         iv0, isa,
                         ie1,
                         nfound,
                         timestamp,
                         nfrag, 
                         frags,
                         nedge, 
                         edges,
                         dfs_predecessor,
                         fragment_rank,
                         dfs_discovered,
                         dfs_finished,
                         dfs_color,
                         fragment_direction_is_forward_in_component,
                         sum_of_ahg_in_component,
                         sum_of_bhg_in_component,
                         dfs_stack
                         );
          }
        }
      }
    }
  }

  finished_a_fragment
    ( iv0, nfound, timestamp,
      fragment_rank, dfs_discovered, dfs_finished, dfs_color);

}


void as_dfs_graph_traversal
(
 FILE *    fout,
 Tfragment frags[],
 Tedge     edges[], 
 IntRank   fragment_rank[]
 )
{ 
  /* A traversal of the graph. */
  IntRank nconnected=0,nfound=0,ncontained=0,ncircl=0,timestamp=0;
  IntFragment_ID vid;
  
  const IntFragment_ID nfrag = GetNumFragments(frags);
  const IntEdge_ID nedge = GetNumEdges(edges);

  DFSstackobj * dfs_stack = dfs_stack_create(nfrag);

  IntRank * dfs_predecessor = NULL; 
  IntRank * dfs_discovered = NULL;
  IntRank * dfs_finished = NULL; 
  DFScolor * dfs_color = NULL; 
  
  SAFE_MALLOC(dfs_predecessor, IntRank, nfrag);
  SAFE_MALLOC(dfs_discovered,  IntRank, nfrag);
  SAFE_MALLOC(dfs_finished,    IntRank, nfrag);
  SAFE_MALLOC(dfs_color,      DFScolor, nfrag);
  
  
  for(vid=0; vid<nfrag; vid++) {
    dfs_predecessor[vid] = AS_DFS_NOT_RANKED;
    fragment_rank[vid] = AS_DFS_NOT_RANKED;
    dfs_discovered[vid] = AS_DFS_NOT_RANKED;
    dfs_finished[vid] = AS_DFS_NOT_RANKED;
    dfs_color[vid] = AS_DFS_WHITE;
  }

  for(vid=0;vid<nfrag;vid++) {
    const int con = get_con_fragment(frags,vid);
    const int lab = get_lab_fragment(frags,vid);
    // We prefer to start each DFS connected component with a
    // non-contained fragment.
    if(
#ifndef FRAGMENT_SUBSET
       (con == FALSE ) &&
#else
       ((AS_CGB_INTERCHUNK_FRAG == lab) ||
        (AS_CGB_INTRACHUNK_FRAG == lab) ||
        (AS_CGB_THRU_FRAG == lab) ||
        (AS_CGB_HANGING_CHUNK_FRAG == lab) ||
        (AS_CGB_BRANCHMULTICONT_FRAG == lab) ) &&
#endif       
       (dfs_predecessor[vid] == AS_DFS_NOT_RANKED)){
      IntFragment_ID iid = get_iid_fragment(frags,vid);
      int sum_of_ahg_in_component = 0;
      int sum_of_bhg_in_component = 0;
      int fragment_direction_is_forward_in_component = TRUE;

      assert(fragment_rank[vid] == AS_DFS_NOT_RANKED);
      assert(dfs_discovered[vid] == AS_DFS_NOT_RANKED);
      assert(dfs_finished[vid] == AS_DFS_NOT_RANKED);
      assert(dfs_color[vid] == AS_DFS_WHITE);
      nconnected++;
      fprintf(fout, 
	      "DFS connected component starts at rank = " F_SIZE_T " and iid=" F_IID "\n",
	      nfound, iid);

        
      dfs_predecessor[vid] = AS_DFS_ROOT;
      as_dfs_recursive
	(fout,
	 vid, FALSE,
	 &nfound, &timestamp,
	 nfrag, frags, nedge, edges, 
	 dfs_predecessor,fragment_rank,
	 dfs_discovered, dfs_finished, dfs_color,
         fragment_direction_is_forward_in_component,
         sum_of_ahg_in_component,
         sum_of_bhg_in_component,
         dfs_stack
         );
      
    }
  }
  
  for(vid=0;vid<nfrag;vid++) {
    if(
       (dfs_predecessor[vid] == AS_DFS_NOT_RANKED) ||
       (fragment_rank[vid] == AS_DFS_NOT_RANKED) ||
       (dfs_discovered[vid] == AS_DFS_NOT_RANKED) ||
       (dfs_finished[vid] == AS_DFS_NOT_RANKED) ||
       (dfs_color[vid] != AS_DFS_BLACK) ) {
      const int lab = get_lab_fragment(frags,vid);
      if( (AS_CGB_UNPLACEDCONT_FRAG != lab) &&
          (AS_CGB_SINGLECONT_FRAG != lab) &&
          (AS_CGB_MULTICONT_FRAG  != lab) ) {
        // The AS_CGB_SINGLECONT_FRAG and AS_CGB_MULTICONT_FRAG
        // fragments do not have dovetail overlaps to them in the
        // visible graph.
        const Fragment_ID uid = get_uid_fragment(frags,vid);
        const IntFragment_ID iid = get_iid_fragment(frags,vid);
        const int con = get_con_fragment(frags,vid);
        const int lab = get_lab_fragment(frags,vid);
        fprintf(stderr,
                "DFS missed fragment (" F_UID "," F_IID "," F_IID ")"
                " %d %d \n",
                uid,iid,vid, con, lab);
      }
    }
    //ncircl++;
    //nfound++;
  }

#if 0
  assert( nfrag == nfound );
  assert( nfrag*2 == timestamp );
#endif
  
  fprintf(fout,"FRAGMENT OVERLAP GRAPH TRAVERSAL STATISTICS\n\n"
	  "%15" F_IIDP " : number of fragments in overlap graph\n"
	  "%15" F_SIZE_TP " : number of fragments found by search\n"
	  "%15" F_SIZE_TP " : number of weakly connected graph components\n"
	  "%15" F_SIZE_TP " : number of contained fragments\n"
	  "%15" F_SIZE_TP " : number of circular chunks\n\n\n",
	  nfrag,nfound,nconnected,ncontained,ncircl);

  /* 
     Some of the search_method()s assume that circular chunks do not occur.
     If this assert ever fails than look for fragments with
     get_lab_fragment() == AS_CGB_INTRACHUNK_FRAG, but have not been visited.
  */
  SAFE_FREE(dfs_predecessor);
  SAFE_FREE(dfs_discovered);
  SAFE_FREE(dfs_finished);
  SAFE_FREE(dfs_color);
  dfs_stack_destroy(dfs_stack);
  
}
