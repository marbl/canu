
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
= "$Id: AS_CGB_walk.c,v 1.1.1.1 2004-04-14 13:50:05 catmandew Exp $";
/*********************************************************************
 *
 * Module: AS_CGB_walk.c
 * Description: Walks the skeleton fragment overlap graph to determine
 * if an overlap is transitively inferable from the skeleton fragment 
 * overlap graph.
 * Assumptions:
 * TBD
 *
 * Author: Clark Mobarry
 *********************************************************************/

/*************************************************************************/
/* System include files */

/*************************************************************************/
/* Local include files */
#include "AS_CGB_all.h"

/* Conditional compilation */
#define JUST_OVERHANGS
#define USE_MARKED_EDGES
#undef PRINT_WALK
#undef DEBUG

#define O2CI(a)    ((a) >> 1) 
// From the chunk-end, return the chunk index.
#define O2S(a)     ((a) & 1)
// From the chunk-end, return the chunk suffix.
#define CIS2O(a,b) (((a) << 1) + ((b) & 1))
// From the chunk index and suffix, return the chunk-end.
#define O2O(a)     ((a) ^ 1)
// From the chunk-end, return the mate chunk-end.

int is_there_an_overlap_path
(
  const Tfragment *frags,
  const Tedge     *edges,
  const IntFragment_ID target_avx,
  const int        target_asx,
#ifdef MATCH_TARGET_EDGE
  const IntFragment_ID target_bvx,
  const int        target_bsx,
  const int        target_ahg,
  const int        target_bhg,
  const Tnes       target_nes,
  const int        target_is_dovetail,
  const int        target_is_from_contained,
  const int        target_is_to_contained,
  const int        target_is_dgn_contained,
#endif // MATCH_TARGET_EDGE
  /* recursion variables: */
  const int        search_depth, // Use zero at the top level.
  const IntFragment_ID current_avx,  // Use target_avx at the top level.
  const int        current_asx,  // Use target_asx at the top level.
  const int        current_ahg,  // Use zero at the top level.
  const int        current_bhg,  // Use zero at the top level.
  const int        last_edge_was_containment, // From target_nes at the top level.
  /* search path limiting: */
  const int        walk_depth,   // The maximum depth of the stack.
  const int        tolerance,    // For the overlaps
  IntFragment_ID visited_a[],
#ifdef MATCH_TARGET_EDGE
  IntFragment_ID visited_b[],
#endif // MATCH_TARGET_EDGE
  // Was this fragment visited from the target overlap edge before? That is
  // by the same (target_avx, target_asx, target_bvx, target_bsx).
  const int        work_limit_per_candidate_edge,
  // Maximum number of edges to explore per candidate edge.
  int *            work_tally_per_candidate_edge,
  // Current number of edges explored per candidate edge.
  /* diagnostic variables: */
#ifdef WALK_DEPTH_DIAGNOSTICS
  int64      search_depth_histogram[], // How many times this level is visited.
  int64      search_path_histogram[],  // How many paths were this length.
#endif // WALK_DEPTH_DIAGNOSTICS
  int64      *ntrans_test_fail // How many paths failed the last test.
  )
{
  // Given an undirected (directed) graph G = (V,E) with n verticies
  // and an array visited(n) initially set to FALSE, this algorithm
  // visits all vertices reachable from v.  G and visited are global.
  // We modifiy this to be a depth first search for each target
  // overlap edge.

  // The arrays visited_a[] and visited_b[] must be initialized out of
  // the range [0,2*nfrag-1] before the rooted depth first search
  // begins.

  // The return value is an non-negative integer. If the query overlap
  // is(is not) transitively inferable from a path in the skeleton
  // graph, then the return value is positive(zero).  A return value
  // of one is the normal status of an inferable edge. A return value
  // of two is stronger in that this edge is removable.


  const int tolerance_ahg = tolerance;
  const int tolerance_bhg = tolerance;

  const IntEdge_ID ir0 = get_segstart_vertex(frags,current_avx,current_asx);
  const int nnode  = get_seglen_vertex(frags,current_avx,current_asx);

#ifndef JUST_OVERHANGS
  const int target_aln = get_length_fragment(frags,target_avx);
  const int target_bln = get_length_fragment(frags,target_bvx);
#endif // JUST_OVERHANGS
  int in1;


  assert(tolerance >= 0);
  if(!(target_is_dovetail || target_is_from_contained ||
       target_is_to_contained || target_is_dgn_contained)) return FALSE;
  // This should be outside of the recursion stack.

  if((*work_tally_per_candidate_edge >=
      work_limit_per_candidate_edge)) return FALSE;
  
  if(search_depth >= walk_depth) {
    //search_path_histogram[search_depth]++;
    return FALSE; // For backward compatability...
  }
  
  // Set up a locally rooted depth first search.  We are depending on
  // visited_a[] and visited_b[] to be initialized out of the range
  // [0,2*nfrag-1] before the rooted depth first search begins.
  
  visited_a[CIS2O(current_avx,current_asx)]=CIS2O(target_avx,target_asx);
  visited_b[CIS2O(current_avx,current_asx)]=CIS2O(target_bvx,target_bsx);

  // for each edge ir1 adjacent to v 
  for(in1=0;in1<nnode;in1++) { 
    // Assume dovetail and to-contained edges sorted for increasing ahg.
    // What about FROM_CONTAINED edges??
    const IntEdge_ID  ir1 = ir0+in1;

    /* Dovetail Overlap types (symmetric edges) */
    const int ir1_is_dovetail = is_a_dvt_edge(edges,ir1);

    /* Containment Overlap types (asymmetric edges) */
    const int ir1_is_from_contained = is_a_frc_edge(edges,ir1);
    // A containment edge from a locally contained fragment to the
    // containing fragment.
    const int ir1_is_to_contained = is_a_toc_edge(edges,ir1);
    // A containment edge to a locally contained fragment from the
    // containing fragment.

    
    if((*work_tally_per_candidate_edge
        >= work_limit_per_candidate_edge)) return FALSE;
    (*work_tally_per_candidate_edge) ++;
    // This appears to be the correct place, but hey.
    
    if( // Restrictions on overlap type
       ((ir1_is_dovetail && target_is_dovetail) ||
        (ir1_is_to_contained && target_is_to_contained) ||
        (ir1_is_from_contained && target_is_from_contained))
      ) // Restrictions on overlap type
    {
    
      const int new_last_edge_was_containment = 
	last_edge_was_containment &&  ir1_is_from_contained;
        // The path is restricted to be (FRC)*(DVT)* .
        // This logic depends on "last_edge_was_containment"
        // initialized to TRUE before a inferring path is searched.

      /* Using the Overlap Record notation: */
      const IntFragment_ID ir1avx = get_avx_edge(edges,ir1);
      const int        ir1asx = get_asx_edge(edges,ir1);
      const int        ir1ahg = get_ahg_edge(edges,ir1);
      const IntFragment_ID ir1bvx = get_bvx_edge(edges,ir1);
      const int        ir1bsx = get_bsx_edge(edges,ir1);
      const int        ir1bhg = get_bhg_edge(edges,ir1);
      
      const IntFragment_ID new_avx = ir1bvx;
      const int            new_asx = !ir1bsx;
      const int new_ahg = current_ahg + ir1ahg;
      const int new_bhg = current_bhg + ir1bhg;
      
#ifndef JUST_OVERHANGS
      const int new_aln = target_aln;
      const int new_bln = get_length_fragment(frags,ir1bvx);
      // Note that since we are constructing the path from the
      // A-side that new_aln and target_aln are equal.
      
      const int target_alp = target_aln - target_ahg;
      const int target_blp = target_bln - target_bhg;
      const int new_alp = new_aln - new_ahg;
      const int new_blp = new_bln - new_bhg;
      // Note that we are implicitly using Granger^s extended
      // fragments when new_ahg < 0 or new_bhg < 0.
#endif // JUST_OVERHANGS
      
#ifdef WALK_DEPTH_DIAGNOSTICS
      search_depth_histogram[search_depth]++;
#endif      
      assert(current_avx == ir1avx);
      assert(current_asx == ir1asx);
      
      // Each middle fragment must include the target overlap sequence. 
      if(
        /* If a middle fragment does not include the target overlap,
           then exit the path. */
#ifndef JUST_OVERHANGS
        (new_alp < target_alp - tolerance_ahg) ||
        (new_blp < target_blp - tolerance_bhg)
#else  // JUST_OVERHANGS
        !(((target_ahg-new_ahg)*(target_ahg >= 0 ? 1 : -1)
           + tolerance_ahg > 0)
          &&
          ((target_bhg-new_bhg)*(target_bhg >= 0 ? 1 : -1)
           + tolerance_bhg > 0)
           )
#endif // JUST_OVERHANGS
        ) continue;
      
      if(search_depth > 0) {
        if( (ir1bvx == target_bvx) &&
            (ir1bsx == target_bsx) ) {
          if( (ABS(new_ahg - target_ahg) <= tolerance_ahg) &&
              (ABS(new_bhg - target_bhg) <= tolerance_bhg) ) {
#ifdef WALK_DEPTH_DIAGNOSTICS
	    search_path_histogram[search_depth]++;
	    // Record the length of a succesful path.
#endif // WALK_DEPTH_DIAGNOSTICS

#ifdef PRINT_WALK
         {
           FILE *fwalk = stderr;
           fprintf(fwalk,
                   "Walk1: " F_IID ":%d " F_IID ":%d %d %d %d : %d " F_IID ":%d " F_IID ":%d %d %d %d\n",
                   get_iid_fragment(frags,target_avx),
                   target_asx,
                   get_iid_fragment(frags,target_bvx),
                   target_bsx,
                   target_ahg,
                   target_bhg,
                   target_nes,
                   search_depth, // Use zero at the top level.
                   get_iid_fragment(frags,current_avx),  // Use target_avx at the top level.
                   current_asx,  // Use target_asx at the top level.
                   get_iid_fragment(frags,ir1bvx),
                   ir1bsx,
                   current_ahg,  // Use zero at the top level.
                   current_bhg,  // Use zero at the top level.
                   new_last_edge_was_containment); // From target_nes at the top level.
         }   
#endif // PRINT_WALK
            return TRUE; // A chord overlap!!
          } else {
            (*ntrans_test_fail)++;
          }
        }
      }

      // If this vertex has been processed from the (target_avx,
      // target_asx) before, then terminate the search.  This makes
      // the graph search into a local DFS where the root of the DFS
      // is (target_avx,target_asx).
      
      // Depth first search only: if (visited(new_avx,new_asx) ==
      // FALSE) for this local search, then do not follow a redundant
      // sub-path.

      if(
        (visited_a[CIS2O(new_avx,new_asx)]
         ==CIS2O(target_avx,target_asx))
        &&
        (visited_b[CIS2O(new_avx,new_asx)]
         ==CIS2O(target_bvx,target_bsx)) 
        ) continue;
      
      {
       int iflag = is_there_an_overlap_path
        ( frags, edges,
          target_avx, target_asx,
#ifdef MATCH_TARGET_EDGE
          target_bvx, target_bsx, 
          target_ahg, target_bhg,
	  target_nes,
          target_is_dovetail,
          target_is_from_contained,
          target_is_to_contained,
          target_is_dgn_contained,
#endif // MATCH_TARGET_EDGE
          /* recursion variables: */
          (search_depth+1), // Use zero at the top level.
          new_avx, new_asx, 
          new_ahg, new_bhg,
          new_last_edge_was_containment,
          /* search path limiting: */
          walk_depth,
          tolerance,
          visited_a, // Was this fragment seen from the target overlap before?
#ifdef MATCH_TARGET_EDGE
          visited_b, 
#endif // MATCH_TARGET_EDGE
          work_limit_per_candidate_edge,
          // Maximum number of edges to explore per candidate edge.
          work_tally_per_candidate_edge,
          // Current number of edges explored per candidate edge.
          /* diagnostics: */
#ifdef WALK_DEPTH_DIAGNOSTICS
	  search_depth_histogram,
	  search_path_histogram,
#endif // WALK_DEPTH_DIAGNOSTICS
          ntrans_test_fail
	  );
       if(iflag) {
#ifdef PRINT_WALK
         {
           FILE *fwalk = stderr;
           fprintf(fwalk,
                   "Walk2: " F_IID ":%d " F_IID ":%d %d %d %d : %d " F_IID ":%d " F_IID ":%d %d %d %d\n",
                   get_iid_fragment(frags,target_avx),
                   target_asx,
                   get_iid_fragment(frags,target_bvx),
                   target_bsx,
                   target_ahg,
                   target_bhg,
                   target_nes,          
                   search_depth, // Use zero at the top level.
                   get_iid_fragment(frags,current_avx),  // Use target_avx at the top level.
                   current_asx,  // Use target_asx at the top level.
                   get_iid_fragment(frags,ir1bvx), 
                   ir1bsx,  
                   current_ahg,  // Use zero at the top level.
                   current_bhg,  // Use zero at the top level.
                   new_last_edge_was_containment); // From target_nes at the top level.
         }   
#endif // PRINT_WALK
         return iflag; // propagate a eureka
       }
      }
    } // Restrictions on overlap types.
  } // in1

  return FALSE;
}
