
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
= "$Id: AS_FGB_buddy.c,v 1.2 2004-09-23 20:25:01 mcschatz Exp $";
/* *******************************************************************
 *
 * Module: AS_CGB_buddy_dechorder.c
 * 
 * Description:
 *
 * Assumptions: 
 *
 * Author: Clark Mobarry
 *********************************************************************/

/*************************************************************************/
/* System include files */

/*************************************************************************/
/* Local include files */
#include "AS_CGB_all.h"
#define MAX(a,b) ((a)>(b)?(a):(b))

/****************************************************************************/
/* File Scope Globals */

static int TIMINGS = TRUE;

/*************************************************************************/
/* Conditional compilation */
#undef DEBUG

/*************************************************************************/
/* Global Defines */
#define EDGE_NOT_FOUND (~((IntEdge_ID)(0)))

/*************************************************************************/
/*************************************************************************/


static int is_dovetail_for_buddy_classification
(const Tnes nes0) 
{
  return
    (AS_CGB_INTERCHUNK_EDGE == nes0) ||
    (AS_CGB_INTRACHUNK_EDGE == nes0) ||
    (AS_CGB_TOUCHES_CONTAINED_EDGE == nes0) ||
    (AS_CGB_BETWEEN_CONTAINED_EDGE == nes0) ||
    (AS_CGB_TOUCHES_CRAPPY_DVT == nes0 ) ||
    (AS_CGB_BETWEEN_CRAPPY_DVT == nes0 )
    ;
}

#if 0
static int dovetail_degree_for_blizzard_pattern
(
 const Tfragment * frags,
 const Tedge     * edges,
 const IntFragment_ID vid,
 const int suf
 ) {
  const IntEdge_ID ie0 = get_segstart_vertex(frags,vid,suf);
  const int in = get_seglen_vertex(frags,vid,suf);
  IntEdge_ID ie1 = ie0;
  int count = 0;
  {
    // Count dovetail overlaps to starting from ie1:
    IntEdge_ID ie2;
    for(ie2 = ie1; ie2 < ie1 + in; ie2++) {
      const Tnes nes2 = get_nes_edge(edges,ie2);
      if( is_dovetail_for_blizzard_pattern(nes2) ) count ++;
    }
  }
  return count;
}
#endif

static IntEdge_ID find_thickest_dovetail_edge_for_buddy_classification
(
 const Tfragment * frags,
 const Tedge     * edges,
 const IntFragment_ID vid,
 const int suf
 ) {
  // Assumes that the fragment-end adjacency list is already sorted
  // from thickest to thinnest.
  
  const IntEdge_ID ie0 = get_segstart_vertex(frags,vid,suf);
  const int in0 = get_seglen_vertex(frags,vid,suf);
  IntEdge_ID ie1 = ie0;
  IntEdge_ID thickest_dovetail_edge = EDGE_NOT_FOUND;
  
#if 0
  if( in0 == 0) {
    fprintf(stdout,
	    "ie0=" F_IID " in0=%d iid=" F_IID " vid=" F_IID " suf=%d\n",
	    ie0, in0, get_iid_fragment(frags,vid), vid, suf);
  }
  // assert( in0 > 0 );
#endif
  // We only invoke this routine from a fragment-end that has a
  // dovetail overlap.

  {
    // Skip over non-dovetail overlaps to find edge ie1:

    for(ie1 = ie0; ie1 < ie0 + in0; ie1++) {
#if 1
      const Tnes nes1 = get_nes_edge(edges,ie1);
      if( is_dovetail_for_buddy_classification(nes1) ) {
	thickest_dovetail_edge = ie1;
	break;
      }
#else
      const int ahg1 = get_ahg_edge(edges,ie1);
      const int bhg1 = get_bhg_edge(edges,ie1);
      if( (ahg1 > 0) && (bhg1 > 0) ) {
	// Is dovetail for FGB buddy dechorder....
	thickest_dovetail_edge = ie1;
	break;
      }
#endif
    }
  } 
  return thickest_dovetail_edge;
}

static IntEdge_ID count_buddy_overlaps
( Tfragment * frags,
  Tedge     * edges
  )
{
  //const IntFragment_ID nfrag = GetNumFragments(frags);
  const IntEdge_ID nedge = GetNumEdges(edges);
  IntEdge_ID ie0;
  IntEdge_ID buddy_edges_count = 0;
#ifdef DEBUG
  IntEdge_ID buddy_invalid_edges_count = 0;
#endif  
  time_t tp1,tp2;

  if(TIMINGS) {
    tp1 = time(NULL); fprintf(stderr,"Begin buddy edge classification.\n");
    system_top();
  }

  for(ie0=0; ie0 < nedge; ie0++) {
    // For each candidate dovetail overlap edge.
    const IntFragment_ID b0_vid = get_avx_edge(edges,ie0);
    const IntFragment_ID c0_vid = get_bvx_edge(edges,ie0);
    const int b0_suf = get_asx_edge(edges,ie0);
    const int c0_suf = get_bsx_edge(edges,ie0);
    const Tnes nes0 = get_nes_edge(edges,ie0);
#if 0
    fprintf(stdout,
	    "ie0=" F_IID " b0_iid=" F_IID " b0_suf=%d c0_iid=" F_IID " c0_suf=%d nes=%d\n",
	    ie0,
	    get_iid_fragment(frags,b0_vid),
	    b0_suf,
	    get_iid_fragment(frags,c0_vid),
	    c0_suf,
	    nes0);
#endif
    if( is_dovetail_for_buddy_classification(nes0) ) {
      const IntEdge_ID ie1 =
        find_thickest_dovetail_edge_for_buddy_classification
	(frags,edges,b0_vid,b0_suf);
      // This is edge B0-A0.
      const IntEdge_ID ie2 =
        find_thickest_dovetail_edge_for_buddy_classification
	(frags,edges,c0_vid,c0_suf);
      // This is edge C0-D0.

      if( (ie1 != EDGE_NOT_FOUND) &&
	  (ie2 != EDGE_NOT_FOUND) ) {
	const IntFragment_ID a0_vid = get_bvx_edge(edges,ie1);
	const IntFragment_ID d0_vid = get_bvx_edge(edges,ie2);
	const int a0_suf = get_bsx_edge(edges,ie1);
	const int d0_suf = get_bsx_edge(edges,ie2);
	
	//  --- A1-A0 ---- B0-B1 ---
	//               /
	//              /          
	//             /
	//            /
	//  --- C1-C0 ---- D0-D1 ---
	
	if( (a0_vid == c0_vid) && (a0_suf == c0_suf) &&
	    (d0_vid == b0_vid) && (d0_suf == b0_suf)
	    ) {
	  // The B0-C0 overlap is a mutually thickest overlap.
	  
	  buddy_edges_count ++;
	  set_nes_edge(edges,ie0,AS_CGB_INTRACHUNK_EDGE);
	  fix_overlap_edge_mate(frags,edges,ie0);
	  
#ifdef DEBUG
	  if( get_inv_edge(edges,ie0) ) {
	    buddy_invalid_edges_count ++;
	    fprintf(stdout,"BUDDYINV nes=%d\n",nes0);
	    fprintf(stdout,
		    "a0_iid=" F_IID " a0_suf=%d "
		    "b0_iid=" F_IID " b0_suf=%d "
		    "c0_iid=" F_IID " c0_suf=%d "
		    "d0_iid=" F_IID " d0_suf=%d "
		    "\n",
		    get_iid_fragment(frags,a0_vid), a0_suf,
		    get_iid_fragment(frags,b0_vid), b0_suf,
		    get_iid_fragment(frags,c0_vid), c0_suf,
		    get_iid_fragment(frags,d0_vid), d0_suf
		    );
	    
	  } else {
	    fprintf(stdout,"BUDDYVAL nes=%d\n",nes0);
	  }
#endif            
	} else if( AS_CGB_INTRACHUNK_EDGE == nes0 ) {
	  set_nes_edge(edges,ie0,AS_CGB_INTERCHUNK_EDGE);
	}
      }
    }
  }

  if(TIMINGS) {
    tp1 = time(NULL); fprintf(stderr,"Finished buddy edge classification.\n");
    system_top();
  }

  {
    fprintf(stderr,"buddy_edges_count=" F_IID "\n",
            buddy_edges_count);
#ifdef DEBUG
    fprintf(stderr,"buddy_invalid_edges_count=" F_IID "\n",
            buddy_invalid_edges_count);
#endif            
  } 

  return buddy_edges_count/2;
}

//////////////////////////////////////////////////////////////////////

void buddy_dechorder
(
 Tfragment frags[],
 Tedge edges[]
 ) {
  count_buddy_overlaps( frags, edges);
}

#if 0
  // This routine depends on the edges are in an adjacency list
  // representation which is sorted by ahg. The operation count is
  // O(nedges) of stride one operations and O(nfrags) of indirect
  // operations?
  
  const IntFragment_ID nfrag=GetNumFragments(frags);
  const IntEdge_ID nedge=GetNumEdges(edges);
  IntEdge_ID count_dvt=0, count_con=0;

  {
    // Initialize the edges for marking the thickest dovetail overlaps
    // for each fragment-end.

    IntEdge_ID ie0=0;
    {
      const int ahg0 = get_ahg_edge(edges,ie0);
      const int bhg0 = get_bhg_edge(edges,ie0);
      if( ahg0 > 0 && bhg0 > 0 ) {
        set_nes_edge(edges,ie0,AS_CGB_INTERCHUNK_EDGE);
        // Set the dovetail edges to a uniform status.
      }
    }
    for(ie0=1; ie0 < nedge; ie0++) {
      const IntFragment_ID avx0 = get_avx_edge(edges,ie0);
      const int asx0 = get_asx_edge(edges,ie0);
      const int ahg0 = get_ahg_edge(edges,ie0);
      const int bhg0 = get_bhg_edge(edges,ie0);
      
      IntEdge_ID ie1 = ie0 - 1;
      const IntFragment_ID avx1 = get_avx_edge(edges,ie1);
      const int asx1 = get_asx_edge(edges,ie1);
      const int ahg1 = get_ahg_edge(edges,ie1);
      const Tnes nes1 = get_nes_edge(edges,ie1);
      
      assert( !((avx0 == avx1)&&(asx0 == asx1)) || (ahg0 >= ahg1) );
      // The O(nfrags) of indirection algorithm will use assumption that
      // the adjacency list is sorted for increasing ahg.
      
      if(is_a_dvt_simple(ahg0,bhg0)) {
        set_nes_edge(edges,ie0,AS_CGB_INTERCHUNK_EDGE);
        // Set the dovetail edges to a uniform status.
      }
      
    }
  }
    
  {
    // Mark the thickest dovetail overlap for each fragment end.

    IntFragment_ID iv0;
    int is0;
    for(iv0=0; iv0 < nfrag; iv0++) {
      for(is0=0; is0 < 1; is0++) {
        const IntEdge_ID ie0 = get_segstart_vertex(frags,iv0,is0);
        const IntEdge_ID in0 = get_seglen_vertex(frags,iv0,is0);
        IntEdge_ID ie1 = ie0;
        for(ie1=ie0; ie1 < ie0 + in0; ie1++) {
          const int ahg1 = get_ahg_edge(edges,ie1);
          const int bhg1 = get_bhg_edge(edges,ie1);
          assert( ahg1 >= 0 );
          if(is_a_dvt_simple(ahg1,bhg1)) {
            set_nes_edge(edges,ie1,AS_CGB_INTRACHUNK_EDGE);
            break;
            // We have the thickest dovetail overlap for this fragment
            // end.
          }
        } 
      }
    }
  }

  {
    // Only keep the marked overlaps that were the thickest from both
    // fragment ends as AS_CGB_INTRACHUNK. What remains are the buddy
    // dovetail overlaps.

    IntEdge_ID ie0;
    for(ie0=0; ie0 < nedge; ie0++) {
      const Tnes nes0 = get_nes_edge(edges,ie0);
      if( nes0 == AS_CGB_INTRACHUNK_EDGE ) {
        IntEdge_ID ie1 = find_overlap_edge_mate(frags,edges,ie0);
        const Tnes nes1 = get_nes_edge(edges,ie1);
        if( nes1 != AS_CGB_INTRACHUNK_EDGE ) {
          set_nes_edge(edges,ie0,AS_CGB_INTERCHUNK_EDGE);
        }
      }
    }
  }

  // Find the fragments at the ends of the chunks, and stuff.
  chunk_classification_of_fragments( frags, edges);


  {
    // Make the light-weight chunks with coordinates for each
    // fragment.

    const int pass = 0;
    // pass==0 is to form the light-chunks, that is chunks ignoring
    // all of the contained fragments.
    
    /* Initialize a flag for chunk following. */
    { int vid; for(vid=0;vid<GetNumFragments(frags);vid++) { 
      contained_timesinchunks[vid] = 0;
    }}
    
    { IntFragment_ID vid; for(vid=0;vid<GetNumFragments(frags);vid++) { 
      fragment_timesinchunks[vid] = 0;
      set_cid_fragment(frags,vid,CHUNK_NOT_VISITED);
      // Just a sentinal chunk id to designate that a fragment has not
      // been placed in any chunks yet.  The contained fragments will
      // be check to see if they can be placed in multiple chunks.
      set_o3p_fragment(frags,vid,0);
      set_o5p_fragment(frags,vid,0);
      // We need to check if a contained fragment is uniquely placed
      // in a light weight chunk by its coordinates as well.
    }}
    
    ResetVA_AChunkMesg(chunkfrags);
    ResetVA_AChunkMesg(thechunks);
    /* Initialize flags for chunk membership. */

    // Assign chunks to essential fragments at the ends of chunks.
    { IntFragment_ID vid; for(vid=0;vid<GetNumFragments(frags);vid++) 
      /* beginning loop through fragment */ {
      const Tlab ilabv = get_lab_fragment(frags,vid);
      // const IntFragment_ID iafr = get_iid_fragment(frags,vid);
      // const IntFragment_ID iafr = get_iid_fragment(frags,vid);
      if(
         (ilabv == AS_CGB_SOLO_FRAG) ||
         (ilabv == AS_CGB_HANGING_FRAG) ||
         (ilabv == AS_CGB_THRU_FRAG) ||
         ((ilabv == AS_CGB_INTERCHUNK_FRAG)
          &&(fragment_timesinchunks[vid] == 0)) ||
         ((ilabv == AS_CGB_HANGING_CHUNK_FRAG)
          &&(fragment_timesinchunks[vid] == 0))
         ) {
        fill_a_chunk_starting_at
          ( pass, vid,
            frags, edges,
            fragment_timesinchunks, contained_timesinchunks,
            chunkfrags, thechunks);
      }
    }
    /* End of a loop over all fragments. */
    }
  }


  // Now if a dovetail overlap (1) has both fragment-ends in the same
  // chunk, (2) is AS_CGB_INTERCHUNK, and (3) is consistent with the
  // chunk coordinates, then it is a chord overlap. Mark it and its
  // mate as AS_CGB_REMOVED_BY_TRANSITIVITY_DVT.

  // Pack the edge data structure.
}
#endif


