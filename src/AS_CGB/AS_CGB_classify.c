
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
= "$Id: AS_CGB_classify.c,v 1.3 2005-03-22 19:02:06 jason_miller Exp $";
/* *******************************************************************
 *
 * Module: AS_CGB_classify.c
 * 
 * Description: These routines classify the fragments and overlaps in
 * a local manner as to their status in chunks.
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
//#include "AS_CGB_blizzard.h"

/*************************************************************************/
/* Conditional compilation */
#undef DEBUGVIEW
#undef USE_BUDDYCHUNK_EGDES

/*************************************************************************/
/* Global Defines */

/*************************************************************************/
/* File Scope Global Variables */
static int TIMINGS = TRUE;

/*************************************************************************/
#if 0
static void remap_edge_labels
(
 Tfragment frags[],
 Tedge edges[])
{
  /* Remap the edge labels to a near FGB status.  The exception are
    edges marked by spurs, crappy fragments, branchpoints, ....  */
  
  const IntEdge_ID nedge = GetNumEdges(edges);

  { IntEdge_ID ie; for(ie=0; ie<nedge; ie++) {
#if 0
    const IntFragment_ID iavx = get_avx_edge(edges,ie);
    const IntFragment_ID ibvx = get_bvx_edge(edges,ie);
    const IntFragment_ID aid = get_iid_fragment(frags,iavx);
    const IntFragment_ID bid = get_iid_fragment(frags,ibvx);
    const int iahg = get_ahg_edge(edges,ie);
    const int ibhg = get_bhg_edge(edges,ie);
#endif    
    const Tnes nes = get_nes_edge(edges,ie);

    switch(nes) {

      /* The containment overlaps... 'C' */
    case AS_CGB_TO_CONTAINED_EDGE:
      // do nothing ...
      break;
      
    case AS_CGB_FROM_CONTAINED_EDGE:
      // do nothing ...
      break;

      /* The dovetail overlaps... 'D' */
    case AS_CGB_DOVETAIL_EDGE:
    case AS_CGB_THICKEST_EDGE:
      break;

    case AS_CGB_BETWEEN_CONTAINED_EDGE:
      set_nes_edge(edges,ie,AS_CGB_DOVETAIL_EDGE); 
      fix_overlap_edge_mate(frags,edges,ie);
      break;
      
    case AS_CGB_INTERCHUNK_EDGE:
    case AS_CGB_BUDDYCHUNK_EDGE:
    case AS_CGB_INTRACHUNK_EDGE:
    case AS_CGB_TOUCHES_CONTAINED_EDGE:
      set_nes_edge(edges,ie,AS_CGB_INTERCHUNK_EDGE); 
      fix_overlap_edge_mate(frags,edges,ie);
      break;

    case AS_CGB_TOUCHES_CRAPPY_DVT:
    case AS_CGB_TOUCHES_CRAPPY_TOC:
    case AS_CGB_BETWEEN_CRAPPY_DVT:
    case AS_CGB_BETWEEN_CRAPPY_TOC:
    case AS_CGB_MARKED_BY_BRANCH_DVT:
    case AS_CGB_MARKED_BY_BREAKER:
    case AS_CGB_MARKED_BY_DELETED_DVT:
    case AS_CGB_MARKED_BY_DELETED_FRC:
    case AS_CGB_MARKED_BY_DELETED_TOC:
      // Leave as is ....
      break;
    default:
      fprintf(stderr,"Unsupported edge type %d\n",nes);
      assert(FALSE);
    }
  }}
}
#endif

static void maskout_overlaps_touching_crappy_fragments
( Tfragment * frags,
  Tedge     *edges
)
{
  /*
   * Label each dovetail edge to and from every crappy fragment.
   */
  const IntEdge_ID nedge = GetNumEdges(edges);
  { IntEdge_ID ie; for(ie=0; ie<nedge; ie++) {
    const IntFragment_ID iavx = get_avx_edge(edges,ie);
    const IntFragment_ID ibvx = get_bvx_edge(edges,ie); 
    const Tnes ines = get_nes_edge(edges,ie);
    Tnes jnes = ines;
    const int ikeep
      =  (AS_CGB_HANGING_CRAPPY_FRAG != get_lab_fragment(frags,iavx))
      && (AS_CGB_HANGING_CRAPPY_FRAG != get_lab_fragment(frags,ibvx));
    const int iremove
      =  (AS_CGB_HANGING_CRAPPY_FRAG == get_lab_fragment(frags,iavx))
      && (AS_CGB_HANGING_CRAPPY_FRAG == get_lab_fragment(frags,ibvx));
    
    switch(ines) {
    case AS_CGB_INTERCHUNK_EDGE:
    case AS_CGB_TOUCHES_CRAPPY_DVT:
    case AS_CGB_BETWEEN_CRAPPY_DVT:
      {
	jnes = ines;
	
	if(!ikeep && !iremove) { jnes = AS_CGB_TOUCHES_CRAPPY_DVT;}
	// This is a dovetail overlap between a globally crappy fragment
	// and a non-crappy fragment.
	
	if(iremove) { jnes = AS_CGB_BETWEEN_CRAPPY_DVT;}
	// This is a dovetail overlap between two globally crappy
	// fragments.
	
	if( ines != jnes ) {
	  set_nes_edge(edges,ie,jnes);
	  fix_overlap_edge_mate(frags, edges, ie);
	}
      }
      break;
    case AS_CGB_CONTAINED_EDGE:
    case AS_CGB_TOUCHES_CRAPPY_CON:
    case AS_CGB_BETWEEN_CRAPPY_CON:
      {
	jnes = ines;
	
	if(!ikeep && !iremove) { jnes = AS_CGB_TOUCHES_CRAPPY_CON;}
	// This is a containment overlap between a globally crappy fragment
	// and a non-crappy fragment.
	
	if(iremove) { jnes = AS_CGB_BETWEEN_CRAPPY_CON;}
	// This is a containment overlap between two globally crappy
	// fragments.
	
	if( ines != jnes ) {
	  set_nes_edge(edges,ie,jnes);
	  fix_overlap_edge_mate(frags, edges, ie);
	}
      }
      break;
    case AS_CGB_DOVETAIL_EDGE:
    case AS_CGB_THICKEST_EDGE:
      // Do nothing ....
      break;	
    default:
      fprintf(stderr,"Unexpected overlap label ines=%d\n",ines);
      assert(FALSE);
      break;
    }
  }}
}


// Essential overlaps do not involve spur or contained fragments.

// All dovetail overlaps are a superset of
// non-chord  dovetail edges are a superset of
// essential (non-chord not involving spur or contained fragments)

// dovetail edges are a superset of
// thickest edges are a superset of
// interchunk dovetail edges is a superset of
// thickchunk dovetail edges is a superset of
// buddychunk dovetail edges is a superset of
// intrachunk dovetail edges.

// all - non_chord = chord overlaps.

// non_chord - essential = touches and between contained, touches and
// between spurs, ...

// essential - interchunk = 

static int is_dovetail_for_chunk_classification
( int ines
  )
{
  int iret = FALSE;
  switch(ines) {
  case AS_CGB_INTERCHUNK_EDGE:
  case AS_CGB_BUDDYCHUNK_EDGE:
  case AS_CGB_INTRACHUNK_EDGE:
    iret = TRUE;
    break;

  case AS_CGB_CONTAINED_EDGE:
    /* The containment overlaps... 'C' */
    break;
    
  case AS_CGB_DOVETAIL_EDGE:
  case AS_CGB_THICKEST_EDGE:
  case AS_CGB_TOUCHES_CRAPPY_DVT:
  case AS_CGB_BETWEEN_CRAPPY_DVT:
  case AS_CGB_TOUCHES_CRAPPY_CON:
  case AS_CGB_BETWEEN_CRAPPY_CON:
    // ignored for the chunk intruder indegree.
    break;
    
  case AS_CGB_MARKED_BY_DELETED_DVT:
  case AS_CGB_MARKED_BY_DELETED_CON:
    // Leave as is ....
    break;
    
  default:
    fprintf(stderr,"Unsupported edge type %d\n",ines);
    assert(FALSE);
    break;
  }
  return iret;
}

#ifdef USE_BUDDYCHUNK_EGDES
static void buddychunk_classification_of_overlaps
(
  Tfragment frags[],
  Tedge edges[]
)
{
  //const IntFragment_ID nfrag = GetNumFragments(frags);
  const IntEdge_ID nedge = GetNumEdges(edges);
  IntEdge_ID ie0;
#ifdef DEBUG
  IntEdge_ID buddy_edges_count = 0;
  IntEdge_ID buddy_invalid_edges_count = 0;
#endif  
  time_t tp1 = 0,tp2;

  if(TIMINGS) {
    tp1 = time(NULL); fprintf(stderr,"Begin buddychunk edge classification.\n");
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
    if( is_dovetail_for_chunk_classification(nes0) ) {
      const IntEdge_ID ie1 =
        find_thickest_dovetail_edge_for_buddy_classification
	(frags,edges,b0_vid,b0_suf);
      // This is edge B0-A0.
      const IntEdge_ID ie2 =
        find_thickest_dovetail_edge_for_buddy_classification
	(frags,edges,c0_vid,c0_suf);
      // This is edge C0-D0.

      if( (ie1 != AS_CGB_EDGE_NOT_FOUND) &&
	  (ie2 != AS_CGB_EDGE_NOT_FOUND) ) {
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
	  set_nes_edge(edges,ie0,AS_CGB_BUDDYCHUNK_EDGE);
	  set_nes_edge(edges,ie1,AS_CGB_BUDDYCHUNK_EDGE);
	  //fix_overlap_edge_mate(frags,edges,ie0);
	  
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
	} else if( AS_CGB_BUDDYCHUNK_EDGE == nes0 ) {
          assert(FALSE);
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

  //return buddy_edges_count/2;
}
#endif // USE_BUDDYCHUNK_EGDES

static void find_in_degree_of_essential_dovetails
(
  /* input only */
  Tfragment *frags,
  Tedge     *edges,
  const int walk_depth,
  /* output only */
  int        fragment_npx[],
  // The dovetail edge in-degree for the prefix of the fragment read.
  int        fragment_nsx[]
  // The dovetail edge in-degree for the suffix of the fragment read.
  )
{
  // We want to count the number of in-coming essential dovetail
  // overlaps to each fragment-end in the transitively reduced graph
  // for the purpose of determining INTRACHUNK edges.
  
  const IntFragment_ID nfrag = GetNumFragments(frags);
  const IntEdge_ID nedge = GetNumEdges(edges);

  {
    IntFragment_ID vid;
    // #pragma omp parallel for
    for(vid=0; vid<nfrag; vid++){
      fragment_npx[vid] = 0;
      // The number of incoming essential dovetail edges that overlap
      // the prefix of the read.
      fragment_nsx[vid] = 0;
      // The number of incoming essential dovetail edges that overlap
      // the suffix of the read.
    }
  }
  
  /* Loop over edges add the number of each type of overlaps
   with each vertex. */
  /* OMP ISSUE:  This is a scatter-add. */
  assert(((!0) == 1) && ((!1) == 0));
  {
    IntEdge_ID ie; 
    for(ie=0; ie<nedge; ie++) {
      //const IntFragment_ID iavx = get_avx_edge(edges,ie); 
      const IntFragment_ID ibvx = get_bvx_edge(edges,ie); 
      const int iasx = get_asx_edge(edges,ie);
      const int ibsx = get_bsx_edge(edges,ie);
      const Tnes ines = get_nes_edge(edges,ie);
      const int edge_flag
        = is_dovetail_for_chunk_classification(ines);
      assert( (iasx == 0) || (iasx == 1));
      assert( (ibsx == 0) || (ibsx == 1));

      if(edge_flag) {
        /* Count the number of essential edges that overlap the
           prefix and the suffix respectively of each fragment. */
        // For the in-coming degree use the B-read.
        fragment_npx[ibvx] += 1-ibsx;
        fragment_nsx[ibvx] +=   ibsx;
      }
    }
  }
  
#if 0
  {
    IntFragment_ID vid;
    FILE *fp = fopen("Q1.txt","w");
    // #pragma omp parallel for
    for(vid=0; vid<nfrag; vid++){
      fprintf(fp,"% 10" F_IIDP " % 10" F_IIDP " | % 5d % 5d\n",
      vid,
      get_iid_fragment(frags,vid),
      fragment_npx[vid],
      // The number of dovetail edges that overlap the prefix of the
      // read.
      fragment_nsx[vid]);
      // The number of dovetail edges that overlap the suffix of the
      // read.
    }
    fclose(fp);
  }
#endif
#if 0
  {
    view_fgb_chkpnt("Q2", frags, edges);
  }
#endif  
}

static void intrachunk_classification_of_overlaps
(
  Tfragment frags[],
  Tedge edges[],
  const int walk_depth
)
{
  const IntFragment_ID nfrag = GetNumFragments(frags);
  const IntEdge_ID nedge = GetNumEdges(edges);
  int *fragment_npx = NULL; 
  int *fragment_nsx = NULL;
  SAFE_MALLOC(fragment_npx, int, nfrag);
  SAFE_MALLOC(fragment_nsx, int, nfrag);

#ifdef DEBUGVIEW
  view_fgb_chkpnt("Q0", frags, edges);
#endif    

  find_in_degree_of_essential_dovetails
    ( frags, edges,
      walk_depth,
      fragment_npx, fragment_nsx);
 
  /* #pragma omp for */  /* There is a scatter-add in this loop. */
  {
    IntEdge_ID ie;
    for(ie=0; ie<nedge; ie++) {
      const Tnes ines = get_nes_edge(edges,ie);
      const IntFragment_ID iavx = get_avx_edge(edges,ie);
      const IntFragment_ID ibvx = get_bvx_edge(edges,ie); 
      const int iasx = get_asx_edge(edges,ie);
      const int ibsx = get_bsx_edge(edges,ie); 
      assert( (ibsx == 0) || (ibsx == 1));
      assert( (iasx == 0) || (iasx == 1));

      // If this buddychunk edge is the only essential dovetail
      // overlap between two fragments, then it is an intra-chunk
      // edge.

      switch(ines) {
      case AS_CGB_INTERCHUNK_EDGE:
      case AS_CGB_INTRACHUNK_EDGE:
#ifdef USE_BUDDYCHUNK_EGDES
        break;
#endif // USE_BUDDYCHUNK_EGDES
      case AS_CGB_BUDDYCHUNK_EDGE:
        {
          /* When the range of asx and bsx is known to be {0,1}, then
             we can simplify this calculation. */
          /* Recall the number of essential edges that share the same
             terminal of the A-fragment as this edge. */
          const int naov = (iasx ? fragment_nsx[iavx] : fragment_npx[iavx]);
          /* Recall the number of essential edges that share the same
             terminal of the B-fragment as this edge. */
          const int nbov = (ibsx ? fragment_nsx[ibvx] : fragment_npx[ibvx]);
          int iconvert = (((naov == 1)&&(nbov == 1)));
          
          // REMEMBER THIS DETAIL!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          // Lastly, we protect dovetail overlaps to globally contained
          // fragments from being labeled as AS_CGB_INTRACHUNK_EDGE.
          const int iacn = get_con_fragment(frags,iavx);
          const int ibcn = get_con_fragment(frags,ibvx); 
          // assert((! iacn) && (! ibcn));
          iconvert = iconvert && ((! iacn) && (! ibcn));
          // REMEMBER THIS DETAIL!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          
          set_nes_edge
            (edges,ie, (iconvert ? AS_CGB_INTRACHUNK_EDGE : ines));
	  fix_overlap_edge_mate(frags,edges,ie);
        }
        break;
        
      case AS_CGB_DOVETAIL_EDGE:
      case AS_CGB_THICKEST_EDGE:
      case AS_CGB_TOUCHES_CONTAINED_EDGE:
      case AS_CGB_BETWEEN_CONTAINED_EDGE:
        break;

        /* The containment overlaps... 'C' */
      case AS_CGB_CONTAINED_EDGE:
        // do nothing ...
        break;

      case AS_CGB_TOUCHES_CRAPPY_DVT:
      case AS_CGB_BETWEEN_CRAPPY_DVT: 
      case AS_CGB_MARKED_BY_BRANCH_DVT:
      case AS_CGB_MARKED_BY_BREAKER:
      case AS_CGB_MARKED_BY_DELETED_DVT:

      case AS_CGB_TOUCHES_CRAPPY_CON:
      case AS_CGB_BETWEEN_CRAPPY_CON:
      case AS_CGB_MARKED_BY_DELETED_CON:
        // Leave as is ....
        break;
      default:
        fprintf(stderr,"Unsupported edge type %d\n",ines);
        assert(FALSE);
        break;
      } // End switch statement.
    }
  }
  SAFE_FREE(fragment_npx);
  SAFE_FREE(fragment_nsx);
}

static void chunk_classification_of_fragments
( Tfragment frags[], Tedge edges[]
)
{
  // We assume that the if the dovetail out-degree is zero, then the
  // dovetail in-degree must be zero.
  
  // AS_CGB_SOLO_FRAG, AS_HANGING_CRAPPY_FRAG (spur), contained, and deleted
  // fragments are invariant.

  // The fragments with the labels: AS_CGB_THRU_FRAG,
  // AS_CGB_INTERCHUNK_FRAG, AS_CGB_INTRACHUNK_FRAG,
  // AS_CGB_HANGING_FRAG, and AS_CGB_HANGING_CHUNK_FRAG are
  // re-labelled.

  // They receive labels from the set: AS_CGB_THRU_FRAG,
  // AS_CGB_INTERCHUNK_FRAG, AS_CGB_INTRACHUNK_FRAG,
  // AS_CGB_HANGING_FRAG, AS_CGB_HANGING_CHUNK_FRAG, and the death
  // sentence AS_CGB_HANGING_CRAPPY_FRAG.

  const IntFragment_ID nfrag = GetNumFragments(frags);
  const IntEdge_ID nedge = GetNumEdges(edges);
  char *intrachunk_npx = NULL; // boolean
  char *intrachunk_nsx = NULL;
  char *interchunk_npx = NULL; // boolean
  char *interchunk_nsx = NULL;

  SAFE_MALLOC(intrachunk_npx, char, nfrag); // boolean
  SAFE_MALLOC(intrachunk_nsx, char, nfrag);
  SAFE_MALLOC(interchunk_npx, char, nfrag); // boolean
  SAFE_MALLOC(interchunk_nsx, char, nfrag);

  /* Label each singleton chunk with a label consistent with its
     fragment graph status.  */

  /* At this point we expect the following information at the reads:
     fragment_nsx   fragment_npx
     0            0           this should not happen if there is
     at least one essential overlap.
     0            *           termination of a contig
     *            0           termination of a contig
     1            >1          fan-out 
     >1           1           fan-in
     >1           >1          "junction box"
  */

  {
    IntFragment_ID vid;
    // #pragma omp parallel for
    for(vid=0; vid<nfrag; vid++){
      intrachunk_nsx[vid] = 0; /* Will be the number of intrachunk edges 
				 that overlap the suffix of the read. */
      intrachunk_npx[vid] = 0; /* Will be the number of intrachunk edges
				 that overlap the prefix of the read. */
      interchunk_nsx[vid] = 0; /* Will be the number of interchunk edges 
				 that overlap the suffix of the read. */
      interchunk_npx[vid] = 0; /* Will be the number of interchunk edges
				 that overlap the prefix of the read. */
    }
  }
  
  {
    IntEdge_ID ie;
    for(ie=0; ie<nedge; ie++) {
      const IntFragment_ID iavx = get_avx_edge(edges,ie); 
      //const IntFragment_ID ibvx = get_bvx_edge(edges,ie); 
      const int iasx = get_asx_edge(edges,ie);
      const int ibsx = get_bsx_edge(edges,ie); 
      const Tnes ines = get_nes_edge(edges,ie);
      assert( (iasx == 0) || (iasx == 1));
      assert( (ibsx == 0) || (ibsx == 1));

      switch(ines) {
      case AS_CGB_INTRACHUNK_EDGE:
        {
          /* Count the number of intrachunk edges that overlap the
             prefix and the suffix respectively. */
          intrachunk_npx[iavx] += 1-iasx;
          intrachunk_nsx[iavx] +=   iasx;
          // Compute the intrachunk out-degree.
        }
        break;
      case AS_CGB_INTERCHUNK_EDGE:
      case AS_CGB_BUDDYCHUNK_EDGE:
        {
          interchunk_npx[iavx] += 1-iasx;
          interchunk_nsx[iavx] +=   iasx;
          // Compute the non-intrachunk dovetail out-degree.
        }
        break;
      case AS_CGB_TOUCHES_CONTAINED_EDGE:
        assert(FALSE);
        break;

      case AS_CGB_DOVETAIL_EDGE:
      case AS_CGB_THICKEST_EDGE:
      case AS_CGB_BETWEEN_CONTAINED_EDGE:
      case AS_CGB_TOUCHES_CRAPPY_DVT:
      case AS_CGB_BETWEEN_CRAPPY_DVT:
      case AS_CGB_MARKED_BY_BRANCH_DVT:

      case AS_CGB_CONTAINED_EDGE:
      case AS_CGB_TOUCHES_CRAPPY_CON:
      case AS_CGB_BETWEEN_CRAPPY_CON:
        // Do nothing.
        break;
      default:
        fprintf(stderr,"Unsupported edge type %d\n",ines);
        assert(FALSE);
        break;
        
      }
    }
  }

  {
    // Now examine the intrachunk and interchunk degree of each
    // fragment-end to find the chunk status of a fragment.
    
    IntFragment_ID vid;
    // #pragma omp parallel for
    for(vid=0; vid<nfrag; vid++) {
      const int deleted = get_del_fragment(frags,vid);
      const int contained = get_con_fragment(frags,vid);
      const int inpx = intrachunk_npx[vid]; 
      /* The number of intrachunk edges that overlap the prefix of the
         read. */
      const int insx = intrachunk_nsx[vid]; 
      /* The number of intrachunk edges that overlap the suffix of the
	 read. */
      const int jnpx = interchunk_npx[vid]; 
      /* The number of interchunk edges that overlap the prefix of the
         read. */
      const int jnsx = interchunk_nsx[vid]; 
      /* The number of interchunk edges that overlap the suffix of the
	 read. */

      const Tlab ilab = get_lab_fragment(frags,vid);

      // const IntFragment_ID iid = get_iid_fragment(frags,vid);
      
      /* There can not be more than one intrachunk edge into each port 
	   of the fragment. */
      assert( (inpx == 0) || (inpx == 1) ); 
      assert( (insx == 0) || (insx == 1) );
      /* if( inpx != 0 || insx != 0 ) then the fragment must be essential. */
      /* if( !(inpx == 0 && insx == 0 ) ) then the fragment must be essential. */

      /* Reset the fragment label: */
      set_container_fragment(frags,vid,0);
      if( deleted ) {
	set_lab_fragment(frags,vid,AS_CGB_DELETED_FRAG);
      } else if( (AS_CGB_BRANCHMULTICONT_FRAG == ilab) ) {
        // a protected fragment type.
      } else if ( contained ) {
	Tlab ilabel = AS_CGB_UNPLACEDCONT_FRAG;
        set_lab_fragment(frags,vid,ilabel);
      } else if ( AS_CGB_HANGING_CRAPPY_FRAG == ilab ) {
        // a protected fragment type in aggressive spur removal.
      } else {
        
	Tlab ilabel = ilab;

        int prefix_blunt = (inpx == 0) && (jnpx == 0);
        int suffix_blunt = (insx == 0) && (jnsx == 0);

        if( (inpx != 0) && (jnpx != 0) ) {
          const IntFragment_ID iid = get_iid_fragment(frags,vid);
          fprintf(stderr,"Woops: iid=" F_IID " vid=" F_IID " inpx=%d jnpx=%d\n",
                  iid,vid,inpx,jnpx);
        }
        if( (insx != 0) && (jnsx != 0) ) {
          const IntFragment_ID iid = get_iid_fragment(frags,vid);
          fprintf(stderr,"Woops: iid=" F_IID " vid=" F_IID " insx=%d jnsx=%d\n",
                  iid,vid,insx,jnsx);
        }

        assert( (inpx == 0) || (jnpx == 0) );
        assert( (insx == 0) || (jnsx == 0) );
        // If there is a intrachunk edge then there must not be an
        // interchunk edge.

	ilabel = ( (prefix_blunt && suffix_blunt) ?
		   AS_CGB_SOLO_FRAG : ilabel );
        
        if( (prefix_blunt && (! suffix_blunt)) ||
            (suffix_blunt && (! prefix_blunt)) ) {
	  ilabel = ( ((inpx == 1) && (insx == 0)) ||
		     ((inpx == 0) && (insx == 1)) ? 
		     AS_CGB_HANGING_CHUNK_FRAG : ilabel);
	  ilabel = ( ((inpx == 0) && (insx == 0)) ?
		     AS_CGB_HANGING_FRAG : ilabel);
        }

        if( (! prefix_blunt) && (! suffix_blunt) ) {

	  ilabel = AS_CGB_THRU_FRAG;
	  ilabel = ( ((inpx == 1) && (insx == 1)) ?
		     AS_CGB_INTRACHUNK_FRAG : ilabel);
	  ilabel = ( ((inpx == 0) && (insx == 1)) ||
		     ((inpx == 1) && (insx == 0)) ? 
		     AS_CGB_INTERCHUNK_FRAG : ilabel);
	  ilabel = ( ((inpx == 0) && (insx == 0)) ?
		     AS_CGB_THRU_FRAG : ilabel);
	}

	set_lab_fragment(frags,vid,ilabel);
      }
    }
  }

  // check_edges(frags, edges);
  SAFE_FREE(intrachunk_npx);
  SAFE_FREE(intrachunk_nsx);
  SAFE_FREE(interchunk_npx);
  SAFE_FREE(interchunk_nsx);
}

void chunk_classification_dvt
(
 Tfragment frags[],
 Tedge edges[],
 const int walk_depth,
 const int remove_blizzard_overlaps
 )
{
  time_t tp1 = 0, tp2;
  fprintf(stderr,"** Entered chunk_classification_dvt\n");
  fprintf(stderr,"**   remove_blizzard_overlaps=%d\n",remove_blizzard_overlaps);
  
    if(TIMINGS) {
      tp1 = time(NULL); 
      fprintf(stderr,
	      "Begin making the global fragment assignment\n");
    } 
#if 0
    remap_edge_labels( frags, edges);
    count_fragment_and_edge_labels( frags, edges,
                                    "after remap_edge_labels");
    // Remap the edges to a near FGB state.  The exceptions being
    // AS_TOUCHES_CRAPPY_DVT, AS_TOUCHES_CRAPPY_TOC, and
    // AS_CGB_MARKED_BY_BRANCH_DVT.
#endif

    // Needs a global thread synchronization here.
    // DANGER: I still need to make this incremental data safe.
    if(TIMINGS) {
      tp2 = time(NULL); 
      fprintf(stderr,
	      "%10" F_TIME_TP " sec: Finished making global "
	      "fragment assignment.\n",(tp2-tp1));
    }

#ifdef USE_BUDDYCHUNK_EGDES
    buddychunk_classification_of_overlaps
      (frags, edges, walk_depth);
    count_fragment_and_edge_labels( frags, edges,
                                    "after buddychunk_classification_of_overlaps");
#endif // USE_BUDDYCHUNK_EGDES
    
    intrachunk_classification_of_overlaps
      (frags, edges,
       walk_depth);
    count_fragment_and_edge_labels( frags, edges,
                                    "after intrachunk_classification_of_overlaps");
    if(TIMINGS) {
      tp2 = time(NULL);
      fprintf(stderr,
	      "%10" F_TIME_TP " sec: Finished chunk_classification_of_overlaps\n",
	      (tp2-tp1));
    }
#ifdef DEBUGVIEW
    view_fgb_chkpnt( "chunk_classification_of_overlaps",
                     frags, edges); 
#endif
    if(TIMINGS) {
      tp1 = time(NULL); 
      fprintf(stderr,"Begin chunk_classification_of_fragments\n");
    }
    chunk_classification_of_fragments(frags, edges);
    if(TIMINGS) {
      tp2 = time(NULL);
      fprintf(stderr,
	      "%10" F_TIME_TP " sec: Finished chunk_classification_of_fragments\n",
	      (tp2-tp1));
    }
    count_fragment_and_edge_labels( frags, edges,
                                    "after chunk_classification_of_fragments");
#ifdef DEBUGVIEW
    view_fgb_chkpnt( "chunk_classification_of_fragments",
                     frags, edges); 
#endif

  fprintf(stderr,"** Exit chunk_classification_dvt\n");
}
