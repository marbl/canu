
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

static char CM_ID[] = "$Id: AS_CGB_classify.c,v 1.8 2007-07-19 09:50:27 brianwalenz Exp $";

//  Description: These routines classify the fragments and overlaps in
//  a local manner as to their status in chunks.
//
//  Author: Clark Mobarry

#include "AS_CGB_all.h"

#undef DEBUGVIEW


static void maskout_overlaps_touching_crappy_fragments(Tfragment * frags,
                                                       Tedge     *edges) {
  //  Label each dovetail edge to and from every crappy fragment.

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
  fragment_npx = safe_malloc(sizeof(int) * nfrag);
  fragment_nsx = safe_malloc(sizeof(int) * nfrag);

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


      switch(ines) {
      case AS_CGB_INTERCHUNK_EDGE:
      case AS_CGB_INTRACHUNK_EDGE:
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
  safe_free(fragment_npx);
  safe_free(fragment_nsx);
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

  intrachunk_npx = safe_malloc(sizeof(char) * nfrag); // boolean
  intrachunk_nsx = safe_malloc(sizeof(char) * nfrag);
  interchunk_npx = safe_malloc(sizeof(char) * nfrag); // boolean
  interchunk_nsx = safe_malloc(sizeof(char) * nfrag);

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

  safe_free(intrachunk_npx);
  safe_free(intrachunk_nsx);
  safe_free(interchunk_npx);
  safe_free(interchunk_nsx);
}

void chunk_classification_dvt(Tfragment frags[],
                              Tedge edges[],
                              const int walk_depth) {

  intrachunk_classification_of_overlaps(frags, edges, walk_depth);
  //count_fragment_and_edge_labels( frags, edges, "after intrachunk_classification_of_overlaps");
  //view_fgb_chkpnt( "chunk_classification_of_overlaps", frags, edges); 

  chunk_classification_of_fragments(frags, edges);
  //count_fragment_and_edge_labels( frags, edges, "after chunk_classification_of_fragments");
  //view_fgb_chkpnt( "chunk_classification_of_fragments", frags, edges); 
}
