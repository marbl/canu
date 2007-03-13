
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
= "$Id: AS_CGB_fgb.c,v 1.9 2007-03-13 03:03:46 brianwalenz Exp $";
/*********************************************************************
 *
 * Module: AS_CGB_fgb.c
 * Description: The fragment overlap graph builder.
 * Assumptions:
 * 
 * (1) A unique OVL record appears only once in the input stream.
 * 
 * (2) A unique OFG record appears only once in the input stream.
 * 
 * (3) If an edge is in a chunk, then it is in no other chunks.
 * 
 * (4) If an essential fragment is in a chunk, then it is in no other chunks.
 * 
 * (5) The graph traversal subroutine, as_graph_traversal(), assumes
 * that circular chunks do not occur.
 *
 * Author: Clark Mobarry
 *********************************************************************/

/*************************************************************************/
/* System include files */

/*************************************************************************/
/* Local include files */
#include "AS_CGB_all.h"
//#include "AS_FGB_buddy.h"
//#include "AS_CGB_blizzard.h"      
#include "AS_FGB_hanging_fragment.h"
#include "AS_FGB_contained.h"

/****************************************************************************/
/* File Scope Globals */

static int TIMINGS = TRUE;
static int TIMINGS1 = FALSE;

/* Conditional compilation */
#define IN_PLACE_PERMUTATION
#undef SKIP_DECORDING

/* Some run time optimizations... */
#undef SMALLER_ADJACENCY_LIST

#undef REORDERING1

#undef SELF_CHECKING_CODE
#undef DEBUGGING

#define CHECK92
#ifdef SELF_CHECKING_CODE
#define CHECK80
#define CHECK90
#define CHECK96
#define CHECK97
#define CHECK98
#define CHECK99
#endif // SELF_CHECKING_CODE

#ifdef DEBUGGING
#undef DEBUG_SHORTCUT1
#undef DEBUG_SHORTCUT2
#undef DEBUG_SHORTCUT3
#define DEBUG_SHORTCUT4
#define DEBUG_SHORTCUT5
#undef DEBUG_INPUT
#undef DEBUG21
#undef DEBUG22
#undef DEBUG23
#undef DEBUG24
#undef DEBUG25
#define DEBUG26
#undef DEBUG27
#endif /*DEBUGGING*/

/*************************************************************************/
/* Constant definitions; Macro definitions; type definitions */

#define FRAGMENT_NOT_FOUND (~((IntFragment_ID)(0)))

/*************************************************************************/
/* Static Globals */

/*************************************************************************/
/* Function prototypes for internal static functions */

/**************************************************************/
/* Utility functions */


static void in_place_permute
(
 const size_t  ndata,  // The number of data records.
 const size_t  size,   // The size of the data records in bytes.
 const IntRank rank[], // The destination rank of each data record
 // (zero based) for a scatter-permutation.
 void  *       data    // The data.
)
{ // Inspired by Gene^s in-place permute ...
  char * done = NULL;
  char  saved_source[size], saved_target[size];
  size_t ii, jj=0;

  done = safe_malloc(sizeof(char) * ndata);
  
  fprintf(stderr,"Permutation in-place " F_SIZE_T " items of " F_SIZE_T " bytes\n",
          ndata, size);
#if 1
  memset(done,FALSE,ndata);
#else
  for( ii = 0; ii < ndata; ii++ ) {
    done[ii] = FALSE;
  }
#endif  
  for( ii = 0; ii < ndata; ii++ ) {
    if( ! done[ii] ) { 
      size_t it;
      // fprintf(stderr,"Start a permutation cycle at ii=" F_SIZE_T "\n",ii);
      it = ii;
      memcpy(saved_source,((char *)data)+ii*size,size);
      do {
        
	it = rank[it];
	// assert(it >= 0);
	assert(it < ndata);
	// fprintf(stderr," it=" F_SIZE_T "\n",it);
	assert( FALSE == done[it] );
	memcpy(saved_target,((char *)data)+it*size,size);
	// save the data at the target.
	memcpy(((char *)data)+it*size,saved_source,size); 
	// push the data.
	done[it] = TRUE;
	jj++;
	memcpy(saved_source,saved_target,size);
      } while( ii != it );
    }
  }
  assert(jj == ndata); // Was this a permutation??
  safe_free(done);
}


#if 0
static void check_edges1(/* Input Only */
		 Tfragment frags[], 
		 Tedge edges[],
		 int scrub)
{
  /* This routine is for checks overlaps not assuming adjacency
     list sorted edges.
  */
  //const IntFragment_ID nfrag = GetNumFragments(frags);
  const IntEdge_ID nedge = GetNumEdges(edges);
  //IntFragment_ID iv0;
  //int is0;
  IntEdge_ID ie0;
  time_t tp1 = 0,tp2;

  if(TIMINGS) {
    tp1 = time(NULL); fprintf(stderr,"Begin check_edges1\n");
    system_top();
  }

  for(ie0=0; ie0 < nedge; ie0++) {
	const IntFragment_ID avx = get_avx_edge(edges,ie0);
	const IntFragment_ID bvx = get_bvx_edge(edges,ie0);
	const IntFragment_ID afr = get_iid_fragment(frags,avx);
	const IntFragment_ID bfr = get_iid_fragment(frags,bvx);
	const int asx = get_asx_edge(edges,ie0);
	const int bsx = get_bsx_edge(edges,ie0);
	const int ahg = get_ahg_edge(edges,ie0);
	const int bhg = get_bhg_edge(edges,ie0);
	const Tnes ines = get_nes_edge(edges,ie0);
	const uint32 qua = get_qua_edge(edges,ie0);
	//const int inv = get_inv_edge(edges,ie0);
	int bad_edge = FALSE;

	if( (avx == bvx)&&(asx == bsx) ) { 
	  fprintf(stderr,F_IID ": " F_IID " " F_IID " " F_IID " " F_IID " %d %d %d %d %d\n",
		  ie0,
		  afr,bfr,avx,bvx,
		  asx,bsx,
		  ahg,bhg,
		  ines);
	  fprintf(stderr,"check_edges1 found an impossible overlap:\n");
	  assert(FALSE);
	}

	switch(ines) {
	  /* The dovetail overlaps */
	case AS_CGB_INTERCHUNK_EDGE:
	case AS_CGB_INTRACHUNK_EDGE:
	case AS_CGB_TOUCHES_CONTAINED_EDGE:
	case AS_CGB_BETWEEN_CONTAINED_EDGE:
	case AS_CGB_REMOVED_BY_TRANSITIVITY_DVT:
	case AS_CGB_REMOVED_BY_THRESHOLD_DVT:
	case AS_CGB_MARKED_BY_BRANCH_DVT:
          // get_dvt_edge()
	  if(!is_a_dvt_simple(ahg,bhg)) {
	    bad_edge = TRUE;
	  }
	  break;
	case AS_CGB_CONTAINED_EDGE:
	case AS_CGB_REMOVED_BY_TRANSITIVITY_CON:
	case AS_CGB_REMOVED_BY_THRESHOLD_CON:
	  if(
#if 1
             is_a_dvt_simple(ahg,bhg)
#else             
             ! (is_a_frc_simple(ahg,bhg) || is_a_toc_simple(ahg,bhg) || is_a_dgn_simple(ahg,bhg))
#endif             
             ) {
	    bad_edge = TRUE;
	  }
	  break;

	case AS_CGB_REMOVED_BY_DUPLICATE_DVT:
	case AS_CGB_REMOVED_BY_DUPLICATE_CON:
	  break;
	default:
	  break;
	}

	if( bad_edge ) {
	  fprintf(stderr,
		  F_IID ": " F_IID " " F_IID " " F_IID " " F_IID " %d %d %d %d %d %d\n"
		  ie0,
		  afr,bfr,avx,bvx,
		  asx,bsx,
		  ahg,bhg,
		  ines,qua);
	}
  }
  if(TIMINGS) {
    tp2 = time(NULL); 
    fprintf(stderr,"%10" F_TIME_TP " sec: Finished check_edges1\n",(tp2-tp1));
    system_top();
  }
}
#endif

#if 0
static void reorder_fragments
(
 Tfragment frags[],
 IntRank   fragment_rank[],
 Tedge     edges[])
{
  const IntFragment_ID nfrag = GetNumFragments(frags);
  const IntEdge_ID nedge = GetNumEdges(edges);

#ifdef CHECK96
  {
    IntFragment_ID ivold;
    IntRank * fragment_tmp = NULL;
    fragment_tmp = safe_malloc(sizeof(IntRank) * nfrag);

    /* Begin validating the permutation. */
    // #pragma omp parallel for
    for(ivold=0; ivold<nfrag; ivold++){
      fragment_tmp[ivold] = FRAGMENT_NOT_VISITED;
    }
    // #pragma omp parallel for
    for(ivold=0; ivold<nfrag; ivold++){
      fragment_tmp[fragment_rank[ivold]] = ivold;
    }
    // #pragma omp parallel for
    for(ivold=0; ivold<nfrag; ivold++){
      assert(fragment_tmp[ivold] != FRAGMENT_NOT_VISITED);
    }
    /* End validating the permutation. */
    safe_free(fragment_tmp);
  }
#endif // CHECK96

#ifdef CHECK80
  {
    const IntFragment_ID max_frag_iid = 1000000;
    check_for_fragment_corruption( frags, edges, max_frag_iid );
  }
#endif

  /* re-order the fragments. */
  /* move the fragment data */
  in_place_permute
    ( nfrag, sizeof(Afragment), 
      fragment_rank, GetVA_Afragment(frags,0));

#ifdef CHECK80
  {
    const IntFragment_ID max_frag_iid = 1000000;
    check_for_fragment_corruption( frags, edges, max_frag_iid );
  }
#endif

  {
    /* Remap the fragment references from IIDs to VIDs (fragment
       identifiers). */
    // #pragma omp parallel for
    IntEdge_ID iedge;
    for(iedge=0;iedge<nedge;iedge++){
      set_avx_edge(edges,iedge, 
		   (IntFragment_ID)fragment_rank[get_avx_edge(edges,iedge)]);
      set_bvx_edge(edges,iedge,
		   (IntFragment_ID)fragment_rank[get_bvx_edge(edges,iedge)]);
    }
  }
  
  /* The edge data is ready to be resorted for memory locality. */
}
#endif

static void setup_segments
( /* Modify */
 Tfragment frags[],
 Tedge edges[])
{
  const IntFragment_ID nfrag = GetNumFragments(frags);
  const IntEdge_ID nedge = GetNumEdges(edges);

  /* Set up the segmentation information as vertex data.  This routine
     assumes that the edge records are already sorted by
     (axv,asx). Thus if the fragments have been reordered then the
     edges need to have been reordered before this routine. */
  int nrec_used;
  IntFragment_ID min_frag_vid; /* The minimum fragment iid used. */
  IntFragment_ID max_frag_vid; /* The maximum fragment iid used. */
  IntFragment_ID num_frag_used; /* The number of fragments used. */
  int max_nnode; /* The maximum number of edges adjacent to a vertex. */

  assert(nfrag > 0);
  // #pragma omp for
  { IntFragment_ID ifrag;
  for(ifrag=0;ifrag<nfrag;ifrag++) {
      set_seglen_vertex(frags,ifrag,FALSE,0);
      set_seglen_vertex(frags,ifrag,TRUE,0);
  }}
  max_frag_vid = min_frag_vid = get_avx_edge(edges,0);
  num_frag_used = 0; /* the number of vertices used by the edges. */
  max_nnode = 0;

#ifndef USE_A_SCAN
  assert(nedge > 0);
  { IntEdge_ID iedge;
  for(iedge=0;iedge<nedge;iedge++){
    const IntFragment_ID iavx = get_avx_edge(edges,iedge);
    const int iasx = get_asx_edge(edges,iedge);
    // const IntFragment_ID ibvx = get_bvx_edge(edges,iedge);
    // const int ibsx = get_bsx_edge(edges,iedge);
    const IntEdge_ID iseglen = get_seglen_vertex(frags,iavx,iasx) + 1;

    // printf(F_IID " : " F_IID " %d " F_IID " %d \n",iedge,iavx,iasx,ibvx,ibsx);
    // assert( iavx >= 0 );
    // assert( ibvx >= 0 );
    assert( iavx < nfrag );
    //assert( ibvx < nfrag );

    min_frag_vid = MIN(min_frag_vid,iavx);
    max_frag_vid = MAX(max_frag_vid,iavx);
    set_seglen_vertex(frags,iavx,iasx,iseglen);
  }}

  // Now assume that the edges are sorted by the (avx,asx) vertex.
  {
    IntFragment_ID ifrag = 0;
    IntEdge_ID isegstart = 0;

    for(ifrag=0;ifrag<nfrag;ifrag++) { /* Note the minus one. */
      int isuffix;
      if((get_seglen_vertex(frags,ifrag,FALSE) > 0 ) ||
	 (get_seglen_vertex(frags,ifrag,TRUE) > 0 )) { 
	num_frag_used ++;
      }
      for(isuffix=0;isuffix<2;isuffix++) {
	/* The number of fragments used in the overlap graph. */
	max_nnode = MAX(max_nnode,get_seglen_vertex(frags,ifrag,isuffix));
	/* The prefix and suffix segments are interleaved. */
	set_segstart_vertex(frags,ifrag,isuffix,isegstart);
	isegstart += get_seglen_vertex(frags,ifrag,isuffix);
      }
    }
    assert(isegstart == nedge); 
  }

#else /* USE_A_SCAN */
  /* A scan through the list of records. */
  iavx = FRAGMENT_NOT_VISITED;  /* put the initial value out of range */
  ibvx = FRAGMENT_NOT_VISITED;
  /* #pragma omp for private(iedge)  */
  for(iedge=0;iedge<nedge;iedge++){
    IntFragment_ID last_avx,last_bvx;
    int nnode;
    last_avx = iavx;
    last_bvx = ibvx;
    iavx = get_avx_edge(edges,iedge);
    ibvx = get_bvx_edge(edges,iedge);
    assert( iavx >= 0 );
    assert( ibvx >= 0 );
    assert( iavx < nfrag );
    assert( ibvx < nfrag );
    if( iavx != last_avx || iasx != last_asx ) {
      fprintf(stderr,"last_avx,iavx=" F_IID "," F_IID "\n",last_avx,iavx);
      assert( iavx == last_avx+1);
      set_segstart_vertex(frags,iavx,iasx,iedge);
      (*nfrag_used)++;
      (*min_ifrag) = MIN((*min_ifrag),iavx);
      (*max_ifrag) = MAX((*max_ifrag),iavx);
      last_bvx = FRAGMENT_NOT_VISITED;
      nnode = 0;
    }
    if( ibvx != last_bvx) {
      nnode++;
      set_seglen_vertex(frags,iavx,iasx,nnode);
    }
    (*max_nnode) = MAX((*max_nnode),nnode);
  }
  set_segstart_vertex(frags,(*nfrag_used),nedge,iasx);
#endif /* USE_A_SCAN */

  // #pragma omp barrier

  { /* Check the validity of the segments. */
    IntFragment_ID ifrag;
    nrec_used = 0;
    for(ifrag=0;ifrag<nfrag;ifrag++) 
      {
	nrec_used += get_seglen_vertex(frags,ifrag,FALSE);
	nrec_used += get_seglen_vertex(frags,ifrag,TRUE);
      }
    // fprintf(stderr,"nrec_used=%d,nedge=" F_IID "\n",nrec_used,nedge);
  }
  fprintf(stderr,
	  "Information on the fragments refered to by overlaps.\n"
	  "nfrag=" F_IID ",min_frag_vid=" F_IID ",max_frag_vid=" F_IID ","
	  "nfrag_used=" F_IID ",max_nnode=%d\n",
	  nfrag,min_frag_vid,max_frag_vid,num_frag_used,max_nnode);
  assert(nrec_used == nedge);
  // assert(min_frag_vid >= 0);
  assert((nfrag == 0) || (max_frag_vid < nfrag));
} // static void setup_segments( Tfragment frags[], Tedge edges[])


static void pack_the_edges
(
 Tfragment *frags,
 Tedge     *edges,
 TIntEdge_ID *next_edge_obj
)
{
  // Compress the removed edges from the graph.

  const IntEdge_ID nedge = GetNumEdges(edges);
  IntEdge_ID iold, idup=0;
  IntEdge_ID
    n_AS_CGB_UNUSED_EDGE=0,
    n_AS_CGB_REMOVED_BY_TRANSITIVITY_DVT=0,
    n_AS_CGB_REMOVED_BY_THRESHOLD_DVT=0,
    n_AS_CGB_REMOVED_BY_DUPLICATE_DVT=0,
    n_AS_CGB_REMOVED_BY_TRANSITIVITY_CON=0,
    n_AS_CGB_REMOVED_BY_THRESHOLD_CON=0,
    n_AS_CGB_REMOVED_BY_DUPLICATE_CON=0;
  int old_dovetail_degree = 0, new_dovetail_degree = 0;
  int old_containment_degree = 0, new_containment_degree = 0;

  IntFragment_ID old_avx = 0;
  IntFragment_ID old_asx = 3; // Out of range
  IntFragment_ID new_avx = 0;
  IntFragment_ID new_asx = 0;
  
  verify_that_the_edges_are_in_order(edges);

  for(iold=0;iold<nedge;iold++) {
    const Tnes ines = get_nes_edge(edges,iold);
    
#if 1
    new_avx = get_avx_edge(edges,iold);
    new_asx = get_asx_edge(edges,iold);
    if( (old_avx != new_avx) || (old_asx != new_asx) ) {
      // When the avx,asx change the segment changes.
      
      if( new_dovetail_degree > old_dovetail_degree ) {
        assert(FALSE);
      }

      if( (new_dovetail_degree == 0) && (old_dovetail_degree > 0) ) {
        fprintf(stderr,
                "WARNING: fragment-end became dovetail disconnected iid=" F_IID " vid=" F_IID " suf=%d\n",
                get_iid_fragment(frags,old_avx), old_avx, old_asx);
        //assert(FALSE);
      }
      // We are inspecting to see if a fragment-end's dovetail degree
      // changes from positive to zero.
      
      if( (new_containment_degree == 0) && (old_containment_degree > 0) ) {
        fprintf(stderr,
                "WARNING: fragment-end became containment disconnected iid=" F_IID " vid=" F_IID " suf=%d\n",
                get_iid_fragment(frags,old_avx), old_avx, old_asx);
        //assert(FALSE);
      }

      old_dovetail_degree = 0; new_dovetail_degree = 0;
      old_containment_degree = 0; new_containment_degree = 0;
      // Reset the counters because we are in a new segment.
    }
    old_avx = new_avx; old_asx = new_asx;
#endif
      
    // Pack the edge to the list.
    switch(ines) {
    case AS_CGB_UNUSED_EDGE:
      n_AS_CGB_UNUSED_EDGE++; break;
      
    case AS_CGB_REMOVED_BY_TRANSITIVITY_DVT:
      n_AS_CGB_REMOVED_BY_TRANSITIVITY_DVT++;
      old_dovetail_degree++; break;
    case AS_CGB_REMOVED_BY_THRESHOLD_DVT:
      n_AS_CGB_REMOVED_BY_THRESHOLD_DVT++;
      old_dovetail_degree++; break;
    case AS_CGB_REMOVED_BY_DUPLICATE_DVT:
      n_AS_CGB_REMOVED_BY_DUPLICATE_DVT++;
      old_dovetail_degree++; break;

    case AS_CGB_REMOVED_BY_TRANSITIVITY_CON:
      n_AS_CGB_REMOVED_BY_TRANSITIVITY_CON++;
      old_containment_degree++; break;
    case AS_CGB_REMOVED_BY_THRESHOLD_CON:
      n_AS_CGB_REMOVED_BY_THRESHOLD_CON++;
      old_containment_degree++; break;
    case AS_CGB_REMOVED_BY_DUPLICATE_CON:
      n_AS_CGB_REMOVED_BY_DUPLICATE_CON++;
      old_containment_degree++; break;
      
    case AS_CGB_DOVETAIL_EDGE:
    case AS_CGB_THICKEST_EDGE:
    case AS_CGB_INTERCHUNK_EDGE:
    case AS_CGB_BUDDYCHUNK_EDGE:
    case AS_CGB_INTRACHUNK_EDGE:
    case AS_CGB_TOUCHES_CONTAINED_EDGE:
    case AS_CGB_BETWEEN_CONTAINED_EDGE:
    case AS_CGB_TOUCHES_CRAPPY_DVT:
      {
        old_dovetail_degree++;
        new_dovetail_degree++;
        if (idup != iold)
          memcpy( GetVA_Aedge(edges,idup),
                  GetVA_Aedge(edges,iold),
                  sizeof(Aedge));
	idup++;
      }
      break;

    case AS_CGB_CONTAINED_EDGE:
    case AS_CGB_TOUCHES_CRAPPY_CON:
      {
        old_containment_degree++;
        new_containment_degree++;
        if (idup != iold)
          memcpy( GetVA_Aedge(edges,idup),
                  GetVA_Aedge(edges,iold),
                  sizeof(Aedge));
	idup++;
      }
      break;

    default:
      fprintf(stderr,"Unexpected edge label = %d\n", ines);
      assert(FALSE);
      break;
    }
  }
  ResetToRange_Aedge(edges,idup);
  if(nedge != idup) {
    ResetToRangeVA_IntEdge_ID(next_edge_obj,0);
    // The next_edge_obj data is invalidated.
  }
  setup_segments(/* Modify */ frags, edges);

  fprintf(stderr,
	  "Pack the edges:\n"
	  "raw_nedge=" F_IID " reduced_nedge=" F_IID "\n"
          "n_AS_CGB_UNUSED_EDGE                =%15" F_IIDP "\n"
	  "n_AS_CGB_REMOVED_BY_TRANSITIVITY_DVT=%15" F_IIDP "\n"
	  "n_AS_CGB_REMOVED_BY_THRESHOLD_DVT   =%15" F_IIDP "\n"
	  "n_AS_CGB_REMOVED_BY_DUPLICATE_DVT   =%15" F_IIDP "\n"
	  "n_AS_CGB_REMOVED_BY_TRANSITIVITY_CON=%15" F_IIDP "\n"
	  "n_AS_CGB_REMOVED_BY_THRESHOLD_CON   =%15" F_IIDP "\n"
	  "n_AS_CGB_REMOVED_BY_DUPLICATE_CON   =%15" F_IIDP "\n"
          ,
	  nedge, idup,
          n_AS_CGB_UNUSED_EDGE,
	  n_AS_CGB_REMOVED_BY_TRANSITIVITY_DVT,
	  n_AS_CGB_REMOVED_BY_THRESHOLD_DVT,
	  n_AS_CGB_REMOVED_BY_DUPLICATE_DVT,
	  n_AS_CGB_REMOVED_BY_TRANSITIVITY_CON,
	  n_AS_CGB_REMOVED_BY_THRESHOLD_CON,
	  n_AS_CGB_REMOVED_BY_DUPLICATE_CON
          );
}

static void rebuild_next_edge_obj
(
  Tedge       * edges,
  TIntEdge_ID * next_edge_obj
 )
{
  fprintf(stderr, "WARNING "__FILE__ " line %d : rebuild_next_edge_obj() is only a stub.\n", __LINE__ );
  //assert(FALSE);
}

#if 0
static void reorder_edges_qsort_all
( Tfragment frags[],
  Tedge edges[] 
  )
{ // Sort the edges....
  const IntEdge_ID nedge = GetNumEdges(edges);
  time_t tp1 = 0,tp2;
  QsortCompare compare_edge = get_compare_edge_function();

fprintf(stderr,"Sort the edges by (avx,asx,ahg,bvx,bsx,bhg).\n");
  if(TIMINGS) {
    tp1 = time(NULL); fprintf(stderr,"Begin qsort\n");
    system_top();
  }
  qsort(GetVA_Aedge(edges,0),nedge,sizeof(Aedge),
	(int (*)(const void *,const void *))compare_edge); /* FIX this */
  /* The internal structure of the edges object of type Tedge is
     used here.  We are using the fact that the first address of the
     edges object is a pointer to an array of type Aedge and length
     nedges.*/
  if(TIMINGS) {
    tp2 = time(NULL); 
    fprintf(stderr,"%10" F_TIME_TP " sec: Finished qsort\n",(tp2-tp1));
    system_top();
  }
  setup_segments(/* Modify */ frags, edges);
}
#endif

static void reorder_edges_bin_and_qsort
( Tfragment *frags,
  Tedge *edges,
  TIntEdge_ID *next_edge_obj
)
{
  QsortCompare compare_edge = get_compare_edge_function();
  const IntFragment_ID nfrag = GetNumFragments(frags);
  const IntEdge_ID nedge = GetNumEdges(edges);

  ResetToRangeVA_IntEdge_ID(next_edge_obj,0);
  // The next_edge_obj data is invalidated.
  
  if( nedge > 0) { // Sort the edges ......
    const int max_frag_len = 2048;
    const IntEdge_ID max_nbins = MAX(max_frag_len,2*nfrag);
    time_t tp1 = 0,tp2;
 
    IntEdge_ID * seglen = NULL;
    IntEdge_ID * segstart = NULL;
    IntRank    * edge_rank = NULL;

    assert(nedge != 0);
    assert(max_nbins != 0);

    seglen    = safe_malloc(sizeof(IntEdge_ID) * max_nbins);
    segstart  = safe_malloc(sizeof(IntEdge_ID) * max_nbins);
    edge_rank = safe_malloc(sizeof(IntRank) * nedge);

    { // FRAGMENT-END SORT
      if(TIMINGS) {
	fprintf(stderr,"Sort the edges by avx,asx\n");
	tp1 = time(NULL); fprintf(stderr,"Begin bin rank\n");
        system_top();
      }
      { // Sort all of the edges ....
	const IntEdge_ID ie_start  = 0;
	const IntEdge_ID ie_finish = nedge + ie_start;

	{ IntFragment_ID iv0; int is0;
	for(iv0=0;iv0<nfrag;iv0++) for(is0=0;is0<2;is0++) {
	  const IntFragment_ID ivert = 2*iv0+is0;
	  seglen[ivert] = 0;
	  segstart[ivert] = 0;
	  set_seglen_vertex(frags,iv0,is0,0);
	  set_segstart_vertex(frags,iv0,is0,0);
	}}
	
	if(TIMINGS1) {
	  tp2 = time(NULL); 
	  fprintf(stderr,"%10" F_TIME_TP " sec: so far at 1\n",(tp2-tp1));
          system_top();
	}

	// #pragma omp parallel for
	// Count the number of edges at each fragment-end.
	{ IntEdge_ID ie;
	for(ie=ie_start;ie<ie_finish;ie++) {
	  const IntFragment_ID iavx = get_avx_edge(edges,ie);
	  const int iasx = get_asx_edge(edges,ie);
	  const IntFragment_ID ivert = 2*iavx+iasx;
	  const int nnode = seglen[ivert] + 1;
	  seglen[ivert] = nnode;
	  set_seglen_vertex(frags,iavx,iasx,nnode);
	}}
	
	if(TIMINGS1) {
	  tp2 = time(NULL); 
	  fprintf(stderr,"%10" F_TIME_TP " sec: so far at 2\n",(tp2-tp1));
          system_top();
	}

	// Scan the number of edges at each fragment-end to find the
	// segment start for each fragment end.
	{ 
	  IntFragment_ID iv0; int is0; IntEdge_ID isum=0;
	  for(iv0=0;iv0<nfrag;iv0++) for(is0=0;is0<2;is0++) {
	    const IntFragment_ID ivert = 2*iv0 + is0;
	    const int nnode = seglen[ivert];
	    segstart[ivert] = isum;
	    set_segstart_vertex(frags,iv0,is0,isum);
	    isum += nnode;
	    seglen[ivert] = 0;
	  }
	  assert( isum == (ie_finish - ie_start) );
	}
	
	if(TIMINGS1) {
	  tp2 = time(NULL); 
	  fprintf(stderr,"%10" F_TIME_TP " sec: so far at 3\n",(tp2-tp1));
          system_top();
        }

	// #pragma omp parallel for
	// Assign a rank to each edge.
	{ IntEdge_ID ie;
	for(ie=ie_start;ie<ie_finish;ie++) {
	  const IntFragment_ID iavx = get_avx_edge(edges,ie);
	  const int iasx = get_asx_edge(edges,ie);
	  const IntFragment_ID ivert = 2*iavx + iasx;
	  const int nnode = seglen[ivert];
	  const IntEdge_ID start = segstart[ivert];
	  edge_rank[ie] = ie_start + start + nnode;
	  seglen[ivert] = nnode + 1;
	}}
      } // Sort all of the edges...

      if(TIMINGS) {
	tp2 = time(NULL); 
	fprintf(stderr,"%10" F_TIME_TP " sec: Finished bin rank\n",(tp2-tp1));
        system_top();
      }
      
#ifdef CHECK99
      {
	IntEdge_ID * edge_index = NULL;
	IntEdge_ID ie;

	safe_malloc(edge_index, IntEdge_ID, nedge);
	for(ie=0;ie<nedge;ie++) {
	  edge_index[ie] = AS_FGB_EDGE_NOT_VISITED;
	}
	for(ie=0;ie<nedge;ie++) {
	  edge_index[edge_rank[ie]] = ie;
	}
	for(ie=0;ie<nedge;ie++) {
	  assert(edge_index[ie] != AS_FGB_EDGE_NOT_VISITED);
	  // assert(edge_index[ie] >= 0 );
	  assert(edge_index[ie] < nedge );
	}
	safe_free(edge_index);
      }
#endif // CHECK99

#ifdef IN_PLACE_PERMUTATION
      if(TIMINGS) {
	tp1 = time(NULL); fprintf(stderr,"Begin in-place permutation\n");
        system_top();
      }
      // Now permute the edges by the rank.
      in_place_permute
	( nedge, sizeof(Aedge), 
	  edge_rank, GetVA_Aedge(edges,0));
      if(TIMINGS) {
	tp2 = time(NULL); 
	fprintf(stderr,
		"%10" F_TIME_TP " sec: Finished in-place permutation\n",(tp2-tp1));
        system_top();
      }
#else
      if(TIMINGS) {
	tp1 = time(NULL); fprintf(stderr,"Begin VA permutation\n");
        system_top();
      }
      // Now permute the edges by the rank.
      ScatterInPlace_VA( edges, nedge, edge_rank);
      if(TIMINGS) {
	tp2 = time(NULL); 
	fprintf(stderr,
		"%10" F_TIME_TP " sec: Finished VA permutation\n",(tp2-tp1));
        system_top();
      }
#endif

    } // FRAGMENT-END SORT

    safe_free(edge_rank);
    safe_free(seglen);
    safe_free(segstart);
  }
  
   { // Sort the edges (in fragment-end segments)
     time_t tp1 = 0,tp2;
    
#ifdef DEBUGGING
     fprintf(stderr,"Sort the segment of edges by avx\n");
#endif
     if(TIMINGS) {
       tp1 = time(NULL); 
       fprintf(stderr,"Begin qsort of the fragment-end segments\n");
       system_top();
     }
     { IntFragment_ID iv0; int is0;
     for(iv0=0;iv0<nfrag;iv0++) for(is0=0;is0<2;is0++) {
       const IntEdge_ID ie_start 
	 = get_segstart_vertex(frags,iv0,is0);
       const IntEdge_ID nnode 
	 = get_seglen_vertex(frags,iv0,is0);
       qsort(GetVA_Aedge(edges,ie_start),nnode,sizeof(Aedge),
	     (int (*)(const void *,const void *))compare_edge); /* FIX this */
       /* The internal structure of the edges object of type Tedge is
	  used here.  We are using the fact that the first address of the
	  edges object is a pointer to an array of type Aedge and length
	  nedges.*/
     }}
     if(TIMINGS) { 
       tp2 = time(NULL); 
       fprintf(stderr,
	       "%10" F_TIME_TP " sec: Finished qsort of the fragment-end segments\n",
	       (tp2-tp1));
       system_top();
     }
   }

  rebuild_next_edge_obj( edges, next_edge_obj);
}

void reorder_edges
( Tfragment *frags,
  Tedge *edges,
  TIntEdge_ID *next_edge_obj
)
{
  reorder_edges_bin_and_qsort( frags, edges, next_edge_obj);
  setup_segments(/* Modify */ frags, edges);
}

/////////////////////////////////////////////////////////////////////

#define NBINS 100
void graph_locality_diagnostic
( Tfragment *frags, 
  Tedge     *edges,
  char      *named,
  char      *namec
  )
{
  time_t tp1 = 0, tp2;
  if(TIMINGS) {
    tp1 = time(NULL); fprintf(stderr,"Begin writing locality diagnostic.\n");
  }
  { 
    FILE *ffga = stderr;
    IntEdge_ID ie1;
    const IntFragment_ID nfrag = GetNumFragments(frags);
    const IntEdge_ID nedge = GetNumEdges(edges);
    const int nsample=500;
    const int nbucket=500;
    int twod[NBINS][NBINS] = {{0}};
    int twoc[NBINS][NBINS] = {{0}};
    
    Histogram_t 
      *edges_locality_histogram
      = create_histogram(nsample,nbucket,0,TRUE);
    for(ie1=0; ie1 < nedge; ie1++) {
      if((get_nes_edge(edges,ie1) != AS_CGB_REMOVED_BY_DUPLICATE_DVT) &&
         (get_nes_edge(edges,ie1) != AS_CGB_REMOVED_BY_DUPLICATE_CON)
         ) {
	const IntFragment_ID iv0 = get_avx_edge(edges,ie1);
        const IntFragment_ID iv1 = get_bvx_edge(edges,ie1);
        const int ahg = get_ahg_edge(edges,ie1);
        const int bhg = get_bhg_edge(edges,ie1);
        const int iv1_iv0 = iv1 - iv0;
        const int vdiff = (iv1_iv0 > 0 ? iv1_iv0 : -iv1_iv0);
        {
          int i0 = (NBINS * iv0) / nfrag;
          int i1 = (NBINS * iv1) / nfrag;
          i0 = MIN(i0,NBINS);
          i1 = MIN(i1,NBINS);
          if( is_a_dvt_simple(ahg,bhg) ) { 
            twod[i0][i1] ++;
          } else {
            twoc[i0][i1] ++;
          }
        }
        // fprintf(fdiag,F_IID " " F_IID "\n", iv0, iv1);
        add_to_histogram(edges_locality_histogram, vdiff, NULL);
      }
    }
    fprintf(ffga,"\n\nHistogram of "
	    "the locality of graph edges\n");
    print_histogram(ffga,edges_locality_histogram, 0, 1);
    free_histogram(edges_locality_histogram);

    {
      int i0,i1; 
      FILE * fdiagd = fopen(named,"w");
      assert(fdiagd != NULL);
      for(i1=0;i1<NBINS;i1++) {
        for(i0=0;i0<NBINS;i0++) {
          const int density = twod[i0][i1];
          if( density > 0) {
            fprintf(fdiagd,"%d %d %d\n", i0, i1, density);
          }
        }
      }
      fclose(fdiagd);
    }

    {
      int i0,i1; 
      FILE * fdiagc = fopen(namec,"w");
      assert(fdiagc != NULL);
      for(i1=0;i1<NBINS;i1++) {
        for(i0=0;i0<NBINS;i0++) {
          const int density = twoc[i0][i1];
          if( density > 0) {
            fprintf(fdiagc,"%d %d %d\n", i0, i1, density);
          }
        }
      }
      fclose(fdiagc);
    }
  }
  if(TIMINGS) {
    tp2 = time(NULL); fprintf(stderr,"%10" F_TIME_TP " sec: Finished "
                              "locality diagnostic.\n",
                              (tp2-tp1));
  }
}

#if 0
static void reorder_graph
(
 Tfragment *frags, 
 Tedge     *edges,
 TIntEdge_ID * next_edge_obj
 )
{
  const int output_graph_locality_diagnostic = FALSE;
  time_t tp1 = 0, tp2;
  const IntFragment_ID nfrag = GetNumFragments(frags);
  IntRank * fragment_rank = NULL;

  fragment_rank = safe_malloc((sizeof(IntRank) * nfrag); // for graph traversal
  
  if(output_graph_locality_diagnostic){
    graph_locality_diagnostic( frags, edges, "Before.diagd", "Before.diagc");
  }

  if(TIMINGS) { 
    tp1 = time(NULL); 
    fprintf(stderr,"Begin: Find a new ordering of the fragments in fragment_rank.\n");
    system_top();
  }
  as_graph_traversal(stderr, frags, edges, fragment_rank);
  if(TIMINGS) { 
    tp2 = time(NULL); 
    fprintf(stderr,"%10" F_TIME_TP " sec: "
            "Finished finding a new ordering of the fragments in fragment_rank.\n",
            (tp2-tp1));
    system_top();
  }

  if(TIMINGS) { 
    tp1 = time(NULL); 
    fprintf(stderr,"Begin: reorder_fragments.\n");
    system_top();
  }
  reorder_fragments( frags, fragment_rank, edges);
  if(TIMINGS) {
    tp2 = time(NULL); 
    fprintf(stderr,"%10" F_TIME_TP " sec: "
            "Finished reorder_fragments.\n",
            (tp2-tp1));
    system_top();
  }
  
  if(TIMINGS) {
    tp1 = time(NULL); 
    fprintf(stderr,"Begin: reorder_edges.\n");
    system_top();
  }
  reorder_edges(frags, edges, next_edge_obj);
  if(TIMINGS) {
    tp2 = time(NULL);
    fprintf(stderr,"%10" F_TIME_TP " sec: "
            "Finished reorder_edges.\n",
            (tp2-tp1));
    system_top();
  }
  
  if(output_graph_locality_diagnostic){
    graph_locality_diagnostic( frags, edges, "After.diagd", "After.diagc");
  }

#ifdef DEBUG_GRAPH_DIAGNOSTICS
  {
    char strtmp[2048] = {0};
    strcat(strtmp,"graph0p");
    fprintf(stderr,"graph_diagnostics\n");
    graph_diagnostics
      (strtmp, frags, edges, AS_CGB_INTERCHUNK,AS_CGB_THRU,0);
  }
#endif /*DEBUG_GRAPH_DIAGNOSTICS*/
  
  safe_free(fragment_rank);
}
#endif

/****************************************************************
 *
 * The main graph transformation routines
 *
 ****************************************************************/

#if 0
static void check_transitive_overlap_marking
(
 Tfragment frags[],
 Tedge     edges[])
{
  /* Check for overly aggressive transitive overlap removal. That is,
   * transitive overlap removal should NEVER remove all the overlaps
   * on a fragment end.  
   */

  const IntFragment_ID nfrag = GetNumFragments(frags);
  IntEdge_ID ie;
  IntFragment_ID ifrag;
  int is; 

  for(ifrag=0; ifrag<nfrag; ifrag++) {
    for(is=0; is<2; is++) {
      const IntEdge_ID segstart = get_segstart_vertex(frags,ifrag,is);
      const int seglength = get_seglen_vertex(frags,ifrag,is);
      int count = 0;
      for(ie=segstart; ie<segstart+seglength; ie++) {
	const Tnes ines = get_nes_edge(edges,ie);
	switch(ines) {
	  /* The containment overlaps... */
	case AS_CGB_CONTAINED_EDGE:
	  count ++;
	  break;
	  
	  /* The dovetail overlaps... */
	case AS_CGB_DOVETAIL_EDGE:
	case AS_CGB_THICKEST_EDGE:
	case AS_CGB_INTERCHUNK_EDGE:
	case AS_CGB_TOUCHES_CONTAINED_EDGE:
	case AS_CGB_BETWEEN_CONTAINED_EDGE:
	  count ++;
	  break;
	case AS_CGB_REMOVED_BY_TRANSITIVITY_DVT:
	case AS_CGB_REMOVED_BY_DUPLICATE_DVT:
	case AS_CGB_REMOVED_BY_TRANSITIVITY_CON:
	case AS_CGB_REMOVED_BY_DUPLICATE_CON:
	  break;
	default:
          fprintf(stderr,"Unexpected overlap edge type nes=%d\n",ines);
	  assert(FALSE);
	}
      }

      if( (seglength > 0) && 
	  (count == 0) &&
	  (! get_del_fragment(frags,ifrag))
	  //(AS_CGB_DELETED_FRAG != get_lab_fragment(frags,ifrag))
          ) {
	fprintf(stderr,"Disconnected fragment end " F_IID " %d\n",
		ifrag,is);
      }
    }
  }
}
#endif


#if 0
static void find_median_fragment_end_degree
(
 Tfragment frags[],
 Tedge edges[],
 const char Output_Graph_Store[],
 int *median_fragment_end_degree
 ) {
  const IntFragment_ID nfrag=GetNumFragments(frags);
  const IntEdge_ID nedge=GetNumEdges(edges);

  { // Find the median adjacency degree...
    int * bincount = NULL;
    int maxcount = 0;
    IntEdge_ID totcount = 0;
    { IntFragment_ID iv0; int is0;
    for(iv0=0;iv0<nfrag;iv0++) for(is0=0;is0<2;is0++) {
      const int nnode = get_seglen_vertex(frags,iv0,is0);
      maxcount = ( maxcount > nnode ? maxcount : nnode );
      totcount += nnode;
      assert(nnode >= 0);
    }}
    assert(totcount == nedge);

    bincount = safe_malloc(sizeof(int) * (maxcount+1));

    { int ii; for(ii=0;ii<maxcount+1;ii++) {
      bincount[ii] = 0;
    }}
    // #pragma omp parallel for
    { IntFragment_ID iv0; int is0;
    for(iv0=0;iv0<nfrag;iv0++) for(is0=0;is0<2;is0++) {
      const int nnode = get_seglen_vertex(frags,iv0,is0);
      bincount[nnode] ++;
    }}

#ifdef CHECK97
    totcount = 0;
    { int ii; 
    for(ii=0;(ii<maxcount+1);ii++) {
      totcount += ii*bincount[ii];
    }}
    assert(totcount == nedge);
#endif // CHECK97
    
    { 
      int ii; IntFragment_ID sumcount=0;
      for(ii=0;(ii<maxcount+1);ii++) {
	sumcount += bincount[ii];
	if(sumcount >= nfrag) { break;} 
	// There are 2*nfrag fragment ends.
      }
      (*median_fragment_end_degree) = ii;
    }
    safe_free(bincount);
  }
}
#endif

/////////////////////////////////////////////////////////////////////

typedef struct {
  int min_bin;
  int max_bin;
  int allocated_number_of_bins;
  int * bins;
} QuickHistogram;

static void histogram_clear( QuickHistogram * histogram )
{
  assert(NULL != histogram);
  {
    int ibin;
    int min_bin = MAX(histogram->min_bin,0);
    int max_bin = MIN(histogram->max_bin,histogram->allocated_number_of_bins-1);
    /* Clear the used part of the histogram. */
    for(ibin=min_bin; ibin <= max_bin; ibin++) {
      histogram->bins[ibin] = 0;
    }
  }
  histogram->min_bin = histogram->allocated_number_of_bins;
  histogram->max_bin = 0;
}

static void histogram_initialize
( QuickHistogram * histogram, int allocated_number_of_bins)
{
  assert(NULL != histogram);
  assert(allocated_number_of_bins > 0);
  histogram->allocated_number_of_bins = allocated_number_of_bins;
  histogram->min_bin = 0;
  histogram->max_bin = allocated_number_of_bins - 1;
  histogram->bins = NULL;
  histogram->bins = safe_malloc(sizeof(int) * allocated_number_of_bins);
  histogram_clear( histogram );
}

void transitive_edge_marking
(
 TStateGlobals * gstate, // For time interval check pointing
 THeapGlobals  * heapva, // For time interval check pointing
 Tfragment     * frags,
 Tedge         * edges,
 TIntEdge_ID   * next_edge_obj,
 const int walk_depth,
 const int cutoff_fragment_end_degree,
 const int work_limit_per_candidate_edge,
 const IntFragment_ID iv_start,
 const int analysis_flag,
 const char Output_Graph_Store[]
 ) {
  const IntFragment_ID nfrag=GetNumFragments(frags);

  IntFragment_ID * visited_a = NULL;
  IntFragment_ID * visited_b = NULL;
  
#ifdef WALK_DEPTH_DIAGNOSTICS
  int64 search_depth_histogram[walk_depth];
  // What was the search depths visited?
  int64 search_path_histogram[walk_depth];
  // What was the search path depth when successful?

  int64 successful_search_depth_histogram[walk_depth];
  // What was the search depths visited?
  int64 successful_search_path_histogram[walk_depth];
  // What was the search path depth when successful?

  int64 failed_search_depth_histogram[walk_depth];
  // What was the search depths visited?
  //int64 failed_search_path_histogram[walk_depth];
  // What was the search path depth when successful?
#endif // WALK_DEPTH_DIAGNOSTICS

  int64 successful_searches = 0, failed_searches = 0;

  // Maximum number of edges to explore per candidate edge.
  int64 work_tally_per_candidate_edge_histogram
    [work_limit_per_candidate_edge+1];
  // The last entry represents a failure.
  
  /* Transitive Edge Removal */
  int64 num_of_triangles_visited = 0;
  int64 num_of_quads_visited = 0;
  int64 ntrans_test_fail = 0;
  FILE *ftrans = NULL;
  
  const int report_interval = 60; // one minute
  const int check_point_interval = 4*60*60; // four hours
  time_t tp3,next_report_time,next_check_point_time;

  visited_a = safe_malloc(sizeof(IntFragment_ID) * 2 * nfrag);
  visited_b = safe_malloc(sizeof(IntFragment_ID) * 2 * nfrag);

  // Was this fragment seen from the target overlap edge before?
  
#ifdef DONT_RUN_IN_SYMMETRIC_MODE    
  if(GetNumVA_IntEdge_ID(next_edge_obj) > 0) {
    // We do not know if the graph is in adjacency list representation.
    reorder_edges( frags, edges, next_edge_obj );
  }
#endif // DONT_RUN_IN_SYMMETRIC_MODE    

#if 1
  check_symmetry_of_the_edge_mates( frags, edges, next_edge_obj);
#ifdef DEBUG_RISM
  view_fgb_chkpnt( "RISM_before_trans_before_packing", frags, edges);
#endif // DEBUG_RISM
  {
    time_t tp1 = 0, tp2;
    if(TIMINGS) {
      tp1 = time(NULL);
      fprintf(stderr,"Begin packing the edges before transitive edge marking.\n");
      system_top();
    }
    /* Reduce the amount of memory used for the graph. */
    pack_the_edges( frags, edges, next_edge_obj);
    if(TIMINGS) {
      tp2 = time(NULL); 
      fprintf(stderr,"%10" F_TIME_TP " sec: Finished packing the edges before transitive edge marking.\n",
              (tp2-tp1));
      system_top();
    }
  }
#ifdef DEBUG_RISM
  view_fgb_chkpnt( "RISM_before_trans_after_packing", frags, edges);
#endif // DEBUG_RISM
  check_symmetry_of_the_edge_mates( frags, edges, next_edge_obj);
#endif


  
#if 0 // Time interval triggered diagnostics
  if(NULL != Output_Graph_Store) {
    char fname[CMD_BUFFER_SIZE-1] = {0};
    ftrans = fopen(strcat(strcpy(fname,Output_Graph_Store),
                 "/fgb.trans"),"a");
  }
#endif // Time interval triggered diagnostics
  
  { IntFragment_ID iv0; for(iv0=iv_start;iv0<nfrag;iv0++) {
    int is0; for(is0=0;is0<2;is0++) {
      visited_a[2*iv0+is0] = 2*nfrag;
      visited_b[2*iv0+is0] = 2*nfrag;
    }
  }}

#ifdef WALK_DEPTH_DIAGNOSTICS
  { int64 i; for(i=0;i<walk_depth;i++) {
    successful_search_depth_histogram[i] = 0;
    successful_search_path_histogram[i] = 0;
    failed_search_depth_histogram[i] = 0;
    //failed_search_path_histogram[i] = 0;
  }}
#endif // WALK_DEPTH_DIAGNOSTICS

  { int i; for(i=0;i<work_limit_per_candidate_edge+1;i++) {
    work_tally_per_candidate_edge_histogram[i]=0;
  }}
  
  // Should these histograms be stored in the check-point or remain a
  // batch quanitity?
  
  tp3 = time(NULL);
  next_report_time = tp3;
  next_check_point_time = tp3+check_point_interval;

  fprintf(stderr,"transitively inferable edge marking\n");
  fprintf(stderr,"Cutoff fragment-end degree=%d\n",
	  cutoff_fragment_end_degree);

  if(ftrans != NULL) {
    fprintf(ftrans,"transitively inferable edge marking\n");
    fprintf(ftrans,"Cutoff fragment-end degree=%d\n",
	    cutoff_fragment_end_degree);
    fprintf(ftrans,
	    "trans: seconds, fragments, edges, triangles, quads\n");
  }

  /* Begin: Check each vertex in the fragment overlap graph. */
  {
#ifndef _OPENMP
    IntFragment_ID iv0;
#else // _OPENMP
    int iv0; // OpenMP requires a signed integer type
#pragma omp parallel for schedule(dynamic,10000)
#endif // _OPENMP
    for(iv0=iv_start;iv0<nfrag;iv0++) {
      int is0; for(is0=0;is0<2;is0++) {
        const IntEdge_ID ir0 = get_segstart_vertex(frags,iv0,is0);
        const int nnode  = get_seglen_vertex(frags,iv0,is0);
        int in2;
        
#ifdef _OPENMP
#pragma omp master
#endif // _OPENMP
	if(NULL != ftrans) {
	  time_t tp4=time(NULL);
	  if( (tp4 >= next_report_time) ||
	      (tp4 >= next_check_point_time) ) { 
	    fprintf(ftrans,
		    "trans: %8" F_TIME_TP ", %15" F_IIDP ",%1d, %15" F_IIDP ", %20" F_S64P ", %20" F_S64P "\n", 
		    tp4-tp3,
		    iv0,is0,ir0,
		    num_of_triangles_visited,
		    num_of_quads_visited);
	    fflush(ftrans);
	    next_report_time = tp4 + report_interval;
	  }
	}
      
#ifdef _OPENMP
#pragma omp critical
#endif // _OPENMP
      if(NULL != ftrans) {
        if( time(NULL) >= next_check_point_time ) {
#if 0
          time_t tp1 = 0, tp2;
          if(TIMINGS) {
            tp1 = time(NULL);
            fprintf(stderr,"Begin packing the edges for checkpoint.\n");
            system_top();
          }
          /* Reduce the amount of memory used for the graph. */
          pack_the_edges( frags, edges, next_edge_obj);
          if(TIMINGS) {
            tp2 = time(NULL); 
            fprintf(stderr,"%10" F_TIME_TP " sec: Finished packing the edges\n",
                    (tp2-tp1));
            system_top();
          }
#endif
          if(NULL != ftrans && NULL != Output_Graph_Store) {
            fprintf(ftrans,"trans: %8" F_TIME_TP " starting check point\n", time(NULL)-tp3);
            { 
              {	  
                char thePath1[CMD_BUFFER_SIZE-1]={0};
                char thePath2[CMD_BUFFER_SIZE-1]={0};
                int ierr=0;
                sprintf(thePath1,"%s/%s",Output_Graph_Store,"fgb.ckp_tmp");
                sprintf(thePath2,"%s/%s",Output_Graph_Store,"fgb.ckp_trans");
                write_fgb_store(thePath1, gstate, heapva);
                ierr = accept_tmp_as_final_file( thePath1, thePath2);
                assert(ierr == 0);
              }
              if(analysis_flag) {
                FILE *ffga=NULL;
                char thePath3[CMD_BUFFER_SIZE-1]={0};
                char thePath4[CMD_BUFFER_SIZE-1]={0};
                const int ProcessFragmentAnnotationsForSimulatorCoordinates
                  = (analysis_flag > 1);
                int ierr=0;
                sprintf(thePath3,"%s/%s",Output_Graph_Store,"fga.ckp_tmp");
                sprintf(thePath4,"%s/%s",Output_Graph_Store,"fga.ckp_trans");
                ffga = fopen(thePath3,"w");
                fragment_graph_analysis
                  (/* Input Only */
                   (gstate->max_frag_iid),
                   (heapva->frags),
                   (heapva->edges),
                   (heapva->frag_annotations),
                   ProcessFragmentAnnotationsForSimulatorCoordinates,
                   /* Output only */
                   ffga
                   );
                fclose(ffga);
                ierr = accept_tmp_as_final_file( thePath3, thePath4);
                assert(ierr == 0);
              }
            }
            fprintf(ftrans,"trans: %8" F_TIME_TP " finished check point\n", time(NULL)-tp3);
            fflush(ftrans);
            next_check_point_time = time(NULL) + check_point_interval;
          }
        }
      }

      //if(nnode < cutoff_fragment_end_degree) // Limit the run time.
      for(in2=0;in2<nnode;in2++) { 
	const IntEdge_ID ir2 = ir0+in2;
	const Tnes ir2nes = get_nes_edge(edges,ir2);
	const IntFragment_ID ir2bvx = get_bvx_edge(edges,ir2);
	const int ir2bsx = get_bsx_edge(edges,ir2);
#ifdef SMALLER_ADJACENCY_LIST
	const int nnodeb = get_seglen_vertex(frags,ir2bvx,ir2bsx);
#endif
        
	// num_of_edges_visited ++;

        const int ir2_is_dovetail = is_a_dvt_edge(edges,ir2);
        const int ir2_is_from_contained = is_a_frc_edge(edges,ir2);
        const int ir2_is_to_contained = is_a_toc_edge(edges,ir2);
        const int ir2_is_dgn_contained = is_a_dgn_edge(edges,ir2);

#ifdef SMALLER_ADJACENCY_LIST
	if((nnode > nnodeb) && ( ir2_is_dovetail )) continue;
	// Limit the run time by pursuing the smaller adjacency list
	// for transitively inferable overlap removal.
#endif // SMALLER_ADJACENCY_LIST
	
	if(!(
             ir2_is_dovetail ||
             ir2_is_from_contained ||
             ir2_is_to_contained  ||
             ir2_is_dgn_contained  ||
             (AS_CGB_TOUCHES_CRAPPY_DVT == ir2nes) ||
             (AS_CGB_TOUCHES_CRAPPY_CON == ir2nes) ||
             /* The removed overlaps .... */
             (AS_CGB_REMOVED_BY_TRANSITIVITY_DVT == ir2nes) ||
             (AS_CGB_REMOVED_BY_THRESHOLD_DVT == ir2nes) ||
             (AS_CGB_REMOVED_BY_DUPLICATE_DVT == ir2nes) ||
             (AS_CGB_REMOVED_BY_TRANSITIVITY_CON == ir2nes) ||
             (AS_CGB_REMOVED_BY_THRESHOLD_CON == ir2nes) ||
             (AS_CGB_REMOVED_BY_DUPLICATE_CON == ir2nes) 
           )) {
          const IntFragment_ID ir2avx = get_avx_edge(edges,ir2);
          const int ir2asx = get_asx_edge(edges,ir2);
          const IntFragment_ID ir2afr = get_iid_fragment(frags,ir2avx);
          const IntFragment_ID ir2bfr = get_iid_fragment(frags,ir2bvx);
          fprintf(stderr,"BUUBA " F_IID " %d " F_IID " %d %d\n",
                  ir2afr,ir2asx,ir2bfr,ir2bsx,ir2nes);
        }


        assert(
          ir2_is_dovetail ||
          ir2_is_from_contained ||
          ir2_is_to_contained  ||
          ir2_is_dgn_contained  ||
          (AS_CGB_TOUCHES_CRAPPY_DVT == ir2nes) ||
          (AS_CGB_TOUCHES_CRAPPY_CON == ir2nes) ||
          (AS_CGB_MARKED_BY_BRANCH_DVT == ir2nes) ||
          /* The removed overlaps .... */
          (AS_CGB_REMOVED_BY_TRANSITIVITY_DVT == ir2nes) ||
          (AS_CGB_REMOVED_BY_THRESHOLD_DVT == ir2nes) ||
          (AS_CGB_REMOVED_BY_DUPLICATE_DVT == ir2nes) ||
          (AS_CGB_REMOVED_BY_TRANSITIVITY_CON == ir2nes) ||
          (AS_CGB_REMOVED_BY_THRESHOLD_CON == ir2nes) ||
          (AS_CGB_REMOVED_BY_DUPLICATE_CON == ir2nes)
          );

	if(
           (
            // Do not waste our time on already marked overlaps!
            // But RULE 6 can remove these marked edges too.
            /* The marked overlaps .... */
            (AS_CGB_MARKED_BY_BRANCH_DVT != ir2nes) &&
            (AS_CGB_MARKED_BY_DELETED_DVT != ir2nes) &&
            (AS_CGB_MARKED_BY_DELETED_CON != ir2nes) 
            )
            &&
           /* Do not waste time on removed overlaps .... */
	   (ir2nes != AS_CGB_REMOVED_BY_TRANSITIVITY_DVT) &&
	   (ir2nes != AS_CGB_REMOVED_BY_THRESHOLD_DVT) && 
	   (ir2nes != AS_CGB_REMOVED_BY_DUPLICATE_DVT) &&
	   (ir2nes != AS_CGB_REMOVED_BY_TRANSITIVITY_CON) &&
	   (ir2nes != AS_CGB_REMOVED_BY_THRESHOLD_CON) && 
	   (ir2nes != AS_CGB_REMOVED_BY_DUPLICATE_CON)
	   ) // Filter the candidate edge for removal.
	  {
	    int iremove = FALSE;

	    const IntFragment_ID ir2avx = get_avx_edge(edges,ir2);
	    const int        ir2asx = get_asx_edge(edges,ir2);
	    const IntFragment_ID ir2bvx = get_bvx_edge(edges,ir2);
	    const int        ir2bsx = get_bsx_edge(edges,ir2);
	    const int        ir2ahg = get_ahg_edge(edges,ir2);
	    const int        ir2bhg = get_bhg_edge(edges,ir2);
	    
	    /* Implement Gene^s original acceptable error criterion. */
	    const int alpha      = AS_CGB_TRANSITIVE_SLOP_ALPHA;
	    const int epsilon256 = (int)(256*AS_CGB_TRANSITIVE_SLOP_EPSILON);
	    // 2^8 = 256
	    const int ir2aln = get_length_fragment(frags,ir2avx);
	    const int ir2bln = get_length_fragment(frags,ir2bvx);
	    const int tolerance
	      = alpha + ((epsilon256*(ir2aln-ir2ahg+ir2bln-ir2bhg)) >> 9);
	    // = alpha + epsilon*0.5*(ir2aln-ir2ahg+ir2bln-ir2bhg);
            // For reproducibilty using integer arithematic.
            
	    assert(tolerance >= 0);
	    assert(ir2aln > ir2ahg);
	    assert(ir2bln > ir2bhg);
	    assert(alpha >= 0);
	    assert(epsilon256 >= 0);

	    if(walk_depth == 0) {
              assert(FALSE);
#if 0
	      iremove = mark_edge_if_there_is_a_path
		(ir2,
		 frags,
		 edges,
		 cutoff_fragment_end_degree,
		 tolerance,
		 /* diagnostics */
                 &ntrans_test_fail,
                 &num_of_triangles_visited,
                 &num_of_quads_visited
		 );
#endif              
	    } else {
               int work_tally_per_candidate_edge = 0;
               // Current number of edges explored per candidate edge.
              
#ifdef WALK_DEPTH_DIAGNOSTICS
              { int64 i; for(i=0;i<walk_depth;i++) {
                search_depth_histogram[i] = 0;
                search_path_histogram[i] = 0;
              }}
#endif // WALK_DEPTH_DIAGNOSTICS
              {
                
                const int target_is_dovetail = is_a_dvt_edge(edges,ir2);
                const int target_is_from_contained = is_a_frc_edge(edges,ir2);
                const int target_is_to_contained = is_a_toc_edge(edges,ir2);
                const int target_is_dgn_contained = is_a_dgn_edge(edges,ir2);
                
                const int last_edge_was_containment =
#if 1
                  FALSE  // 
#else                  
                  TRUE   // Allow containment and dovetail overlaps to infer a dovetail
#endif
                  ;
                  

                      iremove = is_there_an_overlap_path
                  ( frags, edges,
                    ir2avx,ir2asx,
#ifdef MATCH_TARGET_EDGE
                    ir2bvx,ir2bsx,ir2ahg,ir2bhg,ir2nes,
                    target_is_dovetail,
                    target_is_from_contained,
                    target_is_to_contained,
                    target_is_dgn_contained,
#endif // MATCH_TARGET_EDGE
                    /* recursion variables: */
                    0, // Use zero at the top level.
                    ir2avx,
                    ir2asx,
                    0, // Use zero at the top level.
                    0, // Use zero at the top level.
                    last_edge_was_containment,
                    // The path is restricted to be (FRC)*(DVT)* .
                    // This logic depends on "last_edge_was_containment"
                    // initialized to TRUE before a inferring path is searched.
                    /* search path limiting: */
                    walk_depth,   // The maximum depth of the stack.
                    tolerance,    // For the overlaps
                    visited_a,
                    visited_b,
                    // Was this fragment seen from the target overlap
                    // edge before?
                    work_limit_per_candidate_edge,
                    &work_tally_per_candidate_edge,
                    /* diagnostic variables: */
#ifdef WALK_DEPTH_DIAGNOSTICS
                    search_depth_histogram,
                    // How many times this depth is visited.
                    search_path_histogram,
                    // How many paths were this length.
#endif // WALK_DEPTH_DIAGNOSTICS
                    &ntrans_test_fail
                    // How many paths failed the last test.
                    // frpt
                    );
              }
              
              work_tally_per_candidate_edge_histogram[work_tally_per_candidate_edge]++;

	    }

	    if(iremove) {
	      switch(ir2nes) {
              case AS_CGB_REMOVED_BY_TRANSITIVITY_CON:
              case AS_CGB_REMOVED_BY_TRANSITIVITY_DVT:
                assert(FALSE);
                break;
              case AS_CGB_CONTAINED_EDGE:
                set_nes_edge(edges,ir2,AS_CGB_REMOVED_BY_TRANSITIVITY_CON);
                break;
	      default:
		// Dovetail overlaps
                set_nes_edge(edges,ir2, AS_CGB_REMOVED_BY_TRANSITIVITY_DVT);
		break;
	      }
	      fix_overlap_edge_mate(frags, edges, ir2);
              successful_searches ++;
#ifdef WALK_DEPTH_DIAGNOSTICS
              { int64 i; for(i=0;i<walk_depth;i++) {
                successful_search_depth_histogram[i] += 
                  search_depth_histogram[i];
                successful_search_path_histogram[i] +=
                  search_path_histogram[i];
              }}
#endif // WALK_DEPTH_DIAGNOSTICS
	    } else {
              failed_searches ++;
#ifdef WALK_DEPTH_DIAGNOSTICS
              { int64 i; for(i=0;i<walk_depth;i++) {
                failed_search_depth_histogram[i] += 
                  search_depth_histogram[i];
                //failed_search_path_histogram[i] += search_path_histogram[i];
              }}
#endif // WALK_DEPTH_DIAGNOSTICS
            }
	  }
      }
    }}
  } /* End: Check each vertex in the fragment overlap graph. */

#ifdef DEBUGGING
  fprintf(stderr,"The transitive edge removal test failed " F_S64
	  " times at the length comparison.\n",ntrans_test_fail);
#endif

  if(NULL != ftrans) {
    time_t tp4=time(NULL);
    int64 successful_search_path_total = 0;
    int64 successful_search_depth_total = 0;
    int64 failed_search_depth_total = 0;
    const IntEdge_ID nedge=GetNumEdges(edges);
    fprintf(ftrans,
	    "trans: %8" F_TIME_TP ",%10" F_IIDP ",%1d,%10" F_IIDP "\n",
	    tp4-tp3,
	    nfrag, 1, nedge);
    fprintf(ftrans,"transx:     %20" F_S64P ", %20" F_S64P "\n", 
            num_of_triangles_visited,
	    num_of_quads_visited);
    fprintf(ftrans,
	    "trans: successful_searches=" F_S64 ",failed_searches=" F_S64 "\n",
              successful_searches,failed_searches);
#ifdef WALK_DEPTH_DIAGNOSTICS
    fprintf(ftrans,
	    "trans: walk_depth, suc.paths, suc.visits, unsuc.visits\n");
    { int i; for(i=0;i<walk_depth;i++) {
      successful_search_path_total += successful_search_path_histogram[i];
      successful_search_depth_total += successful_search_depth_histogram[i];
      failed_search_depth_total += failed_search_depth_histogram[i];
      fprintf(ftrans,
              "%6d %15" F_S64P " %15" F_S64P " %15" F_S64P "\n",
              i,
              successful_search_path_histogram[i],
              successful_search_depth_histogram[i],
              failed_search_depth_histogram[i]
        );
    }}
#endif // WALK_DEPTH_DIAGNOSTICS

    fprintf(ftrans,
	    "trans: walk_tally_per_candidate_edge, counts\n");
    { int i; for(i=0;i<work_limit_per_candidate_edge+1;i++) {
      if(work_tally_per_candidate_edge_histogram[i] > 0) {
        fprintf(ftrans,"%10d %15" F_S64P "\n",
                i, work_tally_per_candidate_edge_histogram[i]);
          }
    }}
    
    // Should these histograms be stored in the check-point or remain a
    // batch quanitity?
    
    fprintf(ftrans,
            "totals:%15" F_S64P " %15" F_S64P " %15" F_S64P "\n",
            successful_search_path_total,
            successful_search_depth_total,
            failed_search_depth_total);
    fflush(ftrans);
  }

#if 1
  check_symmetry_of_the_edge_mates( frags, edges, next_edge_obj);
#ifdef DEBUG_RISM
  view_fgb_chkpnt( "RISM_after_trans_before_packing", frags, edges);
#endif // DEBUG_RISM
  {
    time_t tp1 = 0, tp2;
    if(TIMINGS) {
      tp1 = time(NULL);
      fprintf(stderr,"Begin packing the edges after transitive edge marking.\n");
      system_top();
    }
    /* Reduce the amount of memory used for the graph. */
    pack_the_edges( frags, edges, next_edge_obj);
    if(TIMINGS) {
      tp2 = time(NULL); 
      fprintf(stderr,"%10" F_TIME_TP " sec: Finished packing the edges after transitive edge marking.\n",
              (tp2-tp1));
      system_top();
    }
  }
#ifdef DEBUG_RISM
  view_fgb_chkpnt( "RISM_after_trans_after_packing", frags, edges);
#endif // DEBUG_RISM
  check_symmetry_of_the_edge_mates( frags, edges, next_edge_obj);
#endif
  
  if(NULL != ftrans && NULL != Output_Graph_Store) {
    fprintf(ftrans,"trans: %8" F_TIME_TP " starting check point\n", time(NULL)-tp3);
    { 
      int ierr=0;
      { 
        char thePath1[CMD_BUFFER_SIZE-1]={0};
        char thePath2[CMD_BUFFER_SIZE-1]={0};
        sprintf(thePath1,"%s/%s",Output_Graph_Store,"fgb.ckp_tmp");
        sprintf(thePath2,"%s/%s",Output_Graph_Store,"fgb.ckp_trans");
        write_fgb_store(thePath1, gstate, heapva);
        ierr = accept_tmp_as_final_file( thePath1, thePath2);
        assert(ierr == 0);
      }
      if(analysis_flag) {
        FILE *ffga=NULL;
        char thePath3[CMD_BUFFER_SIZE-1]={0};
        char thePath4[CMD_BUFFER_SIZE-1]={0};
        const int ProcessFragmentAnnotationsForSimulatorCoordinates
          = (analysis_flag > 1);
        sprintf(thePath3,"%s/%s",Output_Graph_Store,"fga.ckp_tmp");
        sprintf(thePath4,"%s/%s",Output_Graph_Store,"fga.ckp_trans");
        ffga = fopen(thePath3,"w");
        fragment_graph_analysis
          (/* Input Only */
           gstate->max_frag_iid,
           heapva->frags,
           heapva->edges,
           heapva->frag_annotations,
           ProcessFragmentAnnotationsForSimulatorCoordinates,
           /* Output only */
           ffga
           );
        fclose(ffga);
        ierr = accept_tmp_as_final_file( thePath3, thePath4);
        assert(ierr == 0);
      }
    }
    fprintf(ftrans,"trans: %8" F_TIME_TP " finished check point\n", time(NULL)-tp3);
    fflush(ftrans);
  }
  if(NULL != ftrans) { fclose(ftrans);}
  safe_free(visited_a);
  safe_free(visited_b);
}
