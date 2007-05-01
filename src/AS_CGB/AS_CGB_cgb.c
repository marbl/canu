
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
= "$Id: AS_CGB_cgb.c,v 1.10 2007-05-01 14:41:43 granger_sutton Exp $";
/* *******************************************************************
 *
 * Module: AS_CGB_cgb.c
 * 
 * Description: This module builds the chunk graph from the fragment
 * essential overlap graph with contained fragments as an
 * augmentation, and writes the chunk graph using the assembler
 * message library.
 *
 * Author: Clark Mobarry
 *********************************************************************/

#include "AS_CGB_all.h"

#define GNAT
#undef DEBUG30
#undef DEBUGGING
#undef DEBUG71

#undef USE_SUM_OF_BHG_FOR_CONTAINMENTS

#define PROCESS_CHAFF

#undef HANGING_BY_REDUCED_EDGES

#undef DEBUG21

#undef DEBUGGING
#ifdef DEBUGGING
#undef DEBUG_SHORTCUT1
#undef DEBUG_SHORTCUT4
#define DEBUG_SHORTCUT1
#define DEBUG_SHORTCUT4
#define DEBUG08
#define DEBUG09
#define DEBUG10
#define DEBUG11
#define DEBUG22
#define DEBUG91
#endif /*DEBUGGING*/

#undef COORDS_FROM_A_END
#define COORDS_FROM_A_END

#define AS_CGB_EDGE_NOT_VISITED     CDS_INT32_MAX


static int comparefloats(const void * const a, const void * const b) 
{
  /* This comparison function is to be used with ANSI qsort() for
     sorting local arrival rates. */
  int icom;
  float aa = *((float*)a);
  float bb = *((float *)b);
  icom = (int)(100000.0 * (aa - bb));
  return icom ;
}

static void check_edge_trimming
( 
 Tfragment frags[],
 Tedge     edges[])
{
  /* Check for overly aggressive transitive overlap trimming. That is,
   * inform when edge trimming removes all the overlaps on a fragment
   * end.  
   */

  const IntFragment_ID nfrag = GetNumFragments(frags);
  IntEdge_ID ie;
  IntFragment_ID ifrag;
  int is; 

  for(ifrag=0; ifrag<nfrag; ifrag++) {
    const Tlab ilab = get_lab_fragment(frags,ifrag);
    for(is=0; is<2; is++) {
      const IntEdge_ID segstart = get_segstart_vertex(frags,ifrag,is);
      const int seglength = get_seglen_vertex(frags,ifrag,is);
      int new_count_dvt = 0;
      int old_count_dvt = 0;
      int new_count_con = 0;
      int old_count_con = 0;
      for(ie=segstart; ie<segstart+seglength; ie++) {
	const Tnes ines = get_nes_edge(edges,ie);
	switch(ines) {
	  /* The containment overlaps... */
	case AS_CGB_CONTAINED_EDGE:
	  new_count_con ++;
	  old_count_con ++;
	  break;
	  /* The dovetail overlaps... */
	case AS_CGB_DOVETAIL_EDGE:
	case AS_CGB_THICKEST_EDGE:
	case AS_CGB_INTERCHUNK_EDGE:
	case AS_CGB_INTRACHUNK_EDGE:
	case AS_CGB_TOUCHES_CONTAINED_EDGE:
	case AS_CGB_BETWEEN_CONTAINED_EDGE:
        case AS_CGB_TOUCHES_CRAPPY_DVT:
	  new_count_dvt ++;
	  old_count_dvt ++;
	  break;
        case AS_CGB_BETWEEN_CRAPPY_DVT:
	case AS_CGB_MARKED_BY_BRANCH_DVT:
	case AS_CGB_REMOVED_BY_TRANSITIVITY_DVT:
	case AS_CGB_REMOVED_BY_THRESHOLD_DVT:
	  old_count_dvt ++;
	  break;
        case AS_CGB_TOUCHES_CRAPPY_CON:
        case AS_CGB_BETWEEN_CRAPPY_CON:
	case AS_CGB_REMOVED_BY_TRANSITIVITY_CON:
	case AS_CGB_REMOVED_BY_THRESHOLD_CON:
          old_count_con ++;
	  break;
	case AS_CGB_REMOVED_BY_DUPLICATE_DVT:
	case AS_CGB_REMOVED_BY_DUPLICATE_CON:
          break;
	default:
          fprintf(stderr,"ie=" F_IID " ines=%d\n", ie, ines);
	  assert(FALSE);
	}
      }

      if( (old_count_dvt > 0) && (new_count_dvt == 0) &&
          (FALSE == get_con_fragment(frags,ifrag)) &&
          (AS_CGB_HANGING_CRAPPY_FRAG != ilab) &&
          (AS_CGB_DELETED_FRAG != ilab) ) {
	fprintf(stderr,"Dovetail Disconnected (non-contained) fragment end"
                " " F_IID " %d lab=%d old_count=%d new_count=%d\n",
		get_iid_fragment(frags,ifrag),
                is, ilab, old_count_dvt, new_count_dvt);
      }
    }
  }
}


/*************************************************************************/

static void add_fragment_to_chunk
( const int pass,
  const int work_limit_placing_contained_fragments,
  int * work_placing_contained_fragments,
  /* input only */
  Tfragment frags[],
  Tedge edges[],
  const IntChunk_ID ichunk,
  const IntFragment_ID vid,
  const int iforward,
  //const Tlab ilabel,
  BPTYPE sum_of_ahg,
  BPTYPE sum_of_bhg,
  /* input/output */
  int fragment_timesinchunks[],
  TChunkFrag chunkfrags[],
  IntFragment_ID *nfrag_essential_in_chunk,
  BPTYPE *nbase_essential_sampled_in_chunk,
  IntFragment_ID *nfrag_contained_in_chunk,
  BPTYPE *nbase_contained_sampled_in_chunk,
  /* output only */
  BPTYPE *max_ahg_of_contained,
  BPTYPE *max_bhg_of_contained
)

{
  /* This is a recursive routine. */
  /* Given a undirected (directed) graph G = (V,E) with n-vertices and
  an array PLACEMENT() initially set to a sentinal value, the algorithm
  visits all vertices reachable from v using dovetail and to_contained
  edges. */

  const Tlab ilabel = get_lab_fragment(frags,vid);
  const int ilen  = get_length_fragment(frags,vid);
  const BPTYPE ioffseta = sum_of_ahg;
#ifdef USE_SUM_OF_BHG_FOR_CONTAINMENTS
  const BPTYPE ioffsetb = sum_of_bhg;
#else 
  const BPTYPE ioffsetb = sum_of_ahg + ilen;
#endif

  {
    if( (AS_CGB_UNPLACEDCONT_FRAG == ilabel) ||
	 (AS_CGB_SINGLECONT_FRAG == ilabel) ||
	 (AS_CGB_MULTICONT_FRAG == ilabel)
	) {
      (*nfrag_contained_in_chunk)++;
      (*nbase_contained_sampled_in_chunk) += ilen;
    } else {
      (*nfrag_essential_in_chunk)++;
      (*nbase_essential_sampled_in_chunk) += ilen;
    }
    {
      //set_lab_fragment(frags,vid,ilabel);
      const BPTYPE offset5p = (iforward ? ioffseta : ioffsetb);
      const BPTYPE offset3p = (iforward ? ioffsetb : ioffseta);
      set_o5p_fragment(frags,vid,offset5p);
      set_o3p_fragment(frags,vid,offset3p);
      set_cid_fragment(frags,vid,ichunk);
    }
  }

  if(!( (pass == 0) && 
	((AS_CGB_UNPLACEDCONT_FRAG == ilabel) ||
	 (AS_CGB_SINGLECONT_FRAG == ilabel) ||
	 (AS_CGB_MULTICONT_FRAG == ilabel)
	 )))
  {
    // Do not place any contained fragments on the first pass.
    AChunkFrag mychunkfrag;
    mychunkfrag.vid  = vid;
    AppendAChunkFrag(chunkfrags,&mychunkfrag);
    fragment_timesinchunks[vid]++;
  }

  if((AS_CGB_UNLABELED_FRAG != ilabel) && 
     (AS_CGB_HANGING_CRAPPY_FRAG != ilabel) &&
     (AS_CGB_ORPHANEDCONT_FRAG != ilabel) &&
     (AS_CGB_MARKED_BREAKER_FRAG != ilabel) &&
     (AS_CGB_REMOVED_BREAKER_FRAG != ilabel) &&
     (AS_CGB_MULTICONT_FRAG != ilabel) &&    // This is a paranoid restriction.
     (AS_CGB_BRANCHMULTICONT_FRAG != ilabel) // This is a paranoid restriction.
     ) {
    // Then do not place any contained fragments into this chunk.
    int ivsx; for(ivsx=0;ivsx<2;ivsx++) {
      const IntFragment_ID iavx=vid;
      const IntFragment_ID iasx=ivsx^(!iforward);
      IntEdge_ID irc;
      BPTYPE max_ahg = 0;
      /* Set up to search from the essential vertex iavx/iasx . */
      const IntEdge_ID ir0 = get_segstart_vertex(frags,iavx,iasx);
      const IntEdge_ID ir1 = ir0 + get_seglen_vertex(frags,iavx,iasx);
      
      for(irc=ir0; irc<ir1; irc++)
        /* Begin examining the edges from this vertex. */ {
	
        /* These are set for the coverage discriminator statistics.*/
        const int iahg = get_ahg_edge(edges,irc); 
        const int ibhg = get_bhg_edge(edges,irc);
	
        // These are set for orienting the fragments to the chunk
        // coordinates.
        const IntFragment_ID ibvx = get_bvx_edge(edges,irc);
        const int ibsx = get_bsx_edge(edges,irc);
        const Tnes inese = get_nes_edge(edges,irc);
        
        
        // The dovetail edges are classified by whether the suffix
        // (3-prime end, graphicly represented by the arrow head) of
        // each fragment is in the overlap region.  Let iasx/ibsx be the
        // flag that says whether the suffix of the A/B fragment is in
        // the overlap. Let iahg/ibhg be the number of base pairs of the
        // A/B fragment not in the overlap region.
        
        // Here iasx==TRUE for both:
        // frag A ------------>       normal-dovetail (ibsx==FALSE), and
        // frag B       ---------->
        // frag A ------------>       innie-dovetail  (ibsx==TRUE).
        // frag B       <----------
        // Here iasx==FALSE for both:
        // frag A <------------        outtie-dovetail (ibsx==FALSE), and
        // frag B       ---------->
        // frag A <------------        antinormal-dovetail (ibsx==TRUE).
        // frag B       <----------
        
        // The overlap detector does not emit antinormal overlaps since
        // they are redundant with normal overlaps if we are allowed to
        // swap the A and B fragments. Thus there are only three
        // independent orientations for dovetail overlaps if we are
        // allowed to swap the A and B fragments.
        
        // There are only two independent orientations for containment overlaps
        // if we are allowed to swap the A and B fragments.
        
        // Here iasx==TRUE for both:
        // frag A ---------------->     normal-contained (ibsx==FALSE), and
        // frag B       ------->...
        // frag A ---------------->     innie-contained  (ibsx==TRUE).
        // frag B       <-------...
        // Here iasx==FALSE for both:
        // frag A <----------------     outtie-contained (ibsx==FALSE), and
        // frag B       ------->...
        // frag A <----------------     antinormal-contained (ibsx==TRUE).
        // frag B       <-------...
        // Those dots (...) are to graphically represent a non-positive ibhg.
        
        if( 
           // Begin checking for the allowed overlap types.
           (
            (AS_CGB_CONTAINED_EDGE == inese)&&
            ((is_a_toc_edge(edges,irc))
             || (is_a_dgn_edge(edges,irc)&&(get_avx_edge(edges,irc) < get_bvx_edge(edges,irc)))
             )
            ) && (TRUE == get_con_fragment(frags,ibvx))
           ) /* Found a containment edge. */ {
          
          const Tlab jlabel = get_lab_fragment(frags,ibvx);
          const int  jlen = get_length_fragment(frags,ibvx);
          BPTYPE joffseta, joffsetb;
          BPTYPE offset_error;
          
          const IntFragment_ID iafr = get_iid_fragment(frags,iavx);
          
          /* iforward indicates if the essential fragment read (that is
             fragment A) is parallel (TRUE) or antiparallel (FALSE) to 
             the orientation of the chunk. */
          
          const int jforward = !((iforward ^ iasx) ^ ibsx);
          /* jforward indicates if the contained fragment is parallel (or
             anti-parallel) to the orientation of the chunk. */
          
          joffseta = ioffseta + (iforward^iasx ? -ibhg : iahg);
          joffsetb = ioffsetb + (iforward^iasx ? -iahg : ibhg);
          
          offset_error = AS_CGB_TRANSITIVE_SLOP_ALPHA 
            + (int)(AS_CGB_TRANSITIVE_SLOP_EPSILON*jlen);
          
          if( 
             // During the first pass discover the number of times a
             // fragment can be contained in different chunks.
             (!((get_cid_fragment(frags,ibvx) == ichunk) &&
               (abs(
                    (jforward ?  joffseta : joffsetb ) // offset5p
                    -get_o5p_fragment(frags,ibvx)) <= offset_error) &&
               (abs( 
                    (jforward ?  joffsetb : joffseta ) // offset3p
                    -get_o3p_fragment(frags,ibvx)) <= offset_error)))
             // Remember that the contained fragment must be fit into
             // only one location within the chunk.
	     &&
	     ((AS_CGB_UNPLACEDCONT_FRAG == jlabel) ||
	      (AS_CGB_SINGLECONT_FRAG == jlabel))
             ) /* Follow the containment edge */ {

	    switch(jlabel) {
	    case AS_CGB_UNPLACEDCONT_FRAG:
	      set_lab_fragment(frags,ibvx,AS_CGB_SINGLECONT_FRAG);
	      set_container_fragment(frags,ibvx,iafr);
	      break;
	    case AS_CGB_SINGLECONT_FRAG:
	      set_lab_fragment(frags,ibvx,AS_CGB_MULTICONT_FRAG);
	      set_container_fragment(frags,ibvx,iafr);
	      //set_container_fragment(frags,ibvx,0);
	      break;
	    default:
	      assert(FALSE);
	      break;
	    }

	    max_ahg = MAX(max_ahg,iahg);
	    
	    add_fragment_to_chunk
	      ( pass,
		work_limit_placing_contained_fragments,
		work_placing_contained_fragments,
		/* input only */
		frags,
		edges,
		ichunk,
		ibvx,
		jforward,
		//jlabel,
		joffseta,
		joffsetb,
		/* input/output */
		fragment_timesinchunks,
		chunkfrags,
		nfrag_essential_in_chunk,
		nbase_essential_sampled_in_chunk,
		nfrag_contained_in_chunk,
		nbase_contained_sampled_in_chunk,
		/* output only */
		max_ahg_of_contained,
		max_bhg_of_contained
		);

	  }  /* Follow the containment edge */
        } // Found a containment edge.
      }
      if(ivsx) {
        *max_bhg_of_contained = max_ahg;
      } else {
        *max_ahg_of_contained = max_ahg;
      }
    }
  }
}
  

static IntEdge_ID find_intrachunk_edge
(/*Input Only*/
 Tfragment frags[],
 Tedge edges[], 
 IntFragment_ID iavx,
 int iasx,
 /*Output Only*/
 IntFragment_ID *ibvx,
 int *ibsx)
{ /* Find the AS_CGB_INTRACHUNK_EDGE edge
     leaving the AS_CGB_INTERCHUNK_FRAG or AS_CGB_HANGING_CHUNK_FRAG fragment. 
     Input:  iavx, iasx, frags, edges
     Input/Output: ibvx, ibsx is the new b-fragment otherwise unchanged.
     Return: iedge if found, AS_CGB_EDGE_NOT_VISITED if not found
  */
  
  /* Search all vertices adjacent from vertex iavx/iasx
   for a AS_CGB_INTRACHUNK_EDGE edge. */
  const IntEdge_ID ie0 = get_segstart_vertex(frags,iavx,iasx);
  const int ne0 = get_seglen_vertex(frags,iavx,iasx);
  IntEdge_ID iedge = AS_CGB_EDGE_NOT_VISITED;
  IntEdge_ID ie1;

  for(ie1=ie0; ie1 < ie0+ne0; ie1++) {
    const IntFragment_ID iv1 = get_bvx_edge(edges,ie1);
    const int is1 = get_bsx_edge(edges,ie1);
    const Tnes ines1 = get_nes_edge(edges,ie1);
    if((ines1 == AS_CGB_INTRACHUNK_EDGE)) {
      iedge = ie1;
      *ibvx = iv1;
      *ibsx = is1;
      break;
    }
  }
  /* Did not find an expected intra-chunk edge if iedge == AS_CGB_EDGE_NOT_VISITED. */
  return iedge;
}


static void make_a_chunk
(/* Input only */
 const int         pass,
 const int         work_limit_placing_contained_fragments,
 const IntChunk_ID ichunk,
 Tfragment         frags[], /* read-only EXCEPT lab and cid are modified */
 Tedge             edges[],
 const IntFragment_ID  vid,
 /* Input and output */
 int        fragment_timesinchunks[],
 TChunkFrag chunkfrags[],
 /* Output only */
 IntFragment_ID *chunk_avx,
 int        *chunk_asx,
 IntFragment_ID *chunk_bvx,
 int        *chunk_bsx,
 IntFragment_ID *nfrag_essential_in_chunk,
 BPTYPE     *nbase_essential_in_chunk,
 BPTYPE     *nbase_essential_sampled_in_chunk,
 IntFragment_ID *nfrag_contained_in_chunk,
 BPTYPE     *nbase_contained_sampled_in_chunk,
 BPTYPE     *rho)
{
  /* 
     ichunk    is the label to give each fragment in the chunk.
     chunk_avx/chunk_asx is the fragment and port at 
               the A-tip of the chunk.
     chunk_bvx/chunk_bsx is the fragment and port at 
               the B-tip of the chunk.

  */
  IntFragment_ID iavx = vid;
  int iasx;
  IntFragment_ID ibvx;
  int ibsx;
  IntEdge_ID iedge; /* The intra-chunk edge. */
  BPTYPE sum_of_ahg=0, sum_of_bhg=get_length_fragment(frags,vid);
  const Tlab ilabel = get_lab_fragment(frags,vid);
  BPTYPE max_ahg_of_contained_first = 0;
  BPTYPE max_bhg_of_contained_last = 0;
  BPTYPE max_ahg_of_contained = 0;
  BPTYPE max_bhg_of_contained = 0;
  int work_placing_contained_fragments = 0;

  *nfrag_essential_in_chunk = 0;
  *nbase_essential_in_chunk = 0;
  *nbase_essential_sampled_in_chunk = 0;
  *nfrag_contained_in_chunk = 0;
  *nbase_contained_sampled_in_chunk = 0;
  
#ifndef GNAT
  assert(fragment_timesinchunks[vid] == 0);
#else // GNAT  
  if(!(fragment_timesinchunks[vid] == 0)) {
    fprintf( stderr,
             "GNAT1: iid=" F_IID " vid=" F_IID " lab=%d fragment_timesinchunks=%d pass=%d \n"
             // "ibfr=%lu ibvx=%lu\n"
             ,
             get_iid_fragment(frags,vid), vid,
             ilabel, fragment_timesinchunks[vid],
             pass
             // ,get_iid_fragment(frags,vid), vid,
             );
  }
#endif // GNAT  
  
  if(
     (AS_CGB_SOLO_FRAG == ilabel) ||
     (AS_CGB_HANGING_FRAG == ilabel) ||
     (AS_CGB_HANGING_CRAPPY_FRAG == ilabel) ||
     (AS_CGB_THRU_FRAG == ilabel) ||
     (AS_CGB_BRANCHMULTICONT_FRAG == ilabel) ||
     (AS_CGB_ORPHANEDCONT_FRAG == ilabel) ||
     ((AS_CGB_MULTICONT_FRAG == ilabel) && (pass > 0))
     ) {
    //const Tlab jlabel=ilabel;
    int iforward=TRUE;
    
    work_placing_contained_fragments = 0;
    add_fragment_to_chunk
      (/* input only */
       pass,
       work_limit_placing_contained_fragments,
       &work_placing_contained_fragments,
       frags, edges, ichunk, 
       vid, iforward,
       //jlabel,
       sum_of_ahg,
       sum_of_bhg,
       /* input/output */
       fragment_timesinchunks,
       chunkfrags, 
       nfrag_essential_in_chunk, 
       nbase_essential_sampled_in_chunk,
       nfrag_contained_in_chunk, 
       nbase_contained_sampled_in_chunk,
       /* output only */
       &max_ahg_of_contained_first,
       &max_bhg_of_contained_last
       );

    *rho  = ( max_ahg_of_contained_first + max_bhg_of_contained_last )/2;
    *chunk_avx = vid;
    *chunk_asx = FALSE;
    /* The 5-prime end of the fragment is the A-end of the chunk. */
    *chunk_bvx = vid; 
    /* Both ends of a spanned chunk are from the same fragment. */
    *chunk_bsx = TRUE; 
    /* The 3-prime end of the fragment is the B-end of the chunk. */
    (*nbase_essential_in_chunk) = get_length_fragment(frags,vid);
    return;
  }
  
  assert((AS_CGB_INTERCHUNK_FRAG == ilabel) ||
	 (AS_CGB_HANGING_CHUNK_FRAG == ilabel)
         || ((AS_CGB_INTRACHUNK_FRAG == ilabel)
	     &&(fragment_timesinchunks[vid] == 0))
         // For circular chunks
         );
  {
    //IntFragment_ID iavx, ibvx; /* The vertices of the intra-chunk edges. */
    //int        iasx, ibsx; /* The vertices of the intra-chunk edges. */
    //IntEdge_ID iedge_prefix,iedge_suffix;

    /* Find the first intra-chunk edge of the chunk. */
    
    /* First try iasx == FALSE, */
    IntFragment_ID iedge_prefix =
      find_intrachunk_edge(frags,edges,iavx,FALSE,&ibvx,&ibsx);

    /* Second try iasx == TRUE, */
    IntFragment_ID iedge_suffix =
      find_intrachunk_edge(frags,edges,iavx,TRUE,&ibvx,&ibsx);
    
    const Tlab alab = get_lab_fragment(frags,iavx);

    /* If the first fragment is not a AS_CGB_INTRACHUNK_FRAG fragment,
       then there must be one and only one intra-chunk edge leaving a
       chunk-end. */
    if(!(
         (AS_CGB_INTRACHUNK_FRAG == alab) ||
         ((((iedge_prefix == AS_CGB_EDGE_NOT_VISITED) &&
            (iedge_suffix != AS_CGB_EDGE_NOT_VISITED)) ||
           ((iedge_prefix != AS_CGB_EDGE_NOT_VISITED) && 
            (iedge_suffix == AS_CGB_EDGE_NOT_VISITED)) )))) {
      fprintf(stderr,"Make A Chunk pass=%d\n",pass);
      fprintf(stderr,"ichunk,iid,iavx,lab=" F_IID "," F_IID "," F_IID ",%d\n",
	      ichunk,get_iid_fragment(frags,iavx),iavx,alab);
      fprintf(stderr,"prefix segstart=" F_IID ",seglen=" F_S32"\n",
	      get_segstart_vertex(frags,iavx,FALSE),
	      get_seglen_vertex(frags,iavx,FALSE));
      fprintf(stderr,"suffix segstart=" F_IID ",seglen=" F_S32 "\n",
	      get_segstart_vertex(frags,iavx,TRUE),
	      get_seglen_vertex(frags,iavx,TRUE));
      fprintf(stderr,"iedge_prefix=" F_IID ", iedge_suffix=" F_IID "\n",
	      iedge_prefix,iedge_suffix);
      // view_fgb_chkpnt("MakeAChunk", frags, edges);
    }


    assert((AS_CGB_INTRACHUNK_FRAG == alab) ||
           (((iedge_prefix == AS_CGB_EDGE_NOT_VISITED) && 
             (iedge_suffix != AS_CGB_EDGE_NOT_VISITED)) ||
            ((iedge_prefix != AS_CGB_EDGE_NOT_VISITED) && 
             (iedge_suffix == AS_CGB_EDGE_NOT_VISITED))));

    iasx  = (iedge_prefix == AS_CGB_EDGE_NOT_VISITED);
    iedge = (iasx ? iedge_suffix : iedge_prefix);
    
    iavx  = get_avx_edge(edges,iedge);
    iasx  = get_asx_edge(edges,iedge); /* True if the suffix of the first 
					  fragment is in the chunk. */ 
    ibvx  = get_bvx_edge(edges,iedge);
    ibsx  = get_bsx_edge(edges,iedge);
    assert(vid == iavx);

    {
      const Tlab alab = get_lab_fragment(frags,iavx);
      const Tlab blab = get_lab_fragment(frags,ibvx);
      if(!(((AS_CGB_INTRACHUNK_FRAG == blab) ||
            (AS_CGB_INTERCHUNK_FRAG == blab) ||
            (AS_CGB_HANGING_CHUNK_FRAG == blab) ))) {

        { 
          IntFragment_ID iafr = get_iid_fragment(frags,iavx);
          IntFragment_ID ibfr = get_iid_fragment(frags,ibvx);
          fprintf(stdout,"FIRST INTRACHUNK EDGE "
                  "iedge, iafr,iavx,iasx, ibfr,ibvx,ibsx ="
                  F_IID ", " F_IID "," F_IID ",%d, " F_IID "," F_IID ",%d\n", 
                  iedge, iafr,iavx,iasx, ibfr,ibvx,ibsx);
        }

        fprintf(stderr,"alab=%d blab=%d\n",
                alab, blab);
      }
      
      assert((AS_CGB_INTRACHUNK_FRAG == blab) ||
             (AS_CGB_INTERCHUNK_FRAG == blab) ||
             (AS_CGB_HANGING_CHUNK_FRAG == blab) );
    }
    
    assert(fragment_timesinchunks[iavx] == 0);
    assert(fragment_timesinchunks[ibvx] == 0);
  }

  {
    const int iforward = iasx;
    assert( (ilabel == AS_CGB_INTERCHUNK_FRAG) ||
            (ilabel == AS_CGB_INTRACHUNK_FRAG) ||
	    // For circular chunks
	    (ilabel == AS_CGB_HANGING_CHUNK_FRAG) );

    /* Process the first AS_CGB_INTERCHUNK fragment specially. */
    
    assert(fragment_timesinchunks[iavx] == 0);
    *chunk_avx = iavx;
    *chunk_asx = ! iasx;

    work_placing_contained_fragments = 0;
    add_fragment_to_chunk
      (	/* input only */
       pass,
       work_limit_placing_contained_fragments,
       &work_placing_contained_fragments,
       frags, edges, ichunk, 
       vid, iforward,
       //ilabel,
       sum_of_ahg,
       sum_of_bhg,
       /* input/output */
       fragment_timesinchunks,
       chunkfrags, 
       nfrag_essential_in_chunk, 
       nbase_essential_sampled_in_chunk,
       nfrag_contained_in_chunk, 
       nbase_contained_sampled_in_chunk,
       /* output only */
       &max_ahg_of_contained_first,
       &max_bhg_of_contained
       );
  }

  /* Follow unvisited INTRACHUNK edges to the other end of a chunk and
     process vertices along the way.  A circular chunk will stop at a
     previously visited INTRACHUNK fragment. */

  while(
	(fragment_timesinchunks[ibvx] == 0) &&
	(get_nes_edge(edges,iedge) == AS_CGB_INTRACHUNK_EDGE)) {
    /* ibvx will always be a valid fragment in the chunk. */

    { /* Process the INTRACHUNK edge */
      /* These are set for the coverage discriminator statistics.*/
      const int iahg = get_ahg_edge(edges,iedge); 
      const int ibhg = get_bhg_edge(edges,iedge);
      assert(fragment_timesinchunks[ibvx] == 0);
      
      sum_of_ahg += iahg;
      sum_of_bhg += ibhg;
    }

    switch(get_lab_fragment(frags,ibvx)) {
    case AS_CGB_INTRACHUNK_FRAG:
      { /* We process the INTRACHUNK fragments here. */
	//const Tlab ilabel = get_lab_fragment(frags,ibvx);

	const int iforward = ! ibsx;

	assert(fragment_timesinchunks[ibvx] == 0);
	
        work_placing_contained_fragments = 0;
	add_fragment_to_chunk
	  ( /* input only */
	   pass,
	   work_limit_placing_contained_fragments,
	   &work_placing_contained_fragments,
	   frags, edges, ichunk, 
	   ibvx, iforward,
           //ilabel,
           sum_of_ahg,
	   sum_of_bhg,
	   /* input/output */
	   fragment_timesinchunks,
	   chunkfrags, 
	   nfrag_essential_in_chunk, 
	   nbase_essential_sampled_in_chunk,
	   nfrag_contained_in_chunk, 
	   nbase_contained_sampled_in_chunk,
	   /* output only */
	   &max_ahg_of_contained,
	   &max_bhg_of_contained
	   );
      }
      /* Find another intra-chunk edge from the other side of the fragment. */
      iavx = ibvx; iasx = ! ibsx;
      iedge = find_intrachunk_edge(frags,edges,iavx,iasx,&ibvx,&ibsx);
#ifdef DEBUG21
      if( iedge == AS_CGB_EDGE_NOT_VISITED ) {
	const IntFragment_ID iafr = get_iid_fragment(frags,iavx);
	const IntFragment_ID ibfr = get_iid_fragment(frags,ibvx);
	fprintf(stdout,"AN INTRACHUNK EDGE "
		"iedge, iafr,iavx,iasx, ibfr,ibvx,ibsx ="
		F_IID ", " F_IID "," F_IID ",%d, " F_IID "," F_IID ",%d\n", 
		iedge, iafr,iavx,iasx, ibfr,ibvx,ibsx);
      }
#endif
      assert(iedge != AS_CGB_EDGE_NOT_VISITED);
      break;
    case AS_CGB_HANGING_CHUNK_FRAG:
    case AS_CGB_INTERCHUNK_FRAG:
      { /* We process the last AS_CGB_INTERCHUNK_FRAG fragment here. */
	//const Tlab ilabel = get_lab_fragment(frags,ibvx);
	const int iforward = ! ibsx;

	assert(fragment_timesinchunks[ibvx] == 0);

        work_placing_contained_fragments = 0;
	add_fragment_to_chunk
	  ( /* input only */
	   pass,
	   work_limit_placing_contained_fragments,
	   &work_placing_contained_fragments,
	   frags, edges, ichunk, 
	   ibvx, iforward,
           //ilabel,
           sum_of_ahg,
	   sum_of_bhg,
	   /* input/output */
	   fragment_timesinchunks,
	   chunkfrags, 
	   nfrag_essential_in_chunk, 
	   nbase_essential_sampled_in_chunk,
	   nfrag_contained_in_chunk, 
	   nbase_contained_sampled_in_chunk,
	   /* output only */
	   &max_ahg_of_contained,
	   &max_bhg_of_contained_last
	   );
      }
      /* End shortcut to the end of a chunk. */
      /* Check again that the chunk ended! */
      iavx = ibvx; iasx = ! ibsx;
      assert( AS_CGB_EDGE_NOT_VISITED == 
	      find_intrachunk_edge(frags,edges,iavx,iasx,&ibvx,&ibsx));
      break;
    default:
      {
	const IntFragment_ID iafr = get_iid_fragment(frags,iavx);
	const IntFragment_ID ibfr = get_iid_fragment(frags,ibvx);
	fprintf(stderr,"Clark's missing fragment label case: "
		"iedge, iafr,iavx,iasx, ibfr,ibvx,ibsx, nes ="
		F_IID ", " F_IID "," F_IID ",%d, " F_IID "," F_IID ",%d, %d\n", 
		iedge, iafr,iavx,iasx, ibfr,ibvx,ibsx,
                get_nes_edge(edges,iedge));
        fprintf(stderr,"Label on fragment was: %d\n", 
		get_lab_fragment(frags,ibvx));
      }
      assert(FALSE);
    }
  }

  *chunk_bvx = ibvx;
  *chunk_bsx = ! ibsx;
#ifndef COORDS_FROM_A_END
  (*rho) = (sum_of_ahg+sum_of_bhg
            +max_ahg_of_contained+max_bhg_of_contained)/2;
  (*nbase_essential_in_chunk) = 
    (*rho) + get_length_fragment(frags,(*chunk_bvx));
#else
  (*rho) = sum_of_ahg+max_ahg_of_contained;
  (*nbase_essential_in_chunk) = 
    (*rho) + get_length_fragment(frags,(*chunk_bvx));
#endif
}

/*************************************************************************/


static void fill_a_chunk_starting_at
  (
   const int pass,
   const int work_limit_placing_contained_fragments,
   const IntFragment_ID vid,
   Tfragment frags[],
   /*Input Only*/
   Tedge edges[],
   /*Input/Output*/
   /*Output Only*/
   int fragment_timesinchunks[],
   TChunkFrag chunkfrags[],
   TChunkMesg thechunks[])
{
  // Get the next available chunk id.
  const IntChunk_ID ichunk 
    = (IntChunk_ID)GetNumVA_AChunkMesg(thechunks);
  // Get the next available index to put chunk fragments.
  const IntFragment_ID irec_start_of_chunk 
    = (IntFragment_ID)GetNumVA_AChunkFrag(chunkfrags);
  const Tlab lab = get_lab_fragment(frags,vid);
  
  IntFragment_ID chunk_avx,chunk_bvx;
  int chunk_asx,chunk_bsx;
  IntFragment_ID 
    nfrag_in_chunk=0,
    nfrag_essential_in_chunk=0,
    nfrag_contained_in_chunk=0;
  BPTYPE nbase_essential_in_chunk=0;
  BPTYPE nbase_essential_sampled_in_chunk=0;
  BPTYPE nbase_sampled_in_chunk=0;
  BPTYPE nbase_contained_sampled_in_chunk=0;
  BPTYPE rho;
  AChunkMesg mychunk;

#if 1
  if(!((AS_CGB_UNPLACEDCONT_FRAG != lab) &&
         (fragment_timesinchunks[vid] == 0)))
  {
    IntFragment_ID iid = get_iid_fragment(frags,vid);
    {
      fprintf(stderr,
              "PIKACU: pass=%d iid=" F_IID " vid=" F_IID " cid=" F_IID " "
              "con=%d lab=%d container=" F_IID " ftic=%d\n",
              pass, iid, vid, ichunk,
              get_con_fragment(frags,vid),
              get_lab_fragment(frags,vid),
              get_container_fragment(frags,vid),
	      fragment_timesinchunks[vid]
              );
    }
  }
#endif

  assert((!get_con_fragment(frags,vid)) ||
         (AS_CGB_UNPLACEDCONT_FRAG == lab) ||
         (AS_CGB_SINGLECONT_FRAG == lab) ||
         (AS_CGB_MULTICONT_FRAG == lab) ||
         (AS_CGB_BRANCHMULTICONT_FRAG == lab) ||
         (AS_CGB_ORPHANEDCONT_FRAG == lab) );

#ifndef GNAT  
  assert(
         (AS_CGB_UNPLACEDCONT_FRAG != lab) &&
         (fragment_timesinchunks[vid] == 0)
         );
#else // GNAT  
  if(!(
       (AS_CGB_UNPLACEDCONT_FRAG != lab) &&
       (fragment_timesinchunks[vid] == 0)
       )) {
    fprintf( stderr,
             "GNAT2: iid=" F_IID " vid=" F_IID " lab=%d fragment_timesinchunks=%d"
             " pass=%d\n",
             get_iid_fragment(frags,vid), vid, lab, fragment_timesinchunks[vid],
             pass );
  }
#endif // GNAT  
         
  make_a_chunk
    (/* Input only */
     pass,
     work_limit_placing_contained_fragments,
     ichunk,
     frags,
     edges,
     vid,
     /* Input and output */
     fragment_timesinchunks,
     chunkfrags,
     /* Output only */
     &chunk_avx,
     &chunk_asx,
     &chunk_bvx,
     &chunk_bsx,
     &nfrag_essential_in_chunk,
     &nbase_essential_in_chunk,
     &nbase_essential_sampled_in_chunk,
     &nfrag_contained_in_chunk,
     &nbase_contained_sampled_in_chunk,
     &rho);

  if(pass == 0) {
    // Build light chunks
    nfrag_in_chunk = nfrag_essential_in_chunk;
    nbase_sampled_in_chunk = nbase_essential_sampled_in_chunk;
  } else {
    // pass > 0
    // Build heavy chunks
    nfrag_in_chunk = nfrag_essential_in_chunk 
      + nfrag_contained_in_chunk;
    nbase_sampled_in_chunk = nbase_essential_sampled_in_chunk
      + nbase_contained_sampled_in_chunk;
  }
  assert(GetNumVA_AChunkFrag(chunkfrags) == 
	 irec_start_of_chunk+nfrag_in_chunk);
    
  mychunk.iaccession   = ichunk;
  mychunk.bp_length    = nbase_essential_in_chunk;
  mychunk.rho          = rho;
  mychunk.a_branch_type = AS_NO_BPOINT;
  mychunk.b_branch_type = AS_NO_BPOINT;
  mychunk.a_branch_point = 0;
  mychunk.b_branch_point = 0;
  mychunk.num_frags    = nfrag_in_chunk;
  mychunk.f_list       = irec_start_of_chunk;
  mychunk.isrc = 0; 

  mychunk.chunk_avx    = chunk_avx;
  mychunk.chunk_asx    = chunk_asx;
  mychunk.chunk_bvx    = chunk_bvx;
  mychunk.chunk_bsx    = chunk_bsx;

  mychunk.a_list_raw   = get_segstart_vertex(frags,chunk_avx,chunk_asx);
  mychunk.a_degree_raw = get_seglen_vertex(frags,chunk_avx,chunk_asx);
  mychunk.b_list_raw   = get_segstart_vertex(frags,chunk_bvx,chunk_bsx);
  mychunk.b_degree_raw = get_seglen_vertex(frags,chunk_bvx,chunk_bsx);

  assert(nfrag_in_chunk>0);

#if 0
  if(nbase_in_genome > 0) {
    // Rho is the number of bases in the chunk for the purposes of the
    // coverage discriminator statistic. It is the sum of the fragment
    // overhangs in the chunk. Thus a singleton chunk has a rho equal
    // to zero.

    // A singleton chunk provides no information as to its local
    // fragment arrival rate. We need at least two closely spaced
    // fragments to get a local estimate of the fragment arrival rate.

    // The arrival rate of fragments in the chunk is:
    // {arrival_rate_local = ((float)(nfrag_in_chunk-1))/(float)rho}.
    // The problem with this formula is that a singleton chunk has a
    // coverage discriminator statistic of 0/0.
	
    // Now formula for the coverage discriminator statistic for the
    // chunk is {(arrival_rate_global/arrival_rate_local - ln(2))
    // *(nfrag_in_chunk-1)}.  The division by zero singularity cancels
    // out to give the formula before a hack for guides:
    // {(arrival_rate_global*rho - ln(2)*(nfrag_in_chunk-1)}.


    // The hack for guides recognizes that guides are not randomly
    // spaced fragments in the genome.  Therefore they can not
    // contribute to the fragment count in the local arrival rate
    // discriminator statistic.

    const float ln_2 = 0.693147; // Natural logarithm of 2.
    const float arrival_rate_global 
      = ((float)(GetNumFragments(frags) - num_of_guides_total))
      /((float)nbase_in_genome);
    const float coverage_statistic 
      = rho*arrival_rate_global
      - ((float)(nfrag_in_chunk-1)-(float)(num_of_guides_in_chunk))
	 *ln_2;
    // Be careful of creating a negative number while using unsigned
    // integer arthimetic.  The behavior is not what we want here.

    /* The coverage discriminator statistic should be positive for
       single coverage, negative for multiple coverage, and near zero
       for indecisive. */
    mychunk.coverage_stat= coverage_statistic;
  } else {
    mychunk.coverage_stat = 0.;
  }
#endif
  
  SetVA_AChunkMesg(thechunks,ichunk,&mychunk);
}




static void make_the_chunks
( /*Input Only*/
 const int work_limit_placing_contained_fragments,
 Tfragment frags[],
 Tedge edges[],
 /*Input/Output*/
 /*Output Only*/
 int fragment_timesinchunks[],
 TChunkFrag chunkfrags[],
 TChunkMesg thechunks[]
 )
{
  const int num_passes=2;
  // To implement the logic to eject multiply containable fragments.

  count_fragment_and_edge_labels( frags, edges,
                                  "before chunking passes");

  { int pass; for(pass=0; pass<num_passes; pass++) {

    fprintf(stderr," before pass=%d\n",pass);
    // pass==0 is to form the light-chunks, that is chunks ignoring
    // all of the contained fragments.

    // pass==1 is to form the heavy-chunks, that is the light chunks
    // with the uniquely contained fragments accreted onto the light
    // chunks.

    ResetVA_AChunkMesg(chunkfrags);
    ResetVA_AChunkMesg(thechunks);
    /* Initialize flags for chunk membership. */

    { IntFragment_ID vid; for(vid=0;vid<GetNumFragments(frags);vid++) { 
      Tlab lab = get_lab_fragment(frags,vid);
      fragment_timesinchunks[vid] = 0;
      set_cid_fragment(frags,vid,CHUNK_NOT_VISITED);
      // Just a sentinal chunk id to designate that a fragment has not
      // been placed in any chunks yet.  The contained fragments will
      // be check to see if they can be placed in multiple chunks.
      set_o3p_fragment(frags,vid,0);
      set_o5p_fragment(frags,vid,0);
      // We need to check if a contained fragment is uniquely placed
      // in a light weight chunk by its coordinates as well.

      if((AS_CGB_UNPLACEDCONT_FRAG == lab) ||
	 (AS_CGB_SINGLECONT_FRAG == lab) ) {
	set_lab_fragment(frags,vid,AS_CGB_UNPLACEDCONT_FRAG);
      }
      set_container_fragment(frags,vid,0); 
      // Zero is the sentinal fragment iid for no container fragment.
    }}

    // Assign chunks to essential fragments at the ends of chunks.
    { IntFragment_ID vid; for(vid=0;vid<GetNumFragments(frags);vid++) 
      /* beginning loop through fragment */ {
      const Tlab ilabv = get_lab_fragment(frags,vid);
      // const IntFragment_ID iafr = get_iid_fragment(frags,vid);
      if(
	 (ilabv == AS_CGB_SOLO_FRAG) ||
	 (ilabv == AS_CGB_HANGING_FRAG) ||
	 (ilabv == AS_CGB_HANGING_CRAPPY_FRAG) ||
         (ilabv == AS_CGB_BRANCHMULTICONT_FRAG) ||
	 (ilabv == AS_CGB_THRU_FRAG) ||
	 ((ilabv == AS_CGB_INTERCHUNK_FRAG)
	  &&(fragment_timesinchunks[vid] == 0)) ||
	 ((ilabv == AS_CGB_HANGING_CHUNK_FRAG)
	  &&(fragment_timesinchunks[vid] == 0)) ||
         ((ilabv == AS_CGB_MULTICONT_FRAG) && (pass > 0))
	 ) {
	fill_a_chunk_starting_at
	  ( pass, 
            work_limit_placing_contained_fragments, vid,
	    frags, edges,
	    fragment_timesinchunks,
	    chunkfrags, thechunks);
      }
    }} /* End of a loop over all fragments. */
    
    count_fragment_and_edge_labels( frags, edges,
                                    "In a chunking pass at A");

    // Assign chunks to Circular chunks.
    { IntFragment_ID vid; for(vid=0;vid<GetNumFragments(frags);vid++) 
      /* beginning loop through fragment */ {
      const Tlab ilabv = get_lab_fragment(frags,vid);
      if( (fragment_timesinchunks[vid] == 0)
	  && (ilabv == AS_CGB_INTRACHUNK_FRAG) ) {
        IntFragment_ID iid0 = get_iid_fragment(frags,vid);
        fprintf(stderr,"Starting a possibly circular chunk with fragment IID="
                F_IID "\n",
                iid0);
	fill_a_chunk_starting_at
	  ( pass, 
            work_limit_placing_contained_fragments, vid,
	    frags, edges,
	    fragment_timesinchunks,
	    chunkfrags, thechunks);
        fprintf(stderr,"Finished a possibly circular chunk with CID="
                F_IID "\n",
                get_cid_fragment(frags,vid));
      }
    }} /* End of a loop over all fragments. */

    count_fragment_and_edge_labels( frags, edges,
                                    "In a chunking pass at B");

  { IntFragment_ID vid; for(vid=0;vid<GetNumFragments(frags);vid++) { 
      const Tlab ilabv = get_lab_fragment(frags,vid);
      if( (fragment_timesinchunks[vid] > 1) ) {

        assert(TRUE == get_con_fragment(frags,vid));
	// We have a multiply contained fragment.

	if(!((ilabv == AS_CGB_UNPLACEDCONT_FRAG) ||
             (ilabv == AS_CGB_MULTICONT_FRAG) )) {
          const IntFragment_ID iid = get_iid_fragment(frags,vid);
          fprintf(stderr,"ERROR: mis-labeled MULTICONT fragment iid=" F_IID " lab=%d \n",
                  iid, ilabv);
        }
	assert( (ilabv == AS_CGB_UNPLACEDCONT_FRAG) ||
	        (ilabv == AS_CGB_MULTICONT_FRAG)
                );

	set_lab_fragment(frags,vid,AS_CGB_MULTICONT_FRAG);

        // This field is used by the consensus phase to walk the
        // fragment overlap layout tree for each unitig.
      }
  }}

  fprintf(stderr," after pass=%d\n",pass);
    // view_fgb_chkpnt( frags, edges);
    count_fragment_and_edge_labels( frags, edges,
                                    "after a chunking pass");
  }} // End of chunk building passes.

  fprintf(stderr,"End of chunk building passes.\n");

#ifdef PROCESS_CHAFF
  {  // Process the chaff fragments...
    
    // Assign a singleton chunk for each MULTICONT and ORPHAN contained
    // fragment.  In addition if there is a circular chunk, do not
    // crash.

    IntFragment_ID vid;
    IntFragment_ID nfrag=GetNumFragments(frags);

    for(vid=0;vid<nfrag;vid++) {
      /* beginning loop through fragment */
      const Tlab lab = get_lab_fragment(frags,vid);
      if(lab == AS_CGB_UNPLACEDCONT_FRAG) {
        set_lab_fragment(frags,vid,AS_CGB_ORPHANEDCONT_FRAG);
      }
    }
    
    for(vid=0;vid<nfrag;vid++) {
      /* beginning loop through fragment */
      const Tlab lab = get_lab_fragment(frags,vid);

      if( (fragment_timesinchunks[vid] == 0) &&
          (lab != AS_CGB_DELETED_FRAG) ) {
        const int pass = num_passes; // CMM This is a HACK!!!!
        // The remaining AS_CGB_UNPLACEDCONT_FRAG fragments are placed in a
        // greedy fashion.
	fill_a_chunk_starting_at
	  ( pass, 
            work_limit_placing_contained_fragments, vid,
	    frags, edges,
	    fragment_timesinchunks,
	    chunkfrags, thechunks);
      }
    }} /* End of a loop over all fragments. */
  fprintf(stderr," after pass=%d\n",num_passes);
#endif // PROCESS_CHAFF
  
  // Quality control!!!
  // Check that each fragment is refered to once and only once in all 
  // of the unitigs.
  fprintf(stderr,"make_the_chunks: Unitig quality control\n");

  { IntFragment_ID vid; for(vid=0;vid<GetNumFragments(frags);vid++) { 
    fragment_timesinchunks[vid] = 0;
  }}

  { 
    const IntChunk_ID nchunks 
      = (IntChunk_ID)GetNumVA_AChunkMesg(thechunks);
    IntChunk_ID ichunk;
    for(ichunk=0;ichunk<nchunks;ichunk++){
      IntFragment_ID iv, irec_start_of_chunk;
      AChunkMesg * const mychunkptr = GetAChunkMesg(thechunks,ichunk);
      assert(mychunkptr != NULL);
      irec_start_of_chunk = mychunkptr->f_list;

      for(iv=0; iv<(mychunkptr->num_frags); iv++) {
	IntFragment_ID ivc,iid,vid;
	AChunkFrag * mychunkfrag = NULL;
	ivc = irec_start_of_chunk+iv;
	mychunkfrag = GetAChunkFrag(chunkfrags,ivc);
	vid = mychunkfrag->vid;
	iid = get_iid_fragment(frags,vid);
	
        if(0 == iv) {
          const int container = get_container_fragment(frags,vid);
          if( container != 0 ) {
            fprintf(stderr,
                    "The first fragment in a chunk with container != 0 "
                    " uid=" F_UID
                    " iid=" F_IID " vid=" F_IID " cid=" F_IID " type=%d label=%d"
                    " con=%d container=" F_IID "\n",
                    get_uid_fragment(frags,vid),
                    get_iid_fragment(frags,vid),
                    vid,
                    get_cid_fragment(frags,vid),
                    get_typ_fragment(frags,vid),
                    get_lab_fragment(frags,vid),
                    get_con_fragment(frags,vid),
                    get_container_fragment(frags,vid));

#if 0
            // Attempt a fix ....
            set_container_fragment(frags,vid,0);
            set_o5p_fragment(frags,vid,0);
            set_o3p_fragment(frags,vid,get_length_fragment(frags,vid));
            set_lab_fragment(frags,vid,AS_CGB_ORPHANEDCONT_FRAG);
#endif
            
          }
        }
        fragment_timesinchunks[vid] ++;
      }
    }
  }

  { IntFragment_ID vid; for(vid=0;vid<GetNumFragments(frags);vid++) { 
    if( (fragment_timesinchunks[vid] != 1) &&
        (AS_CGB_DELETED_FRAG != get_lab_fragment(frags,vid))) {
      Fragment_ID uid = get_uid_fragment(frags,vid);
      IntFragment_ID iid = get_iid_fragment(frags,vid);
      IntChunk_ID cid = get_cid_fragment(frags,vid);
      Tlab lab = get_lab_fragment(frags,vid);
      fprintf(stderr,"make_the_chunks: Quality control\n"
	      "uid=" F_UID " iid=" F_IID " vid=" F_IID " cid=" F_IID " lab=%d "
              "fragment_timesinchunks=%d\n",
	      uid,iid,vid,cid,lab,fragment_timesinchunks[vid]);
    }
#if 0
    assert((fragment_timesinchunks[vid] == 0)|| 
	   (fragment_timesinchunks[vid] == 1));
#endif    
    if( (fragment_timesinchunks[vid] == 0) &&
        (AS_CGB_DELETED_FRAG != get_lab_fragment(frags,vid))) {
      fprintf(stderr,"This fragment is not in any chunk -- Orphaned fragment "
              "uid=" F_UID " iid=" F_IID " vid=" F_IID " cid=" F_IID " type=%d\n", 
	      get_uid_fragment(frags,vid),
	      get_iid_fragment(frags,vid),
	      vid,
	      get_cid_fragment(frags,vid),
	      get_typ_fragment(frags,vid)
	      );
    }
  }}

}

/*************************************************************************/

float compute_the_global_fragment_arrival_rate
( 
 const int           recalibrate, /* Boolean flag to recalibrate global arrival rate to max
				     unique local arrival rate */
 const float         cgb_unique_cutoff, /* threshold for unique chunks */
 FILE               *fout,
 /* Input Only */
 const BPTYPE        nbase_in_genome,
 /* Input/Output */
 const Tfragment     frags[],     /* The internal representation of
				     the fragment reads. I have one
				     extra for the segstart field. */
 const Tedge         edges[],     /* The internal representation of the
				     overlaps. */
 const float         estimated_global_fragment_arrival_rate,
 /* Output Only */
 const TChunkFrag   *chunkfrags,
 const TChunkMesg   *thechunks
 )
{ 
  IntChunk_ID ichunk = 0;
  BPTYPE total_rho = 0;
  IntFragment_ID total_nfrags = 0;
  IntFragment_ID total_randomly_sampled_fragments_in_genome = 0;
  float computed_global_fragment_arrival_rate_tmp = 0.f;
  float best_global_fragment_arrival_rate;
  const IntChunk_ID nchunks = (IntChunk_ID)GetNumVA_AChunkMesg(thechunks);
  
  if(NULL != fout){
    fprintf(fout,"compute_the_global_fragment_arrival_rate(%d,%f,fout,%d,ptr,ptr,%f,ptr,ptr)\n",
	    recalibrate,cgb_unique_cutoff,nbase_in_genome,estimated_global_fragment_arrival_rate);
  }
  for(ichunk=0;ichunk<nchunks;ichunk++) {
    const BPTYPE rho 
      = GetAChunkMesg(thechunks,ichunk)->rho; // The sum of overhangs ...
    const int number_of_randomly_sampled_fragments_in_chunk
      = count_the_randomly_sampled_fragments_in_a_chunk
      ( frags, chunkfrags, thechunks, ichunk);
    
    total_rho += rho;
    total_nfrags += ( number_of_randomly_sampled_fragments_in_chunk > 0 
		      ? number_of_randomly_sampled_fragments_in_chunk - 1 
		      : 0 );
    total_randomly_sampled_fragments_in_genome +=
      number_of_randomly_sampled_fragments_in_chunk;
  }

  if(NULL != fout) {
    fprintf(fout,"Total rho    = " BPFORMAT "\n", total_rho);
    fprintf(fout,"Total nfrags = " F_IID "\n", total_nfrags);
    computed_global_fragment_arrival_rate_tmp
      = ( total_rho > 0 ? ((float)total_nfrags)/((float)total_rho) : 0.f );
  
    fprintf(fout,"Estimated genome length = " BPFORMAT "\n",
            nbase_in_genome);
    fprintf(fout,"Estimated global_fragment_arrival_rate=%f\n",
            (estimated_global_fragment_arrival_rate));
    fprintf(fout,"Computed global_fragment_arrival_rate =%f\n",
            computed_global_fragment_arrival_rate_tmp);
    fprintf(fout,"Total number of randomly sampled fragments in genome = " F_IID "\n",
            total_randomly_sampled_fragments_in_genome);
    fprintf(fout,"Computed genome length  = %f\n",
            ( computed_global_fragment_arrival_rate_tmp > 0.f
              ? (total_randomly_sampled_fragments_in_genome)
              / ( computed_global_fragment_arrival_rate_tmp)
              : 0.f ));
  }

  best_global_fragment_arrival_rate =
    ( nbase_in_genome == 0
      ? computed_global_fragment_arrival_rate_tmp
      : estimated_global_fragment_arrival_rate );

  if(NULL != fout) {
    fprintf(fout,"Used global_fragment_arrival_rate=%f\n",
            (best_global_fragment_arrival_rate));
    fprintf(fout,"Used global_fragment_arrival_distance=%f\n",
            ((best_global_fragment_arrival_rate) > 0.
             ? 1./(best_global_fragment_arrival_rate)
             : 0.));
  }

  if(recalibrate && (nbase_in_genome == 0)){
    int i;
    float min_local_arrival_rate = best_global_fragment_arrival_rate;
    float max_local_arrival_rate = best_global_fragment_arrival_rate;
    float *arrival_rate_array = NULL, *arrival_rate_ptr = NULL;
    size_t num_arrival_rates, arrival_rate_array_size = 0;

    for(ichunk=0;ichunk<nchunks;ichunk++) {
      const BPTYPE rho 
	= GetAChunkMesg(thechunks,ichunk)->rho; // The sum of overhangs ...
      if((int)rho > 10000){
	arrival_rate_array_size += (int)rho / 10000;
      }
    }
    arrival_rate_ptr = arrival_rate_array = safe_malloc(sizeof(*arrival_rate_array) * arrival_rate_array_size);
    for(ichunk=0;ichunk<nchunks;ichunk++){
      const BPTYPE rho
	= GetAChunkMesg(thechunks,ichunk)->rho; // The sum of overhangs ...
      /*      const int number_of_randomly_sampled_fragments_in_chunk
	= count_the_randomly_sampled_fragments_in_a_chunk
	( frags, chunkfrags, thechunks, ichunk);*/

      if((int)rho > 10000)/* &&
	 (compute_coverage_statistic( rho, number_of_randomly_sampled_fragments_in_chunk,
				      best_global_fragment_arrival_rate ) > cgb_unique_cutoff))*/{
	const int num_10000 = (int)rho / 10000;
	const int number_of_randomly_sampled_fragments_in_chunk
	  = count_the_randomly_sampled_fragments_in_a_chunk
	  (frags, chunkfrags, thechunks, ichunk);
	const float local_arrival_rate = ((float)(number_of_randomly_sampled_fragments_in_chunk-1)) /
	  ((float)rho);
	assert(num_10000 > 0);
	for(i=0;i<num_10000;i++){
	  assert((size_t)(arrival_rate_ptr - arrival_rate_array) < arrival_rate_array_size);
	  *arrival_rate_ptr++ = local_arrival_rate;
	}
      }
    }
    num_arrival_rates = (size_t)(arrival_rate_ptr - arrival_rate_array);
    if(num_arrival_rates > 0){
      float recalibrated_fragment_arrival_rate, tmp_fragment_arrival_rate;
      qsort(arrival_rate_array, num_arrival_rates,
	    sizeof(*arrival_rate_array),
	    &comparefloats);
      min_local_arrival_rate = arrival_rate_array[0];
      max_local_arrival_rate = arrival_rate_array[num_arrival_rates-1];
      recalibrated_fragment_arrival_rate =
	arrival_rate_array[((num_arrival_rates * 19) / 20)];
      if((min_local_arrival_rate * 2.0) > (best_global_fragment_arrival_rate * 1.25)){
	tmp_fragment_arrival_rate = best_global_fragment_arrival_rate * 1.25;
      }else{
	tmp_fragment_arrival_rate = min_local_arrival_rate * 2.0;
      }
      if(tmp_fragment_arrival_rate > recalibrated_fragment_arrival_rate){
	best_global_fragment_arrival_rate = recalibrated_fragment_arrival_rate;
      }else{
	best_global_fragment_arrival_rate = tmp_fragment_arrival_rate;
      }
    }
    if(NULL != fout) {
      fprintf(fout,"Used recalibrated global_fragment_arrival_rate=%f\n",
	      (best_global_fragment_arrival_rate));
      fprintf(fout,"Used recalibrated global_fragment_arrival_distance=%f\n",
	      ((best_global_fragment_arrival_rate) > 0.
	       ? 1./(best_global_fragment_arrival_rate)
	       : 0.));
      fprintf(fout,"Chunk arrival rates sorted at 1/10s\n");
      for(i=0;i<10;i++){
	fprintf(fout,"%f\n",arrival_rate_array[((num_arrival_rates * i) / 10)]);
      }
      fprintf(fout,"%f\n",max_local_arrival_rate);
    }
    safe_free(arrival_rate_array);
  }
  return best_global_fragment_arrival_rate;
}



void chunk_graph_build_1
(
 /* Input Only */
 const char * const Graph_Store_File_Prefix,
 const int work_limit_placing_contained_fragments,		  
 const int walk_depth,
 const IntFragment_ID max_frag_iid,
 const BPTYPE nbase_in_genome,
 const char * chimeras_file,
 const char * spurs_file,
 const int recalibrate_global_arrival_rate,
 const float cgb_unique_cutoff,
 /* Input/Output */
 Tfragment     frags[],     /* The internal representation of
			       the fragment reads. I have one
			       extra for the segstart field. */
 Tedge         edges[],     /* The internal representation of the
			       overlaps. */
 /* Output Only */
 float         *global_fragment_arrival_rate,
 TChunkFrag    *chunkfrags,
 TChunkMesg    *thechunks
 ) 
{
  IntFragment_ID nfrag = GetNumFragments(frags);
  // IntEdge_ID nedge = GetNumEdges(edges);

  int *fragment_timesinchunks = NULL;

  fragment_timesinchunks = safe_malloc(sizeof(int) * nfrag);
  
  assert(chunkfrags != NULL);
  assert(thechunks  != NULL);
  assert(fragment_timesinchunks != NULL);


#ifdef DEBUG71
  view_fgb_chkpnt("AAA", frags, edges);
#endif  

  make_the_chunks(/*Input Only*/
		  work_limit_placing_contained_fragments,		  
		  frags,
		  edges,
		  /*Input/Output*/
		  /*Output Only*/
		  fragment_timesinchunks,
		  chunkfrags,
		  thechunks
		  );
  
  check_edge_trimming( frags, edges);

  if( nbase_in_genome == 0) {
    (*global_fragment_arrival_rate) 
      = compute_the_global_fragment_arrival_rate
      ( recalibrate_global_arrival_rate, cgb_unique_cutoff, stderr, nbase_in_genome, frags, edges, *global_fragment_arrival_rate,
        chunkfrags, thechunks );
  }

  AssertPtr(frags);
  AssertPtr(edges);
  AssertPtr(thechunks);
  AssertPtr(chunkfrags);
    
  {
    // Now that the global fragment arrival rate is available, we set
    // the coverage statistic.
    const IntChunk_ID nchunks = (IntChunk_ID)GetNumVA_AChunkMesg(thechunks);
    IntChunk_ID chunk_index;

    for(chunk_index=0; chunk_index < nchunks; chunk_index++) {
      AChunkMesg * const mychunkptr = GetAChunkMesg( thechunks, chunk_index);
      const int number_of_randomly_sampled_fragments_in_chunk
        = count_the_randomly_sampled_fragments_in_a_chunk
        ( frags, chunkfrags, thechunks, chunk_index);
      const BPTYPE rho = mychunkptr->rho;
      const float coverage_stat = compute_coverage_statistic
        ( rho,
          number_of_randomly_sampled_fragments_in_chunk,
          (*global_fragment_arrival_rate) );
      mychunkptr->coverage_stat = coverage_stat;
      // SetVA_AChunkMesg(thechunks,chunk_index,mychunkptr);
    }
  }

  if( chimeras_file ) {
    cds_uint32 num_chimeras = 0;
    num_chimeras = count_chimeras
      (chimeras_file, cgb_unique_cutoff, frags, edges, chunkfrags, thechunks);
    check_edge_trimming( frags, edges);
  }
    
  if( spurs_file ) {
    cds_uint32 num_crappies = 0;
    num_crappies = count_crappies
      (spurs_file, cgb_unique_cutoff, frags, edges, chunkfrags, thechunks);
    check_edge_trimming( frags, edges);
  }
    
  safe_free(fragment_timesinchunks);
}

void chunk_graph_build_2
(
 /* Input Only */
 const int use_consensus,
 const int dont_find_branch_points,
 const float cgb_unique_cutoff,
 const float global_fragment_arrival_rate,
 GateKeeperStore *gkpStore,
 /* Input/Output */
 Tfragment     frags[],     /* The internal representation of
			       the fragment reads. I have one
			       extra for the segstart field. */
 Tedge         edges[],     /* The internal representation of the
			       overlaps. */
 /* Output Only */
 TChunkFrag    *chunkfrags,
 TChunkMesg    *thechunks,
 VA_TYPE(char) *chunkseqs,
 VA_TYPE(char) *chunkquas,
 VA_TYPE(char) *chunksrcs,
 FILE          *fbpts1,
 FILE          *fbpts2
 ) 
{

  /* In case the analysis flag is not set, just make sure that 
     a few variable are initialized so the output routine
     does not bomb. */
  {
    /* All the isrc fields of the chunks are initialized to zero. So
       the simulator annotation at location zero must be
       initialized. */
    EnableRangeVA_char(chunksrcs,1);
    *(Getchar(chunksrcs,0)) = '\0'; 
    // Set the default chunk annotations to the empty string.
  }

#ifdef BRANCHPOINTS
#error branchpoints defined
  if( ! dont_find_branch_points ) {
    find_the_branch_points
      ( use_consensus, cgb_unique_cutoff, global_fragment_arrival_rate,
        gkpStore, frags, edges, chunkfrags, 
	chunkseqs, thechunks, fbpts1);

#ifdef STORE_BRANCH_POINTS_AT_FRAGMENT
    {
      FILE * fout = fbpts2;
      if(NULL != fout) {
        const IntFragment_ID nfrag = GetNumFragments(frags);
        IntFragment_ID vid;
        fprintf(fout, "# fragment-iid fragment-suffix branch-point-type branch-point-offset\n");
        for( vid=0; vid<nfrag; vid++) {
          const IntFragment_ID iid = get_iid_fragment(frags,vid);
          int vsx;
          for( vsx=0; vsx<2; vsx++) {
            int bpt = get_bpt_vertex(frags,vid,vsx);
            if( bpt >0 ) {
              fprintf(fout, F_IID " %d %d %d\n", iid, vsx, AS_INTO_UNIQUE, bpt);
            }
          }
        }
      }
    }
#endif // STORE_BRANCH_POINTS_AT_FRAGMENT
    
#ifndef STORE_BRANCH_POINTS_AT_FRAGMENT
    chunk_end_edge_trimmer
      ( /* input only */
       cgb_unique_cutoff,
       global_fragment_arrival_rate,
       frags,
       chunkfrags,
       /* modify */
       edges, thechunks);
#else // STORE_BRANCH_POINTS_AT_FRAGMENT
    fragment_end_edge_trimmer
      (
       /* input only */
       frags,
       /* modify */
       edges);
#endif // STORE_BRANCH_POINTS_AT_FRAGMENT
    check_edge_trimming( frags, edges);

#ifdef UNDIRECTED_DEBUG_2
    check_edges1(frags, edges, FALSE);
    check_edges2(frags, edges, FALSE);
#endif /*UNDIRECTED_DEBUG_2*/
  } // ! dont_find_branch_points
#endif // BRANCHPOINTS

  /* END CHUNK GRAPH BUILDING */
}

////////////////////////////////////////////////////////////////

#ifdef DEBUG30
    {
      IntFragment_ID iid = get_iid_fragment(frags,vid);
#if 1
  const BPTYPE offset5p = (iforward ? ioffseta : ioffsetb);
  const BPTYPE offset3p = (iforward ? ioffsetb : ioffseta);
#endif
      if( iid == 9334 ) {
	fprintf(stderr,
		"DEBUG xx iid=" F_IID " ichunk=" F_IID " " 
		" ftic=%d\n"
		"sumahg=" BPFORMAT " sumbhg=" BPFORMAT "\n"
		"offset5p=" BPFORMAT " offset3p=" BPFORMAT "\n"
		,
		iid, ichunk, 
		fragment_timesinchunks[vid],
		sum_of_ahg, sum_of_bhg,
		offset5p, offset3p
		);
      }
    }
#endif // DEBUG30

#ifdef DEBUG30
          {
            IntFragment_ID aid = get_iid_fragment(frags,iavx);
            IntFragment_ID bid = get_iid_fragment(frags,ibvx);
            int asx = get_asx_edge(edges,irc);
            int bsx = get_bsx_edge(edges,irc);
            if( aid == 9334 || bid == 9334 ) {
              fprintf(stderr,
                      "DEBUG xa bid=" F_IID " ichunk=" F_IID " " 
                      " %d " BPFORMAT " " BPFORMAT 
                      " ftic=%d asx=%d bsx=%d ahg=%d bhg=%d\n"
                      " aid=" F_IID " acid=" F_IID 
                      " ao5p=" BPFORMAT " ao3p=" BPFORMAT "\n"
                      " bid=" F_IID " bcid=" F_IID
                      " bo5p=" BPFORMAT " bo3p=" BPFORMAT "\n"
                      " ofe=" BPFORMAT "\n",
                      bid, ichunk,
                      jforward, joffseta, joffsetb,
                      fragment_timesinchunks[ibvx],
                      asx, bsx, iahg, ibhg,
                      aid, get_cid_fragment(frags,iavx),
                      get_o5p_fragment(frags,iavx),
                      get_o3p_fragment(frags,iavx),
                      bid, get_cid_fragment(frags,ibvx),
                      get_o5p_fragment(frags,ibvx),
                      get_o3p_fragment(frags,ibvx),
                      offset_error);
            }
          }
#endif
          
#ifdef DEBUG30
          printf("Contained iafr,iavx,ibfr,ibvx,iforward,jbackward,\n"
                 "ioffseta,ioffsetb,iahg,ibhg,joffseta,joffsetb=\n");
          printf("%5" F_IIDP " %5" F_IIDP " %5" F_IIDP " %5" F_IIDP " %2d %2d \n " BPFORMAT5 " " BPFORMAT5 " %5d %5d " BPFORMAT5 " " BPFORMAT5 "\n",
                 iafr,iavx,ibfr,ibvx,iforward,jforward,
                 ioffseta,ioffsetb,iahg,ibhg,joffseta,joffsetb);
#endif
