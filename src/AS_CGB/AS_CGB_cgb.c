
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

static char CM_ID[] = "$Id: AS_CGB_cgb.c,v 1.21 2007-07-27 20:25:45 brianwalenz Exp $";

//  This module builds the chunk graph from the fragment essential
//  overlap graph with contained fragments as an augmentation, and
//  writes the chunk graph using the assembler message library.
// 
//  Author: Clark Mobarry
//
//  The basic chunk data structures are:
// 
//  (1) a segmented array of type "TChunkFrag" called chunkfrags[] used
//  to store the fragments of each chunk,
//
//  (2) a segmented array of type "TChunkOverlap" called chunkedges[]
//  used to store the overlaps between chunks,
//
//  (3) an array of type "TChunkMesg" thechunks[].
//
//  The first three arrays are segmented to store variable length
//  information for the chunks.  Each member of the last array has a
//  segment header index into each of the three previous arrays. 

#include "AS_CGB_all.h"

#define AS_CGB_EDGE_NOT_VISITED     INT32_MAX



//  AS_CGB_chimeras.c
uint32 count_chimeras(const char  chimeras_report_filename[],
                      const float cgb_unique_cutoff,
                      const Tfragment     frags[],
                      const Tedge         edges[],
                      TChunkFrag          chunk_frags[],
                      TChunkMesg          chunks[]);

//  AS_CGB_chimeras.c
uint32 count_crappies(const char  crappies_report_filename[],
                      const float cgb_unique_cutoff,
                      const Tfragment     frags[],
                      const Tedge         edges[],
                      TChunkFrag          chunk_frags[],
                      TChunkMesg          chunks[]);


static
int
comparefloats(const void * const a, const void * const b) {
  return((int)(100000.0 * (*((float *)a) - *((float *)b))));
}


//  Check for overly aggressive transitive overlap trimming. That is,
//  inform when edge trimming removes all the overlaps on a fragment
//  end.
//
static
void
check_edge_trimming(Tfragment *frags, Tedge *edges) {
  IntFragment_ID nfrag = GetNumFragments(frags);
  IntEdge_ID     ie;
  IntFragment_ID ifrag;
  int            is; 

  for(ifrag=0; ifrag<nfrag; ifrag++) {
    Tlab ilab = get_lab_fragment(frags,ifrag);
    for(is=0; is<2; is++) {
      IntEdge_ID segstart = get_segstart_vertex(frags,ifrag,is);
      int seglength = get_seglen_vertex(frags,ifrag,is);
      int new_count_dvt = 0;
      int old_count_dvt = 0;
      int new_count_con = 0;
      int old_count_con = 0;
      for(ie=segstart; ie<segstart+seglength; ie++) {
        Tnes ines = get_nes_edge(edges,ie);
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
	fprintf(stderr,"Dovetail Disconnected (non-contained) fragment end "F_IID" %d lab=%d old_count=%d new_count=%d\n",
		get_iid_fragment(frags,ifrag),
                is, ilab, old_count_dvt, new_count_dvt);
      }
    }
  }
}




static
void
add_fragment_to_chunk(int pass,
                      Tfragment frags[],
                      Tedge edges[],
                      IntChunk_ID ichunk,
                      IntFragment_ID vid,
                      int iforward,
                      int64  sum_of_ahg,
                      int64  sum_of_bhg,
                      int ftic[],
                      TChunkFrag chunkfrags[],
                      IntFragment_ID *nfrag_essential_in_chunk,
                      int64  *nbase_essential_sampled_in_chunk,
                      IntFragment_ID *nfrag_contained_in_chunk,
                      int64  *nbase_contained_sampled_in_chunk,
                      int64  *max_ahg_of_contained,
                      int64  *max_bhg_of_contained) {

  //  This is a recursive routine.

  //  Given a undirected (directed) graph G = (V,E) with n-vertices
  //  and an array PLACEMENT() initially set to a sentinal value, the
  //  algorithm visits all vertices reachable from v using dovetail
  //  and to_contained edges.

  Tlab   ilabel   = get_lab_fragment(frags,vid);
  int64  ioffseta = sum_of_ahg;
  int64  ioffsetb = sum_of_ahg + get_length_fragment(frags,vid);
  int    ivsx;

  //  If, instead, ioffsetb = sum_of_bhg, that would have been the
  //  USE_SUM_OF_BHG_FOR_CONTAINMENTS conditional compilation.

  if ((AS_CGB_UNPLACEDCONT_FRAG == ilabel) ||
      (AS_CGB_SINGLECONT_FRAG   == ilabel) ||
      (AS_CGB_MULTICONT_FRAG    == ilabel)) {
    (*nfrag_contained_in_chunk)++;
    (*nbase_contained_sampled_in_chunk) += get_length_fragment(frags,vid);
  } else {
    (*nfrag_essential_in_chunk)++;
    (*nbase_essential_sampled_in_chunk) += get_length_fragment(frags,vid);
  }

  set_o5p_fragment(frags,vid,(iforward ? ioffseta : ioffsetb));
  set_o3p_fragment(frags,vid,(iforward ? ioffsetb : ioffseta));
  set_cid_fragment(frags,vid,ichunk);

  //  Add this fragment to the chunk, unless:
  //    it is pass 0 and we're contained
  //    it is pass 1 and we're multiply contained
  //
  int  skipadd = FALSE;
  if ((pass == 0) && ((AS_CGB_UNPLACEDCONT_FRAG == ilabel) ||
                      (AS_CGB_SINGLECONT_FRAG   == ilabel) ||
                      (AS_CGB_MULTICONT_FRAG    == ilabel)))
    skipadd = TRUE;
  if ((pass == 1) && (AS_CGB_MULTICONT_FRAG    == ilabel))
    skipadd = TRUE;

  if (!skipadd) {
    AppendAChunkFrag(chunkfrags,&vid);
    ftic[vid]++;
  }

  //  If we're a bad fragment type, don't continue
  //
  if ((AS_CGB_UNLABELED_FRAG        == ilabel) || 
      (AS_CGB_HANGING_CRAPPY_FRAG   == ilabel) ||
      (AS_CGB_ORPHANEDCONT_FRAG     == ilabel) ||
      (AS_CGB_MARKED_BREAKER_FRAG   == ilabel) ||
      (AS_CGB_REMOVED_BREAKER_FRAG  == ilabel) ||
      (AS_CGB_MULTICONT_FRAG        == ilabel) ||
      (AS_CGB_BRANCHMULTICONT_FRAG  == ilabel))
    return;


  for(ivsx=0;ivsx<2;ivsx++) {
    IntFragment_ID iavx=vid;
    IntFragment_ID iasx=ivsx^(!iforward);
    IntEdge_ID irc;
    int64  max_ahg = 0;

    //  Set up to search from the essential vertex iavx/iasx
    IntEdge_ID ir0 = get_segstart_vertex(frags,iavx,iasx);
    IntEdge_ID ir1 = ir0 + get_seglen_vertex(frags,iavx,iasx);
      
    //  Begin examining the edges from this vertex
    for(irc=ir0; irc<ir1; irc++) {

      //  These are set for the coverage discriminator statistics
      int iahg = get_ahg_edge(edges,irc); 
      int ibhg = get_bhg_edge(edges,irc);
	
      // These are set for orienting the fragments to the chunk
      // coordinates.
      IntFragment_ID ibvx = get_bvx_edge(edges,irc);
      int ibsx            = get_bsx_edge(edges,irc);
      Tnes inese          = get_nes_edge(edges,irc);

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
        
      // Begin checking for the allowed overlap types.

      if (((AS_CGB_CONTAINED_EDGE == inese) &&
           ((is_a_toc_edge(edges,irc)) ||
            (is_a_dgn_edge(edges,irc) && (get_avx_edge(edges,irc) < get_bvx_edge(edges,irc))))) &&
          (TRUE == get_con_fragment(frags,ibvx)))  {

        //  Found a containment edge

        int jlabel = get_lab_fragment(frags,ibvx);
          
        if ((AS_CGB_UNPLACEDCONT_FRAG == jlabel) ||
            (AS_CGB_SINGLECONT_FRAG   == jlabel)) {

          //  iforward indicates if the essential fragment read (that
          //  is fragment A) is parallel (TRUE) or antiparallel
          //  (FALSE) to the orientation of the chunk.
          //
          int   jforward = !((iforward ^ iasx) ^ ibsx);
          int64 joffseta = ioffseta + (iforward^iasx ? -ibhg : iahg);
          int64 joffsetb = ioffsetb + (iforward^iasx ? -iahg : ibhg);
          int64 offset_error = AS_CGB_TRANSITIVE_SLOP_ALPHA + (int)(AS_CGB_TRANSITIVE_SLOP_EPSILON * get_length_fragment(frags,ibvx));
          
          // During the first pass discover the number of times a
          // fragment can be contained in different chunks.
          //
          // Remember that the contained fragment must be fit into
          // only one location within the chunk.

          int o5p = abs((jforward ? joffseta : joffsetb) - get_o5p_fragment(frags,ibvx));
          int o3p = abs((jforward ? joffsetb : joffseta) - get_o3p_fragment(frags,ibvx));

          if (!((get_cid_fragment(frags,ibvx) == ichunk) &&
                (o5p <= offset_error) &&
                (o3p <= offset_error))) {

            //  Follow the containment edge

            switch (jlabel) {
              case AS_CGB_UNPLACEDCONT_FRAG:
                set_lab_fragment(frags,ibvx,AS_CGB_SINGLECONT_FRAG);
                set_container_fragment(frags,ibvx, get_iid_fragment(frags,iavx));
                break;
              case AS_CGB_SINGLECONT_FRAG:
                set_lab_fragment(frags,ibvx,AS_CGB_MULTICONT_FRAG);
                set_container_fragment(frags,ibvx, get_iid_fragment(frags,iavx));
                break;
              default:
                //  Do nothing.
                break;
            }

            max_ahg = MAX(max_ahg,iahg);

            add_fragment_to_chunk(pass,
                                  frags,
                                  edges,
                                  ichunk,
                                  ibvx,
                                  jforward,
                                  joffseta,
                                  joffsetb,
                                  ftic,
                                  chunkfrags,
                                  nfrag_essential_in_chunk,
                                  nbase_essential_sampled_in_chunk,
                                  nfrag_contained_in_chunk,
                                  nbase_contained_sampled_in_chunk,
                                  max_ahg_of_contained,
                                  max_bhg_of_contained);
          }  //  follow the edge
        }  //  label is contained
      } //  Found a containment edge
    }

    if(ivsx)
      *max_bhg_of_contained = max_ahg;
    else
      *max_ahg_of_contained = max_ahg;
  }  //  over both ivsx (ends?)
}
  



//  Find the AS_CGB_INTRACHUNK_EDGE edge leaving the
//  AS_CGB_INTERCHUNK_FRAG or AS_CGB_HANGING_CHUNK_FRAG fragment.
//
//  Input:  iavx, iasx, frags, edges
//  Input/Output: ibvx, ibsx is the new b-fragment otherwise unchanged.
//  Return: iedge if found, AS_CGB_EDGE_NOT_VISITED if not found
static
IntEdge_ID
find_intrachunk_edge(Tfragment frags[],
                     Tedge edges[], 
                     IntFragment_ID  iavx, int  iasx,
                     IntFragment_ID *ibvx, int *ibsx) {

  //  Search all vertices adjacent from vertex iavx/iasx for a
  //  AS_CGB_INTRACHUNK_EDGE edge. */

  IntEdge_ID ie0 = get_segstart_vertex(frags,iavx,iasx);
  int        ne0 = ie0 + get_seglen_vertex(frags,iavx,iasx);

  for (; ie0 < ne0; ie0++)
    if (get_nes_edge(edges,ie0) == AS_CGB_INTRACHUNK_EDGE) {
      *ibvx = get_bvx_edge(edges,ie0);
      *ibsx = get_bsx_edge(edges,ie0);
      return(ie0);
    }

  return(AS_CGB_EDGE_NOT_VISITED);
}




static
void
make_a_chunk(const int         pass,
             const IntChunk_ID ichunk,
             Tfragment         frags[],
             Tedge             edges[],
             const IntFragment_ID  vid,
             int        ftic[],
             TChunkFrag chunkfrags[],
             IntFragment_ID *chunk_avx,
             int        *chunk_asx,
             IntFragment_ID *chunk_bvx,
             int        *chunk_bsx,
             IntFragment_ID *nfrag_essential_in_chunk,
             int64      *nbase_essential_in_chunk,
             int64      *nbase_essential_sampled_in_chunk,
             IntFragment_ID *nfrag_contained_in_chunk,
             int64      *nbase_contained_sampled_in_chunk,
             int64      *rho) {

  //  ichunk              is the label to give each fragment in the chunk.
  //  chunk_avx/chunk_asx is the fragment and port at the A-tip of the chunk.
  //  chunk_bvx/chunk_bsx is the fragment and port at the B-tip of the chunk.

  IntFragment_ID iavx = vid;
  int            iasx;
  IntFragment_ID ibvx;
  int            ibsx;
  IntEdge_ID     iedge;     // The intra-chunk edge
  int64          sum_of_ahg=0;
  int64          sum_of_bhg=get_length_fragment(frags,vid);
  Tlab           ilabel = get_lab_fragment(frags,vid);
  int64          max_ahg_of_contained_first = 0;
  int64          max_bhg_of_contained_last = 0;
  int64          max_ahg_of_contained = 0;
  int64          max_bhg_of_contained = 0;
  Tlab           alab, blab;

  *nfrag_essential_in_chunk = 0;
  *nbase_essential_in_chunk = 0;
  *nbase_essential_sampled_in_chunk = 0;
  *nfrag_contained_in_chunk = 0;
  *nbase_contained_sampled_in_chunk = 0;

  if (ftic[vid] != 0)
    fprintf(stderr, "GNAT1: iid="F_IID" vid="F_IID" lab=%d ftic=%d pass=%d\n",
            get_iid_fragment(frags,vid), vid, ilabel, ftic[vid], pass);
  assert(ftic[vid] == 0);

  if((AS_CGB_SOLO_FRAG            == ilabel) ||
     (AS_CGB_HANGING_FRAG         == ilabel) ||
     (AS_CGB_HANGING_CRAPPY_FRAG  == ilabel) ||
     (AS_CGB_THRU_FRAG            == ilabel) ||
     (AS_CGB_BRANCHMULTICONT_FRAG == ilabel) ||
     (AS_CGB_ORPHANEDCONT_FRAG    == ilabel) ||
     ((AS_CGB_SINGLECONT_FRAG     == ilabel) && (pass >= 1)) ||
     ((AS_CGB_MULTICONT_FRAG      == ilabel) && (pass >= 2))) {

    add_fragment_to_chunk(pass,
                          frags, edges, ichunk, 
                          vid,
                          TRUE,
                          sum_of_ahg,
                          sum_of_bhg,
                          ftic,
                          chunkfrags, 
                          nfrag_essential_in_chunk, 
                          nbase_essential_sampled_in_chunk,
                          nfrag_contained_in_chunk, 
                          nbase_contained_sampled_in_chunk,
                          &max_ahg_of_contained_first,
                          &max_bhg_of_contained_last);

    *rho  = ( max_ahg_of_contained_first + max_bhg_of_contained_last )/2;
    *chunk_avx = vid;
    *chunk_asx = FALSE;    //  The 5-prime end of the fragment is the A-end of the chunk
    *chunk_bvx = vid;      //  Both ends of a spanned chunk are from the same fragment
    *chunk_bsx = TRUE;     //  The 3-prime end of the fragment is the B-end of the chunk
    (*nbase_essential_in_chunk) = get_length_fragment(frags,vid);
    return;
  }
  
  //  INTRACHUNK && times==0 --> for circular chunks
  //
  assert((AS_CGB_INTERCHUNK_FRAG    == ilabel) ||
         (AS_CGB_HANGING_CHUNK_FRAG == ilabel) ||
         ((AS_CGB_INTRACHUNK_FRAG   == ilabel) && (ftic[vid] == 0)));

  {
    //  Find the first intra-chunk edge of the chunk
    
    IntFragment_ID iedge_prefix = find_intrachunk_edge(frags,edges,iavx,FALSE,&ibvx,&ibsx);
    IntFragment_ID iedge_suffix = find_intrachunk_edge(frags,edges,iavx,TRUE,&ibvx,&ibsx);
    
    alab = get_lab_fragment(frags,iavx);

    //  If the first fragment is not a AS_CGB_INTRACHUNK_FRAG fragment,
    //  then there must be one and only one intra-chunk edge leaving a
    //  chunk-end
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
    }
    assert((AS_CGB_INTRACHUNK_FRAG == alab) ||
           (((iedge_prefix == AS_CGB_EDGE_NOT_VISITED) && 
             (iedge_suffix != AS_CGB_EDGE_NOT_VISITED)) ||
            ((iedge_prefix != AS_CGB_EDGE_NOT_VISITED) && 
             (iedge_suffix == AS_CGB_EDGE_NOT_VISITED))));

    iasx  = (iedge_prefix == AS_CGB_EDGE_NOT_VISITED);
    iedge = (iasx ? iedge_suffix : iedge_prefix);
    
    iavx  = get_avx_edge(edges,iedge);
    iasx  = get_asx_edge(edges,iedge); /* True if the suffix of the first fragment is in the chunk. */ 
    ibvx  = get_bvx_edge(edges,iedge);
    ibsx  = get_bsx_edge(edges,iedge);

    assert(vid == iavx);

    alab = get_lab_fragment(frags,iavx);
    blab = get_lab_fragment(frags,ibvx);
    if (!(((AS_CGB_INTRACHUNK_FRAG    == blab) ||
           (AS_CGB_INTERCHUNK_FRAG    == blab) ||
           (AS_CGB_HANGING_CHUNK_FRAG == blab))))
      fprintf(stdout,"FIRST INTRACHUNK EDGE iedge,iafr,iavx,iasx,ibfr,ibvx,ibsx,alab,blab "F_IID","F_IID","F_IID",%d,"F_IID","F_IID",%d,%d,%d\n", 
              iedge, get_iid_fragment(frags,iavx), iavx, iasx, get_iid_fragment(frags,ibvx), ibvx, ibsx, alab, blab);
    assert((AS_CGB_INTRACHUNK_FRAG    == blab) ||
           (AS_CGB_INTERCHUNK_FRAG    == blab) ||
           (AS_CGB_HANGING_CHUNK_FRAG == blab) );
  }

  assert(ftic[iavx] == 0);
  assert(ftic[ibvx] == 0);

  assert((ilabel == AS_CGB_INTERCHUNK_FRAG) ||
         (ilabel == AS_CGB_INTRACHUNK_FRAG) ||
         (ilabel == AS_CGB_HANGING_CHUNK_FRAG));

  //  Process the first AS_CGB_INTERCHUNK fragment specially

  *chunk_avx = iavx;
  *chunk_asx = !iasx;

  add_fragment_to_chunk(pass,
                        frags, edges, ichunk, 
                        vid, iasx,
                        sum_of_ahg,
                        sum_of_bhg,
                        ftic,
                        chunkfrags, 
                        nfrag_essential_in_chunk, 
                        nbase_essential_sampled_in_chunk,
                        nfrag_contained_in_chunk, 
                        nbase_contained_sampled_in_chunk,
                        &max_ahg_of_contained_first,
                        &max_bhg_of_contained);

  //  Follow unvisited INTRACHUNK edges to the other end of a chunk
  //  and process vertices along the way.  A circular chunk will stop
  //  at a previously visited INTRACHUNK fragment.

  while((ftic[ibvx] == 0) &&
        (iedge != AS_CGB_EDGE_NOT_VISITED) &&
	(get_nes_edge(edges,iedge) == AS_CGB_INTRACHUNK_EDGE)) {

    // ibvx will always be a valid fragment in the chunk

    // These are set for the coverage discriminator statistics

    sum_of_ahg += get_ahg_edge(edges,iedge);
    sum_of_bhg += get_bhg_edge(edges,iedge);

    switch(get_lab_fragment(frags,ibvx)) {
      case AS_CGB_INTRACHUNK_FRAG:
        //  CMM says "Process all INTRACHUNK frags."
        add_fragment_to_chunk(pass,
                              frags, edges, ichunk, 
                              ibvx, !ibsx,
                              sum_of_ahg,
                              sum_of_bhg,
                              ftic,
                              chunkfrags, 
                              nfrag_essential_in_chunk, 
                              nbase_essential_sampled_in_chunk,
                              nfrag_contained_in_chunk, 
                              nbase_contained_sampled_in_chunk,
                              &max_ahg_of_contained,
                              &max_bhg_of_contained);
        iavx = ibvx;
        iasx = !ibsx;
        iedge = find_intrachunk_edge(frags,edges,iavx,iasx,&ibvx,&ibsx);

        //  Check that the chunk didn't end.
        if (iedge == AS_CGB_EDGE_NOT_VISITED)
          fprintf(stdout,"AN INTRACHUNK EDGE iedge,iafr,iavx,iasx,ibfr,ibvx,ibsx = "F_IID","F_IID","F_IID",%d,"F_IID","F_IID",%d\n", 
                  iedge,get_iid_fragment(frags,iavx),iavx,iasx,get_iid_fragment(frags,ibvx),ibvx,ibsx);
        assert(iedge != AS_CGB_EDGE_NOT_VISITED);
        break;
      case AS_CGB_HANGING_CHUNK_FRAG:
      case AS_CGB_INTERCHUNK_FRAG:
        //  CMM says "This is the last INTERCHUNK frag we're going to process."
        add_fragment_to_chunk(pass,
                              frags, edges, ichunk, 
                              ibvx, !ibsx,
                              sum_of_ahg,
                              sum_of_bhg,
                              ftic,
                              chunkfrags, 
                              nfrag_essential_in_chunk, 
                              nbase_essential_sampled_in_chunk,
                              nfrag_contained_in_chunk, 
                              nbase_contained_sampled_in_chunk,
                              &max_ahg_of_contained,
                              &max_bhg_of_contained_last);
        iavx = ibvx;
        iasx = !ibsx;

        //  Check that the chunk ended.
        iedge = find_intrachunk_edge(frags,edges,iavx,iasx,&ibvx,&ibsx);
        if (iedge != AS_CGB_EDGE_NOT_VISITED)
          fprintf(stdout,"AN INTERCHUNK EDGE iedge,iafr,iavx,iasx,ibfr,ibvx,ibsx = "F_IID","F_IID","F_IID",%d,"F_IID","F_IID",%d\n", 
                  iedge,get_iid_fragment(frags,iavx),iavx,iasx,get_iid_fragment(frags,ibvx),ibvx,ibsx);
        assert(iedge == AS_CGB_EDGE_NOT_VISITED);

        //  Really, just BPW looking for bugs.  If we hit this, code
        //  previous to 2007-07-25 would have kept inserting frags
        //  into this chunk.  The (iedge != NOT_VISITED) in the while
        //  loop was added at that time.
        assert(ftic[ibvx] != 0);
        break;
      default:
        fprintf(stderr,"Clark's missing fragment label case: iedge,iafr,iavx,iasx,ibfr,ibvx,ibsx,nes = "F_IID","F_IID","F_IID",%d,"F_IID","F_IID",%d,%d\n",
                iedge,get_iid_fragment(frags,iavx),iavx,iasx,get_iid_fragment(frags,ibvx),ibvx,ibsx,get_nes_edge(edges,iedge));
        fprintf(stderr,"Label on fragment was: %d\n", get_lab_fragment(frags,ibvx));
        assert(FALSE);
    }
  }

  *chunk_bvx = ibvx;
  *chunk_bsx = !ibsx;

  //  Undefining COORDS_FROM_A_END would have set rho to:
  //    (sum_of_ahg + sum_of_bhg + max_ahg_of_contained + max_bhg_of_contained) / 2
  //
  (*rho)                      = sum_of_ahg + max_ahg_of_contained;
  (*nbase_essential_in_chunk) = (*rho) + get_length_fragment(frags,(*chunk_bvx));
}



static
void
fill_a_chunk_starting_at(const int pass,
                         const IntFragment_ID vid,
                         Tfragment frags[],
                         Tedge edges[],
                         int ftic[],
                         TChunkFrag chunkfrags[],
                         TChunkMesg thechunks[]) {

  // Get the next available chunk id.
  IntChunk_ID ichunk = (IntChunk_ID)GetNumVA_AChunkMesg(thechunks);

  // Get the next available index to put chunk fragments.
  IntFragment_ID irec_start_of_chunk = (IntFragment_ID)GetNumVA_AChunkFrag(chunkfrags);
  Tlab lab = get_lab_fragment(frags,vid);

  IntFragment_ID  chunk_avx,chunk_bvx;
  int             chunk_asx,chunk_bsx;
  IntFragment_ID  nfrag_in_chunk=0;
  IntFragment_ID  nfrag_essential_in_chunk=0;
  IntFragment_ID  nfrag_contained_in_chunk=0;
  int64           nbase_essential_in_chunk=0;
  int64           nbase_essential_sampled_in_chunk=0;
  int64           nbase_sampled_in_chunk=0;
  int64           nbase_contained_sampled_in_chunk=0;
  int64           rho;

  //  If we're unplaced, it's an error.
  //
  if ((AS_CGB_UNPLACEDCONT_FRAG == lab) || (ftic[vid] != 0))
    fprintf(stderr, "PIKACU: pass=%d iid="F_IID" vid="F_IID" cid="F_IID" con=%d lab=%d container="F_IID" ftic=%d\n",
            pass,
            get_iid_fragment(frags,vid),
            vid,
            ichunk,
            get_con_fragment(frags,vid),
            get_lab_fragment(frags,vid),
            get_container_fragment(frags,vid),
            ftic[vid]);
  assert((AS_CGB_UNPLACEDCONT_FRAG != lab) && (ftic[vid] == 0));

  assert((get_con_fragment(frags,vid) == 0) ||
         (AS_CGB_UNPLACEDCONT_FRAG    == lab) ||
         (AS_CGB_SINGLECONT_FRAG      == lab) ||
         (AS_CGB_MULTICONT_FRAG       == lab) ||
         (AS_CGB_BRANCHMULTICONT_FRAG == lab) ||
         (AS_CGB_ORPHANEDCONT_FRAG    == lab));

  make_a_chunk(pass,
               ichunk,
               frags,
               edges,
               vid,
               ftic,
               chunkfrags,
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

  // Build light chunks on pass 0, otherwise, build heavy chunks.
  if(pass == 0) {
    nfrag_in_chunk         = nfrag_essential_in_chunk;
    nbase_sampled_in_chunk = nbase_essential_sampled_in_chunk;
  } else {
    nfrag_in_chunk         = nfrag_essential_in_chunk         + nfrag_contained_in_chunk;
    nbase_sampled_in_chunk = nbase_essential_sampled_in_chunk + nbase_contained_sampled_in_chunk;
  }

  assert(GetNumVA_AChunkFrag(chunkfrags) == irec_start_of_chunk+nfrag_in_chunk);

  assert(nfrag_in_chunk > 0);
    
  AChunkMesg ch = {0};

  ch.iaccession   = ichunk;
  ch.bp_length    = nbase_essential_in_chunk;
  ch.rho          = rho;

  ch.num_frags    = nfrag_in_chunk;
  ch.f_list       = irec_start_of_chunk;

  ch.chunk_avx    = chunk_avx;
  ch.chunk_asx    = chunk_asx;
  ch.chunk_bvx    = chunk_bvx;
  ch.chunk_bsx    = chunk_bsx;

  ch.a_list_raw   = get_segstart_vertex(frags,chunk_avx,chunk_asx);
  ch.a_degree_raw = get_seglen_vertex(frags,chunk_avx,chunk_asx);
  ch.b_list_raw   = get_segstart_vertex(frags,chunk_bvx,chunk_bsx);
  ch.b_degree_raw = get_seglen_vertex(frags,chunk_bvx,chunk_bsx);

  SetVA_AChunkMesg(thechunks, ichunk, &ch);
}




static
void
make_the_chunks(Tfragment frags[],
                Tedge edges[],
                TChunkFrag chunkfrags[],
                TChunkMesg thechunks[]) {

  IntFragment_ID vid;
  int            pass;


  //  Counts the number of times a fragment occurs in the chunks.
  //
  int *ftic = safe_malloc(sizeof(int) * GetNumFragments(frags));


  // pass==0 is to form the light-chunks, that is chunks ignoring
  // all of the contained fragments.
  //
  // pass==1 is to form the heavy-chunks, that is the light chunks
  // with the uniquely contained fragments accreted onto the light
  // chunks.
  //
  for(pass=0; pass<2; pass++) {
    ResetVA_AChunkMesg(chunkfrags);
    ResetVA_AChunkMesg(thechunks);

    // Initialize flags for chunk membership:
    //
    // Assign a sentinal chunk id to designate that the fragment has
    // not been placed in any chunks yet.  The contained fragments
    // will be check to see if they can be placed in multiple chunks.
    //
    // We need to check if a contained fragment is uniquely placed in
    // a light weight chunk by its coordinates as well; we clear the
    // coordinates.
    //
    // Clear the container fragment IID.
    //
    for(vid=0;vid<GetNumFragments(frags);vid++) { 
      ftic[vid] = 0;

      set_cid_fragment(frags,vid,CHUNK_NOT_VISITED);
      set_o3p_fragment(frags,vid,0);
      set_o5p_fragment(frags,vid,0);

      set_container_fragment(frags,vid,0); 

      if (AS_CGB_SINGLECONT_FRAG == get_lab_fragment(frags,vid))
	set_lab_fragment(frags,vid,AS_CGB_UNPLACEDCONT_FRAG);
    }


    // Assign chunks to essential fragments at the ends of chunks.
    //
    for(vid=0;vid<GetNumFragments(frags);vid++) {
      Tlab lab = get_lab_fragment(frags,vid);

      if ((lab == AS_CGB_SOLO_FRAG) ||
          (lab == AS_CGB_HANGING_FRAG) ||
          (lab == AS_CGB_HANGING_CRAPPY_FRAG) ||
          (lab == AS_CGB_BRANCHMULTICONT_FRAG) ||
          (lab == AS_CGB_THRU_FRAG) ||
          ((lab == AS_CGB_INTERCHUNK_FRAG)    && (ftic[vid] == 0)) ||
          ((lab == AS_CGB_HANGING_CHUNK_FRAG) && (ftic[vid] == 0)))
	fill_a_chunk_starting_at(pass, vid,
                                 frags, edges,
                                 ftic,
                                 chunkfrags, thechunks);
    }

    //count_fragment_and_edge_labels( frags, edges, "In a chunking pass at A");


    // Assign chunks to Circular chunks.
    //
    for(vid=0;vid<GetNumFragments(frags);vid++)  {
      if ((ftic[vid] == 0) &&
          (get_lab_fragment(frags,vid) == AS_CGB_INTRACHUNK_FRAG)) {
        fprintf(stderr,"Starting a possibly circular chunk with fragment IID="F_IID"\n", get_iid_fragment(frags,vid));
	fill_a_chunk_starting_at(pass, vid,
                                 frags, edges,
                                 ftic,
                                 chunkfrags, thechunks);
        fprintf(stderr,"Finished a possibly circular chunk with CID="F_IID"\n", get_cid_fragment(frags,vid));
      }
    }

    //count_fragment_and_edge_labels( frags, edges, "In a chunking pass at B");


    //  Check that fragments used more than once are labeled properly.
    //
    for(vid=0;vid<GetNumFragments(frags);vid++) { 
      if ((ftic[vid] > 1)) {
        Tlab lab = get_lab_fragment(frags,vid);
        if ((lab != AS_CGB_UNPLACEDCONT_FRAG) &&
            (lab != AS_CGB_MULTICONT_FRAG))
          fprintf(stderr,"ERROR: mis-labeled MULTICONT fragment iid=" F_IID " lab=%d \n", get_iid_fragment(frags,vid), lab);

	assert((lab == AS_CGB_UNPLACEDCONT_FRAG) || (lab == AS_CGB_MULTICONT_FRAG));
        assert(TRUE == get_con_fragment(frags,vid));

	set_lab_fragment(frags,vid,AS_CGB_MULTICONT_FRAG);
      }
    }

    //view_fgb_chkpnt( frags, edges);
    //count_fragment_and_edge_labels( frags, edges, "after a chunking pass");
  }


  fprintf(stderr,"End of chunk building passes.\n");



  // Process the chaff fragments...  Assign a singleton chunk for each
  // MULTICONT and ORPHAN contained fragment.  In addition if there is
  // a circular chunk, do not crash.
  //
  // The remaining AS_CGB_UNPLACEDCONT_FRAG fragments are placed in a
  // greedy fashion.

  IntFragment_ID nfrag=GetNumFragments(frags);

  for(vid=0;vid<nfrag;vid++)
    if (get_lab_fragment(frags,vid) == AS_CGB_UNPLACEDCONT_FRAG)
      set_lab_fragment(frags,vid,AS_CGB_ORPHANEDCONT_FRAG);
    
  for(vid=0;vid<nfrag;vid++)
    if ((ftic[vid] == 0) && (get_lab_fragment(frags,vid) != AS_CGB_DELETED_FRAG))
      fill_a_chunk_starting_at(2, vid,
                               frags, edges,
                               ftic,
                               chunkfrags, thechunks);
  


  // Quality control!!!  Check that each fragment is refered to once
  // and only once in all of the unitigs.

  for(vid=0;vid<GetNumFragments(frags);vid++)
    ftic[vid] = 0;

  IntChunk_ID nchunks = (IntChunk_ID)GetNumVA_AChunkMesg(thechunks);
  IntChunk_ID ichunk;

  for(ichunk=0;ichunk<nchunks;ichunk++){
    AChunkMesg      *ch = GetAChunkMesg(thechunks,ichunk);
    IntFragment_ID   iv;

    for (iv=0; iv<ch->num_frags; iv++) {
      IntFragment_ID vid = *GetAChunkFrag(chunkfrags, ch->f_list + iv);

      if ((iv == 0) && (get_container_fragment(frags,vid) != 0))
        fprintf(stderr, "The first fragment in a chunk is contained!   iid="F_IID" vid="F_IID" cid="F_IID" type=%d label=%d con=%d container="F_IID"\n",
                get_iid_fragment(frags,vid),
                vid,
                get_cid_fragment(frags,vid),
                get_typ_fragment(frags,vid),
                get_lab_fragment(frags,vid),
                get_con_fragment(frags,vid),
                get_container_fragment(frags,vid));
      assert(!((iv == 0) && (get_container_fragment(frags,vid) != 0)));

      ftic[vid]++;
    }
  }

  for(vid=0;vid<GetNumFragments(frags);vid++) {
    if (AS_CGB_DELETED_FRAG != get_lab_fragment(frags,vid)) {
      if (ftic[vid] == 0)
        fprintf(stderr,"make_the_chunks QC -- iid="F_IID" vid="F_IID" cid="F_IID" lab=%d ORPHANED FRAGMENT\n",
                get_iid_fragment(frags,vid),
                vid,
                get_cid_fragment(frags,vid),
                get_lab_fragment(frags,vid));
      if (ftic[vid] > 1)
        fprintf(stderr,"make_the_chunks QC -- iid="F_IID" vid="F_IID" cid="F_IID" lab=%d ftic=%d\n",
                get_iid_fragment(frags,vid),
                vid,
                get_cid_fragment(frags,vid),
                get_lab_fragment(frags,vid),
                ftic[vid]);

      assert((ftic[vid] == 0) || (ftic[vid] == 1));
    }
  }

  safe_free(ftic);
}



int
count_the_randomly_sampled_fragments_in_a_chunk(Tfragment   frags[],
                                                TChunkFrag  chunkfrags[],
                                                TChunkMesg  thechunks[],
                                                IntChunk_ID chunk_index) { 
  IntFragment_ID   nf = 0;
  AChunkMesg      *ch = GetAChunkMesg(thechunks,chunk_index);
  int              num_frags = ch->num_frags;
  IntFragment_ID   ii = 0;
  
  for(ii=0;ii<num_frags;ii++) {
    FragType type = get_typ_fragment(frags, *GetVA_AChunkFrag(chunkfrags,ch->f_list+ii));

    // Only AS_READ, and AS_EXTR fragments are to be used in Gene
    // Myers coverage discriminator A-statistic.
    //
    if (type == AS_READ || type == AS_EXTR)
      nf++;
  }
  return(nf);
}


float
compute_the_global_fragment_arrival_rate(int           recalibrate,
                                         float         cgb_unique_cutoff,
                                         FILE         *fout,
                                         int64         nbase_in_genome,
                                         Tfragment    *frags,
                                         Tedge        *edges,
                                         float         estimated_global_fragment_arrival_rate,
                                         TChunkFrag   *chunkfrags,
                                         TChunkMesg   *thechunks) { 
  IntChunk_ID ichunk = 0;
  int64  total_rho = 0;
  IntFragment_ID total_nfrags = 0;
  IntFragment_ID total_randomly_sampled_fragments_in_genome = 0;
  float computed_global_fragment_arrival_rate_tmp = 0.f;
  float best_global_fragment_arrival_rate;
  const IntChunk_ID nchunks = (IntChunk_ID)GetNumVA_AChunkMesg(thechunks);
  
  for(ichunk=0;ichunk<nchunks;ichunk++) {
    int64  rho = GetAChunkMesg(thechunks,ichunk)->rho; // The sum of overhangs ...
    int    nf  = count_the_randomly_sampled_fragments_in_a_chunk (frags, chunkfrags, thechunks, ichunk);

    total_rho       += rho;
    total_nfrags    += ( nf > 0 ? nf - 1 : 0 );
    total_randomly_sampled_fragments_in_genome += nf;
  }

  if(NULL != fout) {
    fprintf(fout,"Total rho    = "F_S64"\n", total_rho);
    fprintf(fout,"Total nfrags = "F_IID"\n", total_nfrags);
    computed_global_fragment_arrival_rate_tmp = ( total_rho > 0 ? ((float)total_nfrags)/((float)total_rho) : 0.f );
  
    fprintf(fout,"Estimated genome length = "F_S64"\n", nbase_in_genome);
    fprintf(fout,"Estimated global_fragment_arrival_rate=%f\n", (estimated_global_fragment_arrival_rate));
    fprintf(fout,"Computed global_fragment_arrival_rate =%f\n", computed_global_fragment_arrival_rate_tmp);
    fprintf(fout,"Total number of randomly sampled fragments in genome = " F_IID "\n", total_randomly_sampled_fragments_in_genome);
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

  //  recalibrate: boolean flag to recalibrate global arrival rate to
  //  max unique local arrival rate

  if(recalibrate && (nbase_in_genome == 0)){
    int i;
    float *arrival_rate_array = NULL, *arrival_rate_ptr = NULL;
    size_t arrival_rate_array_size = 0;

    for(ichunk=0;ichunk<nchunks;ichunk++) {
      const int64  rho 
	= GetAChunkMesg(thechunks,ichunk)->rho; // The sum of overhangs ...
      if((int)rho > 10000){
	arrival_rate_array_size += (int)rho / 10000;
      }
    }
    if((float)(arrival_rate_array_size * 20000) > (float)total_rho){// Check there are enough large unitigs
      float min_10_local_arrival_rate = best_global_fragment_arrival_rate;
      float median_local_arrival_rate = best_global_fragment_arrival_rate;
      float max_local_arrival_rate = best_global_fragment_arrival_rate;
      float recalibrated_fragment_arrival_rate = best_global_fragment_arrival_rate;
      size_t num_arrival_rates;
      int median_index;

      arrival_rate_ptr = arrival_rate_array = safe_malloc(sizeof(*arrival_rate_array) * arrival_rate_array_size);
      for(ichunk=0;ichunk<nchunks;ichunk++){
	int64  rho = GetAChunkMesg(thechunks,ichunk)->rho; // The sum of overhangs ...

	if ((int)rho > 10000) {
          int num_10000 = (int)rho / 10000;
          int nf = count_the_randomly_sampled_fragments_in_a_chunk (frags, chunkfrags, thechunks, ichunk);
          float local_arrival_rate = ((float)(nf-1)) / ((float)rho);

	  assert(num_10000 > 0);

	  for(i=0;i<num_10000;i++){
	    assert((size_t)(arrival_rate_ptr - arrival_rate_array) < arrival_rate_array_size);
	    *arrival_rate_ptr++ = local_arrival_rate;
	  }
	}
      }
      num_arrival_rates = (size_t)(arrival_rate_ptr - arrival_rate_array);
      if(num_arrival_rates > 0){
	float tmp_fragment_arrival_rate, max_diff_arrival_rate;
	float prev_arrival_rate, cur_arrival_rate, diff_arrival_rate;
	int max_diff_index;
	qsort(arrival_rate_array, num_arrival_rates, sizeof(*arrival_rate_array), comparefloats);
	min_10_local_arrival_rate = arrival_rate_array[num_arrival_rates / 10];
	median_index = (num_arrival_rates * 5) / 10;
	median_local_arrival_rate = arrival_rate_array[median_index];
	max_local_arrival_rate = arrival_rate_array[num_arrival_rates-1];
	recalibrated_fragment_arrival_rate =
	  arrival_rate_array[(num_arrival_rates * 19) / 20];
	prev_arrival_rate = min_10_local_arrival_rate;
	max_diff_arrival_rate = 0.0;
	for(i=num_arrival_rates / 10;i<median_index;i++){
	  cur_arrival_rate = arrival_rate_array[i];
	  diff_arrival_rate = cur_arrival_rate - prev_arrival_rate;
	  prev_arrival_rate = cur_arrival_rate;
	  if(diff_arrival_rate > max_diff_arrival_rate){
	    max_diff_arrival_rate = diff_arrival_rate;
	  }
	}
	max_diff_arrival_rate *= 2.0;
	max_diff_index = num_arrival_rates - 1;
	for(i=median_index;i<num_arrival_rates;i++){
	  cur_arrival_rate = arrival_rate_array[i];
	  diff_arrival_rate = cur_arrival_rate - prev_arrival_rate;
	  prev_arrival_rate = cur_arrival_rate;
	  if(diff_arrival_rate > max_diff_arrival_rate){
	    max_diff_arrival_rate = diff_arrival_rate;
	    max_diff_index = i - 1;
	    break;
	  }
	}
	max_diff_arrival_rate = arrival_rate_array[max_diff_index];
	if((min_10_local_arrival_rate * 2.0) > (median_local_arrival_rate * 1.25)){
	  tmp_fragment_arrival_rate = median_local_arrival_rate * 1.25;
	}else{
	  tmp_fragment_arrival_rate = min_10_local_arrival_rate * 2.0;
	}
	if(tmp_fragment_arrival_rate < recalibrated_fragment_arrival_rate){
	  recalibrated_fragment_arrival_rate = tmp_fragment_arrival_rate;
	}
	if(max_diff_arrival_rate < recalibrated_fragment_arrival_rate){
	  recalibrated_fragment_arrival_rate = max_diff_arrival_rate;
	}
      }
      if(recalibrated_fragment_arrival_rate > best_global_fragment_arrival_rate){
	best_global_fragment_arrival_rate = recalibrated_fragment_arrival_rate;
#if 0
	if(NULL != fout) {
	  fprintf(fout,"Used recalibrated global_fragment_arrival_rate=%f\n",
		  (best_global_fragment_arrival_rate));
	  fprintf(fout,"Used recalibrated global_fragment_arrival_distance=%f\n",
		  ((best_global_fragment_arrival_rate) > 0.
		   ? 1./(best_global_fragment_arrival_rate)
		   : 0.));
	  fprintf(fout,"Chunk arrival rates sorted at 1/100s\n");
	  for(i=0;i<100;i++){
	    fprintf(fout,"%f\n",arrival_rate_array[((num_arrival_rates * i) / 100)]);
	  }
	  fprintf(fout,"%f\n",max_local_arrival_rate);
	}
#endif
      }
      safe_free(arrival_rate_array);
    }
  }
  return best_global_fragment_arrival_rate;
}



void chunk_graph_build_1(const char * const Graph_Store_File_Prefix,
                         const int walk_depth,
                         const int64  genome_length,
                         const char * chimeras_file,
                         const char * spurs_file,
                         const int recalibrate_global_arrival_rate,
                         const float cgb_unique_cutoff,
                         Tfragment     frags[],
                         Tedge         edges[],
                         float         *global_fragment_arrival_rate,
                         TChunkFrag    *chunkfrags,
                         TChunkMesg    *thechunks) {

  make_the_chunks(frags, edges, chunkfrags, thechunks);

  check_edge_trimming( frags, edges);

  if (genome_length == 0)
    (*global_fragment_arrival_rate) = compute_the_global_fragment_arrival_rate(recalibrate_global_arrival_rate,
                                                                               cgb_unique_cutoff,
                                                                               stderr,
                                                                               genome_length,
                                                                               frags,
                                                                               edges,
                                                                               *global_fragment_arrival_rate,
                                                                               chunkfrags,
                                                                               thechunks);
   
  // Now that the global fragment arrival rate is available, we set
  // the coverage statistic.

  const IntChunk_ID nchunks = (IntChunk_ID)GetNumVA_AChunkMesg(thechunks);
  IntChunk_ID chunk_index;

  for(chunk_index=0; chunk_index < nchunks; chunk_index++) {
    AChunkMesg *ch = GetAChunkMesg( thechunks, chunk_index);

    ch->coverage_stat = compute_coverage_statistic(ch->rho,
                                                   count_the_randomly_sampled_fragments_in_a_chunk(frags,
                                                                                                   chunkfrags,
                                                                                                   thechunks,
                                                                                                   chunk_index),
                                                   (*global_fragment_arrival_rate) );
  }

  if( chimeras_file ) {
    uint32 num_chimeras = 0;
    num_chimeras = count_chimeras(chimeras_file, cgb_unique_cutoff, frags, edges, chunkfrags, thechunks);
    check_edge_trimming( frags, edges);
  }
    
  if( spurs_file ) {
    uint32 num_crappies = 0;
    num_crappies = count_crappies(spurs_file, cgb_unique_cutoff, frags, edges, chunkfrags, thechunks);
    check_edge_trimming( frags, edges);
  }
}
