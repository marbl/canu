
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
= "$Id: AS_FGB_io.c,v 1.20 2007-05-14 13:40:55 brianwalenz Exp $";
/* *******************************************************************
 *
 * Module: AS_FGB_io.c
 * 
 * Description: Fragment Overlap Graph Builder file input and output.
 * This functional unit reads a *.ovl prototype i/o file
 * an massages it for the data structures in the chunk graph builder.
 *
 * Assumptions: 
 * Author: Clark M. Mobarry
 *********************************************************************/

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "AS_global.h"
#include "AS_UTL_version.h"
#include "AS_MSG_pmesg.h"
#include "AS_CGB_all.h"
#include "AS_CGB_util.h"
#include "AS_OVS_overlapStore.h"


int REAPER_VALIDATION = FALSE;

#undef REPORT_THIS_BUG


//  Define this to silently ignore overlaps with no valid fragment.
//  e.g., if you're throwing out fragments after overlaps have been
//  computed, and don't remove those overlaps.  You'll probably have
//  trouble in CGW, so good luck!
//
#undef USE_FRAGMENT_SUBSET

//  Define to print a message whenever we read an overlap
//  that refers to a non-existent frag.
//
#undef REPORT_DEGENERATE_OVERLAPS

//  Define to print a message whenever we delete an
//  undefined fragment
//
#undef REPORT_DELETED_UNDEF_FRAGS



static void keep_a_list_of_extreme_elements
( 
 void const * const new_element, // A new edge record
 IntEdge_ID const * const next, // beginning of the array of next_edge pointers.
 size_t const start_idx,
 void * const base, // array of previously stored edge records.
 size_t const size, // size of an edge record
 int (*compar)(const void *, const void *),
 const int direction
 )
{
  // Using UNIX qsort conventions for the comparison function.
  // The last element in each list must have (ip == next(ip)).
  // The last element in the list is the most extreme element.

  size_t ip = start_idx; 
  // The index of the first element saved in the list.

  assert((direction == 1) || // save the maximum elements
	 (direction == -1)); // save the minimum elements

  if(!(direction*compar(new_element,((char*)base+ip*size)) >= 0)) {
    // The equals case is to make sure that new records appended on to
    // the list are processed since they may be there simply because
    // the list was not full yet.
    return;
  }

  while( ip != next[ip] ) {
    // While not the last element in the list.
    size_t ipp = next[ip]; 
    // ipp = ip +/- 1 in an array
    if(direction*compar(new_element,((char*)base+ipp*size)) > 0) {
      memcpy(((char*)base+ip*size),((char*)base+ipp*size), size);
      ip = ipp;
    } else {
      break;
    }
  }

  memcpy(((char*)base+ip*size), new_element, size);
  return;
}


static IntEdge_ID Insert_Aedge_into_the_edge_array_simple
(
 Tfragment   * frags,
 Tedge       * edges,
 IntEdge_ID  * nedges_delta,
 TIntEdge_ID * next_edge_obj,
 Aedge       * the_edge
 )
{
  const IntEdge_ID iedge = GetNumEdges(edges); // The next available location.
  (*nedges_delta)++;
  EnableRangeVA_Aedge(edges,iedge+1);
  EnableRangeVA_IntEdge_ID(next_edge_obj,iedge+1);
  // Why do not we set the value?
  memcpy( GetVA_Aedge(edges,iedge), the_edge, sizeof(Aedge));
  return iedge;
}

static void Insert_Aedge_into_the_edge_array_wrapper
(
 Tfragment   * frags,
 Tedge       * edges,
 IntEdge_ID  * nedges_delta,
 TIntEdge_ID * next_edge_obj,
 Aedge       * the_edge,
 int           dvt_double_sided_threshold_fragment_end_degree,
 int           con_double_sided_threshold_fragment_end_degree
 )
{
  const QsortCompare compare_edge = get_compare_edge_function();
  const IntFragment_ID avx = the_edge->avx;
  const int asx = the_edge->asx;
  const int ahg = the_edge->ahg;
  const int bhg = the_edge->bhg;
  
  if(is_a_dvt_simple(ahg,bhg)) {  // A dovetail overlap. See get_dvt_edge().
    const int degree=get_seglen_dvt_vertex(frags,avx,asx);
    assert( degree >= 0 );
    assert( degree <= dvt_double_sided_threshold_fragment_end_degree);
    if(degree < dvt_double_sided_threshold_fragment_end_degree) {
      // Not a full list. Add the new record to the array.
      // What about duplicates?
      const IntEdge_ID iedge = Insert_Aedge_into_the_edge_array_simple
        ( frags, edges, nedges_delta, next_edge_obj, the_edge);
      // Added the edge to the edge array.
      set_seglen_dvt_vertex(frags,avx,asx,(degree+1));
      *GetVA_IntEdge_ID(next_edge_obj,iedge)
        = ( degree > 1
            ? get_segend_vertex(frags,avx,asx)
            // A reference to the old top of the list.
            : iedge
            // A reference to self is a sentinal for the end of the list.
            );
      set_segend_vertex(frags,avx,asx,iedge);
    }
    keep_a_list_of_extreme_elements
      ( the_edge,
        GetVA_IntEdge_ID(next_edge_obj,0),
        get_segend_vertex(frags,avx,asx),
        GetVA_Aedge(edges,0),
        sizeof(Aedge),
        compare_edge,
        -1);  // save the minimum ahg elements.
  }else { // A containment overlap a.k.a. not a dovetail overlap.
    const int degree=get_seglen_frc_vertex(frags,avx,asx);
    assert( degree >= 0 );
    assert( degree <= con_double_sided_threshold_fragment_end_degree);
    if(degree < con_double_sided_threshold_fragment_end_degree) {
      // Not a full list. Add the new record to the array.
      // What about duplicates?
      const IntEdge_ID iedge = Insert_Aedge_into_the_edge_array_simple
        ( frags, edges, nedges_delta, next_edge_obj, the_edge);
      // Added the edge to the edge array.
      set_seglen_frc_vertex(frags,avx,asx,(degree+1));
      *GetVA_IntEdge_ID(next_edge_obj,iedge)
        = ( degree > 1
            ? get_segstart_vertex(frags,avx,asx)
            // A reference to the old top of the list.
            : iedge
            // A reference to self is a sentinal for the end of the list.
            );
      set_segstart_vertex(frags,avx,asx,iedge);
    }
    keep_a_list_of_extreme_elements
      ( the_edge,
        GetVA_IntEdge_ID(next_edge_obj,0),
        get_segstart_vertex(frags,avx,asx),
        GetVA_Aedge(edges,0),
        sizeof(Aedge),
        compare_edge,
#ifndef CONTAINMENT_STACKING
        -1  // save the minimum ahg elements.
#else // CONTAINMENT_STACKING
        1   // save the minimum abs(ahg) elements.
#endif // CONTAINMENT_STACKING
        );

    {
      const int degree=get_seglen_frc_vertex(frags,avx,asx);
      assert( degree >= 0 );
      assert( degree <= con_double_sided_threshold_fragment_end_degree);
    }
  }
}


static void insert_dovetail_into_the_edge_array_Aedge
(
 Aedge       * the_edge,
 Tfragment   * frags,
 Tedge       * edges,
 IntEdge_ID  * nedges_delta,
 TIntEdge_ID * next_edge_obj,
 const int dvt_double_sided_threshold_fragment_end_degree,
 const int con_double_sided_threshold_fragment_end_degree,
 const int intrude_with_non_blessed_overlaps_flag
 )
{
  if(
     ( (get_spur_fragment(frags,the_edge->bvx) == FALSE)
       // Only keep this dovetail edge if the distal fragment is
       // non-spur. This isolates the spur fragments.
       ) &&
     ( intrude_with_non_blessed_overlaps_flag
       || the_edge->blessed
       || ! get_blessed_vertex(frags, the_edge->bvx, the_edge->bsx)
       // If the edge is blessed than keep it.
       // If the edge is not blessed and the target fragment-end is not blessed, then keep it.
       )
    //( keep_dvt_overlaps_to_contained_fragments_in_reaper_pass ||
    //(get_con_fragment(frags,the_edge->bvx) == FALSE)
    // Only keep this dovetail edge if the distal fragment is
    // non-contained.
    
     ) {
          Insert_Aedge_into_the_edge_array_wrapper
            ( frags, edges, nedges_delta, next_edge_obj,
              the_edge,
              dvt_double_sided_threshold_fragment_end_degree,
              con_double_sided_threshold_fragment_end_degree
              );
#ifndef PRIOR_INC_COUNT_BLESSED
        inc_raw_dvt_count_vertex(frags,the_edge->avx,the_edge->asx);
#endif // PRIOR_INC_COUNT_BLESSED
  }
}

static void add_overlap_to_graph
(
 Aedge  an_edge,
 Tfragment  frags[],
 Tedge      edges[],
 FragmentHashObject   *afr_to_avx,         // Part of a hash table replacement
 TIntEdge_ID    *next_edge_obj,
 const int dvt_double_sided_threshold_fragment_end_degree,
 const int con_double_sided_threshold_fragment_end_degree,
 const int intrude_with_non_blessed_overlaps_flag,
 IntEdge_ID * novl_dovetail,
 IntEdge_ID * novl_containment,
 IntEdge_ID * nedges_delta
)
{
  const IntFragment_ID iafr = an_edge.avx;
  const int iasx = an_edge.asx;
  const int iahg = an_edge.ahg;

  const IntFragment_ID ibfr = an_edge.bvx;
  const int ibsx = an_edge.bsx;
  const int ibhg = an_edge.bhg;

  const Tnes ines = AS_CGB_SafeCast_cdsInt8_to_Tnes(an_edge.nes);
  const uint32 qua = an_edge.quality;
  const int iinv = an_edge.invalid;
  //const int grangered = an_edge.grangered;
  //const int reflected = an_edge.reflected;
  
  const int is_dovetail = is_a_dvt_simple(iahg,ibhg) ;

  const IntFragment_ID iavx = get_vid_FragmentHash(afr_to_avx,iafr);
  const IntFragment_ID ibvx = get_vid_FragmentHash(afr_to_avx,ibfr);

#ifdef USE_FRAGMENT_SUBSET
  if((iavx == AS_CGB_NOT_SEEN_YET) ||
     (ibvx == AS_CGB_NOT_SEEN_YET))
    return;
#else
  if ((iavx == AS_CGB_NOT_SEEN_YET))
    fprintf(stderr, "Unseen fragment iid="F_IID" is referred to in an overlap.  I assert!\n", iafr);
  if ((ibvx == AS_CGB_NOT_SEEN_YET))
    fprintf(stderr, "Unseen fragment iid="F_IID" is referred to in an overlap.  I assert!\n", ibfr);

  assert(iavx != AS_CGB_NOT_SEEN_YET);
  assert(ibvx != AS_CGB_NOT_SEEN_YET);
#endif

  const int ialn = get_length_fragment(frags,iavx);
  const int ibln = get_length_fragment(frags,ibvx);
  
  switch(ines) {
  case AS_CGB_DOVETAIL_EDGE: (*novl_dovetail)++; break;
  case AS_CGB_CONTAINED_EDGE: (*novl_containment)++; break;
  default:
    fprintf(stderr,"an_edge->nes=%d\n",
            an_edge.nes);
    assert(FALSE);
  }
  
#ifdef REPORT_THIS_BUG
  if(ibmn <= 0 && (AS_CGB_DOVETAIL_EDGE == ines)) {
    fprintf(stderr,"BUG: bmn <= 0 AS_CGB_DOVETAIL_EDGE \n");
    fprintf(stderr,
            "afr,bfr,ahg,bhg,amn,amx=" F_IID "," F_IID ",%d,%d,%d,%d\n",
            iafr,ibfr,iahg,ibhg,iamn,iamx);
  }
  if(ibmn > 0 && (AS_CGB_TO_CONTAINED == ines)) {
    fprintf("BUG: bmn > 0 for AS_CGB_TO_CONTAINED \n");
    fprintf(stderr,
            "afr,bfr,ahg,bhg,amn,amx=" F_IID "," F_IID ",%d,%d,%d,%d\n",
            iafr,ibfr,iahg,ibhg,iamn,iamx);
  }
#endif

  // Are iafr and ibfr in the range of previously read records?
  assert(iafr != ibfr); // No self-overlaps.
  assert(iafr > 0);     // Pointing to an undefined vertex?
  assert(ibfr > 0);
  
  // Are the over-hang distances compatible with the graph topology?
  // If this is true, then the transitive edge marking routine can
  // make time saving assumptions about overlaps that do not need to 
  // be considered.
  
  assert(!((iahg<0)&&(ibhg<0)));  // Art claims this.

  assert( is_a_dvt_simple(iahg,ibhg) ||
          is_a_frc_simple(iahg,ibhg) ||
          is_a_toc_simple(iahg,ibhg) ||
          is_a_dgn_simple(iahg,ibhg) );
  
  assert((!(AS_CGB_DOVETAIL_EDGE == ines)) || is_a_dvt_simple(iahg,ibhg));
  // This is the definition of a dovetail overlap.  That is a dovetail
  // overlap has a positive a-hang and a positive b-hang.
  assert((!(AS_CGB_CONTAINED_EDGE == ines))
         || is_a_toc_simple(iahg,ibhg) || is_a_frc_simple(iahg,ibhg)
         || is_a_dgn_simple(iahg,ibhg) );
  // This is the definition of containment overlaps in
  // the input.

  assert(iafr == get_iid_fragment(frags,iavx));
  assert(ibfr == get_iid_fragment(frags,ibvx));

  if((ialn <= iahg) || (ibln <= ibhg) || (ialn <= -ibhg) || (ibln <= -iahg) ) {
    if(ialn < iahg)      fprintf(stderr,"INPUT ERROR: aln < ahg\n");
    if(ialn < -ibhg)     fprintf(stderr,"INPUT ERROR: aln < -bhg\n");
    if(ialn == iahg)     fprintf(stderr,"INPUT ERROR: aln == ahg\n");
    if(ialn == -ibhg)    fprintf(stderr,"INPUT ERROR: aln == -bhg\n");

    if(ibln < ibhg)      fprintf(stderr,"INPUT ERROR: bln < bhg\n");
    if(ibln < -iahg)     fprintf(stderr,"INPUT ERROR: bln < -ahg\n");
    if(ibln == ibhg)     fprintf(stderr,"INPUT ERROR: bln == bhg\n");
    if(ibln == -iahg)    fprintf(stderr,"INPUT ERROR: bln == -ahg\n");

    fprintf(stderr," afr=" F_IID " bfr=" F_IID " aln=%d bln=%d\n", iafr, ibfr, ialn, ibln);
    fprintf(stderr," avx=" F_IID " bvx=" F_IID " ahg=%d bhg=%d\n", iavx, ibvx, iahg, ibhg);
  }

  assert(ialn > iahg);
  assert(ibln > ibhg);
  assert(ialn > -ibhg);
  assert(ibln > -iahg);
  
  if((get_del_fragment(frags,iavx) == TRUE) ||
     (get_del_fragment(frags,ibvx) == TRUE) ) {
    return;
    // If either of the fragments referenced in the overlap have
    // been deleted, then ignore this overlap during input.
  }
  
  if((get_lab_fragment(frags,iavx) == AS_CGB_REMOVED_BREAKER_FRAG) ||
     (get_lab_fragment(frags,ibvx) == AS_CGB_REMOVED_BREAKER_FRAG) ) {
    // ines=AS_CGB_REMOVED_BY_BREAKER;
    return;  // Remove the overlap from the graph.
  // If either of the fragments referenced in the overlap has
  // been marked as removed by breaker, set the edge label
  }
  
  // Note that we can not count overlaps to deleted and removed
  // fragments.  Contained and spur fragments still have potential
  // mate link information if they are singly placeable.

#ifdef REPORT_DEGENERATE_OVERLAPS
  if( (iahg == 0) && (ibhg == 0) ) {
    fprintf(stdout,
            "Degenerate Overlap " F_IID " %d %d " F_IID " %d %d %d %d\n",
            get_iid_fragment(frags,iavx), iasx, iahg,
            get_iid_fragment(frags,ibvx), ibsx, ibhg,
            qua,
            ines
            );
  }
#endif
  
  {
    Aedge  the_raw_new_edge = {0};

    {
      the_raw_new_edge.avx = iavx;
      the_raw_new_edge.asx = iasx;
      the_raw_new_edge.ahg = iahg;
      
      the_raw_new_edge.bvx = ibvx;
      the_raw_new_edge.bsx = ibsx;
      the_raw_new_edge.bhg = ibhg;
      
      the_raw_new_edge.nes = ines;
      the_raw_new_edge.quality = qua;
      the_raw_new_edge.invalid = iinv;
      the_raw_new_edge.grangered = FALSE;
      the_raw_new_edge.reflected = FALSE;
    }
    
    if( is_dovetail ) {
      {
        Aedge the_edge = the_raw_new_edge;

        insert_dovetail_into_the_edge_array_Aedge
          (
           & the_edge,
           frags,
           edges,
           nedges_delta,
           next_edge_obj,
           dvt_double_sided_threshold_fragment_end_degree,
           con_double_sided_threshold_fragment_end_degree,
           intrude_with_non_blessed_overlaps_flag
           );
#ifdef PRIOR_INC_COUNT_BLESSED
        inc_raw_dvt_count_vertex(frags,the_edge.avx,the_edge.asx);
#endif // PRIOR_INC_COUNT_BLESSED

        if( ! REAPER_VALIDATION) {
          /* 
             The overlapper connects two fragment-ends in an overlap
             relationship:
             
             A    ---------------->
             B          -------------->
             
             The direction mate edge preserves the which fragment-ends are
             in the overlap:
             
             B^c  <----------------
             A^c       <---------------- 
             
          */
          reflect_Aedge( &the_edge, &the_raw_new_edge);
          the_edge.reflected = FALSE;
          // We do not set the reflected bit during the initial file input
          // because we assume that a bug-free ovlStore has both overlaps.
          // However the standard ovl file has only one OVL record per
          // overlap.
          insert_dovetail_into_the_edge_array_Aedge
            (
             & the_edge,
             frags,
             edges,
             nedges_delta,
             next_edge_obj,
             dvt_double_sided_threshold_fragment_end_degree,
             con_double_sided_threshold_fragment_end_degree,
             intrude_with_non_blessed_overlaps_flag
             );
#ifdef PRIOR_INC_COUNT_BLESSED
          inc_raw_dvt_count_vertex(frags,the_edge.bvx,the_edge.bsx);
#endif // PRIOR_INC_COUNT_BLESSED
        }
      }
    } else { // A containment overlap 
      {
        Aedge the_edge_1, the_edge_2;

        assert(AS_CGB_CONTAINED_EDGE == the_raw_new_edge.nes);
        
        if(is_a_frc_simple(the_raw_new_edge.ahg,the_raw_new_edge.bhg)
           || (is_a_dgn_simple(the_raw_new_edge.ahg,the_raw_new_edge.bhg)&&
               (the_raw_new_edge.avx > the_raw_new_edge.bvx))
           ) {
          the_edge_1 = the_raw_new_edge; 
          // Accept as a from-contained overlap.
        } else if(
                  is_a_toc_simple(the_raw_new_edge.ahg,the_raw_new_edge.bhg)
                  || (is_a_dgn_simple(the_raw_new_edge.ahg,the_raw_new_edge.bhg)&&
                      (the_raw_new_edge.avx < the_raw_new_edge.bvx))
                  ) {
          reflect_Aedge( &the_edge_1, &the_raw_new_edge);
          // Convert in-place a to-contained overlap into a from-contained.
          the_edge_1.reflected = FALSE;
          // We do not set the reflected bit during the initial file input
          // because we assume that a bug-free ovlStore has both overlaps.
          // However the standard ovl file has only one OVL record per
          // overlap.
        } else {
          fprintf(stderr,
                  "the_raw_new_edge.ahg,the_raw_new_edge.bhg,the_raw_new_edge.avx,the_raw_new_edge.bvx=\n"
                  "    %d, %d, " F_IID ", " F_IID "\n",
                  the_raw_new_edge.ahg,the_raw_new_edge.bhg,the_raw_new_edge.avx,the_raw_new_edge.bvx
                  );
          assert(FALSE);
        }

        inc_raw_frc_count_fragment(frags,the_edge_1.avx);
        inc_raw_toc_count_fragment(frags,the_edge_1.bvx);

        assert(AS_CGB_CONTAINED_EDGE == the_edge_1.nes);

        if(
           ( intrude_with_non_blessed_overlaps_flag
             || the_edge_1.blessed
             || ! get_blessed_vertex(frags, the_edge_1.bvx, the_edge_1.bsx)
             // If the edge is blessed than keep it.
             // If the edge is not blessed and the target fragment-end is not blessed,
             // then keep it.
             )
           ) {
          Insert_Aedge_into_the_edge_array_wrapper
            ( frags, edges, nedges_delta, next_edge_obj,
              &the_edge_1,
              dvt_double_sided_threshold_fragment_end_degree,
              con_double_sided_threshold_fragment_end_degree
              );
        }
        
        /* The overlapper guarantees that ahg>0 in the reported
           overlaps.
           
           A    ---------------->
           B          ------->...
           
           So create the granger mate edge:
           A^c  <----------------
           B^c     <-------......  
        */
        granger_Aedge( &the_edge_2, &the_edge_1);
        if(
           ( intrude_with_non_blessed_overlaps_flag
             || the_edge_2.blessed
             || ! get_blessed_vertex(frags, the_edge_2.bvx, the_edge_2.bsx)
             // If the edge is blessed than keep it.
             // If the edge is not blessed and the target fragment-end is not blessed,
             // then keep it.
             )
           ) {
          Insert_Aedge_into_the_edge_array_wrapper
            ( frags, edges, nedges_delta, next_edge_obj,
              &the_edge_2,
              dvt_double_sided_threshold_fragment_end_degree,
              con_double_sided_threshold_fragment_end_degree
            );
        }
      }
    }
  }
}

/****************************************************************************/

void
process_gkp_store_for_fragments(char *gkpStoreName,
                                Tfragment   *frags,
                                Tedge       *edges,
                                FragmentHashObject *afr_to_avx,
                                IntFragment_ID    *min_frag_iid,
                                IntFragment_ID    *max_frag_iid) {

  fragRecord       *fr = new_fragRecord();
  FragStream       *fs = NULL;

  IntFragment_ID    iid = 0;
  IntFragment_ID    vid = 0;

  GateKeeperStore   *gkp = openGateKeeperStore(gkpStoreName, FALSE);
  if (gkp == NULL) {
    fprintf(stderr, "Failed to open %s\n", gkpStoreName);
    exit(1);
  }

  assert(0 == GetNumFragments(frags));

  *min_frag_iid = getLastElemFragStore(gkp) + 1;
  *max_frag_iid = 0;

  fs = openFragStream(gkp, FRAG_S_INF);

  //unsigned int firstElem = getFirstElemFragStore(gkp);
  ///unsigned int lastElem  = getLastElemFragStore(gkp) + 1;
  //resetFragStream(fs, firstElem, lastElem);

  while (nextFragStream(fs, fr)) {
    if (getFragRecordIsDeleted(fr) == FALSE) {
      iid = getFragRecordIID(fr);
  
      *min_frag_iid = MIN(*min_frag_iid, iid);
      *max_frag_iid = MAX(*max_frag_iid, iid);

      //  Argh!  This needs to be here, other code depends on the
      //  range of the VA being the number of fragments.
      //
      EnableRangeVA_Afragment(frags, vid + 1);

      set_vid_FragmentHash(afr_to_avx, iid, vid);

      set_iid_fragment(frags, vid, iid);
      set_cid_fragment(frags, vid, iid);
      set_typ_fragment(frags, vid, AS_READ);
      set_del_fragment(frags, vid, FALSE);
      set_length_fragment(frags, vid,
                          getFragRecordClearRegionEnd  (fr, AS_READ_CLEAR_OBT) -
                          getFragRecordClearRegionBegin(fr, AS_READ_CLEAR_OBT));

      // Assume that each fragment spans a chunk.
      set_lab_fragment(frags, vid, AS_CGB_UNLABELED_FRAG);

      // A flag specifying if this fragment is known to be contained.
      // Set default to non-contained.
      set_con_fragment(frags, vid, FALSE);

      // Zero if this flag is not contained, but equal to the
      // fragment containing this one in a unitig layout graph.
      set_container_fragment(frags, vid,0); 

#ifdef STORE_BRANCH_POINTS_AT_FRAGMENT
      set_bpt_vertex(frags, vid, FALSE, 0); // Sentinal for no branch point.
      set_bpt_vertex(frags, vid, TRUE, 0);
#endif // STORE_BRANCH_POINTS_AT_FRAGMENT

      // Set the counts of raw overlaps seen to zero.
      set_raw_dvt_count_vertex(frags, vid, FALSE, 0);
      set_raw_dvt_count_vertex(frags, vid, TRUE, 0);
      set_raw_toc_count_fragment(frags, vid, 0);
      set_raw_frc_count_fragment(frags, vid, 0);

      // Initialize the lists for the edge trimming.
      set_seglen_dvt_vertex(frags, vid, FALSE, 0);
      set_seglen_frc_vertex(frags, vid, FALSE, 0);
      set_seglen_dvt_vertex(frags, vid, TRUE, 0);
      set_seglen_frc_vertex(frags, vid, TRUE, 0);

      set_blessed_vertex(frags, vid, FALSE, FALSE);
      set_blessed_vertex(frags, vid, TRUE, FALSE);

      set_src_fragment(frags, vid, 0);

      vid++;
    }
  }

  del_fragRecord(fr);
  closeFragStream(fs);
  closeGateKeeperStore(gkp);
}




void process_ovl_store(char * OVL_Store_Path,
                       Tfragment  frags[],
                       Tedge      edges[],
                       FragmentHashObject *afr_to_avx,
                       TIntEdge_ID     next_edge_obj[],
                       const int dvt_double_sided_threshold_fragment_end_degree,
                       const int con_double_sided_threshold_fragment_end_degree,
                       const int intrude_with_non_blessed_overlaps_flag,
                       const uint32 overlap_error_threshold) {
  OverlapStore  *ovs;
  OVSoverlap     olap;

  IntEdge_ID novl_dovetail = 0;
  IntEdge_ID novl_containment = 0;
  IntEdge_ID nedges_delta = 0;

  ovs = AS_OVS_openOverlapStore(OVL_Store_Path);

  //  Copy the information in  (* olap)  into  (* an_edge)  with
  //  appropriate conversions.
  //
  //  This improper dovetail overlap (a_hang<0)&&(b_hang<0)
  //  A_frag    >>>>>>>>>>>
  //  B_frag  >>>>>>>>
  //  becomes a proper dovetail overlap (ahg>0)&&(bhg>0)
  //  A_frag  <<<<<<<<<<<
  //  B_frag       <<<<<<<<
  //
  //  This improper to-contained overlap (a_hang==0)&&(b_hang<0)
  //  A_frag  >>>>>>>>>>
  //  B_frag  >>>>>>>...
  //  becomes a proper to-contained overlap (ahg>0)&&(bhg==0)
  //  A_frag  <<<<<<<<<<
  //  B_frag     <<<<<<<
  //
  //  This improper from-contained overlap (a_hang<0)&&(b_hang==0)
  //  A_frag  ...>>>>>>>
  //  B_frag  >>>>>>>>>>
  //  becomes a proper from-contained overlap (ahg==0)&&(bhg>0)
  //  A_frag  <<<<<<<
  //  B_frag  <<<<<<<<<<
  //
  //  A degenerate overlap (a_hang==0)&&(b_hang==0)
  //  A_frag  >>>>>>>>>>
  //  B_frag  >>>>>>>>>>
  while  (AS_OVS_readOverlapFromStore(ovs, &olap, AS_OVS_TYPE_OVL)) {
    Aedge  e = {0};

    int improper = (((olap.dat.ovl.a_hang <  0) && (olap.dat.ovl.b_hang <  0)) ||
                    ((olap.dat.ovl.a_hang == 0) && (olap.dat.ovl.b_hang <  0)) ||
                    ((olap.dat.ovl.a_hang <  0) && (olap.dat.ovl.b_hang == 0)));

    e.avx = olap.a_iid;
    e.asx = !improper;
    e.ahg = (improper ? -olap.dat.ovl.b_hang : olap.dat.ovl.a_hang);
  
    e.bvx = olap.b_iid;
    e.bsx = (!improper) ^ (!olap.dat.ovl.flipped);
    e.bhg = (improper ? -olap.dat.ovl.a_hang : olap.dat.ovl.b_hang);
  
    e.nes       = (is_a_dvt_simple(e.ahg, e.bhg) ? AS_CGB_DOVETAIL_EDGE : AS_CGB_CONTAINED_EDGE);
    e.quality   = olap.dat.ovl.corr_erate;
    e.invalid   = FALSE;
    e.reflected = FALSE;
    e.grangered = FALSE;
    e.blessed   = FALSE;

    assert((e.ahg > 0) || (e.bhg > 0) || ((e.ahg == 0) && (e.bhg == 0)));

    // Avoid entering the containment overlap twice.
    if(((AS_CGB_CONTAINED_EDGE == e.nes) &&
        (is_a_frc_simple(e.ahg,e.bhg))))
      continue;

    if (e.quality <= overlap_error_threshold)
      add_overlap_to_graph(e,
                           frags,
                           edges,
                           afr_to_avx,
                           next_edge_obj,
                           dvt_double_sided_threshold_fragment_end_degree,
                           con_double_sided_threshold_fragment_end_degree,
                           intrude_with_non_blessed_overlaps_flag,
                           &novl_dovetail,
                           &novl_containment,
                           &nedges_delta);
  }
 
  AS_OVS_closeOverlapStore(ovs);

  fprintf(stderr,"novl_dovetail    = " F_IID "\n", novl_dovetail);
  fprintf(stderr,"novl_containment = " F_IID "\n", novl_containment);
  fprintf(stderr,"nedges_delta     = " F_IID "\n", nedges_delta);
}

/****************************************************************************/

void input_messages_from_a_file(FILE       *fovl,
                                Tfragment  frags[],
                                Tedge      edges[],
                                FragmentHashObject *afr_to_avx,
                                TIntEdge_ID       *next_edge_obj,
                                const int dvt_double_sided_threshold_fragment_end_degree,
                                const int con_double_sided_threshold_fragment_end_degree,
                                const int intrude_with_non_blessed_overlaps_flag,
                                const uint32 overlap_error_threshold) {
  /* Input a batch of fragment read and overlap records from a stream. */
  /* Keep a copy of the number of fragments and edges before
     the new data is read in. */
  const IntEdge_ID nedge_old = GetNumEdges(edges);
  
  IntEdge_ID nedge_new = nedge_old;
  
  IntEdge_ID novl_dovetail=0;  /* The number of overlap records read. */
  IntEdge_ID novl_containment=0;
  IntEdge_ID novl_degenerate=0;

  IntEdge_ID nedge_delta=0;   
  
  /* It is assumed that in the overlap records that new fragments
     point to old fragments.  */

  GenericMesg  *pmesg;

  while( EOF != ReadProtoMesg_AS(fovl, &pmesg)) {
    const MessageType imesgtype = pmesg->t;

    //
    //  bubble popper writes this input, otherwise, it's unused
    //

    if (pmesg->t == MESG_OVL) {

      //  This improper dovetail overlap (a_hang<0)&&(b_hang<0)
      //  A_frag    >>>>>>>>>>>
      //  B_frag  >>>>>>>>
      //  becomes a proper dovetail overlap (ahg>0)&&(bhg>0)
      //  A_frag  <<<<<<<<<<<
      //  B_frag       <<<<<<<<
      //
      //  This improper to-contained overlap (a_hang==0)&&(b_hang<0)
      //  A_frag  >>>>>>>>>>
      //  B_frag  >>>>>>>...
      //  becomes a proper to-contained overlap (ahg>0)&&(bhg==0)
      //  A_frag  <<<<<<<<<<
      //  B_frag     <<<<<<<
      //
      //  This improper from-contained overlap (a_hang<0)&&(b_hang==0)
      //  A_frag  ...>>>>>>>
      //  B_frag  >>>>>>>>>>
      //  becomes a proper from-contained overlap (ahg==0)&&(bhg>0)
      //  A_frag  <<<<<<<
      //  B_frag  <<<<<<<<<<
      //
      //  A degenerate overlap (a_hang==0)&&(b_hang==0)
      //  A_frag  >>>>>>>>>>
      //  B_frag  >>>>>>>>>>

      OverlapMesg *o = (OverlapMesg *)pmesg->m;
      Aedge        e = {0};

      int improper = ((o->ahg <  0) && (o->bhg <  0)) || ((o->ahg == 0) && (o->bhg <  0)) || ((o->ahg < 0) && (o->bhg == 0)) ;

      e.avx = o->aifrag;
      e.asx = improper ^ ((o->orientation == AS_NORMAL) || (o->orientation == AS_INNIE));
      e.ahg = (improper) ? -o->bhg : o->ahg;
  
      e.bvx = o->bifrag;
      e.bsx = e.asx ^ !((o->orientation == AS_INNIE)  || (o->orientation == AS_OUTTIE));
      e.bhg = (improper) ? -o->ahg : o->bhg;
  
      e.nes       = (is_a_dvt_simple(e.ahg, e.bhg)) ? AS_CGB_DOVETAIL_EDGE : AS_CGB_CONTAINED_EDGE;
      e.quality   = o->quality;
      e.invalid   = FALSE;
      e.grangered = FALSE;
      e.reflected = FALSE;
      e.blessed   = FALSE;

      assert( (e.ahg>0) || (e.bhg>0) || ((e.ahg == 0) && (e.bhg == 0)) );

      if (e.quality < overlap_error_threshold)
        add_overlap_to_graph(e,
                             frags,
                             edges,
                             afr_to_avx,
                             next_edge_obj,
                             dvt_double_sided_threshold_fragment_end_degree,
                             con_double_sided_threshold_fragment_end_degree,
                             intrude_with_non_blessed_overlaps_flag,
                             &novl_dovetail,
                             &novl_containment,
                             &nedge_delta);
    } else {
      fprintf(stderr,"Unexpected message type %d (%s)\n",imesgtype, MessageTypeName[imesgtype]);
      assert(FALSE);
    }
  }


  fprintf(stderr,"Input %10" F_IIDP " OVL records (skipped %10"F_IIDP" degenerate).\n",novl_dovetail+novl_containment, novl_degenerate);
  fprintf(stderr,"      %10" F_IIDP " OVL dovetail records.\n",novl_dovetail);
  fprintf(stderr,"      %10" F_IIDP " OVL containment records.\n",novl_containment);

  nedge_new = nedge_old + nedge_delta;
  assert(nedge_new == GetNumEdges(edges));
}
