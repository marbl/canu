
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
= "$Id: AS_FGB_io.c,v 1.13 2007-01-29 20:40:58 brianwalenz Exp $";
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

/*********************************************************************/
/* System include files */
#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/*********************************************************************/
/* Local include files */
#include "AS_global.h"
#include "AS_UTL_version.h"
#include "AS_MSG_pmesg.h"
#include "AS_CGB_all.h"
#include "OlapStoreOVL.h"  // In cds/AS/src/AS_OVL
#include "AS_CGB_util.h"

/****************************************************************************/
/* File Scope Globals */

static int TIMINGS = TRUE;
static int TIMINGS1 = FALSE;
int REAPER_VALIDATION = FALSE;

/************************************************************************/
/* File Scope Preprocessor Globals */

#undef REPORT_THIS_BUG
#undef ZAPZAP
#undef DEBUG113
#undef DEBUG280
#undef DEBUG_MOSQ_01

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


static void add_OFGMesg_to_graph
( OFGMesg      * ofg_mesg,
  VA_TYPE(char) * the_ofg_source,
  Tfragment    * frags,
  IntFragment_ID nfrag_base,
  FragmentHashObject   *afr_to_avx,  // part of a hash table replacement
  IntFragment_ID *Pnofg,
  IntFragment_ID *Pmin_frag_iid,
  IntFragment_ID *Pmax_frag_iid,
  VA_TYPE(char)  *frag_annotations
  )
{
  const Fragment_ID    uid = ofg_mesg->eaccession;
  const IntFragment_ID iid = ofg_mesg->iaccession;
  IntFragment_ID vid = AS_CGB_NOT_SEEN_YET;
  int a_new_fragment = FALSE;
  
  assert(ofg_mesg->action == AS_ADD ||
         ofg_mesg->action == AS_DELETE );

  vid = get_vid_FragmentHash(afr_to_avx,iid);
  
  if(ofg_mesg->action == AS_DELETE) {
    // Has this fragment been processed before???
    if( vid != AS_CGB_NOT_SEEN_YET ) {
      assert( vid >=0 );
      assert( vid < nfrag_base + (*Pnofg));
      set_del_fragment(frags,vid,TRUE);
      //set_con_fragment(frags,vid,FALSE); // Set default to non-contained
    } else {
#ifdef REPORT_DELETED_UNDEF_FRAGS
      fprintf(stderr,"Deleting an undefined fragment (uid,iid)=("
              F_UID "," F_IID ")\n", uid,iid);
#endif
    }
  }
  
  // Has this fragment been processed before???
  a_new_fragment = (vid == AS_CGB_NOT_SEEN_YET);
#ifdef ZAPZAP
  fprintf(stderr,"AS_CGB_NOT_SEEN_YET=" F_S32 "\n",AS_CGB_NOT_SEEN_YET);
  fprintf(stderr,"vid=" F_IID "\n",vid);
  fprintf(stderr,"nofg=" F_IID "\n",(*Pnofg));
  fprintf(stderr,"a_new_fragment=%d\n",a_new_fragment);
  fprintf(stderr,"get_vid_FragmentHash(afr_to_avx,iid)=" F_IID "\n",get_vid_FragmentHash(afr_to_avx,iid));
#endif // ZAPZAP

  assert( (a_new_fragment == TRUE) ||
          (get_vid_FragmentHash(afr_to_avx,iid) <= nfrag_base+(*Pnofg)) );
  
  if( ofg_mesg->action == AS_ADD) {
    
    /* Put the record where it belongs in the array.
       This array is indexed by the overlaps. */
    
    FragType typ = ofg_mesg->type;
    CDS_COORD_t length = ofg_mesg->clear_rng.end - ofg_mesg->clear_rng.bgn;
    
    // imate  = ofg_mesg->u.read.imate;
    // distance = ofg_mesg->u.read.idistance;
    if((*Pmin_frag_iid) == 0 &&
       (*Pmax_frag_iid) == 0 ) {
      (*Pmin_frag_iid) = iid;
      (*Pmax_frag_iid) = iid;
    } else {
      (*Pmin_frag_iid) = MIN((*Pmin_frag_iid),iid);
      (*Pmax_frag_iid) = MAX((*Pmax_frag_iid),iid);
    }
    
    
    // This is a new fragment!
    
    if( TRUE == a_new_fragment ) {
      vid = nfrag_base+(*Pnofg);
      set_vid_FragmentHash(afr_to_avx,iid,vid);
      (*Pnofg)++;
      EnableRangeVA_Afragment(frags,nfrag_base+(*Pnofg));
      
      set_uid_fragment(frags,vid,uid);
      set_iid_fragment(frags,vid,iid);
      set_cid_fragment(frags,vid,iid);
      set_typ_fragment(frags,vid,typ);
      set_del_fragment(frags,vid,FALSE); // Default to good.
      set_length_fragment(frags,vid,length);
      // set_imate_fragment(frags,vid,imate);
      // set_distance_fragment(frags,vid,distance);
      
      set_lab_fragment(frags,vid,AS_CGB_UNLABELED_FRAG);
      // Assume that each fragment spans a chunk.
      set_con_fragment(frags,vid,FALSE); // Set default to non-contained.
      // A flag specifying if this fragment is known to be contained.
      set_container_fragment(frags,vid,0); 
      // Zero if this flag is not contained, but equal to the
      // fragment containing this one in a unitig layout graph.
#ifdef STORE_BRANCH_POINTS_AT_FRAGMENT
      set_bpt_vertex(frags,vid,FALSE,0); // Sentinal for no branch point.
      set_bpt_vertex(frags,vid,TRUE,0);
#endif // STORE_BRANCH_POINTS_AT_FRAGMENT
      
      set_raw_dvt_count_vertex(frags,vid,FALSE,0);
      set_raw_dvt_count_vertex(frags,vid,TRUE,0);
      set_raw_toc_count_fragment(frags,vid,0);
      set_raw_frc_count_fragment(frags,vid,0);
      // Set the counts of raw overlaps seen to zero.
      
      set_seglen_dvt_vertex(frags,vid,FALSE,0);
      set_seglen_frc_vertex(frags,vid,FALSE,0);
      set_seglen_dvt_vertex(frags,vid,TRUE,0);
      set_seglen_frc_vertex(frags,vid,TRUE,0);
      // Initialize the lists for the edge trimming.

      set_blessed_vertex(frags,vid,FALSE,FALSE);
      set_blessed_vertex(frags,vid,TRUE,FALSE);
#if 0
      set_frozen_dvt_vertex(frags,vid,FALSE,0);
      set_frozen_dvt_vertex(frags,vid,TRUE,0);
      set_frozen_frc_vertex(frags,vid,FALSE,0);
      set_frozen_frc_vertex(frags,vid,TRUE,0);
#endif

      if(NULL != frag_annotations) {
        char * the_text =
          ( //ofg_mesg_source_as_offsets_in_frgsrc &&
           (NULL != the_ofg_source)
           ? GetVA_char(the_ofg_source, (size_t)ofg_mesg->source)
           : ofg_mesg->source );
        
        /* Save the simulator annotations. */
        size_t nsrc = GetNumVA_char(frag_annotations);
        size_t nsource = nsrc+(NULL != the_text ? strlen(the_text) : 0)+1; 
        // Remember the terminal null character.
        // fprintf(fp_frgsrc,"IO: " F_SIZE_T " <%s>\n",nsrc,the_text);
        EnableRangeVA_char(frag_annotations,nsource);
        if(NULL != the_text) {
          strcpy(Getchar(frag_annotations,nsrc),the_text);
        }
        set_src_fragment(frags,vid,nsrc);
        
        // ofg_mesg->source = nsrc; ofg_mesg_source_as_offsets_in_frgsrc = TRUE;
        // Replace the C-pointer to absolute address with an index into a VA.
      } else {
        set_src_fragment(frags,vid,0);
      }

#ifdef SAVE_FRGSRC_IN_VA
      if(NULL != ofg_mesg->source ) {
        /* Save the simulator annotations. */
        size_t nsrc = GetNumVA_char(frag_annotations);
        size_t nsource = nsrc+strlen(ofg_mesg->source)+1; 
        // fprintf(fp_frgsrc,"IO: " F_SIZE_T " <%s>\n",nsrc,ofg_mesg->source);
        /* remember the terminal null character. */
        EnableRangeVA_char(frag_annotations,nsource);
        strcpy(Getchar(frag_annotations,nsrc),ofg_mesg->source);
        set_src_fragment(frags,vid,nsrc);

        // ofg_mesg->source = nsrc;
        // Replace the C-pointer to absolute address with an index into a VA.
      } else {
        set_src_fragment(frags,vid,0);
        // Zero is a special offset meaning no source field for this
        // fragment.
      }
#endif // SAVE_FRGSRC_IN_VA        
    } else {
      // Demand that it is a duplicate definition!!
      const int duplicate =
        (uid == get_uid_fragment(frags,vid)) &&
        (iid == get_iid_fragment(frags,vid)) &&
        (typ == get_typ_fragment(frags,vid)) &&
        (length == get_length_fragment(frags,vid)); // a duplicate fragment?
        // get_del_fragment(frags,vid)
      if( duplicate ) {
        fprintf(stderr,"IGNORED DUPLICATE FRAG (" F_UID "," F_IID ")\n",
                uid, iid);
      } else {
        fprintf(stderr,"IGNORED INCONSISENT FRAG (" F_UID "," F_IID ")\n",
                uid, iid);
      }
      assert(TRUE == duplicate);
    }
  }
}

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
  //const IntFragment_ID bvx = the_edge->bvx;
  //const int bsx = the_edge->bsx;
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
    {
      const IntEdge_ID iedge = get_segend_vertex(frags,avx,asx);
      keep_a_list_of_extreme_elements
        ( the_edge,
          GetVA_IntEdge_ID(next_edge_obj,0),
          iedge,
          GetVA_Aedge(edges,0),
          sizeof(Aedge),
          compare_edge,
          -1);  // save the minimum ahg elements.
    }
  }else { // A containment overlap a.k.a. not a dovetail overlap.
    {
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
    }
    {
      const IntEdge_ID iedge = get_segstart_vertex(frags,avx,asx);
      keep_a_list_of_extreme_elements
        ( the_edge,
          GetVA_IntEdge_ID(next_edge_obj,0),
          iedge,
          GetVA_Aedge(edges,0),
          sizeof(Aedge),
          compare_edge,
#ifndef CONTAINMENT_STACKING
          -1  // save the minimum ahg elements.
#else // CONTAINMENT_STACKING
          1   // save the minimum abs(ahg) elements.
#endif // CONTAINMENT_STACKING
          );
    }
#if 1
    {
      const int degree=get_seglen_frc_vertex(frags,avx,asx);
      assert( degree >= 0 );
      assert( degree <= con_double_sided_threshold_fragment_end_degree);
    }
#endif    
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
#ifdef STORE_OVERLAP_EXTREMES
  const int iamn = an_edge.amn;
  const int iamx = an_edge.amx;
#endif // STORE_OVERLAP_EXTREMES

  const IntFragment_ID ibfr = an_edge.bvx;
  const int ibsx = an_edge.bsx;
  const int ibhg = an_edge.bhg;
#ifdef STORE_OVERLAP_EXTREMES
  const int ibmn = an_edge.bmn;
  const int ibmx = an_edge.bmx;
#endif // STORE_OVERLAP_EXTREMES

  const Tnes ines = AS_CGB_SafeCast_cdsInt8_to_Tnes(an_edge.nes);
  const CGB_ERATE_TYPE qua = an_edge.quality;
  const int iinv = an_edge.invalid;
  //const int grangered = an_edge.grangered;
  //const int reflected = an_edge.reflected;
  
  const int is_dovetail =  is_a_dvt_simple(iahg,ibhg);

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
  
#ifdef STORE_OVERLAP_EXTREMES
#ifdef REPORT_THIS_BUG
  if(iamn < 0) {
    fprintf(stderr,"BUG: amn < 0, adjusting to zero....\n");
    fprintf(stderr,"afr,bfr,ahg,bhg,amn,amx=" F_IID "," F_IID ",%d,%d,%d,%d\n",
            iafr,ibfr,iahg,ibhg,iamn,iamx);
  }
#endif // REPORT_THIS_BUG
  assert( iamn >= 0);
  assert( iamn <= iahg);
  assert( iamx >= iahg);
#endif // STORE_OVERLAP_EXTREMES

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

  {
    if((ialn <= iahg) || (ibln <= ibhg) || (ialn <= -ibhg) || (ibln <= -iahg) ) {
      FILE *ferror = stdout;
      if(ialn < iahg) {
        fprintf(ferror,"INPUT ERROR: aln < ahg\n");
      }
      if(ialn < -ibhg) {
        fprintf(ferror,"INPUT ERROR: aln < -ahg\n");
      }
      if(ialn == iahg) {
        fprintf(ferror,"INPUT ERROR: aln == ahg\n");
      }
      if(ialn == -ibhg) {
        fprintf(ferror,"INPUT ERROR: aln == -bhg\n");
      }
      if(ibln < ibhg) {
        fprintf(ferror,"INPUT ERROR: bln < bhg\n");
      }
      if(ibln < -iahg) {
        fprintf(ferror,"INPUT ERROR: bln < -bhg\n");
      }
      if(ibln == ibhg) {
        fprintf(ferror,"INPUT ERROR: bln == bhg\n");
      }
      if(ibln == -iahg) {
        fprintf(ferror,"INPUT ERROR: bln == bhg\n");
      }
      fprintf(ferror," afr=" F_IID " bfr=" F_IID " aln=%d bln=%d\n",
              iafr, ibfr, ialn, ibln);
      fprintf(ferror," avx=" F_IID " bvx=" F_IID " ahg=%d bhg=%d\n",
              iavx, ibvx, iahg, ibhg);
    }
  }
  assert( ialn > iahg);
  assert( ibln > ibhg);
  assert( ialn > -ibhg);
  assert( ibln > -iahg);
  
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
            "Degenerate Overlap " F_IID " %d %d " F_IID " %d %d " CGB_ERATE_FORMAT " %d\n",
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

#ifdef STORE_OVERLAP_EXTREMES
      the_raw_new_edge.amn = iamn;
      the_raw_new_edge.amx = iamx;
      the_raw_new_edge.bmn = ibmn;
      the_raw_new_edge.bmx = ibmx;
#endif // STORE_OVERLAP_EXTREMES
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

static void Convert_OverlapMesg_to_Aedge
(
 const OverlapMesg * const ovl_mesg,
 // The raw overlap input from proto-IO.
 Aedge * const the_new_raw_edge
 // The raw edge converted to the internal overlap edge.
 )
{
  const IntFragment_ID a_frag = ovl_mesg->aifrag; 
  const IntFragment_ID b_frag = ovl_mesg->bifrag;
  /* Assembler internal Fragment ids. */
  const int a_hang = ovl_mesg->ahg;
  const int b_hang = ovl_mesg->bhg;
  const int regular =
    (ovl_mesg->orientation == AS_NORMAL) ||
    (ovl_mesg->orientation == AS_INNIE); 
  const int flipped =
    (ovl_mesg->orientation == AS_INNIE) ||
    (ovl_mesg->orientation == AS_OUTTIE);
  
  //const int iamn = MAX(0,ovl_mesg->min_offset);
  //const int iamx = ovl_mesg->max_offset;
  // This is an approximation.
  //const int ibmn = (ibhg - (iamx-iahg));
  //const int ibmx = (ibhg + (iahg-iamn));

  const CGB_ERATE_TYPE qua = ovl_mesg->quality;
  // Assume that the overlap is valid unless proven otherwise.

  const int improper = ((a_hang <  0) && (b_hang <  0)) 
    || ((a_hang == 0) && (b_hang <  0)) || ((a_hang <  0) && (b_hang == 0)) ;
  const IntFragment_ID iafr = a_frag;
  const IntFragment_ID ibfr = b_frag;
  const int iahg = (improper ? -b_hang : a_hang);
  const int ibhg = (improper ? -a_hang : b_hang);
  const int iasx = improper ^ regular;
  const int ibsx = iasx ^ (!flipped);

  //  This improper dovetail overlap (a_hang<0)&&(b_hang<0)
  //  A_frag    >>>>>>>>>>>
  //  B_frag  >>>>>>>>
  //  becomes a proper dovetail overlap (ahg>0)&&(bhg>0)
  //  A_frag  <<<<<<<<<<<
  //  B_frag       <<<<<<<<

  //  This improper to-contained overlap (a_hang==0)&&(b_hang<0)
  //  A_frag  >>>>>>>>>>
  //  B_frag  >>>>>>>...
  //  becomes a proper to-contained overlap (ahg>0)&&(bhg==0)
  //  A_frag  <<<<<<<<<<
  //  B_frag     <<<<<<<

  //  This improper from-contained overlap (a_hang<0)&&(b_hang==0)
  //  A_frag  ...>>>>>>>
  //  B_frag  >>>>>>>>>>
  //  becomes a proper from-contained overlap (ahg==0)&&(bhg>0)
  //  A_frag  <<<<<<<
  //  B_frag  <<<<<<<<<<

  //  A degenerate overlap (a_hang==0)&&(b_hang==0)
  //  A_frag  >>>>>>>>>>
  //  B_frag  >>>>>>>>>>

  const int is_dovetail = is_a_dvt_simple(iahg,ibhg);
#ifdef STORE_OVERLAP_EXTREMES
  const int iamn = iahg;
  const int iamx = iahg;
  const int ibmn = ibhg;
  const int ibmx = ibhg;
#endif
  const int iinv = FALSE;
  const Tnes ines = ( is_dovetail
		      ? AS_CGB_DOVETAIL_EDGE
                      : AS_CGB_CONTAINED_EDGE );

  assert( (iahg>0) || (ibhg>0) || ((iahg == 0) && (ibhg == 0)) );
  
  the_new_raw_edge->avx = iafr;
  the_new_raw_edge->asx = iasx;
  the_new_raw_edge->ahg = iahg;
  
  the_new_raw_edge->bvx = ibfr;
  the_new_raw_edge->bsx = ibsx;
  the_new_raw_edge->bhg = ibhg;
  
  the_new_raw_edge->nes = ines;
  the_new_raw_edge->quality = qua;
  the_new_raw_edge->invalid = iinv;
  the_new_raw_edge->grangered = FALSE;
  the_new_raw_edge->reflected = FALSE;
  the_new_raw_edge->blessed = FALSE;

#ifdef STORE_OVERLAP_EXTREMES
  the_new_raw_edge->amn = iamn;
  the_new_raw_edge->amx = iamx;
  the_new_raw_edge->bmn = ibmn;
  the_new_raw_edge->bmx = ibmx;
#endif // STORE_OVERLAP_EXTREMES
}


static void add_OverlapMesg_to_graph
(
 OverlapMesg * ovl_mesg,
 VA_TYPE(char) * the_ovl_source,
 Tfragment  frags[],
 Tedge      edges[],
 FragmentHashObject   *afr_to_avx,   // part of a hash table replacement
 TIntEdge_ID     next_edge_obj[],
 const int dvt_double_sided_threshold_fragment_end_degree,
 const int con_double_sided_threshold_fragment_end_degree,
 const int intrude_with_non_blessed_overlaps_flag,
 const CGB_ERATE_TYPE overlap_error_threshold,
 IntEdge_ID * novl_dovetail,
 IntEdge_ID * novl_containment,
 IntEdge_ID * novl_degenerate,
 IntEdge_ID * nedges_delta
)
{
 Aedge  an_edge = {0};
 Convert_OverlapMesg_to_Aedge( ovl_mesg, &an_edge);

 if(overlap_error_threshold < an_edge.quality ) { return;}

 add_overlap_to_graph
   (
    an_edge,
    frags,
    edges,
    afr_to_avx, // Part of a hash table replacement
    next_edge_obj,
    dvt_double_sided_threshold_fragment_end_degree,
    con_double_sided_threshold_fragment_end_degree,
    intrude_with_non_blessed_overlaps_flag,
    novl_dovetail,
    novl_containment,
    nedges_delta
    );
}

/****************************************************************************/

static void Convert_Olap_To_Aedge
(
 const Long_Olap_Data_t * const olap,
 // The raw edge input from the overlap store.
 Aedge * const the_new_raw_edge
 // The raw edge as an Unitigger internal overlap edge.
 )
  
  //  Copy the information in  (* olap)  into  (* the_new_raw_edge)  with
  //  appropriate conversions.  Both must have been allocated by
  //  the calling routine.
  
{
  // Beginning of the input data:
  const IntFragment_ID a_frag = (olap -> a_iid);
  const IntFragment_ID b_frag = (olap -> b_iid);
  const int a_hang = (olap -> a_hang);
  const int b_hang = (olap -> b_hang);
  const int regular = TRUE;
  //(ovl_mesg->orientation == AS_NORMAL) ||
  //(ovl_mesg->orientation == AS_INNIE); 
  const int flipped = (olap -> flipped);

#ifndef USE_FLOAT_ERATE
# ifndef USE_CORRECTED_ERROR_RATE
  const CGB_ERATE_TYPE qua = olap -> orig_erate;
# else // USE_CORRECTED_ERROR_RATE
  const CGB_ERATE_TYPE qua = olap -> corr_erate;
# endif // USE_CORRECTED_ERROR_RATE
#else // USE_FLOAT_ERATE
# ifndef USE_CORRECTED_ERROR_RATE
  const CGB_ERATE_TYPE qua = Expand_Quality (olap -> orig_erate);
# else // USE_CORRECTED_ERROR_RATE
  const CGB_ERATE_TYPE qua = Expand_Quality (olap -> corr_erate);
# endif // USE_CORRECTED_ERROR_RATE
#endif // USE_FLOAT_ERATE
    
  // Ending of the input data.

  const int improper = ((a_hang <  0) && (b_hang <  0)) 
    || ((a_hang == 0) && (b_hang <  0)) || ((a_hang <  0) && (b_hang == 0)) ;
  const IntFragment_ID iafr = a_frag;
  const IntFragment_ID ibfr = b_frag;
  const int iahg = (improper ? -b_hang : a_hang);
  const int ibhg = (improper ? -a_hang : b_hang);
  const int iasx = improper ^ regular;
  const int ibsx = iasx ^ (!flipped);

  //  This improper dovetail overlap (a_hang<0)&&(b_hang<0)
  //  A_frag    >>>>>>>>>>>
  //  B_frag  >>>>>>>>
  //  becomes a proper dovetail overlap (ahg>0)&&(bhg>0)
  //  A_frag  <<<<<<<<<<<
  //  B_frag       <<<<<<<<

  //  This improper to-contained overlap (a_hang==0)&&(b_hang<0)
  //  A_frag  >>>>>>>>>>
  //  B_frag  >>>>>>>...
  //  becomes a proper to-contained overlap (ahg>0)&&(bhg==0)
  //  A_frag  <<<<<<<<<<
  //  B_frag     <<<<<<<

  //  This improper from-contained overlap (a_hang<0)&&(b_hang==0)
  //  A_frag  ...>>>>>>>
  //  B_frag  >>>>>>>>>>
  //  becomes a proper from-contained overlap (ahg==0)&&(bhg>0)
  //  A_frag  <<<<<<<
  //  B_frag  <<<<<<<<<<

  //  A degenerate overlap (a_hang==0)&&(b_hang==0)
  //  A_frag  >>>>>>>>>>
  //  B_frag  >>>>>>>>>>

#ifdef STORE_OVERLAP_EXTREMES
  const int iamn = iahg;
  const int iamx = iahg;
  const int ibmn = ibhg;
  const int ibmx = ibhg;
#endif
  const int iinv = FALSE;

  const int is_dovetail = is_a_dvt_simple(iahg,ibhg);
  const Tnes ines = ( is_dovetail
		      ? AS_CGB_DOVETAIL_EDGE
                      : AS_CGB_CONTAINED_EDGE );

  assert( (iahg>0) || (ibhg>0) || ((iahg == 0) && (ibhg == 0)) );
  
  the_new_raw_edge->avx = iafr;
  the_new_raw_edge->asx = iasx;
  the_new_raw_edge->ahg = iahg;
  
  the_new_raw_edge->bvx = ibfr;
  the_new_raw_edge->bsx = ibsx;
  the_new_raw_edge->bhg = ibhg;
  
  the_new_raw_edge->nes = ines;
  the_new_raw_edge->quality = qua;
  the_new_raw_edge->invalid = iinv;
  the_new_raw_edge->reflected = FALSE;
  the_new_raw_edge->grangered = FALSE;
  the_new_raw_edge->blessed = FALSE;

#ifdef STORE_OVERLAP_EXTREMES
  the_new_raw_edge->amn = iamn;
  the_new_raw_edge->amx = iamx;
  the_new_raw_edge->bmn = ibmn;
  the_new_raw_edge->bmx = ibmx;
#endif // STORE_OVERLAP_EXTREMES

  return;
}

static void add_Long_Olap_Data_t_to_graph
(
 Long_Olap_Data_t  * olap,
 Tfragment  frags[],
 Tedge      edges[],
 FragmentHashObject  *afr_to_avx, // Part of a hash table replacement
 TIntEdge_ID     next_edge_obj[],
 const int dvt_double_sided_threshold_fragment_end_degree,
 const int con_double_sided_threshold_fragment_end_degree,
 const int intrude_with_non_blessed_overlaps_flag,
 const CGB_ERATE_TYPE overlap_error_threshold,
 IntEdge_ID * novl_dovetail,
 IntEdge_ID * novl_containment,
 IntEdge_ID * nedges_delta
)
{
  Aedge  an_edge = {0};
  Convert_Olap_To_Aedge ( olap, & an_edge);
  if(
     ((AS_CGB_CONTAINED_EDGE == an_edge.nes) &&
      (is_a_frc_simple(an_edge.ahg,an_edge.bhg)))
      // Avoid entering the containment overlap twice.
      ){
    return;
  }
  if(overlap_error_threshold < an_edge.quality ) { return;}
  add_overlap_to_graph
    ( an_edge,
      frags,
      edges,
      afr_to_avx, // Part of a hash table replacement
      next_edge_obj,
      dvt_double_sided_threshold_fragment_end_degree,
      con_double_sided_threshold_fragment_end_degree,
      intrude_with_non_blessed_overlaps_flag,
      novl_dovetail,
      novl_containment,
      nedges_delta
      );
}

void InsertFragmentsIntoGraph
(
  THeapGlobals * heapva,
  VA_TYPE(OFGMesg)  * the_ofg_messages, // Additional fragments in a VA.
  VA_TYPE(char) * the_ofg_source,
  FragmentHashObject   * afr_to_avx,  // Part of a hash table replacement
  IntFragment_ID * Pnofg,
  IntFragment_ID * Pmin_frag_iid,
  IntFragment_ID * Pmax_frag_iid,
  VA_TYPE(char)  * frag_annotations
  )
{
  assert(NULL != heapva);
  assert(NULL != heapva->frags);
  assert(NULL != the_ofg_messages);
  assert(NULL != afr_to_avx);
  assert(NULL != Pnofg);
  assert(NULL != Pmin_frag_iid);
  assert(NULL != Pmax_frag_iid);
  //assert(NULL != frag_annotations);
  {
    // Add the fragments of a Variable Length Array into the fragment
    // overlap graph.
    const IntFragment_ID number_of_new_fragments = GetNumVA_OFGMesg(the_ofg_messages);
    const IntFragment_ID nfrag_base = GetNumFragments(heapva->frags);
    IntFragment_ID ii;

#if 1
    fprintf(stderr,"number_of_new_fragments = " F_IID " nfrag_base = " F_IID "\n",
            number_of_new_fragments, nfrag_base );
#endif    

    for( ii = 0; ii < number_of_new_fragments; ii++ ) {
      OFGMesg * ofg_mesg = GetVA_OFGMesg(the_ofg_messages, ii);
      assert(NULL != ofg_mesg);
#if 0
      fprintf(stderr,"  ii = " F_IID "\n", ii );
#endif    
      add_OFGMesg_to_graph
        ( ofg_mesg,
          the_ofg_source,
          heapva->frags,
          nfrag_base,
          afr_to_avx,
          Pnofg,
          Pmin_frag_iid,
          Pmax_frag_iid,
          frag_annotations
          );
    }
#if 0
    fprintf(stderr,"  ii = " F_IID "\n", ii );
#endif    
  }
  // return 0;
}

void InsertOverlapsIntoGraph
(
  THeapGlobals * heapva,
  VA_TYPE(OverlapMesg) * the_ovl_messages, // Additional overlaps in a VA.
  VA_TYPE(char)        * the_ovl_source,
  FragmentHashObject * afr_to_avx,        // Part of a hash table replacement
  const int dvt_double_sided_threshold_fragment_end_degree,
  const int con_double_sided_threshold_fragment_end_degree,
  const int intrude_with_non_blessed_overlaps_flag,
  const CGB_ERATE_TYPE overlap_error_threshold
  )
{
  assert(NULL != heapva);
  assert(NULL != heapva->edges);
  assert(NULL != the_ovl_messages);
  assert(NULL != afr_to_avx);
  {
    // Add the overlaps of a Variable Length Array into the fragment
    // overlap graph.
    
    const IntEdge_ID number_of_new_overlaps = GetNumVA_OverlapMesg(the_ovl_messages);
    IntEdge_ID novl_dovetail = 0;
    IntEdge_ID novl_containment = 0;
    IntEdge_ID novl_degenerate = 0;
    IntEdge_ID nedges_delta = 0;
    IntEdge_ID ii;
    
    for( ii = 0; ii < number_of_new_overlaps; ii++ ) {
      OverlapMesg * ovl_mesg = GetVA_OverlapMesg(the_ovl_messages, ii);
      assert(NULL != ovl_mesg);
      add_OverlapMesg_to_graph
        ( ovl_mesg,
          the_ovl_source,
          heapva->frags,
          heapva->edges,
          afr_to_avx,  // Part of a hash table replacement
          heapva->next_edge_obj,
          dvt_double_sided_threshold_fragment_end_degree,
          con_double_sided_threshold_fragment_end_degree,
          intrude_with_non_blessed_overlaps_flag,
          overlap_error_threshold,
          &novl_dovetail,
          &novl_containment,
          &novl_degenerate,
          &nedges_delta );
    }
    fprintf(stderr,"Added %10" F_IIDP " OVL records.\n",novl_dovetail+novl_containment);
    fprintf(stderr,"      %10" F_IIDP " OVL dovetail records.\n",novl_dovetail);
    fprintf(stderr,"      %10" F_IIDP " OVL containment records.\n",novl_containment);
    fprintf(stderr,"Found %10" F_IIDP " OVL records with no valid frag (degenerate).\n",novl_degenerate);
  }
  //return 0;
}

/****************************************************************************/

void process_ovl_store
(
 char * OVL_Store_Path,
 Tfragment  frags[],
 // The internal representation of the fragment reads. 
 Tedge      edges[],
 // The internal representation of the overlaps.
 VA_TYPE(char)   frag_annotations[],
 FragmentHashObject *afr_to_avx, // Part of a hash table replacement
 IntFragment_ID *min_frag_iid,
 IntFragment_ID *max_frag_iid,
 TIntEdge_ID     next_edge_obj[],
 BPTYPE         *nbase_in_genome,
 const int dvt_double_sided_threshold_fragment_end_degree,
 const int con_double_sided_threshold_fragment_end_degree,
 const int intrude_with_non_blessed_overlaps_flag,
 const CGB_ERATE_TYPE overlap_error_threshold
 ) 
{
 Long_Olap_Data_t  olap;
 IntFragment_ID last=0;
 IntFragment_ID Start = (*min_frag_iid);
 IntFragment_ID Stop  = (*max_frag_iid);
 
 IntEdge_ID novl_dovetail = 0;
 IntEdge_ID novl_containment = 0;
 IntEdge_ID nedges_delta = 0;

 OVL_Store_t  * my_store = New_OVL_Store ();
 OVL_Stream_t  * my_stream = New_OVL_Stream ();

 Open_OVL_Store (my_store, OVL_Store_Path);

 last = Last_Frag_In_OVL_Store (my_store);
 Start = 1;
 Stop  = last;
 
 if(Stop > last) { Stop = last;}
 
 fprintf(stderr," Start = " F_IID " Stop =" F_IID "\n",  Start, Stop);
 Init_OVL_Stream (my_stream, Start, Stop, my_store);
 while  (Next_From_OVL_Stream (& olap, my_stream))
   {
     add_Long_Olap_Data_t_to_graph
       ( &olap,
         frags, edges, afr_to_avx, next_edge_obj,
         dvt_double_sided_threshold_fragment_end_degree,
         con_double_sided_threshold_fragment_end_degree,
         intrude_with_non_blessed_overlaps_flag,
         overlap_error_threshold,
         &novl_dovetail,
         &novl_containment,
         &nedges_delta
         );
   }
 
 Free_OVL_Stream (my_stream);
 Free_OVL_Store (my_store);

 fprintf(stderr,"novl_dovetail    = " F_IID "\n", novl_dovetail);
 fprintf(stderr,"novl_containment = " F_IID "\n", novl_containment);
 fprintf(stderr,"nedges_delta     = " F_IID "\n", nedges_delta);
}

/****************************************************************************/
/****************************************************************************/

static void input_mesgs_internal
(int             argc, 
 char           *argv[],
 FILE           *fovl, 
 FILE           *filk,
 int            *Pnadt,
 int            *Pnidt, 
 int            *Pnilk,
 IntFragment_ID *Pnofg,
 IntEdge_ID     *Pnedges_delta,
 IntEdge_ID     *Pnovl_dovetail, 
 IntEdge_ID     *Pnovl_containment,
 IntEdge_ID     *Pnovl_degenerate,
 BPTYPE         *nbase_in_genome,
 IntFragment_ID  nfrag_base,
 Tfragment       frags[],
 IntEdge_ID      nedge_base,
 Tedge           edges[],
 VA_TYPE(char)   frag_annotations[],
 FragmentHashObject *afr_to_avx, // Part of a hash table replacement
 IntFragment_ID *Pmin_frag_iid,
 IntFragment_ID *Pmax_frag_iid,
 TIntEdge_ID     next_edge_obj[],
 const int dvt_double_sided_threshold_fragment_end_degree,
 const int con_double_sided_threshold_fragment_end_degree,
 const int intrude_with_non_blessed_overlaps_flag,
 const CGB_ERATE_TYPE overlap_error_threshold
)
{
  /* It is assumed that in the overlap records that new fragments
    point to old fragments.  */

  int nadt=0,nidt=0,nilk=0;
  IntEdge_ID nedges_delta=0,novl_dovetail=0,novl_containment=0,novl_degenerate=0;
  GenericMesg *pmesg;
  VA_TYPE(char) * the_ofg_source = NULL;
  VA_TYPE(char) * the_ovl_source = NULL;

  while( EOF != ReadProtoMesg_AS(fovl, &pmesg)) {
    const MessageType imesgtype = pmesg->t;

    // printf("XXX pmesg->t = %d\n", imesgtype);
    switch(imesgtype) {
    case MESG_ADT: 
      {
	AuditMesg  *adt_mesg = (AuditMesg  *)(pmesg->m);
	AuditLine  *adl;
#if 0
        AuditLine auditLine;
#endif
	/*  Audit Message record-- for now scavenge for genome length */
	nadt++;
	if(NULL != filk) {
          static const char teststr[] = "Genome Length is";
	  for(adl = adt_mesg->list; adl != NULL; adl = adl->next) {
	    const size_t len1 = strlen(teststr);
	    const size_t len2 = strlen(adl->comment);
	    const size_t len3 = (len1 < len2 ? len1 : len2);
	    if(0 == strncmp(teststr,adl->comment,len3))
	      { 
		sscanf(  &((adl->comment)[strlen(teststr)]), BPFORMAT,
			 nbase_in_genome);}
#if 0
	    fprintf(stderr,
		    "AS_FGB_io.c: nbase_in_genome=" BPFORMAT "\n",
		    (*nbase_in_genome));
#endif          
	  }
#if 1
#if 0
	  AppendAuditLine_AS(adt_mesg, &auditLine, time(NULL), "CGB", 
			     CM_ID, "(empty)");
#else
	  VersionStampADT(adt_mesg, argc, argv);
#endif
#else
	  VersionStampADTWithCommentAndVersion
	    (adt_mesg, argc, argv,
	     (GlobalParamText != NULL
	      ? GlobalParamText:""),"(blank)");
	  AppendAuditLine_AS
	    (adt_mesg, &auditLine, time(NULL),
	     argv[0], "",params);
	  
#endif        
	  WriteProtoMesg_AS(filk,pmesg);
	}
      }
      break;
    case MESG_IDT: 
      {
	InternalDistMesg  *idt_mesg;
	/*  Distance record--skip for now */
	idt_mesg = (InternalDistMesg  *)(pmesg->m);
	nidt++;
	if(NULL != filk) {
	  WriteProtoMesg_AS(filk,pmesg);
	}
      }
      break;
    case MESG_ILK: 
      {
	nilk++;
	if(NULL != filk) {
	  WriteProtoMesg_AS(filk,pmesg);
	}
      }
      break;
    case MESG_OFG: 
      {
        add_OFGMesg_to_graph
          ( (OFGMesg *) pmesg->m,
            the_ofg_source,
            frags,
            nfrag_base,
            afr_to_avx,  // part of a hash table replacement
            Pnofg,
            Pmin_frag_iid,
            Pmax_frag_iid,
            frag_annotations
            );
      }
      break;
    case MESG_OVL:
      {
        add_OverlapMesg_to_graph
          ( (OverlapMesg *) pmesg->m,
            the_ovl_source,
            frags, edges, afr_to_avx, next_edge_obj,
            dvt_double_sided_threshold_fragment_end_degree,
            con_double_sided_threshold_fragment_end_degree,
            intrude_with_non_blessed_overlaps_flag,
            overlap_error_threshold,
            &novl_dovetail,
            &novl_containment,
            &novl_degenerate,
            &nedges_delta
            );
      }
      break;
    case MESG_IBC:
      {
	// Just swallow these...
	if(NULL != filk) {
	  // WriteProtoMesg_AS(filk,pmesg);
	}
      }
      break;
    case MESG_IBA:
      {
	if(NULL != filk) {
	  WriteProtoMesg_AS(filk,pmesg);
	}
      }
      break;
    default:
      {
	fprintf(stderr,"Unexpected message type %d\n",imesgtype);
	fprintf(stderr,"Message typename %s\n",
		MessageTypeName[imesgtype]);
	assert(FALSE);
      }
      break;
    }
  }
  *Pnadt = nadt;
  *Pnidt = nidt;
  *Pnilk = nilk;
  *Pnedges_delta = nedges_delta;
  *Pnovl_dovetail = novl_dovetail;
  *Pnovl_containment = novl_containment;
  *Pnovl_degenerate = novl_degenerate;
}


void input_messages_from_a_file
(int        argc, 
 char       *argv[],
 FILE       *fovl,
 FILE       *filk,
 Tfragment  frags[],
 // The internal representation of the fragment reads. 
 Tedge      edges[],
 // The internal representation of the overlaps.
 VA_TYPE(char) frag_annotations[],
 FragmentHashObject *afr_to_avx, // Part of a hash table replacement
 IntFragment_ID    *min_frag_iid,
 IntFragment_ID    *max_frag_iid,
 TIntEdge_ID       *next_edge_obj,
 BPTYPE *nbase_in_genome,
 const int dvt_double_sided_threshold_fragment_end_degree,
 const int con_double_sided_threshold_fragment_end_degree,
 const int intrude_with_non_blessed_overlaps_flag,
 const CGB_ERATE_TYPE overlap_error_threshold
 ) 
{
  /* Input a batch of fragment read and overlap records from a stream. */
  /* Keep a copy of the number of fragments and edges before
     the new data is read in. */
  const IntFragment_ID nfrag_old = GetNumFragments(frags); 
  const IntEdge_ID nedge_old = GetNumEdges(edges);
  
  IntFragment_ID nfrag_new = nfrag_old; // Note that this is NOT const.
  IntEdge_ID nedge_new = nedge_old;
  
  int nadt=0;   /* The number of audit messages read. */
  int nidt=0;   /* The number of distance records read. */
  int nilk=0;   /* The number of internal link messages. */
  IntFragment_ID nofg=0;   /* The number of fragment records read. */
  IntEdge_ID novl_dovetail=0,novl_containment=0,novl_degenerate=0; /* The number of overlap records read. */
  IntEdge_ID nedge_delta=0;   
  
  time_t tp1,tp2; // int32 seconds from the beginning of 1970.

  if(TIMINGS) { 
    time(&tp1); fprintf(stderr,"Begin input of OVL file\n");
    system_top();
  }

  if(nfrag_old > 0){
    const IntFragment_ID iid = get_iid_fragment(frags,0);
    (*max_frag_iid) = MAX((*max_frag_iid),iid);
    (*min_frag_iid) = MIN((*min_frag_iid),iid);
  }
  if(nfrag_old > 1) {
    IntFragment_ID vid = 0;
    for(vid=1;vid<nfrag_old;vid++) {
      const IntFragment_ID iid = get_iid_fragment(frags,vid);
      (*max_frag_iid) = MAX((*max_frag_iid),iid);
      (*min_frag_iid) = MIN((*min_frag_iid),iid);
    }
  }
  assert((*min_frag_iid) <= (*max_frag_iid));
  
  
#ifdef DEBUG21
  {
    IntFragment_ID iv0;
    printf("All fragments before input_messages\n");
    printf("nfrag_old=" F_IID "\n",nfrag_old);
    for(iv0=0; iv0<nfrag_old; iv0++) {
      const IntFragment_ID iid = get_iid_fragment(frags,iv0);

      printf("%5" F_IIDP ": iid " F_IID ", src %10" F_SIZE_TP ", "
	     "prefix %8" F_U32P ", %5" F_S32P ", suffix %8" F_S32P ", %5 "F_S32P " : %5d \n",
	     iv0,
	     get_iid_fragment(frags,iv0),
	     get_src_fragment(frags,iv0),
	     get_segstart_vertex(frags,iv0,FALSE),
	     get_seglen_vertex(frags,iv0,FALSE),
	     get_seglen_frc_vertex(frags,iv0,FALSE),
	     get_segstart_vertex(frags,iv0,TRUE),
	     get_seglen_vertex(frags,iv0,TRUE),
	     get_seglen_frc_vertex(frags,iv0,TRUE),
	     get_lab_fragment(frags,iv0)
	     );
    }
    printf("min_frag_iid=" F_IID " max_frag_iid=" F_IID "\n",
	   (*min_frag_iid),(*max_frag_iid));
  }
#endif /*DEBUG21*/
  input_mesgs_internal
    (argc,argv,
     fovl,filk,
     &nadt,&nidt,&nilk,&nofg,&nedge_delta,
     &novl_dovetail,&novl_containment,&novl_degenerate,
     nbase_in_genome,
     nfrag_old, frags,
     nedge_old, edges,
     frag_annotations,
     afr_to_avx,  // Part of a hash table replacement
     min_frag_iid,max_frag_iid,
     next_edge_obj,
     dvt_double_sided_threshold_fragment_end_degree,
     con_double_sided_threshold_fragment_end_degree,
     intrude_with_non_blessed_overlaps_flag,
     overlap_error_threshold
     );

  fprintf(stderr,"Input %10d ADT records.\n",nadt);
  fprintf(stderr,"Input %10d IDT records.\n",nidt);
  fprintf(stderr,"Input %10d ILK records.\n",nilk);
  fprintf(stderr,"Input %10d OFG records.\n",nofg);
  fprintf(stderr,"Input %10" F_IIDP " OVL records (skipped %10"F_IIDP" degenerate).\n",novl_dovetail+novl_containment, novl_degenerate);
  fprintf(stderr,"min_frag_iid=" F_IID " max_frag_iid=" F_IID "\n",
	  (*min_frag_iid),(*max_frag_iid));

  fprintf(stderr,"      %10" F_IIDP " OVL dovetail records.\n",novl_dovetail);
  fprintf(stderr,"      %10" F_IIDP " OVL containment records.\n",novl_containment);

  if(TIMINGS) {
    time(&tp2); 
    fprintf(stderr,"%10" F_TIME_TP " sec: Finished reading input file\n",
	    (tp2-tp1));
    system_top();
  }

  fprintf(stderr,"nfrag_old=" F_IID ",nofg=" F_IID ",nedge_old=" F_IID ",nedge_delta=" F_IID "\n",
          nfrag_old,nofg,nedge_old,nedge_delta);

  nfrag_new = nfrag_old + nofg;
  assert(nfrag_new == GetNumFragments(frags));
  
  if(TIMINGS) { 
    time(&tp1); fprintf(stderr,"Begin input phase check\n");
    system_top();
  }
#ifdef DEBUG_INPUT
  {
    FILE *fchk;
    char strtmp[1024]={0};
    strcpy(strtmp,File_Prefix);
    strcat(strtmp,".chk_raw");
    fprintf(stderr,"Openning %s for chunk info.\n",strtmp);
    if((fchk = fopen(strtmp,"w")) == NULL)
      assert(0);
    fprintf(stderr,"Openned %s for chunk info.\n",strtmp);
    output_chk_raw(fchk,nidt,nofg,novl,frags,edges);
    fclose(fchk);
  }
#endif /*DEBUG_INPUT*/
  

  if(TIMINGS) {
    time(&tp2); 
    fprintf(stderr,"%10" F_TIME_TP " sec: Finished input phase check\n",(tp2-tp1));
    system_top();
  }

  nedge_new = nedge_old + nedge_delta;
  assert(nedge_new == GetNumEdges(edges));

  fprintf(stderr,
          "AS_FGB_io.c: nbase_in_genome=" BPFORMAT "\n",
          (*nbase_in_genome));
}

