
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
= "$Id: AS_CGB_main.c,v 1.4 2005-03-22 19:48:29 jason_miller Exp $";
/*********************************************************************
 *
 * Module:  AS_CGB_main.c
 * Description:
 *    Unitigger main
 * 
 *    Reference: ChunkGraphBuilder.rtf
 *
 *    Command Line Interface:
 *        $ cgb
 *
 * -p <parameterPath>    If specified, parameters are loaded from this path
 * -X Enable developer's options.
 *    Any option not stated in the tool's SOP should require this flag
 * -P Output is ASCII (default is binary)
 * -D <level>  Specify debug level (only valid with -X)
 * -v <level>  Specify verbosity level (only valid with -X)
 * -A <level> Run the graph analyzer
 * -E <int>  The number of errors allowed.
 * -C <string> The check-pointing information.
 * -R <string> The restart information.
 *
 *        -l the length of the genome
 *        -C Use a consensus sequence to compute the branch points
 *        -R <int> Which global containment rule to use.
 *        -W <int> Limit of search depth for fragment graph walking in 
 *           transitively inferable overlap removal.
 *        -Y Do not bother to count chimeras
 *
 *       FragStorePath:   Used to find/create a directory that contains the
 *       fragment store 
 *
 *       GraphStorePath:  Used to read/write a fragment graph store.
 *      
 *       OverlapInputFile: The file with new OVL records to process. 
 *
 *       Unitigger processes the input file, reporting errors and
 *       warnings.  Any message that causes an error warning is output
 *       to stderr, along with the associated warning/error messages.
 *       The fragment graph stage produces a <inputFile>.fgb file 
 *       and the UnitiggerStore is updated. The chunk graph builder stage
 *       reads the UnitiggerStore and produces a <inputFile>.cgb 
 *       In addition an optional <inputFile>.cga file is produced
 *       with diagnostic information about the fragment overlap graph
 *       and the unitig graph.
 *
 *       See ChunkGraphBuilder.rtf for descriptions of the checks
 *       performed.
 *
 *    Programmer:  C. M. Mobarry
 *       Written:  Jan 1999
 * 
 *********************************************************************/

/*********************************************************************/

/* Local include files */
#include "AS_UTL_version.h"  
#include "AS_CGB_all.h"
#include "AS_CGB_unitigger_globals.h"

/****************************************************************************/
#define DEBUGGING
#define USE_IBA_FILE
#undef LOAD_FRAGSTORE_INTO_MEMORY
#undef SWITCH_CONTAINMENT_DIRECTION_CGB
#undef NOT_PACKED
#undef DEBUG_RISM

/* These are optimization parameters. */

/****************************************************************************/
/* Globals */

static int TIMINGS = TRUE;
extern int REAPER_VALIDATION;

/****************************************************************************/

static void reflect_an_edge_in_place
(
 Tedge * edges,
 IntEdge_ID iedge
) 
{
  Aedge * the_edge = GetVA_Aedge(edges,iedge);
  reflect_Aedge(the_edge,the_edge);
}

static IntEdge_ID get_the_thickest_dvt_overlap_from_vertex
( Tfragment * frags,
  Tedge     * edges,
  TIntEdge_ID * next_edge_obj, // If non-NULL do not assume an adjaceny representation.
  IntFragment_ID avx,
  int            asx
  ) {
  // Does the graph have to be in adjacency format?
  const IntEdge_ID nedge = GetNumEdges(edges);
  IntEdge_ID ie_thickest = AS_CGB_EDGE_NOT_FOUND;
  
  if(GetNumVA_IntEdge_ID(next_edge_obj) == 0) {
    // Assume adjacency mode
    const IntEdge_ID segstart = get_segstart_vertex(frags,avx,asx);
    const int        seglen   = get_seglen_vertex(frags,avx,asx);
    IntEdge_ID ie;
    for( ie=segstart; ie < segstart+seglen; ie++) {
      const Tnes nes = get_nes_edge(edges,ie);
      assert(
             (AS_CGB_DOVETAIL_EDGE == nes) || 
             (AS_CGB_THICKEST_EDGE == nes) || 
             (AS_CGB_INTERCHUNK_EDGE == nes) ||
             (AS_CGB_TOUCHES_CRAPPY_DVT == nes) ||
             (AS_CGB_BETWEEN_CRAPPY_DVT == nes) ||
             (AS_CGB_CONTAINED_EDGE == nes) ||
             (AS_CGB_TOUCHES_CRAPPY_CON == nes) ||
             (AS_CGB_BETWEEN_CRAPPY_CON == nes) );
      if(
         (AS_CGB_DOVETAIL_EDGE == nes) || 
         (AS_CGB_THICKEST_EDGE == nes) || 
         (AS_CGB_INTERCHUNK_EDGE == nes) 
         // get_dvt_edge(edges,ie)
         ) {
        ie_thickest = ie;
        break;
      }
    }
  } else {
    // Assume hash mode
    const int degree = get_seglen_dvt_vertex(frags,avx,asx);
    assert(degree <= 1);
    assert(FALSE);
    // assert that "unitigger -x 1" was used.
    if(degree == 1 ) { 
      ie_thickest = get_segend_vertex(frags,avx,asx);
    }
  }

  assert((AS_CGB_EDGE_NOT_FOUND == ie_thickest) || (nedge > ie_thickest));

#ifdef ZIGGY
  if(nedge > ie_thickest) {
    IntFragment_ID bvx = get_bvx_edge(edges,ie_thickest);
    int            bsx = get_bsx_edge(edges,ie_thickest);
    int           acon = get_con_fragment(frags,avx);
    int           bcon = get_con_fragment(frags,bvx);
    Tlab          alab = get_lab_fragment(frags,avx);
    Tlab          blab = get_lab_fragment(frags,bvx);
    fprintf(stdout,"ZIGGY path "
            "afr= " F_IID " asx= %d acon= %d alab= %d bfr= " F_IID " bsx= %d bcon= %d blab= %d\n",
            get_iid_fragment(frags, avx), asx, acon, alab, 
            get_iid_fragment(frags, bvx), bsx, bcon, blab );
  }
#endif // ZIGGY
  
  return ie_thickest;
}

static void identify_thickest_overlaps
( Tfragment * frags,
  Tedge     * edges,
  TIntEdge_ID * next_edge_obj 
  ) 
{
  const IntEdge_ID nedge = GetNumEdges(edges);
  const IntFragment_ID nfrag = GetNumFragments(frags);
  
  if(GetNumVA_IntEdge_ID(next_edge_obj) > 0) {
    // Assume hash mode
    IntEdge_ID ie;
    assert(FALSE);
    for(ie=0; ie < nedge; ie++) {
      const Tnes nes = get_nes_edge(edges,ie);
      if((AS_CGB_DOVETAIL_EDGE == nes) ||
         (AS_CGB_THICKEST_EDGE == nes)
         ) {
        const IntFragment_ID avx = get_avx_edge(edges,ie);
        const IntFragment_ID asx = get_asx_edge(edges,ie);
        const int degree = get_seglen_dvt_vertex(frags,avx,asx);
        if( degree != 1 ) {
          const IntFragment_ID bvx = get_bvx_edge(edges,ie);
          const IntFragment_ID bsx = get_bsx_edge(edges,ie);
          fprintf(stderr,"ERROR: afr= " F_IID " asx= %d bfr= " F_IID " bsx= %d nes= %d degree= %d\n",
                  get_iid_fragment(frags,avx), asx,
                  get_iid_fragment(frags,bvx), bsx,
                  nes, degree);
        }
        assert(degree == 1); // assert that "unitigger -x 1" was used.
        if((AS_CGB_DOVETAIL_EDGE == nes)) {
          set_nes_edge(edges,ie,AS_CGB_THICKEST_EDGE);
          fix_overlap_edge_mate(frags,edges,ie);
        }
      }
    }
  } else {
    // Assume an adjaceny representation.
    IntFragment_ID avx;
    int asx;
    for(avx=0; avx < nfrag; avx++) {
      for(asx=0; asx < 2; asx++) {
        const IntEdge_ID ie = get_the_thickest_dvt_overlap_from_vertex
          ( frags, edges, next_edge_obj, // If non-NULL do not assume an adjacency representation.
            avx, asx );
        if(AS_CGB_EDGE_NOT_FOUND != ie) {
          const Tnes nes = get_nes_edge(edges,ie);
          //assert((AS_CGB_DOVETAIL_EDGE == nes) || (AS_CGB_THICKEST_EDGE == nes));
          if((AS_CGB_DOVETAIL_EDGE == nes)) {
            set_nes_edge(edges,ie,AS_CGB_THICKEST_EDGE);
            fix_overlap_edge_mate(frags,edges,ie);
          }
        }
      }
    }
  }
}

static void identify_bmpc_paths
( Tfragment * frags,
  Tedge     * edges,
  TIntEdge_ID * next_edge_obj, // If non-NULL do not assume an adjaceny representation.
  IntEdge_ID ie_start
  ) {

  // This BMPC path marking implementation restricted to following the
  // thickest overlap graph edges connecting contained fragments. The
  // index "ie_start" is assumed to be an outgoing non-chordal
  // dovetail overlap from a known essential fragment.

  // The side effects are: (1) the starting edge is marked essential.
  // (2) that every thickest edge visited will be marked essential.
  // (3) every contained fragment visited will be marked essential.
  
  IntEdge_ID ie_now;
  IntEdge_ID ie_next;
  IntEdge_ID nedge = GetNumEdges(edges);
  
  for( ie_now = ie_start;
       (ie_now < nedge) && (AS_CGB_THICKEST_EDGE == get_nes_edge(edges,ie_now));
       // Note that nedge < AS_CGB_EDGE_NOT_FOUND.
       ie_now = ie_next
       ) {
    const IntFragment_ID bvx = get_bvx_edge(edges,ie_now);
    const int bcon = get_con_fragment(frags,bvx);
    set_nes_edge(edges,ie_now,AS_CGB_INTERCHUNK_EDGE);
    fix_overlap_edge_mate(frags,edges,ie_now);

    if(bcon) {
      const int   bsx = get_bsx_edge(edges,ie_now);
      const Tlab  lab = get_lab_fragment(frags,bvx);
      assert((AS_CGB_UNPLACEDCONT_FRAG == lab) ||
             (AS_CGB_BRANCHMULTICONT_FRAG == lab));
      // What if this fragment was labeled a spur?
      set_lab_fragment( frags, bvx, AS_CGB_BRANCHMULTICONT_FRAG);
#ifdef ZIGGY
      {
        IntFragment_ID avx = get_avx_edge(edges,ie_now);
        int            asx = get_asx_edge(edges,ie_now);
        int           acon = get_con_fragment(frags,avx);
        int           bcon = get_con_fragment(frags,bvx);
        Tlab          alab = get_lab_fragment(frags,avx);
        Tlab          blab = get_lab_fragment(frags,bvx);
        fprintf(stdout,"ZIGGY start "
                "afr= " F_IID " asx= %d acon= %d alab= %d bfr= " F_IID " bsx= %d bcon= %d blab= %d\n",
                get_iid_fragment(frags, avx), asx, acon, alab, 
                get_iid_fragment(frags, bvx), bsx, bcon, blab );
      }
#endif // ZIGGY

      ie_next =
        get_the_thickest_dvt_overlap_from_vertex
        (frags,edges,next_edge_obj,bvx,!bsx);
      // Get the thickest dvt overlap from the fragment's other end.

      if( AS_CGB_EDGE_NOT_FOUND != ie_next ) {
        Tnes nes = get_nes_edge(edges,ie_now);
        assert( (AS_CGB_INTERCHUNK_EDGE == nes) || (AS_CGB_THICKEST_EDGE == nes) );
        assert(bvx == get_avx_edge(edges,ie_next));
        assert(!bsx == get_asx_edge(edges,ie_next));
      }
      
    } else {
      // ie_next = nedge; // terminate
      break;
    }
  }
}

static void identify_essential_components
( Tfragment * frags,
  Tedge     * edges,
  TIntEdge_ID * next_edge_obj 
  ) {
  const IntEdge_ID nfrag = GetNumFragments(frags);

  { // Assume an adjaceny representation.
    IntFragment_ID avx;
    int asx;
    for(avx=0; avx < nfrag; avx++) {
      const IntFragment_ID con = get_con_fragment(frags,avx);
      const IntFragment_ID lab = get_lab_fragment(frags,avx);
      if( (!con) && ((AS_CGB_HANGING_FRAG == lab) || (AS_CGB_THRU_FRAG == lab)) ) {
        for(asx=0; asx < 2; asx++) {
          const IntEdge_ID ie = get_the_thickest_dvt_overlap_from_vertex
            ( frags, edges, next_edge_obj, // If non-NULL do not assume an adjaceny representation.
              avx, asx );
          if(AS_CGB_EDGE_NOT_FOUND != ie) {
            const Tnes nes = get_nes_edge(edges,ie);
            assert((AS_CGB_THICKEST_EDGE == nes) || (AS_CGB_INTERCHUNK_EDGE == nes));
            if(REAPER_VALIDATION) {
              if((AS_CGB_THICKEST_EDGE == nes)) {
                set_nes_edge(edges,ie,AS_CGB_INTERCHUNK_EDGE);
                fix_overlap_edge_mate(frags,edges,ie);
              }
            } else {
              identify_bmpc_paths( frags, edges, next_edge_obj, ie);
            }
          }
        }
      }
    }
  }
}


/****************************************************************************/


#ifdef  SWITCH_CONTAINMENT_DIRECTION_CGB

static void reflect_containment_direction_in_place
(
 Tedge * edges,
 int     become_to_contained
) 
{
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
  const IntEdge_ID nedge = GetNumEdges(edges);
  IntEdge_ID iedge;
  time_t tp1, tp2;
  
  if(TIMINGS) {
    tp1 = time(NULL); fprintf(stderr,"Begin reflect_containment_direction_in_place\n");
    system_top();
  }

  for(iedge=0; iedge<nedge; iedge++) {
    const Tnes nes = get_nes_edge(edges, iedge);
    const int ahg = get_ahg_edge(edges, iedge);
    const int bhg = get_bhg_edge(edges, iedge);
    if((AS_CGB_CONTAINED_EDGE == nes)&&
       (
        (become_to_contained &&
         (is_a_frc_edge(edges,iedge) ||
          (is_a_dgn_edge(edges,iedge)&&(get_avx_edge(edges,iedge) > get_bvx_edge(edges,iedge))))
         )
        ||
        (!become_to_contained &&
         (is_a_toc_edge(edges,iedge) ||
          (is_a_dgn_edge(edges,iedge)&&(get_avx_edge(edges,iedge) < get_bvx_edge(edges,iedge))))
         )
        )
       ) {
      reflect_an_edge_in_place( edges, iedge);
    }
  }
  if(TIMINGS) { 
    tp2 = time(NULL); 
    fprintf(stderr,"%10" F_TIME_TP " sec: Finished reflect_containment_direction_in_place\n",(tp2-tp1));
    system_top();
  }
}

static void convert_all_containment_overlaps_direction_TOC
( Tfragment * frags,
  Tedge     * edges,
  TIntEdge_ID * next_edge_obj, // If non-NULL do not assume an adjaceny representation.
  const int become_to_contained
 )
{
  time_t tp1, tp2;
  if(TIMINGS) {
    tp1 = time(NULL);
    fprintf(stderr,"Begin switching the containment direction in graph store.\n");
    system_date();
  }
  reflect_containment_direction_in_place( edges, become_to_contained);
  reorder_edges( frags, edges, next_edge_obj);
  if(TIMINGS) {
    tp2 = time(NULL); 
    fprintf(stderr,"%10" F_TIME_TP " sec: Finished switching the containment direction in graph store.\n",
            (tp2-tp1));
    system_date();
  }
}
#endif // SWITCH_CONTAINMENT_DIRECTION_CGB

static void maskout_overlaps_touching_crappy_fragments
( Tfragment * frags,
  Tedge     *edges,
  TIntEdge_ID * next_edge_obj // If non-NULL do not assume an adjaceny representation.
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
    case AS_CGB_DOVETAIL_EDGE:
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
#ifdef NOT_PACKED
    case AS_CGB_REMOVED_BY_TRANSITIVITY_DVT:
    case AS_CGB_REMOVED_BY_DUPLICATE_DVT:
    case AS_CGB_REMOVED_BY_TRANSITIVITY_CON:
    case AS_CGB_REMOVED_BY_DUPLICATE_CON:
      // Do nothing ....
      break;	
#endif      
    default:
      fprintf(stderr,"Unexpected overlap label ines=%d\n",ines);
      assert(FALSE);
      break;
    }
  }}
}


int main_cgb
(
 int argc,
 char * argv [],
 TStateGlobals * gstate,
 THeapGlobals  * heapva,
 UnitiggerGlobals * rg
 )
{
  FragStoreHandle TheFragStore = NULLFRAGSTOREHANDLE;

  int status = 0; // Main program exit status.
  int ierr = 0;

  time_t tp1,tp2; 
  //The time() function returns the time in seconds since the
  //Epoch. The Epoch is referenced to 00:00:00 CUT (Coordinated
  //Universal Time) 1 Jan 1970.

#ifdef _OPENMP
  fprintf(stderr,"Compiled by an OpenMP compilant compiler\n");
#endif

#ifdef DEBUGGING
  fprintf(stderr,__FILE__ " compiled at "  __DATE__ " " __TIME__ "\n");
  fprintf(stderr,"%s\n",CM_ID);
#endif

  system_date();

  
  if( NULL != rg->frag_store ) { 
    char *theFragStorePath = rg->frag_store;
    // Used to squirrel away a clean copy of the store.
  fprintf(stderr, __FILE__ " theFragStorePath = <%s>\n", theFragStorePath);
#ifndef LOAD_FRAGSTORE_INTO_MEMORY
    // Open the fragstore, readonly
    TheFragStore = openFragStore(theFragStorePath,"r");
#else // LOAD_FRAGSTORE_INTO_MEMORY
    TheFragStore = loadFragStore(theFragStorePath);
#endif // LOAD_FRAGSTORE_INTO_MEMORY
    assert(TheFragStore != NULLFRAGSTOREHANDLE);
  }
  
  /* The remaining command line arguments are input Assembler
     message files. */
  
  count_fragment_and_edge_labels( heapva->frags, heapva->edges,
                                  "After reading the fragment graph store");
  check_symmetry_of_the_edge_mates( heapva->frags, heapva->edges, heapva->next_edge_obj);

  if(rg->aggressive_spur_fragment_marking) {
    maskout_overlaps_touching_crappy_fragments
      ( heapva->frags, heapva->edges, heapva->next_edge_obj );
    count_fragment_and_edge_labels
      ( heapva->frags, heapva->edges, "after maskout_overlaps_touching_crappy_fragments");
    // Do we need to patch up the graph here??
#ifdef DEBUG_RISM
    view_fgb_chkpnt( "maskout_overlaps_touching_crappy_fragments", heapva->frags, heapva->edges);
#endif // DEBUG_RISM
    check_symmetry_of_the_edge_mates( heapva->frags, heapva->edges, heapva->next_edge_obj);
  }


  ///////////////////////////////////////////////////////////////////////////////////
  
  identify_thickest_overlaps( heapva->frags, heapva->edges, heapva->next_edge_obj );
  count_fragment_and_edge_labels( heapva->frags, heapva->edges,
                                  "After identify_thickest_edges");
#ifdef DEBUG_RISM
  view_fgb_chkpnt( "identify_thickest_overlaps", heapva->frags, heapva->edges);
  check_symmetry_of_the_edge_mates( heapva->frags, heapva->edges, heapva->next_edge_obj);
#endif    

  identify_essential_components( heapva->frags, heapva->edges, heapva->next_edge_obj );
  count_fragment_and_edge_labels( heapva->frags, heapva->edges,
                                  "After identify_essential_components");
#ifdef DEBUG_RISM
  view_fgb_chkpnt( "identify_essential_components", heapva->frags, heapva->edges);
  check_symmetry_of_the_edge_mates( heapva->frags, heapva->edges, heapva->next_edge_obj);
#endif    

#if 0
  identify_buddy_overlaps(heapva->frags,heapva->edges);
  count_fragment_and_edge_labels
    ( heapva->frags, heapva->edges,
      "After identify_buddy_overlaps");
#endif
  
  chunk_classification_dvt
    ( heapva->frags, heapva->edges,
      rg->walk_depth,
      rg->remove_blizzard_overlaps
      );
  
  count_fragment_and_edge_labels( heapva->frags, heapva->edges,
                                  "In main after chunk_classification_dvt");
  check_symmetry_of_the_edge_mates( heapva->frags, heapva->edges, heapva->next_edge_obj);


  //////////////////////////////////////////////////////////////////////
    
#ifdef SWITCH_CONTAINMENT_DIRECTION_CGB
  convert_all_containment_overlaps_direction_TOC
    ( heapva->frags, heapva->edges, heapva->next_edge_obj, TRUE);
  count_fragment_and_edge_labels( heapva->frags, heapva->edges,
                                  "After switching the containment direction");
#endif // SWITCH_CONTAINMENT_DIRECTION_CGB


  { // beginning of the non-incremental phase
    
    const IntFragment_ID nfrag = GetNumFragments(heapva->frags);

    {
      /* Count the total amount of guide fragments. */
      IntFragment_ID ifrag;
      IntFragment_ID num_of_guides_total = 0;
      for(ifrag=0;ifrag<nfrag;ifrag++) { 
	const FragType type = get_typ_fragment(heapva->frags,ifrag);
	if((type != AS_READ) && 
	   (type != AS_EBAC) &&
           (type != AS_EXTR)) {
          num_of_guides_total++;
        }
	// Only AS_READ, AS_EBAC, & AS_EXTR fragments are to be used in Gene
	// Myers coverage statistic.

	set_cid_fragment(heapva->frags,ifrag,ifrag); // While we are here ....
      }
      fprintf(stderr,
	      "Total number of guides counted: " F_IID " of " F_IID " fragments\n",
	      num_of_guides_total, nfrag);
      gstate->nfrag_randomly_sampled_in_genome
        = (int)nfrag - (int)num_of_guides_total;

    } 

    if(rg->nbase_in_genome > 0) { 
      gstate->nbase_in_genome = rg->nbase_in_genome;
      gstate->global_fragment_arrival_rate 
        = ((float)gstate->nfrag_randomly_sampled_in_genome)
        / ((float)gstate->nbase_in_genome);
      fprintf(stderr,"gstate->nbase_in_genome=" BPFORMAT "\n",
            gstate->nbase_in_genome);
    }

    // Check for a common error.  The error is trying to run CGB with
    // a fragment graph store that has no overlaps.
    { 
      IntEdge_ID numedge = GetNumEdges(heapva->edges);
      if( numedge == 0 ) {
	fprintf(stderr,
		"CGB Error! You silly rabbit.\n"
		"This fragment graph store has no fragment overlaps.\n"
		"All unitigs would have been be singletons.\n");
      }
      assert(numedge > 0);
    }


    if(TIMINGS) {
      time(&tp1); fprintf(stderr,"Begin marking graph\n");
    }

#ifdef STORE_BRANCH_POINTS_AT_FRAGMENT
    if(branch_points_input_file_name != NULL &&
       branch_points_input_file_name[0] != '\0') {
      fprintf(stderr,"CGB input known branch-points <%s>\n",
              branch_points_input_file_name );
      {
        FragmentHash * afr_to_avx = build_FragmentHash( heapva->frags, 0);
        // STORE_BRANCH_POINTS_AT_FRAGMENT
        char input_line[CMD_BUFFER_SIZE] = {0};
        
        FILE * fbpts = NULL;
        SAFE_FOPEN( fbpts, branch_points_input_file_name,"r");
#if 0
        fopen(branch_points_input_file_name,"r");
        if(NULL == fbpts) {
          fprintf(stderr,"Could not open <%s>\n",
                  branch_points_input_file_name);
          exit(1);
        }
#endif
        assert(NULL != afr_to_avx);

        // Input the known fragment-end branch-points.
        while(NULL != fgets(input_line,CMD_BUFFER_SIZE-2,fbpts)) {
          IntFragment_ID iid, vid;
          int suf, bpt_type, bpt_offset, old_bpt_offset=0;
          if(input_line[0] == '#') continue; // ignore comment lines
          sscanf(input_line,F_IID " %d %d %d",&iid,&suf,&bpt_type,&bpt_offset);
          vid = get_vid_FragmentHash(afr_to_avx,iid);
          // There is no bounds checking here.  This is
          // dangerous. What if the external file has bad data?
          old_bpt_offset = get_bpt_vertex(heapva->frags,vid,suf);
          bpt_offset = max(bpt_offset,old_bpt_offset);
          set_bpt_vertex(heapva->frags,vid,suf,bpt_offset);
        }
        // STORE_BRANCH_POINTS_AT_FRAGMENT
        destroy_FragmentHash(afr_to_avx);
        SAFE_FCLOSE(fbpts);

      fragment_end_edge_trimmer
        (
         /* input only */
         heapva->frags,
         /* modify */
         heapva->edges);
      }
    }
#endif // STORE_BRANCH_POINTS_AT_FRAGMENT

#if 0
    view_fgb_chkpnt( rg->Output_Graph_Store_Prefix, heapva->frags, heapva->edges);
#endif    
#ifdef DEBUG_RISM
  view_fgb_chkpnt( "finished_marking_graph", heapva->frags, heapva->edges);
  check_symmetry_of_the_edge_mates( heapva->frags, heapva->edges, heapva->next_edge_obj);
#endif    
    
    if(TIMINGS) {
      time(&tp2); 
      fprintf(stderr,"%10" F_TIME_TP " sec: Finished marking graph\n",(tp2-tp1));
    }
    
    count_fragment_and_edge_labels( heapva->frags, heapva->edges,
                                    "In main before build 1");

    if(TIMINGS) {
      tp1 = time(NULL); fprintf(stderr,"Begin chunk graph build 1\n");
      system_date();
    }
    chunk_graph_build_1
      (
       /* Input Only */
       rg->Output_Graph_Store_Prefix,
       rg->work_limit_placing_contained_fragments,		  
       rg->walk_depth,
       gstate->max_frag_iid,
       rg->nbase_in_genome,
       rg->chimeras_file,
       rg->spurs_file,
       rg->cgb_unique_cutoff,
       heapva->frags,     /* The internal representation of
                              the fragment reads. I have one
                              extra for the segstart field. */
       heapva->edges,     /* The internal representation of the
                              overlaps. */
       &(gstate->global_fragment_arrival_rate),
       heapva->chunkfrags,
       heapva->thechunks
       );
    if(TIMINGS) {
      tp2 = time(NULL); 
      fprintf(stderr,"%10" F_TIME_TP " sec: Finished chunk graph build 1\n",
              (tp2-tp1));
      system_date();
    }

    count_fragment_and_edge_labels( heapva->frags, heapva->edges,
                                    "In main after build 1");


  // BEGIN CGB ITERATION PHASE
    {
      int iterator;
      for(iterator=0;iterator< rg->num_cgb_passes;iterator++) {

        if(TIMINGS) {
          tp1 = time(NULL); fprintf(stderr,"Begin chunk graph build 2\n");
          system_date();
        }
	{
	  char strtmp1[CMD_BUFFER_SIZE-1];
	  FILE * fbpts1 = fopen
	    (strcat(strcpy(strtmp1,rg->Output_Graph_Store_Prefix),
                    ".cgb_bpts1"),"a");
	  FILE * fbpts2 = fopen
	    (strcat(strcpy(strtmp1,rg->Output_Graph_Store_Prefix),
                    ".cgb_bpts2"),"a");
          assert(NULL != fbpts1);
          assert(NULL != fbpts2);
          
	  chunk_graph_build_2
	    (
	     /* Input Only */
	     rg->use_consensus,
	     rg->dont_find_branch_points,
             rg->cgb_unique_cutoff,
             gstate->global_fragment_arrival_rate,
             TheFragStore,
             /* Input/Output */
	     heapva->frags,     /* The internal representation of
				    the fragment reads. I have one
				    extra for the segstart field. */
	     heapva->edges,     /* The internal representation of the
				    overlaps. */
	     heapva->chunkfrags,
	     heapva->thechunks,
             /* Output Only */
	     heapva->chunkseqs,
	     heapva->chunkquas,
	     heapva->chunksrcs,
	     fbpts1,
             fbpts2
	     );
          SAFE_FCLOSE(fbpts1);
          SAFE_FCLOSE(fbpts2);
        }
        if(TIMINGS) {
          tp2 = time(NULL); 
          fprintf(stderr,"%10" F_TIME_TP " sec: Finished chunk graph build 2\n",
                  (tp2-tp1));
          system_date();
        }


        count_fragment_and_edge_labels( heapva->frags, heapva->edges,
                                        "Before chunk_classification_dvt");

        if(TIMINGS) {
          time(&tp1); fprintf(stderr,"Begin re-marking graph\n");
        }
        chunk_classification_dvt
          ( heapva->frags, heapva->edges,
            rg->walk_depth,
            rg->remove_blizzard_overlaps
            );
        if(TIMINGS) {
          time(&tp2); 
          fprintf(stderr,"%10" F_TIME_TP " sec: Finished re-marking graph\n",(tp2-tp1));
        }
    
        count_fragment_and_edge_labels( heapva->frags, heapva->edges,
                                        "After chunk_classification_dvt");

        if(TIMINGS) {
          tp1 = time(NULL); fprintf(stderr,"Begin chunk graph build 1\n");
          system_date();
        }
        chunk_graph_build_1
          (
           /* Input Only */
	   rg->Output_Graph_Store_Prefix,
	   rg->work_limit_placing_contained_fragments,		  
           rg->walk_depth,
           gstate->max_frag_iid,
           rg->nbase_in_genome,
           rg->chimeras_file,
           rg->spurs_file,
           rg->cgb_unique_cutoff,
           heapva->frags,     /* The internal representation of
                                  the fragment reads. I have one
                                  extra for the segstart field. */
           heapva->edges,     /* The internal representation of the
                                  overlaps. */
           &(gstate->global_fragment_arrival_rate),
           heapva->chunkfrags,
           heapva->thechunks
           );
        if(TIMINGS) {
          tp2 = time(NULL); 
          fprintf(stderr,"%10" F_TIME_TP " sec: Finished chunk graph build 1\n",
                  (tp2-tp1));
          system_date();
        }
      }
    } // END CGB ITERATION PHASE


    count_fragment_and_edge_labels( heapva->frags, heapva->edges,
                                  "Before save_the_chunk_graph");

  } // end of non-incremental phase


#ifdef SWITCH_CONTAINMENT_DIRECTION_CGB
  convert_all_containment_overlaps_direction_TOC
    ( heapva->frags, heapva->edges, heapva->next_edge_obj, FALSE);
  count_fragment_and_edge_labels( heapva->frags, heapva->edges,
                                  "After switching the containment direction");
#endif // SWITCH_CONTAINMENT_DIRECTION_CGB

#if 0
  {
    char strtmp[CMD_BUFFER_SIZE-1];
    sprintf(strtmp,"%s.fgb",rg->Graph_Store_File_Prefix);
    fprintf(stderr,"Last graph is %s\n", strtmp);
    view_fgb_chkpnt( strtmp, heapva->frags, heapva->edges);
  }
#endif    

  /**************** Finish Process Input  *********************/
  system_date();
  
  fprintf(stderr,"close the fragStore\n");

  if(TheFragStore != NULLFRAGSTOREHANDLE) {
    ierr = closeFragStore(TheFragStore); assert(ierr == 0);
    TheFragStore = NULLFRAGSTOREHANDLE;
  }
  
  return( status);
}

#if 0
int main(int argc, char * argv [])
{
  TStateGlobals gstate = {0}
  THeapGlobals  heapva = {0};
  UnitiggerGlobals rg  = {0};
  int status;
  status = main_cgb( argc, argv, &gstate, &heapva, &rg);
  exit(status != EXIT_SUCCESS);
}
#endif
