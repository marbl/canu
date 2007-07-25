
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

static char CM_ID[] = "$Id: AS_CGB_main.c,v 1.15 2007-07-25 10:29:50 brianwalenz Exp $";

#include "AS_UTL_version.h"  
#include "AS_CGB_all.h"
#include "AS_CGB_unitigger_globals.h"

#undef SWITCH_CONTAINMENT_DIRECTION_CGB

extern int REAPER_VALIDATION;

//  AS_CGB_classify.c
void chunk_classification_dvt(Tfragment frags[],
                              Tedge edges[],
                              const int walk_depth);

//  AS_CGB_cgb.c
void chunk_graph_build_1(const char * const Graph_Store_File_Prefix,
                         const int walk_depth,
                         const int64  genome_length,
                         const char * chimeras_file,
                         const char * spurs_file,
                         const int recalibrate_global_arrival_rate,
                         const float cgb_unique_cutoff,
                         Tfragment frags[],
                         Tedge     edges[],
                         float         *global_fragment_arrival_rate,
                         TChunkFrag    *chunkfrags,
                         TChunkMesg    *thechunks);


static IntEdge_ID get_the_thickest_dvt_overlap_from_vertex
( Tfragment * frags,
  Tedge     * edges,
  IntFragment_ID avx,
  int            asx
  ) {
  // Does the graph have to be in adjacency format?
  const IntEdge_ID nedge = GetNumEdges(edges);
  IntEdge_ID ie_thickest = AS_CGB_EDGE_NOT_FOUND;
  
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

  assert((AS_CGB_EDGE_NOT_FOUND == ie_thickest) || (nedge > ie_thickest));

  return ie_thickest;
}

static void identify_thickest_overlaps
( Tfragment * frags,
  Tedge     * edges
  ) 
{
  const IntEdge_ID nedge = GetNumEdges(edges);
  const IntFragment_ID nfrag = GetNumFragments(frags);
  
  // Assume an adjaceny representation.
  IntFragment_ID avx;
  int asx;
  for(avx=0; avx < nfrag; avx++) {
    for(asx=0; asx < 2; asx++) {
      const IntEdge_ID ie = get_the_thickest_dvt_overlap_from_vertex ( frags, edges, avx, asx );
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

static void identify_bmpc_paths
( Tfragment * frags,
  Tedge     * edges,
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

      ie_next = get_the_thickest_dvt_overlap_from_vertex (frags,edges,bvx,!bsx);
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
  Tedge     * edges
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
          const IntEdge_ID ie = get_the_thickest_dvt_overlap_from_vertex ( frags, edges, avx, asx );
          if(AS_CGB_EDGE_NOT_FOUND != ie) {
            const Tnes nes = get_nes_edge(edges,ie);
            assert((AS_CGB_THICKEST_EDGE == nes) || (AS_CGB_INTERCHUNK_EDGE == nes));
            if(REAPER_VALIDATION) {
              if((AS_CGB_THICKEST_EDGE == nes)) {
                set_nes_edge(edges,ie,AS_CGB_INTERCHUNK_EDGE);
                fix_overlap_edge_mate(frags,edges,ie);
              }
            } else {
              identify_bmpc_paths( frags, edges, ie);
            }
          }
        }
      }
    }
  }
}




#ifdef  SWITCH_CONTAINMENT_DIRECTION_CGB

static void reflect_an_edge_in_place(Tedge * edges,
                                     IntEdge_ID iedge) {
  Aedge * the_edge = GetVA_Aedge(edges,iedge);
  reflect_Aedge(the_edge,the_edge);
}


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
}

static void convert_all_containment_overlaps_direction_TOC
( Tfragment * frags,
  Tedge     * edges,
  const int become_to_contained
 )
{
  time_t tp1, tp2;
  reflect_containment_direction_in_place( edges, become_to_contained);
  reorder_edges( frags, edges);
}
#endif // SWITCH_CONTAINMENT_DIRECTION_CGB

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
    default:
      fprintf(stderr,"Unexpected overlap label ines=%d\n",ines);
      assert(FALSE);
      break;
    }
  }}
}


int main_cgb(THeapGlobals  * heapva,
             UnitiggerGlobals * rg) {

  GateKeeperStore *gkpStore = openGateKeeperStore(rg->frag_store, FALSE);
  
  //count_fragment_and_edge_labels( heapva->frags, heapva->edges, "After reading the fragment graph store");
  //check_symmetry_of_the_edge_mates( heapva->frags, heapva->edges);

  if(rg->aggressive_spur_fragment_marking) {
    maskout_overlaps_touching_crappy_fragments ( heapva->frags, heapva->edges);

    // Do we need to patch up the graph here??

    //count_fragment_and_edge_labels ( heapva->frags, heapva->edges, "after maskout_overlaps_touching_crappy_fragments");
    //view_fgb_chkpnt( "maskout_overlaps_touching_crappy_fragments", heapva->frags, heapva->edges);
    //check_symmetry_of_the_edge_mates( heapva->frags, heapva->edges);
  }

  
  identify_thickest_overlaps( heapva->frags, heapva->edges);
  //count_fragment_and_edge_labels( heapva->frags, heapva->edges, "After identify_thickest_edges");
  //view_fgb_chkpnt( "identify_thickest_overlaps", heapva->frags, heapva->edges);
  //check_symmetry_of_the_edge_mates( heapva->frags, heapva->edges);

  identify_essential_components( heapva->frags, heapva->edges );
  //count_fragment_and_edge_labels( heapva->frags, heapva->edges, "After identify_essential_components");
  //view_fgb_chkpnt( "identify_essential_components", heapva->frags, heapva->edges);
  //check_symmetry_of_the_edge_mates( heapva->frags, heapva->edges);

  chunk_classification_dvt(heapva->frags, heapva->edges, rg->walk_depth);
  
  //count_fragment_and_edge_labels( heapva->frags, heapva->edges, "In main after chunk_classification_dvt");
  //check_symmetry_of_the_edge_mates( heapva->frags, heapva->edges);
   
#ifdef SWITCH_CONTAINMENT_DIRECTION_CGB
  convert_all_containment_overlaps_direction_TOC ( heapva->frags, heapva->edges, TRUE);
  //count_fragment_and_edge_labels( heapva->frags, heapva->edges, "After switching the containment direction");
#endif // SWITCH_CONTAINMENT_DIRECTION_CGB


  { // beginning of the non-incremental phase
    
    IntFragment_ID nfrag = GetNumFragments(heapva->frags);
    IntFragment_ID ifrag;
    IntFragment_ID num_of_guides_total = 0;

    /* Count the total amount of guide fragments. */

    for(ifrag=0;ifrag<nfrag;ifrag++) { 
      const FragType type = get_typ_fragment(heapva->frags,ifrag);

      // Only AS_READ & AS_EXTR fragments are to be used in Gene Myers
      // coverage statistic.
      //
      if((type != AS_READ) && (type != AS_EXTR))
        num_of_guides_total++;

      set_cid_fragment(heapva->frags,ifrag,ifrag); // While we are here ....
    }

    fprintf(stderr, "Total number of guides counted: " F_IID " of " F_IID " fragments\n",
            num_of_guides_total, nfrag);

    if(rg->genome_length > 0) { 
      heapva->nbase_in_genome = rg->genome_length;
      heapva->global_fragment_arrival_rate = ((float)(nfrag - num_of_guides_total) /
                                              (float)heapva->nbase_in_genome);
    }

    // Check for a common error.  The error is trying to run CGB with
    // a fragment graph store that has no overlaps.
    if (GetNumEdges(heapva->edges) == 0)
      fprintf(stderr,
              "CGB Error! You silly rabbit.\n"
              "This fragment graph store has no fragment overlaps.\n"
              "All unitigs would have been be singletons.\n");
    assert(GetNumEdges(heapva->edges) > 0);

    //view_fgb_chkpnt( "finished_marking_graph", heapva->frags, heapva->edges);
    //check_symmetry_of_the_edge_mates( heapva->frags, heapva->edges);
    //count_fragment_and_edge_labels( heapva->frags, heapva->edges, "In main before build 1");

    chunk_graph_build_1(rg->Output_Graph_Store_Prefix,
                        rg->walk_depth,
                        rg->genome_length,
                        rg->chimeras_file,
                        rg->spurs_file,
                        rg->recalibrate_global_arrival_rate,
                        rg->cgb_unique_cutoff,
                        heapva->frags,
                        heapva->edges,
                        &(heapva->global_fragment_arrival_rate),
                        heapva->chunkfrags,
                        heapva->thechunks);

    //count_fragment_and_edge_labels( heapva->frags, heapva->edges, "In main after build 1");

  } // end of non-incremental phase


#ifdef SWITCH_CONTAINMENT_DIRECTION_CGB
  convert_all_containment_overlaps_direction_TOC ( heapva->frags, heapva->edges, FALSE);
  //count_fragment_and_edge_labels( heapva->frags, heapva->edges, "After switching the containment direction");
#endif // SWITCH_CONTAINMENT_DIRECTION_CGB

  closeGateKeeperStore(gkpStore);

  return(0);
}
