
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

static char *rcsid = "$Id: AS_FGB_main.c,v 1.25 2008-10-08 22:02:54 brianwalenz Exp $";

#include "AS_CGB_all.h"

int (*compare_edge_function)(const void *a, const void *b);

VA_DEF(IntEdge_ID);

//  AS_CGB_fga.c
void fragment_graph_analysis(Tfragment frags[],
                             Tedge     edges[],
                             FILE      *ffga);


//  AS_CGB_fgb.c
void transitive_edge_marking(Tfragment     * frags,
                             Tedge         * edges,
                             const int walk_depth,
                             const int cutoff_fragment_end_degree,
                             const int work_limit_per_candidate_edge);

//  AS_CGB_edgemate.c
void append_the_edge_mates (Tfragment frags[],
                            Tedge edges[]);



void
input_messages_from_a_file(FILE       *fovl,
                           Tfragment  frags[],
                           Tedge      edges[],
                           IntFragment_ID *afr_to_avx,
                           VA_TYPE(IntEdge_ID)  *next_edge,
                           const int dvt_double_sided_threshold_fragment_end_degree,
                           const int con_double_sided_threshold_fragment_end_degree,
                           const int intrude_with_non_blessed_overlaps_flag,
                           const uint32 overlap_error_threshold);

void
process_gkp_store_for_fragments(char *gkpStoreName,
                                Tfragment   *frags,
                                Tedge       *edges);

void
process_ovl_store(char * OVL_Store_Path,
                  Tfragment  frags[],
                  Tedge      edges[],
                  IntFragment_ID *afr_to_avx,
                  VA_TYPE(IntEdge_ID)  *next_edge,
                  const int dvt_double_sided_threshold_fragment_end_degree,
                  const int con_double_sided_threshold_fragment_end_degree,
                  const int intrude_with_non_blessed_overlaps_flag,
                  const uint32 overlap_error_threshold);




static void output_mesgs(Tfragment frags[],
                         Tedge     edges[],
                         FILE *fcgb) {

  // Output the OVL messages:

  const IntEdge_ID  nedge = GetNumEdges(edges);
  IntEdge_ID ie;
  for(ie=0;ie<nedge;ie++){
    OverlapMesg ovl_mesg;

    const IntFragment_ID avx = get_avx_edge(edges,ie);
    const int asx = get_asx_edge(edges,ie);
    const int ahg = get_ahg_edge(edges,ie);

    const IntFragment_ID bvx = get_bvx_edge(edges,ie);
    const int bsx = get_bsx_edge(edges,ie);
    const int bhg = get_bhg_edge(edges,ie);

    const Tnes    nes = get_nes_edge(edges,ie);
    const uint32  qua = get_qua_edge(edges,ie);

    const IntFragment_ID aid = get_iid_fragment(frags,avx);
    const IntFragment_ID bid = get_iid_fragment(frags,bvx);
    // Assembler internal Fragment ids.

    ovl_mesg.aifrag = aid;
    ovl_mesg.bifrag = bid;

    ovl_mesg.ahg = ahg;
    ovl_mesg.bhg = bhg;
    ovl_mesg.min_offset = ahg;
    ovl_mesg.max_offset = ahg;

    ovl_mesg.orientation =
      ( asx ?
        ( bsx ? AS_INNIE : AS_NORMAL ) :
        ( bsx ? AS_ANTI  : AS_OUTTIE ) );

    switch(nes) {
      case AS_CGB_DOVETAIL_EDGE:
      case AS_CGB_INTERCHUNK_EDGE:
      case AS_CGB_INTRACHUNK_EDGE:
      case AS_CGB_TOUCHES_CONTAINED_EDGE:
      case AS_CGB_BETWEEN_CONTAINED_EDGE:
        ovl_mesg.overlap_type = AS_DOVETAIL; break;
      case AS_CGB_REMOVED_BY_TRANSITIVITY_DVT:
        ovl_mesg.overlap_type = AS_DOVETAIL_CHORD; break;
      case AS_CGB_CONTAINED_EDGE:
        ovl_mesg.overlap_type = AS_CONTAINMENT; break;
      case AS_CGB_REMOVED_BY_TRANSITIVITY_CON:
        ovl_mesg.overlap_type = AS_CONTAINMENT_CHORD; break;
      default:
        fprintf(stderr,"Unexpected overlap edge type: nes=%d\n", nes);
        assert(FALSE);
    }

    ovl_mesg.quality = AS_OVS_decodeQuality(qua);
    ovl_mesg.polymorph_ct = 0;
    ovl_mesg.alignment_trace = NULL;

    if(
       ((is_a_dvt_simple(ahg,bhg))&&(aid < bid))
       // Output only one dovetail overlap edge record per overlap.
       ||
       ((!is_a_dvt_simple(ahg,bhg))&&(asx))
       // Output only one containment overlap edge record per overlap.
       // Choose the NORMAL and INNIE orientations.
       ){
      GenericMesg   pmesg;
      pmesg.t = MESG_OVL;
      pmesg.m = &ovl_mesg;
      WriteProtoMesg_AS(fcgb,&pmesg);
    }
  }
}





static void process_one_ovl_file(const char Batch_File_Name[],
                                 THeapGlobals   * heapva,
                                 IntFragment_ID * afr_to_avx,
                                 VA_TYPE(IntEdge_ID)  * next_edge,
                                 const int dvt_double_sided_threshold_fragment_end_degree,
                                 const int con_double_sided_threshold_fragment_end_degree,
                                 const int intrude_with_non_blessed_overlaps_flag,
                                 const uint32 overlap_error_threshold) {

  FILE *fovl = fopen(Batch_File_Name,"r");
  if(NULL == fovl){
    fprintf(stderr,"* Can not open input file %s\n",Batch_File_Name);
    exit(1);
  }

  input_messages_from_a_file(fovl,
                             heapva->frags,
                             heapva->edges,
                             afr_to_avx,
                             next_edge,
                             dvt_double_sided_threshold_fragment_end_degree,
                             con_double_sided_threshold_fragment_end_degree,
                             intrude_with_non_blessed_overlaps_flag,
                             overlap_error_threshold);

  fclose(fovl);
}




int compare_edge_weak(const void * const aa, const void * const bb)
{
  // This comparison function is used with ANSI qsort() for sorting
  // the edges to form contiguous segments.

  // The lesser edge is the one we keep in the Reaper.
  int icom;
  Aedge *a = (Aedge *)aa;
  Aedge *b = (Aedge *)bb;

  icom = ((a->avx) - (b->avx));
  if( icom == 0 ) {
    icom = ((a->asx) - (b->asx));
    //if( icom == 0 ) {
    //icom = (b->blessed - a->blessed);
      if( icom == 0 ) {
        icom = ((a->ahg) - (b->ahg));
        // Favor the minimum ahg.
        if( icom == 0 ) {
          icom = ((b->bhg) - (a->bhg));
          // Favor the maximum bhg.

          if( icom == 0 ) {
            // The following is unnecessary, but useful for the binary
            // search in the adjaceny lists and regression output.
            icom = ((a->bvx) - (b->bvx));
            if( icom == 0 )
              icom = ((a->bsx) - (b->bsx));
          }
        } // End of regression stuff.
      }
      //}
  }
  return icom ;
}


int compare_edge_strong(const void * const aa, const void * const bb)
{
  // This comparison function is used with ANSI qsort() for sorting
  // the edges to form contiguous segments.

  // The lesser edge is the one we keep in the Reaper.
  int icom;
  Aedge *a = (Aedge *)aa;
  Aedge *b = (Aedge *)bb;

  icom = ((a->avx) - (b->avx));
  if( icom == 0 ) {
    icom = ((a->asx) - (b->asx));
    if( icom == 0 ) {
      icom = (b->blessed - a->blessed);
      if( icom == 0 ) {
        icom = ((a->ahg) - (b->ahg));
        // Favor the minimum ahg.
        if( icom == 0 ) {
          icom = ((b->bhg) - (a->bhg));
          // Favor the maximum bhg.

          if( icom == 0 ) {
            // The following is unnecessary, but useful for the binary
            // search in the adjaceny lists and regression output.
            icom = ((a->bvx) - (b->bvx));
            if( icom == 0 ) {
              icom = ((a->bsx) - (b->bsx));
              if( icom == 0 )
                icom = (a->reflected - b->reflected);
            }
          }
        } // End of regression stuff.
      }
    }
  }
  return icom ;
}



int main_fgb(THeapGlobals  * heapva,
             UnitiggerGlobals * rg) {
  int status = 0;
  int ierr = 0;

  IntFragment_ID *afr_to_avx = NULL;

  int did_processing_phase_2 = FALSE;

  compare_edge_function = compare_edge_strong;

  if((rg->frag_store) && (rg->frag_store[0] != '\0'))
    process_gkp_store_for_fragments(rg->frag_store,
                                    heapva->frags,
                                    heapva->edges);


  // Re-hash the fragment IID to fragment VID mapping using the
  // fragments in the store.  (duplicated, search for BUILD_AFR_TO_AVX
  {
    IntFragment_ID  iv    = 0;
    IntFragment_ID  im    = 0;
    IntFragment_ID  nfrag = GetNumFragments(heapva->frags);

    for (iv=0; iv<nfrag; iv++) {
      IntFragment_ID iid = get_iid_fragment(heapva->frags,iv);
      im = MAX(im, iid);
    }

    assert(im < AS_CGB_NOT_SEEN_YET);

    afr_to_avx = safe_calloc(im + 1, sizeof(IntFragment_ID));

    for(iv=0; iv<nfrag; iv++)
      afr_to_avx[get_iid_fragment(heapva->frags,iv)] = iv;
  }

  VA_TYPE(IntEdge_ID)  *next_edge = CreateVA_IntEdge_ID(rg->maxedges);

  if (rg->ovl_files_list_fname != NULL)
    process_one_ovl_file(rg->ovl_files_list_fname,
                         heapva,
                         afr_to_avx,
                         next_edge,
                         rg->dvt_double_sided_threshold_fragment_end_degree,
                         rg->con_double_sided_threshold_fragment_end_degree,
                         rg->intrude_with_non_blessed_overlaps_flag,
                         rg->overlap_error_threshold);



  // Process the blessed ovl file.
  if(rg->blessed_overlaps_input_filename){
    assert(0 == GetNumEdges(heapva->edges));

    process_one_ovl_file(rg->blessed_overlaps_input_filename,
                         heapva,
                         afr_to_avx,
                         next_edge,
                         rg->dvt_double_sided_threshold_fragment_end_degree,
                         rg->con_double_sided_threshold_fragment_end_degree,
                         rg->intrude_with_non_blessed_overlaps_flag,
                         AS_OVS_encodeQuality(1.0));

    {
      const IntFragment_ID nfrag = GetNumFragments(heapva->frags);
      const IntEdge_ID nedge = GetNumEdges(heapva->edges);
      IntEdge_ID ie;
      for( ie=0; ie < nedge; ie ++) {
        const IntFragment_ID avx = get_avx_edge(heapva->edges,ie);
        const IntFragment_ID asx = get_asx_edge(heapva->edges,ie);
        // assert(nedge > ie);
        assert(nfrag > avx);
        set_blessed_edge(heapva->edges, ie, TRUE);
        set_blessed_vertex(heapva->frags, avx, asx, TRUE);
      }
    }
  }


  // Process the bubble smoothing ovl file.
  if(rg->bubble_overlaps_filename[0])
    process_one_ovl_file(rg->bubble_overlaps_filename,
                         heapva,
                         afr_to_avx,
                         next_edge,
                         rg->dvt_double_sided_threshold_fragment_end_degree,
                         rg->con_double_sided_threshold_fragment_end_degree,
                         rg->intrude_with_non_blessed_overlaps_flag,
                         rg->overlap_error_threshold);


  if ((rg->OVL_Store_Path) && (rg->OVL_Store_Path[0] != '\0'))
    process_ovl_store(rg->OVL_Store_Path,
                      heapva->frags,
                      heapva->edges,
                      afr_to_avx,
                      next_edge,
                      rg->dvt_double_sided_threshold_fragment_end_degree,
                      rg->con_double_sided_threshold_fragment_end_degree,
                      rg->intrude_with_non_blessed_overlaps_flag,
                      rg->overlap_error_threshold);

  safe_free(afr_to_avx);

  Delete_VA(next_edge);

  /////////////////////////////////////////////////////////////////

  //count_fragment_and_edge_labels ( heapva->frags, heapva->edges, "Before separate_fragments_as_solo_hanging_thru");

  separate_fragments_as_solo_hanging_thru(heapva->frags,heapva->edges);

  compare_edge_function = compare_edge_weak;

  // Now do not distinguish the blessed overlaps.

  reorder_edges( heapva->frags, heapva->edges);
  //count_fragment_and_edge_labels( heapva->frags, heapva->edges, "RISM_reorder_edges");
  //check_symmetry_of_the_edge_mates( heapva->frags, heapva->edges);

  append_the_edge_mates( heapva->frags, heapva->edges);

  // Currently the spur finding and micro-bubble smoothing code
  // needs the dovetail edge mates to exist.

  //count_fragment_and_edge_labels( heapva->frags, heapva->edges, "RISM_append_the_edge_mates");
  //check_symmetry_of_the_edge_mates( heapva->frags, heapva->edges);


  //  This routine marks duplicated edges for deletion.
  //  It is only effective when the edges are fully sorted.
  //
  {
    IntEdge_ID nedge = GetNumEdges(heapva->edges);
    IntEdge_ID ie1   = 0;
    int64      count = 0;
    int        icompare;

    for (ie1=1; ie1<nedge; ie1++) {
      Aedge edge0 = *GetVA_Aedge(heapva->edges,ie1-1);
      Aedge edge1 = *GetVA_Aedge(heapva->edges,ie1);

      edge0.blessed   = edge1.blessed   = FALSE;
      edge0.reflected = edge1.reflected = FALSE;

      icompare = (*compare_edge_function)(&edge1,&edge0);
      assert(icompare >= 0);

      if (icompare == 0) {
        count++;

        if (is_a_dvt_edge(heapva->edges,ie1) == TRUE)
          set_nes_edge(heapva->edges,ie1,AS_CGB_REMOVED_BY_DUPLICATE_DVT);
        else
          set_nes_edge(heapva->edges,ie1,AS_CGB_REMOVED_BY_DUPLICATE_CON);
      }
    }

    fprintf(stderr, F_S64" duplicate edges found.\n", count);
  }


  //count_fragment_and_edge_labels( heapva->frags, heapva->edges, "RISM_delete_duplicate_edges");
  //check_symmetry_of_the_edge_mates( heapva->frags, heapva->edges);

  // Currently the spur finding code needs the dovetail edge mates to
  // exist.

  if( (GetNumFragments(heapva->frags) > 0) &&
      (GetNumEdges(heapva->edges) > 0) &&
      (rg->dechord_the_graph) &&
      (!did_processing_phase_2)
      ) /* The FGB computations need a populated graph store. */ {

    transitive_edge_marking(heapva->frags,
                            heapva->edges,
                            rg->walk_depth,
                            rg->cutoff_fragment_end_degree,
                            rg->work_limit_per_candidate_edge);
  }

  identify_early_spur_fragments( heapva->frags, heapva->edges);

  //count_fragment_and_edge_labels ( heapva->frags, heapva->edges, "Before contained_fragment_marking_frc");
  contained_fragment_marking_frc( heapva->frags, heapva->edges);
  //count_fragment_and_edge_labels ( heapva->frags, heapva->edges, "After contained_fragment_marking_frc");


  char  thePath[FILENAME_MAX];
  sprintf(thePath,"%s.fga.ckp",rg->Output_Graph_Store_Prefix);
  FILE *ffga = fopen(thePath,"w");

  fragment_graph_analysis(heapva->frags,
                          heapva->edges,
                          ffga);
  fclose(ffga);


  if( rg->create_dump_file ) {
    FILE *folp = NULL;
    char strtmp[FILENAME_MAX];

    fprintf(stderr,"Opening dump file to write a batch of "
	    "ADT+IDT+OFG+OVL messages.\n");

    if(NULL == (folp = fopen(rg->Dump_File_Name,"w"))){
      fprintf(stderr,"* Can not open output file %s\n",rg->Dump_File_Name);
      exit(1);
    }

    output_mesgs (heapva->frags,
                  heapva->edges,
                  folp);
    fclose(folp);
  }

  return( status);
}
