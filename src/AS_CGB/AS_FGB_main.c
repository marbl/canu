
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
= "$Id: AS_FGB_main.c,v 1.11 2007-04-12 18:54:44 brianwalenz Exp $";
/*********************************************************************
 *
 * Module:  AS_FGB_main.c
 * Description:
 *    Unitigger main
 * 
 *    Reference: FragmentGraphBuilder.rtf
 *
 *    Command Line Interface:
 *        $ fgb
 * -p <parameterPath>    If specified, parameters are loaded from this path
 * -i <inputStorePath>    The path of the input store
 * -o <outputStorePath>    The path of the output store
 * -c There is no input store, create a new store for output.
 *    No -i can be specified.
 * -a Append to the specified input store
 *    No -o can be specified).  Upon failure, the store should not be altered.
 * -f Overwrite the output store, if it exists
 *    (default is to fail if the output exists)
 * -X Enable developer's options.
 *    Any option not stated in the tool's SOP should require this flag
 * -P Output is ASCII (default is binary)
 * -D <level>  Specify debug level (only valid with -X)
 * -v <level>  Specify verbosity level (only valid with -X)
 * -A <level> Run the graph analyzer
 * -E <int>  The number of errors allowed.
 * -C <string> The check-pointing information.
 * -R <string> The restart information.
 * -T Output the chord overlaps as well.
 * -n <int> The expected number of fragment reads.
 * -m <int> The expected number of fragment overlaps.
 * -x <int> The dovetail double-sided degree threshold.
 * -y <int> The dovetail single-sided degree threshold.
 * -z <int> The containment degree threshold.
 * -r <int> Which global containment rule to use.
 * -w <int> Limit of search depth for fragment graph walking in 
 *          transitively inferable overlap removal.
 * <OverlapInputFile>*
 * 
 * Examples:
 * -c -o <outputStorePath>  creates a store to collect the output
 * -c -i   is an error condition
 * -a -i <inputStorePath>  appends the output to the input store
 * -a -o   is an error condition
 *
 *       GraphStorePath:  Used to read/write a fragment graph store.
 *       It must be named *.fgb.
 *      
 *       OverlapInputFile: The file with new OVL records to process. 
 *       It must be named <inputFile>.ovl.
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
 *       See documentation for descriptions of the checks
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
#include "AS_FGB_hanging_fragment.h"
#include "AS_FGB_contained.h"
#include "AS_CGB_unitigger_globals.h"
#include "AS_FGB_buildFragmentHash.h"

/****************************************************************************/
#define DEBUGGING
#define PREPEND_IBA_FILE
#undef DEBUG_VIEW1
#undef DEBUG_VIEW2
#undef DEBUG_VIEW3
#undef DEBUG_VIEW4
#undef DEBUG_VIEW5
#undef DEBUG_VIEW6
#undef DEBUG_RISM

/* These are optimization parameters. */

/****************************************************************************/
/* File Scope Globals */

static int TIMINGS = TRUE;
static int TIMINGS1 = FALSE;

/****************************************************************************/

static void processing_phase_3
(
  TStateGlobals * gstate,
  THeapGlobals  * heapva,
  const char Output_Graph_Store[],
  const int analysis_flag
)
{
  system_date();
  ReportHeapUsage_CGB( gstate, heapva, stderr);

#if 0
  view_fgb_chkpnt( Output_Graph_Store, heapva->frags, heapva->edges);
  system_date();
#endif

  if(NULL != Output_Graph_Store) {
    { 
      int ierr=0;
      char thePath1[CMD_BUFFER_SIZE-1]={0};
      char thePath2[CMD_BUFFER_SIZE-1]={0};
      sprintf(thePath1,"%s/%s",Output_Graph_Store,"fgb.ckp_tmp");
      sprintf(thePath2,"%s/%s",Output_Graph_Store,"fgb.ckp");
      write_fgb_store(thePath1, gstate, heapva);
      system_date();
      ierr = accept_tmp_as_final_file( thePath1, thePath2);
      assert(ierr == 0);
      // The temporary checkpoint file was moved into its final
      // position.
    }
    if(analysis_flag) {
      const int ProcessFragmentAnnotationsForSimulatorCoordinates
        = (analysis_flag > 1);
      FILE *ffga=NULL;
      int ierr=0;
      char thePath3[CMD_BUFFER_SIZE-1]={0};
      char thePath4[CMD_BUFFER_SIZE-1]={0};
      sprintf(thePath3,"%s/%s",Output_Graph_Store,"fga.ckp_tmp");
      sprintf(thePath4,"%s/%s",Output_Graph_Store,"fga.ckp");
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
      system_date();
      ierr = accept_tmp_as_final_file( thePath3, thePath4);
      assert(ierr == 0);
      // The temporary analysis file could not be moved into its final
      // position.
    }

    {
      char thePath5[CMD_BUFFER_SIZE-1]={0};
      char thePath6[CMD_BUFFER_SIZE-1]={0};
      int ierr=0;
      sprintf(thePath5,"%s/%s",Output_Graph_Store,"fgb.iba_tmp");
      sprintf(thePath6,"%s/%s",Output_Graph_Store,"fgb.iba");
      fprintf(stderr,"thePath5=<%s>\n", thePath5);
      fprintf(stderr,"thePath6=<%s>\n", thePath6);
      system_date();
      ierr = accept_tmp_as_final_file( thePath5, thePath6);
      assert(ierr == 0);
      // The temporary batch info file could not be moved into its final
      // position.
    }
    system_date();
  }
}

static void output_mesgs
(/* Input Only*/
 const Tfragment frags[],
 const Tedge     edges[],
 const VA_TYPE(char) fragsrc[],
 /* Read Only */
 /* Append Only*/
 FILE *fcgb)
{

  // Output the OFG messages:
  {
    const IntFragment_ID nfrag = GetNumFragments(frags);
    IntFragment_ID iv;
    for(iv=0;iv<nfrag;iv++){
      OFGMesg ofg_mesg;

      ofg_mesg.action     = (get_del_fragment(frags,iv)
                             ? AS_DELETE : AS_ADD);
      ofg_mesg.eaccession = get_uid_fragment(frags,iv);
      ofg_mesg.iaccession = get_iid_fragment(frags,iv);
      ofg_mesg.type       = get_typ_fragment(frags,iv);
      //Locale_ID    		elocale;
      //SeqInterval  		locale_pos;
      ofg_mesg.entry_time = 0;
      ofg_mesg.clear_rng.bgn = 0;
      ofg_mesg.clear_rng.end = get_length_fragment(frags,iv);
      ofg_mesg.source = NULL;
  
      if(fragsrc != NULL) {
	const size_t isrc = get_src_fragment(frags,iv);
	ofg_mesg.source = Getchar(fragsrc,isrc);
      } else {
        ofg_mesg.source = "";
      }
      
      {
        assert(0);        //  MESG_OFR was removed, so we now crash
        GenericMesg pmesg;
        //pmesg.t = MESG_OFR;
        pmesg.m = &ofg_mesg;
        WriteProtoMesg_AS(fcgb,&pmesg);
      }
    }
  }
  
  // Output the OVL messages:
  {
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

      signed char delta[1] = {AS_ENDOF_DELTA_CODE};
      
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
      ovl_mesg.delta = delta;

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
  
#if 0  
  if(analysis_flag) {
    fprintf(stderr,"\n\nHistogram of the OVL types\n");
    print_histogram(stderr, OVL_types_histogram, 0, 1);
  }
  if(NULL != ovl_types_histogram) {
    free_histogram(fom_types_histogram);
  }
#endif
}

static void dump_next_edge_obj
( FILE * fout,
  VA_TYPE(IntEdge_ID) * next_edge_obj
  )
{
  const size_t length = GetNumElements_VA(next_edge_obj);
  size_t ii;
  for(ii=0; ii< length; ii++) {
    const size_t jj = *GetVA_IntEdge_ID(next_edge_obj,ii);
    if( ii != jj) { 
      fprintf(fout,F_SIZE_T " " F_SIZE_T "\n", ii, jj );
    }
  }
}

/****************************************************************************/


static void process_one_ovl_file
(int        argc, 
 char       *argv[],
 const char Batch_File_Name[],
 const char Output_Graph_Store[],
 TStateGlobals  * gstate,
 THeapGlobals   * heapva,
 FragmentHashObject * afr_to_avx,    // Part of a hash table replacement
 const int dvt_double_sided_threshold_fragment_end_degree,
 const int con_double_sided_threshold_fragment_end_degree,
 const int intrude_with_non_blessed_overlaps_flag,
 const uint32 overlap_error_threshold,
 const int check_point_level
)
{ // Process the ovl files.
  
  if( check_point_level == 0 ) { // Process the ovl file....

    FILE *fovl = fopen(Batch_File_Name,"r");
    if(NULL == fovl){
      fprintf(stderr,"* Can not open input file %s\n",Batch_File_Name);
      exit(1);
    }
    
    input_messages_from_a_file
      (argc, argv, // For ADT version stamp.
       fovl,
       (heapva->frags), // The internal representation of the fragments.
       (heapva->edges), // The internal representation of the overlaps.
       (heapva->frag_annotations),
       afr_to_avx, 
       &(gstate->min_frag_iid),
       &(gstate->max_frag_iid),
       (heapva->next_edge_obj),
       &(gstate->nbase_in_genome),
       dvt_double_sided_threshold_fragment_end_degree,
       con_double_sided_threshold_fragment_end_degree,
       intrude_with_non_blessed_overlaps_flag,
       overlap_error_threshold
       );
    
    fclose(fovl);
    
  }
} // Process the ovl files.

static void delete_duplicate_edges
(/* Input Only */
 Tfragment frags[], 
 /* Modify */
 Tedge edges[],
 TIntEdge_ID *next_edge_obj
 )
{
  /* 
     This routine marks duplicated edges for deletion.
     It is only effective when the edges are fully sorted.
  */
  //const IntFragment_ID nfrag = GetNumFragments(frags);
  const IntEdge_ID nedge = GetNumEdges(edges);
  IntEdge_ID ie1=0;
  time_t tp1 = 0,tp2;

  const QsortCompare compare_edge = get_compare_edge_function();

  IntEdge_ID count=0;

  if(TIMINGS) {
    tp1 = time(NULL); fprintf(stderr,"Begin delete duplicate edges\n");
    system_top();
  }

  for(ie1=1; ie1 < nedge; ie1++) {
    const IntEdge_ID ie0 = ie1 - 1;
    Aedge edge0 = *GetVA_Aedge(edges,ie0);
    Aedge edge1 = *GetVA_Aedge(edges,ie1);
    edge0.blessed = FALSE;
    edge1.blessed = FALSE;
    edge0.reflected = FALSE;
    edge1.reflected = FALSE;
    {
      const int icompare =
        compare_edge(&edge1,&edge0);
      // This comparison function must satisfy the UNIX qsort semantics.
      
      // The comparison function will be called with two parameters that
      // point to the two elements to be compared. The comparison
      // function must return an integer less than, equal to, or greater
      // than zero, depending on whether the first element in the
      // comparison is considered less than, equal to, or greater than
      // the second element.
      
      assert(icompare >= 0);
      // The edges must be sorted before running this routine.
      
      if( icompare == 0 ) {
        
        count++;
        if( is_a_dvt_edge(edges,ie1) == TRUE ) {
          set_nes_edge(edges,ie1,AS_CGB_REMOVED_BY_DUPLICATE_DVT);
        } else {
          // const int ahg = get_ahg_edge(edges,ie1);
          // const int bhg = get_bhg_edge(edges,ie1);
          set_nes_edge(edges,ie1,AS_CGB_REMOVED_BY_DUPLICATE_CON);
        }
        //fix_overlap_edge_mate(frags,edges,ie1);
        {
          const IntFragment_ID avx = get_avx_edge(edges,ie1);
          const IntFragment_ID bvx = get_bvx_edge(edges,ie1);
          
          fprintf(stdout,
                  "DELETED DUPLICATE EDGE " F_IID ":"
                  "(" F_U64 "," F_IID "," F_IID ") %d %d, "
                  "(" F_U64 "," F_IID "," F_IID ") %d %d : %d "
                  " %d %d "
                  "\n",
                  ie1,
                  get_uid_fragment(frags,avx),
                  get_iid_fragment(frags,avx),
                  get_avx_edge(edges,ie1),
                  get_asx_edge(edges,ie1),
                  get_ahg_edge(edges,ie1),
                  get_uid_fragment(frags,bvx),
                  get_iid_fragment(frags,bvx),
                  get_bvx_edge(edges,ie1),
                  get_bsx_edge(edges,ie1),
                  get_bhg_edge(edges,ie1),
                  get_nes_edge(edges,ie1),
                  get_qua_edge(edges,ie1),
                  get_blessed_edge(edges,ie1)
                  );
        }
      }
    }
  }
  fprintf(stderr, " " F_IID " duplicate edges found.\n",count);
  if(TIMINGS) {
    tp2 = time(NULL); 
    fprintf(stderr,"%10" F_TIME_TP " sec: Finished duplicate edge deletion\n",(tp2-tp1));
    system_top();
  }
  if(count > 0) {
    ResetToRangeVA_IntEdge_ID(next_edge_obj,0);
    // The next_edge_obj data is invalidated.
  }
}

static void check_edges3( Tfragment frags[], Tedge edges[]) {
  const IntEdge_ID nedge = GetNumEdges(edges);
  IntEdge_ID ir2;
  for(ir2=0;ir2<nedge; ir2++) {
    const IntFragment_ID ir2avx = get_avx_edge(edges,ir2);
    const IntFragment_ID ir2bvx = get_bvx_edge(edges,ir2);
    const int ir2ahg = get_ahg_edge(edges,ir2);
    const int ir2bhg = get_bhg_edge(edges,ir2);
    const int ir2aln = get_length_fragment(frags,ir2avx);
    const int ir2bln = get_length_fragment(frags,ir2bvx);
    if( (ir2aln <= ir2ahg) || (ir2bln <= ir2bhg) ||
        (ir2aln <= -ir2bhg) || (ir2bln <= -ir2ahg) ) {
      const IntFragment_ID ir2afr = get_iid_fragment(frags,ir2avx);
      const IntFragment_ID ir2bfr = get_iid_fragment(frags,ir2bvx);
      const Tnes ir2nes = get_nes_edge(edges,ir2);
      fprintf(stderr,"ERROR: Bad edge "
              " ir2afr=" F_IID " ir2bfr=" F_IID " ir2avx=" F_IID " ir2bvx=" F_IID " ir2ahg=%d ir2bhg=%d ir2nes=%d\n",
              ir2afr, ir2bfr, ir2avx, ir2bvx, ir2ahg, ir2bhg, ir2nes);
      
    }
  }
}

/****************************************************************************/

int main_fgb
(
 int argc,
 char * argv [],
 TStateGlobals * gstate,
 THeapGlobals  * heapva,
 UnitiggerGlobals * rg
 )
{
  int status = 0;
  int ierr = 0;

  time_t tp1 = 0,tp2; 

  /* The Store/Checkpoint */
  FragmentHashObject *afr_to_avx = NULL;

  int did_processing_phase_2 = FALSE;

  set_compare_edge_function(compare_edge_strong);
  // This is for the blessed overlaps.

  
  // Re-hash the fragment IID to fragment VID mapping using the
  // fragments in the store.
  assert(afr_to_avx == NULL);
  afr_to_avx = build_FragmentHash(heapva->frags, rg->as_cgb_max_frag_iid );
  assert(afr_to_avx != NULL);

  // WARNING do we need to rebuild the next_edge_obj here?

  if(NULL != heapva->the_ofg_messages ) {
    IntFragment_ID nofg = 0;
    
    fprintf(stderr,"Process the ofg VA.\n");

    InsertFragmentsIntoGraph
      (
       heapva,
       heapva->the_ofg_messages, // Additional fragments in a VA.
       heapva->the_ofg_source,
       afr_to_avx,        // Part of a hash table replacement
       &nofg,
       &(gstate->min_frag_iid),
       &(gstate->max_frag_iid),
       heapva->frag_annotations
       );
  }
  

  if (rg->ovl_files_list_fname != NULL)
      process_one_ovl_file
        ( argc, 
          argv,
          rg->ovl_files_list_fname,
          rg->Output_Graph_Store,
          gstate,
          heapva,
          afr_to_avx,
          rg->dvt_double_sided_threshold_fragment_end_degree,
          rg->con_double_sided_threshold_fragment_end_degree,
          rg->intrude_with_non_blessed_overlaps_flag,
          rg->overlap_error_threshold,
          rg->check_point_level
          );
  

  if(NULL != rg->blessed_overlaps_input_filename){
    // Process the blessed ovl file.

    assert(0 == GetNumEdges(heapva->edges));

    process_one_ovl_file
      ( argc, 
        argv,
        rg->blessed_overlaps_input_filename,
        rg->Output_Graph_Store,
        gstate,
        heapva,
        afr_to_avx,
        rg->dvt_double_sided_threshold_fragment_end_degree,
        rg->con_double_sided_threshold_fragment_end_degree,
        rg->intrude_with_non_blessed_overlaps_flag,
        AS_OVS_encodeQuality(1.0),
        rg->check_point_level
        );
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
    // The blessed overlaps are loaded.
  }
  

  if(NULL != rg->bubble_overlaps_filename){
    // Process the bubble smoothing ovl file.

    process_one_ovl_file
      ( argc, 
        argv,
        rg->bubble_overlaps_filename,
        rg->Output_Graph_Store,
        gstate,
        heapva,
        afr_to_avx,
        rg->dvt_double_sided_threshold_fragment_end_degree,
        rg->con_double_sided_threshold_fragment_end_degree,
        rg->intrude_with_non_blessed_overlaps_flag,
        rg->overlap_error_threshold,
        rg->check_point_level
        );
  }

  if( NULL != heapva->the_ovl_messages ) {
    
    fprintf(stderr,"Process the ovl VA.\n");

    InsertOverlapsIntoGraph
      (
       heapva,
       heapva->the_ovl_messages, // Additional overlaps in a VA.
       heapva->the_ovl_source,
       afr_to_avx,        // Part of a hash table replacement
       rg->dvt_double_sided_threshold_fragment_end_degree,
       rg->con_double_sided_threshold_fragment_end_degree,
       rg->intrude_with_non_blessed_overlaps_flag,
       rg->overlap_error_threshold
       );

  }

  if((NULL != rg->OVL_Store_Path) &&
     (rg->OVL_Store_Path[0] != '\0')) {
    if(TIMINGS) {
      tp1 = time(NULL);
      fprintf(stderr,"Begin reading overlap store\n");
      system_top();
    }
    process_ovl_store
      ( rg->OVL_Store_Path,
        heapva->frags, // The internal representation of the fragment reads. 
        heapva->edges, // The internal representation of the overlaps.
        heapva->frag_annotations,
        afr_to_avx,
        &(gstate->min_frag_iid),
        &(gstate->max_frag_iid),
        heapva->next_edge_obj,
        &(gstate->nbase_in_genome),
        rg->dvt_double_sided_threshold_fragment_end_degree,
        rg->con_double_sided_threshold_fragment_end_degree,
        rg->intrude_with_non_blessed_overlaps_flag,
        rg->overlap_error_threshold
        );
    if(TIMINGS) {
      time(&tp2); 
      fprintf(stderr,"%10" F_TIME_TP " sec: Finished reading input file\n",
              (tp2-tp1));
      system_top();
    }
  }

  destroy_FragmentHash(afr_to_avx);

  /////////////////////////////////////////////////////////////////
  fprintf(stderr,"All data is in.\n");
  ReportHeapUsage_CGB( gstate, heapva, stderr); 

  if(NULL != rg->Output_Graph_Store) {
    fprintf(stderr,"Check point after all the data is in.\n");
    { 
      {	  
        char thePath1[CMD_BUFFER_SIZE-1]={0};
        char thePath2[CMD_BUFFER_SIZE-1]={0};
        int ierr=0;
        sprintf(thePath1,"%s/%s",rg->Output_Graph_Store,"fgb.ckp_tmp");
        sprintf(thePath2,"%s/%s",rg->Output_Graph_Store,"fgb.ckp_all_data_in");
        write_fgb_store(thePath1, gstate, heapva);
        ierr = accept_tmp_as_final_file( thePath1, thePath2);
        assert(ierr == 0);
      }
      if(rg->analysis_flag) {
        FILE *ffga=NULL;
        char thePath3[CMD_BUFFER_SIZE-1]={0};
        char thePath4[CMD_BUFFER_SIZE-1]={0};
        const int ProcessFragmentAnnotationsForSimulatorCoordinates
          = (rg->analysis_flag > 1);
        int ierr=0;
        sprintf(thePath3,"%s/%s",rg->Output_Graph_Store,"fga.ckp_tmp");
        sprintf(thePath4,"%s/%s",rg->Output_Graph_Store,"fga.ckp_all_data_in");
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
  }
  /////////////////////////////////////////////////////////////////
  

#ifdef DEBUG_VIEW1
    if( NULL != rg->Output_Graph_Store ) {
      char strtmp[CMD_BUFFER_SIZE-1];
      sprintf(strtmp,"%s.all_data_is_in",
	      rg->Output_Graph_Store);
      view_fgb_chkpnt( strtmp, heapva->frags, heapva->edges);
    }
#endif    
  
#if 0
  {
    FILE *before = fopen("next_edge_obj.before","w");
    FILE *after  = fopen("next_edge_obj.after","w");
    // Re-establish the next_edge[] pointers.
    fprintf(stderr,"Re-establish the next_edge[] pointers.\n");
    dump_next_edge_obj(before, heapva->next_edge_obj);
    // clear next_edge_obj
    Reset_VA(heapva->next_edge_obj);
    rebuild_next_edge_obj( heapva->edges, heapva->next_edge_obj);
    dump_next_edge_obj(after, heapva->next_edge_obj);
    fclose(before);
    fclose(after);
  }
#endif

  count_fragment_and_edge_labels
    ( heapva->frags, heapva->edges,
      "Before separate_fragments_as_solo_hanging_thru");
  separate_fragments_as_solo_hanging_thru(heapva->frags,heapva->edges);

#ifdef DEBUG_VIEW2
    if( NULL != rg->Output_Graph_Store ) {
      char strtmp[CMD_BUFFER_SIZE-1];
      sprintf(strtmp,"%s.separate_solo_hanging_thru",
	      rg->Output_Graph_Store);
      view_fgb_chkpnt( strtmp, heapva->frags, heapva->edges);
    }
#endif    

    set_compare_edge_function(compare_edge_weak);
    // Now do not distinguish the blessed overlaps.

    reorder_edges( heapva->frags, heapva->edges, heapva->next_edge_obj);
    count_fragment_and_edge_labels( heapva->frags, heapva->edges, "RISM_reorder_edges");
#ifdef DEBUG_RISM
    view_fgb_chkpnt( "RISM_reorder_edges", heapva->frags, heapva->edges);
#endif// DEBUG_RISM
    check_symmetry_of_the_edge_mates( heapva->frags, heapva->edges, heapva->next_edge_obj);
    
    append_the_edge_mates( heapva->frags, heapva->edges, heapva->next_edge_obj);
    // Currently the spur finding and micro-bubble smoothing code
    // needs the dovetail edge mates to exist.
    count_fragment_and_edge_labels( heapva->frags, heapva->edges, "RISM_append_the_edge_mates");
#ifdef DEBUG_RISM
    view_fgb_chkpnt( "RISM_append_the_edge_mates", heapva->frags, heapva->edges);
#endif // DEBUG_RISM
    check_symmetry_of_the_edge_mates( heapva->frags, heapva->edges, heapva->next_edge_obj);

    delete_duplicate_edges( heapva->frags, heapva->edges, heapva->next_edge_obj);
    count_fragment_and_edge_labels( heapva->frags, heapva->edges, "RISM_delete_duplicate_edges");
#ifdef DEBUG_RISM
    view_fgb_chkpnt( "RISM_delete_duplicate_edges", heapva->frags, heapva->edges);
#endif // DEBUG_RISM
    check_symmetry_of_the_edge_mates( heapva->frags, heapva->edges, heapva->next_edge_obj);
    // Currently the spur finding code needs the dovetail edge mates to
    // exist.

#ifdef DEBUG_VIEW3
  if( NULL != rg->Output_Graph_Store ) {
    char strtmp[CMD_BUFFER_SIZE-1];
    sprintf(strtmp,"%s.jake_3.fgb",
            rg->Output_Graph_Store);
    fprintf(stderr,"View graph is %s\n", strtmp);
    view_fgb_chkpnt( strtmp, heapva->frags, heapva->edges);
    check_symmetry_of_the_edge_mates( heapva->frags, heapva->edges, heapva->next_edge_obj);
  }
#endif    

  if( (GetNumFragments(heapva->frags) > 0) &&
      (GetNumEdges(heapva->edges) > 0) &&
      (rg->dechord_the_graph) &&
      (!did_processing_phase_2) 
      ) /* The FGB computations need a populated graph store. */ {

    if(TIMINGS) {
      tp1 = time(NULL); 
      fprintf(stderr,"Begin transitively inferable edge marking\n");
      system_top();
    }
    transitive_edge_marking
      ( gstate, heapva, heapva->frags, heapva->edges, heapva->next_edge_obj,
        rg->walk_depth,
        rg->cutoff_fragment_end_degree,
        rg->work_limit_per_candidate_edge,
        rg->iv_start,
        rg->analysis_flag,
        rg->Output_Graph_Store);
    
    // Needs a global thread synchronization here.
    if(TIMINGS) {
      tp2 = time(NULL);
      fprintf(stderr,
              "%10" F_TIME_TP " sec: Finished transitively inferable edge marking\n",
              (tp2-tp1));
      system_top();
    }
  } // Finish processing the input.
  
  identify_early_spur_fragments( heapva->frags, heapva->edges);
#ifdef DEBUG_VIEW4
  if( NULL != rg->Output_Graph_Store ) {
    char strtmp[CMD_BUFFER_SIZE-1];
    sprintf(strtmp,"%s.identify_early_spur_fragments",
            rg->Output_Graph_Store);
    fprintf(stderr,"Last graph is %s\n", strtmp);
    view_fgb_chkpnt( strtmp, heapva->frags, heapva->edges);
    check_symmetry_of_the_edge_mates( heapva->frags, heapva->edges, heapva->next_edge_obj);
  }
#endif
  
  count_fragment_and_edge_labels
    ( heapva->frags, heapva->edges,
      "Before contained_fragment_marking_frc");
  contained_fragment_marking_frc( heapva->frags, heapva->edges);
  count_fragment_and_edge_labels
    ( heapva->frags, heapva->edges,
      "After contained_fragment_marking_frc");
  

  if( ! rg->create_dump_file ) {
    processing_phase_3
      (
       gstate, heapva,
       rg->Output_Graph_Store,
       rg->analysis_flag
       );
  }
    
  /**************** Finish Process *********************/
  
#ifdef DEBUG_VIEW5
  if(NULL != rg->Output_Graph_Store ) {
    char strtmp[CMD_BUFFER_SIZE-1];
    sprintf(strtmp,"%s.reaper",rg->Output_Graph_Store);
    fprintf(stderr,"Last graph is %s\n", strtmp);
    view_fgb_chkpnt( strtmp, heapva->frags, heapva->edges);
    check_symmetry_of_the_edge_mates( heapva->frags, heapva->edges, heapva->next_edge_obj);
  }
#endif    

  if( rg->create_dump_file ) {
    FILE *folp = NULL;
    char strtmp[CMD_BUFFER_SIZE-1];
    
    fprintf(stderr,"Opening dump file to write a batch of "
	    "ADT+IDT+OFG+OVL messages.\n");

    if(NULL == (folp = fopen(rg->Dump_File_Name,"w"))){
      fprintf(stderr,"* Can not open output file %s\n",rg->Dump_File_Name);
      exit(1);
    }

    output_mesgs
      (/* Input Only*/
       (heapva->frags), (heapva->edges),
       (heapva->frag_annotations),
       folp);
    ierr = fclose(folp); assert(ierr == 0);
  }
  
#ifdef DEBUG_VIEW6
    if( NULL != rg->Output_Graph_Store ) {
      char strtmp[CMD_BUFFER_SIZE-1];
      sprintf(strtmp,"%s.fgb_main",
	      rg->Output_Graph_Store);
      view_fgb_chkpnt( strtmp, heapva->frags, heapva->edges);
      check_symmetry_of_the_edge_mates( heapva->frags, heapva->edges, heapva->next_edge_obj);
    }
#endif    

  system_date();
  return( status);
}
