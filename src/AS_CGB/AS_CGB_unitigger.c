
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

static char CM_ID[] = "$Id: AS_CGB_unitigger.c,v 1.24 2007-08-04 22:26:58 brianwalenz Exp $";

#include "AS_UTL_version.h"
#include "AS_CGB_all.h"
#include "AS_CGB_Bubble.h"

extern int REAPER_VALIDATION;

void
chunk_graph_analysis(THeapGlobals *heapva,
                     UnitiggerGlobals *rg);


void
output_the_chunks(Tfragment     *frags,
                  Tedge         *edges,
                  TChunkFrag    *chunkfrags,
                  TChunkMesg    *thechunks,
                  float          global_fragment_arrival_rate,
                  int            fragment_count_target,
                  char          *Graph_Store_File_Prefix) {

  IntChunk_ID       chunk_index;
  IntChunk_ID       nchunks = (IntChunk_ID)GetNumVA_AChunkMesg(thechunks);
  IntFragment_ID    nfrag   = GetNumFragments(frags);
  
  char filename[FILENAME_MAX];
  
  FILE *fcgb = NULL;
  int fragment_count = 0;
  int file_count = 0;

  VA_TYPE(IntMultiPos) * the_imps = CreateVA_IntMultiPos(0);
  
  for(chunk_index=0;chunk_index < nchunks; chunk_index++) /* a */ {
    AChunkMesg     *ch   = GetVA_AChunkMesg(thechunks,chunk_index);
    IntFragment_ID  num_frags = ch->num_frags;
    GenericMesg     pmesg;
    IntUnitigMesg   achunk;
    IntFragment_ID  ivc;

    for(ivc=0; ivc<num_frags; ivc++) {
      IntFragment_ID vid    = *GetVA_AChunkFrag(chunkfrags, ch->f_list + ivc);
      IntMultiPos    a_frag = {0};

      a_frag.type         = get_typ_fragment(frags,vid);
      a_frag.ident        = get_iid_fragment(frags,vid);
      a_frag.contained  = get_container_fragment(frags,vid);
      a_frag.position.bgn = get_o5p_fragment(frags,vid);
      a_frag.position.end = get_o3p_fragment(frags,vid);
      a_frag.delta_length = 0;
      a_frag.delta        = NULL;

      SetVA_IntMultiPos(the_imps,ivc,&a_frag);
    }

    achunk.consensus = "";
    achunk.quality   = "";
      
    achunk.iaccession     = ch->iaccession;
    achunk.source         = "gen> @@ [0,0]";
    achunk.coverage_stat  = compute_coverage_statistic(ch->rho,
                                                       count_the_randomly_sampled_fragments_in_a_chunk(frags,
                                                                                                       chunkfrags,
                                                                                                       thechunks,
                                                                                                       chunk_index),
                                                       global_fragment_arrival_rate);
    achunk.status         = AS_UNASSIGNED;

    achunk.consensus      = "";
    achunk.quality        = "";
    achunk.length         = ch->bp_length;

    achunk.forced         = FALSE;
    achunk.num_frags      = ch->num_frags;
    achunk.f_list         = GetVA_IntMultiPos(the_imps,0);
    achunk.num_vars       = 0;

    fragment_count += ch->num_frags;

    if(fragment_count_target > 0) {
      if((fragment_count >= fragment_count_target)) {
        if(NULL != fcgb) { fclose(fcgb); fcgb = NULL;}
        fragment_count = 0;
      }
      if(NULL == fcgb) {
        file_count ++;
        sprintf(filename,"%s_%03d.cgb",Graph_Store_File_Prefix,file_count);
        fcgb = fopen(filename,"w");
        assert(NULL != fcgb);
      }
    } else {
      if(NULL == fcgb) {
        sprintf(filename,"%s.cgb",Graph_Store_File_Prefix);
        fcgb = fopen(filename,"w");
        assert(NULL != fcgb);
      }
    }
    
    pmesg.t = MESG_IUM;
    pmesg.m = &achunk;
    WriteProtoMesg_AS(fcgb,&pmesg);
  }

  DeleteVA_IntMultiPos(the_imps);

  fclose(fcgb);
}






int
ParseCommandLine(UnitiggerGlobals * rg,
                 int argc,
                 char * argv[]) {

  int illegal=FALSE;
  int blessed_overlaps_output_flag = FALSE;
  
  int ch,errflg=0;
  optarg = NULL;

  argc = AS_configure(argc, argv);

  while (!errflg && 
         ((ch = getopt(argc, argv, 
                       "B:F:H:I:L:S:U:W:Y:"
                       "d:e:h:j:kl:m:n:o:p:su:w:x:y:z:"
                       "56:7"
                       )) != EOF)) {

    switch(ch) {
      /* The required command line options: */
      case 'B':
        // -B <int> : Specifies a target number of fragments in IMP
        // records per cgb file.
        rg->fragment_count_target = atoi(optarg);
        break;
      case 'F':
        // -F <path> : identify fragment store
        rg->frag_store = optarg;
        break;
      case 'H':
        // -H <filename> : identify chimeras file
        rg->chimeras_file = optarg;
        fprintf( stderr, " * chimeras file is %s\n", rg->chimeras_file );
        break;
      case 'I':
        // -I <file> : The input overlapper overlap store.
        rg->OVL_Store_Path = optarg;
        fprintf(stderr,"* The input overlapper overlap store is <%s>.\n",
                rg->OVL_Store_Path);
        break;
      case 'L':
        rg->ovl_files_list_fname = optarg;
        fprintf(stderr,"*  ovl_files_list_fname <%s>\n",
                rg->ovl_files_list_fname);
        break;
      case 'S':
        // -S <filename> : identify spurs file
        rg->spurs_file = optarg;
        fprintf( stderr, " * spurs file is %s\n", rg->spurs_file );
        break;
      case 'U':
        // -U 
        rg->bubble_smoothing_flag = atoi(optarg);
        fprintf(stderr, "* Bubble smoothing is now %s.\n", (rg->bubble_smoothing_flag) ? "ON" : "OFF");
        break;
      case 'W':
        // -W int : The maximum path walk depth to allow during
        // transitive edge graph reduction.
        rg->walk_depth = atoi(optarg);
        fprintf(stderr,"* walk_depth set to %d\n",
                rg->walk_depth);
        break;
      case 'Y':
        // -Y 
        rg->dont_count_chimeras = atoi(optarg);
        break;

      case 'd':
        // -d : De-chord the fragment overlap graph.
        rg->dechord_the_graph = atoi(optarg);
        fprintf(stderr,"* Graph dechording is %s.\n",
                (rg->dechord_the_graph ? "on" : "off"));
        break;
      case 'e':
        // -e <float or int> : Overlaps with error rates above this
        // value will be ignored on input.
        {
          int float_erate_value = (NULL != strstr(optarg,"."));
          if(float_erate_value) {
            rg->overlap_error_threshold = AS_OVS_encodeQuality(atof(optarg) / 100.0);
          } else {
            rg->overlap_error_threshold = AS_OVS_encodeQuality(atof(optarg) / 1000.0);
          }
          fprintf(stderr,"The overlap_error_threshold = %f%%\n",
                  AS_OVS_decodeQuality(rg->overlap_error_threshold) * 100.0);
        }
        break;
      case 'h':
        // help
        goto UsageStatement;
        // Usage( argv[0], NULL );
        break;
      case 'j':
        // -j <float> : Astat cut-off threshold
        rg->cgb_unique_cutoff = atof(optarg);
        fprintf(stderr,"* cgb_unique_cutoff set to %f\n", rg->cgb_unique_cutoff);
        break;
      case 'k':
        // -k : Recalibrate the global arrival rate to be the max unique local arrival rate
        rg->recalibrate_global_arrival_rate = TRUE;
        fprintf(stderr,
                "* -k: Recalibrate the global arrival rate to be the max unique local arrival rate.\n");
        break;
      case 'l':
        // -l <int> : the length of the genome
        rg->genome_length = STR_TO_INT64(optarg, NULL, 10);
        break;
      case 'm':
        // -m <int> : Pre-allocate space to process this many additional
        // overlaps.
        rg->maxedges = atol(optarg);
        break;
      case 'n':
        // -n <int> : Pre-allocate space to process this many additional
        // fragments.
        rg->maxfrags = atol(optarg);
        break;

      case 'o':
        rg->Output_Graph_Store_Prefix = (char *)optarg;
        break;

      case 's':
        // -s 
        rg->aggressive_spur_fragment_marking = FALSE;
        // Aggressive removal of the early spurs.
        break;

      case 'u':
        // -u <file> : Create a OVL file compatible dump of the fragment
        // graph store.
        rg->Dump_File_Name = optarg;
        fprintf(stderr,"* OVL dump file set to <%s>.\n", rg->Dump_File_Name);
        rg->create_dump_file = (rg->Dump_File_Name[0] != '\0');
        break;

      case 'w':
        // -w <int> : The work limit per candidate edge for
        // de-chording.
        rg->work_limit_per_candidate_edge = atoi(optarg);
        fprintf(stderr,
                "* work_limit_per_candidate_edge set to %d\n",
                rg->work_limit_per_candidate_edge);
        break;
      case 'x':
        // -x <int> : The threshold adjaceny list position for "Double
        // sided thresholding".
        rg->dvt_double_sided_threshold_fragment_end_degree = atoi(optarg);
        fprintf(stderr,
                "* dvt_double_sided_threshold_fragment_end_degree set to %d\n",
                rg->dvt_double_sided_threshold_fragment_end_degree);
        break;
      case 'y':
        // -y <int> : Intrude with non-blessed overlap edges to blessed fragment ends.
        rg->intrude_with_non_blessed_overlaps_flag = atoi(optarg);
        fprintf(stderr,
                "* Intrude with non-blessed overlap edges to blessed fragment ends %d\n",
                rg->intrude_with_non_blessed_overlaps_flag);
        break;
      case 'z':
        // -z <int> : The threshold adjaceny list position for "Containment degree
        // trimming".
        rg->con_double_sided_threshold_fragment_end_degree = atoi(optarg);
        fprintf(stderr,
                "* con_double_sided_threshold_fragment_end_degree set to %d\n",
                rg->con_double_sided_threshold_fragment_end_degree);
        break;

      case '5':
        REAPER_VALIDATION = TRUE;
        break;
      case '6':
        rg->blessed_overlaps_input_filename = optarg;
        break;
      case '7':
        blessed_overlaps_output_flag = TRUE;
        break;
      default :
        fprintf(stderr,"Unrecognized option -%c\n",optopt);
        errflg++;
        illegal = TRUE;
        goto UsageStatement;
    }
  }

  // Get the list of ovl files on the command line.
  rg->num_ovl_files = argc - optind;
  rg->the_ovl_files = &(argv[optind]);


  if(rg->bubble_smoothing_flag && (NULL == rg->frag_store)) {
    fprintf(stderr,"Error: bubble smoothing needs a fragment store.\n");
    exit(1);
  }

  if(! rg->dont_count_chimeras) {
    static char chimeras_report_filename[128]={0};
    static char crappies_report_filename[128]={0};

    sprintf(crappies_report_filename,"%s.cgb_crappies", rg->Output_Graph_Store_Prefix);
    rg->spurs_file = crappies_report_filename;

    sprintf(chimeras_report_filename,"%s.cgb_chimeras", rg->Output_Graph_Store_Prefix);
    rg->chimeras_file = chimeras_report_filename;
  }

  if(blessed_overlaps_output_flag) {
    static char blessed_overlaps_output_filename[128]={0};
    sprintf(blessed_overlaps_output_filename,"%s.blovl", rg->Output_Graph_Store_Prefix);
    rg->blessed_overlaps_output_filename = blessed_overlaps_output_filename;
  }

 UsageStatement:
  if((illegal == TRUE)) {
    fprintf(stderr, "USAGE: %s <option>*\n"
            "\t-B <int>        Specifies the target number of fragments per partition.\n"
            "\t-F <directory>  The fragment store name.\n"
            "\t-H <filename>   chimeras file.\n"
            "\t-I <directory>  Read the OVL store.\n"
            "\t-L <filename>   The input OverlapFragMesgs; asm.ofg.\n"
            "\t-S <filename>   Spurs file.\n"
            "\t-U <boolean>    Find bubble smoothing overlaps.\n"
            "\t-W <int>        Limit in path length for graph walking.\n"
            "\t-Y <boolean>    Do not count chimera fragments.\n"
            "\t-d <int>        Enable/Disable de-chording of the fragment overlap graph.\n"
            "\t-e <n>          Overlaps with error rate about this are ignored on input.\n"
            "\t\t                An integer value is in parts per thousand.\n"
            "\t-h              Help.\n"
            "\t-j <int>        Unique unitig cut-off\n"
            "\t-k              Recalibrate the global arrival rate to be the max unique local arrival rate\n"
            "\t-l <int>        Specify length of the genome.\n"
            "\t-m <nEdge>      Pre-allocate memory\n"
            "\t-n <nFrag>      Pre-allocate memory\n"
            "\t-o <pfx>        output to this prefix.\n"
            "\t-s              Disable early spur fragment removal.\n"
            "\t-u <filename>   Create a OVL compatible dump of the graph.\n"
            "\t-w <int>        The work limit per candidate edge for de-chording.\n"
            "\t-x <int>        Dovetail outgoing degree threshold per fragment-end.\n"
            "\t-y <int>        Ignore non-blessed overlap edges to blessed fragment ends.\n"
            "\t-z <int>        Containment outgoing degree threshold per fragment-end.\n"
            ,
            argv[0]);
    exit (EXIT_FAILURE);
  }

  /* The remaining command line arguments are input Assembler
     message files. */

  return 0;
}




int
main(int argc, char **argv) {
  THeapGlobals     *heapva = (THeapGlobals      *)safe_calloc(sizeof(THeapGlobals), 1);
  UnitiggerGlobals *rg     = (UnitiggerGlobals  *)safe_calloc(sizeof(UnitiggerGlobals), 1);

  rg->work_limit_per_candidate_edge = 1000;
  rg->dvt_double_sided_threshold_fragment_end_degree = 1;
  rg->con_double_sided_threshold_fragment_end_degree = 1;
  rg->intrude_with_non_blessed_overlaps_flag = FALSE;
  rg->cutoff_fragment_end_degree = 1000000;
  rg->dechord_the_graph = TRUE;
  rg->cgb_unique_cutoff = CGB_UNIQUE_CUTOFF;
  rg->walk_depth = 100;
  rg->overlap_error_threshold = AS_OVS_encodeQuality(1.0);
  rg->recalibrate_global_arrival_rate = FALSE;
  rg->walk_depth=100;
  rg->output_iterations_flag = TRUE;
  rg->aggressive_spur_fragment_marking = TRUE;

  rg->maxfrags = 40000;
  rg->maxedges = 40000;
  rg->maxtext  = 40000;

  ParseCommandLine( rg, argc, argv);

  //BasicUnitigger( argc, argv, gstate, heapva, rg);

 again:
  heapva->frags             = CreateVA_Afragment (rg->maxfrags);
  heapva->edges             = CreateVA_Aedge     (rg->maxedges); 
  heapva->chunkfrags        = CreateVA_AChunkFrag(0);
  heapva->thechunks         = CreateVA_AChunkMesg(0);
  heapva->nbase_in_genome              = 0;
  heapva->global_fragment_arrival_rate = 0;

  main_fgb(heapva, rg);
  main_cgb(heapva, rg);

  //  End of BasicUnitigger

  //  Do bubble smoothing if enabled.
  if(rg->bubble_smoothing_flag) {

    //  Turn this on, else we'll nuke the crapppies files when we redo
    //  work.
    rg->dont_count_chimeras = TRUE;

    //  Turn this off so we don't do smoothing again.
    rg->bubble_smoothing_flag = FALSE;

    /* BUBBLE SMOOTHING CONFIGURATION PARAMETERS - See also
       AS_CGB_Bubble.h           - Flag to enable debugging output, 
       AS_CGB_Bubble_Popper.h    - Alignment parameters. 
       AS_CGB_Bubble_VertexSet.h - Default parameters for bubble finding. */

    AS_CGB_Bubble_List_t bubbles = NULL;

    // The OVL records needed to remove the bubbles
    sprintf(rg->bubble_overlaps_filename, "%s.bubble_edges.ovl", rg->Output_Graph_Store_Prefix);
    FILE *bfp = fopen(rg->bubble_overlaps_filename, "w");

    GateKeeperStore *gkpStore = openGateKeeperStore(rg->frag_store, FALSE);
    AS_CGB_Bubble_find_and_remove_bubbles(gkpStore,
                                          heapva->frags, heapva->edges,
                                          heapva->thechunks, heapva->chunkfrags, 
                                          heapva->global_fragment_arrival_rate,
                                          bfp,
                                          stderr,
                                          rg->Output_Graph_Store_Prefix);
    closeGateKeeperStore(gkpStore);
    fclose(bfp);

    /* NOTE: 0's in following call indicate use of defaults. */
  
    bubbles = AS_CGB_Bubble_find_bubbles(heapva->frags, heapva->edges, 0, 0, 0);    

    {
      char bfn[500] = {0};
      sprintf(bfn,"%s.cgb_bubbles_txt",rg->Output_Graph_Store_Prefix);
      FILE *bfp = fopen(bfn, "w");
      AS_CGB_Bubble_List_t bubbles_tmp = bubbles;

      while (bubbles_tmp != NULL) {
        fprintf(bfp, F_IID" %d "F_IID" %d\n", 
                get_iid_fragment(heapva->frags, bubbles_tmp->start),
                bubbles_tmp->start_sx,
                get_iid_fragment(heapva->frags, bubbles_tmp->end),
                bubbles_tmp->end_sx);
        bubbles_tmp = bubbles_tmp->next;
      }
      fclose(bfp);
    }

    //  Like I said, redo all our work.
    Delete_VA(heapva->frags);
    Delete_VA(heapva->edges);
    Delete_VA(heapva->chunkfrags);
    Delete_VA(heapva->thechunks);

    goto again;
  }

  chunk_graph_analysis(heapva, rg);

  output_the_chunks(heapva->frags,
                    heapva->edges,
                    heapva->chunkfrags,
                    heapva->thechunks,
                    heapva->global_fragment_arrival_rate,
                    rg->fragment_count_target,
                    rg->Output_Graph_Store_Prefix);



  // Determine the blessed overlaps.  They are the overlaps
  // that are interior to u-unitigs.  In particular, intrachunk overlaps 
  //
  IntFragment_ID nfrag = GetNumFragments(heapva->frags);
  IntEdge_ID     nedge = GetNumEdges(heapva->edges);
  IntEdge_ID     ie;

  FILE          *fcgb = NULL;

  if (rg->blessed_overlaps_output_filename)
    fcgb = fopen(rg->blessed_overlaps_output_filename,"w");

  GenericMesg    pmesg;
  OverlapMesg    omesg;

  pmesg.t = MESG_OVL;
  pmesg.m = &omesg;

  for (ie=0; ie < nedge; ie ++) {
    Tnes nes = get_nes_edge(heapva->edges,ie);
    IntFragment_ID avx = get_avx_edge(heapva->edges,ie);

    assert(nfrag > avx);

    AChunkMesg *ch = GetVA_AChunkMesg(heapva->thechunks, get_cid_fragment(heapva->frags,avx));

    if(ch->coverage_stat >= rg->cgb_unique_cutoff) {
      if (AS_CGB_INTRACHUNK_EDGE == nes )
        set_blessed_edge(heapva->edges, ie, TRUE);

      if ((AS_CGB_SINGLECONT_FRAG == get_lab_fragment(heapva->frags,avx)) &&
          (AS_CGB_CONTAINED_EDGE == nes) &&
          (is_a_frc_edge(heapva->edges,ie)))
        set_blessed_edge(heapva->edges, ie, TRUE);
    }

    // output latest set of blessed overlap edges.

    if ((fcgb) && (get_blessed_edge(heapva->edges,ie))) {
      IntFragment_ID avx = get_avx_edge(heapva->edges,ie);
      int asx = get_asx_edge(heapva->edges,ie);
      int ahg = get_ahg_edge(heapva->edges,ie);
        
      IntFragment_ID bvx = get_bvx_edge(heapva->edges,ie);
      int bsx = get_bsx_edge(heapva->edges,ie);
      int bhg = get_bhg_edge(heapva->edges,ie);
        
      uint32 qua = get_qua_edge(heapva->edges,ie);

      IntFragment_ID aid = get_iid_fragment(heapva->frags,avx);
      IntFragment_ID bid = get_iid_fragment(heapva->frags,bvx);
        
      signed char delta[1] = {AS_ENDOF_DELTA_CODE};

      omesg.aifrag = aid;
      omesg.bifrag = bid;
          
      omesg.ahg        = ahg;
      omesg.bhg        = bhg;
      omesg.min_offset = ahg;
      omesg.max_offset = ahg;
          
      omesg.orientation = asx ? (bsx ? AS_INNIE : AS_NORMAL) :
        (bsx ? AS_ANTI  : AS_OUTTIE);

      switch (get_nes_edge(heapva->edges,ie)) {
        case AS_CGB_DOVETAIL_EDGE:
        case AS_CGB_THICKEST_EDGE:
        case AS_CGB_INTERCHUNK_EDGE:
        case AS_CGB_INTRACHUNK_EDGE:
        case AS_CGB_TOUCHES_CONTAINED_EDGE:
        case AS_CGB_BETWEEN_CONTAINED_EDGE:
        case AS_CGB_TOUCHES_CRAPPY_DVT:
          omesg.overlap_type = AS_DOVETAIL;
          break;
        case AS_CGB_CONTAINED_EDGE:
        case AS_CGB_TOUCHES_CRAPPY_CON:
        case AS_CGB_BETWEEN_CRAPPY_CON:
          omesg.overlap_type = AS_CONTAINMENT;
          break;
        default:
          fprintf(stderr,"Unexpected overlap edge type: nes=%d\n", nes);
          assert(FALSE);
      }

      omesg.quality = AS_OVS_decodeQuality(qua);
      omesg.polymorph_ct = 0;
      omesg.delta = delta;
          
      WriteProtoMesg_AS(fcgb,&pmesg);
    }
  }

  if (fcgb)
    fclose(fcgb);

  //  This writes the .fgv and .fge files.
  view_fgb_chkpnt(rg->Output_Graph_Store_Prefix,
                  heapva->frags,
                  heapva->edges);

  
  Delete_VA(heapva->frags);
  Delete_VA(heapva->edges);
  Delete_VA(heapva->chunkfrags);
  Delete_VA(heapva->thechunks);

  fprintf(stderr, "All done.  Bye.\n");

  exit(0);
}
