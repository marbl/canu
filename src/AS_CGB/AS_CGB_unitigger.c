
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
= "$Id: AS_CGB_unitigger.c,v 1.8 2007-01-28 21:52:24 brianwalenz Exp $";
/*********************************************************************
 *
 * Module: AS_CGB_unitigger.c
 * Description:
 *    Unitigger main
 * 
 *    Reference:
 *
 *    Command Line Interface:
 *
 *    Programmer:  C. M. Mobarry
 *       Written:  Dec 2001
 * 
 *********************************************************************/

// "Mini-unitigger"
// unitigger -P -A 2 -Y 1 -N 1000 -x 1 -z 10 -j 5 -c -a
// ../Test11/a006l.ovl
// 1> nohup.unitigger.1.wb 2> nohup.unitigger.2.wb

// cd /work/assembly/cmobarry/trunk/data/a006/Test12
// prefix=ax; rm -rf $prefix.* ;

// unitigger -P -A 2 -Y 1 -N 1000 -x 1 -z 10 -j 5 -U 1 -F
// ../Test11/a006l.frgStore -c -o $prefix.fgbStore ../Test11/a006l.ovl

// unitigger -P -A 2 -N 1000 -x 1 -z 10 -j 5 -s -S
// ${prefix}.cgb_crappies -F ../Test11/a006l.frgStore -O
// ${prefix}.breaker_overlaps.ovl -c -o $prefix.fgbStore
// ../Test11/a006l.ovl

// unitigger -P -A 2 -Y 1 -N 1000 -n 1000 -j 5
//   -F xy.frgStore -f -c -o xy4.fgbStore -I xy.ovlStore -L xy.ofgList
//   1> xy4.stdout 2> xy4.stderr

#include "AS_CGB_all.h"
#include "AS_UTL_version.h"
#include "AS_CGB_unitigger_globals.h"
#include "AS_CGB_Bubble.h"
#ifdef REPAIR_BREAKERS
#include "AS_CGB_repair_breakers.h"
#endif // REPAIR_BREAKERS


extern int REAPER_VALIDATION;
static int TIMINGS = TRUE;
static int skip_fgb = FALSE;
static int skip_cgb = FALSE;

/****************************************************************************/

static void output_OVL_mesgs
(/* Input Only*/
 MesgWriter WriteMesg_AS,
 const Tfragment frags[],
 const Tedge     edges[],
 const VA_TYPE(char) fragsrc[],
 /* Read Only */
 FILE *filk,
 /* Append Only*/
 FILE *fcgb)
{
  
  assert(NULL != frags);
  assert(NULL != edges);
  
  // Output the OVL messages:
  if(NULL != fcgb) {
    const IntEdge_ID  nedge = GetNumEdges(edges);
    IntEdge_ID ie;
    for(ie=0;ie<nedge;ie++){
      const int blessed = get_blessed_edge(edges,ie);
      if( blessed ) {

        const IntFragment_ID avx = get_avx_edge(edges,ie);
        const int asx = get_asx_edge(edges,ie);
        const int ahg = get_ahg_edge(edges,ie);
#ifdef STORE_OVERLAP_EXTREMES
        const int amn = get_amn_edge(edges,ie);
        const int amx = get_amx_edge(edges,ie);
#endif // STORE_OVERLAP_EXTREMES
        
        const IntFragment_ID bvx = get_bvx_edge(edges,ie);
        const int bsx = get_bsx_edge(edges,ie);
        const int bhg = get_bhg_edge(edges,ie);
#ifdef STORE_OVERLAP_EXTREMES
        const int bmn = get_bmn_edge(edges,ie);
        const int bmx = get_bmx_edge(edges,ie);
#endif // STORE_OVERLAP_EXTREMES
        
        const Tnes  nes = get_nes_edge(edges,ie);
        const CGB_ERATE_TYPE qua = get_qua_edge(edges,ie);
        
        const IntFragment_ID aid = get_iid_fragment(frags,avx);
        const IntFragment_ID bid = get_iid_fragment(frags,bvx);
        // Assembler internal Fragment ids.
        
        signed char delta[1] = {AS_ENDOF_DELTA_CODE};

        {
          OverlapMesg ovl_mesg;
          ovl_mesg.aifrag = aid;
          ovl_mesg.bifrag = bid;
          
          ovl_mesg.ahg = ahg;
          ovl_mesg.bhg = bhg;
#ifndef STORE_OVERLAP_EXTREMES
          ovl_mesg.min_offset = ahg;
          ovl_mesg.max_offset = ahg;
#else // STORE_OVERLAP_EXTREMES
          ovl_mesg.min_offset = amn;
          ovl_mesg.max_offset = amx;
#endif // STORE_OVERLAP_EXTREMES
          
          ovl_mesg.orientation =
            ( asx ?
              ( bsx ? AS_INNIE : AS_NORMAL ) :
              ( bsx ? AS_ANTI  : AS_OUTTIE ) );
          
          switch(nes) {
          case AS_CGB_DOVETAIL_EDGE:
          case AS_CGB_THICKEST_EDGE:
          case AS_CGB_INTERCHUNK_EDGE:
          case AS_CGB_INTRACHUNK_EDGE:
          case AS_CGB_TOUCHES_CONTAINED_EDGE:
          case AS_CGB_BETWEEN_CONTAINED_EDGE:
          case AS_CGB_TOUCHES_CRAPPY_DVT:
            ovl_mesg.overlap_type = AS_DOVETAIL; break;
          case AS_CGB_CONTAINED_EDGE:
          case AS_CGB_TOUCHES_CRAPPY_CON:
          case AS_CGB_BETWEEN_CRAPPY_CON:
            ovl_mesg.overlap_type = AS_CONTAINMENT; break;
          default:
            fprintf(stderr,"Unexpected overlap edge type: nes=%d\n", nes);
            assert(FALSE);
          }

          ovl_mesg.quality = CGB_ERATE_TYPE_to_cds_float32(qua);
          ovl_mesg.polymorph_ct = 0;
          ovl_mesg.delta = delta;
          
          if(blessed){
            GenericMesg   pmesg;
            pmesg.t = MESG_OVL;
            pmesg.m = &ovl_mesg;
            WriteMesg_AS(fcgb,&pmesg);
          }
        }
      }
    }
  }
  
#if 0  
  if(analysis_flag) {
    fprintf(stderr,"\n\nHistogram of the OVL types\n");
    print_histogram(stderr, ovl_types_histogram, 0, 1);
  }
  if(NULL != ovl_types_histogram) {
    free_histogram(fom_types_histogram);
  }
#endif
}

#ifdef NOT_USED

static FILE *  File_Open
(const char * Filename, const char * Mode)
     /* Open  Filename  in  Mode  and return a pointer to its control
      *  block.  If fail, print a message and exit. */
{
  FILE  *  fp;
  fp = fopen (Filename, Mode);
  if  (fp == NULL)
    {
      fprintf (stderr, "ERROR:  Could not open file  %s \n", Filename);
      exit (EXIT_FAILURE);
    }
  return  fp;
}

static int File_Exists (const char * Filename){
  /* Test for filename existence */
  FILE  *  fp;
  fp = fopen (Filename, "r+");
  if(fp) {
    fclose(fp);
    return 1;
  }
  return  0;
}

#endif // NOT_USED

#if 0
int nerrs = 0;   // Number of errors in current run
int maxerrs = 10; // Number of errors allowed before we punt

static int incrementErrors(int num, FILE *msgFile){
  //  fprintf(msgFile,"* incrementErrors by %d -- (%d,%d)\n",
  //  num, nerrs, maxerrs);
  nerrs += num;
  if(nerrs >= maxerrs){
    fprintf(msgFile, "Maximum allowed errors reached %d > %d...bye\n",
	    nerrs, maxerrs);
    return(EXIT_FAILURE);
  }
  //fprintf(stderr,"* incrementErrors returning SUCCESS\n");
  return(EXIT_SUCCESS);
}
#endif


/****************************************************************************/
#define ANALYSIS_EXTENSION     ".cga"

static void save_the_chunk_graph
(
 MesgWriter WriteMesg_AS,
 const char * const Graph_Store_File_Prefix,
 const char * const Graph_Store_File_Name,
 const int analysis_level,
 const char * bubble_boundaries_filename,
 const int output_fom_messages,
 const float cgb_unique_cutoff,
 const int iterator,
 const int fragment_count_target,
 int argc,
 char * argv [],
 TStateGlobals * gstate,
 THeapGlobals  * heapva
 )
{
  FILE *fcga = NULL;
  FILE *fcam = NULL;
  FILE *fwrn = stderr;
  FILE *fp_unitig_statistics = NULL;

  char strtmp2[CMD_BUFFER_SIZE-1];
  int ierr;
  time_t tp1 = 0, tp2;

  if(analysis_level > 0) {
    if(NULL != Graph_Store_File_Prefix) {
      sprintf(strtmp2,"%s%s.%d",Graph_Store_File_Prefix,
	      ANALYSIS_EXTENSION,iterator);
      SAFE_FOPEN( fcga, strtmp2, "w");
      sprintf(strtmp2,"%s%s.%d",Graph_Store_File_Prefix,".cam",iterator);
      SAFE_FOPEN( fcam, strtmp2, "w");
      sprintf(strtmp2,"%s%s.%d",Graph_Store_File_Prefix,".cus",iterator);
      SAFE_FOPEN( fp_unitig_statistics, strtmp2, "w");
    }

    if(TIMINGS) {
      tp1 = time(NULL); fprintf(stderr,"Begin chunk graph analysis\n");
      system_date();
    }
    chunk_graph_analysis
      (
       /* Input Only */
       analysis_level,
       gstate->max_frag_iid,
       heapva->frags,  /* The internal representation of
                           the fragment reads. */
       heapva->edges,  /* The internal representation of the  overlaps. */
       heapva->frag_annotations,
       gstate->nbase_in_genome,
       cgb_unique_cutoff,
       gstate->global_fragment_arrival_rate,
       bubble_boundaries_filename,
       heapva->chunkfrags,
               
       /* Modify the chunk annotation, by setting the "isrc" index into 
          chunksrc array. */
       heapva->thechunks,
       heapva->chunksrcs, 
       // The character handle used to store the annotations.
       fcga, fcam, fp_unitig_statistics, fwrn );

    SAFE_FCLOSE(fcga);
    SAFE_FCLOSE(fcam);
    SAFE_FCLOSE(fp_unitig_statistics);
    
    if(TIMINGS) {
      tp2 = time(NULL); 
      fprintf(stderr,"%10" F_TIME_TP " sec: Finished chunk graph analysis\n",
              (tp2-tp1));
      system_date();
    }
  }

  /* unitigger output */ {  

    if(TIMINGS) { 
      tp1 = time(NULL); fprintf(stderr,"Begin chunk graph message I/O\n");
      system_date();
    }

    {
      FILE *fcgb = NULL;
      char strtmp2[CMD_BUFFER_SIZE-1];
      
      fprintf(stderr,"Openning output file to append the "
              "IUM+UOM messages.\n");
      strcpy(strtmp2,Graph_Store_File_Prefix);

      if(fragment_count_target == 0) {
          sprintf(strtmp2,"%s.cgb_tmp",Graph_Store_File_Prefix);
      } else {
        int file_count = 0;
          sprintf(strtmp2,"%s_%03d.cgb_tmp",
		  Graph_Store_File_Prefix,file_count);
      }

      SAFE_FOPEN( fcgb, strtmp2, "w");

      {
        // Output the snapshot ADT message.
        GenericMesg mesg;
        AuditMesg  adt_mesg;
        
        /*  Audit Message record */
        mesg.t = MESG_ADT;
        mesg.m = &adt_mesg;
        mesg.s = sizeof(AuditMesg);
        adt_mesg.list = NULL;
        VersionStampADT(&adt_mesg, argc, argv);
        WriteMesg_AS(fcgb,&mesg);
      }
      
#ifdef USE_IBA_FILE
      {
        FILE *fiba = NULL;
        fprintf(stderr,"Openning intermediate file to read the "
                "IBA+ADT(+ADL)+IDT messages.\n");
        sprintf(strtmp2,"%s/fgb.iba",Graph_Store_File_Name);
        if(NULL == (fiba = fopen(strtmp2,"r"))){
          fprintf(stderr,"*WARNING: Can not open fgb.iba file %s\n",strtmp2);
          // exit(1);
        }
        
        if(NULL != fiba) {
          // Output the ADT messages from all the incremental processing.
          MesgReader ReadMesg_AS = (MesgReader)InputFileType_AS(fiba);
          GenericMesg *pmesg = NULL;
          while(EOF != ReadMesg_AS(fiba,&pmesg)) {
            WriteMesg_AS(fcgb,pmesg);
          }
        }
        if(NULL != fiba) { ierr = fclose(fiba); assert(ierr == 0); fiba = NULL;}
      }
#endif        
      if(NULL != fcgb) { ierr = fclose(fcgb); assert(ierr == 0); fcgb = NULL;}
    }
    
    output_the_chunks
      (/*Input Only*/
       WriteMesg_AS,
       heapva->frags,
       heapva->edges,
       heapva->frag_annotations,
       heapva->chunkfrags,
       heapva->thechunks,
       heapva->chunkseqs,
       heapva->chunkquas,
       heapva->chunksrcs,
       analysis_level,
       output_fom_messages,
       gstate->global_fragment_arrival_rate,
       fragment_count_target,
       Graph_Store_File_Prefix
       // heapva->the_imp_messages,
       // heapva->the_ium_messages
       );


    if(TIMINGS) {
      tp2 = time(NULL); 
      fprintf(stderr,
              "%10" F_TIME_TP " sec: Finished chunk graph I/O\n",(tp2-tp1));
      system_date();
    }
  }
}


static
void do_bubble_smoothing
(
 FragStoreHandle TheFragStore,
 THeapGlobals *heapva,
 float global_arrival_rate,
 char *file_prefix,
 char ** bubble_overlaps_filename,
 char ** bubble_boundaries_filename
 )
{
  /* BUBBLE SMOOTHING CONFIGURATION PARAMETERS - See also
     AS_CGB_Bubble.h           - Flag to enable debugging output, 
     AS_CGB_Bubble_Popper.h    - Alignment parameters. 
     AS_CGB_Bubble_VertexSet.h - Default parameters for bubble finding. */

  const char * const ovl_file_suffix = ".bubble_edges.ovl";
  /* The OVL records needed to remove the
     bubbles is written to a filename 
     generated by adding the prefix parameter
     to this suffix. */
  
  int enable_aux_file_output = 1;
  // Set to 1 to produce the two auxillary files below.  0 disables.
  // File to which a textual list of the bubbles is written.

  /* Other local variables: */
  time_t start_time, end_time;
  AS_CGB_Bubble_List_t bubbles = NULL, bubbles_tmp = NULL;

  fprintf(stderr, "*** BEGIN BUBBLE SMOOTHING\n");
  start_time = time(NULL);

  {
    FILE *bubble_overlaps_file = NULL;
    static char bubble_boundaries_filename_tmp[500] = {0};
    static char bubble_overlaps_filename_tmp[500] = {0};
    sprintf(bubble_boundaries_filename_tmp,"%s.cgb_bubbles_txt",file_prefix);
    *bubble_boundaries_filename = bubble_boundaries_filename_tmp;
    sprintf(bubble_overlaps_filename_tmp, "%s%s", file_prefix, ovl_file_suffix);
    *bubble_overlaps_filename = bubble_overlaps_filename_tmp;
    bubble_overlaps_file = fopen(bubble_overlaps_filename_tmp, "w");
    if (!bubble_overlaps_file) {
      fprintf(stderr, "!!! ERROR: Could not open %s for OVL output.\n", 
              bubble_overlaps_filename_tmp);
    } else {
      AS_CGB_Bubble_find_and_remove_bubbles
        ( TheFragStore,
          heapva->frags, heapva->edges,
          heapva->thechunks, heapva->chunkfrags, 
          global_arrival_rate,
          bubble_overlaps_file, stderr,
          file_prefix);
      SAFE_FCLOSE(bubble_overlaps_file);
    }
  }
  end_time = time(NULL);
  fprintf(stderr, "  * Bubble smoothing took " F_TIME_T " m, " F_TIME_T " s (including output)\n",
	  ((end_time - start_time) / 60), ((end_time - start_time) % 60));
  {
    char command_line[1000]={0};
    sprintf(command_line,"rm -f %s",*bubble_boundaries_filename);
    system(command_line);
    // Eliminate spurious files...;)
  }
  if (enable_aux_file_output) {
    FILE * bubble_boundaries_file = NULL;
    fprintf(stderr, "  * Outputting Auxillary Bubble Information\n");
    fprintf(stderr, "  * %s ... \n", *bubble_boundaries_filename);
    SAFE_FOPEN(bubble_boundaries_file, *bubble_boundaries_filename, "w");

    if (!bubble_boundaries_file) {
      fprintf(stderr, "  * ... FAILED.  Could not open %s.\n",
              *bubble_boundaries_filename) ;
    } 
    else {
      /* NOTE: 0's in following call indicate use of defaults. */
      bubbles = AS_CGB_Bubble_find_bubbles(heapva->frags, heapva->edges, 
					   0, 0, 0);
      bubbles_tmp = bubbles;
      while (bubbles_tmp != NULL) {
	fprintf(bubble_boundaries_file, F_IID " %d " F_IID " %d\n", 
		get_iid_fragment(heapva->frags, bubbles_tmp->start),
		bubbles_tmp->start_sx,
		get_iid_fragment(heapva->frags, bubbles_tmp->end),
		bubbles_tmp->end_sx);
	bubbles_tmp = bubbles_tmp->next;
      }
      SAFE_FCLOSE(bubble_boundaries_file);
      fprintf(stderr, "  * ... Done.\n");
    }
  }
    
  fprintf(stderr, "*** END BUBBLE SMOOTHING\n");
}


/****************************************************************************/
/****************************************************************************/

static void FinalizeGlobals( UnitiggerGlobals * rg)
{
#ifdef USE_STATIC_STRINGS_2
  if(NULL != rg->program_name) { free(rg->program_name);}
  if(NULL != rg->OVL_Store_Path) { free(rg->OVL_Store_Path);}
  if(NULL != rg->Input_Graph_Store) { free(rg->Input_Graph_Store);}
  if(NULL != rg->Output_Graph_Store) { free(rg->Output_Graph_Store);}
  if(NULL != rg->Dump_File_Name) { free(rg->Dump_File_Name);}
  if(NULL != rg->ovl_files_list_fname) { free(rg->ovl_files_list_fname);}
  if(NULL != rg->Fragment_Subset_IIDs_File_Name) { free(rg->Fragment_Subset_IIDs_File_Name);}
#endif // USE_STATIC_STRINGS_2
}

// set default or uninitialized values for important variables
static void InitializeGlobals
(
  UnitiggerGlobals * rg,
  char * program_name
  )
{
  memset( rg, 0, sizeof( UnitiggerGlobals ) );
#ifndef USE_STATIC_CMD_STRINGS
  rg->program_name = program_name;
#else // USE_STATIC_CMD_TRINGS
  strcpy(rg->program_name,program_name);
#endif //USE_STATIC_CMD_STRINGS
  rg->input_store_flag = FALSE;
  rg->output_store_flag = FALSE;
  rg->create_store_flag = FALSE;
  rg->append_store_flag = FALSE;
  rg->clobber_store_flag = FALSE;
  rg->developer_mode_flag = FALSE;
  rg->num_threads = 4;
  rg->work_limit_per_candidate_edge = 1000;
  rg->remove_blizzard_overlaps = FALSE;
  rg->use_all_overlaps_in_reaper_pass = TRUE;
  rg->as_cgb_max_frag_iid = 100000 * CGB_MULTIPLIER;
  rg->dvt_double_sided_threshold_fragment_end_degree = 1;
  rg->con_double_sided_threshold_fragment_end_degree = 1;
  rg->intrude_with_non_blessed_overlaps_flag = FALSE;
  rg->cutoff_fragment_end_degree = 1000000;
#ifdef REPAIR_BREAKERS
  rg->breaker_fix = 3;
#endif // REPAIR_BREAKERS
  rg->dechord_the_graph = TRUE;
  rg->check_point_level = 0;
  rg->compress_the_graph = TRUE;
  rg->cgb_unique_cutoff = CGB_UNIQUE_CUTOFF;
  rg->walk_depth = 100;
  rg->iv_start = 0;
#if 0
  rg->overlap_error_threshold = AS_READ_ERROR_RATE; // cds/AS/inc/AS_global.h
#else
  rg->overlap_error_threshold = PerMil_to_CGB_ERATE_TYPE(1000);
#endif
  rg->work_limit_placing_contained_fragments = 20;
  rg->walk_depth=100;
  rg->iv_start=0;
  rg->num_threads = 4;
  rg->output_iterations_flag = TRUE;
  rg->aggressive_spur_fragment_marking = TRUE;
}

static int ParseCommandLine
(
  UnitiggerGlobals * rg,
  int argc,
  char * argv[]
  )
{
  int ierr = 0;
  int illegal=FALSE;
  int blessed_overlaps_output_flag = FALSE;
  
  { /* Parse the argument list using "man 3 getopt". */ 
    int ch,errflg=0;
    optarg = NULL;
    while (!errflg && 
	   ((ch = getopt(argc, argv, 
			 "A:B:C:D:E:F:G:H:I:J:K:L:M:N:O:PQ:R:S:V:U:W:XY:Z:"
			 "ab:cd:e:fg:hi:j:l:m:n:o:p:qr:st:u:v:w:x:y:z:"
                         "156:789"
                         )) != EOF)) {
      // Used command line options:
      // Supp. Assembler:               P          a c  f hi     op
      // Unsup. Assem.  :  CDE            R     X                     t v
      // FGB CGB common :A                  
      // FGB specific   :      G IJ L N  Q    VW      de       mn        wxyz
      // CGB specific   : B        K         U   YZ b       j l    q s
#ifdef REPAIR_BREAKERS
      // repair_cgb_spur:     F H      O   S             g          r
#endif // REPAIR_BREAKERS
      // left-overs     :            M      T                k             
      switch(ch) {
        /* The required command line options: */
      case 'A':
        // -A <int> : Specify an analysis level.
        rg->analysis_level = atoi(optarg);
        //analysis_flag = (analysis_level > 0);
        rg->analysis_flag = rg->analysis_level;
        break;
      case 'B':
        // -B <int> : Specifies a target number of fragments in IMP
        // records per cgb file.
        rg->fragment_count_target = atoi(optarg);
        break;
      case 'C':
#if 0
        // -C <string> : The check pointing instructions.
        char check_point_instructions[CMD_BUFFER_SIZE-1] = "";
#ifndef USE_STATIC_CMD_STRINGS
        check_point_instructions = optarg;
#else // USE_STATIC_CMD_STRINGS
        strcpy(check_point_instructions,optarg);
#endif // USE_STATIC_CMD_STRINGS
        fprintf(stderr,"* Check-pointing instructions <%s>\n",
                check_point_instructions);
#else
        // -C <int> : The check pointing level.
        rg->check_point_level = atoi(optarg);
        assert(rg->check_point_level > 0);
        fprintf(stderr,"* check_point_level <%d>\n",
                rg->check_point_level);
#endif        
        break;
      case 'D':
        // -D <int> : Specify a debugging level.
        rg->debug_level = atoi(optarg);
        break;
#if 0
      case 'E':
        // -E <int> : Set the allowed number of input errors.
        rg->maxerrs = atoi(optarg);
        fprintf(stderr,"* maxerrs set to %d\n", rg->maxerrs );
        break;
#endif        
      case 'F':
        // -F <path> : identify fragment store
        rg->frag_store = optarg;
        fprintf( stderr, " * frag store is %s\n", rg->frag_store );
        break;
      case 'G':
        // -G <file> : 
#ifndef USE_STATIC_CMD_STRINGS
        rg->Fragment_Subset_IIDs_File_Name = optarg;
#else // USE_STATIC_CMD_STRINGS
        strcpy(rg->Fragment_Subset_IIDs_File_Name,optarg);
#endif // USE_STATIC_CMD_STRINGS
        fprintf(stderr,"* Fragment subset IIDs file set to <%s>.\n",
                rg->Fragment_Subset_IIDs_File_Name);
        break;
      case 'H':
        // -H <filename> : identify chimeras file
        rg->chimeras_file = optarg;
        fprintf( stderr, " * chimeras file is %s\n", rg->chimeras_file );
        break;
      case 'I':
        // -I <file> : The input overlapper overlap store.
#ifndef USE_STATIC_CMD_STRINGS
        rg->OVL_Store_Path = optarg;
#else // USE_STATIC_CMD_STRINGS
        strcpy(rg->OVL_Store_Path,optarg);
#endif // USE_STATIC_CMD_STRINGS
        fprintf(stderr,"* The input overlapper overlap store is <%s>.\n",
                rg->OVL_Store_Path);
        break;
      case 'J':
        // -J <boolean> : Remove the blizzard overlaps.
        rg->remove_blizzard_overlaps = atoi(optarg);
        break;
      case 'K':
        // -K Input a file of known fragment-end based
        // branch-points.
        fprintf(stderr,"* -N %s\n", optarg);
#ifndef USE_STATIC_CMD_STRINGS
        rg->branch_points_input_file_name = optarg;
#else // USE_STATIC_CMD_STRINGS
        strcpy(rg->branch_points_input_file_name, optarg);
#endif // USE_STATIC_CMD_STRINGS
        break;
      case 'L':
#ifndef USE_STATIC_CMD_STRINGS
        rg->ovl_files_list_fname = optarg;
#else // USE_STATIC_CMD_STRINGS
        strcpy(rg->ovl_files_list_fname,optarg);
#endif // USE_STATIC_CMD_STRINGS
        fprintf(stderr,"*  ovl_files_list_fname <%s>\n",
                rg->ovl_files_list_fname);
        break;
      case 'M':
        // -M <int> : //un-used
        break;
      case 'N':
        // -N <int> : max IID of the fragments.
        rg->as_cgb_max_frag_iid = atol(optarg);
	fprintf(stderr,"** as_cgb_max_frag_iid = " F_IID "\n",
                rg->as_cgb_max_frag_iid);
        break;
      case 'O':
        // -O <filename> : identify overlap file
        rg->ovl_file = optarg;
        fprintf( stderr, " * overlap file is %s\n", rg->ovl_file );
        break;
      case 'P':
        // -P : Any "protoIO" messages will be in ASCII rather than
        // binary.
	fprintf(stderr,"** ASCII mode\n");
        rg->as_proto_output = TRUE;
        break;
      case 'Q':
        // -Q <int>: Specify Reaper pass.
        rg->reaper_pass = atoi(optarg);
        rg->use_all_overlaps_in_reaper_pass = (rg->reaper_pass == 0);
        break;
#if 0
      case 'R':
        // -R <string> : The restart instructions.
        char restart_instructions[CMD_BUFFER_SIZE-1] = "";
#ifndef USE_STATIC_CMD_STRINGS
        restart_instructions = optarg;
#else // USE_STATIC_CMD_STRINGS
        strcpy(restart_instructions,optarg);
#endif // USE_STATIC_CMD_STRINGS
        fprintf(stderr,"* Restart instructions <%s>\n",
                restart_instructions);
        break;
#else
      case 'R':
        // -R int : The particular transitive edge graph reduction
        // scheme to use.
        fprintf(stderr,"* -R is obsolete\n");
        break;
#endif        
      case 'S':
        // -S <filename> : identify spurs file
        rg->spurs_file = optarg;
        fprintf( stderr, " * spurs file is %s\n", rg->spurs_file );
        break;
      case 'T':
        // -T :
        fprintf(stderr,"* -T is no longer used\n");
        break;
      case 'V':
        // -V int : The transitive edge reduction will restart at this
        // fragment^s adjacency list.
        rg->iv_start = atoi(optarg);
        fprintf(stderr,"* iv_start set to " F_IID "\n",
                rg->iv_start);
        break;
      case 'U':
        // -U 
	rg->bubble_smoothing_flag = atoi(optarg);
	fprintf(stderr, "* Bubble smoothing is now ");
	if (rg->bubble_smoothing_flag)
	  fprintf(stderr, "ON.\n");
	else
	  fprintf(stderr, "OFF.\n");
	break;
      case 'W':
        // -W int : The maximum path walk depth to allow during
        // transitive edge graph reduction.
        rg->walk_depth = atoi(optarg);
        fprintf(stderr,"* walk_depth set to %d\n",
                rg->walk_depth);
        break;
      case 'X':
        // -X : Turn on developer mode.
        rg->developer_mode_flag = TRUE;
        fprintf(stderr,
                "* -X: Turn on developer mode.\n");
        break;
      case 'Y':
        // -Y 
	rg->dont_count_chimeras = atoi(optarg);
	break;
      case 'Z':
        // -Z <boolean> : Remove the blizzard overlaps
        rg->remove_blizzard_overlaps = atoi( optarg );
        break;


      case 'a':
        // -a : Append to the input store.
        rg->append_store_flag = TRUE;
        fprintf(stderr,"* Appending the input store.\n");
        break;
      case 'b':
        // -b <int>
	rg->num_cgb_passes = atoi(optarg);
	assert(rg->num_cgb_passes >= 0);
	break;
      case 'c':
        // -c : Create a new store.
        rg->create_store_flag = TRUE;
        fprintf(stderr,
                "* -c: There is no input store, creating one.\n");
        fprintf(stderr,"* Creating a new fragment graph store.\n");
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
          const int float_erate_value = (NULL != strstr(optarg,"."));
          if(float_erate_value) {
            rg->overlap_error_threshold = cds_float32_to_CGB_ERATE_TYPE(atof(optarg));
          } else {
            rg->overlap_error_threshold = PerMil_to_CGB_ERATE_TYPE(atoi(optarg));
          }
          fprintf(stderr,"The overlap_error_threshold = %d per thousand\n",
                  rg->overlap_error_threshold);
        }
        break;
      case 'f':
        // -f : Over-write the output store if it already exists.
        rg->clobber_store_flag = TRUE;
        fprintf(stderr,
                "* -f: Over-write the output store if it already exists.\n");
        break;
#ifdef REPAIR_BREAKERS
      case 'g':
        // -g <int> : specify modifications to be made
        rg->breaker_fix = atoi( optarg );
        switch( rg->breaker_fix )
          {
          case 1:
            fprintf( stderr,
                     " * modifications: marking all\n" );
            break;
          case 2:
            fprintf( stderr,
                     " * modifications: removing all\n" );
            break;
          case 3:
            fprintf( stderr,
                     " * modifications: adding overlaps & marking\n" );
            break;
          case 4:
            fprintf( stderr,
                     " * modifications: adding overlaps & removing\n" );
            break;
          default:
            fprintf(stderr, "ERROR: Invalid -g repair mode\n" );
            goto UsageStatement;
            break;
          }
        break;
#endif // REPAIR_BREAKERS
      case 'h':
        // help
        goto UsageStatement;
        // Usage( argv[0], NULL );
        break;
      case 'i':
        // -i <file> : The input fragment overlap graph store.
#ifndef USE_STATIC_CMD_STRINGS
        rg->Input_Graph_Store = optarg;
#else // USE_STATIC_CMD_STRINGS
        strcpy(rg->Input_Graph_Store,optarg);
#endif // USE_STATIC_CMD_STRINGS
        fprintf(stderr,"* The input fragment graph store is <%s>.\n",
                rg->Input_Graph_Store);
	assert(NULL != rg->Input_Graph_Store);
        rg->input_store_flag = (rg->Input_Graph_Store[0] != '\0');
        break;
      case 'j':
        // -j <float> : Astat cut-off threshold
	rg->cgb_unique_cutoff = atof(optarg);
	fprintf(stderr,"* cgb_unique_cutoff set to %f\n", rg->cgb_unique_cutoff);
	break;
      case 'l':
        // -l <int> : the length of the genome
	rg->nbase_in_genome = STR_TO_INT64(optarg, NULL, 10);
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
        // -o <file> : The output fragment overlap graph store.
#ifndef USE_STATIC_CMD_STRINGS
        rg->Output_Graph_Store = optarg;
#else // USE_STATIC_CMD_STRINGS
        strcpy(rg->Output_Graph_Store,optarg);
#endif // USE_STATIC_CMD_STRINGS
        fprintf(stderr,"* The output fragment graph store is <%s>.\n",
                rg->Output_Graph_Store);
	assert(NULL != rg->Output_Graph_Store);
        rg->output_store_flag = (rg->Output_Graph_Store[0] != '\0');
        break;
#if 0
      case 'p':
        // -p <file> : The parameters file.
        {
          char The_Parameters_File_Name[CMD_BUFFER_SIZE-1] = "";
          // The file for default parameters.
#ifndef USE_STATIC_CMD_STRINGS
          The_Parameters_File_Name = optarg;
#else // USE_STATIC_CMD_STRINGS
          strcpy(The_Parameters_File_Name,optarg);
#endif // USE_STATIC_CMD_STRINGS
          fprintf(stderr,"* The parameters file is <%s>.\n",
                  The_Parameters_File_Name);
          loadParams(The_Parameters_File_Name);
        }
        break;
#endif


#if 0
      case 'p':
        // -p <file> : The parameters file.
        {
          char *allParams = NULL;
          char Parameter_File_Name[CMD_BUFFER_SIZE-1] = "";
          // The file for default parameters.
          int iret = 0;
          
          fprintf(stderr,"* -p %s\n", optarg);
          strcpy(Parameter_File_Name, optarg);
#if 0
          rg->params_file = optarg;
          fprintf( stderr, "* Parameters file is <%s>.\n", rg->params_file );
          loadParams( rg->params_file );
#endif          
          /* Load parameters file */
          iret = loadParams(Parameter_File_Name);
          if(iret != PARAM_PROC_SUCCESS) {
            fprintf(stderr,
                    "Parameter file <%s> not successfully processed\n",
                    Parameter_File_Name);
          }
          getAllParams("cgb", &allParams);
          if(allParams != NULL) {
            // char *GlobalParamText = NULL;
            // Now a global variable to be placed into a ADL record.
            char *param_string = NULL;
            if(GlobalParamText != NULL) free(GlobalParamText);
            GlobalParamText = (char *)safe_malloc(strlen(allParams)+100);
            sprintf(GlobalParamText,
                    "Parameters from parameter file are:\n%s\n",
                    allParams);
            fprintf(stderr,"* Loaded parameters %s\n", allParams);
            
            /* Load all command line parameters from parameters file */
            param_string = getParam("cgb.estimated_genome_length");
            if(param_string != NULL){
              fprintf(stderr,
                      "* Read estimated genome length: %s "
                      "from parameter file\n",
                      param_string);
              nbase_in_genome = STR_TO_INT64(param_string, NULL, 10);
            }
            free(param_string);
          }
          free(allParams);
        }
        break;
#endif        
      case 'q':
	rg->dont_find_branch_points = TRUE;
        break;
      case 'r':
        // -r <filename> : identify iid file
        rg->iid_file = optarg;
        fprintf( stderr, " * iid file is %s\n", rg->iid_file );
        break;
      case 's':
        // -s 
        rg->aggressive_spur_fragment_marking = FALSE;
        // Aggressive removal of the early spurs.
        break;
      case 't':
        // -t <int> : Specify the number of threads.
        rg->num_threads = atoi(optarg);
        break;

      case 'u':
        // -u <file> : Create a OVL file compatible dump of the fragment
        // graph store.
#ifndef USE_STATIC_CMD_STRINGS
        rg->Dump_File_Name = optarg;
#else // USE_STATIC_CMD_STRINGS
        strcpy(rg->Dump_File_Name,optarg);
#endif // USE_STATIC_CMD_STRINGS
        fprintf(stderr,"* OVL dump file set to <%s>.\n", rg->Dump_File_Name);
        rg->create_dump_file = (rg->Dump_File_Name[0] != '\0');
        break;
      case 'v':
        // -v <int> : Specify a verbosity level.
      rg->verbosity_level = atoi(optarg);
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
        //rg->blessed_overlaps_output_filename = optarg;
        break;
      case '8':
        // -8 reserved for debugging.
        skip_fgb = TRUE;
        break;
      case '9':
        // -9 reserved for debugging.
        skip_cgb = TRUE;
        break;
      case '?':
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

#ifdef REPAIR_BREAKERS
    // get the list of cgb files for repair_breakers.
    rg->num_cgb_files = argc - optind;
    if( rg->num_cgb_files <= 0 )
      Usage( argv[0], "ERROR: Need one or more .cgb files\n" );
    rg->the_cgb_files = &(argv[optind]);
#endif // REPAIR_BREAKERS

    assert(rg->as_cgb_max_frag_iid < AS_CGB_NOT_SEEN_YET);
    if(rg->as_cgb_max_frag_iid == 0) {
      fprintf(stderr,
              "UNITIGGER ERROR: Until I get around to implementing a hash table"
              " for the encountered IIDs, the \"unitigger -N <max_iid>\""
              " command line option is necessary.\n");
      exit(1);
    }
    // The -c option means that the input store string is empty.
    if(rg->input_store_flag && rg->create_store_flag ) {
      fprintf
        (stderr,"ERROR: The -c and -i options incompatible.\n");
      illegal = TRUE;
      goto UsageStatement;
    }
    if((! rg->input_store_flag) && (! rg->create_store_flag)) {
      fprintf
        (stderr,"ERROR: The -c or -i command line option required.\n");
      illegal = TRUE;
    }
    
    // The -a option means that the output store string is empty.
    if(rg->output_store_flag && rg->append_store_flag && rg->create_dump_file) {
      fprintf
        (stderr,"ERROR: The -a and -o and -O options incompatible.\n");
      illegal = TRUE;
      goto UsageStatement;
    }
    if((! rg->output_store_flag) && (! rg->append_store_flag)
       && (! rg->create_dump_file) ) {
      fprintf
        (stderr,"ERROR: One of the -a or -o or -O command line options required.\n");
      illegal = TRUE;
    }
    if(rg->output_store_flag + rg->append_store_flag + rg->create_dump_file != 1) {
      fprintf
        (stderr,"ERROR: Only one of -a or -o or -O command line options can be specified.\n");
      fprintf(stderr,"%d %d %d\n",
              rg->output_store_flag, rg->append_store_flag, rg->create_dump_file);
      illegal = TRUE;
      goto UsageStatement;
    }

    assert( !illegal || (rg->create_store_flag ^ rg->input_store_flag));
    assert( !illegal || (rg->append_store_flag ^ rg->output_store_flag));
    
    if(rg->create_store_flag && rg->append_store_flag) {
#if 0
      fprintf
        (stderr,"ERROR: Cannot append to a null input store.\n");
      illegal = TRUE;
      goto UsageStatement;
#else
      fprintf
        (stderr,"WARNING: Appending to a null input store.\n");
#endif
    }

    if(rg->input_store_flag) {
      char buffer[CMD_BUFFER_SIZE-1];
#if 0
      if(rg->append_store_flag) {
        sprintf(buffer,"chmod -R u+w %s", rg->Input_Graph_Store);
        ierr = system(buffer);
        if(ierr) {
          fprintf(stderr,
                  "ERROR: Trouble getting write access to the input graph "
                  "store directory.\n");
        illegal = TRUE;
        goto UsageStatement;
        }
      }
#endif      
      sprintf(buffer,"(cd %s)" , rg->Input_Graph_Store);
      ierr = system(buffer);
      if(ierr) {
        fprintf(stderr,
                "ERROR: Trouble validating "
                "the input graph store directory.\n");
        illegal = TRUE;
        goto UsageStatement;
      }
    }
      
    if(rg->append_store_flag) {
#ifndef USE_STATIC_CMD_STRINGS
      rg->Output_Graph_Store = rg->Input_Graph_Store;
#else // USE_STATIC_CMD_STRINGS
      strcpy(rg->Output_Graph_Store,rg->Input_Graph_Store);
#endif // USE_STATIC_CMD_STRINGS
    }

    if(rg->output_store_flag) {
      char buffer[CMD_BUFFER_SIZE-1];
      if(rg->clobber_store_flag) {
        sprintf(buffer,"rm -rf %s", rg->Output_Graph_Store);
        ierr = system(buffer);
        if(ierr) {
          fprintf(stderr,
                  "ERROR: Trouble clobbering the output graph "
                  "store directory.\n");
          illegal = TRUE;
          goto UsageStatement;
        }
      }
      sprintf(buffer,"mkdir %s" , rg->Output_Graph_Store);
      ierr = system(buffer);
      if(ierr) {
        fprintf(stderr,
                "ERROR: Trouble creating the output graph store directory.\n"
                "       Did the directory already exist?\n");
        illegal = TRUE;
        goto UsageStatement;
      }
    }
    
    if(rg->bubble_smoothing_flag && (NULL == rg->frag_store)) {
      fprintf(stderr,"Error: bubble smoothing needs a fragment store.\n");
      assert(NULL != rg->frag_store);
    }

#if 0
    if(argc - optind == 0) {
      fprintf
        (stderr,"ERROR: The -i command line option was not used "
         "and no input was supplied.\n");
      exit(1);
    }
#endif

    {
      char *Graph_Store_File_Name   = NULL;
      char *Graph_Store_File_Suffix = NULL;
      
      if(NULL != rg->Output_Graph_Store) {
        char *theGraphStorePath = rg->Output_Graph_Store;
        static char Graph_Store_File_Name_Original[CMD_BUFFER_SIZE];
        static char Graph_Store_File_Name_Tmp[CMD_BUFFER_SIZE];
        fprintf(stderr, __FILE__ " theGraphStorePath = <%s>\n", theGraphStorePath);
        strcpy(Graph_Store_File_Name_Original,theGraphStorePath);
        Graph_Store_File_Name = Graph_Store_File_Name_Original;
        fprintf(stderr,"Graph_Store_File_Name = %s\n",Graph_Store_File_Name);
        strcpy(Graph_Store_File_Name_Tmp,Graph_Store_File_Name);
        rg->Output_Graph_Store_Prefix = Graph_Store_File_Name_Tmp;
        Graph_Store_File_Suffix = strrchr(Graph_Store_File_Name_Tmp,'.');
        if( Graph_Store_File_Suffix != NULL) {
          *Graph_Store_File_Suffix = '\0'; Graph_Store_File_Suffix++;
        }
      }
      
      fprintf(stderr, __FILE__ " %d: rg->Output_Graph_Store   = %s\n",
              __LINE__, rg->Output_Graph_Store);
      fprintf(stderr, __FILE__ " %d: rg->Output_Graph_Store_Prefix = %s\n",
              __LINE__, rg->Output_Graph_Store_Prefix);
    }
    
    if(! rg->dont_count_chimeras) {
      static char chimeras_report_filename[128]={0};
      static char crappies_report_filename[128]={0};
      assert(NULL != rg->Output_Graph_Store_Prefix);

      sprintf(crappies_report_filename,"%s.cgb_crappies",
	      rg->Output_Graph_Store_Prefix);
      rg->spurs_file = crappies_report_filename;

      sprintf(chimeras_report_filename,"%s.cgb_chimeras",
	      rg->Output_Graph_Store_Prefix);
      rg->chimeras_file = chimeras_report_filename;
    }

    if(blessed_overlaps_output_flag) {
      static char blessed_overlaps_output_filename[128]={0};
      assert(NULL != rg->Output_Graph_Store_Prefix);
      sprintf(blessed_overlaps_output_filename,"%s.blovl",
	      rg->Output_Graph_Store_Prefix);
      rg->blessed_overlaps_output_filename = blessed_overlaps_output_filename;
    }

  UsageStatement:
    if((illegal == TRUE)) {
	fprintf
          (stderr, "USAGE: %s <option>*\n"
           "\tThe available options are:\n"
           "\t[-A <int>] Run the fragment overlap graph analyzer.\n"
	   "\t[-B <int>] Specifies the target number of fragments per partition.\n"
           "\t[-C <int>] Check point level for restart. (not implemented)\n"
           "\t[-D <int>] Specify a debugging level. (not implemented)\n"
           "\t[-E <int>] Specify the maximum number of soft errors.\n"
	   "\t[-F <directory>] The fragment store name.\n"
	   "\t[-G <file>] A file of fragment IIDs.\n"
	   "\t[-H <filename>] chimeras file.\n"
           "\t[-I <directory>] Read the OVL store.\n"
	   "\t[-J <boolean>] Remove blizzard overlaps in FGB.\n"
	   "\t[-K <filename>] File of known fragment-end based branch-points.\n"
           "\t[-L <filename>] A file that lists the input overlap files.\n"
           "\t[-M <boolean>] Unused\n"
           "\t[-N <max IID of the fragments> ]\n"
#ifdef REPAIR_BREAKERS
	   "\t[-O <filename>] Repair breakers output file.\n"
#endif // REPAIR_BREAKERS
           "\t[-P] Specify ASCII output for protoIO messages.\n"
	   "\t[-Q <int>] Specify Reaper Pass.\n"
           "\t[-R <int>] Which containment rule to use. \n"
	   "\t[-S <filename>] Spurs file.\n"
           "\t[-T] Output transitively removeable overlaps.\n"
	   "\t[-U <boolean>] Find bubble smoothing overlaps.\n"
	   "\t[-V] Transitive edge reduction restarts at this VID.\n"
           "\t[-W <int>] Limit in path length for graph walking.\n"
           "\t[-X] Enable the expert options.\n"
	   "\t[-Y <boolean>] Do not count chimera fragments.\n"
	   "\t[-Z <boolean>] Remove the blizzard overlaps in CGB.\n"
           "\t(-a|-o <GraphStore>) Append or specify the output GraphStore name.\n"
	   "\t[-b <int>] Number of cgb passes for finding branch points.\n"
           "\t(-c|-i <GraphStore>) Create or specify the input GraphStore name.\n"
           "\t[-d <int>] Enable/Disable de-chording of the fragment overlap graph.\n"
           "\t[-e <float or int>] Overlaps with error rate about this are ignored on input.\n"
           "\t\t An integer value is in parts per thousand.\n"
           "\t[-f] Allow the output graph store to be over-written.\n"
#ifdef REPAIR_BREAKERS
	   "\t[-g <int>] Spur and chimera repair approach.\n"
#endif // REPAIR_BREAKERS
	   "\t[-h] Help.\n"
           "\t(-c|-i <GraphStore>) Create or specify the input GraphStore name.\n"
	   "\t[-j <int>] Unique unitig cut-off\n"
	   "\t[-k] Unused\n"
	   "\t[-l <int>] Specify length of the genome.\n"
           "\t[-m <additional number of edge records>] Pre-allocate memory\n"
           "\t[-n <additional number of fragments>] Pre-allocate memory\n"
           "\t(-a|-o <GraphStore>) Append or specify the output GraphStore name.\n"
           "\t[-p <file>] Specify the parameters file.\n"
	   "\t[-q] Do not find branch points.\n"
#ifdef REPAIR_BREAKERS
	   "\t[-r <filename>] Affected fragment IID file for repair breakers.\n"
#endif // REPAIR_BREAKERS
	   "\t[-s] Disable early spur fragment removal.\n"
           "\t[-t <int>] Specify the number of threads.\n"
           "\t[-u <filename>] Create a OVL compatible dump of the graph.\n"
           "\t[-v <int>] Specify a verbosity level.\n"
           "\t[-w <int>] The work limit per candidate edge for de-chording.\n"
           "\t[-x <int>] Dovetail outgoing degree threshold per fragment-end.\n"
           "\t[-y <int>] Ignore non-blessed overlap edges to blessed fragment ends.\n"
           "\t[-z <int>] Containment outgoing degree threshold per fragment-end.\n"
           ,
           argv[0]);
        exit (EXIT_FAILURE);
      }

    /* The remaining command line arguments are input Assembler
       message files. */
  }
  return 0;
}

static void InitializeUnitiggerGlobals
(
 int argc,
 char * argv [],
 //TStateGlobals    * gstate,
 //THeapGlobals     * heapva,
 UnitiggerGlobals * rg
 )
{
  InitializeGlobals( rg, argv[0] );
  ParseCommandLine( rg, argc, argv);
}


// static size_t
// ReportMemorySize_VA( VarArrayType *va, char *name, FILE *stream)

static void BasicUnitigger
(
 int                argc,
 char             * argv[],
 TStateGlobals    * gstate,
 THeapGlobals     * heapva,
 UnitiggerGlobals * rg
 )
{
  int status = 0;
  
  ReportHeapUsage_CGB( gstate, heapva, stderr); 
  
  /* Create globals */ {
    // read the fgb store to initialize the heap.
    time_t tp1 = time(NULL);
    time_t tp2;
    
    fprintf(stderr,"Begin reading the fragment graph store.\n");
    system_date();
    read_fgb_store(
                    rg->Input_Graph_Store,
                    gstate, heapva,
                    rg->maxfrags, 
                    rg->maxedges, 
                    rg->maxtext);
    tp2 = time(NULL); 
    fprintf(stderr,"%10" F_TIME_TP " sec: Finished reading the fragment graph store.\n",
            (tp2-tp1));
    system_date();
  }

  if( ! skip_fgb ) {
  ReportHeapUsage_CGB( gstate, heapva, stderr); 
  status = main_fgb( argc, argv, gstate, heapva, rg);
  fprintf(stderr,"Basic Unitigger main_fgb status = %d\n",
          status);
  skip_fgb = FALSE;
  }
  
  if( ! skip_cgb ) {
  ReportHeapUsage_CGB( gstate, heapva, stderr); 
  status = main_cgb( argc, argv, gstate, heapva, rg);
  fprintf(stderr,"Basic Unitigger main_cgb status = %d\n",
          status);
  skip_cgb = FALSE;
  }
  
  ReportHeapUsage_CGB( gstate, heapva, stderr);

}

static void StandardUnitigger
(
 int                argc,
 char             * argv[],
 TStateGlobals    * gstate,
 THeapGlobals     * heapva,
 UnitiggerGlobals * rg
 )
{
  open_fgb_store( gstate, heapva );

  BasicUnitigger( argc, argv, gstate, heapva, rg);

#ifdef REPAIR_BREAKERS
  if(FALSE) {
    int status = 0;
    status = main_repair_globals( argc, argv, gstate, heapva, rg);
    fprintf(stderr,"main_repair_globals status = %d\n", status);
    // Actually should be before bubble smoothing
    exit(1);
    
    //clear_fgb_store( gstate, heapva, TRUE, TRUE, TRUE, TRUE); 
    close_fgb_store( gstate, heapva );
    open_fgb_store( gstate, heapva );
    
    BasicUnitigger( argc, argv, gstate, heapva, rg);
  } // main_repair_globals
#endif // REPAIR_BREAKERS

  if(rg->bubble_smoothing_flag) {
    FragStoreHandle TheFragStore = NULLFRAGSTOREHANDLE;
    const int dont_count_chimeras = rg->dont_count_chimeras;   // Save flag.
    rg->dont_count_chimeras = TRUE;

    if( NULL != rg->frag_store ) { 
      char *theFragStorePath = rg->frag_store;
      fprintf(stderr, __FILE__ " theFragStorePath = <%s>\n", theFragStorePath);
#ifndef LOAD_FRAGSTORE_INTO_MEMORY
      // Open the fragstore, readonly
      TheFragStore = openFragStore(theFragStorePath,"r");
#else // LOAD_FRAGSTORE_INTO_MEMORY
      TheFragStore = loadFragStore(theFragStorePath);
#endif // LOAD_FRAGSTORE_INTO_MEMORY
      assert(TheFragStore != NULLFRAGSTOREHANDLE);
    }
    
#if 0
    {
      const int iterator=99;

      // Checkpoint before Bubble Smoothing
      {
        int ierr;
        char * checkpoint99 = NULL;
        checkpoint99 = (char *)safe_malloc(100+2*sizeof(rg->Output_Graph_Store));
        sprintf(checkpoint99,"cp -r %s %s.99", rg->Output_Graph_Store, rg->Output_Graph_Store);
        // Using "mv" causes the following error during the final output.
        // AS_FGB_main.c Jan 16 2002 09:22:37
        // Creating file to append the new IBA+ADT(+ADL) messages.
        // ERROR: Trouble copying the IBA+ADT file.

        processing_phase_3
          ( gstate,
            heapva,
            rg->Output_Graph_Store,
            rg->analysis_level
            );
        ierr = system(checkpoint99);
        free(checkpoint99);
      }

      // Checkpoint before Bubble Smoothing
      {
        MesgWriter WriteMesg_AS = NULL;
        MesgWriter ErrorWriter_AS = NULL;
        
        // VersionStamp(argc,argv);
        WriteMesg_AS = OutputFileType_AS((rg->as_proto_output
                                          ? AS_PROTO_OUTPUT : AS_BINARY_OUTPUT));
        ErrorWriter_AS = OutputFileType_AS(AS_PROTO_OUTPUT);
        
        if(NULL != rg->Output_Graph_Store_Prefix) {
          save_the_chunk_graph
            (
             WriteMesg_AS,
             rg->Output_Graph_Store_Prefix,
             rg->Output_Graph_Store,
             rg->analysis_level,
             rg->bubble_boundaries_filename,
             rg->output_fom_messages,
             rg->cgb_unique_cutoff,
             iterator, // rg->num_cgb_passes,
             rg->fragment_count_target,
             argc,
             argv,
             gstate,
             heapva
             );
        }
      }
    }
#endif        

    ReportHeapUsage_CGB( gstate, heapva, stderr); 

    // The bubble smoothing code belongs here rather than in main_cgb().
    /* Do bubble smoothing if enabled. */
    if (rg->bubble_smoothing_flag) {
      do_bubble_smoothing
        ( TheFragStore,
          heapva,
          gstate->global_fragment_arrival_rate, 
          rg->Output_Graph_Store_Prefix,
          &(rg->bubble_overlaps_filename),
          &(rg->bubble_boundaries_filename)
          );

      fprintf(stderr,"close the fragStore\n");
      
      if(TheFragStore != NULLFRAGSTOREHANDLE) {
        int ierr;
        ierr = closeFragStore(TheFragStore); assert(ierr == 0);
        TheFragStore = NULLFRAGSTOREHANDLE;
      }
      
      ReportHeapUsage_CGB( gstate, heapva, stderr); 

    }

    // clear_fgb_store operation:
    close_fgb_store( gstate, heapva );
    open_fgb_store( gstate, heapva );

    BasicUnitigger( argc, argv, gstate, heapva, rg);
    
    rg->dont_count_chimeras = dont_count_chimeras; // Restore flag.
  }
  
#if 1
  /////////////////////////////////////////////////////////////////
  fprintf(stderr,"Starting checkpoint before producing Assembler message output.\n");
  ReportHeapUsage_CGB( gstate, heapva, stderr); 

  if(NULL != rg->Output_Graph_Store) {
    fprintf(stderr,"Check point after all the data is in.\n");
    { 
      {	  
        char thePath1[CMD_BUFFER_SIZE-1]={0};
        char thePath2[CMD_BUFFER_SIZE-1]={0};
        int ierr=0;
        sprintf(thePath1,"%s/%s",rg->Output_Graph_Store,"fgb.ckp_tmp");
        sprintf(thePath2,"%s/%s",rg->Output_Graph_Store,"fgb.ckp");
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
        sprintf(thePath4,"%s/%s",rg->Output_Graph_Store,"fga.ckp");
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
  fprintf(stderr,"Finished checkpoint before producing Assembler message output.\n");
  /////////////////////////////////////////////////////////////////
#endif

  {
    MesgWriter WriteMesg_AS = NULL;
    MesgWriter ErrorWriter_AS = NULL;

    // VersionStamp(argc,argv);
    WriteMesg_AS = (MesgWriter)OutputFileType_AS((rg->as_proto_output
                                      ? AS_PROTO_OUTPUT : AS_BINARY_OUTPUT));
    ErrorWriter_AS = (MesgWriter)OutputFileType_AS(AS_PROTO_OUTPUT);
    
    if(NULL != rg->Output_Graph_Store_Prefix) {
      save_the_chunk_graph
        (
         WriteMesg_AS,
         rg->Output_Graph_Store_Prefix,
         rg->Output_Graph_Store,
         rg->analysis_level,
         rg->bubble_boundaries_filename,
         rg->output_fom_messages,
         rg->cgb_unique_cutoff,
         rg->num_cgb_passes,
         rg->fragment_count_target,
         argc,
         argv,
         gstate,
         heapva
         );
    }
  }

  ///////////////////////////////
  // bless some overlaps  

  view_fgb_chkpnt( rg->Output_Graph_Store_Prefix, heapva->frags, heapva->edges);

  {
    {
      // Determine the blessed overlaps.  They are the overlaps
      // that are interior to u-unitigs.  In particular, intrachunk overlaps 
      
      const IntFragment_ID nfrag = GetNumFragments(heapva->frags);
      const IntEdge_ID nedge = GetNumEdges(heapva->edges);
      IntEdge_ID ie;


      for( ie=0; ie < nedge; ie ++) {
        const Tnes nes = get_nes_edge(heapva->edges,ie);
        const IntFragment_ID avx = get_avx_edge(heapva->edges,ie);
        const Tlab alab = get_lab_fragment(heapva->frags,avx);
        const IntChunk_ID chunk_index = get_cid_fragment(heapva->frags,avx);
        // assert(nedge > ie);
        assert(nfrag > avx);
        {
          const AChunkMesg * const mychunk = 
	    GetVA_AChunkMesg(heapva->thechunks,chunk_index);
          const cds_float32 coverage_statistic = mychunk->coverage_stat;
          if(coverage_statistic >= rg->cgb_unique_cutoff) {
              if( AS_CGB_INTRACHUNK_EDGE == nes ) {
                set_blessed_edge(heapva->edges, ie, TRUE);
              }
              if( (AS_CGB_SINGLECONT_FRAG == alab) &&
		  (AS_CGB_CONTAINED_EDGE == nes) && (is_a_frc_edge(heapva->edges,ie))
		  ) {
                set_blessed_edge(heapva->edges, ie, TRUE);
              }
            }
        }
      }
    }

    {
      // Determine the blessed overlaps that are not interior to u-unitigs.
      
      const IntFragment_ID nfrag = GetNumFragments(heapva->frags);
      const IntEdge_ID nedge = GetNumEdges(heapva->edges);
      IntEdge_ID ie;


      for( ie=0; ie < nedge; ie ++) {
        const IntFragment_ID avx = get_avx_edge(heapva->edges,ie);
        const IntFragment_ID asx = get_asx_edge(heapva->edges,ie);
        const IntFragment_ID bvx = get_bvx_edge(heapva->edges,ie);
        const IntFragment_ID bsx = get_bsx_edge(heapva->edges,ie);
        const IntChunk_ID chunk_index = get_cid_fragment(heapva->frags,avx);
        // assert(nedge > ie);
        assert(nfrag > avx);
        if( get_blessed_edge(heapva->edges,ie) ){
          const AChunkMesg * const mychunk = 
	    GetVA_AChunkMesg(heapva->thechunks,chunk_index);
          const cds_float32 coverage_statistic = mychunk->coverage_stat;
          if(coverage_statistic < rg->cgb_unique_cutoff) {
	    fprintf(stdout,"BONIUU: "
		    "afr= " F_IID " asx= %d bfr= " F_IID " bsx= %d olt= %c\n",
		    get_iid_fragment(heapva->frags,avx),asx,
		    get_iid_fragment(heapva->frags,bvx),bsx,
		    (is_a_dvt_edge(heapva->edges,ie) ? 'D' : 'C')
		    );
            }
        }
      }
    }

    view_fgb_chkpnt( rg->Output_Graph_Store_Prefix, heapva->frags, heapva->edges);

    {
      // output latest set of blessed overlap edges.

      // VersionStamp(argc,argv);
      MesgWriter
        WriteMesg_AS = (MesgWriter)OutputFileType_AS
        ((rg->as_proto_output
          ? AS_PROTO_OUTPUT : AS_BINARY_OUTPUT));
      
      FILE *filk = NULL;
      FILE *fcgb = (NULL == rg->blessed_overlaps_output_filename
                    ? NULL
                    : fopen(rg->blessed_overlaps_output_filename,"w")
                    );
      
      if(NULL != fcgb) {
        output_OVL_mesgs
          (/* Input Only*/
           WriteMesg_AS,
           heapva->frags, 
           heapva->edges,
           heapva->frag_annotations,
           /* Read Only */
           filk,
           /* Append Only*/
           fcgb);
        fclose(fcgb);
      }
    }
  }
  
  close_fgb_store( gstate, heapva );
}

int main(int argc, char * argv [])
{
  TStateGlobals gstate = {0};
  THeapGlobals  heapva = {0};
  UnitiggerGlobals rg  = {0};
  int status = EXIT_SUCCESS;
  
  InitializeUnitiggerGlobals( argc, argv, &rg);

  StandardUnitigger( argc, argv, &gstate, &heapva, &rg );

  exit(status);
}






