
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
static char CM_ID[] = "$Id: AS_CGB_unitigger.c,v 1.17 2007-05-01 14:41:43 granger_sutton Exp $";
/*********************************************************************
 *
 * Module: AS_CGB_unitigger.c
 * Description:
 *    Unitigger main
 * 
 *    Programmer:  C. M. Mobarry
 *       Written:  Dec 2001
 *     ReWritten:  Apr 2007, B.Walenz
 * 
 *********************************************************************/

#include "AS_CGB_all.h"
#include "AS_UTL_version.h"
#include "AS_CGB_unitigger_globals.h"
#include "AS_CGB_Bubble.h"

extern int REAPER_VALIDATION;





void
output_the_chunks(const Tfragment frags[],
                  const Tedge     edges[],
                  const VA_TYPE(char) * const fragsrc,
                  const TChunkFrag    chunkfrags[],
                  const TChunkMesg    thechunks[],
                  const VA_TYPE(char) chunkseqs[],
                  const VA_TYPE(char) chunkquas[],
                  const VA_TYPE(char) chunksrc[],
                  const int analysis_level,
                  const float global_fragment_arrival_rate,
                  const int fragment_count_target,
                  const char * const Graph_Store_File_Prefix) {

  IntFragment_ID max_num_frags_per_chunk=0;
  IntEdge_ID max_num_ovlps_per_chunk=0;

  IntChunk_ID       chunk_index;
  const IntChunk_ID nchunks = (IntChunk_ID)GetNumVA_AChunkMesg(thechunks);
  const IntFragment_ID  nfrag   = GetNumFragments(frags);
  
  size_t  * frag_source_index = NULL; 
  VA_TYPE(char)  * newfragsrc = NULL;

  const int nsample=500;
  const int nbucket=500;

  char filename[FILENAME_MAX];
  
  FILE *fcgb = NULL;
  int fragment_count = 0;
  int file_count = 0;

  VA_TYPE(IntMultiPos) * the_imps = CreateVA_IntMultiPos(0);
  
  // Determine the maximum sizes of the variant parts of a chunk message.
  for(chunk_index=0;chunk_index < nchunks; chunk_index++){
    const AChunkMesg * const mychunk = GetVA_AChunkMesg(thechunks,chunk_index);

    max_num_frags_per_chunk = MAX(max_num_frags_per_chunk,
				  mychunk->num_frags);
    max_num_ovlps_per_chunk = MAX(max_num_ovlps_per_chunk,
				  mychunk->a_degree_raw 
				  + mychunk->b_degree_raw);
  }

  if(analysis_level > 1 && (fragsrc != NULL)) {

    frag_source_index = safe_malloc(sizeof(size_t) * nfrag);
    newfragsrc = CreateVA_char(2*GetNumVA_char(fragsrc));

    // Append to the fragment source field.  This needs to be a separate
    // pass since the Assembler I/O routines use pointers AND the
    // variable length arrays are free to realloc the space for the
    // character data.

    for(chunk_index=0;chunk_index < nchunks; chunk_index++){
      const AChunkMesg * const mychunk = GetVA_AChunkMesg(thechunks,chunk_index);
      const IntFragment_ID num_frags = mychunk->num_frags;
      IntFragment_ID ivc;
      // assert(mychunk->f_list >= 0);
      for(ivc=0; ivc<num_frags; ivc++) {
	const IntFragment_ID ivn = mychunk->f_list + ivc; 
	// Get the ivc-th fragment of the chunk.
	// assert(ivn < nfrag);
	// Get the next fragment in the chunk.
	const IntFragment_ID vid     = GetVA_AChunkFrag(chunkfrags,ivn)->vid; 
	const size_t isrc = get_src_fragment(frags,vid);

	size_t offset, nsource;
	char source[1024];
	char clabel;
	int iret = 0;

	switch(get_lab_fragment(frags,vid)) {
	case AS_CGB_SOLO_FRAG:            clabel = 'L'; break;
	case AS_CGB_HANGING_FRAG:         clabel = 'H'; break;
	case AS_CGB_HANGING_CHUNK_FRAG:   clabel = 'B'; break;
	case AS_CGB_HANGING_CRAPPY_FRAG:  clabel = 'F'; break;
	case AS_CGB_THRU_FRAG:            clabel = 'T'; break;
	case AS_CGB_ORPHANEDCONT_FRAG:    clabel = 'O'; break;
	case AS_CGB_MULTICONT_FRAG:       clabel = 'M'; break;
	case AS_CGB_BRANCHMULTICONT_FRAG: clabel = 'P'; break;
	case AS_CGB_INTERCHUNK_FRAG:      clabel = 'E'; break;
	case AS_CGB_INTRACHUNK_FRAG:      clabel = 'I'; break;
	case AS_CGB_SINGLECONT_FRAG:      clabel = 'C'; break;
	case AS_CGB_DELETED_FRAG:         clabel = 'D'; break;
	default:
	  assert(FALSE);
	}
	iret = sprintf(source,"%slab>%c%c\n", 
		       GetVA_char(fragsrc,isrc), clabel,
		       (get_con_fragment(frags,vid) ? 'C' : 'E' ));
	assert(iret > 0);
	nsource = strlen(source);
	offset  = GetNumVA_char(newfragsrc);
	EnableRangeVA_char(newfragsrc,offset+nsource+1);
	/* Remember the terminal null character for a char string. */
	strcpy(GetVA_char(newfragsrc,offset),source);
	frag_source_index[vid] = offset;
      }
    }
  }

  for(chunk_index=0;chunk_index < nchunks; chunk_index++) /* a */ {
    const AChunkMesg * const mychunk = GetVA_AChunkMesg(thechunks,chunk_index);
    const IntFragment_ID num_frags = mychunk->num_frags;
    GenericMesg   pmesg;
    IntUnitigMesg achunk;
    IntFragment_ID ivc;
    int forced  = FALSE;

    const int number_of_randomly_sampled_fragments_in_chunk
      = count_the_randomly_sampled_fragments_in_a_chunk
      ( frags, chunkfrags, thechunks, chunk_index);
    const float coverage_statistic 
      = compute_coverage_statistic
      ( mychunk->rho,
        number_of_randomly_sampled_fragments_in_chunk,
        global_fragment_arrival_rate );

    for(ivc=0; ivc<num_frags; ivc++) {
      const IntFragment_ID ivn = mychunk->f_list + ivc; 
      // Get the ivc-th fragment of the chunk.
      // assert(ivn < nfrag);
      // Get the next fragment in the chunk.
      const IntFragment_ID vid = GetVA_AChunkFrag(chunkfrags,ivn)->vid; 
      // const IntFragment_ID iid  = get_iid_fragment(frags,vid);
      IntMultiPos a_frag;

      a_frag.type   = get_typ_fragment(frags,vid);
      a_frag.ident  = get_iid_fragment(frags,vid);
      a_frag.contained  = get_container_fragment(frags,vid);
      a_frag.position.bgn = get_o5p_fragment(frags,vid);
      a_frag.position.end = get_o3p_fragment(frags,vid);
      a_frag.delta_length = 0;
      a_frag.delta        = NULL;

      SetVA_IntMultiPos(the_imps,ivc,&a_frag);
    }

    achunk.consensus = "";
    achunk.quality   = "";
      
    // The ProtoSpec specifies that the first chunk id is ZERO.
    achunk.iaccession     = mychunk->iaccession;
#ifdef AS_ENABLE_SOURCE
    achunk.source         = GetVA_char(chunksrc,mychunk->isrc);
#endif
    achunk.coverage_stat  = coverage_statistic;
    achunk.status         = AS_UNASSIGNED;
    achunk.a_branch_point = mychunk->a_branch_point;
    achunk.b_branch_point = mychunk->b_branch_point;

    achunk.consensus      = "";
    achunk.quality        = "";
    achunk.length         = mychunk->bp_length;

    achunk.forced         = forced;
    achunk.num_frags      = mychunk->num_frags;
    achunk.f_list         = GetVA_IntMultiPos(the_imps,0);
    achunk.num_vars       = 0;

    fragment_count += mychunk->num_frags;
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

  if( NULL != frag_source_index ) {
    assert(NULL != newfragsrc);
    DeleteVA_char(newfragsrc); newfragsrc = NULL;
    safe_free(frag_source_index);
  }

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

  while (!errflg && 
         ((ch = getopt(argc, argv, 
                       "A:B:C:D:E:F:G:H:I:J:K:L:M:N:O:PQ:R:S:V:U:W:XY:Z:"
                       "ab:cd:e:fg:hi:j:kl:m:n:o:p:qr:st:u:v:w:x:y:z:"
                       "156:789"
                       )) != EOF)) {

    // Used command line options:
    // Supp. Assembler:               P          a c  f hi     op
    // Unsup. Assem.  :  CDE            R     X                     t v
    // FGB CGB common :A                  
    // FGB specific   :      G IJ L N  Q    VW      de       mn        wxyz
    // CGB specific   : B        K         U   YZ b       jkl    q s
    // left-overs     :            M      T                             

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
        // -C <int> : The check pointing level.
        rg->check_point_level = atoi(optarg);
        assert(rg->check_point_level > 0);
        fprintf(stderr,"* check_point_level <%d>\n", rg->check_point_level);
        break;
      case 'D':
        // -D <int> : Specify a debugging level.
        rg->debug_level = atoi(optarg);
        break;
      case 'F':
        // -F <path> : identify fragment store
        rg->frag_store = optarg;
        break;
      case 'G':
        // -G <file> : 
        rg->Fragment_Subset_IIDs_File_Name = optarg;
        fprintf(stderr,"* Fragment subset IIDs file set to <%s>.\n", rg->Fragment_Subset_IIDs_File_Name);
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
      case 'J':
        // -J <boolean> : Remove the blizzard overlaps.
        rg->remove_blizzard_overlaps = atoi(optarg);
        break;
      case 'L':
        rg->ovl_files_list_fname = optarg;
        fprintf(stderr,"*  ovl_files_list_fname <%s>\n",
                rg->ovl_files_list_fname);
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
      case 'Z':
        // -Z <boolean> : Remove the blizzard overlaps
        rg->remove_blizzard_overlaps = atoi( optarg );
        break;

      case 'b':
        // -b <int>
        rg->num_cgb_passes = atoi(optarg);
        assert(rg->num_cgb_passes >= 0);
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
        rg->Output_Graph_Store_Prefix = (char *)optarg;
        break;

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

      case 'u':
        // -u <file> : Create a OVL file compatible dump of the fragment
        // graph store.
        rg->Dump_File_Name = optarg;
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

  assert(rg->as_cgb_max_frag_iid < AS_CGB_NOT_SEEN_YET);
  if(rg->as_cgb_max_frag_iid == 0) {
    fprintf(stderr,
            "UNITIGGER ERROR: Until I get around to implementing a hash table"
            " for the encountered IIDs, the \"unitigger -N <max_iid>\""
            " command line option is necessary.\n");
    exit(1);
  }


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
            "\t-A <int>        Run the fragment overlap graph analyzer.\n"
            "\t-B <int>        Specifies the target number of fragments per partition.\n"
            "\t-C <int>        Check point level for restart. (not implemented)\n"
            "\t-D <int>        Specify a debugging level. (not implemented)\n"
            "\t-E <int>        Specify the maximum number of soft errors.\n"
            "\t-F <directory>  The fragment store name.\n"
            "\t-G <file>       A file of fragment IIDs.\n"
            "\t-H <filename>   chimeras file.\n"
            "\t-I <directory>  Read the OVL store.\n"
            "\t-J <boolean>    Remove blizzard overlaps in FGB.\n"
            "\t-K <filename>   File of known fragment-end based branch-points.\n"
            "\t-L <filename>   The input OverlapFragMesgs; asm.ofg.\n"
            "\t-N <maxIID>\n"
            "\t-P              Specify ASCII output for protoIO messages.\n"
            "\t-Q <int>        Specify Reaper Pass.\n"
            "\t-R <int>        Which containment rule to use. \n"
            "\t-S <filename>   Spurs file.\n"
            "\t-T              Output transitively removeable overlaps.\n"
            "\t-U <boolean>    Find bubble smoothing overlaps.\n"
            "\t-V              Transitive edge reduction restarts at this VID.\n"
            "\t-W <int>        Limit in path length for graph walking.\n"
            "\t-Y <boolean>    Do not count chimera fragments.\n"
            "\t-Z <boolean>    Remove the blizzard overlaps in CGB.\n"
            "\t-b <int>        Number of cgb passes for finding branch points.\n"
            "\t-d <int>        Enable/Disable de-chording of the fragment overlap graph.\n"
            "\t-e <n>          Overlaps with error rate about this are ignored on input.\n"
            "\t\t                An integer value is in parts per thousand.\n"
            "\t-h              Help.\n"
            "\t-j <int>        Unique unitig cut-off\n"
            "\t-l <int>        Specify length of the genome.\n"
            "\t-m <nEdge>      Pre-allocate memory\n"
            "\t-n <nFrag>      Pre-allocate memory\n"
            "\t-o <pfx>        output to this prefix.\n"
            "\t-p <file>       Specify the parameters file.\n"
            "\t-q              Do not find branch points.\n"
            "\t-s              Disable early spur fragment removal.\n"
            "\t-u <filename>   Create a OVL compatible dump of the graph.\n"
            "\t-v <int>        Specify a verbosity level.\n"
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
  TStateGlobals    *gstate = (TStateGlobals     *)safe_calloc(sizeof(TStateGlobals), 1);
  THeapGlobals     *heapva = (THeapGlobals      *)safe_calloc(sizeof(THeapGlobals), 1);
  UnitiggerGlobals *rg     = (UnitiggerGlobals  *)safe_calloc(sizeof(UnitiggerGlobals), 1);

  rg->program_name = argv[0];
  rg->work_limit_per_candidate_edge = 1000;
  rg->remove_blizzard_overlaps = FALSE;
  rg->use_all_overlaps_in_reaper_pass = TRUE;
  rg->as_cgb_max_frag_iid = 100000 * CGB_MULTIPLIER;
  rg->dvt_double_sided_threshold_fragment_end_degree = 1;
  rg->con_double_sided_threshold_fragment_end_degree = 1;
  rg->intrude_with_non_blessed_overlaps_flag = FALSE;
  rg->cutoff_fragment_end_degree = 1000000;
  rg->dechord_the_graph = TRUE;
  rg->check_point_level = 0;
  rg->compress_the_graph = TRUE;
  rg->cgb_unique_cutoff = CGB_UNIQUE_CUTOFF;
  rg->walk_depth = 100;
  rg->iv_start = 0;
  rg->overlap_error_threshold = AS_OVS_encodeQuality(1.0);
  rg->recalibrate_global_arrival_rate = FALSE;
  rg->work_limit_placing_contained_fragments = 20;
  rg->walk_depth=100;
  rg->iv_start=0;
  rg->output_iterations_flag = TRUE;
  rg->aggressive_spur_fragment_marking = TRUE;

  rg->maxfrags = 40000;
  rg->maxedges = 40000;
  rg->maxtext  = 40000;

  ParseCommandLine( rg, argc, argv);

  //BasicUnitigger( argc, argv, gstate, heapva, rg);

 again:
  gstate->store_version       = 1;
  gstate->state_of_the_store  = 1; // Just initialized.
  gstate->min_frag_iid        = 0;
  gstate->max_frag_iid        = 0;
  gstate->nbase_in_genome     = 0;

  heapva->frags             = CreateVA_Afragment (rg->maxfrags);
  heapva->edges             = CreateVA_Aedge     (rg->maxedges); 
  heapva->frag_annotations  = CreateVA_char      (rg->maxtext);
  heapva->next_edge_obj     = CreateVA_IntEdge_ID(rg->maxedges); 

  heapva->chunkfrags = CreateVA_AChunkFrag(0);
  heapva->thechunks  = CreateVA_AChunkMesg(0);

  heapva->chunkseqs = CreateVA_char(0);
  heapva->chunkquas = CreateVA_char(0);
  heapva->chunksrcs = CreateVA_char(0);

  main_fgb(gstate, heapva, rg);
  main_cgb(gstate, heapva, rg);

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

    AS_CGB_Bubble_List_t bubbles = NULL, bubbles_tmp = NULL;

    static char bubble_boundaries_filename_tmp[500] = {0};
    static char bubble_overlaps_filename_tmp[500] = {0};

    rg->bubble_boundaries_filename = bubble_boundaries_filename_tmp;
    rg->bubble_overlaps_filename   = bubble_overlaps_filename_tmp;

    sprintf(rg->bubble_boundaries_filename,"%s.cgb_bubbles_txt",rg->Output_Graph_Store_Prefix);
    sprintf(rg->bubble_overlaps_filename, "%s.bubble_edges.ovl", rg->Output_Graph_Store_Prefix);

    // The OVL records needed to remove the bubbles
    GateKeeperStore *gkpStore = openGateKeeperStore(rg->frag_store, FALSE);
    FILE *bubble_overlaps_file = fopen(rg->bubble_overlaps_filename, "w");
    AS_CGB_Bubble_find_and_remove_bubbles(gkpStore,
                                          heapva->frags, heapva->edges,
                                          heapva->thechunks, heapva->chunkfrags, 
                                          gstate->global_fragment_arrival_rate,
                                          bubble_overlaps_file,
                                          stderr,
                                          rg->Output_Graph_Store_Prefix);
    fclose(bubble_overlaps_file);
    closeGateKeeperStore(gkpStore);

    /* NOTE: 0's in following call indicate use of defaults. */
  
    bubbles_tmp = bubbles = AS_CGB_Bubble_find_bubbles(heapva->frags, heapva->edges, 0, 0, 0);    

    FILE *bubble_boundaries_file = fopen(rg->bubble_boundaries_filename, "w");
    while (bubbles_tmp != NULL) {
      fprintf(bubble_boundaries_file, F_IID " %d " F_IID " %d\n", 
              get_iid_fragment(heapva->frags, bubbles_tmp->start),
              bubbles_tmp->start_sx,
              get_iid_fragment(heapva->frags, bubbles_tmp->end),
              bubbles_tmp->end_sx);
      bubbles_tmp = bubbles_tmp->next;
    }
    fclose(bubble_boundaries_file);


    //  Like I said, redo all our work.
    Delete_VA(heapva->frags);
    Delete_VA(heapva->edges);
    Delete_VA(heapva->next_edge_obj);
    Delete_VA(heapva->chunkfrags);
    Delete_VA(heapva->thechunks);
    Delete_VA(heapva->chunkseqs);
    Delete_VA(heapva->chunkquas); 
    Delete_VA(heapva->frag_annotations);
    Delete_VA(heapva->chunksrcs);

    goto again;
  }
  

  if(rg->analysis_level > 0) {
    char strtmp2[FILENAME_MAX];

    sprintf(strtmp2,"%s%s.%d",rg->Output_Graph_Store_Prefix, ".cga",rg->num_cgb_passes);
    FILE *fcga = fopen(strtmp2, "w");

    sprintf(strtmp2,"%s%s.%d",rg->Output_Graph_Store_Prefix,".cam",rg->num_cgb_passes);
    FILE *fcam = fopen(strtmp2, "w");

    sprintf(strtmp2,"%s%s.%d",rg->Output_Graph_Store_Prefix,".cus",rg->num_cgb_passes);
    FILE *fstat = fopen(strtmp2, "w");

    chunk_graph_analysis(rg->analysis_level,
                         gstate->max_frag_iid,
                         heapva->frags,  /* The internal representation of the fragment reads. */
                         heapva->edges,  /* The internal representation of the  overlaps. */
                         heapva->frag_annotations,
                         gstate->nbase_in_genome,
                         rg->recalibrate_global_arrival_rate,
                         rg->cgb_unique_cutoff,
                         gstate->global_fragment_arrival_rate,
                         rg->bubble_boundaries_filename,
                         heapva->chunkfrags,
                         /* Modify the chunk annotation, by setting the "isrc" index into 
                            chunksrc array. */
                         heapva->thechunks,
                         heapva->chunksrcs, 
                         fcga, fcam, fstat, stderr );

    fclose(fcga);
    fclose(fcam);
    fclose(fstat);
  }

  output_the_chunks(heapva->frags,
                    heapva->edges,
                    heapva->frag_annotations,
                    heapva->chunkfrags,
                    heapva->thechunks,
                    heapva->chunkseqs,
                    heapva->chunkquas,
                    heapva->chunksrcs,
                    rg->analysis_level,
                    gstate->global_fragment_arrival_rate,
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

    AChunkMesg *mychunk = GetVA_AChunkMesg(heapva->thechunks,
                                           get_cid_fragment(heapva->frags,avx));

    assert(nfrag > avx);

    if(mychunk->coverage_stat >= rg->cgb_unique_cutoff) {
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
  Delete_VA(heapva->next_edge_obj);
  Delete_VA(heapva->chunkfrags);
  Delete_VA(heapva->thechunks);
  Delete_VA(heapva->chunkseqs);
  Delete_VA(heapva->chunkquas); 
  Delete_VA(heapva->frag_annotations);
  Delete_VA(heapva->chunksrcs);

  fprintf(stderr, "All done.  Bye.\n");

  exit(0);
}
