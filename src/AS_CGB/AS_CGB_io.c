
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
= "$Id: AS_CGB_io.c,v 1.1.1.1 2004-04-14 13:49:43 catmandew Exp $";
/* *******************************************************************
 *
 * Module: AS_CGB_io.c
 * 
 * Description: Chunk Graph Builder file input and output.
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

/************************************************************************/

#ifdef DEBUGGING
#define DEBUG_SHORTCUT
#undef DEBUG_SHORTCUT
#endif

/*************************************************************************/
#define FILENAME_MAX_LEN 1024
#define IGNORE (AS_NO_OVERLAP)
/*************************************************************************/

void output_the_chunks
(/* Input Only*/
 MesgWriter WriteMesg_AS,
 const Tfragment frags[],
 const Tedge     edges[],
 const VA_TYPE(char) * const fragsrc,
 const TChunkFrag    chunkfrags[],
 const TChunkMesg    thechunks[],
 const VA_TYPE(char) chunkseqs[],
 const VA_TYPE(char) chunkquas[],
 const VA_TYPE(char) chunksrc[],
 const int analysis_level,
 const int output_fom_messages,
 const float global_fragment_arrival_rate,
 const int fragment_count_target,
 const char * const Graph_Store_File_Prefix
#if 0
 /* Modified */
 VA_TYPE(IntMultiPos)   * the_imps,
 VA_TYPE(IntUnitigMesg) * the_iums
#endif
 )
{
  IntFragment_ID max_num_frags_per_chunk=0;
  IntEdge_ID max_num_ovlps_per_chunk=0;

  IntChunk_ID       chunk_index;
  const IntChunk_ID nchunks = (IntChunk_ID)GetNumVA_AChunkMesg(thechunks);
  const IntFragment_ID  nfrag   = GetNumFragments(frags);
  
  size_t  * frag_source_index = NULL; 
  VA_TYPE(char)  * newfragsrc = NULL;

  const int nsample=500;
  const int nbucket=500;

  char filename[FILENAME_MAX_LEN];
  
  FILE *fcgb = NULL;
  int fragment_count = 0;
  int file_count = 0;
  
#if 0
  ResetToRange_VA(the_imps,0);
  ResetToRange_VA(the_iums,0);
#else
 VA_TYPE(IntMultiPos) * the_imps = CreateVA_IntMultiPos(0);
#endif
  
  // Determine the maximum sizes of the variant parts of a chunk message.
  for(chunk_index=0;chunk_index < nchunks; chunk_index++){
    const AChunkMesg * const mychunk = GetVA_AChunkMesg(thechunks,chunk_index);

    max_num_frags_per_chunk = max(max_num_frags_per_chunk,
				  mychunk->num_frags);
    max_num_ovlps_per_chunk = max(max_num_ovlps_per_chunk,
				  mychunk->a_degree_raw 
				  + mychunk->b_degree_raw);
  }

  if(analysis_level > 1 && (fragsrc != NULL)) {

    SAFE_MALLOC(frag_source_index, size_t, nfrag);
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

#if 0
    fprintf(stdout, __FILE__ " coverage stat " F_IID " %f %f\n",
            chunk_index, mychunk->coverage_stat, coverage_statistic);
#endif
    
    //assert(mychunk->f_list >= 0);
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
      a_frag.source =   ( (analysis_level == 0)
			  ? NULL
			  :  ( (NULL != frag_source_index)
			       ? GetVA_char(newfragsrc,frag_source_index[vid])
			       : GetVA_char(fragsrc, get_src_fragment(frags,vid)) 
			       ));
      a_frag.position.bgn = get_o5p_fragment(frags,vid);
      a_frag.position.end = get_o3p_fragment(frags,vid);
      a_frag.delta_length = 0;
      a_frag.delta        = NULL;

      SetVA_IntMultiPos(the_imps,ivc,&a_frag);
    }

    {
      //achunk.consensus      = GetVA_char(chunkseqs,mychunk->seq_loc);
      //achunk.quality        = GetVA_char(chunkquas,mychunk->qua_loc);
      achunk.consensus = "";
      achunk.quality   = "";
    }
      
    // The ProtoSpec specifies that the first chunk id is ZERO.
    achunk.iaccession     = mychunk->iaccession;
    achunk.source         = GetVA_char(chunksrc,mychunk->isrc);
    achunk.coverage_stat  = coverage_statistic;
    achunk.status         = AS_UNASSIGNED;
    achunk.a_branch_point = mychunk->a_branch_point;
    achunk.b_branch_point = mychunk->b_branch_point;

    {
      achunk.consensus      = "";
      achunk.quality        = "";
      achunk.length         = mychunk->bp_length;
      //achunk.length         = strlen(achunk.consensus);
    }
    achunk.forced         = forced;
    achunk.num_frags      = mychunk->num_frags;
    achunk.f_list         = GetVA_IntMultiPos(the_imps,0);

    fragment_count += mychunk->num_frags;
    if(fragment_count_target > 0) {
      if((fragment_count >= fragment_count_target)) {
        if(NULL != fcgb) { fclose(fcgb); fcgb = NULL;}
        fragment_count = 0;
      }
      if(NULL == fcgb) {
        file_count ++;
        sprintf(filename,"%s_%03d.cgb_tmp",Graph_Store_File_Prefix,file_count);
        fcgb = fopen(filename,"w");
        assert(NULL != fcgb);
      }
    } else {
      if(NULL == fcgb) {
        sprintf(filename,"%s.cgb_tmp",Graph_Store_File_Prefix);
        fcgb = fopen(filename,"a");
	// We are depending on the ADT messages having just been
	// written to a new file with a file name of filename.
        assert(NULL != fcgb);
      }
    }
    
    pmesg.t = MESG_IUM;
    pmesg.m = &achunk;
    WriteMesg_AS(fcgb,&pmesg);
    //fprintf(fpar,"% 5d % 10" F_IIDP "\n", file_count, iid);

  }

  DeleteVA_IntMultiPos(the_imps);

  if(fragment_count_target > 0) {
    if(NULL != fcgb) { fclose(fcgb); fcgb = NULL;}
    fragment_count = 0;
    if(NULL == fcgb) {
      file_count ++;
      sprintf(filename,"%s_%03d.cgb_tmp",Graph_Store_File_Prefix,file_count);
      fcgb = fopen(filename,"w");
      assert(NULL != fcgb);
    }
  }

  if( NULL != frag_source_index ) {
    assert(NULL != newfragsrc);
    DeleteVA_char(newfragsrc); newfragsrc = NULL;
    SAFE_FREE(frag_source_index);
  }

  if( output_fom_messages ) {
    /* Output all the fragment overlaps between unitigs in a form that
       the CGW wants. */

    Histogram_t  * fom_types_histogram = NULL;
    IntFragment_ID vid;
    int vsx;

    if(analysis_level > 0) {
      fom_types_histogram = create_histogram(nsample,nbucket,TRUE,FALSE);
    }

    for( vid=0; vid < nfrag; vid++ ) {
      for( vsx=0; vsx<2 ; vsx ++ ) {

	const IntEdge_ID frag_edge_list = get_segstart_vertex(frags,vid,vsx);
	const int frag_edge_degree = get_seglen_vertex(frags,vid,vsx);
	GenericMesg   pmesg;
	
	int ifoe;
	for(ifoe=0;ifoe < frag_edge_degree; ifoe++){
	
	  const IntEdge_ID ie = frag_edge_list + ifoe;
	  /* A raw fragment overlap edge */
	  const Tnes inese = get_nes_edge(edges,ie);
	
	  const IntFragment_ID iavx = get_avx_edge(edges,ie);
	  // An index into the Tfragment array of the fragment at the
	  // proximal vertex of the edge.
	  const int iasx = get_asx_edge(edges,ie);
	  // A flag indicating whether the suffix of the fragment at
	  // the proximal vertex is in the overlap.
	  const IntFragment_ID ibvx = get_bvx_edge(edges,ie);
	  // An index into the Tfragment array of the fragment at the
	  // distal vertex of the edge.
	  const int ibsx = get_bsx_edge(edges,ie);
	  // A flag indicating whether the suffix of the fragment at
	  // the distant vertex is in the overlap.
	  const int iahg = get_ahg_edge(edges,ie); /* overhang in bp. */
#ifndef STORE_OVERLAP_EXTREMES
          const int iamn = get_ahg_edge(edges,ie);
	  const int iamx = get_ahg_edge(edges,ie);
#else // STORE_OVERLAP_EXTREMES
          const int iamn = get_amn_edge(edges,ie);
	  const int iamx = get_amx_edge(edges,ie);
#endif // STORE_OVERLAP_EXTREMES
	  const int ibhg = get_bhg_edge(edges,ie); /* overhang in bp. */
	  
	  const int ilen = get_length_fragment(frags,iavx);
	  const int best_overlap_length = (ilen - iahg);
	  const int min_overlap_length  = (ilen - iamx);
	  const int max_overlap_length  = (ilen - iamn);
	  const float32 quality         = get_qua_edge(edges,ie);
	  
          const int containment = ! is_a_dvt_simple(iahg,ibhg);
          //const Tlab alab = get_lab_fragment(frags,iavx);
          //const Tlab blab = get_lab_fragment(frags,ibvx);
          const IntChunk_ID cavx = get_cid_fragment(frags,iavx);
          const IntChunk_ID cbvx = get_cid_fragment(frags,ibvx);
          
          assert( iavx == vid );
	  assert( iasx == vsx );
	  
          if( (cavx != cbvx ) &&
              // Do not output overlaps between fragments in the same unitig.
              
              ( containment ? (iasx == TRUE) : (iavx > ibvx) )
              // Output only one edge per overlap. We use the
              // knowledge that each contaiment overlap has two edges
              // from the container fragment to the contained
              // fragment.
              ) {
	    
	    FragOverlapMesg cea;
	    UnitigOverlapType overlap_type;
#if 0
	    const int casx = (iasx ^ (! get_forward_fragment(frags,iavx)));
	    const int cbsx = (ibsx ^ (! get_forward_fragment(frags,ibvx)));
#endif
	    
	    switch(inese) {
	    case AS_CGB_INTERCHUNK_EDGE:
	      overlap_type = AS_OVERLAP; break;
	    case AS_CGB_INTRACHUNK_EDGE:
	      assert(FALSE); break;
	    case AS_CGB_TOUCHES_CONTAINED_EDGE:
	      overlap_type = AS_TOUCHES_CONTAINED_OVERLAP; break;
	    case AS_CGB_BETWEEN_CONTAINED_EDGE:
	      overlap_type = IGNORE; break;
	      // overlap_type = AS_BETWEEN_CONTAINED_OVERLAP; break;
	    case AS_CGB_TOUCHES_CRAPPY_DVT:
	    case AS_CGB_BETWEEN_CRAPPY_DVT:
	      overlap_type = AS_OVERLAP; break;
#if 1
	    case AS_CGB_MARKED_BY_BRANCH_DVT:
	    case AS_CGB_MARKED_BY_BREAKER:
	      overlap_type = IGNORE; break;
#endif                

	    case AS_CGB_CONTAINED_EDGE:
	    case AS_CGB_TOUCHES_CRAPPY_CON:
	    case AS_CGB_BETWEEN_CRAPPY_CON:
              if(is_a_toc_simple(iahg,ibhg)) {
                overlap_type = AS_1_CONTAINS_2_OVERLAP;
              } else if ( is_a_frc_simple(iahg,ibhg)) {
                overlap_type = AS_2_CONTAINS_1_OVERLAP;
              } else {
                assert(FALSE);
              }
              break;
	    case AS_CGB_REMOVED_BY_TRANSITIVITY_CON: // A chord in the graph.
              if(is_a_toc_simple(iahg,ibhg)) {
                overlap_type = AS_1_CONTAINS_2_CHORD_OVERLAP;
              } else if ( is_a_frc_simple(iahg,ibhg)) {
                overlap_type = AS_2_CONTAINS_1_CHORD_OVERLAP;
              } else {
                assert(FALSE);
              }
              break;
	      
	    case AS_CGB_MARKED_BY_DELETED_DVT:
	    case AS_CGB_MARKED_BY_DELETED_CON:
	      overlap_type = IGNORE; break;
	      
	    default:
	      fprintf(stderr,"Unexpected overlap edge type = %d\n",
		      inese);
	      assert(FALSE);
	    }
	    
	    if( overlap_type != IGNORE) {
	      
	      { // The following are fragment overlaps.
#if 0
		const ChunkOrientationType  orient
		  = ( ( iasx) && ( ibsx) ? AB_BA : 0)
		  + ( ( iasx) && (!ibsx) ? AB_AB : 0)
		  + ( (!iasx) && ( ibsx) ? BA_BA : 0)
		  + ( (!iasx) && (!ibsx) ? BA_AB : 0);
#else
		const ChunkOrientationType  orient
		  = ( iasx 
		      ? (ibsx ? AB_BA : AB_AB )
		      : (ibsx ? BA_BA : BA_AB )
		      );
#endif
		char source[80]= {0};
#if 0
		sprintf(source,
                        ">>(" F_IID ",%d)(" F_IID ",%d) %d",
                        cavx,casx,cbvx,cbsx,inese);
#endif
		cea.afrag = get_iid_fragment(frags,iavx);
		cea.bfrag = get_iid_fragment(frags,ibvx);
		cea.orient = orient;
		cea.overlap_type = overlap_type;
		cea.best_overlap_length = best_overlap_length;
		cea.min_overlap_length  = min_overlap_length;
		cea.max_overlap_length  = max_overlap_length;
		cea.quality             = quality;
		cea.source = source;
		
	      }
	      
	      if( NULL != fom_types_histogram ) {
		add_to_histogram(fom_types_histogram, (int)overlap_type, NULL);
	      }
	      
	      pmesg.t = MESG_FOM;
	      pmesg.m = &cea;
	      WriteMesg_AS(fcgb,&pmesg);
	      
	    }
          }
        }
      }
    }

    if(NULL != fom_types_histogram) {
      fprintf(stderr,"\n\nHistogram of the FOM types\n");
      print_histogram(stderr, fom_types_histogram, 0, 1);
      free_histogram(fom_types_histogram);
    }
  }

  fclose(fcgb);

  {
    int iii = 0;
    int ierr;

    if(fragment_count_target > 0) {
      for( iii=0; iii <= file_count ; iii++) {
        char thePath1[CMD_BUFFER_SIZE-1]={0};
	sprintf(filename,"%s_%03d.cgb",Graph_Store_File_Prefix,iii);
        sprintf(thePath1,"%s_tmp",filename);
        ierr = accept_tmp_as_final_file( thePath1, filename);
        assert(ierr == 0);
        // The temporary batch info file could not be moved into its final
        // position.
      }
    } else {
      char thePath1[CMD_BUFFER_SIZE-1]={0};
      sprintf(filename,"%s.cgb", Graph_Store_File_Prefix);
      sprintf(thePath1,"%s_tmp",filename);
      ierr = accept_tmp_as_final_file( thePath1, filename);
      assert(ierr == 0);
      // The temporary batch info file could not be moved into its final
      // position.
    }
  }
}


void convert_the_chunks_to_IUM
(/* Input Only*/
 const Tfragment        * frags,
 const Tedge            * edges,
 const VA_TYPE(char)    * fragsrc,
 const TChunkFrag       * chunkfrags,
 const TChunkMesg       * thechunks,
 const VA_TYPE(char)    * chunkseqs,
 const VA_TYPE(char)    * chunkquas,
 const VA_TYPE(char)    * chunksrc,
 const int                analysis_level,
 const float              global_fragment_arrival_rate,
 /* Filled as output only */
 VA_TYPE(IntMultiPos)   * the_imps,
 VA_TYPE(IntUnitigMesg) * the_iums,
 VA_TYPE(char)   * the_imp_source,
 VA_TYPE(char)   * the_ium_source
)
{
  const IntChunk_ID nchunks = (IntChunk_ID)GetNumVA_AChunkMesg(thechunks);
  const IntFragment_ID  nfrag   = GetNumFragments(frags);
  
  int           fragment_count = 0;
  IntChunk_ID   chunk_index = 0;
  
  EnableRangeVA_IntMultiPos( the_imps, nfrag);
  EnableRangeVA_IntUnitigMesg( the_iums, nchunks);
  EnableRangeVA_char(the_imp_source,0);
  EnableRangeVA_char(the_ium_source,0);
                     
  for(chunk_index=0;chunk_index < nchunks; chunk_index++) {
    const AChunkMesg * const mychunk = GetVA_AChunkMesg(thechunks,chunk_index);
    const IntFragment_ID num_frags = mychunk->num_frags;
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

    //assert(mychunk->f_list >= 0);
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
      a_frag.source =   ( (analysis_level == 0)
			  ? NULL
			  : GetVA_char(fragsrc, get_src_fragment(frags,vid)) );
      a_frag.position.bgn = get_o5p_fragment(frags,vid);
      a_frag.position.end = get_o3p_fragment(frags,vid);
      a_frag.delta_length = 0;
      a_frag.delta        = NULL;

      SetVA_IntMultiPos(the_imps,fragment_count+ivc,&a_frag);
    }

    {
      //achunk.consensus      = GetVA_char(chunkseqs,mychunk->seq_loc);
      //achunk.quality        = GetVA_char(chunkquas,mychunk->qua_loc);
      achunk.consensus = "";
      achunk.quality   = "";
    }
      
    // The ProtoSpec specifies that the first chunk id is ZERO.
    achunk.iaccession     = mychunk->iaccession;
    achunk.source         = GetVA_char(chunksrc,mychunk->isrc);
    achunk.coverage_stat  = coverage_statistic;
    achunk.status         = AS_UNASSIGNED;
    achunk.a_branch_point = mychunk->a_branch_point;
    achunk.b_branch_point = mychunk->b_branch_point;

    {
      achunk.consensus      = "";
      achunk.quality        = "";
      achunk.length         = mychunk->bp_length;
      //achunk.length         = strlen(achunk.consensus);
    }
    achunk.forced         = forced;
    achunk.f_list         = GetVA_IntMultiPos(the_imps,fragment_count);
    achunk.num_frags      = mychunk->num_frags;

    SetVA_IntUnitigMesg( the_iums, chunk_index, &achunk);

    fragment_count += mychunk->num_frags;
  }

  // Warning: Convert the offsets into the_imps into C-pointers.
}


void output_the_IUM_to_file
(/* Input Only*/
 MesgWriter                  WriteMesg_AS,
 const VA_TYPE(char)    *    fragsrc,
 const VA_TYPE(char)    *    chunksrc,
 VA_TYPE(IntMultiPos)   *    the_imps,
 VA_TYPE(IntUnitigMesg) *    the_iums,
 const int                   fragment_count_target,
 const char * const          Graph_Store_File_Prefix
)
{
  const IntChunk_ID nchunks = (IntChunk_ID)GetNumVA_IntUnitigMesg(the_iums);
  char          filename[FILENAME_MAX_LEN] = {0};
  FILE *        fcgb = NULL;
  int           fragment_count = 0;
  int           file_count = 0;
  IntChunk_ID   chunk_index = 0;

  for(chunk_index=0;chunk_index < nchunks; chunk_index++) {
    IntUnitigMesg * mychunk
      = GetVA_IntUnitigMesg(the_iums,chunk_index);

    fragment_count += GetVA_IntUnitigMesg(the_iums,chunk_index)->num_frags;

    if(fragment_count_target > 0) {
      if((fragment_count >= fragment_count_target)) {
        if(NULL != fcgb) { fclose(fcgb); fcgb = NULL;}
        fragment_count = 0;
      }
      if(NULL == fcgb) {
        file_count ++;
        sprintf(filename,"%s_%03d.cgb_tmp",Graph_Store_File_Prefix,file_count);
        fcgb = fopen(filename,"w");
        assert(NULL != fcgb);
      }
    } else {
      if(NULL == fcgb) {
        sprintf(filename,"%s.cgb_tmp",Graph_Store_File_Prefix);
        fcgb = fopen(filename,"a");
	// We are depending on the ADT messages having just been
	// written to a new file with a file name of filename.
        assert(NULL != fcgb);
      }
    }

    {
      GenericMesg   pmesg;
      pmesg.t = MESG_IUM;
      pmesg.m = mychunk;
      WriteMesg_AS(fcgb,&pmesg);
      //fprintf(fpar,"% 5d % 10" F_IIDP "\n", file_count, iid);
    }
  }

  if(fragment_count_target > 0) {
    if(NULL != fcgb) { fclose(fcgb); fcgb = NULL;}
    fragment_count = 0;
    if(NULL == fcgb) {
      file_count ++;
      sprintf(filename,"%s_%03d.cgb_tmp",Graph_Store_File_Prefix,file_count);
      fcgb = fopen(filename,"w");
      assert(NULL != fcgb);
    }
  }

  fclose(fcgb);

  {
    int iii = 0;
    int ierr;

    if(fragment_count_target > 0) {
      for( iii=0; iii <= file_count ; iii++) {
        char thePath1[CMD_BUFFER_SIZE-1]={0};
	sprintf(filename,"%s_%03d.cgb",Graph_Store_File_Prefix,iii);
        sprintf(thePath1,"%s_tmp",filename);
        ierr = accept_tmp_as_final_file( thePath1, filename);
        assert(ierr == 0);
        // The temporary batch info file could not be moved into its final
        // position.
      }
    } else {
      char thePath1[CMD_BUFFER_SIZE-1]={0};
      sprintf(filename,"%s.cgb", Graph_Store_File_Prefix);
      sprintf(thePath1,"%s_tmp",filename);
      ierr = accept_tmp_as_final_file( thePath1, filename);
      assert(ierr == 0);
      // The temporary batch info file could not be moved into its final
      // position.
    }
  }

}


