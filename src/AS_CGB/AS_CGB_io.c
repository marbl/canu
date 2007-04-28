
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
= "$Id: AS_CGB_io.c,v 1.12 2007-04-28 08:46:21 brianwalenz Exp $";
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


void output_the_chunks
(/* Input Only*/
 const Tfragment frags[],
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
 const char * const Graph_Store_File_Prefix
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
