
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
/* *******************************************************************
 * $Id: AS_CGB_edgetrimmer.c,v 1.2 2004-09-23 20:25:01 mcschatz Exp $
 *
 * Module: AS_CGB_edgetrimmer.c
 * 
 * Description: 
 *
 * Assumptions: 
 *
 * Authors: Clark Mobarry
 *
 * This module will trim spurious edges for ALU spanning.
 *
 *********************************************************************/

/*************************************************************************/
/* System include files */

/*************************************************************************/
/* Local include files */
#include "AS_CGB_all.h"

/*************************************************************************/
static char CM_ID[] 
= "$Id: AS_CGB_edgetrimmer.c,v 1.2 2004-09-23 20:25:01 mcschatz Exp $";
/*************************************************************************/

#define DEBUG_WITH_HISTO
#undef DEBUG_WITH_HISTO

#define DEBUGGING
#undef DEBUGGING

#ifdef DEBUGGING
#endif // DEBUGGING

#define AS_CGB_BPT_SLOP 10
//#define AS_CGB_BPT_IN_UNIQUE 40
#define AS_CGB_BPT_IN_UNIQUE 30 
//#define AS_CGB_BPT_IN_UNIQUE 20

void chunk_end_edge_trimmer
(
 /* input only */
 const float         cgb_unique_cutoff,
 const float         global_fragment_arrival_rate,
 const Tfragment     frags[],
 const TChunkFrag    chunkfrags[],
 /* modify */
 Tedge               edges[],
 TChunkMesg          thechunks[]
)
{
#ifdef DEBUG_WITH_HISTO
  FILE * fout = stdout;
  const int nsample = 1000, nbucket = 1000;
  HISTOGRAM 
    *histo_unique_region_overlaps = create_histogram(nsample,nbucket,0,TRUE),
    *histo_repeat_region_overlaps = create_histogram(nsample,nbucket,0,TRUE),
    *histo_unknown_region_overlaps = create_histogram(nsample,nbucket,0,TRUE),
    *histo_removeable_overlaps = create_histogram(nsample,nbucket,0,TRUE);
#endif /*DEBUG_WITH_HISTO*/

  const IntChunk_ID nchunks = (IntChunk_ID)GetNumVA_AChunkMesg(thechunks);
  IntChunk_ID seed_chunk_index;
  int seed_chunk_suffix;
  IntChunk_ID count_of_u_unitig_ends = 0;
  // int iter;

  /* Loop over chunk ends as vertices of the chunk graph. */
  for(seed_chunk_index=0;
      seed_chunk_index < nchunks; seed_chunk_index++) 
    for(seed_chunk_suffix=0; seed_chunk_suffix<2; seed_chunk_suffix++) {
      // const IntChunk_ID seed_chunk_key = 2*seed_chunk_index+seed_chunk_suffix;
      
      /* The return value is the number of branch partners including the
	 seed chunk-end. */
      AChunkMesg * the_seed_chunk = GetVA_AChunkMesg(thechunks,seed_chunk_index);
  
      /* From the seed chunk in the gang, walk the adjacent vertices. */
      
      const IntEdge_ID seed_chunkedgelist
	= ( seed_chunk_suffix == 0 ? 
	    the_seed_chunk->a_list_raw :
	    the_seed_chunk->b_list_raw );
      const int seed_chunkedgedegree 
	= ( seed_chunk_suffix == 0 ? 
	    the_seed_chunk->a_degree_raw :
	    the_seed_chunk->b_degree_raw );
      const int seed_chunk_bpt 
	= ( seed_chunk_suffix == 0 ?
	    the_seed_chunk->a_branch_point :
	    the_seed_chunk->b_branch_point );
      const int number_of_randomly_sampled_fragments_in_chunk
        = count_the_randomly_sampled_fragments_in_a_chunk
        ( frags, chunkfrags, thechunks, seed_chunk_index);
      const BPTYPE rho = GetVA_AChunkMesg(thechunks,seed_chunk_index)->rho;
      const float coverage_stat = compute_coverage_statistic
        ( rho,
          number_of_randomly_sampled_fragments_in_chunk,
          global_fragment_arrival_rate );

#ifdef DEBUG33
      printf("seed_chunk_index,seed_chunk_suffix=" F_IID",%d"
	     "a_branch_type=%d,a_branch_point=" F_S32 ","
	     "b_branch_type=%d,b_branch_point=" F_S32 "\n",
	     seed_chunk_index,seed_chunk_suffix,
	     GetVA_AChunkMesg(thechunks,seed_chunk_index)->a_branch_type,
	     GetVA_AChunkMesg(thechunks,seed_chunk_index)->a_branch_point,
	     GetVA_AChunkMesg(thechunks,seed_chunk_index)->b_branch_type,
	     GetVA_AChunkMesg(thechunks,seed_chunk_index)->b_branch_point);
#endif

      if( coverage_stat >= cgb_unique_cutoff ) {
	// assert that the adjacency list is already sorted by
	// increasing ahg.
	int count_of_unique_region_overlaps = 0;
	int count_of_repeat_region_overlaps = 0;
	int count_of_unknown_region_overlaps = 0;
	int count_of_removeable_overlaps = 0;
	int imate_index, imate_used;
      
	count_of_u_unitig_ends++;
	
	/* From the seed chunk-end find its mate chunk-ends. */
	for(imate_index=0, imate_used=0;
	    (imate_index < seed_chunkedgedegree);
	    imate_index++) {
	  // The adjacency list is already sorted by overhang when
	  // the original fragment adjacency list, that is from
	  // smallest to largest. We are taking just the MAX_NUM_OF_MATES
	  // 
	  const IntEdge_ID mate_edge = seed_chunkedgelist + imate_index;
	  const Tnes       mate_edge_nes = get_nes_edge(edges,mate_edge);

	  // Process the backbone dovetail overlap types used in
	  // count_number_of_backbone_dovetails().
          if( AS_CGB_INTRACHUNK_EDGE == mate_edge_nes ) {
            IntFragment_ID avx = get_avx_edge(edges,mate_edge);
            IntFragment_ID bvx = get_bvx_edge(edges,mate_edge);
            int asx = get_asx_edge(edges,mate_edge);
            int bsx = get_bsx_edge(edges,mate_edge);
            fprintf(stderr,
                    "AS_CGB_INTRACHUNK_EDGE == mate_edge_nes\n"
                    "aid=" F_IID " asx=%d bid=" F_IID " bsx=%d\n"
                    "seed_chunk_bpt=%d (alen-seed_chunk_bpt)=" F_S32 "\n",
                    get_iid_fragment(frags,avx), asx,
                    get_iid_fragment(frags,bvx), bsx,
                    seed_chunk_bpt,(get_length_fragment(frags,avx)-seed_chunk_bpt)
                    );
          }

	  //assert( AS_CGB_INTRACHUNK_EDGE != mate_edge_nes );
          //I forgot about circular chunks!
          
	  if( (AS_CGB_INTERCHUNK_EDGE == mate_edge_nes) ||
              (AS_CGB_TOUCHES_CONTAINED_EDGE == mate_edge_nes)
              ) {
	    const int mate_chunk_ovl = get_best_ovl(frags,edges,mate_edge);
	    
	    imate_used++;
	    assert(AS_CGB_BPT_IN_UNIQUE > AS_CGB_BPT_SLOP);
	    if(mate_chunk_ovl > seed_chunk_bpt + AS_CGB_BPT_IN_UNIQUE) {
	      // This is an anchor edge.
	      count_of_unique_region_overlaps ++;
	    } else if (mate_chunk_ovl <= seed_chunk_bpt + AS_CGB_BPT_SLOP) {
	      // This is a trimmable edge if there is a anchor edge.
	      count_of_repeat_region_overlaps ++;
	    } else {
	      // This edge is in the uncertain region.
	      count_of_unknown_region_overlaps ++;
	    }
	    
	  }
	}
	
#ifdef DEBUGGING
	if( (count_of_unique_region_overlaps == 1 ) &&
	    (count_of_repeat_region_overlaps  > 0 ) &&
	    (count_of_unknown_region_overlaps > 0 ) 
	    ) {
	  fprintf(stderr,"A foiled edge trimming! \n"
		  "seed_chunk_index=" F_IID ", seed_chunk_suffix=%d,"
		  "count_of_unique_region_overlaps=%d,"
		  "count_of_repeat_region_overlaps=%d,"
		  "count_of_unknown_region_overlaps=%d\n",
		  seed_chunk_index,seed_chunk_suffix,
		  count_of_unique_region_overlaps,
		  count_of_repeat_region_overlaps,
		  count_of_unknown_region_overlaps);
	}
#endif // DEBUGGING
	
#ifdef DEBUG_WITH_HISTO
	add_to_histogram(histo_unique_region_overlaps,
			 count_of_unique_region_overlaps,NULL);
	add_to_histogram(histo_repeat_region_overlaps,
			 count_of_repeat_region_overlaps,NULL);
	add_to_histogram(histo_unknown_region_overlaps,
			 count_of_unknown_region_overlaps,NULL);
#endif /*DEBUG_WITH_HISTO*/
	
	// At this point we have a U-unitig with overlaps in the 
	if( (count_of_unique_region_overlaps > 0) &&
	    (count_of_repeat_region_overlaps > 0) ) 
	  for(imate_index=0, imate_used=0;
	      (imate_index < seed_chunkedgedegree);
	      imate_index++) {
	    const IntEdge_ID mate_edge = seed_chunkedgelist + imate_index;
	    const Tnes mate_edge_nes = get_nes_edge(edges,mate_edge);
	    // The adjacency list is already sorted by overhang when
	    // the original fragment adjacency list, that is from
	    // smallest to largest. We are taking just the MAX_NUM_OF_MATES
	    // 
	    if( (AS_CGB_INTERCHUNK_EDGE == mate_edge_nes) ||
                (AS_CGB_TOUCHES_CONTAINED_EDGE == mate_edge_nes)
                ) {
	      const int mate_chunk_ovl = get_best_ovl(frags,edges,mate_edge);
	      imate_used++;
	      
	      assert(AS_CGB_BPT_IN_UNIQUE > AS_CGB_BPT_SLOP);
	      if(mate_chunk_ovl > seed_chunk_bpt + AS_CGB_BPT_IN_UNIQUE) {
		// Keep this edge.
	      } else if (mate_chunk_ovl <= seed_chunk_bpt + AS_CGB_BPT_SLOP) {
		// Trim this edge.
		count_of_removeable_overlaps ++;
		set_nes_edge(edges,mate_edge,AS_CGB_MARKED_BY_BRANCH_DVT);
		fix_overlap_edge_mate(frags,edges,mate_edge);
	      } else {
		// What to do with this edge?
	      }
	    }
	  }

#ifdef DEBUG_WITH_HISTO
	  add_to_histogram(histo_removeable_overlaps,
			   count_of_removeable_overlaps,NULL);
#endif // DEBUG_WITH_HISTO
	}
    }
#ifdef DEBUG_WITH_HISTO
  fprintf(fout,"\n\nHistogram of unique_region_overlaps\n");
  print_histogram(fout,histo_unique_region_overlaps, 0, 1);

  fprintf(fout,"\n\nHistogram of repeat_region_overlaps\n");
  print_histogram(fout,histo_repeat_region_overlaps, 0, 1);

  fprintf(fout,"\n\nHistogram of histo_unknown_region_overlaps\n");
  print_histogram(fout,histo_unknown_region_overlaps, 0, 1);

  fprintf(fout,"\n\nHistogram of histo_removeable_overlaps\n");
  print_histogram(fout,histo_removeable_overlaps, 0, 1);

  free_histogram(histo_unique_region_overlaps);
  free_histogram(histo_repeat_region_overlaps);
  free_histogram(histo_unknown_region_overlaps);
  free_histogram(histo_removeable_overlaps);
#endif /*DEBUG_WITH_HISTO*/
}
	      

void fragment_end_edge_trimmer
(
 /* input only */
 const Tfragment     frags[],
 /* modify */
 Tedge               edges[]
)
{
#ifdef DEBUG_WITH_HISTO
  FILE * fout = stdout;
  const int nsample = 1000, nbucket = 1000;
  HISTOGRAM 
    *histo_unique_region_overlaps = create_histogram(nsample,nbucket,0,TRUE),
    *histo_repeat_region_overlaps = create_histogram(nsample,nbucket,0,TRUE),
    *histo_unknown_region_overlaps = create_histogram(nsample,nbucket,0,TRUE),
    *histo_removeable_overlaps = create_histogram(nsample,nbucket,0,TRUE);
#endif /*DEBUG_WITH_HISTO*/

  // int iter;
  const IntFragment_ID nfrag = GetNumFragments(frags);
  IntFragment_ID seed_fragment_index;
  int seed_fragment_suffix;
  
  /* Loop over fragment ends as vertices of the fragment overlap graph. */
  for(seed_fragment_index=0; seed_fragment_index < nfrag; seed_fragment_index++) 
    for(seed_fragment_suffix=0; seed_fragment_suffix<2; seed_fragment_suffix++) {
      // const IntFragment_ID seed_fragment_key = 2*seed_fragment_index+seed_fragment_suffix;
      
      /* The return value is the number of branch partners including the
	 seed fragment-end. */
      /* From the seed fragment in the gang, walk the adjacent vertices. */
      
      const IntEdge_ID seed_fragment_edge_list =
        get_segstart_vertex( frags, seed_fragment_index, seed_fragment_suffix);
      const int seed_fragment_edge_degree =
        get_seglen_vertex( frags, seed_fragment_index, seed_fragment_suffix);

      int seed_bpt_vertex =
#ifndef STORE_BRANCH_POINTS_AT_FRAGMENT
	0;
#else // STORE_BRANCH_POINTS_AT_FRAGMENT
        get_bpt_vertex(frags, seed_fragment_index, seed_fragment_suffix);
#endif // STORE_BRANCH_POINTS_AT_FRAGMENT
      // A return value of zero is a sentinal for no branch-point.

      
      if( seed_bpt_vertex > 0 ) {
	// assert that the adjacency list is already sorted by
	// increasing ahg.
	int count_of_unique_region_overlaps = 0;
	int count_of_repeat_region_overlaps = 0;
	int count_of_unknown_region_overlaps = 0;
	int count_of_removeable_overlaps = 0;
	int imate_index, imate_used;
      
	/* From the seed fragment-end find its mate fragment-ends. */
	for(imate_index=0, imate_used=0;
	    (imate_index < seed_fragment_edge_degree); 
	    imate_index++) {
	  // The adjacency list is already sorted by overhang when
	  // the original fragment adjacency list, that is from
	  // smallest to largest. We are taking just the MAX_NUM_OF_MATES
	  // 
	  const IntEdge_ID mate_edge = seed_fragment_edge_list + imate_index;
	  const Tnes       mate_edge_nes = get_nes_edge(edges,mate_edge);

#if 0
	  // Process the backbone dovetail overlap types used in
	  // count_number_of_backbone_dovetails().
          if( AS_CGB_INTRACHUNK_EDGE == mate_edge_nes ) {
            IntFragment_ID avx = get_avx_edge(edges,mate_edge);
            IntFragment_ID bvx = get_bvx_edge(edges,mate_edge);
            int asx = get_asx_edge(edges,mate_edge);
            int bsx = get_bsx_edge(edges,mate_edge);
            fprintf(stderr,
                    "AS_CGB_INTRACHUNK_EDGE == mate_edge_nes\n"
                    "aid=" F_IID " asx=%d bid=" F_IID " bsx=%d\n"
                    "seed_fragment_bpt=%d (alen-seed_fragment_bpt)=" F_S32 "\n",
                    get_iid_fragment(frags,avx), asx,
                    get_iid_fragment(frags,bvx), bsx,
                    seed_bpt_vertex,(get_length_fragment(frags,avx)-seed_bpt_vertex)
                    );
          }
#endif
	  //assert( AS_CGB_INTRACHUNK_EDGE != mate_edge_nes );
          //I forgot about circular chunks!

	  if(
	     // Cut the active edges for chunking.
	     (AS_CGB_INTERCHUNK_EDGE == mate_edge_nes) ||
	     (AS_CGB_TOUCHES_CONTAINED_EDGE == mate_edge_nes) ||
	     
	     // Now cut other dovetail edges too.
	     (AS_CGB_BETWEEN_CONTAINED_EDGE == mate_edge_nes)
	     
	     // Someday we may want to cut containment edges too.
	     //(AS_CGB_TO_CONTAINED_EDGE == mate_edge_nes)
	     
	     ) {
	    const int mate_fragment_ovl = get_best_ovl(frags,edges,mate_edge);
	    
	    imate_used++;
	    assert(AS_CGB_BPT_IN_UNIQUE > AS_CGB_BPT_SLOP);
	    if(mate_fragment_ovl > seed_bpt_vertex + AS_CGB_BPT_IN_UNIQUE) {
	      // This is an anchor edge.
	      count_of_unique_region_overlaps ++;
	    } else if (mate_fragment_ovl <= seed_bpt_vertex + AS_CGB_BPT_SLOP) {
	      // This is a trimmable edge if there is a anchor edge.
	      count_of_repeat_region_overlaps ++;
	    } else {
	      // This edge is in the uncertain region.
	      count_of_unknown_region_overlaps ++;
	    }
	    
	  }
	}
	
#ifdef DEBUGGING
	if( (count_of_unique_region_overlaps == 1 ) &&
	    (count_of_repeat_region_overlaps  > 0 ) &&
	    (count_of_unknown_region_overlaps > 0 ) 
	    ) {
	  fprintf(stderr,"A foiled edge trimming! \n"
		  "seed_fragment_index=" F_IID ", seed_fragment_suffix=%d,"
		  "count_of_unique_region_overlaps=%d,"
		  "count_of_repeat_region_overlaps=%d,"
		  "count_of_unknown_region_overlaps=%d\n",
		  seed_fragment_index,seed_fragment_suffix,
		  count_of_unique_region_overlaps,
		  count_of_repeat_region_overlaps,
		  count_of_unknown_region_overlaps);
	}
#endif // DEBUGGING
	
#ifdef DEBUG_WITH_HISTO
	add_to_histogram(histo_unique_region_overlaps,
			 count_of_unique_region_overlaps,NULL);
	add_to_histogram(histo_repeat_region_overlaps,
			 count_of_repeat_region_overlaps,NULL);
	add_to_histogram(histo_unknown_region_overlaps,
			 count_of_unknown_region_overlaps,NULL);
#endif /*DEBUG_WITH_HISTO*/
	
	// At this point we have a fragment with overlaps.
	if( (count_of_unique_region_overlaps > 0) &&
	    (count_of_repeat_region_overlaps > 0) ) 
	  for(imate_index=0, imate_used=0;
	      (imate_index < seed_fragment_edge_degree); 
	      imate_index++) {
	    const IntEdge_ID mate_edge = seed_fragment_edge_list + imate_index;
	    const Tnes mate_edge_nes = get_nes_edge(edges,mate_edge);
	    // The adjacency list is already sorted by overhang when
	    // the original fragment adjacency list, that is from
	    // smallest to largest. We are taking just the MAX_NUM_OF_MATES
	    // 
	    if(
               // Cut the active edges for chunking.
               (AS_CGB_INTERCHUNK_EDGE == mate_edge_nes) ||
               (AS_CGB_TOUCHES_CONTAINED_EDGE == mate_edge_nes) ||
               
               // Now cut other dovetail edges too.
               (AS_CGB_BETWEEN_CONTAINED_EDGE == mate_edge_nes) 
               
               // Someday we may want to cut containment edges too.
               //(AS_CGB_TO_CONTAINED_EDGE == mate_edge_nes)

               ) {
	      const int mate_fragment_ovl = get_best_ovl(frags,edges,mate_edge);
	      imate_used++;
	      
	      assert(AS_CGB_BPT_IN_UNIQUE > AS_CGB_BPT_SLOP);
	      if(mate_fragment_ovl > seed_bpt_vertex + AS_CGB_BPT_IN_UNIQUE) {
		// Keep this edge.
	      } else if (mate_fragment_ovl <= seed_bpt_vertex + AS_CGB_BPT_SLOP) {
		// Trim this edge.
		count_of_removeable_overlaps ++;
                set_nes_edge(edges,mate_edge,AS_CGB_MARKED_BY_BRANCH_DVT);
                fix_overlap_edge_mate(frags,edges,mate_edge);
	      } else {
		// What to do with this edge?
	      }
	    }
	  }

#ifdef DEBUG_WITH_HISTO
	  add_to_histogram(histo_removeable_overlaps,
			   count_of_removeable_overlaps,NULL);
#endif // DEBUG_WITH_HISTO
	}
    }
#ifdef DEBUG_WITH_HISTO
  fprintf(fout,"\n\nHistogram of unique_region_overlaps\n");
  print_histogram(fout,histo_unique_region_overlaps, 0, 1);

  fprintf(fout,"\n\nHistogram of repeat_region_overlaps\n");
  print_histogram(fout,histo_repeat_region_overlaps, 0, 1);

  fprintf(fout,"\n\nHistogram of histo_unknown_region_overlaps\n");
  print_histogram(fout,histo_unknown_region_overlaps, 0, 1);

  fprintf(fout,"\n\nHistogram of histo_removeable_overlaps\n");
  print_histogram(fout,histo_removeable_overlaps, 0, 1);

  free_histogram(histo_unique_region_overlaps);
  free_histogram(histo_repeat_region_overlaps);
  free_histogram(histo_unknown_region_overlaps);
  free_histogram(histo_removeable_overlaps);
#endif /*DEBUG_WITH_HISTO*/
}
	      
