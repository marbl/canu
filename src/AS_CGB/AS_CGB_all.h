
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

//  $Id: AS_CGB_all.h,v 1.20 2007-07-20 17:17:08 brianwalenz Exp $

#ifndef AS_CGB_ALL_INCLUDE
#define AS_CGB_ALL_INCLUDE

#undef DONT_RUN_IN_SYMMETRIC_MODE      // Allow un-mated dovetail edges in CGB.  (??)

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>

#include <fcntl.h>
#include <assert.h>

#include "AS_global.h"
#include "AS_PER_gkpStore.h"
#include "AS_UTL_Var.h"
#include "AS_UTL_param_proc.h"

typedef uint32 IntEdge_ID;

#define AS_CGB_NO_MATE          0
#define FRAGMENT_NOT_VISITED   INT32_MAX
#define CHUNK_NOT_VISITED      INT32_MAX
#define AS_CGB_EDGE_NOT_FOUND  (~((IntEdge_ID)(0)))



// The transitive overlap inference slop is 
//
// abs((length of the middle fragment) +
//     (length of the opposite overlap) -
//     (length of the righthand overlap) -
//     (length of the lefthand overlap)))
//
//  <=  AS_CGB_TRANSITIVE_SLOP_ALPHA +
//      AS_CGB_TRANSITIVE_SLOP_EPSILON*(overlap in bp)
//
//  These values allow for 40 bp slop in an average 512 bp
//  containment overlap.
//
#define AS_CGB_TRANSITIVE_SLOP_ALPHA   20     // bp
#define AS_CGB_TRANSITIVE_SLOP_EPSILON 0.07  // fraction 


typedef enum {

  AS_CGB_UNLABELED_FRAG=0,
  
  /* Initial non-contained fragment types */
  AS_CGB_SOLO_FRAG=1, 
  // A non-contained fragment with no raw dovetail overlaps at
  // all. This could be because of crappy sequencing, hard repetitive
  // sequence screening on both fragment ends, or a physical
  // sequencing gap.
  AS_CGB_HANGING_FRAG=2, 
  // A non-contained fragment with no raw dovetail overlaps on one
  // side. This could be because of crappy sequencing, hard repetitive
  // sequence screening on one fragment end, or a physical sequencing
  // gap.
  AS_CGB_THRU_FRAG=3, 
  // A non-contained fragment with raw dovetail overlaps on both
  // fragment ends. 

  /* Modified non-contained fragment types */
  AS_CGB_HANGING_CRAPPY_FRAG=4,    
  // A formerly hanging fragment that needs to be treated like a
  // contained unitig so that it does not inhibit chunking.
  AS_CGB_HANGING_CHUNK_FRAG=5, 
  // A formerly hanging fragment that terminates a non-singleton light
  // chunk. This could be due to crappy sequencing or hard repetitive
  // sequence screening.
  AS_CGB_INTERCHUNK_FRAG=6, 
  // A formerly thru fragment that terminates a light chunk and
  // has essentail non-containment overlaps to more than one light
  // chunk.
  AS_CGB_INTRACHUNK_FRAG=7, 
  // A formerly thru fragment that is interior to a light chunk. 

  /* Initial contained fragment type */
  AS_CGB_ORPHANEDCONT_FRAG=8,
  // A contained fragment that needs to be singleton unitig because it
  // was not placed into a unitig by a containment edge.

  /* Modified contained fragment types */
  AS_CGB_MULTICONT_FRAG=9, 
  // A contained fragment that needs to be singleton unitig because it
  // could be placed in multiple unitigs by containment edges.
  AS_CGB_SINGLECONT_FRAG=10, 
  // A contained fragment that has essential containment overlaps to
  // only one light chunk.

  /* spur & chimera & other bad fragments */
  AS_CGB_MARKED_BREAKER_FRAG=11,
  // soft - marked
  AS_CGB_REMOVED_BREAKER_FRAG=12,
  // hard - removed
  
  AS_CGB_UNPLACEDCONT_FRAG=13,
  // A contained fragment that has not yet been placed into a chunk.
  AS_CGB_BRANCHMULTICONT_FRAG=14,
  // A contained fragment that even after transitive overlap removal
  // still has a dovetail overlap to a non-contained fragment. Many of
  // these contained fragments are in repetitive regions near a
  // branch-point.
  AS_CGB_ESSENTIAL_CONT_FRAG=15,
  
  /* Other fragment types */
  AS_CGB_DELETED_FRAG=127

} Tlab;


//           dovetail 
// \superset nonchordal
// \superset nontransitive
// \superset thickest (a directed edge concept)
// \superset interchunk (essential edges for chunking.)
// \superset intrachunk


typedef enum {
  AS_CGB_UNUSED_EDGE=0,

  AS_CGB_DOVETAIL_EDGE=1,
  
  AS_CGB_THICKEST_EDGE=2,
  // A dovetail overlap edge that is the thickest from the proximal
  // fragment-end

  AS_CGB_BETWEEN_CONTAINED_EDGE=4,
  // A dovetail overlap between globally contained fragments.

  AS_CGB_TOUCHES_CONTAINED_EDGE=6,
  // A dovetail overlap touching a globally contained fragment and a
  // non-contained fragment.

  /* Containment Overlap types (asymmetric edges) */
  AS_CGB_CONTAINED_EDGE=12,

  /* Dovetail Overlap types (symmetric edges) */
  AS_CGB_INTERCHUNK_EDGE=21,
  // A dovetail overlap exterior to a chunk.

  AS_CGB_INTRACHUNK_EDGE=23,
  // A dovetail overlap interior to a chunk.


  AS_CGB_TOUCHES_CRAPPY_DVT=32,
  // A dovetail overlap touching a crappy fragment and a non-crappy
  // fragment.
  AS_CGB_BETWEEN_CRAPPY_DVT=33,
  // A dovetail overlap between crappy fragments.
  AS_CGB_TOUCHES_CRAPPY_CON=38,
  AS_CGB_BETWEEN_CRAPPY_CON=39,


  AS_CGB_MARKED_BY_BRANCH_DVT=57,
  // An dovetail overlap removed by being on the repeat side of a
  // branch point.

  AS_CGB_MARKED_BY_BREAKER=58,
  // An edge to spur, chimera & other bad fragment
  
  AS_CGB_MARKED_BY_DELETED_DVT=61,
  AS_CGB_MARKED_BY_DELETED_CON=64,

  AS_CGB_REMOVED_BY_TRANSITIVITY_DVT=101,
  // An dovetail overlap that is inferrable by the dovetail chain overlap transivity
  // rules.
  AS_CGB_REMOVED_BY_TRANSITIVITY_CON=104,
  // An containment overlap that is inferrable by the overlap
  // transivity rules.

  AS_CGB_REMOVED_BY_THRESHOLD_DVT=111,
  // An dovetail overlap that is beyond the adjacency degree
  // threshold.
  AS_CGB_REMOVED_BY_THRESHOLD_CON=114,
  // An containment overlap that is beyond the adjacency degree
  // threshold.

  AS_CGB_REMOVED_BY_BREAKER=114,
  // An edge to spur, chimera & other bad fragment
  
  AS_CGB_REMOVED_BY_DUPLICATE_DVT=116,
  AS_CGB_REMOVED_BY_DUPLICATE_CON=119
  // An edge removed because it was a duplicate.

} Tnes;


#include "AS_CGB_unitigger_globals.h"
#include "AS_CGB_methods.h"


//  AS_CGB_util.c
Tnes AS_CGB_SafeCast_cdsInt8_to_Tnes (int8 dataIn) ;



//  AS_FGB_hanging_fragment.c
void separate_fragments_as_solo_hanging_thru(Tfragment frags[], Tedge edges[]);
void identify_early_spur_fragments(Tfragment frags[], Tedge edges[]);

//  AS_FGB_contained.c
void check_containment_edges(Tfragment frags[], Tedge edges[]);
void contained_fragment_marking_frc(Tfragment frags[], Tedge edges[]);

//  AS_FGB_fragmentHash.h
#define AS_CGB_NOT_SEEN_YET  INT32_MAX


//  AS_CGB_traversal.c
void as_graph_traversal(FILE *fout,
			Tfragment frags[],
			Tedge edges[], 
			size_t fragment_visited[]);

//  AS_CGB_fgb.c
void view_fgb_chkpnt(char *Store_Path_Prefix,
                     Tfragment frags[], 
                     Tedge edges[]);


//  AS_CGB_fgb.c
void reorder_edges(Tfragment *frags,
                   Tedge *edges);

//  AS_CGB_count_fragment_and_edge_labels.c
void count_fragment_and_edge_labels(Tfragment frags[],
                                    Tedge     edges[],
                                    char      comment[]);


////////////////////////////////////////
//  AS_CGB_cgb.c
int count_the_randomly_sampled_fragments_in_a_chunk (const Tfragment   frags[],
                                                     const TChunkFrag  chunkfrags[],
                                                     const TChunkMesg  thechunks[],
                                                     const IntChunk_ID chunk_index);


float compute_the_global_fragment_arrival_rate ( 
 const int           recalibrate, /* Boolean flag to recalibrate global arrival rate to max
				     unique local arrival rate */
 const float         cgb_unique_cutoff, /* threshold for unique chunks */
 FILE               *fout,
 /* Input Only */
 const int64         nbase_in_genome,
 /* Input/Output */
 const Tfragment     frags[],
 const Tedge         edges[],
 const float         estimated_global_fragment_arrival_rate,
 /* Output Only */
 const TChunkFrag   *chunkfrags,
 const TChunkMesg   *thechunks
 );

//  AS_CGB_cgb.c (end)
////////////////////////////////////////


////////////////////////////////////////
//  AS_CGB_edgemate.c
void reflect_Aedge( Aedge *new_edge, Aedge *old_edge);
void granger_Aedge( Aedge *new_edge, Aedge *old_edge);

void fix_overlap_edge_mate (const Tfragment frags[], 
                            Tedge edges[],
                            const IntEdge_ID ie0);

IntEdge_ID check_symmetry_of_the_edge_mates (Tfragment frags[],
                                             Tedge edges[]);

//  AS_CGB_edgemate.c (end)
////////////////////////////////////////


////////////////////////////////////////
//  AS_CGB_walk.c
#undef WALK_DEPTH_DIAGNOSTICS
#define MATCH_TARGET_EDGE

int is_there_an_overlap_path(const Tfragment * const frags,
                             const Tedge     * const edges,
                             const IntFragment_ID target_avx,
                             const int        target_asx,
#ifdef MATCH_TARGET_EDGE
                             const IntFragment_ID target_bvx,
                             const int        target_bsx,
                             const int        target_ahg,
                             const int        target_bhg,
                             const Tnes       target_nes,
                             const int        target_is_dovetail,
                             const int        target_is_from_contained,
                             const int        target_is_to_contained,
                             const int        target_is_dgn_contained,
#endif // MATCH_TARGET_EDGE
                             /* recursion variables: */
                             const int        search_depth, // Use zero at the top level.
                             const IntFragment_ID current_avx,  // Use target_avx at the top level.
                             const int        current_asx,  // Use target_asx at the top level.
                             const int        current_ahg,  // Use zero at the top level.
                             const int        current_bhg,  // Use zero at the top level.
                             const int        last_edge_was_contained, // From target_nes at the top level.
                             /* search path limiting: */
                             const int        walk_depth,   // The maximum depth of the stack.
                             const int        tolerance,    // For the overlaps
                             IntFragment_ID visited_a[],
#ifdef MATCH_TARGET_EDGE
                             IntFragment_ID visited_b[],
#endif // MATCH_TARGET_EDGE
                             // Was this fragment visited from the target overlap before? That is
                             // by the same (target_avx, target_asx, target_bvx, target_bsx).
                             const int        work_limit_per_candidate_edge,
                             // Maximum number of edges to explore per candidate edge.
                             int *            work_tally_per_candidate_edge,
                             // Current number of edges explored per candidate edge.
                             /* diagnostic variables: */
#ifdef WALK_DEPTH_DIAGNOSTICS
                             int64      search_depth_histogram[], // How many times this depth is visited.
                             int64      search_path_histogram[],  // How many paths were this length.
#endif // WALK_DEPTH_DIAGNOSTICS
                             int64      *ntrans_test_fail // How many paths failed the last test.
                             );
//  AS_CGB_walk.c (end)
////////////////////////////////////////







// Rho is the number of bases in the chunk between the first fragment
// arrival and the last fragment arrival.  It is the sum of the
// fragment overhangs in the chunk. For intuitive purposes you can
// think of it as the length of the chunk minus the length of the last
// fragment. Thus a singleton chunk has a rho equal to zero.
//
// A singleton chunk provides no information as to its local fragment
// arrival rate. We need at least two closely spaced fragments that
// are randomly sampled from the chunk to get a local estimate of the
// fragment arrival rate.
//
// The local arrival rate of fragments in the chunk is:
//
// arrival_rate_local = ((float)(nfrag_randomly_sampled_in_chunk-1))/(float)rho
//
// The arrival distance of fragments in the chunk is the reciprocal of
// the last formula:
//
// arrival_distance_local = ((float)rho)/((float)(nfrag_randomly_sampled_in_chunk-1))
//
// Note a problem with this formula is that a singleton chunk has a
// coverage discriminator statistic of 0/0.
//
// The formula for the coverage discriminator statistic for the chunk
// is:
//
// (arrival_rate_global/arrival_rate_local - ln(2))*(nfrag_randomly_sampled_in_chunk-1)
//
// The division by zero singularity cancels out to give the formula:
//
// (arrival_rate_global*rho - ln(2)*(nfrag_randomly_sampled_in_chunk-1)
//
// Call fragments that are not randomly sampled in the genome as
// "guide" fragments.  The modification for guides recognizes that
// guides are not randomly spaced fragments in the genome.  Therefore
// they should not contribute to the fragment count in the local
// arrival rate discriminator statistic.
//
// Be careful of falsely creating a negative number while using
// unsigned integer arthimetic.  The behavior is not what we want
// here.
//
// The coverage discriminator statistic should be positive for single
// coverage, negative for multiple coverage, and near zero for
// indecisive.
//
// ADJUST_FOR_PARTIAL_EXCESS: The standard statistic gives log
// likelihood ratio of expected depth vs.  twice expected depth; but
// when enough fragments are present, we can actually test whether
// depth exceeds expected even fractionally; in deeply sequenced
// datasets (e.g. bacterial genomes), this has been observed for
// repetitive segments.
//
static float compute_coverage_statistic(int64  rho,
                                        int    number_of_randomly_sampled_fragments_in_chunk,
                                        float  global_fragment_arrival_rate) {
  
 float ln2 = 0.693147f; // Logarithm base 2 of "e".
 float sqrt2 = 1.414213f;

 float coverage_statistic = 0.f;
 if (global_fragment_arrival_rate > 0.f)
   coverage_statistic = rho*global_fragment_arrival_rate - ln2*(number_of_randomly_sampled_fragments_in_chunk - 1);

#undef ADJUST_FOR_PARTIAL_EXCESS
#ifdef ADJUST_FOR_PARTIAL_EXCESS
  if(rho>0&&global_fragment_arrival_rate>0.f){
    float lambda = global_fragment_arrival_rate * rho;
    float zscore = ((number_of_randomly_sampled_fragments_in_chunk -1)-lambda) / sqrt(lambda);
    float p = .5 - erf(zscore/sqrt2)*.5;
    if(coverage_statistic>5 && p < .001){
      fprintf(stderr,"Standard unitigger a-stat is %f, but only %e chance of this great an excess of fragments: obs = %d, expect = %g rho = " F_S64 " Will reset a-stat to 1.5\n",
	      coverage_statistic,p,
	      number_of_randomly_sampled_fragments_in_chunk-1,
	      lambda,rho);
      return 1.5;
    }
  }
#endif

  return coverage_statistic;
}








extern int (*compare_edge_function)(const void *a, const void *b);


#endif /*AS_CGB_ALL_INCLUDE*/
