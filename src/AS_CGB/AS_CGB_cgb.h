
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
/*********************************************************************
 * $Id: AS_CGB_cgb.h,v 1.4 2005-03-22 19:48:24 jason_miller Exp $
 *
 * Module: AS_CGB_cgb.h
 *
 * Description: The Chunk Graph Builder header file.
 *
 * The basic chunk data structures are:
 * 
 *   (1) a segmented array of type "TChunkFrag" called chunkfrags[] used
 *   to store the fragments of each chunk,
 *
 *   (2) a segmented array of type "TChunkOverlap" called chunkedges[]
 *   used to store the overlaps between chunks,
 *
 *   (3) an array of type "char" called chunksrc[] used to store
 *   annotation strings, and
 *
 *   (4) an array of type "TChunkMesg" thechunks[].
 *
 *   The first three arrays are segmented to store variable length
 *   information for the chunks.  Each member of the last array has a
 *   segment header index into each of the three previous arrays. 
 *
 *
 * Assumptions:
 * Author: Clark Mobarry
 *******************************************************************/

#ifndef AS_CGB_CGB_INCLUDE
#define AS_CGB_CGB_INCLUDE


/* ******************************************************************** */

static float compute_coverage_statistic
(
 BPTYPE rho,
 int    number_of_randomly_sampled_fragments_in_chunk,
 float  global_fragment_arrival_rate
 )
{
  // Rho is the number of bases in the chunk between the first
  // fragment arrival and the last fragment arrival.  It is the sum of
  // the fragment overhangs in the chunk. For intuitive purposes you
  // can think of it as the length of the chunk minus the length of
  // the last fragment. Thus a singleton chunk has a rho equal to
  // zero.
  
  // A singleton chunk provides no information as to its local
  // fragment arrival rate. We need at least two closely spaced
  // fragments that are randomly sampled from the chunk to get a local
  // estimate of the fragment arrival rate.
  
  // The local arrival rate of fragments in the chunk is:
  // {arrival_rate_local =
  // ((float)(nfrag_randomly_sampled_in_chunk-1))/(float)rho}.

  // The arrival distance of fragments in the chunk is the reciprocal
  // of the last formula:
  // {arrival_distance_local =
  // ((float)rho)/((float)(nfrag_randomly_sampled_in_chunk-1)).

  // Note a problem with this formula is that a singleton chunk has a
  // coverage discriminator statistic of 0/0.
  
  // The formula for the coverage discriminator statistic for the
  // chunk is:
  // {(arrival_rate_global/arrival_rate_local - ln(2))
  // *(nfrag_randomly_sampled_in_chunk-1)}.
  // The division by zero singularity cancels out to give the formula:
  // {(arrival_rate_global*rho -
  // ln(2)*(nfrag_randomly_sampled_in_chunk-1)}.
  
  // Call fragments that are not randomly sampled in the genome as
  // "guide" fragments.  The modification for guides recognizes that
  // guides are not randomly spaced fragments in the genome.
  // Therefore they should not contribute to the fragment count in the
  // local arrival rate discriminator statistic.
  
  // to avoid recomputing ...
  const float ln2 = 0.693147f; // Logarithm base 2 of "e".
  const float sqrt2 = 1.414213f;

  const float coverage_statistic =
    ( global_fragment_arrival_rate > 0.f
      ? rho*global_fragment_arrival_rate
      - ln2*(number_of_randomly_sampled_fragments_in_chunk - 1)
      : 0.f );

  // Be careful of falsely creating a negative number while using
  // unsigned integer arthimetic.  The behavior is not what we want
  // here.
  
  // The coverage discriminator statistic should be positive for
  // single coverage, negative for multiple coverage, and near zero
  // for indecisive.

#undef ADJUST_FOR_PARTIAL_EXCESS
  // the standard statistic gives log likelihood ratio of expected depth vs.
  // twice expected depth; but when enough fragments are present, we can actually
  // test whether depth exceeds expected even fractionally; in deeply sequenced datasets
  // (e.g. bacterial genomes), this has been observed for repetitive segments
  if(rho>0&&global_fragment_arrival_rate>0.f){
    float lambda = global_fragment_arrival_rate * rho;
    float zscore = ((number_of_randomly_sampled_fragments_in_chunk -1)-lambda)/
      sqrt(lambda);
    float p = .5 - erf(zscore/sqrt2)*.5;
    if(coverage_statistic>5 && p < .001){
      fprintf(stderr,"Standard unitigger a-stat is %f , but only %e chance of this great an excess of fragments: obs = %d, expect = %g rho = " F_S64 " Will%s reset a-stat to 1.5\n",
	      coverage_statistic,p,
	      number_of_randomly_sampled_fragments_in_chunk-1,
	      lambda,rho,
#ifdef ADJUST_FOR_PARTIAL_EXCESS
	      ""
#else
	      " not"
#endif
	      );
#ifdef ADJUST_FOR_PARTIAL_EXCESS
      return 1.5;
#endif
    }
  }

  return coverage_statistic;
}


static int count_the_randomly_sampled_fragments_in_a_chunk
(/* Input Only*/
 const Tfragment   frags[],
 const TChunkFrag  chunkfrags[],
 const TChunkMesg  thechunks[],
 const IntChunk_ID chunk_index
 )
{ 
  IntFragment_ID num_of_randomly_sampled_fragments_in_the_chunk = 0;
  const AChunkMesg * const mychunk = GetAChunkMesg(thechunks,chunk_index);
  const int num_frags = mychunk->num_frags;
  IntFragment_ID ii = 0;
  
  for(ii=0;ii<num_frags;ii++) {
    const IntFragment_ID ifrag = 
      GetVA_AChunkFrag(chunkfrags,mychunk->f_list+ii)->vid;
    const FragType type = 
      get_typ_fragment(frags,ifrag);
    if (type == AS_READ || type == AS_EBAC || type == AS_EXTR){
      num_of_randomly_sampled_fragments_in_the_chunk ++;
      // Only AS_READ, AS_EBAC, and AS_EXTR fragments are to be used in Gene
      // Myers coverage discriminator A-statistic.
    }
  }
  return num_of_randomly_sampled_fragments_in_the_chunk;
}

/* ******************************************************************** */

float compute_the_global_fragment_arrival_rate
( 
 FILE *fout,
 /* Input Only */
 const BPTYPE nbase_in_genome,
/* Input/Output */
 const Tfragment     frags[],     /* The internal representation of
			       the fragment reads. I have one
			       extra for the segstart field. */
 const Tedge         edges[],     /* The internal representation of the
			       overlaps. */
 const float         estimated_global_fragment_arrival_rate,
 /* Output Only */
 const TChunkFrag    *chunkfrags,
 const TChunkMesg    *thechunks
 );

void chunk_graph_build_1
(
 /* Input Only */
 const char * const Graph_Store_File_Prefix,
 const int work_limit_placing_contained_fragments,		  
 const int walk_depth,
 const IntFragment_ID max_frag_iid,
 const BPTYPE nbase_in_genome,
 const char * chimeras_file,
 const char * spurs_file,
 const float cgb_unique_cutoff,
 // Don't count chimeras
 /* Input/Output */
 Tfragment frags[],     /* The internal representation of
			   the fragment reads. */
 Tedge     edges[],     /* The internal representation of the
			 overlaps. */
 /* Output Only */
 float         *global_fragment_arrival_rate,
 // This rate is the number of fragments/bp, which in the design space
 // is between zero and one.
 TChunkFrag    *chunkfrags,
 TChunkMesg    *thechunks
 );

void chunk_graph_build_2
(
 /* Input Only */
 const int        use_consensus, 
 // Use the consensus sequence for the chunk to compute the branch points.
 const int        dont_find_branch_points, 
 // Don^t compute branch points at all.
 const float   cgb_unique_cutoff,
 const float   global_fragment_arrival_rate,
 // This rate is the number of fragments/bp, which in the design space
 // is between zero and one.
 FragStoreHandle TheFragStore,
 /* Input/Output */
 Tfragment frags[],     /* The internal representation of
			   the fragment reads. */
 Tedge     edges[],     /* The internal representation of the
			 overlaps. */
 TChunkFrag    *chunkfrags,
 TChunkMesg    *thechunks,
 /* Output Only */
 VA_TYPE(char) *chunkseqs,
 VA_TYPE(char) *chunkquas,
 VA_TYPE(char) *chunksrcs,
 FILE          *fbpts1,
 FILE          *fbpts2
 );

#endif /*AS_CGB_CGB_INCLUDE*/
