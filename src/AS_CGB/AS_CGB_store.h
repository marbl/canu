
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
 * $Id: AS_CGB_store.h,v 1.9 2007-05-14 09:27:10 brianwalenz Exp $
 *
 * Module: AS_CGB_store.h
 * Description: Header file for the code that reads and writes the 
 * check point data.
 * Assumptions:
 * Author: Clark Mobarry
 *********************************************************************/

#ifndef AS_CGB_STORE_INCLUDE
#define AS_CGB_STORE_INCLUDE

#define GENINFO

typedef struct {
  /* fragment name */
#ifdef GENINFO
  BPTYPE  genbgn,genend; /* Coordinates in genome */
#endif
} Afraginfo;


VA_DEF(Afraginfo);
VA_DEF(char);
VA_DEF(OverlapMesg);

typedef VA_TYPE(Afraginfo) Tfraginfo;

#ifdef GENINFO
static void set_genbgn_fraginfo(Tfraginfo * const fraginfo,IntFragment_ID i,BPTYPE value)
{ GetAfraginfo(fraginfo,i)->genbgn = (BPTYPE)value;}
static void set_genend_fraginfo(Tfraginfo * const fraginfo,IntFragment_ID i,BPTYPE value)
{ GetAfraginfo(fraginfo,i)->genend = (BPTYPE)value;}

static BPTYPE get_genbgn_fraginfo(const Tfraginfo * const fraginfo,IntFragment_ID i)
{ return (BPTYPE) GetAfraginfo(fraginfo,i)->genbgn;}
static BPTYPE get_genend_fraginfo(const Tfraginfo * const fraginfo,IntFragment_ID i)
{ return (BPTYPE) GetAfraginfo(fraginfo,i)->genend;}
#endif /*GENINFO*/


// End of fraginfo stuff... 


typedef struct { 
  int32       store_version;
  int32       state_of_the_store;
  int32       unused;

  unsigned int    edges_sorted_by_fragment_end : 1;
  unsigned int    unmated_edges : 1;
  unsigned int    sorted_by_adjacency_lists : 1;
  unsigned int    next_edge_array_valid : 1;
  unsigned int    dechorded : 1;
  unsigned int    transitively_marked : 1;
  unsigned int    using_to_contained_edges : 1;
  unsigned int    unusedbits : 25;

  BPTYPE          nbase_in_genome;
  IntFragment_ID  nfrag_randomly_sampled_in_genome;
  float           global_fragment_arrival_rate;
  // The estimated length of the genome in base pairs.
  IntFragment_ID  min_frag_iid;
  IntFragment_ID  max_frag_iid;

} TStateGlobals;


typedef struct { 

  Tfragment     *frags; /* The current fragment array handle. */
  Tedge         *edges;   /* The current edge array handle. */
  TIntEdge_ID   *next_edge_obj;
  Tfraginfo     *fraginfo;

  TChunkFrag    *chunkfrags;
  TChunkMesg    *thechunks;
  VA_TYPE(char) *chunkseqs;
  VA_TYPE(char) *chunkquas;

  VA_TYPE(char) *frag_annotations;
  VA_TYPE(char) *chunksrcs;

} THeapGlobals;


#endif /* AS_CGB_STORE_INCLUDE */
