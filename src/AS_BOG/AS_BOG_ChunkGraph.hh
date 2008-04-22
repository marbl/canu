
/**************************************************************************
 * This file is part of Celera Assembler, a software program that 
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, The Venter Institute. All rights reserved.
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

#ifndef INCLUDE_AS_BOG_CHUNKGRAPH
#define INCLUDE_AS_BOG_CHUNKGRAPH

#include "AS_BOG_Datatypes.hh"
#include "AS_BOG_BestOverlapGraph.hh"

namespace AS_BOG{

  struct ChunkGraph{

  public:

    // Number of frags edges to cross in 
    static const short FRAG_WALK_NUM = 101;

    ChunkGraph(void);
    ~ChunkGraph(void);

    // Build the ChunkGraph, based on a BOG
    void build(BestOverlapGraph *bovlg);		

    // Chunkability rule
    bool isChunkable( BestEdgeOverlap *beo );

    virtual bool isChunkable( iuid frag_a_id, fragment_end_type which_end);

    // Returns IUID of 5' or 3' end of specified frag_id
    // Since there should only be one out/incoming connection
    iuid getChunking(iuid src_frag_id,
                     fragment_end_type whichEnd);

    void getChunking(iuid src_frag_id, 
                     iuid& five_prime_dst_frag_id, iuid& three_prime_dst_frag_id);

    void setChunking(iuid src_frag_id, 
                     iuid five_prime_dst_frag_id, iuid three_prime_dst_frag_id);

    long getNumFragments(void);
    long countSingletons(void);

    void checkInDegree();

    iuid nextFragByChunkLength();

    // follows the graph path to the next frag end
    FragmentEnd followPath(FragmentEnd);

  protected:
    BestOverlapGraph *bovlg;

  private:

    struct _chunk_unit_struct{
      iuid five_prime;
      iuid three_prime;
    };
    struct _chunk_length {
      iuid fragId;
      short cnt;
    };

    short countChunkWidth(iuid, fragment_end_type );
    iuid countFullWidth(iuid, fragment_end_type );

    static int sortChunkLens(const void*,const void*);

    _chunk_unit_struct *_chunkable_array;
    iuid *_edgePathLen;
    _chunk_length *_chunk_lengths;
    iuid _max_fragments;

  };

  struct PromiscuousChunkGraph : public ChunkGraph {
    bool isChunkable( iuid frag_id, fragment_end_type which_end);
  };

} //AS_BOG namespace


#endif

