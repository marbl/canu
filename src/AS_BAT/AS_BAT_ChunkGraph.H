
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

#ifndef INCLUDE_AS_BAT_CHUNKGRAPH
#define INCLUDE_AS_BAT_CHUNKGRAPH

static const char *rcsid_INCLUDE_AS_BAT_CHUNKGRAPH = "$Id: AS_BAT_ChunkGraph.H,v 1.2 2012-01-05 16:29:26 brianwalenz Exp $";

#include "AS_BAT_Datatypes.H"

class BestOverlapGraph;

class ChunkLength {
public:
  uint32 fragId;
  uint32 cnt;

  bool operator<(ChunkLength const that) const {
    if (cnt == that.cnt)
      return(fragId < that.fragId);
    return(cnt > that.cnt);
  };
};


class ChunkGraph {
public:
  ChunkGraph(const char *output_prefix);
  ChunkGraph(set<AS_IID> *restrict);
  ~ChunkGraph(void) {
    delete [] _chunkLength;
  };

  uint32 nextFragByChunkLength(void) {
    if (_chunkLengthIter >= _maxFragment)
      return(0);
    return(_chunkLength[_chunkLengthIter++].fragId);
  };

private:
  uint64 getIndex(FragmentEnd e);
  uint32 countFullWidth(FragmentEnd firstEnd);

  uint64              _maxFragment;

  //  The usual case, for a chunk graph of all fragments.
  ChunkLength        *_chunkLength;
  uint32              _chunkLengthIter;
  uint32             *_pathLen;

  //  For a chunk graph of a single unitig plus some extra fragments.
  //  This maps the AS_IID to an index in the arrays above.
  map<AS_IID,uint32>  _idMap;
  set<AS_IID>        *_restrict;
};

#endif
