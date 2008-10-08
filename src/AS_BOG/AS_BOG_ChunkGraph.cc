
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

static const char *rcsid = "$Id: AS_BOG_ChunkGraph.cc,v 1.26 2008-10-08 22:02:54 brianwalenz Exp $";

#include "AS_BOG_ChunkGraph.hh"
#include "AS_BOG_BestOverlapGraph.hh"

#include <algorithm>

ChunkGraph::ChunkGraph(FragmentInfo *fi, BestOverlapGraph *BOG){

  _max_fragments   = fi->numFragments();

  iuid *pathLen = new iuid[(_max_fragments+1)*2];
  memset(pathLen, 0, sizeof(iuid)*(_max_fragments+1)*2);

  _chunk_lengths   = new _chunk_length      [_max_fragments];
  for(iuid frag_id=1; frag_id<=_max_fragments; frag_id++){
    _chunk_lengths[frag_id-1].fragId = frag_id;
    _chunk_lengths[frag_id-1].cnt    = (countFullWidth(BOG, pathLen, frag_id, FIVE_PRIME) +
                                        countFullWidth(BOG, pathLen, frag_id, THREE_PRIME));
  }

  std::sort(_chunk_lengths, _chunk_lengths + _max_fragments);

  delete[] pathLen;
}


uint32
ChunkGraph::countFullWidth(BestOverlapGraph *BOG, iuid *pathLen, iuid firstFrag, fragment_end_type end) {
  iuid index = firstFrag * 2 + end;

  if (pathLen[index] > 0)
    return pathLen[index];

  uint32                cnt    = 0;
  uint32                cntMax = 0;
  std::set<FragmentEnd> seen;
  FragmentEnd           lastEnd(firstFrag, end);

  //  Until we run off the chain, or we hit a fragment with a known
  //  length, compute the length FROM THE START.
  //
  while ((lastEnd.fragId() != NULL_FRAG_ID) &&
         (pathLen[index] == 0)) {

    seen.insert(lastEnd);

    pathLen[index] = ++cnt;

    //  Follow the path of lastEnd

    iuid frag = ((lastEnd.fragEnd() == FIVE_PRIME) ?
                 BOG->getBestEdgeOverlap(lastEnd.fragId(), FIVE_PRIME) ->frag_b_id :
                 BOG->getBestEdgeOverlap(lastEnd.fragId(), THREE_PRIME)->frag_b_id);

    if (frag == NULL_FRAG_ID) {
      lastEnd = FragmentEnd();
    } else {
      BOG->followOverlap( &lastEnd );
      assert(frag == lastEnd.fragId());
    }

    index = lastEnd.fragId() * 2 + lastEnd.fragEnd();
  }

  //  Check why we stopped.  If we've seen this fragment before, we
  //  just encountered a cycle in the graph.  Adjust every node in the
  //  cycle to have the correct length.
  //
  //  Otherwise, we've stopped because we've hit a fragment with a
  //  known path length.

  if (seen.find(lastEnd) != seen.end()) {
    //  If we end because of a cycle, mark points in circle same count.
    //  lastEnd will be the first fragment in the path that is in the
    //  circle, the * below (the path follows the top edge around).
    //
    //  >---------*------\
    //            |      |
    //            \------/

    //  pathLen[index] is the length of the path at the *.
    //
    iuid        circleLen = cnt - pathLen[index];
    FragmentEnd currEnd   = lastEnd;

    do {
      pathLen[index] = circleLen;
      BOG->followOverlap(&currEnd);

      index = currEnd.fragId() * 2 + currEnd.fragEnd();
    } while (lastEnd != currEnd);

  } else if (lastEnd.fragId() != NULL_FRAG_ID) {
    //  Otherwise, if lastEnd is not NULL, we must have stopped
    //  because of an existing path.
    cnt += pathLen[index];
  }

  //  Our return value is now whatever count we're at.
  cntMax = cnt;

  //  Traverse the last again, converting path length from the start
  //  to be path length from the end.  Notice that we stop at either
  //  the end of the path -- either NULL or the first frag with known
  //  length -- or the start of the circle.
  //
  FragmentEnd currEnd(firstFrag, end);

  while (currEnd.fragId() != lastEnd.fragId()) {
    pathLen[currEnd.fragId() * 2 + currEnd.fragEnd()] = cnt--;
    BOG->followOverlap( &currEnd );
  }

  return(cntMax);
}
