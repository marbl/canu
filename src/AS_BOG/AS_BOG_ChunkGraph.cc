
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

static const char *rcsid = "$Id: AS_BOG_ChunkGraph.cc,v 1.31 2010-09-25 07:48:45 brianwalenz Exp $";

#include "AS_BOG_Datatypes.hh"
#include "AS_BOG_ChunkGraph.hh"
#include "AS_BOG_BestOverlapGraph.hh"


ChunkGraph::ChunkGraph(FragmentInfo *fi, BestOverlapGraph *bog) {

  _BOG             = bog;

  _maxFragment     = fi->numFragments();

  _pathLen         = new uint32      [_maxFragment * 2 + 2];
  _chunkLength     = new ChunkLength [_maxFragment];
  _chunkLengthIter = 0;

  memset(_pathLen, 0, sizeof(uint32) * (_maxFragment * 2 + 2));

  for (uint32 fid=1; fid <= _maxFragment; fid++) {
    if (bog->isContained(fid))
      continue;

    _chunkLength[fid-1].fragId = fid;
    _chunkLength[fid-1].cnt    = (countFullWidth(FragmentEnd(fid, FIVE_PRIME)) +
                                  countFullWidth(FragmentEnd(fid, THREE_PRIME)));
  }

  delete [] _pathLen;
  _pathLen = NULL;

  std::sort(_chunkLength, _chunkLength + _maxFragment);
}


uint32
ChunkGraph::countFullWidth(FragmentEnd firstEnd) {

  assert(firstEnd.index() < _maxFragment * 2 + 2);

  if (_pathLen[firstEnd.index()] > 0)
    return _pathLen[firstEnd.index()];

  uint32                length = 0;
  std::set<FragmentEnd> seen;
  FragmentEnd           lastEnd = firstEnd;

  //  Until we run off the chain, or we hit a fragment with a known length, compute the length FROM
  //  THE START.
  //
  while ((lastEnd.fragId() != 0) &&
         (_pathLen[lastEnd.index()] == 0)) {

    seen.insert(lastEnd);

    _pathLen[lastEnd.index()] = ++length;

    //  Follow the path of lastEnd

    lastEnd = _BOG->followOverlap(lastEnd);
  }

  //  Check why we stopped.  Three cases:
  //
  //  1)  We ran out of best edges to follow -- lastEnd.fragId() == 0
  //  2)  We encountered a fragment with known length -- _pathLen[lastEnd.index()] > 0
  //  3)  We encountered a self-loop (same condition as case 2)
  //
  //  To distinguish case 2 and 3, we keep a set<> of the fragments we've seen in this construction.
  //  If 'lastEnd' is in that set, then we're case 3.  If so, adjust every node in the cycle to have
  //  the same length, the length of the cycle itself.
  //
  //  'lastEnd' and 'index' are the first fragment in the cycle; we've seen this one before.
  //
  if (lastEnd.fragId() == 0) {
    //  Case 1.  Do nothing.
    ;

  } else if (seen.find(lastEnd) != seen.end()) {
    //  Case 3, a cycle.
    uint32      cycleLen = length - _pathLen[lastEnd.index()] + 1;
    FragmentEnd currEnd  = lastEnd;
    do {
      _pathLen[currEnd.index()] = cycleLen;
      currEnd = _BOG->followOverlap(currEnd);
    } while (lastEnd != currEnd);

  } else {
    //  Case 2, an existing path.
    length += _pathLen[lastEnd.index()];
  }

  //  Our return value is now whatever count we're at.
  uint32 lengthMax = length;

  //  Traverse again, converting "path length from the start" into "path length from the end".  Any
  //  cycle has had its length set correctly already, and we stop at either the start of the cycle,
  //  or at the start of any existing path.
  //
  FragmentEnd currEnd = firstEnd;

  while (currEnd.fragId() != lastEnd.fragId()) {
    _pathLen[currEnd.index()] = length--;
    currEnd = _BOG->followOverlap(currEnd);
  }

#ifdef EMIT_CHUNK_GRAPH_PATH
  seen.clear();

  currEnd = firstEnd;

  fprintf(stderr, "PATH from %d,%d:",
          firstEnd.fragId(),
          (firstEnd.fragEnd() == FIVE_PRIME) ? 5 : 3);

  while ((currEnd.fragId() != 0) &&
         (seen.find(currEnd) == seen.end())) {
    seen.insert(currEnd);

    if (currEnd.fragId() == lastEnd.fragId())
      fprintf(stderr, " LAST");

    fprintf(stderr, " %d,%d(%d)",
            currEnd.fragId(),
            (currEnd.fragEnd() == FIVE_PRIME) ? 5 : 3,
            _pathLen[currEnd.index()]);

    currEnd = _BOG->followOverlap(currEnd);
  }

  if (seen.find(currEnd) != seen.end())
    fprintf(stderr, " CYCLE %d,%d(%d)",
            currEnd.fragId(),
            (currEnd.fragEnd() == FIVE_PRIME) ? 5 : 3,
            _pathLen[currEnd.index()]);

  fprintf(stderr, "\n");
#endif

  if (lengthMax != _pathLen[firstEnd.index()])
    fprintf(stderr, "ERROR: lengthMax %d _pathLen[] %d\n",
            lengthMax, _pathLen[firstEnd.index()]);
  assert(lengthMax == _pathLen[firstEnd.index()]);

  return(_pathLen[firstEnd.index()]);
}
