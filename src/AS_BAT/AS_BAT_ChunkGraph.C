
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

static const char *rcsid = "$Id: AS_BAT_ChunkGraph.C,v 1.2 2012-01-05 16:29:26 brianwalenz Exp $";

#include "AS_BAT_Datatypes.H"
#include "AS_BAT_ChunkGraph.H"
#include "AS_BAT_BestOverlapGraph.H"


ChunkGraph::ChunkGraph(const char *output_prefix) {

  setLogFile(output_prefix, "ChunkGraph");

  _maxFragment = FI->numFragments();
  _restrict    = NULL;

  _pathLen         = new uint32      [_maxFragment * 2 + 2];
  _chunkLength     = new ChunkLength [_maxFragment];
  _chunkLengthIter = 0;

  memset(_pathLen,     0, sizeof(uint32)      * (_maxFragment * 2 + 2));
  memset(_chunkLength, 0, sizeof(ChunkLength) * (_maxFragment));

  for (uint32 fid=1; fid <= _maxFragment; fid++) {
    if (OG->isContained(fid))
      continue;

    if (OG->isSuspicious(fid))
      //  Fragment is suspicious.  We won't seed a BOG from it, and populateUnitig will make only a
      //  singleton.
      continue;

    _chunkLength[fid-1].fragId = fid;
    _chunkLength[fid-1].cnt    = (countFullWidth(FragmentEnd(fid, false)) +
                                  countFullWidth(FragmentEnd(fid, true)));
  }

  delete [] _pathLen;
  _pathLen = NULL;

  std::sort(_chunkLength, _chunkLength + _maxFragment);
}



ChunkGraph::ChunkGraph(set<AS_IID> *restrict) {

  _maxFragment = 0;
  _restrict    = restrict;

  for (set<AS_IID>::iterator it=_restrict->begin(); it != _restrict->end(); it++)
    _idMap[*it] = _maxFragment++;

  _pathLen         = new uint32      [_maxFragment * 2 + 2];
  _chunkLength     = new ChunkLength [_maxFragment];
  _chunkLengthIter = 0;

  memset(_pathLen,     0, sizeof(uint32)      * (_maxFragment * 2 + 2));
  memset(_chunkLength, 0, sizeof(ChunkLength) * (_maxFragment));

  for (set<AS_IID>::iterator it=_restrict->begin(); it != _restrict->end(); it++) {
    uint32  fid = *it;          //  Actual fragment ID
    uint32  fit = _idMap[fid];  //  Local array index

    if (OG->isContained(fid))
      continue;

    _chunkLength[fit].fragId = fid;
    _chunkLength[fit].cnt    = (countFullWidth(FragmentEnd(fid, false)) +
                                countFullWidth(FragmentEnd(fid, true)));
  }

  delete [] _pathLen;
  _pathLen = NULL;

  std::sort(_chunkLength, _chunkLength + _maxFragment);
}




uint64
ChunkGraph::getIndex(FragmentEnd e) {
  if (_restrict == NULL)
    return(e.fragId() * 2 + e.frag3p());

  return(_idMap[e.fragId()] * 2 + e.frag3p());
}


uint32
ChunkGraph::countFullWidth(FragmentEnd firstEnd) {
  uint64   firstIdx = getIndex(firstEnd);

  assert(firstIdx < _maxFragment * 2 + 2);

  if (_pathLen[firstIdx] > 0)
    return _pathLen[firstIdx];

  uint32                length = 0;
  std::set<FragmentEnd> seen;
  FragmentEnd           lastEnd = firstEnd;
  uint64                lastIdx = firstIdx;

  //  Until we run off the chain, or we hit a fragment with a known length, compute the length FROM
  //  THE START.
  //
  while ((lastIdx != 0) &&
         (_pathLen[lastIdx] == 0)) {

    seen.insert(lastEnd);

    _pathLen[lastIdx] = ++length;

    //  Follow the path of lastEnd

    lastEnd = OG->followOverlap(lastEnd);
    lastIdx = getIndex(lastEnd);
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
    uint32      cycleLen = length - _pathLen[lastIdx] + 1;
    FragmentEnd currEnd  = lastEnd;
    uint64      currIdx  = lastIdx;

    do {
      _pathLen[currIdx] = cycleLen;
      currEnd = OG->followOverlap(currEnd);
      currIdx = getIndex(currEnd);
    } while (lastEnd != currEnd);

  } else {
    //  Case 2, an existing path.
    length += _pathLen[lastIdx];
  }

  //  Our return value is now whatever count we're at.
  uint32 lengthMax = length;

  //  Traverse again, converting "path length from the start" into "path length from the end".  Any
  //  cycle has had its length set correctly already, and we stop at either the start of the cycle,
  //  or at the start of any existing path.
  //
  FragmentEnd currEnd = firstEnd;
  uint64      currIdx = firstIdx;

  while (currEnd != lastEnd) {
    _pathLen[currIdx] = length--;
    currEnd = OG->followOverlap(currEnd);
    currIdx = getIndex(currEnd);
  }

  if (logFileFlagSet(LOG_CHUNK_GRAPH)) {
    seen.clear();

    currEnd = firstEnd;
    currIdx = firstIdx;

    fprintf(logFile, "PATH from %d,%d length %d:",
            firstEnd.fragId(),
            (firstEnd.frag3p()) ? 3 : 5,
            _pathLen[firstIdx]);

    while ((currEnd.fragId() != 0) &&
           (seen.find(currEnd) == seen.end())) {
      seen.insert(currEnd);

      if (currEnd == lastEnd)
        fprintf(logFile, " LAST");

      fprintf(logFile, " %d,%d(%d)",
              currEnd.fragId(),
              (currEnd.frag3p()) ? 3 : 5,
              _pathLen[currIdx]);

      currEnd = OG->followOverlap(currEnd);
      currIdx = getIndex(currEnd);
    }

    if (seen.find(currEnd) != seen.end())
      fprintf(logFile, " CYCLE %d,%d(%d)",
              currEnd.fragId(),
              (currEnd.frag3p()) ? 3 : 5,
              _pathLen[currIdx]);

    fprintf(logFile, "\n");
  }

  if (lengthMax != _pathLen[firstIdx])
    fprintf(logFile, "ERROR: lengthMax %d _pathLen[] %d\n",
            lengthMax, _pathLen[firstIdx]);
  assert(lengthMax == _pathLen[firstIdx]);

  return(_pathLen[firstIdx]);
}
