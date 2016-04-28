
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' (http://kmer.sourceforge.net)
 *  both originally distributed by Applera Corporation under the GNU General
 *  Public License, version 2.
 *
 *  Canu branched from Celera Assembler at its revision 4587.
 *  Canu branched from the kmer project at its revision 1994.
 *
 *  This file is derived from:
 *
 *    src/AS_BAT/AS_BAT_ChunkGraph.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2010-NOV-23 to 2013-AUG-01
 *      are Copyright 2010,2012-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-DEC-19 to 2014-DEC-22
 *      are Copyright 2014 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2016-JAN-11
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_BAT_FragmentInfo.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_ChunkGraph.H"

#include "AS_BAT_Logging.H"


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
    if (OG->isContained(fid)) {
      if (logFileFlagSet(LOG_CHUNK_GRAPH))
        writeLog("read %u contained\n", fid);
      continue;
    }

    if (OG->isSuspicious(fid)) {
      if (logFileFlagSet(LOG_CHUNK_GRAPH))
        writeLog("read %u suspicious\n", fid);
      continue;
    }

    uint32  l5 = countFullWidth(FragmentEnd(fid, false));
    uint32  l3 = countFullWidth(FragmentEnd(fid, true));

    _chunkLength[fid-1].fragId = fid;
    _chunkLength[fid-1].cnt    = l5 + l3;
  }

  delete [] _pathLen;
  _pathLen = NULL;

  std::sort(_chunkLength, _chunkLength + _maxFragment);
}



ChunkGraph::ChunkGraph(set<uint32> *restrict) {

  _maxFragment = 0;
  _restrict    = restrict;

  for (set<uint32>::iterator it=_restrict->begin(); it != _restrict->end(); it++)
    _idMap[*it] = _maxFragment++;

  _pathLen         = new uint32      [_maxFragment * 2 + 2];
  _chunkLength     = new ChunkLength [_maxFragment];
  _chunkLengthIter = 0;

  memset(_pathLen,     0, sizeof(uint32)      * (_maxFragment * 2 + 2));
  memset(_chunkLength, 0, sizeof(ChunkLength) * (_maxFragment));

  for (set<uint32>::iterator it=_restrict->begin(); it != _restrict->end(); it++) {
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

    writeLog("path from %d,%d length %d:",
            firstEnd.fragId(),
            (firstEnd.frag3p()) ? 3 : 5,
            _pathLen[firstIdx]);

    while ((currEnd.fragId() != 0) &&
           (seen.find(currEnd) == seen.end())) {
      seen.insert(currEnd);

      if (currEnd == lastEnd)
        writeLog(" LAST");

      writeLog(" %d,%d(%d)",
              currEnd.fragId(),
              (currEnd.frag3p()) ? 3 : 5,
              _pathLen[currIdx]);

      currEnd = OG->followOverlap(currEnd);
      currIdx = getIndex(currEnd);
    }

    if (seen.find(currEnd) != seen.end())
      writeLog(" CYCLE %d,%d(%d)",
              currEnd.fragId(),
              (currEnd.frag3p()) ? 3 : 5,
              _pathLen[currIdx]);

    writeLog("\n");
  }

  if (lengthMax != _pathLen[firstIdx])
    writeLog("ERROR: lengthMax %d _pathLen[] %d\n",
            lengthMax, _pathLen[firstIdx]);
  assert(lengthMax == _pathLen[firstIdx]);

  return(_pathLen[firstIdx]);
}
