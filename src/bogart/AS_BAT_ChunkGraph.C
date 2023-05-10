
/******************************************************************************
 *
 *  This file is part of canu, a software program that assembles whole-genome
 *  sequencing reads into contigs.
 *
 *  This software is based on:
 *    'Celera Assembler' r4587 (http://wgs-assembler.sourceforge.net)
 *    the 'kmer package' r1994 (http://kmer.sourceforge.net)
 *
 *  Except as indicated otherwise, this is a 'United States Government Work',
 *  and is released in the public domain.
 *
 *  File 'README.licenses' in the root directory of this distribution
 *  contains full conditions and disclaimers.
 */

#include "AS_BAT_ReadInfo.H"
#include "AS_BAT_BestOverlapGraph.H"
#include "AS_BAT_ChunkGraph.H"

#include "AS_BAT_Logging.H"

#include <set>
#include <algorithm>



ChunkGraph::ChunkGraph(const char *prefix) {
  uint32   maxID    = RI->numReads();

  //  Allocate space for storing the full path length from each read.  There
  //  is no zeroth read, and this path length will always be zero.  When
  //  _chunkLength is sorted, this zeroth read will be at the end of the
  //  list, which is used to stop the iteration in nextReadByChunkLength().
  //
  _chunkLength     = new ChunkLength [maxID + 1];
  _chunkLengthIter = 0;

  for (uint32 fid=0; fid <= maxID; fid++) {
    _chunkLength[fid].readId  = fid;
    _chunkLength[fid].pathLen = 0;
  }

  //  Temporary storage of the path length from each read end.  Initialize to
  //  bogus value so we can use '0' as a valid path length.

  uint32  *endPathLen = new uint32 [maxID * 2 + 2];

  endPathLen[0] = 0;   //  But the zeroth read must have a zero length path.
  endPathLen[1] = 0;   //  It _shouldn't_ be used, but not verified.

  for (uint32 epl=2; epl < maxID * 2 + 2; epl++)
    endPathLen[epl] = 0;

  //  For each actual read, compute both end path lengths, and save the
  //  total path length in _chunkLength.

  FILE *chunkLog = (logFileFlagSet(LOG_CHUNK_GRAPH)) ? AS_UTL_openOutputFile(prefix, '.', "chunkGraph.log") : NULL;

  for (uint32 fid=1; fid <= maxID; fid++) {
    if ((RI->isValid(fid)       == false) ||     //  Read just doesn't exist.
        (OG->isContained(fid)   == true))        //  Read is contained, not in a path.
      continue;

    if (OG->isCoverageGap(fid)  == true)         //  Read is chimeric.  Explicitly skip, otherwise
      continue;                                  //  assert in countFillWidth() fails.

    _chunkLength[fid].pathLen = (countFullWidth(ReadEnd(fid, false), endPathLen, chunkLog) +
                                 countFullWidth(ReadEnd(fid, true),  endPathLen, chunkLog));
  }

  AS_UTL_closeFile(chunkLog, prefix, '.', "chunkGraph.log");

  delete [] endPathLen;

  //  Sort by decreasing path length.

  auto decreasingReadCount = [](ChunkLength const &a, ChunkLength const &b) {
                               return((a.pathLen > b.pathLen) || ((a.pathLen == b.pathLen) && (a.readId < b.readId)));
                             };

  std::sort(_chunkLength, _chunkLength + maxID + 1, decreasingReadCount);
}



ChunkGraph::~ChunkGraph(void) {
  delete [] _chunkLength;
};



//  Return the ReadEnd we'd get by following the edge out of the supplied
//  ReadEnd.
//
//  If there is no edge, a ReadEnd with readId == 0 is returned.
//
ReadEnd
followOverlap(ReadEnd end) {
  BestEdgeOverlap *edge = OG->getBestEdgeOverlap(end);

  return(ReadEnd(edge->readId(), !edge->read3p()));
}


uint64
getIndex(ReadEnd e) {
  return(e.readId() * 2 + e.read3p());
}


uint32
ChunkGraph::countFullWidth(ReadEnd firstEnd, uint32 *endPathLen, FILE *chunkLog) {
  uint64   firstIdx = getIndex(firstEnd);

  assert(firstEnd.readId() != 0);

  if (endPathLen[firstIdx] != 0) {
    if (chunkLog)
      fprintf(chunkLog, "path from %d,%d'(length=%u)\n",
              firstEnd.readId(),
              (firstEnd.read3p()) ? 3 : 5,
              endPathLen[firstIdx]);
    return(endPathLen[firstIdx]);
  }

  uint32                length = 0;
  std::set<ReadEnd>     seen;
  ReadEnd               lastEnd = firstEnd;
  uint64                lastIdx = firstIdx;

  //  Until we run off the chain, or we hit an end with a known length,
  //  compute the length of the path.

  while ((lastEnd.readId() != 0) &&
         (endPathLen[lastIdx] == 0)) {
    seen.insert(lastEnd);

    //  Should never get to a covergeGap read in a path
    //  (but we can get to lopsided).
    assert(OG->isCoverageGap(lastEnd.readId()) == false);

    endPathLen[lastIdx] = ++length;

    lastEnd = followOverlap(lastEnd);
    lastIdx = getIndex(lastEnd);
  }

  //  Check why we stopped.  Three cases:
  //
  //  Case 1.  We just ran out of overlaps.  Do nothing.
  if (lastEnd.readId() == 0) {
    ;
  }

  //  Case 2.  We encountered a read with a path length already set, and this
  //  read is not in our path.
  else if (seen.find(lastEnd) == seen.end()) {
    length += endPathLen[lastIdx];
  }

  //  Case 3.  Same as 2, but this read IS in our path, thus, we're in a loop.
  //  Set the length of all reads in the loop to be the size of the loop.
  //
  //  'lastEnd' and 'lastIdx' are the read we've seen before, aka, the first
  //  read in the loop.  So, starting from the first read in the loop, set
  //  the path length to the cycle length until we get back to that same
  //  first read.
  else {
    uint32      cycleLen = length - endPathLen[lastIdx] + 1;
    ReadEnd     currEnd  = lastEnd;
    uint64      currIdx  = lastIdx;

    do {
      endPathLen[currIdx] = cycleLen;

      currEnd = followOverlap(currEnd);
      currIdx = getIndex(currEnd);
    } while (lastEnd != currEnd);
  }

  //  Traverse again, converting "path length from the start" into "path
  //  length from the end".
  //
  //  For case 3 above, the reads in the loop have had their path lengths set
  //  already, and 'lastEnd' is the first read in the loop.

  ReadEnd     currEnd = firstEnd;
  uint64      currIdx = firstIdx;

  while (currEnd != lastEnd) {
    endPathLen[currIdx] = length--;

    assert(endPathLen[currIdx] > 0);

    currEnd = followOverlap(currEnd);
    currIdx = getIndex(currEnd);
  }

  //  Now just do some logging.

  if (logFileFlagSet(LOG_CHUNK_GRAPH)) {
    seen.clear();

    currEnd = firstEnd;
    currIdx = firstIdx;

    if (chunkLog)
      fprintf(chunkLog, "path from %d,%d'(length=%u):",
              firstEnd.readId(),
              (firstEnd.read3p()) ? 3 : 5,
              endPathLen[firstIdx]);

    while ((currEnd.readId() != 0) &&
           (seen.find(currEnd) == seen.end())) {
      seen.insert(currEnd);

      if ((chunkLog) && (currEnd == lastEnd))
        fprintf(chunkLog, " LAST");

      if (chunkLog)
        fprintf(chunkLog, " %d,%d'(%u)",
                currEnd.readId(),
                (currEnd.read3p()) ? 3 : 5,
                endPathLen[currIdx]);

      currEnd = followOverlap(currEnd);
      currIdx = getIndex(currEnd);
    }

    if ((chunkLog) && (seen.find(currEnd) != seen.end()))
      fprintf(chunkLog, " CYCLE %d,%d'(%u)",
              currEnd.readId(),
              (currEnd.read3p()) ? 3 : 5,
              endPathLen[currIdx]);

    if (chunkLog)
      fprintf(chunkLog, "\n");
  }

  //  And return the path length from this read end.

  assert(endPathLen[firstIdx] > 0);

  return(endPathLen[firstIdx]);
}
