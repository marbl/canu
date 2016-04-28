
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
 *    src/AS_BAT/AS_BAT_ChunkGraph.H
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

#ifndef INCLUDE_AS_BAT_CHUNKGRAPH
#define INCLUDE_AS_BAT_CHUNKGRAPH

#include "AS_global.H"

#include <set>
#include <map>

using namespace std;


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
  ChunkGraph(set<uint32> *restrict);
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
  //  This maps the uint32 to an index in the arrays above.
  map<uint32,uint32>  _idMap;
  set<uint32>        *_restrict;
};



extern ChunkGraph       *CG;

#endif  //  INCLUDE_AS_BAT_CHUNKGRAPH
