
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

#include "runtime.H"

//  Returns a list of read IDs sorted by the number of reads in a BOG path
//  seeded by that read.

class ChunkGraph {
public:
  ChunkGraph(const char *prefix);
  ~ChunkGraph(void);

  uint32 nextReadByChunkLength(void) {
    if (_chunkLength[_chunkLengthIter].pathLen == 0)   //  By construction, the array
      return(0);                                       //  always ends with pathLen == 0.

    return(_chunkLength[_chunkLengthIter++].readId);
  };

private:
  uint32 countFullWidth(ReadEnd firstEnd, uint32 *pathLen, FILE *chunkLog);

  struct ChunkLength {
    uint32 readId;
    uint32 pathLen;
  };

  ChunkLength        *_chunkLength;
  uint32              _chunkLengthIter;
};

extern ChunkGraph       *CG;

#endif  //  INCLUDE_AS_BAT_CHUNKGRAPH
