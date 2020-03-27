
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
