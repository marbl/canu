
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

#ifndef INCLUDE_AS_BAT_TIGVECTOR
#define INCLUDE_AS_BAT_TIGVECTOR

#include "runtime.H"

class Unitig;

class TigVector {
public:
  TigVector(uint32 nReads);
  ~TigVector();

  Unitig   *newUnitig(bool verbose=false);
  void      deleteUnitig(uint32 i);

  size_t    size(void)            {  return(_totalTigs);  };
  Unitig  *&operator[](uint32 i)  {  return(_blocks[i / _blockSize][i % _blockSize]);  };

  void      optimizePositions(const char *prefix, const char *label);

  void      computeArrivalRate(const char *prefix, const char *label);

  void      computeErrorProfiles(const char *prefix, const char *label);
  void      reportErrorProfiles(const char *prefix, const char *label);

  //  Mapping from read to position in a tig.
public:
  void      registerRead(uint32 readId, uint32 tigid=0, uint32 ufpathidx=UINT32_MAX) {
    _inUnitig[readId]  = tigid;
    _ufpathIdx[readId] = ufpathidx;
  };

  uint32    inUnitig(uint32 readId)         {  return(_inUnitig[readId]);   };
  uint32    ufpathIdx(uint32 readId)        {  return(_ufpathIdx[readId]);  };

private:
  uint32    *_inUnitig;      //  Maps a read iid to a unitig id.
  uint32    *_ufpathIdx;     //  Maps a read iid to an index in ufpath

  //  The actual vector.
private:
  uint64     _blockSize;

  uint64     _numBlocks;
  uint64     _maxBlocks;
  Unitig  ***_blocks;
  uint64     _blockNext;

  uint64     _totalTigs;
};


#endif  //  INCLUDE_AS_BAT_TIGVECTOR
