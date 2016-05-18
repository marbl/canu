
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
 *  Modifications by:
 *
 *    Brian P. Walenz beginning on 2016-APR-06
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef INCLUDE_AS_BAT_UNITIG_VECTOR
#define INCLUDE_AS_BAT_UNITIG_VECTOR

#include "AS_global.H"

class Unitig;

class UnitigVector {
public:
  UnitigVector();
  ~UnitigVector();

  Unitig   *newUnitig(bool verbose);
  void      deleteUnitig(uint32 i);

  size_t    size(void)            {  return(_totalUnitigs);  };
  Unitig  *&operator[](uint32 i)  {  return(_blocks[i / _blockSize][i % _blockSize]);  };

  void      computeArrivalRate(const char *prefix, const char *label);

  void      computeErrorProfiles(const char *prefix, const char *label);
  void      reportErrorProfiles(const char *prefix, const char *label);

private:
  uint64              _blockSize;

  uint64              _numBlocks;
  uint64              _maxBlocks;
  Unitig           ***_blocks;
  uint64              _blockNext;

  uint64              _totalUnitigs;
};


#endif  //  INCLUDE_AS_BAT_UNITIG_VECTOR
