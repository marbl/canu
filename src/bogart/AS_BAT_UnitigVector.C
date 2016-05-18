
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

#include "AS_BAT_Logging.H"

#include "AS_BAT_Unitig.H"
#include "AS_BAT_UnitigVector.H"



UnitigVector::UnitigVector() {
  _blockSize    = 1048576;
  _numBlocks    = 1;
  _maxBlocks    = 1024;
  _blocks       = new Unitig ** [_maxBlocks];
  _blocks[0]    = new Unitig  * [_blockSize];
  _blocks[0][0] = NULL;  //  No first unitig.
  _blockNext    = 1;
  _totalUnitigs = 1;
};



UnitigVector::~UnitigVector() {

  //  Delete the unitigs.
  for (uint32 ii=0; ii<_numBlocks; ii++)
    for (uint32 jj=0; jj<_blockSize; jj++)
      delete _blocks[ii][jj];

  //  Delete the blocks.
  for (uint32 ii=0; ii<_numBlocks; ii++)
    delete [] _blocks[ii];

  //  And the block pointers.
  delete [] _blocks;
};



Unitig *
UnitigVector::newUnitig(bool verbose) {
  Unitig *u = new Unitig();

#pragma omp critical
  {
    u->_id = _totalUnitigs++;

    if (verbose)
      writeLog("Creating Unitig %d\n", u->_id);

    if (_blockNext >= _blockSize) {
      assert(_numBlocks < _maxBlocks);

      _blocks[_numBlocks] = new Unitig * [_blockSize];

      memset(_blocks[_numBlocks], 0, sizeof(Unitig **) * _blockSize);

      _numBlocks++;
      _blockNext = 0;
    }

    _blocks[_numBlocks-1][_blockNext++] = u;

    //  The rest are just sanity checks.

    assert((u->id() / _blockSize) == (_numBlocks - 1));
    assert((u->id() % _blockSize) == (_blockNext - 1));

    assert(operator[](u->id()) == u);
  }

  return(u);
};



void
UnitigVector::deleteUnitig(uint32 i) {
  delete _blocks[i / _blockSize][i % _blockSize];
  _blocks[i / _blockSize][i % _blockSize] = NULL;
}



#ifdef CHECK_UNITIG_ARRAY_INDEXING
Unitig *&operator[](uint32 i) {
  uint32  idx = i / _blockSize;
  uint32  pos = i % _blockSize;

  if (((i    >= _totalUnitigs)) ||
      ((idx  >= _numBlocks)) ||
      (((pos >= _blockNext) && (idx >= _numBlocks - 1)))) {
    fprintf(stderr, "UnitigVector::operator[]()--  i="F_U32" with totalUnitigs="F_U64"\n", i, _totalUnitigs);
    fprintf(stderr, "UnitigVector::operator[]()--  blockSize="F_U64"\n", _blockSize);
    fprintf(stderr, "UnitigVector::operator[]()--  idx="F_U32" numBlocks="F_U64"\n", idx, _numBlocks);
    fprintf(stderr, "UnitigVector::operator[]()--  pos="F_U32" blockNext="F_U64"\n", pos, _blockNext);
  }
  assert(i    < _totalUnitigs);
  assert((idx < _numBlocks));
  assert((pos < _blockNext) || (idx < _numBlocks - 1));

  return(_blocks[idx][pos]);
};
#endif







void
UnitigVector::computeArrivalRate(const char *prefix, const char *label) {
  uint32  tiLimit = size();
  uint32  numThreads = omp_get_max_threads();
  uint32  blockSize = (tiLimit < 100000 * numThreads) ? numThreads : tiLimit / 99999;

  fprintf(stderr, "Computing arrival rates for %u unitigs using %u threads.\n", tiLimit, numThreads);

  vector<int32>  hist[6];

  //#pragma omp parallel for schedule(dynamic, blockSize)
  for (uint32 ti=0; ti<tiLimit; ti++) {
    Unitig  *tig = operator[](ti);

    if (tig == NULL)
      continue;

    if (tig->ufpath.size() == 1)
      continue;

    tig->computeArrivalRate(prefix, label, hist);
  }

  for (uint32 ii=1; ii<6; ii++) {
    char  N[FILENAME_MAX];

    sprintf(N, "%s.arrivalRate.%u.dat", prefix, ii);
    FILE *F = fopen(N, "w");
    for (uint32 jj=0; jj<hist[ii].size(); jj++)
      fprintf(F, "%d\n", hist[ii][jj]);
    fclose(F);
  }
}






void
UnitigVector::computeErrorProfiles(const char *prefix, const char *label) {
  uint32  tiLimit = size();
  uint32  numThreads = omp_get_max_threads();
  uint32  blockSize = (tiLimit < 100000 * numThreads) ? numThreads : tiLimit / 99999;

  fprintf(stderr, "Computing error profiles for %u unitigs using %u threads.\n", tiLimit, numThreads);

  //#pragma omp parallel for schedule(dynamic, blockSize)
  for (uint32 ti=0; ti<tiLimit; ti++) {
    Unitig  *tig = operator[](ti);

    if (tig == NULL)
      continue;

    if (tig->ufpath.size() == 1)
      continue;

    tig->computeErrorProfile(prefix, label);
  }

  fprintf(stderr, "Computing error profiles - FINISHED.\n");
}



void
UnitigVector::reportErrorProfiles(const char *prefix, const char *label) {
  uint32  tiLimit = size();
  uint32  numThreads = omp_get_max_threads();
  uint32  blockSize = (tiLimit < 100000 * numThreads) ? numThreads : tiLimit / 99999;

  for (uint32 ti=0; ti<tiLimit; ti++) {
    Unitig  *tig = operator[](ti);

    if (tig == NULL)
      continue;

    if (tig->ufpath.size() == 1)
      continue;

    tig->reportErrorProfile(prefix, label);
  }
}

