
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

#include "AS_BAT_Logging.H"

#include "AS_BAT_Unitig.H"
#include "AS_BAT_TigVector.H"



TigVector::TigVector(uint32 nReads) {

  //  The read-to-tig map

  _inUnitig  = new uint32 [nReads + 1];
  _ufpathIdx = new uint32 [nReads + 1];

  for (uint32 ii=0; ii<nReads+1; ii++) {
    _inUnitig[ii]   = 0;
    _ufpathIdx[ii]  = UINT32_MAX;
  }

  //  The vector

  _blockSize    = 1048576;

  _numBlocks    = 1;
  _maxBlocks    = 1024;
  _blocks       = new Unitig ** [_maxBlocks];
  _blocks[0]    = new Unitig  * [_blockSize];
  memset(_blocks[0], 0, sizeof(Unitig *) * _blockSize);
  _blockNext    = 1;

  _totalTigs    = 1;
};



TigVector::~TigVector() {

  //  Delete the maps

  delete [] _inUnitig;
  delete [] _ufpathIdx;

  //  Delete the tigs.

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
TigVector::newUnitig(bool verbose) {
  Unitig *u = new Unitig(this);

#pragma omp critical
  {
    u->_id = _totalTigs++;

    if (verbose)
      writeLog("Creating Unitig %d\n", u->_id);

    if (_blockNext >= _blockSize) {
      assert(_numBlocks < _maxBlocks);

      _blocks[_numBlocks] = new Unitig * [_blockSize];

      memset(_blocks[_numBlocks], 0, sizeof(Unitig *) * _blockSize);

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
TigVector::deleteUnitig(uint32 i) {
  delete _blocks[i / _blockSize][i % _blockSize];
  _blocks[i / _blockSize][i % _blockSize] = NULL;
}



#ifdef CHECK_UNITIG_ARRAY_INDEXING
Unitig *&operator[](uint32 i) {
  uint32  idx = i / _blockSize;
  uint32  pos = i % _blockSize;

  if (((i    >= _totalTigs)) ||
      ((idx  >= _numBlocks)) ||
      (((pos >= _blockNext) && (idx >= _numBlocks - 1)))) {
    writeStatus("TigVector::operator[]()--  i=" F_U32 " with totalTigs=" F_U64 "\n", i, _totalTigs);
    writeStatus("TigVector::operator[]()--  blockSize=" F_U64 "\n", _blockSize);
    writeStatus("TigVector::operator[]()--  idx=" F_U32 " numBlocks=" F_U64 "\n", idx, _numBlocks);
    writeStatus("TigVector::operator[]()--  pos=" F_U32 " blockNext=" F_U64 "\n", pos, _blockNext);
  }
  assert(i    < _totalTigs);
  assert((idx < _numBlocks));
  assert((pos < _blockNext) || (idx < _numBlocks - 1));

  return(_blocks[idx][pos]);
};
#endif



void
TigVector::computeArrivalRate(const char *prefix, const char *label) {
  uint32  tiLimit = size();
  uint32  numThreads = omp_get_max_threads();
  uint32  blockSize = (tiLimit < 100000 * numThreads) ? numThreads : tiLimit / 99999;

  writeStatus("computeArrivalRate()-- Computing arrival rates for %u tigs, with %u thread%s.\n", tiLimit, numThreads, (numThreads == 1) ? "" : "s");

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

    snprintf(N, FILENAME_MAX, "%s.arrivalRate.%u.dat", prefix, ii);
    FILE *F = AS_UTL_openOutputFile(N);
    for (uint32 jj=0; jj<hist[ii].size(); jj++)
      fprintf(F, "%d\n", hist[ii][jj]);
    AS_UTL_closeFile(F, N);
  }
}






void
TigVector::computeErrorProfiles(const char *prefix, const char *label) {
  uint32  tiLimit = size();
  uint32  numThreads = omp_get_max_threads();
  uint32  blockSize = (tiLimit < 100000 * numThreads) ? numThreads : tiLimit / 99999;

  writeStatus("computeErrorProfiles()-- Computing error profiles for %u tigs, with %u thread%s.\n", tiLimit, numThreads, (numThreads == 1) ? "" : "s");

#pragma omp parallel for schedule(dynamic, blockSize)
  for (uint32 ti=0; ti<tiLimit; ti++) {
    Unitig  *tig = operator[](ti);

    if (tig == NULL)
      continue;

    if (tig->ufpath.size() == 1)
      continue;

    tig->computeErrorProfile(prefix, label);
  }

  writeStatus("computeErrorProfiles()-- Finished.\n");
}



void
TigVector::reportErrorProfiles(const char *prefix, const char *label) {
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

