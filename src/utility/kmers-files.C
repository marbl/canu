
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
 *    Brian P. Walenz beginning on 2018-JUL-21
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"

#include "kmers.H"
#include "kmers-files.H"
#include "kmers-histogram.H"



merylFileBlockReader::merylFileBlockReader() {
  _data        = NULL;

  _blockPrefix = 0;
  _nKmers      = 0;
  _nKmersMax   = 0;

  _kCode       = 0;
  _unaryBits   = 0;
  _binaryBits  = 0;
  _k1          = 0;

  _cCode       = 0;
  _c1          = 0;
  _c2          = 0;

  _suffixes    = NULL;
  _values      = NULL;
}


merylFileBlockReader::~merylFileBlockReader() {
  delete    _data;
  delete [] _suffixes;
  delete [] _values;
}


bool
merylFileBlockReader::loadBlock(FILE *inFile, uint32 activeFile, uint32 activeIteration) {

  //  If _data exists, we've already loaded the block, but haven't used it yet.

  if (_data)
    return(true);

  //  Otherwise, allocate _data, read the block from disk.  If nothing loaded,
  //  return false.

  _data = new stuffedBits();

  _blockPrefix = 0;
  _nKmers      = 0;

  if (_data->loadFromFile(inFile) == false) {
    delete _data;
    _data = NULL;

    return(false);
  }

  //  Decode the header of _data, but don't process the kmers yet.

  uint64 pos   = _data->getPosition();
  uint64 m1    = _data->getBinary(64);
  uint64 m2    = _data->getBinary(64);

  _blockPrefix = _data->getBinary(64);
  _nKmers      = _data->getBinary(64);

  _kCode       = _data->getBinary(8);
  _unaryBits   = _data->getBinary(32);
  _binaryBits  = _data->getBinary(32);
  _k1          = _data->getBinary(64);

  _cCode       = _data->getBinary(8);
  _c1          = _data->getBinary(64);
  _c2          = _data->getBinary(64);

#ifdef SHOW_LOAD
  fprintf(stderr, "loadBlock()-- file %u iter %u:\n", activeFile, activeIteration);
  fprintf(stderr, "    prefix     0x%016lx\n", _blockPrefix);
  fprintf(stderr, "    nKmers     " F_U64 "\n", _nKmers);
  fprintf(stderr, "    kCode      " F_U32 "\n", _kCode);
  fprintf(stderr, "    unaryBits  " F_U32 "\n", _unaryBits);
  fprintf(stderr, "    binaryBits " F_U32 "\n", _binaryBits);
  fprintf(stderr, "    k1efix     " F_U64 "\n", _k1);
  fprintf(stderr, "    cCode      " F_U32 "\n", _cCode);
  fprintf(stderr, "    c1         " F_U64 "\n", _c1);
  fprintf(stderr, "    c2         " F_U64 "\n", _c2);
#endif

  if ((m1 != 0x7461446c7972656dllu) ||
      (m2 != 0x0a3030656c694661llu)) {
    fprintf(stderr, "merylFileReader::nextMer()-- Magic number mismatch in activeFile " F_U32 " activeIteration " F_U32 " position " F_U64 ".\n", activeFile, activeIteration, pos);
    fprintf(stderr, "merylFileReader::nextMer()-- Expected 0x7461446c7972656d got 0x%016" F_X64P "\n", m1);
    fprintf(stderr, "merylFileReader::nextMer()-- Expected 0x0a3030656c694661 got 0x%016" F_X64P "\n", m2);
    exit(1);
  }

  return(true);
}



void
merylFileBlockReader::decodeBlock(void) {
  if (_data == NULL)
    return;

  //fprintf(stderr, "decodeBlock() nKmersMax %lu nKmers %lu\n", _nKmersMax, _nKmers);

  resizeArrayPair(_suffixes, _values, 0, _nKmersMax, _nKmers, resizeArray_doNothing);

  decodeBlock(_suffixes, _values);
}



void
merylFileBlockReader::decodeBlock(kmdata *suffixes, kmvalu *values) {

  if (_data == NULL)
    return;

  kmdata  thisPrefix = 0;

  //  Decode the suffixes.

  if      (_kCode == 1) {
    for (uint32 kk=0; kk<_nKmers; kk++) {
      thisPrefix += (kmdata)_data->getUnary();

      uint32 ls = (_binaryBits <= 64) ? (0)           : (_binaryBits - 64);
      uint32 rs = (_binaryBits <= 64) ? (_binaryBits) : (64);

      suffixes[kk]   = thisPrefix;
      suffixes[kk] <<= ls;
      suffixes[kk]  |= _data->getBinary(ls);
      suffixes[kk] <<= rs;
      suffixes[kk]  |= _data->getBinary(rs);
    }
  }

  else {
    fprintf(stderr, "ERROR: unknown kCode %u\n", _kCode), exit(1);
  }

  //  Decode the values.

  if      (_cCode == 1) {
    for (uint32 kk=0; kk<_nKmers; kk++)
      values[kk] = _data->getBinary(32);
  }

  else if (_cCode == 2) {
    for (uint32 kk=0; kk<_nKmers; kk++)
      values[kk] = _data->getBinary(64);
  }

  else {
    fprintf(stderr, "ERROR: unknown cCode %u\n", _cCode), exit(1);
  }

  delete _data;
  _data = NULL;
}


