
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

#include "clearRangeFile.H"
#include "files.H"

clearRangeFile::clearRangeFile(char const *filename) {

  //  Clear the filename, because it isn't set if the filename passed in is NULL.

  memset(_filename, 0, sizeof(char) * (FILENAME_MAX+1));

  //  Save the name of the file we want to use.

  setFilename(filename);

  //  Allocate initial space for data.

  _modified = true;

  _minSet   = UINT32_MAX;
  _maxSet   = 0;

  _maxAlloc = 0;

  _bgn      = NULL;
  _end      = NULL;
  _flags    = NULL;

#if 0
  _maxAlloc = 128 * 1024;

  _bgn      = new uint32 [_maxAlloc];
  _end      = new uint32 [_maxAlloc];
  _flags    = new uint8  [_maxAlloc];

  memset(_bgn,   0, sizeof(uint32) * _maxAlloc);
  memset(_end,   0, sizeof(uint32) * _maxAlloc);
  memset(_flags, 0, sizeof(uint8)  * _maxAlloc);
#endif
}



clearRangeFile::~clearRangeFile() {

  saveData();

  delete [] _bgn;
  delete [] _end;
  delete [] _flags;
}



void
clearRangeFile::setFilename(char const *filename) {

  if (filename == NULL)   //  Do nothing if no filename supplied; use whatever is there already.
    return;

  memset(_filename, 0, sizeof(char) * (FILENAME_MAX+1));
  strncpy(_filename, filename, FILENAME_MAX);
}



void
clearRangeFile::loadData(char const *filename) {

  setFilename(filename);

  FILE  *F = AS_UTL_openInputFile(_filename);

  loadFromFile(_minSet, "clearRangeFile::minSet",         F);
  loadFromFile(_maxSet, "clearRangeFile::maxSet",         F);

  reallocData(_maxSet, true);

  loadFromFile(_bgn   + _minSet,  "clearRangeFile::bgn",   _maxSet - _minSet + 1, F);
  loadFromFile(_end   + _minSet,  "clearRangeFile::end",   _maxSet - _minSet + 1, F);
  loadFromFile(_flags + _minSet,  "clearRangeFile::flags", _maxSet - _minSet + 1, F);

  AS_UTL_closeFile(F, _filename);

  fprintf(stderr, "clearRangeFile()-- loaded reads %u-%u from '%s'.\n", _minSet, _maxSet, _filename);

  _modified = false;
};



void
clearRangeFile::saveData(char const *filename, bool force) {

  setFilename(filename);

  if ((_modified == false) && (force == false))
    return;

  FILE  *F = AS_UTL_openOutputFile(_filename);

  writeToFile(_minSet, "clearRangeFile::minSet",         F);
  writeToFile(_maxSet, "clearRangeFile::maxSet",         F);

  writeToFile(_bgn   + _minSet,  "clearRangeFile::bgn",   _maxSet - _minSet + 1, F);
  writeToFile(_end   + _minSet,  "clearRangeFile::end",   _maxSet - _minSet + 1, F);
  writeToFile(_flags + _minSet,  "clearRangeFile::flags", _maxSet - _minSet + 1, F);

  AS_UTL_closeFile(F, _filename);

  fprintf(stderr, "clearRangeFile()-- saved reads %u-%u to '%s'.\n", _minSet, _maxSet, _filename);

  _modified = false;
};



void
clearRangeFile::reallocData(uint64 newmax, bool exact) {
  uint64  unused = 0;

  if (newmax < _maxAlloc)
    return;

  if (exact == false) {
    newmax *= 3;
    newmax /= 2;
  }

  if (newmax < 128 * 1024)
    newmax = 128 * 1024;

  setArraySize(_bgn,   _maxAlloc, unused, newmax, resizeArray_copyData | resizeArray_clearNew);
  setArraySize(_end,   _maxAlloc, unused, newmax, resizeArray_copyData | resizeArray_clearNew);
  setArraySize(_flags, _maxAlloc, unused, newmax, resizeArray_copyData | resizeArray_clearNew);

  _maxAlloc = newmax;
}
