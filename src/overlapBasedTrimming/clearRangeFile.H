
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
 *    Brian P. Walenz on 2015-APR-15
 *      are Copyright 2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2018-JAN-17
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef CLEAR_RANGE_FILE_H
#define CLEAR_RANGE_FILE_H

#include "AS_global.H"

#include "sqStore.H"

#include "files.H"

class clearRangeFile {
public:

  //  Create a clear range file.  If the file doesn't exist, maxID must be set to the number
  //  of reads allowed.  If the file does exist, we don't care and load the number of reads
  //  from the file.

  clearRangeFile(char *fileName, sqStore *seq) {

    strncpy(_fileName, fileName, FILENAME_MAX-1);

    if (AS_UTL_fileExists(_fileName) == false) {
      _modified = true;   //  Always gets written, regardless of changes

      _lastID   = seq->sqStore_getNumReads();
      _bgn      = new uint32 [_lastID + 1];
      _end      = new uint32 [_lastID + 1];

      //  Every range is undefined

      memset(_bgn, 0, sizeof(uint32) * (_lastID + 1));
      memset(_end, 0, sizeof(uint32) * (_lastID + 1));

      //  Nope, every range is set to fully clear.

      reset(seq);
    }

    else {
      _modified = false;

      FILE  *F = AS_UTL_openInputFile(_fileName);

      AS_UTL_safeRead(F, &_lastID, "clearRangeFile::lastID", sizeof(uint32), 1);

      assert(seq->sqStore_getNumReads() == _lastID);  //  Sane?

      _bgn = new uint32 [_lastID + 1];
      _end = new uint32 [_lastID + 1];

      AS_UTL_safeRead(F,  _bgn, "clearRangeFile::bgn", sizeof(uint32), _lastID + 1);
      AS_UTL_safeRead(F,  _end, "clearRangeFile::end", sizeof(uint32), _lastID + 1);

      AS_UTL_closeFile(F, _fileName);
    }
  };


  //  Close the clear range file.

  ~clearRangeFile() {
    if (_modified == true) {
      FILE  *F = AS_UTL_openOutputFile(_fileName);

      AS_UTL_safeWrite(F, &_lastID, "clearRangeFile::lastID", sizeof(uint32), 1);

      AS_UTL_safeWrite(F,  _bgn, "clearRangeFile::bgn", sizeof(uint32), _lastID + 1);
      AS_UTL_safeWrite(F,  _end, "clearRangeFile::end", sizeof(uint32), _lastID + 1);

      AS_UTL_closeFile(F, _fileName);
    }

    delete [] _bgn;
    delete [] _end;
  };


  void      reset(sqStore *seq) {
    for (uint32 fi=1; fi <= _lastID; fi++) {
      _bgn[fi] = 0;
      _end[fi] = seq->sqStore_getRead(fi)->sqRead_sequenceLength();
    }
  };


  uint32    bgn(uint32 id)     { assert(id <= _lastID);  return(_bgn[id]); };
  uint32    end(uint32 id)     { assert(id <= _lastID);  return(_end[id]); };


  uint32   &setbgn(uint32 id)  { assert(id <= _lastID);  _modified = true;  return(_bgn[id]); };
  uint32   &setend(uint32 id)  { assert(id <= _lastID);  _modified = true;  return(_end[id]); };


  bool      isDeleted(uint32 id) {
    assert(id <= _lastID);
    return((_bgn[id] == UINT32_MAX) && (_end[id] == UINT32_MAX));
  };


  void      setDeleted(uint32 id) {
    assert(id <= _lastID);
    _modified = true;
    _bgn[id]  = UINT32_MAX;
    _end[id]  = UINT32_MAX;
  };


  void      copy(clearRangeFile *source) {
    if (source == NULL)
      return;

    assert(_lastID   == source->_lastID);

    memcpy(_bgn, source->_bgn, sizeof(uint32) * (_lastID + 1));
    memcpy(_end, source->_end, sizeof(uint32) * (_lastID + 1));
  };


private:
  bool     _modified;

  char     _fileName[FILENAME_MAX];

  uint32   _lastID;  //  [_lastID] is valid; allocated _lastID+1 spots.
  uint32  *_bgn;
  uint32  *_end;
};


#endif // CLEAR_RANGE_FILE_H
