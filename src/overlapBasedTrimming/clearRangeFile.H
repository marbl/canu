
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

#include "runtime.H"

#include "sqStore.H"

#include "files.H"

//  Create a clear range file.  If the file doesn't exist, maxID must be
//  set to the number of reads allowed.  If the file does exist, we don't
//  care and load the number of reads from the file.
//
//  If no input file, set up for a new clear range:
//    it's modified (so it'll get written out)
//    every range is (first cleared) then set to the whole read
//
//  But if there is a file, load the clear ranges saved.


class clearRangeFile {
public:
  clearRangeFile(char const *fileName, sqStore *seq) {

    _modified = false;

    memset(_fileName, 0, sizeof(char) * (FILENAME_MAX+1));
    strncpy(_fileName, fileName, FILENAME_MAX);

    _lastID   = seq->sqStore_lastReadID();
    _bgn      = new uint32 [_lastID + 1];
    _end      = new uint32 [_lastID + 1];

    memset(_bgn, 0, sizeof(uint32) * (_lastID + 1));
    memset(_end, 0, sizeof(uint32) * (_lastID + 1));


    if (fileExists(_fileName) == false) {
      _modified = true;

     reset(seq);
    }

    else {
      _modified = false;

      FILE  *F = AS_UTL_openInputFile(_fileName);

      loadFromFile(_lastID, "clearRangeFile::lastID",             F);
      loadFromFile(_bgn,    "clearRangeFile::bgn",   _lastID + 1, F);
      loadFromFile(_end,    "clearRangeFile::end",   _lastID + 1, F);

      AS_UTL_closeFile(F, _fileName);

      assert(seq->sqStore_lastReadID() == _lastID);  //  Sane?  And if wrong, we've overwritten _bgn and _end, and crash anyway.
    }
  };


  ~clearRangeFile() {
    if (_modified == true) {
      FILE  *F = AS_UTL_openOutputFile(_fileName);

      writeToFile(_lastID, "clearRangeFile::lastID",             F);
      writeToFile(_bgn,    "clearRangeFile::bgn",   _lastID + 1, F);
      writeToFile(_end,    "clearRangeFile::end",   _lastID + 1, F);

      AS_UTL_closeFile(F, _fileName);
    }

    delete [] _bgn;
    delete [] _end;
  };


  void      reset(sqStore *seq) {
    for (uint32 fi=1; fi <= _lastID; fi++) {
      _bgn[fi] = 0;
      _end[fi] = seq->sqStore_getReadLength(fi);
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

  char     _fileName[FILENAME_MAX+1];

  uint32   _lastID;  //  [_lastID] is valid; allocated _lastID+1 spots.
  uint32  *_bgn;
  uint32  *_end;
};


#endif // CLEAR_RANGE_FILE_H
