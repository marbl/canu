
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
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#ifndef CLEAR_RANGE_FILE_H
#define CLEAR_RANGE_FILE_H

#include "AS_global.H"

#include "gkStore.H"

#include "AS_UTL_fileIO.H"

class clearRangeFile {
public:

  //  Create a clear range file.  If the file doesn't exist, maxID must be set to the number
  //  of reads allowed.  If the file does exist, we don't care and load the number of reads
  //  from the file.

  clearRangeFile(char *fileName, gkStore *gkp) {

    strcpy(_fileName, fileName);

    if (AS_UTL_fileExists(_fileName) == false) {
      _modified = true;   //  Always gets written, regardless of changes

      _lastID   = gkp->gkStore_getNumReads();
      _bgn      = new uint32 [_lastID + 1];
      _end      = new uint32 [_lastID + 1];

      //  Every range is undefined

      memset(_bgn, 0, sizeof(uint32) * (_lastID + 1));
      memset(_end, 0, sizeof(uint32) * (_lastID + 1));

      //  Nope, every range is set to fully clear.

      reset(gkp);
    }

    else {
      _modified = false;

      FILE  *F = fopen(_fileName, "r");
      if (errno)
        fprintf(stderr, "clearRangeFile()--  Failed to open '%s' for loading clear ranges: %s\n", _fileName, strerror(errno)), exit(1);

      AS_UTL_safeRead(F, &_lastID, "clearRangeFile::lastID", sizeof(uint32), 1);

      assert(gkp->gkStore_getNumReads() == _lastID);  //  Sane?

      _bgn = new uint32 [_lastID + 1];
      _end = new uint32 [_lastID + 1];

      AS_UTL_safeRead(F,  _bgn, "clearRangeFile::bgn", sizeof(uint32), _lastID + 1);
      AS_UTL_safeRead(F,  _end, "clearRangeFile::end", sizeof(uint32), _lastID + 1);

      fclose(F);
    }
  };


  //  Close the clear range file.

  ~clearRangeFile() {
    if (_modified == true) {
      FILE  *F = fopen(_fileName, "w");
      if (errno)
        fprintf(stderr, "clearRangeFile()--  Failed to open '%s' for saving clear ranges: %s\n", _fileName, strerror(errno)), exit(1);

      AS_UTL_safeWrite(F, &_lastID, "clearRangeFile::lastID", sizeof(uint32), 1);

      AS_UTL_safeWrite(F,  _bgn, "clearRangeFile::bgn", sizeof(uint32), _lastID + 1);
      AS_UTL_safeWrite(F,  _end, "clearRangeFile::end", sizeof(uint32), _lastID + 1);

      fclose(F);
    }

    delete [] _bgn;
    delete [] _end;
  };


  void      reset(gkStore *gkp) {
    for (uint32 fi=1; fi <= _lastID; fi++) {
      _bgn[fi] = 0;
      _end[fi] = gkp->gkStore_getRead(fi)->gkRead_sequenceLength();
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
