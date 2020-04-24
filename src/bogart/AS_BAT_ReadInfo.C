
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

#include "AS_BAT_ReadInfo.H"
#include "AS_BAT_Logging.H"



ReadInfo::ReadInfo(const char *seqStorePath,
                   const char *prefix,
                   uint32      minReadLen,
                   uint32      maxReadLen) {

  uint32    numNotPresent = 0;
  uint32    numShort      = 0;
  uint32    numLong       = 0;
  uint32    numLoaded     = 0;

  _seqStore     = new sqStore(seqStorePath);

  _numBases     = 0;
  _numReads     = _seqStore->sqStore_lastReadID();
  _numLibraries = _seqStore->sqStore_lastLibraryID();

  _readStatus    = new ReadStatus [_numReads + 1];

  //  Initialize to indicate a read that isn't valid.

  for (uint32 fi=0; fi <= _numReads; fi++) {
    _readStatus[fi].readLength = 0;
    _readStatus[fi].libraryID  = 0;
    _readStatus[fi].isPresent  = false;
    _readStatus[fi].isIgnored  = false;
    _readStatus[fi].unused     = 0;
  }

  //  Scan the store.
  //    Flag any read 'ignored' in the store as 'not present' in the assembly.
  //    Flag any read too short              as 'not present' in the assembly.
  //

  for (uint32 fi=1; fi <= _numReads; fi++) {
    uint32   len  = _seqStore->sqStore_getReadLength(fi);

    if (_seqStore->sqStore_isIgnoredRead(fi)) {
      numNotPresent++;
    }

    else if (len < minReadLen) {
      numShort++;
    }

    else if (maxReadLen < len) {
      numLong++;
    }

    else {
      _readStatus[fi].readLength = len;
      _readStatus[fi].libraryID  = _seqStore->sqStore_getLibraryIDForRead(fi);
      _readStatus[fi].isPresent  = true;
      _readStatus[fi].isIgnored  = false;

      _numBases += len;

      numLoaded++;
    }
  }


  writeStatus("ReadInfo()-- Found   %9u reads.\n", numLoaded);

  if (minReadLen > 0)
    writeStatus("ReadInfo()-- Ignored %9u reads shorter than %6u bp.\n", numShort, minReadLen);

  if (maxReadLen < UINT32_MAX)
    writeStatus("ReadInfo()-- Ignored %9u reads longer  than %6u bp.\n", numLong, maxReadLen);
}



ReadInfo::~ReadInfo() {
  delete    _seqStore;
  delete [] _readStatus;
}
