
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
 *    Brian P. Walenz beginning on 2016-AUG-12
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_BAT_ReadInfo.H"
#include "AS_BAT_Logging.H"



ReadInfo::ReadInfo(const char *seqStorePath,
                   const char *prefix,
                   uint32      minReadLen) {

  sqStore  *seqStore      = new sqStore(seqStorePath);
  uint32    numNotPresent = 0;
  uint32    numShort      = 0;
  uint32    numLoaded     = 0;

  _numBases     = 0;
  _numReads     = seqStore->sqStore_lastReadID();
  _numLibraries = seqStore->sqStore_lastLibraryID();

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
    uint32   len  = seqStore->sqStore_getReadLength(fi);

    if (seqStore->sqStore_isIgnoredRead(fi)) {
      numNotPresent++;
    }

    else if (len < minReadLen) {
      numShort++;
    }

    else {
      _readStatus[fi].readLength = len;
      _readStatus[fi].libraryID  = seqStore->sqStore_getLibraryIDForRead(fi);
      _readStatus[fi].isPresent  = true;
      _readStatus[fi].isIgnored  = false;

      _numBases += len;

      numLoaded++;
    }
  }

  delete seqStore;

  if (minReadLen > 0)
    writeStatus("ReadInfo()-- Using %d reads, ignoring %u reads less than " F_U32 " bp long.\n",
                numLoaded, numShort, minReadLen);

  else
    writeStatus("ReadInfo()-- Using %d reads, no minimum read length used.\n",
                numLoaded);
}



ReadInfo::~ReadInfo() {
  delete [] _readStatus;
}
