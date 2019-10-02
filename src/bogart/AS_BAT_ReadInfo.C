
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

  sqStore  *seqStore = new sqStore(seqStorePath);

  _numBases     = 0;
  _numReads     = seqStore->sqStore_lastReadID();
  _numLibraries = seqStore->sqStore_lastLibraryID();

  _readStatus    = new ReadStatus [_numReads + 1];

  for (uint32 i=0; i<_numReads + 1; i++) {
    _readStatus[i].readLength = 0;
    _readStatus[i].libraryID  = 0;
    _readStatus[i].isBackbone = false;
    _readStatus[i].isUnplaced = false;
    _readStatus[i].isLeftover = false;
    _readStatus[i].unused     = 0;
  }

  uint32 numSkipped = 0;
  uint32 numLoaded  = 0;

  FILE   *F = AS_UTL_openOutputFile(prefix, '.', "ignored.reads");

  for (uint32 fi=1; fi<=_numReads; fi++) {
    uint32   len  = seqStore->sqStore_getReadLength(fi);

    if (seqStore->sqStore_isIgnoredRead(fi)) {
      fprintf(F, "%u ignored\n", fi);
      numSkipped++;
      continue;
    }

    if (len < minReadLen) {
      fprintf(F, "%u length %u too short\n", fi, len);
      numSkipped++;
      continue;
    }

    _numBases += len;

    _readStatus[fi].readLength = len;
    _readStatus[fi].libraryID  = seqStore->sqStore_getLibraryIDForRead(fi);

    numLoaded++;
  }

  fclose(F);

  delete seqStore;

  if (minReadLen > 0)
    writeStatus("ReadInfo()-- Using %d reads, ignoring %u reads less than " F_U32 " bp long.\n",
                numLoaded, numSkipped, minReadLen);
  else
    writeStatus("ReadInfo()-- Using %d reads, no minimum read length used.\n",
                numLoaded);
}



ReadInfo::~ReadInfo() {
  delete [] _readStatus;
}
