
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

  sqStore  *seqStore = sqStore::sqStore_open(seqStorePath);

  _numBases     = 0;
  _numReads     = seqStore->sqStore_getNumReads();
  _numLibraries = seqStore->sqStore_getNumLibraries();

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

  for (uint32 fi=1; fi<=_numReads; fi++) {
    sqRead  *read = seqStore->sqStore_getRead(fi);
    uint32   iid  = read->sqRead_readID();
    uint32   len  = read->sqRead_sequenceLength();

    if (len < minReadLen) {
      numSkipped++;
      continue;
    }

    _numBases += len;

    _readStatus[iid].readLength = len;
    _readStatus[iid].libraryID  = read->sqRead_libraryID();

    numLoaded++;
  }

  seqStore->sqStore_close();

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
