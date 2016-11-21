
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



ReadInfo::ReadInfo(gkStore    *gkp,
                           const char *prefix,
                           uint32      minReadLen) {

  _numBases     = 0;
  _numReads     = gkp->gkStore_getNumReads();
  _numLibraries = gkp->gkStore_getNumLibraries();

  _readLength    = new uint32 [_numReads + 1];
  _libIID        = new uint32 [_numReads + 1];

  for (uint32 i=0; i<_numReads + 1; i++) {
    _readLength[i] = 0;
    _libIID[i] = 0;
  }

  uint32 numSkipped = 0;
  uint32 numLoaded  = 0;

  for (uint32 fi=1; fi<=_numReads; fi++) {
    gkRead  *read = gkp->gkStore_getRead(fi);

    if (read->gkRead_sequenceLength() < minReadLen) {
      numSkipped++;

    } else {
      uint32 iid = read->gkRead_readID();
      uint32 lib = read->gkRead_libraryID();

      _numBases        += read->gkRead_sequenceLength();
      _readLength[iid]  = read->gkRead_sequenceLength();
      _libIID[iid]      = lib;

      numLoaded++;
    }
  }

  if (minReadLen > 0)
    writeStatus("ReadInfo()-- Using %d reads, ignoring %u reads less than " F_U32 " bp long.\n",
                numLoaded, numSkipped, minReadLen);
  else
    writeStatus("ReadInfo()-- Using %d reads, no minimum read length used.\n",
                numLoaded);
}



ReadInfo::~ReadInfo() {
  delete [] _readLength;
  delete [] _libIID;
}
