
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
 *  This file is derived from:
 *
 *    src/AS_BAT/AS_BAT_FragmentInfo.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2010-NOV-23 to 2013-AUG-01
 *      are Copyright 2010-2013 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-DEC-19 to 2015-JUN-16
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2016-JAN-11
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

  _numLibraries = gkp->gkStore_getNumLibraries();
  _numReads     = gkp->gkStore_getNumReads();

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

      _readLength[iid] = read->gkRead_sequenceLength();
      _libIID[iid]     = lib;

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
