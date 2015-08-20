
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
 *    Brian P. Walenz beginning on 2014-DEC-19
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

static const char *rcsid = "$Id$";

#include "AS_BAT_Datatypes.H"

const uint64 fiMagicNumber   = 0x6f666e4967617266llu;  //  'fragInfo' until it gets messed up by endianess.
const uint64 fiVersionNumber = 1;


FragmentInfo::FragmentInfo(gkStore    *gkp,
                           const char *prefix,
                           uint32      minReadLen) {

  if (load(prefix))
    return;

  writeLog("FragmentInfo()-- Loading fragment information\n");

  if (minReadLen > 0)
    writeLog("FragmentInfo()-- Reads shorter than "F_U32" bases are forced to be singleton.\n",
             minReadLen);

  _numLibraries = gkp->gkStore_getNumLibraries();
  _numFragments = gkp->gkStore_getNumReads();

  _fragLength    = new uint32 [_numFragments + 1];
  _mateIID       = new uint32 [_numFragments + 1];
  _libIID        = new uint32 [_numFragments + 1];

  _mean          = new double [_numLibraries + 1];
  _stddev        = new double [_numLibraries + 1];

  _numFragsInLib = new uint32 [_numLibraries + 1];
  _numMatesInLib = new uint32 [_numLibraries + 1];

  for (uint32 i=0; i<_numFragments + 1; i++) {
    _fragLength[i] = 0;
    _mateIID[i] = 0;
    _libIID[i] = 0;
  }

  for (uint32 i=1; i<_numLibraries + 1; i++) {
    _mean[i]          = 0.0;
    _stddev[i]        = 0.0;
    //_mean[i]          = gkp->gkStore_getLibrary(i)->mean;
    //_stddev[i]        = gkp->gkStore_getLibrary(i)->stddev;
    _numFragsInLib[i] = 0;
    _numMatesInLib[i] = 0;
  }

  uint32 numSkipped = 0;
  uint32 numLoaded  = 0;

  for (uint32 fi=1; fi<=_numFragments; fi++) {
    gkRead  *read = gkp->gkStore_getRead(fi);

    if (read->gkRead_sequenceLength() < minReadLen) {
      numSkipped++;

    } else {
      uint32 iid = read->gkRead_readID();
      uint32 lib = read->gkRead_libraryID();

      _fragLength[iid] = read->gkRead_sequenceLength();
      _mateIID[iid]    = 0;  //read->gkRead_mateIID();;
      _libIID[iid]     = lib;

      _numFragsInLib[lib]++;

      if (_mateIID[iid])
        _numMatesInLib[lib]++;

      numLoaded++;
    }

    if (((numSkipped + numLoaded) % 10000000) == 0)
      writeLog("FragmentInfo()-- Loading fragment information: skipped:%9d active:%9d\n",
               numSkipped, numLoaded);
  }

  for (uint32 i=0; i<_numLibraries + 1; i++)
    _numMatesInLib[i] /= 2;

  //  Search for and break (and complain) mates to deleted fragments.
  uint32  numBroken = 0;

  for (uint32 i=0; i<_numFragments + 1; i++) {
    if ((_fragLength[i] == 0) ||
        (_mateIID[i] == 0) ||
        (_fragLength[_mateIID[i]] > 0))
      //  This frag deleted, or this frag unmated, or mate of this frag is alive, all good!
      continue;

    assert(_mateIID[_mateIID[i]] == 0);

    if (numBroken++ < 100)
      writeLog("FragmentInfo()-- WARNING!  Mate of fragment %d (fragment %d) is deleted.\n",
               i, _mateIID[i]);

    _mateIID[i] = 0;
  }

  if (numBroken > 0)
    writeLog("FragmentInfo()-- WARNING!  Removed "F_U32" mate relationships.\n", numBroken);

  writeLog("FragmentInfo()-- Loaded %d alive reads, skipped %d short reads.\n",
           numLoaded, numSkipped);

  save(prefix);
}



FragmentInfo::~FragmentInfo() {
  delete [] _fragLength;
  delete [] _mateIID;
  delete [] _libIID;

  delete [] _mean;
  delete [] _stddev;

  delete [] _numFragsInLib;
  delete [] _numMatesInLib;
}



void
FragmentInfo::save(const char *prefix) {
  char  name[FILENAME_MAX];

  sprintf(name, "%s.fragmentInfo", prefix);

  errno = 0;
  FILE *file = fopen(name, "w");
  if (errno) {
    writeLog("FragmentInfo()-- Failed to open '%s' for writing: %s\n", name, strerror(errno));
    writeLog("FragmentInfo()-- Will not save fragment information to cache.\n");
    return;
  }

  writeLog("FragmentInfo()-- Saving fragment information to cache '%s'\n", name);

  AS_UTL_safeWrite(file, &fiMagicNumber,   "fragmentInformationMagicNumber",  sizeof(uint64), 1);
  AS_UTL_safeWrite(file, &fiVersionNumber, "fragmentInformationMagicNumber",  sizeof(uint64), 1);
  AS_UTL_safeWrite(file, &_numFragments,   "fragmentInformationNumFrgs",      sizeof(uint32), 1);
  AS_UTL_safeWrite(file, &_numLibraries,   "fragmentInformationNumLibs",      sizeof(uint32), 1);

  AS_UTL_safeWrite(file,  _fragLength,     "fragmentInformationFragLen",      sizeof(uint32), _numFragments + 1);
  AS_UTL_safeWrite(file,  _mateIID,        "fragmentInformationMateIID",      sizeof(uint32), _numFragments + 1);
  AS_UTL_safeWrite(file,  _libIID,         "fragmentInformationLibIID",       sizeof(uint32), _numFragments + 1);

  AS_UTL_safeWrite(file,  _mean,           "fragmentInformationMean",         sizeof(double), _numLibraries + 1);
  AS_UTL_safeWrite(file,  _stddev,         "fragmentInformationStddev",       sizeof(double), _numLibraries + 1);
  AS_UTL_safeWrite(file,  _numFragsInLib,  "fragmentInformationNumFrgsInLib", sizeof(uint32), _numLibraries + 1);
  AS_UTL_safeWrite(file,  _numMatesInLib,  "fragmentInformationNumMateInLib", sizeof(uint32), _numLibraries + 1);

  fclose(file);
}


bool
FragmentInfo::load(const char *prefix) {
  char  name[FILENAME_MAX];

  sprintf(name, "%s.fragmentInfo", prefix);

  errno = 0;
  FILE *file = fopen(name, "r");
  if (errno)
    return(false);

  uint64  magicNumber   = 0;
  uint64  versionNumber = 0;

  AS_UTL_safeRead(file, &magicNumber,    "fragmentInformationMagicNumber",   sizeof(uint64), 1);
  AS_UTL_safeRead(file, &versionNumber,  "fragmentInformationVersionNumber", sizeof(uint64), 1);
  AS_UTL_safeRead(file, &_numFragments,  "fragmentInformationNumFrgs",       sizeof(uint32), 1);
  AS_UTL_safeRead(file, &_numLibraries,  "fragmentInformationNumLibs",       sizeof(uint32), 1);

  if (magicNumber != fiMagicNumber) {
    writeLog("FragmentInfo()-- File '%s' is not a fragment info; cannot load.\n", name);
    fclose(file);
    return(false);
  }
  if (versionNumber != fiVersionNumber) {
    writeLog("FragmentInfo()-- File '%s' is version "F_U64", I can only read version "F_U64"; cannot load.\n",
            name, versionNumber, fiVersionNumber);
    fclose(file);
    return(false);
  }

  writeLog("FragmentInfo()-- Loading fragment information for "F_U32" fragments and "F_U32" libraries from cache '%s'\n",
          _numFragments, _numLibraries, name);

  _fragLength    = new uint32 [_numFragments + 1];
  _mateIID       = new uint32 [_numFragments + 1];
  _libIID        = new uint32 [_numFragments + 1];

  _mean          = new double [_numLibraries + 1];
  _stddev        = new double [_numLibraries + 1];

  _numFragsInLib = new uint32 [_numLibraries + 1];
  _numMatesInLib = new uint32 [_numLibraries + 1];

  AS_UTL_safeRead(file,  _fragLength,    "fragmentInformationFragLen",      sizeof(uint32), _numFragments + 1);
  AS_UTL_safeRead(file,  _mateIID,       "fragmentInformationMateIID",      sizeof(uint32), _numFragments + 1);
  AS_UTL_safeRead(file,  _libIID,        "fragmentInformationLibIID",       sizeof(uint32), _numFragments + 1);

  AS_UTL_safeRead(file,  _mean,          "fragmentInformationMean",         sizeof(double), _numLibraries + 1);
  AS_UTL_safeRead(file,  _stddev,        "fragmentInformationStddev",       sizeof(double), _numLibraries + 1);
  AS_UTL_safeRead(file,  _numFragsInLib, "fragmentInformationNumFrgsInLib", sizeof(uint32), _numLibraries + 1);
  AS_UTL_safeRead(file,  _numMatesInLib, "fragmentInformationNumMateInLib", sizeof(uint32), _numLibraries + 1);

  fclose(file);

  return(true);
}
