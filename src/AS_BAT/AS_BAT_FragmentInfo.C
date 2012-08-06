
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, The Venter Institute. All rights reserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received (LICENSE.txt) a copy of the GNU General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *************************************************************************/

static const char *rcsid = "$Id: AS_BAT_FragmentInfo.C,v 1.5 2012-08-06 23:32:32 brianwalenz Exp $";

#include "AS_BAT_Datatypes.H"

const uint64 fiMagicNumber   = 0x6f666e4967617266llu;  //  'fragInfo' until it gets messed up by endianess.
const uint64 fiVersionNumber = 1;


FragmentInfo::FragmentInfo(gkStore *gkpStore, const char *prefix) {

  if (load(prefix))
    return;

  writeLog("FragmentInfo()-- Loading fragment information\n");

  gkStream         *fs = new gkStream(gkpStore, 0, 0, GKFRAGMENT_INF);
  gkFragment        fr;

  _numLibraries = gkpStore->gkStore_getNumLibraries();
  _numFragments = gkpStore->gkStore_getNumFragments();

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

  for (uint32 i=0; i<_numLibraries + 1; i++) {
    _mean[i]          = 0.0;
    _stddev[i]        = 0.0;
    _numFragsInLib[i] = 0;
    _numMatesInLib[i] = 0;
  }

  for (uint32 i=1; i<_numLibraries + 1; i++) {
    _mean[i]          = gkpStore->gkStore_getLibrary(i)->mean;
    _stddev[i]        = gkpStore->gkStore_getLibrary(i)->stddev;
    _numFragsInLib[i] = 0;
    _numMatesInLib[i] = 0;
  }

  uint32 numDeleted = 0;
  uint32 numLoaded  = 0;

  while(fs->next(&fr)) {
    if (fr.gkFragment_getIsDeleted()) {
      numDeleted++;
    } else {
      uint32 iid = fr.gkFragment_getReadIID();
      uint32 lib = fr.gkFragment_getLibraryIID();

      _fragLength[iid] = fr.gkFragment_getClearRegionLength();
      _mateIID[iid]    = fr.gkFragment_getMateIID();;
      _libIID[iid]     = lib;

      _numFragsInLib[lib]++;

      if (_mateIID[iid])
        _numMatesInLib[lib]++;

      numLoaded++;
    }

    if (((numDeleted + numLoaded) % 10000000) == 0)
      writeLog("FragmentInfo()-- Loading fragment information deleted:%9d active:%9d\n", numDeleted, numLoaded);
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
    
  writeLog("FragmentInfo()-- Loaded %d alive fragments, skipped %d dead fragments.\n", numLoaded, numDeleted);

  delete fs;

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
