
/**************************************************************************
 * This file is part of Celera Assembler, a software program that
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 1999-2004, Applera Corporation. All rights reserved.
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

#ifndef TGSTORE_H
#define TGSTORE_H

static const char *rcsid_TGSTORE_H = "$Id$";

#include "AS_global.H"
#include "tgTig.H"
//
//  The tgStore is a disk-resident (with memory cache) database of tgTig structures.
//

class tgStore {
public:

  //  Create a tgStore (first constructor).
  //  Open a tgStore (second constructor).
  //
  tgStore(const char *path);
  tgStore(const char *path,
          uint32      version,
          bool        writable=false,
          bool        inplace=false,
          bool        append=false);
  ~tgStore();

  //  Update to the next version.
  //
  void           nextVersion(void);

  //  Add or update a MA in the store.  If keepInCache, we keep a pointer to the tgTig.  THE
  //  STORE NOW OWNS THE OBJECT.
  //
  void           insertTig(tgTig *ma, bool keepInCache);

  //  delete() removes the tig from the cache, and marks it as deleted in the store.
  //
  void           deleteTig(uint32 tigID);

  //  load() will load and cache the MA.  THE STORE OWNS THIS OBJECT.
  //  copy() will load and copy the MA.  It will not cache.  YOU OWN THIS OBJECT.
  //
  tgTig         *loadTig(uint32 tigID);
  void           unloadTig(uint32 tigID, bool discard=false);

  void           copyTig(uint32 tigID, tgTig *ma);

  //  Flush to disk any cached MAs.  This is called by flushCache().
  //
  void           flushDisk(uint32 tigID);
  void           flushDisk(void);

  //  Flush the cache of loaded MAs.  Be aware that this is expensive in that the flushed things
  //  usually just get loaded back into core.
  //
  void           flushCache(uint32 tigID, bool discard=false) { unloadTig(tigID, discard); };
  void           flushCache(void);

  uint32         numTigs(void) { return(_tigLen); };

  //  Accessors to tig data; these do not load the tig from disk.

  bool           isDeleted(uint32 tigID);

  double         getCoverageStat(uint32 tigID);
  double         getMicroHetProb(uint32 tigID);

  bool           getSuggestRepeat(uint32 tigID);
  bool           getSuggestUnique(uint32 tigID);
  bool           getSuggestCircular(uint32 tigID);
  bool           getSuggestHaploid(uint32 tigID);

  uint32         getNumChildren(uint32 tigID);

  void           setCoverageStat(uint32 tigID, double cs);
  void           setMicroHetProb(uint32 tigID, double mp);

  void           setSuggestRepeat(uint32 tigID, bool enable=true);
  void           setSuggestUnique(uint32 tigID, bool enable=true);
  void           setSuggestCircular(uint32 tigID, bool enable=true);
  void           setSuggestHaploid(uint32 tigID, bool enable=true);

  uint32         getVersion(uint32 tigID);

private:
  struct tgStoreEntry {
    tgTigRecord  tigRecord;
    uint64       unusedFlags : 12;  //  One whole bit for future use.
    uint64       flushNeeded : 1;   //  If true, this MAR and associated tig are NOT saved to disk.
    uint64       isDeleted   : 1;   //  If true, this MAR has been deleted from the assembly.
    uint64       svID        : 10;  //  10 -> 1024 versions
    uint64       fileOffset  : 40;  //  40 -> 1 TB file size; offset in file where MA is stored
  };

  void                    init(const char *path_, uint32 version_, bool writable_, bool inplace_, bool append_);

  void                    writeTigToDisk(tgTig *ma, tgStoreEntry *maRecord);

  uint32                  numTigsInMASRfile(char *name);

  void                    dumpMASR(tgStoreEntry* &R, uint32& L, uint32& M, uint32 V);
  void                    loadMASR(tgStoreEntry* &R, uint32& L, uint32& M, uint32 V, bool onlyThisV);

  void                    purgeVersion(uint32 version);
  void                    purgeCurrentVersion(void);

  friend void operationCompress(char *tigName, int tigVers);

  FILE                   *openDB(uint32 V);

  char                    _path[FILENAME_MAX];
  char                    _name[FILENAME_MAX];

  bool                    _writable;               //  We are able to write
  bool                    _inplace;                //  We read and write to the same version
  bool                    _append;                 //  Do not nuke an existing version

  bool                    _newTigs;                //  internal flag, set if tigs were added

  uint32                  _originalVersion;        //  Version we started from (see newTigs in code)
  uint32                  _currentVersion;         //  Version we are writing to

  uint32                  _tigMax;
  uint32                  _tigLen;
  tgStoreEntry           *_tigEntry;
  tgTig                 **_tigCache;

  struct dataFileT {
    FILE   *FP;
    bool    atEOF;
  };

  dataFileT              *_dataFile;       //  dataFile[version]
};


inline
bool
tgStore::isDeleted(uint32 tigID) {
  assert(tigID >= 0);
  return(_tigEntry[tigID].isDeleted);
}

inline
double
tgStore::getCoverageStat(uint32 tigID) {
  assert(tigID >= 0);
  assert(tigID < _tigLen);
  return(_tigEntry[tigID].tigRecord._coverageStat);
}

inline
double
tgStore::getMicroHetProb(uint32 tigID) {
  assert(tigID >= 0);
  assert(tigID < _tigLen);
  return(_tigEntry[tigID].tigRecord._microhetProb);
}

inline
bool
tgStore::getSuggestRepeat(uint32 tigID) {
  assert(tigID >= 0);
  assert(tigID < _tigLen);
  return(_tigEntry[tigID].tigRecord._suggestRepeat);
}

inline
bool
tgStore::getSuggestUnique(uint32 tigID) {
  assert(tigID >= 0);
  assert(tigID < _tigLen);
  return(_tigEntry[tigID].tigRecord._suggestUnique);
}

inline
bool
tgStore::getSuggestCircular(uint32 tigID) {
  assert(tigID >= 0);
  assert(tigID < _tigLen);
  return(_tigEntry[tigID].tigRecord._suggestCircular);
}

inline
bool
tgStore::getSuggestHaploid(uint32 tigID) {
  assert(tigID >= 0);
  assert(tigID < _tigLen);
  return(_tigEntry[tigID].tigRecord._suggestHaploid);
}

inline
uint32
tgStore::getNumChildren(uint32 tigID) {
  assert(tigID >= 0);
  return(_tigEntry[tigID].tigRecord._childrenLen);
}

inline
void
tgStore::setCoverageStat(uint32 tigID, double cs) {
  assert(tigID >= 0);
  assert(tigID < _tigLen);
  _tigEntry[tigID].tigRecord._coverageStat = cs;
  if (_tigCache[tigID])
    _tigCache[tigID]->_coverageStat = cs;
}

inline
void
tgStore::setMicroHetProb(uint32 tigID, double mp) {
  assert(tigID >= 0);
  assert(tigID < _tigLen);
  _tigEntry[tigID].tigRecord._microhetProb = mp;
  if (_tigCache[tigID])
    _tigCache[tigID]->_microhetProb = mp;
}

inline
void
tgStore::setSuggestRepeat(uint32 tigID, bool enable) {
  assert(tigID >= 0);
  assert(tigID < _tigLen);
  _tigEntry[tigID].tigRecord._suggestRepeat = enable;
  if (_tigCache[tigID])
    _tigCache[tigID]->_suggestRepeat = enable;
}

inline
void
tgStore::setSuggestUnique(uint32 tigID, bool enable) {
  assert(tigID >= 0);
  assert(tigID < _tigLen);
  _tigEntry[tigID].tigRecord._suggestUnique = enable;
  if (_tigCache[tigID])
    _tigCache[tigID]->_suggestUnique = enable;
}

inline
void
tgStore::setSuggestCircular(uint32 tigID, bool enable) {
  assert(tigID >= 0);
  assert(tigID < _tigLen);
  _tigEntry[tigID].tigRecord._suggestCircular = enable;
  if (_tigCache[tigID])
    _tigCache[tigID]->_suggestCircular = enable;
}

inline
void
tgStore::setSuggestHaploid(uint32 tigID, bool enable) {
  assert(tigID >= 0);
  assert(tigID < _tigLen);
  _tigEntry[tigID].tigRecord._suggestHaploid = enable;
  if (_tigCache[tigID])
    _tigCache[tigID]->_suggestHaploid = enable;
}

inline
uint32
tgStore::getVersion(uint32 tigID) {
  assert(tigID >= 0);
  assert(tigID < _tigLen);
  return(_tigEntry[tigID].svID);
}

#endif
