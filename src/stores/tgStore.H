
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

#ifndef TGSTORE_H
#define TGSTORE_H

#include "runtime.H"
#include "tgTig.H"
//
//  The tgStore is a disk-resident (with memory cache) database of tgTig structures.
//
//  There are two basic modes of operation:
//    open a store for reading version v
//    open a store for reading version v, and writing to version v+1, erasing v+1 before starting
//    open a store for reading version v, and writing to version v+1, preserving the contents
//    open a store for reading version v, and writing to version v,   preserving the contents
//

enum tgStoreType {       //  writable  inplace  append
  tgStoreCreate    = 0,  //  Make a new one, then become tgStoreWrite
  tgStoreReadOnly  = 1,  //     false        *       * - open version v   for reading; inplace=append=false in the code
  tgStoreWrite     = 2,  //      true    false   false - open version v+1 for writing, purge contents of v+1; standard open for writing
  tgStoreAppend    = 3,  //      true    false    true - open version v+1 for writing, do not purge contents
  tgStoreModify    = 4,  //      true     true   false - open version v   for writing, do not purge contents
};



class tgStore {
public:
  tgStore(const char *path,
          uint32      version = 0,
          tgStoreType type    = tgStoreReadOnly);
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
  void           unloadTig(uint32 tigID, bool discardChanges=false);

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

  uint32         getSourceID(uint32 tigID);
  uint32         getSourceBgn(uint32 tigID);
  uint32         getSourceEnd(uint32 tigID);

  tgTig_class    getClass(uint32 tigID);
  bool           getSuggestRepeat(uint32 tigID);
  bool           getSuggestCircular(uint32 tigID);

  uint32         getNumChildren(uint32 tigID);

  void           setSourceID(uint32 tigID, uint32 id);
  void           setSourceBgn(uint32 tigID, uint32 bgn);
  void           setSourceEnd(uint32 tigID, uint32 end);

  void           setClass(uint32 tigID, tgTig_class c);
  void           setSuggestRepeat(uint32 tigID, bool enable=true);
  void           setSuggestCircular(uint32 tigID, bool enable=true);

  uint32         getVersion(uint32 tigID);

private:
  struct tgStoreEntry {
    tgTigRecord  tigRecord;
    uint64       unusedFlags : 12;  //  One whole bit for future use.
    uint64       flushNeeded : 1;   //  If true, this MAR and associated tig are NOT saved to disk.
    uint64       isDeleted   : 1;   //  If true, this MAR has been deleted from the assembly.
    uint64       svID        : 10;  //  10 -> 1024 versions (HARDCODED in tgStore.C)
    uint64       fileOffset  : 40;  //  40 -> 1 TB file size; offset in file where MA is stored
  };

  void                    writeTigToDisk(tgTig *ma, tgStoreEntry *maRecord);

  uint32                  numTigsInMASRfile(char *name);

  void                    dumpMASR(tgStoreEntry* &R, uint32& L,            uint32 V);
  void                    loadMASR(tgStoreEntry* &R, uint32& L, uint32& M, uint32 V);

  void                    purgeVersion(uint32 version);
  void                    purgeCurrentVersion(void);

  friend void operationCompress(char *tigName, int tigVers);

  FILE                   *openDB(uint32 V);

  char                    _path[FILENAME_MAX+1];   //  Path to the store.
  char                    _name[FILENAME_MAX+1];   //  Name of the currently opened file, and other uses.

  tgStoreType             _type;

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
  return(_tigEntry[tigID].isDeleted);
}

inline
uint32
tgStore::getSourceID(uint32 tigID) {
  assert(tigID < _tigLen);
  return(_tigEntry[tigID].tigRecord._sourceID);
}

inline
uint32
tgStore::getSourceBgn(uint32 tigID) {
  assert(tigID < _tigLen);
  return(_tigEntry[tigID].tigRecord._sourceBgn);
}

inline
uint32
tgStore::getSourceEnd(uint32 tigID) {
  assert(tigID < _tigLen);
  return(_tigEntry[tigID].tigRecord._sourceEnd);
}

inline
tgTig_class
tgStore::getClass(uint32 tigID) {
  assert(tigID < _tigLen);
  return(_tigEntry[tigID].tigRecord._class);
}

inline
bool
tgStore::getSuggestRepeat(uint32 tigID) {
  assert(tigID < _tigLen);
  return(_tigEntry[tigID].tigRecord._suggestRepeat);
}

inline
bool
tgStore::getSuggestCircular(uint32 tigID) {
  assert(tigID < _tigLen);
  return(_tigEntry[tigID].tigRecord._suggestCircular);
}

inline
uint32
tgStore::getNumChildren(uint32 tigID) {
  return(_tigEntry[tigID].tigRecord._childrenLen);
}



inline
void
tgStore::setSourceID(uint32 tigID, uint32 id) {
  assert(tigID < _tigLen);
  _tigEntry[tigID].tigRecord._sourceID = id;
  if (_tigCache[tigID])
    _tigCache[tigID]->_sourceID = id;
}

inline
void
tgStore::setSourceBgn(uint32 tigID, uint32 bgn) {
  assert(tigID < _tigLen);
  _tigEntry[tigID].tigRecord._sourceBgn = bgn;
  if (_tigCache[tigID])
    _tigCache[tigID]->_sourceBgn = bgn;
}

inline
void
tgStore::setSourceEnd(uint32 tigID, uint32 end) {
  assert(tigID < _tigLen);
  _tigEntry[tigID].tigRecord._sourceEnd = end;
  if (_tigCache[tigID])
    _tigCache[tigID]->_sourceEnd = end;
}

inline
void
tgStore::setClass(uint32 tigID, tgTig_class c) {
  assert(tigID < _tigLen);
  _tigEntry[tigID].tigRecord._class = c;
  if (_tigCache[tigID])
    _tigCache[tigID]->_class = c;
}

inline
void
tgStore::setSuggestRepeat(uint32 tigID, bool enable) {
  assert(tigID < _tigLen);
  _tigEntry[tigID].tigRecord._suggestRepeat = enable;
  if (_tigCache[tigID])
    _tigCache[tigID]->_suggestRepeat = enable;
}

inline
void
tgStore::setSuggestCircular(uint32 tigID, bool enable) {
  assert(tigID < _tigLen);
  _tigEntry[tigID].tigRecord._suggestCircular = enable;
  if (_tigCache[tigID])
    _tigCache[tigID]->_suggestCircular = enable;
}

inline
uint32
tgStore::getVersion(uint32 tigID) {
  assert(tigID < _tigLen);
  return(_tigEntry[tigID].svID);
}

#endif
