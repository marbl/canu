
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

  //  loadTig() will load and cache the tig.  THE STORE OWNS THIS OBJECT.
  //  unloadTig() allows the store to release the cached tig.
  //
  //  copyTig() will load and copy the tig to the supplied object.  YOU OWN
  //  THIS OBJECT, and can delete it at will.  The object passed in is the
  //  return value (so you can create and load in one line:
  //    tgTig *tig = store->copyTig(new tgTig);
  //
  tgTig         *loadTig(uint32 tigID);
  void           unloadTig(uint32 tigID, bool discardChanges=false);

  tgTig         *copyTig(uint32 tigID, tgTig *tig);

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

  bool           isDeleted(uint32 tigID) {
    assert(tigID < _tigLen);
    return(_tigEntry[tigID].isDeleted);
  };

  uint32         getVersion(uint32 tigID) {
    assert(tigID < _tigLen);
    return(_tigEntry[tigID].svID);
  };

private:
  struct tgStoreEntryV1 {
    tgTigRecordV1  tigRecord;
    uint64         unusedFlags : 12;  //  One whole bit for future use.
    uint64         flushNeeded : 1;   //  If true, this MAR and associated tig are NOT saved to disk.
    uint64         isDeleted   : 1;   //  If true, this MAR has been deleted from the assembly.
    uint64         svID        : 10;  //  10 -> 1024 versions (HARDCODED in tgStore.C)
    uint64         fileOffset  : 40;  //  40 -> 1 TB file size; offset in file where MA is stored
  };
  struct tgStoreEntry {
    tgTigRecord    tigRecord;
    uint64         unusedFlags : 12;  //  One whole bit for future use.
    uint64         flushNeeded : 1;   //  If true, this MAR and associated tig are NOT saved to disk.
    uint64         isDeleted   : 1;   //  If true, this MAR has been deleted from the assembly.
    uint64         svID        : 10;  //  10 -> 1024 versions (HARDCODED in tgStore.C)
    uint64         fileOffset  : 40;  //  40 -> 1 TB file size; offset in file where MA is stored

    tgStoreEntry &
    operator=(const tgStoreEntryV1 &that) {
      tigRecord   = that.tigRecord;   //  Itself an overloaded operator=().

      unusedFlags = that.unusedFlags;
      flushNeeded = that.flushNeeded;
      isDeleted   = that.isDeleted;
      svID        = that.svID;
      fileOffset  = that.fileOffset;

      return *this;
    }
  };

  void                    writeTigToDisk(tgTig *ma, tgStoreEntry *maRecord);

  uint32                  numTigsInMASRfile(char *name);

  void                    dumpMASR(tgStoreEntry* &R, uint32& L,            uint32 V);
  void                    loadMASR(tgStoreEntry* &R, uint32& L, uint32& M, uint32 V);

  void                    purgeVersion(uint32 version);
  void                    purgeCurrentVersion(void);

  friend void operationCompress(char *tigName, int tigVers);   //  So it can get to purgeVersion().

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


#endif
