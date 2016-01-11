
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
 *    src/AS_CNS/MultiAlignStore.C
 *
 *  Modifications by:
 *
 *    Brian P. Walenz from 2009-OCT-05 to 2014-MAR-31
 *      are Copyright 2009-2014 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Sergey Koren on 2010-FEB-05
 *      are Copyright 2010 J. Craig Venter Institute, and
 *      are subject to the GNU General Public License version 2
 *
 *    Brian P. Walenz from 2014-DEC-22 to 2015-AUG-11
 *      are Copyright 2014-2015 Battelle National Biodefense Institute, and
 *      are subject to the BSD 3-Clause License
 *
 *    Brian P. Walenz beginning on 2015-OCT-29
 *      are a 'United States Government Work', and
 *      are released in the public domain
 *
 *  File 'README.licenses' in the root directory of this distribution contains
 *  full conditions and disclaimers for each license.
 */

#include "AS_global.H"
#include "AS_UTL_fileIO.H"
#include "tgStore.H"

uint32  MASRmagic   = 0x5253414d;  //  'MASR', as a big endian integer
uint32  MASRversion = 1;

#define MAX_VERS   1024  //  Linked to 10 bits in the header file.


tgStore::tgStore(const char *path_,
                 uint32      version_,
                 tgStoreType type_) {

  //  Handle goofy default parameters.  These let us retain the previous behavior (before tgStoreType):
  //    new tgStore("path")    - to create a new store
  //    new tgStore("path", v) - to open an existing store
  //
  //  And still allow
  //    new tgStore("path", v, tgStoreCreate) - create new store, make v the current version

  if (version_ == 0) {
    version_ = 1;
    type_    = tgStoreCreate;
  }

  //  Initialize the object.

  _type = type_;

  strcpy(_path, path_);

  _newTigs           = false;

  _currentVersion    = version_;
  _originalVersion   = version_;

  _tigMax            = 0;
  _tigLen            = 0;
  _tigEntry          = NULL;
  _tigCache          = NULL;

  _dataFile          = new dataFileT [MAX_VERS];

  for (uint32 i=0; i<MAX_VERS; i++) {
    _dataFile[i].FP = NULL;
    _dataFile[i].atEOF = false;
  }

  //  Create a new one?

  if (type_ == tgStoreCreate) {
    AS_UTL_mkdir(_path);                 //  Create the directory, if needed.
    purgeCurrentVersion();               //  Purge any data there currently.

    for (int32 vv=1; vv<version_; vv++)  //  Move to the requested version.
      nextVersion();

    _type = tgStoreWrite;                //  Created, ready for action!

    return;  //  No tigs to load, se we can't do the rest.
  }

  //  Load the tgStoreEntrys for the current version.

  loadMASR(_tigEntry, _tigLen, _tigMax, _currentVersion);

  //  Allocate the cache to the proper size

  _tigCache = new tgTig * [_tigMax];

  for (uint32 xx=0; xx<_tigMax; xx++)
    _tigCache[xx] = NULL;

  //  Check that nothing is marked for flushing, if so, clear the flag.  This shouldn't ever trigger.

  for (uint32 xx=0; xx<_tigLen; xx++)
    if (_tigEntry[xx].flushNeeded != 0) {
      fprintf(stderr, "WARNING: flushNeeded on tig %u\n", xx);
      _tigEntry[xx].flushNeeded = 0;
    }

  //  Set the version for writing.

  switch (type_) {
    case tgStoreCreate:
      break;

    case tgStoreReadOnly:
      if (_tigLen == 0) {
        fprintf(stderr, "tgStore::tgStore()-- ERROR, didn't find any tigs in the store.\n");
        fprintf(stderr, "tgStore::tgStore()--        asked for store '%s', correct?\n", _path);
        fprintf(stderr, "tgStore::tgStore()--        asked for version '%d', correct?\n", _originalVersion);
        exit(1);
      }
      break;

    case tgStoreWrite:
      _currentVersion++;      //  Writes go to the next version.
      purgeCurrentVersion();  //  And clear it.
      break;

    case tgStoreAppend:
      _currentVersion++;      //  Writes go to the next version.
      break;

    case tgStoreModify:
      break;

    default:
      assert(0);
      break;
  }


  //  Fail (again?) if there are no tigs loaded.

  //if (_tigLen == 0) {
  //  fprintf(stderr, "tgStore::tgStore()-- ERROR, didn't find any tigs in the store.  Correct version?\n");
  //  exit(1);
  //}
}


tgStore::~tgStore() {

  flushCache();

  //  If writable, write the data.

  if ((_type == tgStoreWrite) ||
      (_type == tgStoreAppend) ||
      (_type == tgStoreModify))
    dumpMASR(_tigEntry, _tigLen, _currentVersion);

  //  Now just trash ourself.

  delete [] _tigEntry;
  delete [] _tigCache;

  for (uint32 v=0; v<MAX_VERS; v++)
    if (_dataFile[v].FP)
      fclose(_dataFile[v].FP);

  delete [] _dataFile;
}



void
tgStore::purgeVersion(uint32 version) {

  sprintf(_name, "%s/seqDB.v%03d.dat", _path, version);   AS_UTL_unlink(_name);
  sprintf(_name, "%s/seqDB.v%03d.ctg", _path, version);   AS_UTL_unlink(_name);
  sprintf(_name, "%s/seqDB.v%03d.utg", _path, version);   AS_UTL_unlink(_name);
}


void
tgStore::purgeCurrentVersion(void) {
  purgeVersion(_currentVersion);
}



void
tgStore::nextVersion(void) {

  //  Write out any tigs that are cached

  flushDisk();

  //  Dump the MASR's.

  dumpMASR(_tigEntry, _tigLen, _currentVersion);

  //  Close the current version; we'll reopen on demand.

  if (_dataFile[_currentVersion].FP) {
    errno = 0;
    fclose(_dataFile[_currentVersion].FP);
    if (errno)
      fprintf(stderr, "tgStore::nextVersion()-- Failed to close '%s': %s\n", _name, strerror(errno)), exit(1);

    _dataFile[_currentVersion].FP    = NULL;
    _dataFile[_currentVersion].atEOF = false;
  }

  //  Bump to the next version.

  //fprintf(stderr, "tgStore::tgStore()-- moving from version %d to version %d; _tigLen %d _ctgLen %d\n",
  //        _currentVersion, _currentVersion+1, _tigLen, _ctgLen);

  _currentVersion++;

  if (_currentVersion == MAX_VERS) {
    fprintf(stderr, "Too many version, I can't proceed.\n");
    exit(1);
  }

  //  Remove any existing files at that version level.

  purgeCurrentVersion();
}



void
tgStore::writeTigToDisk(tgTig *tig, tgStoreEntry *te) {

  assert(_type != tgStoreReadOnly);

  FILE *FP = openDB(te->svID);

  //  The atEOF flag allows us to skip a seek when we're already (supposed) to be at the EOF.  This
  //  (hopefully) fixes a problem on one system where the seek() was placing the FP just before EOF
  //  (almost like the last block of data wasn't being flushed), and the tell() would then place the
  //  next tig in the middle of the previous one.
  //
  //  It also should (greatly) improve performance over NFS, espeically during BOG and CNS.  Both of
  //  these only write data, so no repositioning of the stream is needed.
  //
  if (_dataFile[te->svID].atEOF == false) {
    AS_UTL_fseek(FP, 0, SEEK_END);
    _dataFile[te->svID].atEOF = true;
  }

  te->flushNeeded = 0;
  te->fileOffset  = AS_UTL_ftell(FP);

  //fprintf(stderr, "tgStore::writeTigToDisk()-- write tig "F_S32" in store version "F_U64" at file position "F_U64"\n",
  //        tig->_tigID, te->svID, te->fileOffset);

  tig->saveToStream(FP);
}



void
tgStore::insertTig(tgTig *tig, bool keepInCache) {

  //  Check that the components do not exceed the bound.
  //
  if (tig->_gappedLen > 0) {
    uint32  len = tig->_gappedLen;
    uint32  swp = 0;
    uint32  neg = 0;
    uint32  pos = 0;

    for (uint32 i=0; i<tig->_childrenLen; i++) {
      tgPosition *read = tig->_children + i;

      if ((read->_max < read->_min))
        fprintf(stderr, "tgStore::insertTig()-- ERROR:   tig %d read %d at (%d,%d) has swapped min/max coordinates\n",
                tig->_tigID, read->_objID, read->_min, read->_max), swp++;
      //  Could fix, but we currently just fail.  This is an algorithmic problem that should be fixed.

      if ((read->_min < 0) || (read->_max < 0))
        fprintf(stderr, "tgStore::insertTig()-- WARNING: tig %d read %d at (%d,%d) has negative position\n",
                tig->_tigID, read->_objID, read->_min, read->_max), neg++;
      if (read->_min < 0)  read->_min = 0;
      if (read->_max < 0)  read->_max = 0;

      if ((read->_min > len) || (read->_max > len))
        fprintf(stderr, "tgStore::insertTig()-- WARNING: tig %d read %d at (%d,%d) exceeded multialign length %d\n",
                tig->_tigID, read->_objID, read->_min, read->_max, len), pos++;
      if (read->_min > len)  read->_min = len;
      if (read->_max > len)  read->_max = len;
    }

#if 0
    if (swp + neg + pos > 0) {
      tig->dumpLayout(stderr);
      fprintf(stderr, "tgStore::insertTig()-- ERROR: tig %d has invalid layout, exceeds bounds of consensus sequence (length %d) -- neg=%d pos=%d -- swp=%d.\n",
              tig->_tigID, len, neg, pos, swp);
    }
#endif
    assert(swp == 0);
    //assert(neg == 0);
    //assert(pos == 0);
  }

  if (tig->_tigID == UINT32_MAX) {
    tig->_tigID = _tigLen;
    _newTigs  = true;

    fprintf(stderr, "tgStore::insertTig()-- Added new tig %d\n", tig->_tigID);
  }

  if (_tigMax <= tig->_tigID) {
    while (_tigMax <= tig->_tigID)
      _tigMax = (_tigMax == 0) ? (1024) : (2 * _tigMax);
    assert(tig->_tigID < _tigMax);

    tgStoreEntry    *nr = new tgStoreEntry [_tigMax];
    tgTig          **nc = new tgTig *      [_tigMax];

    memcpy(nr, _tigEntry, sizeof(tgStoreEntry) * _tigLen);
    memcpy(nc, _tigCache, sizeof(tgTig *)      * _tigLen);

    memset(nr + _tigLen, 0, sizeof(tgStoreEntry) * (_tigMax - _tigLen));
    memset(nc + _tigLen, 0, sizeof(tgTig *)      * (_tigMax - _tigLen));

    for (uint32 xx=_tigLen; xx<_tigMax; xx++) {
      nr[xx].isDeleted = true;  //  Deleted until it gets added, otherwise we try to load and fail.
      nc[xx]           = NULL;
    }

    delete [] _tigEntry;
    delete [] _tigCache;

    _tigEntry = nr;
    _tigCache = nc;
  }

  _tigLen = MAX(_tigLen, tig->_tigID + 1);

  _tigEntry[tig->_tigID].tigRecord       = *tig;

  _tigEntry[tig->_tigID].unusedFlags     = 0;
  _tigEntry[tig->_tigID].flushNeeded     = true;   //  Mark as needing a flush by default
  _tigEntry[tig->_tigID].isDeleted       = false;  //  Now really here!
  _tigEntry[tig->_tigID].svID            = _currentVersion;
  _tigEntry[tig->_tigID].fileOffset      = 123456789;


  //  Write to disk RIGHT NOW unless we're keeping it in cache.  If it is written, the flushNeeded
  //  flag is cleared.
  //
  if ((keepInCache == false) && (_type != tgStoreReadOnly))
    writeTigToDisk(tig, _tigEntry + tig->_tigID);

  //  If the cache is different from this tig, delete the cache.  Not sure why this happens --
  //  did we copy a tig, muck with it, and then want to replace the one in the store?
  //
  if ((_tigCache[tig->_tigID] != tig) && (_tigCache[tig->_tigID] != NULL)) {
    delete _tigCache[tig->_tigID];
    _tigCache[tig->_tigID] = NULL;
  }

  //  Cache it if requested, otherwise clear the cache.
  //
  _tigCache[tig->_tigID] = (keepInCache) ? tig : NULL;
}



void
tgStore::deleteTig(uint32 tigID) {
  assert(tigID <  _tigLen);

  flushDisk(tigID);

  assert(_tigEntry[tigID].flushNeeded == 0);

  assert(_tigEntry[tigID].isDeleted == 0);

  _tigEntry[tigID].isDeleted = 1;

  delete [] _tigCache[tigID];
  _tigCache[tigID] = NULL;
}



tgTig *
tgStore::loadTig(uint32 tigID) {
  bool              cantLoad = true;

  if (_tigLen <= tigID)
    fprintf(stderr, "tgStore::loadTig()-- WARNING: invalid out-of-range tigID "F_S32", only "F_S32" ma in store; return NULL.\n",
            tigID, _tigLen);
  assert(tigID < _tigLen);

  //fprintf(stderr, "tgStore::loadTig()-- Loading tig %u (out of %u) from version %u at offest %lu\n",
  //        tigID, _tigLen,
  //        _tigEntry[tigID].svID,
  //        _tigEntry[tigID].fileOffset);

  //  This is...and is not...an error.  It does indicate something didn't go according to plan, like
  //  loading a tig that doesn't exist (that should be caught by the above 'tigID < _tigLen'
  //  assert).

  if (_tigEntry[tigID].isDeleted == true)
    return(NULL);

  //  This _is_ an error.  If a tig is in version zero, it isn't in the store at all.
  //  Someone did something stupid when adding tigs.

  if (_tigEntry[tigID].svID == 0)
    return(NULL);

  //  Otherwise, we can load something.

  if (_tigCache[tigID] == NULL) {
    FILE *FP = openDB(_tigEntry[tigID].svID);

    //  Since the tig isn't in the cache, it had better NOT be marked as needing to be flushed!
    assert(_tigEntry[tigID].flushNeeded == false);

    //  Seek to the correct position, and reset the atEOF to indicate we're (with high probability)
    //  not at EOF anymore.
    if (_dataFile[_tigEntry[tigID].svID].atEOF == true) {
      fflush(FP);
      _dataFile[_tigEntry[tigID].svID].atEOF = false;
    }

    AS_UTL_fseek(FP, _tigEntry[tigID].fileOffset, SEEK_SET);

    _tigCache[tigID] = new tgTig;
    _tigCache[tigID]->loadFromStream(FP);

    //  ALWAYS assume the incore record is more up to date
    *_tigCache[tigID] = _tigEntry[tigID].tigRecord;

    //  Since we just loaded, no flush is needed.
    _tigEntry[tigID].flushNeeded = 0;
  }

  return(_tigCache[tigID]);
}



void
tgStore::unloadTig(uint32 tigID, bool discardChanges) {

  if (discardChanges)
    _tigEntry[tigID].flushNeeded = 0;

  flushDisk(tigID);

  assert(_tigEntry[tigID].flushNeeded == 0);

  delete _tigCache[tigID];
  _tigCache[tigID] = NULL;
}


void
tgStore::copyTig(uint32 tigID, tgTig *tigcopy) {

  assert(tigID <  _tigLen);

  //  Deleted?  Clear it and return.

  if (_tigEntry[tigID].isDeleted) {
    tigcopy->clear();
    return;
  }

  //  In the cache?  Deep copy it and return.

  if (_tigCache[tigID]) {
    *tigcopy = *_tigCache[tigID];
    return;
  }

  //  Otherwise, load from disk.

  FILE *FP = openDB(_tigEntry[tigID].svID);

  //  Seek to the correct position, and reset the atEOF to indicate we're (with high probability)
  //  not at EOF anymore.

  if (_dataFile[_tigEntry[tigID].svID].atEOF == true) {
    fflush(FP);
    _dataFile[_tigEntry[tigID].svID].atEOF = false;
  }

  AS_UTL_fseek(FP, _tigEntry[tigID].fileOffset, SEEK_SET);

  tigcopy->clear();
  tigcopy->loadFromStream(FP);

  //  ALWAYS assume the incore record is more up to date
  *tigcopy = _tigEntry[tigID].tigRecord;
}



void
tgStore::flushDisk(uint32 tigID) {

  if (_tigEntry[tigID].flushNeeded == 0)
    return;

  writeTigToDisk(_tigCache[tigID], _tigEntry+tigID);
}



void
tgStore::flushDisk(void) {

  for (uint32 tigID=0; tigID<_tigLen; tigID++)
    if ((_tigCache[tigID]) && (_tigEntry[tigID].flushNeeded))
      flushDisk(tigID);
}



void
tgStore::flushCache(void) {

  flushDisk();

  for (uint32 i=0; i<_tigLen; i++)
    if (_tigCache[i]) {
      delete _tigCache[i];
      _tigCache[i] = NULL;
    }
}



uint32
tgStore::numTigsInMASRfile(char *name) {
  uint32        MASRmagicInFile   = 0;
  uint32        MASRversionInFile = 0;
  uint32        MASRtotalInFile   = 0;

  uint32        indxLen = 0;
  uint32        masrLen = 0;

  if (AS_UTL_fileExists(name, false, false) == false)
    return(0);

  errno = 0;
  FILE *F = fopen(name, "r");
  if (errno)
    fprintf(stderr, "tgStore::numTigsInMASRfile()-- Failed to open '%s': %s\n", name, strerror(errno)), exit(1);

  AS_UTL_safeRead(F, &MASRmagicInFile,   "MASRmagic",   sizeof(uint32), 1);
  AS_UTL_safeRead(F, &MASRversionInFile, "MASRversion", sizeof(uint32), 1);
  AS_UTL_safeRead(F, &MASRtotalInFile,   "MASRtotal",   sizeof(uint32), 1);

  fclose(F);

  if (MASRmagicInFile != MASRmagic) {
    fprintf(stderr, "tgStore::numTigsInMASRfile()-- Failed to open '%s': magic number mismatch; file=0x%08x code=0x%08x\n",
            name, MASRmagicInFile, MASRmagic);
    exit(1);
  }

  if (MASRversionInFile != MASRversion) {
    fprintf(stderr, "tgStore::numTigsInMASRfile()-- Failed to open '%s': version number mismatch; file=%d code=%d\n",
            name, MASRversionInFile, MASRversion);
    exit(1);
  }

  return(MASRtotalInFile);
}

void
tgStore::dumpMASR(tgStoreEntry* &R, uint32& L, uint32 V) {

  sprintf(_name, "%s/seqDB.v%03d.tig", _path, V);

  errno = 0;
  FILE *F = fopen(_name, "w");
  if (errno)
    fprintf(stderr, "tgStore::dumpMASR()-- Failed to create '%s': %s\n", _name, strerror(errno)), exit(1);

  AS_UTL_safeWrite(F, &MASRmagic,   "MASRmagic",   sizeof(uint32), 1);
  AS_UTL_safeWrite(F, &MASRversion, "MASRversion", sizeof(uint32), 1);
  AS_UTL_safeWrite(F, &L,           "MASRtotal",   sizeof(uint32), 1);

  uint32  indxLen = 0;

  //fprintf(stderr, "tgStore::dumpMASR()-- Writing '%s' (indxLen=%d masrLen=%d).\n", _name, indxLen, L);

  //  The max isn't written.  On load, max is set to length.

  AS_UTL_safeWrite(F, &indxLen, "MASRindexLen", sizeof(uint32),       1);
  AS_UTL_safeWrite(F, &L,       "MASRlen",      sizeof(uint32),       1);
  AS_UTL_safeWrite(F,  R,       "MASR",         sizeof(tgStoreEntry), L);

  fclose(F);
}


void
tgStore::loadMASR(tgStoreEntry* &R, uint32& L, uint32& M, uint32 V) {

  //  Allocate space for the data.  Search for the latest version, ask it how many tigs are in the
  //  store.
  //
  //  We don't always need to do this, sometimes we're called to update the data with newer data.
  //
  if (R == NULL) {
    for (int32 i=V; i>0; i--) {
      sprintf(_name, "%s/seqDB.v%03d.tig", _path, i);
      L = numTigsInMASRfile(_name);
      if (L > 0)
        break;
    }

    M = L + 1024;
    R = new tgStoreEntry [M];
    memset(R, 0, sizeof(tgStoreEntry) * M);
  }

  sprintf(_name, "%s/seqDB.v%03d.tig", _path, V);

  while ((AS_UTL_fileExists(_name) == false) && (V > 0))  {
    V--;
    sprintf(_name, "%s/seqDB.v%03d.tig", _path, V);
  }

  if (V == 0)
    fprintf(stderr, "tgStore::loadMASR()-- Failed to find any tigs in store '%s'.\n", _path), exit(1);

  errno = 0;
  FILE *F = fopen(_name, "r");
  if (errno)
    fprintf(stderr, "tgStore::loadMASR()-- Failed to open '%s': %s\n", _name, strerror(errno)), exit(1);

  uint32        MASRmagicInFile   = 0;
  uint32        MASRversionInFile = 0;
  uint32        MASRtotalInFile   = 0;

  uint32        indxLen = 0;
  uint32        masrLen = 0;

  AS_UTL_safeRead(F, &MASRmagicInFile,   "MASRmagic",   sizeof(uint32), 1);
  AS_UTL_safeRead(F, &MASRversionInFile, "MASRversion", sizeof(uint32), 1);
  AS_UTL_safeRead(F, &MASRtotalInFile,   "MASRtotal",   sizeof(uint32), 1);
  AS_UTL_safeRead(F, &indxLen,           "MASRindxLen", sizeof(uint32), 1);
  AS_UTL_safeRead(F, &masrLen,           "MASRmasrLen", sizeof(uint32), 1);

  if (MASRmagicInFile != MASRmagic) {
    fprintf(stderr, "tgStore::loadMASR()-- Failed to open '%s': magic number mismatch; file=0x%08x code=0x%08x\n",
            _name, MASRmagicInFile, MASRmagic);
    exit(1);
  }

  if (MASRversionInFile != MASRversion) {
    fprintf(stderr, "tgStore::loadMASR()-- Failed to open '%s': version number mismatch; file=%d code=%d\n",
            _name, MASRversionInFile, MASRversion);
    exit(1);
  }

  //  Check we're consistent.
  if (L < MASRtotalInFile)
    fprintf(stderr, "tgStore::loadMASR()-- '%s' has more tigs ("F_U32") than expected ("F_U32").\n",
            _name, MASRtotalInFile, L), exit(1);

  AS_UTL_safeRead(F,  R, "MASR", sizeof(tgStoreEntry), masrLen);

  fclose(F);
}




FILE *
tgStore::openDB(uint32 version) {

  if (_dataFile[version].FP)
    return(_dataFile[version].FP);

  //  Load the data

  sprintf(_name, "%s/seqDB.v%03d.dat", _path, version);

  //  If version is the _currentVersion, open for writing if allowed.
  //
  //  "a+" technically writes (always) to the end of file, but this hasn't been tested.

  errno = 0;

  if ((_type != tgStoreReadOnly) && (version == _currentVersion)) {
    _dataFile[version].FP    = fopen(_name, "a+");
    _dataFile[version].atEOF = false;
  } else {
    _dataFile[version].FP    = fopen(_name, "r");
    _dataFile[version].atEOF = false;
  }

  if (errno)
    fprintf(stderr, "tgStore::openDB()-- Failed to open '%s': %s\n", _name, strerror(errno)), exit(1);

  return(_dataFile[version].FP);
}
