
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

static const char *rcsid = "$Id$";

#include "AS_global.H"
#include "AS_UTL_fileIO.H"
#include "tgStore.H"

uint32  MASRmagic   = 0x5253414d;  //  'MASR', as a big endian integer
uint32  MASRversion = 1;

#define MAX_VERS   1024
#define MAX_PART   1024

void
tgStore::init(const char *path_, uint32 version_, bool writable_, bool inplace_, bool append_) {

  strcpy(_path, path_);

  _writable          = writable_;
  _inplace           = inplace_;
  _append            = append_;

  _newTigs           = false;

  _currentVersion    = version_;
  _originalVersion   = version_;

  _part              = 0;
  _partMap           = NULL;
  _partMapLen        = 0;

  _tigMax            = 0;
  _tigLen            = 0;
  _tigEntry          = NULL;
  _tigCache          = NULL;

  //  Could use sysconf(_SC_OPEN_MAX) too.  Should make this dynamic?
  //
  _dataFile          = new dataFileT * [MAX_VERS];
  _dataFile[0]       = new dataFileT   [MAX_PART];

  for (uint32 i=1; i<MAX_VERS; i++)
    _dataFile[i] = _dataFile[0] + i * MAX_PART;
}


tgStore::tgStore(const char *path_) {
  init(path_, 1, true, false, false);
  AS_UTL_mkdir(_path);     //  Create the directory, if needed.
  purgeCurrentVersion();  //  Purge any data there currently.
}


tgStore::tgStore(const char *path_,
                 uint32      version_,
                 uint32      partition_,
                 bool        writable_,
                 bool        inplace_,
                 bool        append_) {

  if (writable_ == false)
    inplace_ = append_ = false;

  init(path_, version_, writable_, inplace_, append_);

  _part = partition_;

  if ((_inplace == true) && (_append == true))
    fprintf(stderr, "tgStore::tgStore()-- ERROR, cannot both append and be inplace.\n"), exit(1);

  //  Load the tgStoreEntrys for the current version.

  loadMASR(_tigEntry, _tigLen, _tigMax, _currentVersion, false);

  //  Check that there are tigs.

  if ((_tigLen == 0) && (_append == false)) {
    fprintf(stderr, "tgStore::tgStore()-- ERROR, didn't find any tigs in the store.\n");
    fprintf(stderr, "tgStore::tgStore()--        asked for store '%s', correct?\n", _path);
    fprintf(stderr, "tgStore::tgStore()--        asked for version '%d', correct?\n", _originalVersion);
    fprintf(stderr, "tgStore::tgStore()--        asked for partition=%d, correct?\n", _part);
    fprintf(stderr, "tgStore::tgStore()--        asked for writable=%d inplace=%d append=%d, correct?\n", _writable, _inplace, _append);
    exit(1);
  }

  //  Allocate the cache to the proper size

  _tigCache = new tgTig * [_tigMax];

  //  Open the next version for writing, and remove what is currently there.

  if ((_writable == true) && (_inplace == false) && (_append == false)) {
    _currentVersion++;

    purgeCurrentVersion();
  }

  //  Open the next version for writing, and keep the data that is currently there.

  if (_append == true)
    _currentVersion++;

  //  Check that nothing is marked for flushing, if so, clear the flag.  This shouldn't ever trigger.

  for (uint32 xx=0; xx<_tigLen; xx++)
    if (_tigEntry[xx].flushNeeded != 0) {
      fprintf(stderr, "WARNING: flushNeeded on tig %u\n", xx);
      _tigEntry[xx].flushNeeded = 0;
    }
  //fprintf(stderr, "tgStore::loadMASRfile()-- flushNeeded check finished on '%s'.\n", _name);

  //  Fail (again?) if there are no tigs loaded.

  if (_tigLen == 0) {
    fprintf(stderr, "tgStore::tgStore()-- ERROR, didn't find any tigs in the store.  Correct version?\n");
    exit(1);
  }
}


tgStore::~tgStore() {

  flushCache();

  //  If writable, write the data.

  if (_writable)
    dumpMASR(_tigEntry, _tigLen, _tigMax, _currentVersion);

  //  If new tigs were added, AND we are partitioned, update the global partition.
  //
  //  Not sure what this was doing.  It triggered if new tigs were added, and we were partitioned on unitigs.
  //  It would load the full set of unitig data, then dump it back out.
  //
  //if ((_newTigs) && (_part != 0)) {
  //  assert(_partMap == NULL);
  //
  //  _part = 0;  //  To load the unpartitioned MASR          v-isunitigs v-onlyThisV
  //  loadMASR(_tigEntry, _tigLen, _tigMax, _originalVersion, true, false);
  //  dumpMASR(_tigEntry, _tigLen, _tigMax, _originalVersion, true);
  //}

  delete [] _tigEntry;
  delete [] _tigCache;

  for (uint32 v=0; v<MAX_VERS; v++)
    for (uint32 p=0; p<MAX_PART; p++)
      if (_dataFile[v][p].FP)
        fclose(_dataFile[v][p].FP);

  delete [] _dataFile[0];
  delete [] _dataFile;
}



void
tgStore::purgeVersion(uint32 version) {
  uint32 part = 1;

  sprintf(_name, "%s/seqDB.v%03d.dat", _path, version);   AS_UTL_unlink(_name);
  sprintf(_name, "%s/seqDB.v%03d.ctg", _path, version);   AS_UTL_unlink(_name);
  sprintf(_name, "%s/seqDB.v%03d.utg", _path, version);   AS_UTL_unlink(_name);

  sprintf(_name, "%s/seqDB.v%03d.p001.dat", _path, version);
  while (AS_UTL_unlink(_name)) {
    sprintf(_name, "%s/seqDB.v%03d.p%03d.ctg", _path, version, part);   AS_UTL_unlink(_name);
    sprintf(_name, "%s/seqDB.v%03d.p%03d.utg", _path, version, part);   AS_UTL_unlink(_name);
    sprintf(_name, "%s/seqDB.v%03d.p%03d.dat", _path, version, ++part);
  }
}


void
tgStore::purgeCurrentVersion(void) {
  purgeVersion(_currentVersion);
}



void
tgStore::nextVersion(void) {

  assert(_writable == true);
  assert(_inplace == false);
  assert(_partMap == NULL);
  assert(_part == 0);

  //  Write out any tigs that are cached

  flushDisk();

  //  Dump the MASR's.

  dumpMASR(_tigEntry, _tigLen, _tigMax, _currentVersion);

  //  Close the current version; we'll reopen on demand.

  for (uint32 p=0; p<MAX_PART; p++) {
    if (_dataFile[_currentVersion][p].FP) {
      errno = 0;
      fclose(_dataFile[_currentVersion][p].FP);
      if (errno)
        fprintf(stderr, "tgStore::nextVersion()-- Failed to close '%s': %s\n", _name, strerror(errno)), exit(1);

      _dataFile[_currentVersion][p].FP    = NULL;
      _dataFile[_currentVersion][p].atEOF = false;
    }
  }

  //  Bump to the next version.

  //fprintf(stderr, "tgStore::tgStore()-- moving from version %d to version %d; _tigLen %d _ctgLen %d\n",
  //        _currentVersion, _currentVersion+1, _tigLen, _ctgLen);

  _currentVersion++;

  //  Remove any existing files at that version level.

  purgeCurrentVersion();
}



void
tgStore::writeToPartitioned(uint32 *partMap_, uint32 partMapLen_) {

  assert(_writable == true);          //  Must be writable to write!
  assert(_inplace == false);
  assert(_partMap == NULL);
  assert(_part == 0);

  //  The _dataFile for the unpartitioned data cannot have data in it.

  if (_dataFile[_currentVersion][0].FP != NULL)
    fprintf(stderr, "tgStore::writeToPartitioned()-- ERROR!  There is already data in the unpartitioned store, cannot convert to a partitioned store.\n");
  assert(_dataFile[_currentVersion][0].FP == NULL);

  _partMap    = partMap_;
  _partMapLen = partMapLen_;
}



void
tgStore::writeTigToDisk(tgTig *tig, tgStoreEntry *te) {

  //fprintf(stderr, "tgStore::writeTigToDisk()-- write ma "F_S32" in store version "F_U64" partition "F_U64" at file position "F_U64"\n", tig->_tigID, maRecord->svID, maRecord->ptID, maRecord->fileOffset);

  FILE *FP = openDB(te->svID, te->ptID);

  //  The atEOF flag allows us to skip a seek when we're already (supposed) to be at the EOF.  This
  //  (hopefully) fixes a problem on one system where the seek() was placing the FP just before EOF
  //  (almost like the last block of data wasn't being flushed), and the tell() would then place the
  //  next tig in the middle of the previous one.
  //
  //  It also should (greatly) improve performance over NFS, espeically during BOG and CNS.  Both of
  //  these only write data, so no repositioning of the stream is needed.
  //
  if (_dataFile[te->svID][te->ptID].atEOF == false) {
    AS_UTL_fseek(FP, 0, SEEK_END);
    _dataFile[te->svID][te->ptID].atEOF = true;
  }

  te->flushNeeded = 0;
  te->fileOffset  = AS_UTL_ftell(FP);

  tig->saveToStream(FP);
}



void
tgStore::insertTig(tgTig *tig, bool keepInCache) {

  //  Check that the components do not exceed the bound.
  //
  if (tig->_gappedBases > 0) {
    uint32  len = tig->_gappedLen;
    uint32  neg = 0;
    uint32  pos = 0;

    for (uint32 i=0; i<tig->_childrenLen; i++) {
      tgPosition *read = tig->_children + i;

      if ((read->_bgn < 0) || (read->_end < 0))
        fprintf(stderr, "tgStore::insertTig()-- ERROR: tig %d read %d at (%d,%d) has negative position\n",
                tig->_tigID, read->_objID, read->_bgn, read->_end), neg++;
      if (read->_bgn < 0)  read->_bgn = 0;
      if (read->_end < 0)  read->_end = 0;

      if ((read->_bgn > len) || (read->_end > len))
        fprintf(stderr, "tgStore::insertTig()-- ERROR: tig %d read %d at (%d,%d) exceeded multialign length %d\n",
                tig->_tigID, read->_objID, read->_bgn, read->_end, len), pos++;
      if (read->_bgn > len)  read->_bgn = len;
      if (read->_end > len)  read->_end = len;
    }

    if (neg + pos > 0) {
      tig->dumpLayout(stderr);
      fprintf(stderr, "tgStore::insertTig()-- ERROR: tig %d has invalid layout, exceeds bounds of consensus sequence (length %d) -- neg=%d pos=%d.\n",
              tig->_tigID, len, neg, pos);
    }

    assert(neg == 0);
    assert(pos == 0);
  }


  if (tig->_tigID == UINT32_MAX) {
    tig->_tigID = _tigLen;
    _newTigs  = true;

    //  Make sure that the UTG line (if present) agrees with our new tigID
    //if (GetNumIntUnitigPoss(tig->u_list) == 1)
    //  GetIntUnitigPos(tig->u_list, 0)->_id = tig->_tigID;

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

    delete [] _tigEntry;
    delete [] _tigCache;

    _tigEntry = nr;
    _tigCache  = nc;
  }

  _tigLen = MAX(_tigLen, tig->_tigID + 1);

  assert(_tigEntry->isDeleted == 0);

#if 0
  if (_tigEntry->svID > 0)
    fprintf(stderr, "tgStore::InsetTig()--  Moving tig tigID %d from svID %d to svID %d\n",
            tig->_tigID, _tigEntry->svID, _currentVersion);
#endif

  _tigEntry->unusedFlags     = 0;
  _tigEntry->flushNeeded     = 1;  //  Mark as needing a flush by default
  _tigEntry->isPresent       = 1;
  _tigEntry->isDeleted       = 0;
  _tigEntry->ptID            = 0;
  _tigEntry->svID            = _currentVersion;

  _tigEntry->tigRecord        = *tig;

  //  Decide on which partition to write to.  If any of the partMaps are set, use that.  Else, use
  //  the _part we are restricted to.

  //  ???

  if (_partMap) {
    if (_partMapLen <= tig->_tigID)
      fprintf(stderr, "tgStore::insertTig()--  ERROR!  Attempt to insert a tig (id=%d) into a store partitioned only for %d tigs.\n",
              tig->_tigID, _partMapLen);
    if (_partMap[tig->_tigID] == 0)
      fprintf(stderr, "tgStore::insertTig()--  ERROR!  tig %d is partitioned to partition 0.\n", tig->_tigID);
    if (_partMap[tig->_tigID] == 0xffffffff)
      fprintf(stderr, "tgStore::insertTig()--  ERROR!  tig %d is partitioned, but not in a partition.\n", tig->_tigID);
    assert(_partMapLen > tig->_tigID);
    assert(_partMap[tig->_tigID] != 0);
    assert(_partMap[tig->_tigID] != 0xffffffff);
  }

  _tigEntry->ptID = (_partMap) ? _partMap[tig->_tigID] : _part;

  //  Write to disk RIGHT NOW unless we're keeping it in cache.  If it is written, the flushNeeded
  //  flag is cleared.
  //
  if (keepInCache == false)
    writeTigToDisk(tig, _tigEntry);

  //  If the cache is different from this tig, delete the cache.  Not sure why this happens --
  //  did we copy a tig, muck with it, and then want to replace the one in the store?
  //
  if (_tigCache[tig->_tigID] != tig) {
    delete _tigCache[tig->_tigID];
    _tigCache[tig->_tigID] = NULL;
  }

  //  Cache it if requested, otherwise clear the cache.
  //
  _tigCache[tig->_tigID] = (keepInCache) ? tig : NULL;
}



void
tgStore::deleteTig(uint32 tigID) {
  assert(tigID >= 0);
  assert(tigID <  _tigLen);

  flushDisk(tigID);

  assert(_tigEntry[tigID].flushNeeded == 0);

  assert(_tigEntry[tigID].isPresent == 1);
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

  //  This is...and is not...an error.  It does indicate something didn't go according to plan, like
  //  loading a tig that doesn't exist (that should be caught by the above 'tigID < _tigLen'
  //  assert).  Unfortunately, the 'isPresent' flag is set to false for all tigs not in our
  //  partition, and so we MUST return NULL here.
  //
  if (_tigEntry[tigID].isPresent == 0)
    return(NULL);

  if (_tigEntry[tigID].isDeleted == 1)
    return(NULL);

  //  If we're loading from a specific partition, and it isn't the correct one, we can't load, so return NULL.
  if ((_part > 0) && (_tigEntry[tigID].ptID != _part))
    return(NULL);

  //  Otherwise, we can load something.

  if (_tigCache[tigID] == NULL) {
    FILE *FP = openDB(_tigEntry[tigID].svID, _tigEntry[tigID].ptID);

    //  Since the tig isn't in the cache, it had better NOT be marked as needing to be flushed!
    assert(_tigEntry[tigID].flushNeeded == 0);

    //  Seek to the correct position, and reset the atEOF to indicate we're (with high probability)
    //  not at EOF anymore.
    if (_dataFile[_tigEntry[tigID].svID][_tigEntry[tigID].ptID].atEOF == true) {
      fflush(FP);
      _dataFile[_tigEntry[tigID].svID][_tigEntry[tigID].ptID].atEOF = false;
    }

    AS_UTL_fseek(FP, _tigEntry[tigID].fileOffset, SEEK_SET);

    _tigCache[tigID] = new tgTig;
    _tigCache[tigID]->loadFromStream(FP);

    if (_tigCache[tigID] == NULL)
      fprintf(stderr,"loadTig()-- FAILED for tig "F_S32" in file "F_U64" at offset "F_U64"\n",
              tigID, _tigEntry[tigID].svID, _tigEntry[tigID].fileOffset);
    assert(_tigCache[tigID] != NULL);

    //  ALWAYS assume the incore record is more up to date
    *_tigCache[tigID] = _tigEntry[tigID].tigRecord;

    //  Since we just loaded, no flush is needed.
    _tigEntry[tigID].flushNeeded = 0;
  }

  return(_tigCache[tigID]);
}



void
tgStore::unloadTig(uint32 tigID, bool discard) {

  if (discard)
    _tigEntry[tigID].flushNeeded = 0;

  flushDisk(tigID);

  assert(_tigEntry[tigID].flushNeeded == 0);

  delete _tigCache[tigID];
  _tigCache[tigID] = NULL;
}


void
tgStore::copyTig(uint32 tigID, tgTig *tigcopy) {

  assert(tigID >= 0);
  assert(tigID <  _tigLen);
  assert(_tigEntry[tigID].isPresent == 1);

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

  FILE *FP = openDB(_tigEntry[tigID].svID, _tigEntry[tigID].ptID);

  //  Seek to the correct position, and reset the atEOF to indicate we're (with high probability)
  //  not at EOF anymore.
  
  if (_dataFile[_tigEntry[tigID].svID][_tigEntry[tigID].ptID].atEOF == true) {
    fflush(FP);
    _dataFile[_tigEntry[tigID].svID][_tigEntry[tigID].ptID].atEOF = false;
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



void
tgStore::dumpMASRfile(char *name, tgStoreEntry *R, uint32 L, uint32 M, uint32 part) {

  errno = 0;
  FILE *F = fopen(name, "w");
  if (errno)
    fprintf(stderr, "tgStore::dumpMASRfile()-- Failed to create '%s': %s\n", name, strerror(errno)), exit(1);

  AS_UTL_safeWrite(F, &MASRmagic,   "MASRmagic",   sizeof(uint32), 1);
  AS_UTL_safeWrite(F, &MASRversion, "MASRversion", sizeof(uint32), 1);
  AS_UTL_safeWrite(F, &L,           "MASRtotal",   sizeof(uint32), 1);

  if (_part != 0) {
    uint32        indxLen = 0;
    uint32        masrLen = 0;

    uint32         *indx    = (uint32      *)safe_malloc(sizeof(uint32)      * L);
    tgStoreEntry  *masr    = (tgStoreEntry *)safe_malloc(sizeof(tgStoreEntry) * L);

    //  Copy all the metadata for this partition into our buffer...
    //
    //  ...and, if this is the first partition, copy all the deleted stuff here too.  No client will
    //  ever get this (see loadTig()), but we should haul the crud along with us, I guess.
    //
    for (uint32 i=0; i<L; i++) {
      if ((R[i].ptID == part) ||
          ((_part == 1) && (R[i].isDeleted == 1))) {
        indx[indxLen++] = i;
        masr[masrLen++] = R[i];
      }
    }

    //fprintf(stderr, "tgStore::dumpMASRfile()-- Writing '%s' partitioned (indxLen=%d masrLen=%d).\n", name, indxLen, masrLen);

    AS_UTL_safeWrite(F, &indxLen, "MASRindxLen", sizeof(uint32),      1);
    AS_UTL_safeWrite(F, &masrLen, "MASRlen",     sizeof(uint32),      1);

    AS_UTL_safeWrite(F,  indx,    "MASRindx",    sizeof(uint32),      indxLen);
    AS_UTL_safeWrite(F,  masr,    "MASR",        sizeof(tgStoreEntry), masrLen);

    safe_free(indx);
    safe_free(masr);

  } else {
    uint32  indxLen = 0;

    //fprintf(stderr, "tgStore::dumpMASRfile()-- Writing '%s' unpartitioned (indxLen=%d masrLen=%d).\n", name, indxLen, L);

    AS_UTL_safeWrite(F, &indxLen, "MASRindexLen", sizeof(uint32),      1);
    AS_UTL_safeWrite(F, &L,       "MASRlen",      sizeof(uint32),      1);
    AS_UTL_safeWrite(F,  R,       "MASR",         sizeof(tgStoreEntry), L);
  }

  fclose(F);
}


bool
tgStore::loadMASRfile(char *name, tgStoreEntry* R, uint32 L, uint32 M, uint32 part, bool onlyThisV) {
  uint32        MASRmagicInFile   = 0;
  uint32        MASRversionInFile = 0;
  uint32        MASRtotalInFile   = 0;

  uint32        indxLen = 0;
  uint32        masrLen = 0;

  if (AS_UTL_fileExists(name, false, false) == false)
    return(false);

  errno = 0;
  FILE *F = fopen(name, "r");
  if (errno)
    fprintf(stderr, "tgStore::loadMASRfile()-- Failed to open '%s': %s\n", name, strerror(errno)), exit(1);

  AS_UTL_safeRead(F, &MASRmagicInFile,   "MASRmagic",   sizeof(uint32), 1);
  AS_UTL_safeRead(F, &MASRversionInFile, "MASRversion", sizeof(uint32), 1);
  AS_UTL_safeRead(F, &MASRtotalInFile,   "MASRtotal",   sizeof(uint32), 1);
  AS_UTL_safeRead(F, &indxLen,           "MASRindxLen", sizeof(uint32), 1);
  AS_UTL_safeRead(F, &masrLen,           "MASRmasrLen", sizeof(uint32), 1);

  if (MASRmagicInFile != MASRmagic) {
    fprintf(stderr, "tgStore::loadMASRfile()-- Failed to open '%s': magic number mismatch; file=0x%08x code=0x%08x\n",
            name, MASRmagicInFile, MASRmagic);
    exit(1);
  }

  if (MASRversionInFile != MASRversion) {
    fprintf(stderr, "tgStore::loadMASRfile()-- Failed to open '%s': version number mismatch; file=%d code=%d\n",
            name, MASRversionInFile, MASRversion);
    exit(1);
  }

  //  Check we're consistent.
  if (L < MASRtotalInFile)
    fprintf(stderr, "tgStore::loadMASRfile()-- '%s' has more tigs ("F_U32") than expected ("F_U32").\n",
            name, MASRtotalInFile, L), exit(1);

  if (indxLen > 0) {
    //  A partitioned file.  Load the index, load the data, then copy the data into the real
    //  array.

    uint32        *indx = (uint32      *)safe_malloc(sizeof(uint32)      * indxLen);
    tgStoreEntry *masr = (tgStoreEntry *)safe_malloc(sizeof(tgStoreEntry) * masrLen);

    AS_UTL_safeRead(F, indx, "indx", sizeof(uint32),      indxLen);
    AS_UTL_safeRead(F, masr, "masr", sizeof(tgStoreEntry), masrLen);

    for (uint32 i=0; i<indxLen; i++)
      if ((onlyThisV == false) || (masr[i].svID == _currentVersion))
        R[indx[i]] = masr[i];

    safe_free(indx);
    safe_free(masr);

  } else {
    //  Not a partitioned file.  Can directly load the data.
    AS_UTL_safeRead(F,  R, "MASR", sizeof(tgStoreEntry), masrLen);
  }

  fclose(F);

  return(true);
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
tgStore::dumpMASR(tgStoreEntry* &R, uint32& L, uint32& M, uint32 V) {

  //  We can't both read from a specific partition and create a partitioned store.
  assert((_part == 0) || (_partMap == NULL));

  //  Not partitioned; dump a single file.  If we started off un partitioned, but then became
  //  partitioned, _dataFile[V][0] will exist, and we still want to dump.  This lets us open the
  //  store unpartitioned.
  //
  if ((_part == 0) && (_partMap == NULL)) {
    sprintf(_name, "%s/seqDB.v%03d.tig", _path, V);
    dumpMASRfile(_name, R, L, M, 0);
    return;
  }

  //  Partitioned, but we are restricted to one partition.
  //
  if (_part != 0) {
    sprintf(_name, "%s/seqDB.v%03d.p%03d.tig", _path, V, _part);
    dumpMASRfile(_name, R, L, M, _part);
    return;
  }

  //  Writing to partitions, dump ALL partition files.  The unpartitioned entry saves all the
  //  pointers for tigs that are not partitioned: deleted tigs, and surrogate instances are the only
  //  two that I know of.  See comments on loading.
  //
  sprintf(_name, "%s/seqDB.v%03d.tig", _path, V);
  dumpMASRfile(_name, R, L, M, 0);

  for (uint32 p=1; p<MAX_PART; p++) {
    if (_dataFile[_currentVersion][p].FP) {
      sprintf(_name, "%s/seqDB.v%03d.p%03d.tig", _path, V, p);
      dumpMASRfile(_name, R, L, M, p);
    }
  }
}


void
tgStore::loadMASR(tgStoreEntry* &R, uint32& L, uint32& M, uint32 V, bool onlyThisV) {

  //  We can't both read from a specific partition and create a partitioned store.
  assert((_part == 0) || (_partMap == NULL));

  //  Allocate space for the data.  We don't always need to do this, sometimes
  //  we're called to update the data with newer data.
  //
  if (R == NULL) {

    //  Search for the first *.utg or *.ctg file, ask it how many tigs are in the store.
    //
    for (int32 i=V; i>0; i--) {
      sprintf(_name, "%s/seqDB.v%03d.tig", _path, i);
      L = numTigsInMASRfile(_name);
      if (L > 0)
        break;
    }

    //  Allocate space for the data.
    //
    M = L + 1024;
    R = (tgStoreEntry *)safe_calloc(M, sizeof(tgStoreEntry));
  }

  //  If partitioned, load just the partition requested.
  //
  if (_part != 0) {
    sprintf(_name, "%s/seqDB.v%03d.p%03d.tig", _path, V, _part);
    loadMASRfile(_name, R, L, M, _part, onlyThisV);
    return;
  }

  //  Try to load the full unpartitioned data.  There are two use cases here:
  //    1) we are loading an unpartitioned store;     seqDB.v###.typ
  //    2) we are loading ALL of a parttioned store;  seqDB.v###.p###.typ
  //
  //  This is the first case.
  //
  if ((_part == 0) && (_partMap == NULL)) {
    sprintf(_name, "%s/seqDB.v%03d.tig", _path, V);
    if (loadMASRfile(_name, R, L, M, 0, onlyThisV))
      return;
  }

  //  Must be case 2.  Load all partitions.
  //
  //  Unfortunately, the unpartitioned data is NOT propagated through versions.  For example, cgw
  //  creates version 25.  It is partitioned.  Partitions 1, 2 and 3 exist, as does the
  //  unpartitioned seqDB.v025.ctg file.  Consensus reads partitioned 1-3 for input, and writes data
  //  to version 26, partitions 1-3.  seqDB.v025.ctg is never read by consensus (nor should it be!),
  //  and it never appears in version 26.  Terminator is told to output version 26, but it needs to
  //  load in seqDB.v025.ctg first to get all the unpartitioned tigs.
  //
  for (int32 i=V; i>0; i--) {
    sprintf(_name, "%s/seqDB.v%03d.tig", _path, i);
    if (loadMASRfile(_name, R, L, M, 0, onlyThisV))
      break;
  }

  sprintf(_name, "%s/seqDB.v%03d.p001.tig", _path, V);
  for (uint32 p=1; loadMASRfile(_name, R, L, M, p, onlyThisV); p++)
    sprintf(_name, "%s/seqDB.v%03d.p%03d.tig", _path, V, p+1);
}




FILE *
tgStore::openDB(uint32 version, uint32 partition) {

  if (_dataFile[version][partition].FP)
    return(_dataFile[version][partition].FP);

  //  If partition is zero, open the unpartitioned store.

  if (_part == 0)
    sprintf(_name, "%s/seqDB.v%03d.dat", _path, version);
  else
    sprintf(_name, "%s/seqDB.v%03d.p%03d.dat", _path, version, partition);

  //  Try again: On some large assemblies (or misconfigured partitioning) we exhaust the number of
  //  open files.  This will close the earlier versions (repoened on demand) when we fail to open a
  //  file.
  //
  //  This came into existence after BPW forgot to pass the desired partition size from runCA to 
  //  CGW, and ended up with an assembly with too many partitions.  Gatekeeper couldn't open enough
  //  files for the gkpStore partitioning.  CGW was called again to repartition, but it too opened
  //  too many files.
  //
  int  tryAgain = 1;
 doTryAgain:

  errno = 0;

  //  If version is the _currentVersion, open for writing if allowed.
  //
  //  "a+" technically writes (always) to the end of file, but this hasn't been tested.

  if ((_writable) && (version == _currentVersion)) {
    _dataFile[version][partition].FP = fopen(_name, "a+");
    _dataFile[version][partition].atEOF = false;
  } else {
    _dataFile[version][partition].FP = fopen(_name, "r");
    _dataFile[version][partition].atEOF = false;
  }

  if ((errno) && (tryAgain)) {
    tryAgain = 0;

    fprintf(stderr, "tgStore::openDB()-- Failed to open '%s': %s\n", _name, strerror(errno));
    fprintf(stderr, "tgStore::openDB()-- Trying again.\n");

    for (uint32 v=0; v<_currentVersion; v++)
      for (uint32 p=0; p<MAX_PART; p++)
        if (_dataFile[v][p].FP) {
          fclose(_dataFile[v][p].FP);
          _dataFile[v][p].FP = NULL;
        }
    goto doTryAgain;
  }

  if (errno)
    fprintf(stderr, "tgStore::openDB()-- Failed to open '%s': %s\n", _name, strerror(errno)), exit(1);

  return(_dataFile[version][partition].FP);
}
