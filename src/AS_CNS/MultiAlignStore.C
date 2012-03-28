
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

static const char *rcsid = "$Id: MultiAlignStore.C,v 1.21 2012-03-28 06:11:25 brianwalenz Exp $";

#include "AS_global.h"
#include "AS_UTL_fileIO.h"
#include "MultiAlignStore.h"

uint32  MASRmagic   = 0x5253414d;  //  'MASR', as a big endian integer
uint32  MASRversion = 1;

#define MAX_VERS   1024
#define MAX_PART   1024

void
MultiAlignStore::init(const char *path_, uint32 version_, bool writable_, bool inplace_, bool append_) {

  strcpy(path, path_);

  writable          = writable_;
  inplace           = inplace_;
  append            = append_;

  newTigs           = false;

  currentVersion    = version_;
  originalVersion   = version_;

  unitigPartMap     = NULL;
  contigPartMap     = NULL;

  unitigPart        = 0;
  contigPart        = 0;

  utgMax            = 0;
  utgLen            = 0;
  utgRecord         = NULL;
  utgCache          = NULL;

  ctgMax            = 0;
  ctgLen            = 0;
  ctgRecord         = NULL;
  ctgCache          = NULL;

  //  Could use sysconf(_SC_OPEN_MAX) too.  Should make this dynamic?
  //
  dataFile          = (dataFileT **)safe_calloc(MAX_VERS,            sizeof(dataFileT *));
  dataFile[0]       = (dataFileT  *)safe_calloc(MAX_VERS * MAX_PART, sizeof(dataFileT));

  for (uint32 i=1; i<MAX_VERS; i++)
    dataFile[i] = dataFile[0] + i * MAX_PART;
}


MultiAlignStore::MultiAlignStore(const char *path_) {
  init(path_, 1, true, false, false);
  AS_UTL_mkdir(path);     //  Create the directory, if needed.
  purgeCurrentVersion();  //  Purge any data there currently.
}


MultiAlignStore::MultiAlignStore(const char *path_,
                                 uint32      version_,
                                 uint32      unitigPartition_,
                                 uint32      contigPartition_,
                                 bool        writable_,
                                 bool        inplace_,
                                 bool        append_) {

  init(path_, version_, writable_, inplace_, append_);

  unitigPart = unitigPartition_;
  contigPart = contigPartition_;

  if ((unitigPart != 0) && (contigPart != 0))
    fprintf(stderr, "MultiAlignStore::MultiAlignStore()-- ERROR, cannot set both unitigPart and contigPart.\n"), exit(1);

  if ((writable == false) && (inplace == true))
    fprintf(stderr, "MultiAlignStore::MultiAlignStore()-- ERROR, cannot operate inplace unless writable.\n"), exit(1);

  if ((writable == false) && (append == true))
    fprintf(stderr, "MultiAlignStore::MultiAlignStore()-- ERROR, cannot append unless writable.\n"), exit(1);

  if ((inplace == true) && (append == true))
    fprintf(stderr, "MultiAlignStore::MultiAlignStore()-- ERROR, cannot both append and be inplace.\n"), exit(1);

  //  Load the MultiAlignRs for the current version.

  loadMASR(utgRecord, utgLen, utgMax, currentVersion, TRUE,  FALSE);
  loadMASR(ctgRecord, ctgLen, ctgMax, currentVersion, FALSE, FALSE);

  if ((utgLen == 0) && (ctgLen == 0) && (append == false)) {
    fprintf(stderr, "MultiAlignStore::MultiAlignStore()-- ERROR, didn't find any unitigs or contigs in the store.\n");
    fprintf(stderr, "MultiAlignStore::MultiAlignStore()--        asked for store '%s', correct?\n", path);
    fprintf(stderr, "MultiAlignStore::MultiAlignStore()--        asked for version '%d', correct?\n", originalVersion);
    fprintf(stderr, "MultiAlignStore::MultiAlignStore()--        asked for partition unitig=%d contig=%d, correct?\n", unitigPart, contigPart);
    fprintf(stderr, "MultiAlignStore::MultiAlignStore()--        asked for writable=%d inplace=%d append=%d, correct?\n", writable, inplace, append);
    exit(1);
  }

  //  Allocate the cache to the proper size

  utgCache = (MultiAlignT **)safe_calloc(utgMax, sizeof(MultiAlignT *));
  ctgCache = (MultiAlignT **)safe_calloc(ctgMax, sizeof(MultiAlignT *));

  //  Open the next version for writing, and remove what is currently there.

  if ((writable == true) && (inplace == false) && (append == false)) {
    currentVersion++;

    purgeCurrentVersion();
  }

  //  Open the next version for writing, and keep the data that is currently there.
  //
  //  Special cases for partitioned stores (consensus).  Do NOT load the next version of the
  //  unpartitioned data.  For contig consensus, we want to update the contig data with the stuff in
  //  the next version.  This lets us resume failed ctgcns runs.  However, we do not want any of the
  //  unitig data in the next version.  Unitigs shouldn't change, and likely cannot change since the
  //  store is currently partitioned on contigs.
  //
  if (append == true) {
    currentVersion++;

    if (contigPart == 0)
      //  Contigs are partitioned, so do NOT load the next version of unitigs.  They shouldn't be changing.
      loadMASR(utgRecord, utgLen, utgMax, currentVersion, TRUE,  TRUE);

    if (unitigPart == 0)
      //  Unitigs are partitioned, so do NOT load the next version of contigs.  They shouldn't be changing.
      loadMASR(ctgRecord, ctgLen, ctgMax, currentVersion, FALSE, TRUE);
  }

  if ((utgLen == 0) && (ctgLen == 0)) {
    fprintf(stderr, "MultiAlignStore::MultiAlignStore()-- ERROR, didn't find any unitigs or contigs in the store.  Correct version?\n");
    exit(1);
  }
}


MultiAlignStore::~MultiAlignStore() {

  flushCache();

  //  If writable, and we aren't working on contig-partitioned data, write the unitig data.  If we
  //  are working on contig-partitioned data, when we load the store the next time, we'll fall back
  //  to the previous version for unitigs.

  if ((writable) && (contigPart == 0))
    dumpMASR(utgRecord, utgLen, utgMax, currentVersion, TRUE);

  if ((writable) && (unitigPart == 0))
    dumpMASR(ctgRecord, ctgLen, ctgMax, currentVersion, FALSE);

  //  If new tigs were added, AND we are partitioned, update the global
  //  partition.
  //
  if ((newTigs) && (unitigPart != 0)) {
    assert(unitigPartMap == NULL);

    unitigPart = 0;  //  To load the unpartitioned MASR
    loadMASR(utgRecord, utgLen, utgMax, originalVersion, TRUE, FALSE);
    dumpMASR(utgRecord, utgLen, utgMax, originalVersion, TRUE);
  }

  if ((newTigs) && (contigPart != 0)) {
    assert(contigPartMap == NULL);

    contigPart = 0;
    loadMASR(ctgRecord, ctgLen, ctgMax, originalVersion, FALSE, FALSE);
    dumpMASR(ctgRecord, ctgLen, ctgMax, originalVersion, FALSE);
  }

  safe_free(utgRecord);
  safe_free(utgCache);

  safe_free(ctgRecord);
  safe_free(ctgCache);

  for (uint32 v=0; v<MAX_VERS; v++)
    for (uint32 p=0; p<MAX_PART; p++)
      if (dataFile[v][p].FP)
        fclose(dataFile[v][p].FP);

  safe_free(dataFile[0]);
  safe_free(dataFile);
}



void
MultiAlignStore::purgeCurrentVersion(void) {
  char   name[FILENAME_MAX];
  uint32 part = 1;

  sprintf(name, "%s/seqDB.v%03d.dat", path, currentVersion);   AS_UTL_unlink(name);
  sprintf(name, "%s/seqDB.v%03d.ctg", path, currentVersion);   AS_UTL_unlink(name);
  sprintf(name, "%s/seqDB.v%03d.utg", path, currentVersion);   AS_UTL_unlink(name);

  sprintf(name, "%s/seqDB.v%03d.p001.dat", path, currentVersion);
  while (AS_UTL_unlink(name)) {
    sprintf(name, "%s/seqDB.v%03d.p%03d.ctg", path, currentVersion, part);   AS_UTL_unlink(name);
    sprintf(name, "%s/seqDB.v%03d.p%03d.utg", path, currentVersion, part);   AS_UTL_unlink(name);
    sprintf(name, "%s/seqDB.v%03d.p%03d.dat", path, currentVersion, ++part);
  }
}



void
MultiAlignStore::nextVersion(void) {

  assert(writable == true);
  assert(inplace == false);
  assert(unitigPartMap == NULL);
  assert(contigPartMap == NULL);
  assert(unitigPart == 0);
  assert(contigPart == 0);

  //  Dump the MASR's.

  dumpMASR(utgRecord, utgLen, utgMax, currentVersion, TRUE);
  dumpMASR(ctgRecord, ctgLen, ctgMax, currentVersion, FALSE);

  //  Close the current version; we'll reopen on demand.

  for (uint32 p=0; p<MAX_PART; p++) {
    if (dataFile[currentVersion][p].FP) {
      errno = 0;
      fclose(dataFile[currentVersion][p].FP);
      if (errno)
        fprintf(stderr, "MultiAlignStore::nextVersion()-- Failed to close '%s': %s\n", name, strerror(errno)), exit(1);

      dataFile[currentVersion][p].FP    = NULL;
      dataFile[currentVersion][p].atEOF = false;
    }
  }

  //  Bump to the next version.

  //fprintf(stderr, "MultiAlignStore()-- moving from version %d to version %d; utgLen %d ctgLen %d\n",
  //        currentVersion, currentVersion+1, utgLen, ctgLen);

  currentVersion++;

  //  Remove any existing files at that version level.

  purgeCurrentVersion();
}



void
MultiAlignStore::writeToPartitioned(uint32 *unitigPartMap_, uint32 unitigPartMapLen_,
                                    uint32 *contigPartMap_, uint32 contigPartMapLen_) {

  assert(writable == true);          //  Must be writable to write!
  assert(inplace == false);
  assert(unitigPartMap == NULL);
  assert(contigPartMap == NULL);
  assert(unitigPart == 0);
  assert(contigPart == 0);

  //  The dataFile for the unpartitioned data cannot have data in it.

  if (dataFile[currentVersion][0].FP != NULL)
    fprintf(stderr, "MultiAlignStore::writeToPartitioned()-- ERROR!  There is already data in the unpartitioned store, cannot convert to a partitioned store.\n");
  assert(dataFile[currentVersion][0].FP == NULL);

  //  We cannot handle partitioning both unitigs and contigs.

  if ((unitigPartMap_ != NULL) && (contigPartMap_ != NULL))
    fprintf(stderr, "MultiAlignStore::writeToPartitiond()--  ERROR!  Attempting to partition both unitigs and contigs.\n");
  assert((unitigPartMap_ == NULL) || (contigPartMap_ == NULL));

  unitigPartMap    = unitigPartMap_;
  unitigPartMapLen = unitigPartMapLen_;

  contigPartMap    = contigPartMap_;
  contigPartMapLen = contigPartMapLen_;
}



void
MultiAlignStore::insertMultiAlign(MultiAlignT *ma, bool isUnitig, bool keepInCache) {

  assert(ma->data.unitig_status      != 0);
  assert(ma->data.unitig_unique_rept != 0);
  assert(ma->data.contig_status      != 0);

  //  We can add unitigs only if we are not contig-partitioned, likewise for contigs.
  if ((isUnitig == 1) && (contigPart != 0))
    fprintf(stderr, "MultiAlignStore::insertMultiAlign()-- ERROR: attempted to add a unitig when the store is partitioned on contigs.\n");
  assert((contigPart == 0) || (isUnitig == 0));

  if ((isUnitig == 0) && (unitigPart != 0))
    fprintf(stderr, "MultiAlignStore::insertMultiAlign()-- ERROR: attempted to add a contig when the store is partitioned on unitigs.\n");
  assert((unitigPart == 0) || (isUnitig == 1));

  if (isUnitig == 1) {
    if (ma->maID < 0) {
      ma->maID = utgLen;
      newTigs  = true;
      assert(contigPart       == 0);
      assert(contigPartMap    == NULL);
      assert(contigPartMapLen == 0);
      fprintf(stderr, "insertMultiAlign()-- Added new unitig %d\n", ma->maID);
    }

    if (utgMax <= ma->maID) {
      utgMax = (utgMax == 0) ? (65536) : (2 * utgMax);
      assert(ma->maID < utgMax);

      utgRecord = (MultiAlignR  *)safe_realloc(utgRecord, utgMax * sizeof(MultiAlignR));
      utgCache  = (MultiAlignT **)safe_realloc(utgCache,  utgMax * sizeof(MultiAlignT *));

      memset(utgRecord + utgLen, 0, sizeof(MultiAlignR)   * (utgMax - utgLen));
      memset(utgCache  + utgLen, 0, sizeof(MultiAlignT *) * (utgMax - utgLen));
    }

    utgLen = MAX(utgLen, ma->maID + 1);
  }

  if (isUnitig == 0) {
    if (ma->maID < 0) {
      ma->maID = ctgLen;
      newTigs  = true;
      assert(unitigPart       == 0);
      assert(unitigPartMap    == NULL);
      assert(unitigPartMapLen == 0);
      fprintf(stderr, "insertMultiAlign()-- Added new contig %d\n", ma->maID);
    }

    if (ctgMax <= ma->maID) {
      ctgMax = (ctgMax == 0) ? (65536) : (2 * ctgMax);
      assert(ma->maID < ctgMax);

      ctgRecord = (MultiAlignR  *)safe_realloc(ctgRecord, ctgMax * sizeof(MultiAlignR));
      ctgCache  = (MultiAlignT **)safe_realloc(ctgCache,  ctgMax * sizeof(MultiAlignT *));

      memset(ctgRecord + ctgLen, 0, sizeof(MultiAlignR)   * (ctgMax - ctgLen));
      memset(ctgCache  + ctgLen, 0, sizeof(MultiAlignT *) * (ctgMax - ctgLen));
    }

    ctgLen = MAX(ctgLen, ma->maID + 1);
  }

  MultiAlignR  *maRecord = (isUnitig) ? (utgRecord + ma->maID) : (ctgRecord + ma->maID);

  assert(maRecord->isDeleted == 0);

  maRecord->unusedFlags     = 0;
  maRecord->isPresent       = 1;
  maRecord->isDeleted       = 0;
  maRecord->ptID            = 0;
  maRecord->svID            = currentVersion;

  //  ALWAYS update mad on insert.
  ma->data.num_frags        = GetNumIntMultiPoss(ma->f_list);
  ma->data.num_unitigs      = GetNumIntUnitigPoss(ma->u_list);
  maRecord->mad             = ma->data;

  //  Decide on which partition to write to.
  //
  //  If any of the partMaps are set, use that.  Else, use the unitigPart/contigPart we are
  //  restricted to.
  //
  //  Do NOT allow tigs of the other flavor to be written when we are partitioned.  This might break
  //  future ctgcns changes, where we want to recompute the unitig as we recompute the contig.
  //
  if (isUnitig == true) {
    if ((contigPart > 0) || (contigPartMap != NULL))
      fprintf(stderr, "MultiAlignStore::insertMultiAlign()--  ERROR!  Attempt to insert a unitig into a store parititioned on contigs.\n");
    assert(contigPart    == 0);
    assert(contigPartMap == NULL);
    if (unitigPartMap) {
      if (unitigPartMapLen <= ma->maID)
        fprintf(stderr, "MultiAlignStore::insertMultiAlign()--  ERROR!  Attempt to insert a unitig (id=%d) into a store partitioned only for %d unitigs.\n",
                ma->maID, unitigPartMapLen);
      if (unitigPartMap[ma->maID] == 0)
        fprintf(stderr, "MultiAlignStore::insertMultiAlign()--  ERROR!  Unitig %d is partitioned to partition 0.\n", ma->maID);
      if (unitigPartMap[ma->maID] == 0xffffffff)
        fprintf(stderr, "MultiAlignStore::insertMultiAlign()--  ERROR!  Unitig %d is partitioned, but not in a partition.\n", ma->maID);
      assert(unitigPartMapLen > ma->maID);
      assert(unitigPartMap[ma->maID] != 0);
      assert(unitigPartMap[ma->maID] != 0xffffffff);
    }
    maRecord->ptID = (unitigPartMap) ? unitigPartMap[ma->maID] : unitigPart;
  }

  if (isUnitig == false) {
    if ((unitigPart > 0) || (unitigPartMap != NULL))
      fprintf(stderr, "MultiAlignStore::insertMultiAlign()--  ERROR!  Attempt to insert a contig into a store parititioned on unitigs.\n");
    assert(unitigPart    == 0);
    assert(unitigPartMap == NULL);
    if (contigPartMap) {
      if (contigPartMapLen <= ma->maID)
        fprintf(stderr, "MultiAlignStore::insertMultiAlign()--  ERROR!  Attempt to insert a contig (id=%d) into a store partitioned only for %d contigs.\n",
                ma->maID, contigPartMapLen);
      if (contigPartMap[ma->maID] == 0)
        fprintf(stderr, "MultiAlignStore::insertMultiAlign()--  ERROR!  Contig %d is partitioned to partition 0.\n", ma->maID);
      if (contigPartMap[ma->maID] == 0xffffffff)
        fprintf(stderr, "MultiAlignStore::insertMultiAlign()--  ERROR!  Contig %d is partitioned, but not in a partition.\n", ma->maID);
      assert(contigPartMapLen > ma->maID);
      assert(contigPartMap[ma->maID] != 0);
      assert(contigPartMap[ma->maID] != 0xffffffff);
    }
    maRecord->ptID = (contigPartMap) ? contigPartMap[ma->maID] : contigPart;
  }


  FILE *FP = openDB(maRecord->svID, maRecord->ptID);

  //  The atEOF flag allows us to skip a seek when we're already (supposed) to be at the EOF.  This
  //  (hopefully) fixes a problem on one system where the seek() was placing the FP just before EOF
  //  (almost like the last block of data wasn't being flushed), and the tell() would then place the
  //  next tig in the middle of the previous one.
  //
  //  It also should (greatly) improve performance over NFS, espeically during BOG and CNS.  Both of
  //  these only write data, so no repositioning of the stream is needed.
  //
  if (dataFile[maRecord->svID][maRecord->ptID].atEOF == false) {
    AS_UTL_fseek(FP, 0, SEEK_END);
    dataFile[maRecord->svID][maRecord->ptID].atEOF = true;
  }

  maRecord->fileOffset = AS_UTL_ftell(FP);

  SaveMultiAlignTToStream(ma, FP);

  //fprintf(stderr, "MA %d in store version %d partition %d at file position %d\n", ma->maID, maRecord->svID, maRecord->ptID, maRecord->fileOffset);

  MultiAlignT   **maCache  = (isUnitig) ? (utgCache) : (ctgCache);

  //  If we want to save this in the cache, delete whatever is there (unless it is us) and save it.
  //  The store now owns this MutliAlignT.
  //
  //  If not, delete whatever is there (unless it is us), but do NOT delete us.  The store doesn NOT
  //  own this MultiAlignT, even if it was initially supplied by loadMultiAlign().
  //
  if (maCache[ma->maID] != ma)
    DeleteMultiAlignT(maCache[ma->maID]);

  maCache[ma->maID] = (keepInCache) ? ma : NULL;
}



void
MultiAlignStore::deleteMultiAlign(int32 maID, bool isUnitig) {
  uint32                  maLen    = (isUnitig) ? utgLen    : ctgLen;
  MultiAlignR            *maRecord = (isUnitig) ? utgRecord : ctgRecord;
  MultiAlignT           **maCache  = (isUnitig) ? utgCache  : ctgCache;

  assert(maID >= 0);
  assert(maID < (int32)maLen);

  assert(maRecord[maID].isPresent == 1);
  assert(maRecord[maID].isDeleted == 0);

  maRecord[maID].isDeleted = 1;

  DeleteMultiAlignT(maCache[maID]);

  maCache[maID] = NULL;
}



MultiAlignT *
MultiAlignStore::loadMultiAlign(int32 maID, bool isUnitig) {
  uint32                  maLen    = (isUnitig) ? utgLen    : ctgLen;
  MultiAlignR            *maRecord = (isUnitig) ? utgRecord : ctgRecord;
  MultiAlignT           **maCache  = (isUnitig) ? utgCache  : ctgCache;
  bool                    cantLoad = true;

  assert(maID >= 0);
  assert(maID < (int32)maLen);

  //  This is...and is not...an error.  It does indicate something didn't go according to plan, like
  //  loading a unitig that doesn't exist (that should be caught by the above 'maID < maLen'
  //  assert).  Unfortunately, the 'isPresent' flag is set to FALSE for all tigs not in our
  //  partition, and so we MUST return NULL here.
  //
  if (maRecord[maID].isPresent == 0)
    return(NULL);

  if (maRecord[maID].isDeleted == 1)
    return(NULL);

  //  If we're not reading a specific partition, load.
  if ((isUnitig == true)  && (unitigPart == 0))
    goto canLoad;
  if ((isUnitig == false) && (contigPart == 0))
    goto canLoad;

  //  If we're loading a unitig, and we're limited to a contig partition, load.
  if ((isUnitig == true) && (contigPart != 0))
    goto canLoad;

  //  If we're loading from a specific partition, and it's the correct one, load.
  if ((isUnitig == true)  && (maRecord[maID].ptID == unitigPart))
    goto canLoad;
  if ((isUnitig == false) && (maRecord[maID].ptID == contigPart))
    goto canLoad;

  //  Otherwise, can't load.
  return(NULL);

 canLoad:
  if (maCache[maID] == NULL) {
    FILE *FP = openDB(maRecord[maID].svID, maRecord[maID].ptID);

    //  Seek to the correct position, and reset the atEOF to indicate we're (with high probability)
    //  not at EOF anymore.
    if (dataFile[maRecord[maID].svID][maRecord[maID].ptID].atEOF == true) {
      fflush(FP);
      dataFile[maRecord[maID].svID][maRecord[maID].ptID].atEOF = false;
    }

    AS_UTL_fseek(FP, maRecord[maID].fileOffset, SEEK_SET);

    maCache[maID] = LoadMultiAlignTFromStream(FP);

    if (maCache[maID] == NULL)
      fprintf(stderr,"MultiAlignStore::loadMultiAlign()-- FAILED for %s "F_S32" in file "F_U64" at offset "F_U64"\n",
              (isUnitig ? "Unitig" : "Contig"), maID, maRecord[maID].svID, maRecord[maID].fileOffset);
    assert(maCache[maID] != NULL);

    //  ALWAYS assume the incore mad is more up to date
    maCache[maID]->data = maRecord[maID].mad;
  }

  return(maCache[maID]);
}



void
MultiAlignStore::unloadMultiAlign(int32 maID, bool isUnitig) {

  if (isUnitig)
    DeleteMultiAlignT(utgCache[maID]);
  else
    DeleteMultiAlignT(ctgCache[maID]);
}


void
MultiAlignStore::copyMultiAlign(int32 maID, bool isUnitig, MultiAlignT *macopy) {
  uint32                  maLen    = (isUnitig) ? utgLen    : ctgLen;
  MultiAlignR            *maRecord = (isUnitig) ? utgRecord : ctgRecord;
  MultiAlignT           **maCache  = (isUnitig) ? utgCache  : ctgCache;

  assert(maID >= 0);
  assert(maID < (int32)maLen);
  assert(maRecord[maID].isPresent == 1);

  if (maRecord[maID].isDeleted) {
    ClearMultiAlignT(macopy);
    return;
  }

  if (maCache[maID]) {
    CopyMultiAlignT(macopy, maCache[maID]);
  } else {
    FILE *FP = openDB(maRecord[maID].svID, maRecord[maID].ptID);

    //  Seek to the correct position, and reset the atEOF to indicate we're (with high probability)
    //  not at EOF anymore.

    if (dataFile[maRecord[maID].svID][maRecord[maID].ptID].atEOF == true) {
      fflush(FP);
      dataFile[maRecord[maID].svID][maRecord[maID].ptID].atEOF = false;
    }

    AS_UTL_fseek(FP, maRecord[maID].fileOffset, SEEK_SET);

    ReLoadMultiAlignTFromStream(FP, macopy);
  }

  //  ALWAYS assume the incore mad is more up to date
  macopy->data = maRecord[maID].mad;
}



void
MultiAlignStore::flushCache(void) {
  uint32  numUFlushed = 0;
  uint32  numCFlushed = 0;

  for (uint32 i=0; i<utgLen; i++) {
    if (utgCache[i]) {
      DeleteMultiAlignT(utgCache[i]);
      utgCache[i] = NULL;
      numUFlushed++;
    }
  }

  for (uint32 i=0; i<ctgLen; i++) {
    if (ctgCache[i]) {
      DeleteMultiAlignT(ctgCache[i]);
      ctgCache[i] = NULL;
      numCFlushed++;
    }
  }

  //fprintf(stderr, "MultiAlignStore::flushCache()-- flushed %u unitigs and %u contigs.\n",
  //        numUFlushed, numCFlushed);

#if 0
  for (uint32 p=0; p<MAX_PART; p++)
    if (dataFile[currentVersion][p].FP)
      fflush(dataFile[currentVersion][p].FP);
#endif
}











void
MultiAlignStore::dumpMASRfile(char *name, MultiAlignR *R, uint32 L, uint32 M, uint32 part) {
  errno = 0;
  FILE *F = fopen(name, "w");
  if (errno)
    fprintf(stderr, "MultiAlignStore::dumpMASRfile()-- Failed to create '%s': %s\n", name, strerror(errno)), exit(1);

  AS_UTL_safeWrite(F, &MASRmagic,   "MASRmagic",   sizeof(uint32), 1);
  AS_UTL_safeWrite(F, &MASRversion, "MASRversion", sizeof(uint32), 1);
  AS_UTL_safeWrite(F, &L,           "MASRtotal",   sizeof(uint32), 1);

  if (part != 0) {
    uint32        indxLen = 0;
    uint32        masrLen = 0;

    uint32       *indx    = (uint32      *)safe_malloc(sizeof(uint32)      * L);
    MultiAlignR  *masr    = (MultiAlignR *)safe_malloc(sizeof(MultiAlignR) * L);

    //  Copy all the metadata for this partition into our buffer...
    //
    //  ...and, if this is the first partition, copy all the deleted stuff here too.  No client will
    //  ever get this (see loadMultiAlign()), but we should haul the crud along with us, I guess.
    //
    for (uint32 i=0; i<L; i++) {
      if ((R[i].ptID == part) ||
          ((part == 1) && (R[i].isDeleted == 1))) {
        indx[indxLen++] = i;
        masr[masrLen++] = R[i];
      }
    }

    //fprintf(stderr, "MultiAlignStore::dumpMASRfile()-- Writing '%s' partitioned (indxLen=%d masrLen=%d).\n", name, indxLen, masrLen);

    AS_UTL_safeWrite(F, &indxLen, "MASRindxLen", sizeof(uint32),      1);
    AS_UTL_safeWrite(F, &masrLen, "MASRlen",     sizeof(uint32),      1);

    AS_UTL_safeWrite(F,  indx,    "MASRindx",    sizeof(uint32),      indxLen);
    AS_UTL_safeWrite(F,  masr,    "MASR",        sizeof(MultiAlignR), masrLen);

    safe_free(indx);
    safe_free(masr);

  } else {
    uint32  indxLen = 0;

    //fprintf(stderr, "MultiAlignStore::dumpMASRfile()-- Writing '%s' unpartitioned (indxLen=%d masrLen=%d).\n", name, indxLen, L);

    AS_UTL_safeWrite(F, &indxLen, "MASRindexLen", sizeof(uint32),      1);
    AS_UTL_safeWrite(F, &L,       "MASRlen",      sizeof(uint32),      1);
    AS_UTL_safeWrite(F,  R,       "MASR",         sizeof(MultiAlignR), L);
  }

  fclose(F);
}


bool
MultiAlignStore::loadMASRfile(char *name, MultiAlignR* R, uint32 L, uint32 M, uint32 part, bool onlyThisV) {
  uint32        MASRmagicInFile   = 0;
  uint32        MASRversionInFile = 0;
  uint32        MASRtotalInFile   = 0;

  uint32        indxLen = 0;
  uint32        masrLen = 0;

  uint32       *indx    = NULL;
  MultiAlignR  *masr    = NULL;

  if (AS_UTL_fileExists(name, false, false) == false)
    return(false);

  errno = 0;
  FILE *F = fopen(name, "r");
  if (errno)
    fprintf(stderr, "MultiAlignStore::loadMASRfile()-- Failed to open '%s': %s\n", name, strerror(errno)), exit(1);

  AS_UTL_safeRead(F, &MASRmagicInFile,   "MASRmagic",   sizeof(uint32), 1);
  AS_UTL_safeRead(F, &MASRversionInFile, "MASRversion", sizeof(uint32), 1);
  AS_UTL_safeRead(F, &MASRtotalInFile,   "MASRtotal",   sizeof(uint32), 1);
  AS_UTL_safeRead(F, &indxLen,           "MASRindxLen", sizeof(uint32), 1);
  AS_UTL_safeRead(F, &masrLen,           "MASRmasrLen", sizeof(uint32), 1);

  if (MASRmagicInFile != MASRmagic) {
    fprintf(stderr, "MultiAlignStore::loadMASRfile()-- Failed to open '%s': magic number mismatch; file=0x%08x code=0x%08x\n",
            name, MASRmagicInFile, MASRmagic);
    exit(1);
  }

  if (MASRversionInFile != MASRversion) {
    fprintf(stderr, "MultiAlignStore::loadMASRfile()-- Failed to open '%s': version number mismatch; file=%d code=%d\n",
            name, MASRversionInFile, MASRversion);
    exit(1);
  }

  //  Check we're consistent.
  if (L < MASRtotalInFile)
    fprintf(stderr, "MultiAlignStore::loadMASRfile()-- '%s' has more tigs ("F_U32") than expected ("F_U32").\n",
            name, MASRtotalInFile, L), exit(1);

  if (indxLen > 0) {
    //  A partitioned file.  Load the index, load the data, then copy the data into the real
    //  array.
    indx = (uint32      *)safe_malloc(sizeof(uint32)      * indxLen);
    masr = (MultiAlignR *)safe_malloc(sizeof(MultiAlignR) * masrLen);

    AS_UTL_safeRead(F, indx, "indx", sizeof(uint32),      indxLen);
    AS_UTL_safeRead(F, masr, "masr", sizeof(MultiAlignR), masrLen);

    for (uint32 i=0; i<indxLen; i++)
      if ((onlyThisV == FALSE) || (masr[i].svID == currentVersion))
        R[indx[i]] = masr[i];

    safe_free(indx);
    safe_free(masr);

  } else {
    //  Not a partitioned file.  Can directly load the data.
    AS_UTL_safeRead(F,  R, "MASR", sizeof(MultiAlignR), masrLen);
  }

  fclose(F);

  return(true);
}


uint32
MultiAlignStore::numTigsInMASRfile(char *name) {
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
    fprintf(stderr, "MultiAlignStore::numTigsInMASRfile()-- Failed to open '%s': %s\n", name, strerror(errno)), exit(1);

  AS_UTL_safeRead(F, &MASRmagicInFile,   "MASRmagic",   sizeof(uint32), 1);
  AS_UTL_safeRead(F, &MASRversionInFile, "MASRversion", sizeof(uint32), 1);
  AS_UTL_safeRead(F, &MASRtotalInFile,   "MASRtotal",   sizeof(uint32), 1);

  fclose(F);

  if (MASRmagicInFile != MASRmagic) {
    fprintf(stderr, "MultiAlignStore::numTigsInMASRfile()-- Failed to open '%s': magic number mismatch; file=0x%08x code=0x%08x\n",
            name, MASRmagicInFile, MASRmagic);
    exit(1);
  }

  if (MASRversionInFile != MASRversion) {
    fprintf(stderr, "MultiAlignStore::numTigsInMASRfile()-- Failed to open '%s': version number mismatch; file=%d code=%d\n",
            name, MASRversionInFile, MASRversion);
    exit(1);
  }

  return(MASRtotalInFile);
}

void
MultiAlignStore::dumpMASR(MultiAlignR* &R, uint32& L, uint32& M, uint32 V, bool isUnitig) {
  uint32   part = (isUnitig) ? unitigPart    : contigPart;
  uint32  *pmap = (isUnitig) ? unitigPartMap : contigPartMap;

  //  We can't both read from a specific partition and create a partitioned store.
  assert((part == 0) || (pmap == NULL));

  //  Only one partitioning allowed, either unitigs or contigs.
  assert((unitigPart == 0) || (contigPart == 0));

  //  If partitioned, and told to write the other type, do nothing.
  //
  if (((isUnitig == TRUE)  && (contigPart != 0)) ||
      ((isUnitig == FALSE) && (unitigPart != 0)))
    return;

  //  Not partitioned; dump a single file.  If we started off un partitioned, but then became
  //  partitioned, dataFile[V][0] will exist, and we still want to dump.  This lets us open the
  //  store unpartitioned.
  //
  if ((part == 0) && (pmap == NULL)) {
    sprintf(name, "%s/seqDB.v%03d.%s", path, V, (isUnitig) ? "utg" : "ctg");
    dumpMASRfile(name, R, L, M, 0);
    return;
  }

  //  Partitioned, but we are restricted to one partition.
  //
  if (part != 0) {
    sprintf(name, "%s/seqDB.v%03d.p%03d.%s", path, V, unitigPart + contigPart, (isUnitig) ? "utg" : "ctg");
    dumpMASRfile(name, R, L, M, unitigPart + contigPart);
    return;
  }

  //  Writing to partitions, dump ALL partition files.  The unpartitioned entry saves all the
  //  pointers for tigs that are not partitioned: deleted tigs, and surrogate instances are the only
  //  two that I know of.  See comments on loading.
  //
  sprintf(name, "%s/seqDB.v%03d.%s", path, V, (isUnitig) ? "utg" : "ctg");
  dumpMASRfile(name, R, L, M, 0);

  for (uint32 p=1; p<MAX_PART; p++) {
    if (dataFile[currentVersion][p].FP) {
      sprintf(name, "%s/seqDB.v%03d.p%03d.%s", path, V, p, (isUnitig) ? "utg" : "ctg");
      dumpMASRfile(name, R, L, M, p);
    }
  }
}


void
MultiAlignStore::loadMASR(MultiAlignR* &R, uint32& L, uint32& M, uint32 V, bool isUnitig, bool onlyThisV) {
  uint32   part = (isUnitig) ? unitigPart    : contigPart;
  uint32  *pmap = (isUnitig) ? unitigPartMap : contigPartMap;

  //  We can't both read from a specific partition and create a partitioned store.
  assert((part == 0) || (pmap == NULL));

  //  Only one partitioning allowed, either unitigs or contigs.
  assert((unitigPart == 0) || (contigPart == 0));

  //  Allocate space for the data.  We don't always need to do this, sometimes
  //  we're called to update the data with newer data.
  //
  if (R == NULL) {

    //  Search for the first *.utg or *.ctg file, ask it how many tigs are in the store.
    //
    for (int32 i=V; i>0; i--) {
      sprintf(name, "%s/seqDB.v%03d.%s", path, i, (isUnitig) ? "utg" : "ctg");
      L = numTigsInMASRfile(name);
      if (L > 0)
        break;
    }

    //  Allocate space for the data.
    //
    M = L + 1024;
    R = (MultiAlignR *)safe_calloc(M, sizeof(MultiAlignR));
  }

  //  If partitioned, load just the partition requested.  In contig consensus, we'll request to load
  //  contig partition N, and unitig partition 0, which will load just the contigs we care about,
  //  and all the unitigs.
  //
  if (part != 0) {
    sprintf(name, "%s/seqDB.v%03d.p%03d.%s", path, V, part, (isUnitig) ? "utg" : "ctg");
    loadMASRfile(name, R, L, M, part, onlyThisV);
    return;
  }

  //  Try to load the full unpartitioned data.  There are two use cases here:
  //    1) we are loading an unpartitioned store;     seqDB.v###.typ
  //    2) we are loading ALL of a parttioned store;  seqDB.v###.p###.typ
  //
  //  This is the first case.
  //
  if ((part == 0) && (pmap == NULL)) {
    sprintf(name, "%s/seqDB.v%03d.%s", path, V, (isUnitig) ? "utg" : "ctg");
    if (loadMASRfile(name, R, L, M, 0, onlyThisV))
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
    sprintf(name, "%s/seqDB.v%03d.%s", path, i, (isUnitig) ? "utg" : "ctg");
    if (loadMASRfile(name, R, L, M, 0, onlyThisV))
      break;
  }

  sprintf(name, "%s/seqDB.v%03d.p001.%s", path, V, (isUnitig) ? "utg" : "ctg");
  for (uint32 p=1; loadMASRfile(name, R, L, M, p, onlyThisV); p++)
    sprintf(name, "%s/seqDB.v%03d.p%03d.%s", path, V, p+1, (isUnitig) ? "utg" : "ctg");
}




FILE *
MultiAlignStore::openDB(uint32 version, uint32 partition) {

  if (dataFile[version][partition].FP)
    return(dataFile[version][partition].FP);

  //  If partition is zero, open the unpartitioned store.

  if (partition == 0)
    sprintf(name, "%s/seqDB.v%03d.dat", path, version);
  else
    sprintf(name, "%s/seqDB.v%03d.p%03d.dat", path, version, partition);

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

  //  If version is the currentVersion, open for writing if allowed.
  //
  //  "a+" technically writes (always) to the end of file, but this hasn't been tested.

  if ((writable) && (version == currentVersion)) {
    dataFile[version][partition].FP = fopen(name, "a+");
    dataFile[version][partition].atEOF = false;
  } else {
    dataFile[version][partition].FP = fopen(name, "r");
    dataFile[version][partition].atEOF = false;
  }

  if ((errno) && (tryAgain)) {
    tryAgain = 0;

    fprintf(stderr, "MultiAlignStore::openDB()-- Failed to open '%s': %s\n", name, strerror(errno));
    fprintf(stderr, "MultiAlignStore::openDB()-- Trying again.\n");

    for (uint32 v=0; v<currentVersion; v++)
      for (uint32 p=0; p<MAX_PART; p++)
        if (dataFile[v][p].FP) {
          fclose(dataFile[v][p].FP);
          dataFile[v][p].FP = NULL;
        }
    goto doTryAgain;
  }

  if (errno)
    fprintf(stderr, "MultiAlignStore::openDB()-- Failed to open '%s': %s\n", name, strerror(errno)), exit(1);

  return(dataFile[version][partition].FP);
}



void
MultiAlignStore::dumpMultiAlignR(int32 maID, bool isUnitig) {
  MultiAlignR  *maRecord = (isUnitig) ? utgRecord : ctgRecord;

  fprintf(stdout, "maRecord.isPresent   = "F_U64"\n", maRecord[maID].isPresent);
  fprintf(stdout, "maRecord.isDeleted   = "F_U64"\n", maRecord[maID].isDeleted);
  fprintf(stdout, "maRecord.ptID        = "F_U64"\n", maRecord[maID].ptID);
  fprintf(stdout, "maRecord.svID        = "F_U64"\n", maRecord[maID].svID);
  if (isUnitig) { fprintf(stdout, "maRecord.fileOffset  = "F_U64"\n", maRecord[maID].fileOffset); }
}



void
MultiAlignStore::dumpMultiAlignRTable(bool isUnitig) {
  MultiAlignR  *maRecord = (isUnitig) ? utgRecord : ctgRecord;
  uint32        len      = (isUnitig) ? utgLen    : ctgLen;

  fprintf(stdout, "maID\tisPresent\tisDeleted\tptID\tsvID\tfileOffset");
  if (isUnitig) fprintf(stdout, "\tcovStat");
  fprintf(stdout, "\n");

  for (uint32 i=0; i<len; i++) {
    fprintf(stdout, ""F_U32"\t", i);
    fprintf(stdout, ""F_U64"\t", maRecord[i].isPresent);
    fprintf(stdout, ""F_U64"\t", maRecord[i].isDeleted);
    fprintf(stdout, ""F_U64"\t", maRecord[i].ptID);
    fprintf(stdout, ""F_U64"\t", maRecord[i].svID);
    fprintf(stdout, ""F_U64, maRecord[i].fileOffset);
    if (isUnitig) fprintf(stdout, "\t%d",  getUnitigCoverageStat(i));
    fprintf(stdout, "\n");
  }
}
