
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

static const char *rcsid = "$Id: MultiAlignStore.C,v 1.7 2009-12-03 01:08:29 brianwalenz Exp $";

#include "AS_global.h"
#include "AS_UTL_fileIO.h"
#include "MultiAlignStore.h"

#define MAX_VERS   1024
#define MAX_PART   1024

void
MultiAlignStore::init(const char *path_, uint32 version_, bool writable_, bool inplace_) {

  strcpy(path, path_);

  writable          = writable_;
  creating          = false;
  inplace           = inplace_;

  currentVersion    = version_;

  unitigPartMap     = NULL;
  contigPartMap     = NULL;

  unitigPart        = 0;
  contigPart        = 0;

  utgMax            = 65536;
  utgLen            = 0;
  utgRecord         = (MultiAlignR  *)safe_calloc(utgMax, sizeof(MultiAlignR));
  utgCache          = (MultiAlignT **)safe_calloc(utgMax, sizeof(MultiAlignT *));

  ctgMax            = 65536;
  ctgLen            = 0;
  ctgRecord         = (MultiAlignR  *)safe_calloc(ctgMax, sizeof(MultiAlignR));
  ctgCache          = (MultiAlignT **)safe_calloc(ctgMax, sizeof(MultiAlignT *));

  //  Could use sysconf(_SC_OPEN_MAX) too.  Should make this dynamic?
  //
  dataFile          = (dataFileT **)safe_calloc(MAX_VERS,            sizeof(dataFileT *));
  dataFile[0]       = (dataFileT  *)safe_calloc(MAX_VERS * MAX_PART, sizeof(dataFileT));

  for (uint32 i=1; i<MAX_VERS; i++)
    dataFile[i] = dataFile[0] + i * MAX_PART;
}


MultiAlignStore::MultiAlignStore(const char *path_) {
  init(path_, 1, true, false);

  creating = true;

  AS_UTL_mkdir(path);
}


MultiAlignStore::MultiAlignStore(const char *path_,
                                 uint32      version_,
                                 uint32      unitigPartition_,
                                 uint32      contigPartition_,
                                 bool        writable_,
                                 bool        inplace_) {

  init(path_, version_, writable_, inplace_);

  unitigPart = unitigPartition_;
  contigPart = contigPartition_;

  if ((unitigPart != 0) && (contigPart != 0))
    fprintf(stderr, "MultiAlignStore::MultiAlignStore()-- ERROR, cannot set both unitigPart and contigPart.\n");
  assert((unitigPart == 0) || (contigPart == 0));

  //  Load the MultiAlignRs for the current version.

  loadMASR(utgRecord, "utg", utgLen, utgMax, currentVersion);
  loadMASR(ctgRecord, "ctg", ctgLen, ctgMax, currentVersion);

  //  Reallocate the cache to the proper size

  safe_free(utgCache);
  safe_free(ctgCache);

  utgCache = (MultiAlignT **)safe_calloc(utgMax, sizeof(MultiAlignT *));
  ctgCache = (MultiAlignT **)safe_calloc(ctgMax, sizeof(MultiAlignT *));

  //  Open the next version for writing.

  if ((writable == true) && (inplace == false))
    currentVersion++;
}


MultiAlignStore::~MultiAlignStore() {

  flushCache();

  if (writable) {
    dumpMASR(utgRecord, "utg", utgLen, utgMax, currentVersion);
    dumpMASR(ctgRecord, "ctg", ctgLen, ctgMax, currentVersion);
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
MultiAlignStore::nextVersion(void) {

  assert(writable == true);
  assert(inplace == false);
  assert(unitigPartMap == NULL);
  assert(contigPartMap == NULL);
  assert(unitigPart == 0);
  assert(contigPart == 0);

  //  Dump the MASR's.

  dumpMASR(utgRecord, "utg", utgLen, utgMax, currentVersion);
  dumpMASR(ctgRecord, "ctg", ctgLen, ctgMax, currentVersion);

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

  currentVersion++;

  //  Remove any existing files at that version level.

  {
    char  name[FILENAME_MAX];
    int32 part = 1;

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
}



void
MultiAlignStore::writeToPartitioned(uint32 *unitigPartMap_, uint32 *contigPartMap_) {

  assert(writable == true);          //  Must be writable to write!
  assert(inplace == false);
  assert(unitigPartMap == NULL);
  assert(contigPartMap == NULL);
  assert(unitigPart == 0);
  assert(contigPart == 0);

  //  We allow both of these to be set, though it usually makes little sense to do that.

  unitigPartMap = unitigPartMap_;
  contigPartMap = contigPartMap_;
}



void
MultiAlignStore::insertMultiAlign(MultiAlignT *ma, bool isUnitig, bool keepInCache) {

  assert(ma->maID >= 0);

  assert(ma->data.unitig_status      != 0);
  assert(ma->data.unitig_unique_rept != 0);
  assert(ma->data.contig_status      != 0);

  if (isUnitig == 1) {
    if (utgMax <= ma->maID) {
      utgMax *= 2;

      utgRecord = (MultiAlignR  *)safe_realloc(utgRecord, utgMax * sizeof(MultiAlignR));
      utgCache  = (MultiAlignT **)safe_realloc(utgCache,  utgMax * sizeof(MultiAlignT *));

      memset(utgRecord + utgLen, 0, sizeof(MultiAlignR)   * (utgMax - utgLen));
      memset(utgCache  + utgLen, 0, sizeof(MultiAlignT *) * (utgMax - utgLen));
    }

    utgLen = MAX(utgLen, ma->maID + 1);
  }

  if (isUnitig == 0) {
    if (ctgMax <= ma->maID) {
      ctgMax *= 2;

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
  //    If any of the partMaps are set, use that.
  //
  //    If unitigPart is set, the unitig must have come from that partition.
  //
  //    If contigPart is set, the contig must have come from that partition.  Unitigs,
  //    however, need to be reset to the same partition.
  //
  if (isUnitig == true) {
    maRecord->ptID = unitigPart;

    if (unitigPartMap)
      maRecord->ptID = unitigPartMap[ma->maID];

    if (contigPart)
      maRecord->ptID = contigPart;
  }

  if (isUnitig == false) {
    maRecord->ptID = contigPart;

    if (contigPartMap)
      maRecord->ptID = contigPartMap[ma->maID];
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
  int32                   maLen    = (isUnitig) ? utgLen    : ctgLen;
  MultiAlignR            *maRecord = (isUnitig) ? utgRecord : ctgRecord;
  MultiAlignT           **maCache  = (isUnitig) ? utgCache  : ctgCache;

  assert(maID >= 0);
  assert(maID < maLen);

  assert(maRecord[maID].isPresent == 1);
  assert(maRecord[maID].isDeleted == 0);

  maRecord[maID].isDeleted = 1;

  DeleteMultiAlignT(maCache[maID]);

  maCache[maID] = NULL;
}



MultiAlignT *
MultiAlignStore::loadMultiAlign(int32 maID, bool isUnitig) {
  int32                   maLen    = (isUnitig) ? utgLen    : ctgLen;
  MultiAlignR            *maRecord = (isUnitig) ? utgRecord : ctgRecord;
  MultiAlignT           **maCache  = (isUnitig) ? utgCache  : ctgCache;
  bool                    cantLoad = true;

  assert(maID < maLen);
  assert(maRecord[maID].isPresent == 1);

  if (maRecord[maID].isDeleted)
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
MultiAlignStore::copyMultiAlign(int32 maID, bool isUnitig, MultiAlignT *macopy) {
  int32                   maLen    = (isUnitig) ? utgLen    : ctgLen;
  MultiAlignR            *maRecord = (isUnitig) ? utgRecord : ctgRecord;
  MultiAlignT           **maCache  = (isUnitig) ? utgCache  : ctgCache;

  assert(maID < maLen);
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

  for (int32 i=0; i<utgLen; i++)
    DeleteMultiAlignT(utgCache[i]);

  memset(utgCache, 0, utgMax * sizeof(MultiAlignT *));

  for (int32 i=0; i<ctgLen; i++)
    DeleteMultiAlignT(ctgCache[i]);

  memset(ctgCache, 0, ctgMax * sizeof(MultiAlignT *));

  for (uint32 p=0; p<MAX_PART; p++)
    if (dataFile[currentVersion][p].FP)
      fflush(dataFile[currentVersion][p].FP);
}












void
MultiAlignStore::dumpMASRfile(char *name, MultiAlignR *R, int32 L, int32 M) {
  errno = 0;
  FILE *F = fopen(name, "w");
  if (errno)
    fprintf(stderr, "MultiAlignStore::dumpMASRfile()-- Failed to create '%s': %s\n", name, strerror(errno)), exit(1);

  AS_UTL_safeWrite(F, &L, "MASRlen", sizeof(uint32),      1);
  AS_UTL_safeWrite(F,  R, "MASR",    sizeof(MultiAlignR), L);

  fclose(F);
}


void
MultiAlignStore::dumpMASR(MultiAlignR* &R, char *T, int32& L, int32& M, uint32 V) {

  //  Not partitioned; dump a single file.  If we started off un partitioned, but then became
  //  partitioned, dataFile[V][0] will exist, and we still want to dump.  This lets us open the
  //  store unpartitioned.
  //
  if (((unitigPart == 0) && (unitigPartMap == NULL) &&
       (contigPart == 0) && (contigPartMap == NULL)) ||
      (dataFile[V][0].FP)) {
    sprintf(name, "%s/seqDB.v%03d.%s", path, V, T);
    dumpMASRfile(name, R, L, M);
  }

  if ((unitigPart == 0) && (unitigPartMap == NULL) &&
      (contigPart == 0) && (contigPartMap == NULL))
    return;

  //  Partitioned, but we are restricted to one partition.
  //
  if ((unitigPart != 0) || (contigPart != 0)) {
    assert((unitigPart == 0) || (contigPart == 0));  //  Checked on object creation, too
    sprintf(name, "%s/seqDB.v%03d.p%03d.%s", path, V, unitigPart + contigPart, T);
    dumpMASRfile(name, R, L, M);
    return;
  }

  //  Writing to partitions, dump ALL partition files (at this point, they're all the same).
  //
  assert((unitigPartMap != NULL) || (contigPartMap != NULL));

  for (uint32 p=1; p<MAX_PART; p++) {
    if (dataFile[currentVersion][p].FP) {
      sprintf(name, "%s/seqDB.v%03d.p%03d.%s", path, V, p, T);
      dumpMASRfile(name, R, L, M);
    }
  }
}





void
MultiAlignStore::loadMASRfile(char *name, MultiAlignR* &R, int32& L, int32& M) {
  errno = 0;
  FILE *F = fopen(name, "r");
  if (errno)
    fprintf(stderr, "MultiAlignStore::loadMASRfile()-- Failed to open '%s': %s\n", name, strerror(errno)), exit(1);

  AS_UTL_safeRead(F, &L, "MASRlen", sizeof(uint32), 1);

  if (M < L) {
    M = L;
    safe_free(R);
    R = (MultiAlignR *)safe_calloc(L, sizeof(MultiAlignR));
  }

  AS_UTL_safeRead(F,  R, "MASR",    sizeof(MultiAlignR), L);

  fclose(F);
}


void
MultiAlignStore::loadMASR(MultiAlignR* &R, char *T, int32& L, int32& M, uint32 V) {
  MultiAlignR  *Rp = NULL;
  int32         Lp = 0;
  int32         Mp = 0;

  sprintf(name, "%s/seqDB.v%03d.%s", path, V, T);

  if (AS_UTL_fileExists(name, false, false)) {
    loadMASRfile(name, R, L, M);
    return;
  }

  assert(creating == false);

  //  Load the first partition; nice side effect is that this allocates the R array to the correct
  //  size.

  sprintf(name, "%s/seqDB.v%03d.p001.%s", path, V, T);

  if (AS_UTL_fileExists(name, false, false)) {
    loadMASRfile(name, R, L, M);
  } else {
    fprintf(stderr, "MultiAlignStore::loadMASR()-- ERROR:  Didn't find version %d, partition %d\n", V, 1);
    fprintf(stderr, "  %s/seqDB.v%03d.%s        not present\n", path, V, T);
    fprintf(stderr, "  %s/seqDB.v%03d.p001.%s   not present\n", path, V, T);
    exit(1);
  }

  //  Load and merge the rest of the partitions.

  sprintf(name, "%s/seqDB.v%03d.p002.%s", path, V, T);

  for (uint32 p=2; AS_UTL_fileExists(name, false, false); p++) {
    loadMASRfile(name, Rp, Lp, Mp);

    for (int32 i=0; i<Lp; i++)
      if (Rp[i].ptID == p)
        R[i] = Rp[i];

    sprintf(name, "%s/seqDB.v%03d.p%03d.%s", path, V, p+1, T);
  }

  safe_free(Rp);
}




FILE *
MultiAlignStore::openDB(uint32 version, uint32 partition) {

  if (dataFile[version][partition].FP)
    return(dataFile[version][partition].FP);

  //  If partition is zero, open the unpartitioned store.

  if (partition == 0) {
    sprintf(name, "%s/seqDB.v%03d.dat", path, version);
  } else {
    sprintf(name, "%s/seqDB.v%03d.p%03d.dat", path, version, partition);
  }

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

  if ((inplace) && (version == currentVersion)) {
    dataFile[version][partition].FP = fopen(name, "a+");
    dataFile[version][partition].atEOF = false;
  } else if ((writable) && (version == currentVersion)) {
    dataFile[version][partition].FP = fopen(name, "w+");
    dataFile[version][partition].atEOF = true;
  } else {
    dataFile[version][partition].FP = fopen(name, "r");
    dataFile[version][partition].atEOF = false;
  }

  if ((errno) && (tryAgain)) {
    tryAgain = 0;

    fprintf(stderr, "MultiAlignStore::openDB()-- Failed to open '%s': %s\n", name, strerror(errno));
    fprintf(stderr, "MultiAlignStore::openDB()-- Trying again.\n");

    for (int v=0; v<currentVersion; v++)
      for (int p=0; p<MAX_PART; p++)
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

  fprintf(stdout, "maRecord.isPresent   = %d\n", maRecord[maID].isPresent);
  fprintf(stdout, "maRecord.isDeleted   = %d\n", maRecord[maID].isDeleted);
  fprintf(stdout, "maRecord.ptID        = %d\n", maRecord[maID].ptID);
  fprintf(stdout, "maRecord.svID        = %d\n", maRecord[maID].svID);
  fprintf(stdout, "maRecord.fileOffset  = %d\n", maRecord[maID].fileOffset);
}



void
MultiAlignStore::dumpMultiAlignRTable(bool isUnitig) {
  MultiAlignR  *maRecord = (isUnitig) ? utgRecord : ctgRecord;
  int32         len      = (isUnitig) ? utgLen    : ctgLen;

  fprintf(stdout, "maID\tisPresent\tisDeleted\tptID\tsvID\tfileOffset\n");

  for (int32 i=0; i<len; i++) {
    fprintf(stdout, "%d\t", i);
    fprintf(stdout, "%d\t", maRecord[i].isPresent);
    fprintf(stdout, "%d\t", maRecord[i].isDeleted);
    fprintf(stdout, "%d\t", maRecord[i].ptID);
    fprintf(stdout, "%d\t", maRecord[i].svID);
    fprintf(stdout, "%d\n", maRecord[i].fileOffset);
  }
}
