
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

static char *rcsid = "$Id: AS_PER_gkStore.C,v 1.12 2009-10-26 13:20:26 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <sys/stat.h>

#include "AS_global.h"
#include "AS_PER_genericStore.h"
#include "AS_PER_gkpStore.h"
#include "AS_PER_encodeSequenceQuality.h"
#include "AS_UTL_fileIO.h"



gkStore::gkStore() {
  uid      = createStringStore(NULL, "uid");
  STRtoUID = CreateStringHashTable_AS();
}



gkStore::gkStore(const char *path, int partition) {
  char       name[FILENAME_MAX];

  gkStore_clear();

  strcpy(storePath, path);

  assert(partition > 0);
  partnum = partition;

  sprintf(name,"%s/fpk.%03d", storePath, partnum);
  partfpk = createIndexStore(name, "partfpk", sizeof(gkPackedFragment), 1);

  sprintf(name,"%s/fnm.%03d", storePath, partnum);
  partfnm = createIndexStore(name, "partfnm", sizeof(gkNormalFragment), 1);
  sprintf(name,"%s/qnm.%03d", storePath, partnum);
  partqnm = createStringStore(name, "partqnm");

  sprintf(name,"%s/fsb.%03d", storePath, partnum);
  partfsb = createIndexStore(name, "partfsb", sizeof(gkStrobeFragment), 1);
  sprintf(name,"%s/qsb.%03d", storePath, partnum);
  partqsb = createStringStore(name, "partqsb");
}


void
gkStore::gkStore_open(int writable) {
  char  name[FILENAME_MAX];
  char   mode[4];
  FILE  *gkpinfo;

  sprintf(name,"%s/inf", storePath);
  errno = 0;
  gkpinfo = fopen(name, "r");
  if (errno) {
    fprintf(stderr, "failed to open gatekeeper store '%s': %s\n", name, strerror(errno));
    exit(1);
  }

  if (1 != AS_UTL_safeRead(gkpinfo, &inf, "gkStore_open:header", sizeof(gkStoreInfo), 1)) {
    fprintf(stderr, "failed to open gatekeeper store '%s': couldn't read the header (%s)\n", name, strerror(errno));
    exit(1);
  }

  fclose(gkpinfo);

  if (inf.gkMagic != 1) {
    fprintf(stderr, "gkStore_open()-- Invalid magic; corrupt %s/inf?\n", name);
    exit(1);
  }

  if (inf.gkVersion != 4) {
    fprintf(stderr, "gkStore_open()-- Invalid version!  Found version %d, code supports version 4.\n", inf.gkVersion);
    exit(1);
  }

  if ((inf.gkLibrarySize        != sizeof(gkLibrary)) ||
      (inf.gkPackedFragmentSize != sizeof(gkPackedFragment)) ||
      (inf.gkNormalFragmentSize != sizeof(gkNormalFragment)) ||
      (inf.gkStrobeFragmentSize != sizeof(gkStrobeFragment))) {
    fprintf(stderr, "gkStore_open()--  ERROR!  Incorrect element sizes; code and store are incompatible.\n");
    fprintf(stderr, "                  gkLibrary:          store %5d   code %5d bytes\n", inf.gkLibrarySize, (int)sizeof(gkLibrary));
    fprintf(stderr, "                  gkPackedFragment:   store %5d   code %5d bytes\n", inf.gkPackedFragmentSize, (int)sizeof(gkPackedFragment));
    fprintf(stderr, "                  gkNormalFragment:   store %5d   code %5d bytes\n", inf.gkNormalFragmentSize, (int)sizeof(gkNormalFragment));
    fprintf(stderr, "                  gkStrobeFragment:   store %5d   code %5d bytes\n", inf.gkStrobeFragmentSize, (int)sizeof(gkStrobeFragment));
    exit(1);
  }

  assert((writable == 0) || (writable == 1));

  if (writable) {
    isReadOnly = 0;
    isCreating = 0;
    strcpy(mode, "r+");
  } else {
    isReadOnly = 1;
    isCreating = 0;
    strcpy(mode, "r");
  }

  sprintf(name,"%s/fpk", storePath);
  fpk   = openStore(name, mode);

  sprintf(name,"%s/fnm", storePath);
  fnm   = openStore(name, mode);
  sprintf(name,"%s/snm", storePath);
  snm   = openStore(name, mode);
  sprintf(name,"%s/qnm", storePath);
  qnm   = openStore(name, mode);

  sprintf(name,"%s/fsb", storePath);
  fsb   = openStore(name, mode);
  sprintf(name,"%s/ssb", storePath);
  ssb   = openStore(name, mode);
  sprintf(name,"%s/qsb", storePath);
  qsb   = openStore(name, mode);

  sprintf(name,"%s/lib", storePath);
  lib   = openStore(name, mode);
  lib   = convertStoreToMemoryStore(lib);

  sprintf(name,"%s/uid", storePath);
  uid    = openStore(name, mode);
  
  sprintf(name, "%s/plc", storePath);
  plc    = openStore(name, mode);
  plc    = convertStoreToMemoryStore(plc); 

  if ((NULL == fpk) ||
      (NULL == fnm) || (NULL == snm) || (NULL == qnm) ||
      (NULL == fsb) || (NULL == ssb) || (NULL == qsb) ||
      (NULL == lib) || (NULL == uid) || (NULL == plc)) {
    fprintf(stderr,"Failed to open gkpStore '%s'.\n", storePath);
    exit(1);
  }
}


void
gkStore::gkStore_create(void) {
  char  name[FILENAME_MAX];
  FILE  *gkpinfo;

  isReadOnly = 0;
  isCreating = 1;

  sprintf(name,"%s/inf", storePath);
  if (AS_UTL_fileExists(name, FALSE, TRUE)) {
    fprintf(stderr, "GateKeeper Store '%s' exists; will not create a new one on top of it.\n", storePath);
    exit(1);
  }

  AS_UTL_mkdir(storePath);

  inf.gkMagic              = 1;
  inf.gkVersion            = 4;
  inf.gkLibrarySize        = sizeof(gkLibrary);
  inf.gkPackedFragmentSize = sizeof(gkPackedFragment);
  inf.gkNormalFragmentSize = sizeof(gkNormalFragment);
  inf.gkStrobeFragmentSize = sizeof(gkStrobeFragment);
  inf.gkPlacementSize      = sizeof(gkPlacement);

  sprintf(name,"%s/inf", storePath);
  errno = 0;
  gkpinfo = fopen(name, "w");
  if (errno) {
    fprintf(stderr, "failed to create gatekeeper store '%s': %s\n", name, strerror(errno));
    exit(1);
  }

  AS_UTL_safeWrite(gkpinfo, &inf, "creategkStore:header", sizeof(gkStoreInfo), 1);
  if (fclose(gkpinfo)) {
    fprintf(stderr, "failed to create gatekeeper store '%s': %s\n", name, strerror(errno));
    exit(1);
  }

  sprintf(name,"%s/fpk", storePath);
  fpk = createIndexStore(name, "fpk", sizeof(gkPackedFragment), 1);

  sprintf(name,"%s/fnm", storePath);
  fnm = createIndexStore(name, "fnm", sizeof(gkNormalFragment), 1);
  sprintf(name,"%s/snm", storePath);
  snm = createStringStore(name, "snm");
  sprintf(name,"%s/qnm", storePath);
  qnm = createStringStore(name, "qnm");

  sprintf(name,"%s/fsb", storePath);
  fsb = createIndexStore(name, "fsb", sizeof(gkStrobeFragment), 1);
  sprintf(name,"%s/ssb", storePath);
  ssb = createStringStore(name, "ssb");
  sprintf(name,"%s/qsb", storePath);
  qsb = createStringStore(name, "qsb");

  sprintf(name,"%s/lib", storePath);
  lib = createIndexStore(name, "lib", sizeof(gkLibrary), 1);
  lib = convertStoreToMemoryStore(lib);

  sprintf(name,"%s/uid", storePath);
  uid = createStringStore(name, "uid");

  sprintf(name, "%s/plc", storePath);
  plc = createIndexStore(name, "plc", sizeof(gkPlacement), 1);
  plc = convertStoreToMemoryStore(plc);
  
  sprintf(name,"%s/f2p", storePath);
  FRGtoPLC = CreateScalarHashTable_AS();
  SaveHashTable_AS(name, FRGtoPLC);
  
  sprintf(name,"%s/u2i", storePath);
  UIDtoIID = CreateScalarHashTable_AS();
  SaveHashTable_AS(name, UIDtoIID);

  IIDmax = 1048576;
  IIDtoTYPE = (uint8  *)safe_malloc(sizeof(uint8)  * IIDmax);
  IIDtoTIID = (uint32 *)safe_malloc(sizeof(uint32) * IIDmax);
}



gkStore::gkStore(const char *path, int creatable, int writable) {
  char   name[FILENAME_MAX];

  gkStore_clear();

  strcpy(storePath, path);

  sprintf(name,"%s/inf", storePath);

  if ((AS_UTL_fileExists(name, FALSE, writable)) && (creatable == 0)) {
    gkStore_open(writable);
  } else if (creatable) {
    gkStore_create();
  } else {
    fprintf(stderr, "gkStore::gkStore()-- GateKeeper Store '%s' doesn't exist.\n", storePath);
    exit(1);
  }

  //  Configure external clear ranges.
  clearRange = new gkClearRange * [AS_READ_CLEAR_NUM];

  for (uint32 i=0; i<AS_READ_CLEAR_NUM; i++)
    clearRange[i] = new gkClearRange(this, i, FALSE);

  AS_UID_setGatekeeper(this);
}



gkStore::~gkStore() {
  char  name[FILENAME_MAX];
  FILE *gkpinfo;

  if ((isReadOnly == 0) && (partnum == 0)) {
    sprintf(name,"%s/inf", storePath);
    errno = 0;
    gkpinfo = fopen(name, "w");
    if (errno) {
      fprintf(stderr, "failed to write gatekeeper store into to '%s': %s\n", name, strerror(errno));
      exit(1);
    }
    AS_UTL_safeWrite(gkpinfo, &inf, "closegkStore:header", sizeof(gkStoreInfo), 1);
    if (fclose(gkpinfo)) {
      fprintf(stderr, "failed to close gatekeeper store '%s': %s\n", name, strerror(errno));
      exit(1);
    }
  }

  closeStore(fpk);

  closeStore(fnm);
  closeStore(snm);
  closeStore(qnm);

  closeStore(fsb);
  closeStore(ssb);
  closeStore(qsb);

  closeStore(lib);

  closeStore(uid);

  closeStore(plc);
  DeleteHashTable_AS(FRGtoPLC);
  
  DeleteHashTable_AS(UIDtoIID);
  DeleteHashTable_AS(STRtoUID);

  safe_free(frgUID);

  if (clearRange)
    for (uint32 i=0; i<AS_READ_CLEAR_NUM; i++)
      delete clearRange[i];
  delete [] clearRange;

  safe_free(IIDtoTYPE);
  safe_free(IIDtoTIID);

  closeStore(partfpk);
  closeStore(partfnm);
  closeStore(partfsb);
  closeStore(partqnm);
  closeStore(partqsb);

  DeleteHashTable_AS(partmap);

  gkStore_clear();
}


void
gkStore::gkStore_clear(void) {
  memset(storePath, 0, sizeof(char) * FILENAME_MAX);

  isReadOnly = 1;

  memset(&inf, 0, sizeof(gkStoreInfo));

  fpk = NULL;

  fnm = NULL;
  snm = NULL;
  qnm = NULL;

  fsb = NULL;
  ssb = NULL;
  qsb = NULL;

  lib = NULL;

  uid = NULL;

  plc      = NULL;
  FRGtoPLC = NULL;
  
  UIDtoIID = NULL;
  STRtoUID = NULL;

  frgUID = NULL;

  IIDmax    = 0;
  IIDtoTYPE = 0;
  IIDtoTIID = NULL;

  clearRange = NULL;

  partnum = 0;

  partfpk = NULL;
  partfnm = partqnm = NULL;
  partfsb = partqsb = NULL;

  partmap = NULL;
}


void
gkStore::gkStore_delete(void) {
  char   name[FILENAME_MAX];

  //  This function does both a ~gkStore (needed to close open files, and etc) and then removes the
  //  files from disk.  It does not handle a partitioned store.

  closeStore(fpk);

  closeStore(fnm);
  closeStore(snm);
  closeStore(qnm);

  closeStore(fsb);
  closeStore(ssb);
  closeStore(qsb);

  closeStore(lib);

  closeStore(uid);

  closeStore(plc);
  DeleteHashTable_AS(FRGtoPLC);
  
  DeleteHashTable_AS(UIDtoIID);
  DeleteHashTable_AS(STRtoUID);

  safe_free(frgUID);

  safe_free(IIDtoTYPE);
  safe_free(IIDtoTIID);

  closeStore(partfpk);
  closeStore(partfnm);
  closeStore(partfsb);
  closeStore(partqnm);
  closeStore(partqsb);

  DeleteHashTable_AS(partmap);

  //  Remove files (and close/purge clear ranges).

  sprintf(name,"%s/inf", storePath);  unlink(name);

  sprintf(name,"%s/fpk", storePath);  unlink(name);

  sprintf(name,"%s/fnm", storePath);  unlink(name);
  sprintf(name,"%s/snm", storePath);  unlink(name);
  sprintf(name,"%s/qnm", storePath);  unlink(name);

  sprintf(name,"%s/fsb", storePath);  unlink(name);
  sprintf(name,"%s/ssb", storePath);  unlink(name);
  sprintf(name,"%s/qsb", storePath);  unlink(name);

  sprintf(name,"%s/lib", storePath);  unlink(name);
  sprintf(name,"%s/uid", storePath);  unlink(name);
  sprintf(name,"%s/plc", storePath);  unlink(name);
  sprintf(name,"%s/f2p", storePath);  unlink(name);
  sprintf(name,"%s/u2i", storePath);  unlink(name);

  for (int32 i=0; i<AS_READ_CLEAR_NUM; i++) {
    gkStore_purgeClearRange(i);
    delete clearRange[i];
  }
  delete [] clearRange;

  rmdir(storePath);

  gkStore_clear();
}



////////////////////////////////////////////////////////////////////////////////


void
gkStore::gkStore_loadPartition(uint32 partition) {
  char       name[FILENAME_MAX];
  int        i, f, e;

  assert(partmap    == NULL);
  assert(isCreating == 0);

  if (isReadOnly == 0)
    fprintf(stderr, "WARNING:  loading a partition from a writable gkpStore.\n");
  assert(isReadOnly == 1);

  partnum = partition;

  sprintf(name,"%s/fpk.%03d", storePath, partnum);
  if (AS_UTL_fileExists(name, FALSE, FALSE) == 0) {
    fprintf(stderr, "gkStore_loadPartition()--  Partition %d doesn't exist; normal store used instead.\n", partnum);
    return;
  }

  //  load all our data

  sprintf(name,"%s/fpk.%03d", storePath, partnum);
  partfpk = loadStorePartial(name, 0, 0);

  sprintf(name,"%s/fnm.%03d", storePath, partnum);
  partfnm = loadStorePartial(name, 0, 0);
  sprintf(name,"%s/qnm.%03d", storePath, partnum);
  partqnm = loadStorePartial(name, 0, 0);

  sprintf(name,"%s/fsb.%03d", storePath, partnum);
  partfsb = loadStorePartial(name, 0, 0);
  sprintf(name,"%s/qsb.%03d", storePath, partnum);
  partqsb = loadStorePartial(name, 0, 0);

  //  zip through the frags and build a map from iid to the frag record

  partmap = CreateScalarHashTable_AS();

  f = getFirstElemStore(partfpk);
  e = getLastElemStore(partfpk);
  for (i=f; i<=e; i++) {
    gkPackedFragment *p = (gkPackedFragment *)getIndexStorePtr(partfpk, i);
    if (InsertInHashTable_AS(partmap,
                             (uint64)p->readIID, 0,
                             (INTPTR)(p), 0) != HASH_SUCCESS)
      assert(0);
  }

  f = getFirstElemStore(partfnm);
  e = getLastElemStore(partfnm);
  for (i=f; i<=e; i++) {
    gkNormalFragment *p = (gkNormalFragment *)getIndexStorePtr(partfnm, i);
    if (InsertInHashTable_AS(partmap,
                             (uint64)p->readIID, 0,
                             (INTPTR)(p), 0) != HASH_SUCCESS)
      assert(0);
  }

  f = getFirstElemStore(partfsb);
  e = getLastElemStore(partfsb);
  for (i=f; i<=e; i++) {
    gkStrobeFragment *p = (gkStrobeFragment *)getIndexStorePtr(partfsb, i);
    if (InsertInHashTable_AS(partmap,
                             (uint64)p->readIID, 0,
                             (INTPTR)(p), 0) != HASH_SUCCESS)
      assert(0);
  }
}

////////////////////////////////////////////////////////////////////////////////


void
gkStore::gkStore_load(AS_IID beginIID, AS_IID endIID, int flags) {
  uint32  stType    = 0;
  uint32  stTiid    = 0;
  uint32  edType    = 0;
  uint32  edTiid    = 0;

  int64  firstsm   = 0, lastsm = 0;
  int64  firstmd   = 0, lastmd = 0;
  int64  firstlg   = 0, lastlg = 0;

  //  Position the metadata stream -- this code is similar to
  //  gkStore_loadPartition.

  assert(partmap    == NULL);
  assert(isReadOnly == 1);
  assert(isCreating == 0);

  if (beginIID == 0)   beginIID = 1;
  if (endIID == 0)     endIID   = gkStore_getNumFragments();

  gkStore_decodeTypeFromIID(beginIID, stType, stTiid);
  gkStore_decodeTypeFromIID(endIID,   edType, edTiid);

  if (stType == edType) {
    switch (stType) {
      case GKFRAGMENT_PACKED:
        firstsm = stTiid;
        lastsm  = edTiid;
        break;
      case GKFRAGMENT_NORMAL:
        firstmd = stTiid;
        lastmd  = edTiid;
        break;
      case GKFRAGMENT_STROBE:
        firstlg = stTiid;
        lastlg  = edTiid;
        break;
    }
  } else if ((stType == GKFRAGMENT_PACKED) && (edType == GKFRAGMENT_NORMAL)) {
    firstsm = stTiid;
    lastsm  = STREAM_UNTILEND;
    firstmd = STREAM_FROMSTART;
    lastmd  = edTiid;
  } else if ((stType == GKFRAGMENT_PACKED) && (edType == GKFRAGMENT_STROBE)) {
    firstsm = stTiid;
    lastsm  = STREAM_UNTILEND;
    firstmd = STREAM_FROMSTART;
    lastmd  = STREAM_UNTILEND;
    firstlg = STREAM_FROMSTART;
    lastlg  = edTiid;
  } else if ((stType == GKFRAGMENT_NORMAL) && (edType == GKFRAGMENT_STROBE)) {
    firstmd = stTiid;
    lastmd  = STREAM_UNTILEND;
    firstlg = STREAM_FROMSTART;
    lastlg  = edTiid;
  } else {
    assert(0);
  }

  //  Load the stores.  If we're loading all the way till the end, the
  //  last+1 fragment doesn't exist.  In this case, we (ab)use the
  //  fact that convertStoreToPartialMemoryStore() treats 0 as meaning
  //  "from the start" or "till the end".

  if (firstsm != lastsm) {
    fpk = convertStoreToPartialMemoryStore(fpk, firstsm, lastsm);
  }

  if (firstmd != lastmd) {
    gkNormalFragment mdbeg;
    gkNormalFragment mdend;

    mdbeg.seqOffset = mdend.seqOffset = 0;  //  Abuse.
    mdbeg.qltOffset = mdend.qltOffset = 0;

    if (firstmd != STREAM_FROMSTART)
      getIndexStore(fnm, firstmd, &mdbeg);

    if ((lastmd != STREAM_UNTILEND) && (lastmd + 1 <= getLastElemStore(fnm)))
      getIndexStore(fnm, lastmd+1, &mdend);

    fnm = convertStoreToPartialMemoryStore(fnm, firstmd, lastmd);

    if (flags == GKFRAGMENT_SEQ)
      snm = convertStoreToPartialMemoryStore(snm, mdbeg.seqOffset, mdend.seqOffset);

    if (flags == GKFRAGMENT_QLT)
      qnm = convertStoreToPartialMemoryStore(qnm, mdbeg.qltOffset, mdend.qltOffset);
  }

  if (firstlg != lastlg) {
    gkStrobeFragment lgbeg;
    gkStrobeFragment lgend;

    lgbeg.seqOffset = lgend.seqOffset = 0;  //  Abuse.
    lgbeg.qltOffset = lgend.qltOffset = 0;

    if (firstlg != STREAM_FROMSTART)
      getIndexStore(fsb, firstlg, &lgbeg);

    if ((lastlg != STREAM_UNTILEND) && (lastlg + 1 <= getLastElemStore(fsb)))
      getIndexStore(fsb, lastlg+1, &lgend);

    fsb = convertStoreToPartialMemoryStore(fsb, firstlg, lastlg);

    if (flags == GKFRAGMENT_SEQ)
      ssb = convertStoreToPartialMemoryStore(ssb, lgbeg.seqOffset, lgend.seqOffset);

    if (flags == GKFRAGMENT_QLT)
      qsb = convertStoreToPartialMemoryStore(qsb, lgbeg.qltOffset, lgend.qltOffset);
  }
}


////////////////////////////////////////////////////////////////////////////////

void
gkStore::gkStore_decodeTypeFromIID(AS_IID iid, uint32& type, uint32& tiid) {
  type = 0;
  tiid = 0;

  if (isCreating) {
    type = IIDtoTYPE[iid];
    tiid = IIDtoTIID[iid];
  } else {
    if (iid <= inf.numShort + inf.numMedium + inf.numLong) {
      type = GKFRAGMENT_STROBE;
      tiid = iid - inf.numShort - inf.numMedium;
    }
    if (iid <= inf.numShort + inf.numMedium) {
      type = GKFRAGMENT_NORMAL;
      tiid = iid - inf.numShort;
    }
    if (iid <= inf.numShort) {
      type = GKFRAGMENT_PACKED;
      tiid = iid;
    }
  }

  if (tiid == 0) {
    fprintf(stderr, "gkStore_decodeTypeFromIID()-- ERROR:  fragment iid %d is out of range.\n", iid);
    fprintf(stderr, "gkStore_decodeTypeFromIID()--         numShort=%d numMedium=%d numLong=%d\n",
            inf.numShort, inf.numMedium, inf.numLong);
    assert(0);
  }
}



void
gkStore::gkStore_getFragmentData(gkStream *gst, gkFragment *fr, uint32 flags) {

  assert((flags == GKFRAGMENT_INF) || (flags == GKFRAGMENT_SEQ) || (flags == GKFRAGMENT_QLT));

  fr->gkp = this;

  //  Get this out of the way first, it's a complete special case.
  //  It's also the easy case.  ;-)
  //
  if (fr->type == GKFRAGMENT_PACKED) {
    fr->hasSEQ = 1;
    fr->hasQLT = 1;

    if (fr->enc == NULL) {
      fr->enc = (char *)safe_malloc(sizeof(char) * AS_READ_MAX_NORMAL_LEN + 1);
      fr->seq = (char *)safe_malloc(sizeof(char) * AS_READ_MAX_NORMAL_LEN + 1);
      fr->qlt = (char *)safe_malloc(sizeof(char) * AS_READ_MAX_NORMAL_LEN + 1);
    }

    decodeSequenceQuality(fr->fr.packed.enc, fr->seq, fr->qlt);

    return;
  }

  uint32  actLen = 0;
  int64   nxtOff = 0;

  uint32  seqLen = fr->gkFragment_getSequenceLength();
  int64   seqOff = fr->gkFragment_getSequenceOffset();

  uint32  qltLen = fr->gkFragment_getQualityLength();
  int64   qltOff = fr->gkFragment_getQualityOffset();

  StoreStruct  *store  = NULL;
  StreamStruct *stream = NULL;


  if (flags == GKFRAGMENT_SEQ) {
    fr->hasSEQ = 1;

    if (fr->enc == NULL) {
      fr->enc = (char *)safe_malloc(sizeof(char) * AS_READ_MAX_NORMAL_LEN + 1);
      fr->seq = (char *)safe_malloc(sizeof(char) * AS_READ_MAX_NORMAL_LEN + 1);
      fr->qlt = (char *)safe_malloc(sizeof(char) * AS_READ_MAX_NORMAL_LEN + 1);
    }

    assert(partmap == NULL);

    if (gst == NULL) {
      store = (fr->type == GKFRAGMENT_NORMAL) ? snm : ssb;
      getStringStore(store, seqOff, fr->enc, AS_READ_MAX_NORMAL_LEN, &actLen, &nxtOff);
    } else {
      stream = (fr->type == GKFRAGMENT_NORMAL) ? gst->snm : gst->ssb;
      nextStream(stream, fr->enc, AS_READ_MAX_NORMAL_LEN, &actLen);
    }

    decodeSequence(fr->enc, fr->seq, seqLen);

    assert(fr->seq[seqLen] == 0);
  }


  if (flags == GKFRAGMENT_QLT) {
    fr->hasSEQ = 1;
    fr->hasQLT = 1;

    if (fr->enc == NULL) {
      fr->enc = (char *)safe_malloc(sizeof(char) * AS_READ_MAX_NORMAL_LEN + 1);
      fr->seq = (char *)safe_malloc(sizeof(char) * AS_READ_MAX_NORMAL_LEN + 1);
      fr->qlt = (char *)safe_malloc(sizeof(char) * AS_READ_MAX_NORMAL_LEN + 1);
    }

    if (partmap) {
      store = (fr->type == GKFRAGMENT_NORMAL) ? partqnm : partqsb;
      getStringStore(store, qltOff, fr->enc, AS_READ_MAX_NORMAL_LEN, &actLen, &nxtOff);
    } else if (gst == NULL) {
      store = (fr->type == GKFRAGMENT_NORMAL) ? qnm : qsb;
      getStringStore(store, qltOff, fr->enc, AS_READ_MAX_NORMAL_LEN, &actLen, &nxtOff);
    } else {
      stream = (fr->type == GKFRAGMENT_NORMAL) ? gst->qnm : gst->qsb;
      nextStream(stream, fr->enc, AS_READ_MAX_NORMAL_LEN, &actLen);
    }

    //fr->enc[actLen] = 0;

    decodeSequenceQuality(fr->enc, fr->seq, fr->qlt);
    assert(fr->seq[seqLen] == 0);
    assert(fr->qlt[seqLen] == 0);
  }
}



void
gkStore::gkStore_getFragment(AS_IID iid, gkFragment *fr, int32 flags) {

  fr->hasSEQ = 0;
  fr->hasQLT = 0;

  fr->gkp = this;

  gkStore_decodeTypeFromIID(iid, fr->type, fr->tiid);

  //fprintf(stderr, "gkStore_getFragment()--  Retrieving IID=%d from %d,%d\n", iid, fr->type, fr->tiid);

  if (partmap) {
    //  If partitioned, we have everything in memory.  This is keyed
    //  off of the global IID.
    void *p = (void *)(INTPTR)LookupValueInHashTable_AS(partmap, iid, 0);

    if (p == NULL)
      fprintf(stderr, "getFrag()-- ERROR!  IID "F_IID" not in partition!\n", iid);
    assert(p != NULL);

    assert(fr->isGKP == 0);

    switch (fr->type) {
      case GKFRAGMENT_PACKED:
        memcpy(&fr->fr.packed, p, sizeof(gkPackedFragment));
        break;
      case GKFRAGMENT_NORMAL:
        memcpy(&fr->fr.normal, p, sizeof(gkNormalFragment));
        break;
      case GKFRAGMENT_STROBE:
        memcpy(&fr->fr.strobe, p, sizeof(gkStrobeFragment));
        break;
    }

  } else {
    //  Not paritioned, load from disk store.
    switch (fr->type) {
      case GKFRAGMENT_PACKED:
        getIndexStore(fpk, fr->tiid, &fr->fr.packed);
        break;
      case GKFRAGMENT_NORMAL:
        getIndexStore(fnm, fr->tiid, &fr->fr.normal);
        break;
      case GKFRAGMENT_STROBE:
        getIndexStore(fsb, fr->tiid, &fr->fr.strobe);
        break;
    }
  }

  //  GatekeeperMode assumes these are set.  Set them.  Currently only used by sffToCA, when it
  //  loads frags from the temporary store before detecting mate pairs.
  if (fr->isGKP) {
    fr->gkFragment_getClearRegion(fr->clrBgn, fr->clrEnd, AS_READ_CLEAR_CLR);
    fr->gkFragment_getClearRegion(fr->vecBgn, fr->vecEnd, AS_READ_CLEAR_VEC);
    fr->gkFragment_getClearRegion(fr->maxBgn, fr->maxEnd, AS_READ_CLEAR_MAX);
    fr->gkFragment_getClearRegion(fr->tntBgn, fr->tntEnd, AS_READ_CLEAR_TNT);
  }

  //fprintf(stderr, "gkStore_getFragment()--  Retrieved IID=%d (asked for %d) from %d,%d\n", fr->gkFragment_getReadIID(), iid, fr->type, fr->tiid);

  gkStore_getFragmentData(NULL, fr, flags);
}




void
gkStore::gkStore_setFragment(gkFragment *fr) {
  assert(partmap    == NULL);
  assert(isReadOnly == 0);

  //  Sanity check that type and type IID agree.  These change while the store is being built.
  {
    uint32  type;
    uint32  tiid;

    gkStore_decodeTypeFromIID(fr->gkFragment_getReadIID(), type, tiid);

    assert(type == fr->type);
    assert(tiid == fr->tiid);
  }

  //fprintf(stderr, "gkStore_setFragment()--  Setting IID=%d to %d,%d\n", fr->gkFragment_getReadIID(), fr->type, fr->tiid);

  switch (fr->type) {
    case GKFRAGMENT_PACKED:
      setIndexStore(fpk, fr->tiid, &fr->fr.packed);
      break;
    case GKFRAGMENT_NORMAL:
      setIndexStore(fnm, fr->tiid, &fr->fr.normal);
      break;
    case GKFRAGMENT_STROBE:
      setIndexStore(fsb, fr->tiid, &fr->fr.strobe);
      break;
  }
}




//  Delete fragment with iid from the store.  If the fragment has a
//  mate, remove the mate relationship from both fragmentss.
//
//  If the mate is supplied, delete it too.
//
void
gkStore::gkStore_delFragment(AS_IID iid, bool deleteMateFrag) {
  gkFragment   fr;
  int32        mid = 0;

  assert(partmap    == NULL);
  assert(isReadOnly == 0);

  gkStore_getFragment(iid, &fr, GKFRAGMENT_INF);
  switch (fr.type) {
    case GKFRAGMENT_PACKED:
      mid = fr.fr.packed.mateIID;
      fr.fr.packed.deleted = 1;
      fr.fr.packed.mateIID = 0;
      break;
    case GKFRAGMENT_NORMAL:
      mid = fr.fr.normal.mateIID;
      fr.fr.normal.deleted = 1;
      fr.fr.normal.mateIID = 0;
      break;
    case GKFRAGMENT_STROBE:
      mid = fr.fr.strobe.mateIID;
      fr.fr.strobe.deleted = 1;
      fr.fr.strobe.mateIID = 0;
      break;
  }
  gkStore_setFragment(&fr);

  //  No mate, we're done.
  if (mid == 0)
    return;

  gkStore_getFragment(mid, &fr, GKFRAGMENT_INF);
  switch (fr.type) {
    case GKFRAGMENT_PACKED:
      fr.fr.packed.deleted = fr.fr.packed.deleted || deleteMateFrag;
      fr.fr.packed.mateIID = 0;
      break;
    case GKFRAGMENT_NORMAL:
      fr.fr.normal.deleted = fr.fr.normal.deleted || deleteMateFrag;
      fr.fr.normal.mateIID = 0;
      break;
    case GKFRAGMENT_STROBE:
      fr.fr.strobe.deleted = fr.fr.strobe.deleted || deleteMateFrag;
      fr.fr.strobe.mateIID = 0;
      break;
  }
  gkStore_setFragment(&fr);
}



void
gkStore::gkStore_addFragment(gkFragment *fr) {
  int encLen;

  assert(partmap    == NULL);
  assert(isReadOnly == 0);
  assert(isCreating == 1);

  assert(fr->type != GKFRAGMENT_ERROR);

  int32 iid = gkStore_getNumFragments() + 1;

  if (fr->clrBgn < fr->clrEnd)  clearRange[AS_READ_CLEAR_CLR]->gkClearRange_enableCreate();
  if (fr->vecBgn < fr->vecEnd)  clearRange[AS_READ_CLEAR_VEC]->gkClearRange_enableCreate();
  if (fr->maxBgn < fr->maxEnd)  clearRange[AS_READ_CLEAR_MAX]->gkClearRange_enableCreate();
  if (fr->tntBgn < fr->tntEnd)  clearRange[AS_READ_CLEAR_TNT]->gkClearRange_enableCreate();

  assert(strlen(fr->gkFragment_getSequence()) == fr->gkFragment_getSequenceLength());
  assert(strlen(fr->gkFragment_getQuality())  == fr->gkFragment_getQualityLength());

  assert((fr->gkFragment_getIsDeleted() == 1) || (fr->clrBgn < fr->clrEnd));

  switch (fr->type) {
    case GKFRAGMENT_PACKED:
      fr->tiid          = ++inf.numShort;
      fr->fr.packed.readIID = iid;

      assert(fr->tiid == getLastElemStore(fpk) + 1);

      clearRange[AS_READ_CLEAR_CLR]->gkClearRange_makeSpaceShort(fr->tiid, fr->clrBgn, fr->clrEnd);
      clearRange[AS_READ_CLEAR_VEC]->gkClearRange_makeSpaceShort(fr->tiid, fr->vecBgn, fr->vecEnd);
      clearRange[AS_READ_CLEAR_MAX]->gkClearRange_makeSpaceShort(fr->tiid, fr->maxBgn, fr->maxEnd);
      clearRange[AS_READ_CLEAR_TNT]->gkClearRange_makeSpaceShort(fr->tiid, fr->tntBgn, fr->tntEnd);

      fr->fr.packed.clearBeg = fr->clrBgn;
      fr->fr.packed.clearEnd = fr->clrEnd;
      break;

    case GKFRAGMENT_NORMAL:
      fr->tiid          = ++inf.numMedium;
      fr->fr.normal.readIID = iid;

      assert(fr->tiid == getLastElemStore(fnm) + 1);

      clearRange[AS_READ_CLEAR_CLR]->gkClearRange_makeSpaceMedium(fr->tiid, fr->clrBgn, fr->clrEnd);
      clearRange[AS_READ_CLEAR_VEC]->gkClearRange_makeSpaceMedium(fr->tiid, fr->vecBgn, fr->vecEnd);
      clearRange[AS_READ_CLEAR_MAX]->gkClearRange_makeSpaceMedium(fr->tiid, fr->maxBgn, fr->maxEnd);
      clearRange[AS_READ_CLEAR_TNT]->gkClearRange_makeSpaceMedium(fr->tiid, fr->tntBgn, fr->tntEnd);

      fr->fr.normal.clearBeg = fr->clrBgn;
      fr->fr.normal.clearEnd = fr->clrEnd;
      break;

    case GKFRAGMENT_STROBE:
      fr->tiid           = ++inf.numLong;
      fr->fr.strobe.readIID = iid;

      assert(fr->tiid == getLastElemStore(fsb) + 1);

      clearRange[AS_READ_CLEAR_CLR]->gkClearRange_makeSpaceLong(fr->tiid, fr->clrBgn, fr->clrEnd);
      clearRange[AS_READ_CLEAR_VEC]->gkClearRange_makeSpaceLong(fr->tiid, fr->vecBgn, fr->vecEnd);
      clearRange[AS_READ_CLEAR_MAX]->gkClearRange_makeSpaceLong(fr->tiid, fr->maxBgn, fr->maxEnd);
      clearRange[AS_READ_CLEAR_TNT]->gkClearRange_makeSpaceLong(fr->tiid, fr->tntBgn, fr->tntEnd);

      fr->fr.strobe.clearBeg = fr->clrBgn;
      fr->fr.strobe.clearEnd = fr->clrEnd;
      break;
  }

  assert(fr->tiid != 0);

  //  Set clear ranges...if they're defined.  !!NOTE!!  The "CLR"
  //  range MUST be set last since it is "the" clear range, and
  //  setClearRange() always updates the latest clear range along with
  //  the named region.  (That is, setting VEC will set both VEC and
  //  LATEST.)
  //
  if (fr->vecBgn < fr->vecEnd)  clearRange[AS_READ_CLEAR_VEC]->gkClearRange_setClearRegion(fr, fr->vecBgn, fr->vecEnd);
  if (fr->maxBgn < fr->maxEnd)  clearRange[AS_READ_CLEAR_MAX]->gkClearRange_setClearRegion(fr, fr->maxBgn, fr->maxEnd);
  if (fr->tntBgn < fr->tntEnd)  clearRange[AS_READ_CLEAR_TNT]->gkClearRange_setClearRegion(fr, fr->tntBgn, fr->tntEnd);
  if (fr->clrBgn < fr->clrEnd)  clearRange[AS_READ_CLEAR_CLR]->gkClearRange_setClearRegion(fr, fr->clrBgn, fr->clrEnd);

  if (IIDmax <= iid) {
    IIDmax *= 2;
    IIDtoTYPE = (uint8  *)safe_realloc(IIDtoTYPE, sizeof(uint8)  * IIDmax);
    IIDtoTIID = (uint32 *)safe_realloc(IIDtoTIID, sizeof(uint32) * IIDmax);
  }

  switch (fr->type) {
    case GKFRAGMENT_PACKED:
      assert(fr->seq[fr->fr.packed.seqLen] == 0);
      assert(fr->qlt[fr->fr.packed.seqLen] == 0);

      encodeSequenceQuality(fr->fr.packed.enc, fr->seq, fr->qlt);

      gkStore_setUIDtoIID(fr->fr.packed.readUID, fr->fr.packed.readIID, AS_IID_FRG);
      appendIndexStore(fpk, &fr->fr.packed);

      IIDtoTYPE[iid] = GKFRAGMENT_PACKED;
      IIDtoTIID[iid] = fr->tiid;
      break;

    case GKFRAGMENT_NORMAL:
      assert(fr->seq[fr->fr.normal.seqLen] == 0);
      assert(fr->qlt[fr->fr.normal.seqLen] == 0);

      fr->fr.normal.seqOffset = getLastElemStore(snm) + 1;
      fr->fr.normal.qltOffset = getLastElemStore(qnm) + 1;

      gkStore_setUIDtoIID(fr->fr.normal.readUID, fr->fr.normal.readIID, AS_IID_FRG);
      appendIndexStore(fnm, &fr->fr.normal);

      encLen = encodeSequence(fr->enc, fr->seq);
      appendStringStore(snm, fr->enc, encLen);

      encodeSequenceQuality(fr->enc, fr->seq, fr->qlt);
      appendStringStore(qnm, fr->enc, fr->fr.normal.seqLen);

      IIDtoTYPE[iid] = GKFRAGMENT_NORMAL;
      IIDtoTIID[iid] = fr->tiid;
      break;

    case GKFRAGMENT_STROBE:
      assert(fr->seq[fr->fr.strobe.seqLen] == 0);
      assert(fr->qlt[fr->fr.strobe.seqLen] == 0);

      fr->fr.strobe.seqOffset = getLastElemStore(ssb) + 1;
      fr->fr.strobe.qltOffset = getLastElemStore(qsb) + 1;

      gkStore_setUIDtoIID(fr->fr.strobe.readUID, fr->fr.strobe.readIID, AS_IID_FRG);
      appendIndexStore(fsb, &fr->fr.strobe);

      encLen = encodeSequence(fr->enc, fr->seq);
      appendStringStore(ssb, fr->enc, encLen);

      encodeSequenceQuality(fr->enc, fr->seq, fr->qlt);
      appendStringStore(qsb, fr->enc, fr->fr.strobe.seqLen);

      IIDtoTYPE[iid] = GKFRAGMENT_STROBE;
      IIDtoTIID[iid] = fr->tiid;
      break;
  }

  //  We loaded a fragment regardless of its deleted status.  This is
  //  just the count of fragments in the store.
  inf.frgLoaded++;

  if (fr->gkFragment_getIsDeleted()) {
    //  Errors are defined as a fragment loaded but marked as deleted.
    inf.frgErrors++;
  } else {
    //  Only living fragments count towards random fragments.
    if (fr->gkFragment_getIsNonRandom() == 0)
      inf.numRandom++;
  }
}


////////////////////////////////////////////////////////////////////////////////

void
gkStore::gkStore_buildPartitions(short *partitionMap, uint32 maxPart) {
  gkFragment        fr;

  assert(partmap    == NULL);
  assert(isReadOnly == 1);
  assert(isCreating == 0);

  //  Create the partitions by opening N copies of the gatekeeper store,
  //  and telling each one to make a partition.
  //
  gkStore **gkpart = new gkStore * [maxPart + 1];

  AS_PER_setBufferSize(512 * 1024);

  for (uint32 i=1; i<=maxPart; i++)
    gkpart[i] = new gkStore (storePath, i);

  //  And, finally, add stuff to each partition.
  //
  for (uint32 iid=1; iid<=gkStore_getNumFragments(); iid++) {
    int p;

    gkStore_getFragment(iid, &fr, GKFRAGMENT_QLT);

    p = partitionMap[iid];

    //  Check it's actually partitioned.  Deleted reads won't get
    //  assigned to a partition.
    if (p < 1)
      continue;


    if (fr.type == GKFRAGMENT_PACKED) {
      appendIndexStore(gkpart[p]->partfpk, &fr.fr.packed);
    }


    if (fr.type == GKFRAGMENT_NORMAL) {
      fr.fr.normal.seqOffset = -1;
      fr.fr.normal.qltOffset = getLastElemStore(gkpart[p]->partqnm) + 1;

      appendIndexStore(gkpart[p]->partfnm, &fr.fr.normal);

      encodeSequenceQuality(fr.enc, fr.seq, fr.qlt);
      appendStringStore(gkpart[p]->partqnm, fr.enc, fr.fr.normal.seqLen);
    }


    if (fr.type == GKFRAGMENT_STROBE) {
      fr.fr.strobe.seqOffset = -1;
      fr.fr.strobe.qltOffset = getLastElemStore(gkpart[p]->partqsb) + 1;

      appendIndexStore(gkpart[p]->partfsb, &fr.fr.strobe);

      encodeSequenceQuality(fr.enc, fr.seq, fr.qlt);
      appendStringStore(gkpart[p]->partqsb, fr.enc, fr.fr.strobe.seqLen);
    }
  }

  //  cleanup -- close all the stores

  for (uint32 i=1; i<=maxPart; i++)
    delete gkpart[i];

  delete gkpart;
}

////////////////////////////////////////////////////////////////////////////////
//
//
//
void
gkStore::gkStore_enableClearRange(uint32 which) {
  assert(partmap == NULL);
  //assert(isReadOnly == 0);  --  if not writable, this will allow us to return
  //  undefined clear ranges in, e.g., gatekeeper.  gkFragment_setClearRegion()
  //  then checks this assert.
  clearRange[which]->gkClearRange_enableCreate();
}


void
gkStore::gkStore_purgeClearRange(uint32 which) {
  assert(partmap == NULL);
  clearRange[which]->gkClearRange_purge();
}
