
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

static char *rcsid = "$Id: AS_PER_gkStore.C,v 1.1 2009-06-10 18:05:14 brianwalenz Exp $";

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

  sprintf(name,"%s/fsm.%03d", storePath, partnum);
  partfsm = createIndexStore(name, "partfsm", sizeof(gkShortFragment), 1);

  sprintf(name,"%s/fmd.%03d", storePath, partnum);
  partfmd = createIndexStore(name, "partfmd", sizeof(gkMediumFragment), 1);
  sprintf(name,"%s/qmd.%03d", storePath, partnum);
  partqmd = createStringStore(name, "partqmd");

  sprintf(name,"%s/flg.%03d", storePath, partnum);
  partflg = createIndexStore(name, "partflg", sizeof(gkLongFragment), 1);
  sprintf(name,"%s/qlg.%03d", storePath, partnum);
  partqlg = createStringStore(name, "partqlg");
}


gkStore::gkStore(const char *path, int creatable, int writable) {
  char   name[FILENAME_MAX];
  char   mode[4];
  FILE  *gkpinfo;

  gkStore_clear();

  strcpy(storePath, path);

  sprintf(name,"%s/inf", storePath);

  if ((AS_UTL_fileExists(name, 0, writable)) && (creatable == 0)) {
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

    if (inf.gkVersion != 2) {
      fprintf(stderr, "gkStore_open()-- Invalid version!  Found version %d, code supports version 2.\n", inf.gkVersion);
      exit(1);
    }

    if ((inf.gkLibrarySize        != sizeof(gkLibrary)) ||
        (inf.gkShortFragmentSize  != sizeof(gkShortFragment)) ||
        (inf.gkMediumFragmentSize != sizeof(gkMediumFragment)) ||
        (inf.gkLongFragmentSize   != sizeof(gkLongFragment))) {
      fprintf(stderr, "gkStore_open()--  ERROR!  Incorrect element sizes; code and store are incompatible.\n");
      fprintf(stderr, "                  gkLibrary:          store %5d   code %5d bytes\n", inf.gkLibrarySize, (int)sizeof(gkLibrary));
      fprintf(stderr, "                  gkShortFragment:    store %5d   code %5d bytes\n", inf.gkShortFragmentSize, (int)sizeof(gkShortFragment));
      fprintf(stderr, "                  gkMediumFragment:   store %5d   code %5d bytes\n", inf.gkMediumFragmentSize, (int)sizeof(gkMediumFragment));
      fprintf(stderr, "                  gkLongFragment:     store %5d   code %5d bytes\n", inf.gkLongFragmentSize, (int)sizeof(gkLongFragment));
      exit(1);
    }

    //  writable is -1 if we are creating a partition
    if (writable == -1)
      return;

    if (writable) {
      isReadOnly = 0;
      strcpy(mode, "r+");
    } else {
      isReadOnly = 1;
      strcpy(mode, "r");
    }

    sprintf(name,"%s/fsm", storePath);
    fsm   = openStore(name, mode);

    sprintf(name,"%s/fmd", storePath);
    fmd   = openStore(name, mode);
    sprintf(name,"%s/smd", storePath);
    smd   = openStore(name, mode);
    sprintf(name,"%s/qmd", storePath);
    qmd   = openStore(name, mode);

    sprintf(name,"%s/flg", storePath);
    flg   = openStore(name, mode);
    sprintf(name,"%s/slg", storePath);
    slg   = openStore(name, mode);
    sprintf(name,"%s/qlg", storePath);
    qlg   = openStore(name, mode);

    sprintf(name,"%s/lib", storePath);
    lib   = openStore(name, mode);
    lib   = convertStoreToMemoryStore(lib);

    sprintf(name,"%s/uid", storePath);
    uid = openStore(name, mode);

    if ((NULL == fsm) ||
        (NULL == fmd) || (NULL == smd) || (NULL == qmd) ||
        (NULL == flg) || (NULL == slg) || (NULL == qlg) ||
        (NULL == lib) || (NULL == uid)) {
      fprintf(stderr,"Failed to open gkpStore '%s'.\n", storePath);
      exit(1);
    }

  } else if (creatable) {
    isReadOnly = 0;

    if (AS_UTL_fileExists(name, 0, writable)) {
      fprintf(stderr, "GateKeeper Store '%s' exists; will not create a new one on top of it.\n", path);
      exit(1);
    }

    AS_UTL_mkdir(path);

    inf.gkMagic              = 1;
    inf.gkVersion            = 2;
    inf.gkLibrarySize        = sizeof(gkLibrary);
    inf.gkShortFragmentSize  = sizeof(gkShortFragment);
    inf.gkMediumFragmentSize = sizeof(gkMediumFragment);
    inf.gkLongFragmentSize   = sizeof(gkLongFragment);

    sprintf(name,"%s/inf", path);
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

    sprintf(name,"%s/fsm", path);
    fsm = createIndexStore(name, "fsm", sizeof(gkShortFragment), 1);

    sprintf(name,"%s/fmd", path);
    fmd = createIndexStore(name, "fmd", sizeof(gkMediumFragment), 1);
    sprintf(name,"%s/smd", path);
    smd = createStringStore(name, "smd");
    sprintf(name,"%s/qmd", path);
    qmd = createStringStore(name, "qmd");

    sprintf(name,"%s/flg", path);
    flg = createIndexStore(name, "flg", sizeof(gkLongFragment), 1);
    sprintf(name,"%s/slg", path);
    slg = createStringStore(name, "slg");
    sprintf(name,"%s/qlg", path);
    qlg = createStringStore(name, "qlg");

    sprintf(name,"%s/lib", path);
    lib = createIndexStore(name, "lib", sizeof(gkLibrary), 1);

    sprintf(name,"%s/uid", path);
    uid = createStringStore(name, "uid");

    sprintf(name,"%s/u2i", path);
    UIDtoIID = CreateScalarHashTable_AS();
    SaveHashTable_AS(name, UIDtoIID);

  } else {
    fprintf(stderr, "gkStore::gkStore()-- GateKeeper Store '%s' doesn't exist.\n", path);
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

  closeStore(fsm);

  closeStore(fmd);
  closeStore(smd);
  closeStore(qmd);

  closeStore(flg);
  closeStore(slg);
  closeStore(qlg);

  closeStore(lib);

  closeStore(uid);

  DeleteHashTable_AS(UIDtoIID);
  DeleteHashTable_AS(STRtoUID);

  safe_free(frgUID);

  if (clearRange)
    for (uint32 i=0; i<AS_READ_CLEAR_NUM; i++)
      delete clearRange[i];
  delete [] clearRange;

  closeStore(partfsm);
  closeStore(partfmd);
  closeStore(partflg);
  closeStore(partqmd);
  closeStore(partqlg);

  DeleteHashTable_AS(partmap);

  gkStore_clear();
}


void
gkStore::gkStore_clear(void) {
  memset(storePath, 0, sizeof(char) * FILENAME_MAX);

  isReadOnly = 1;

  memset(&inf, 0, sizeof(gkStoreInfo));

  fsm = NULL;

  fmd = NULL;
  smd = NULL;
  qmd = NULL;

  flg = NULL;
  slg = NULL;
  qlg = NULL;

  lib = NULL;

  uid = NULL;

  UIDtoIID = NULL;
  STRtoUID = NULL;

  frgUID = NULL;

  clearRange = NULL;

  partnum = 0;

  partfsm = NULL;
  partfmd = partqmd = NULL;
  partflg = partqlg = NULL;

  partmap = NULL;
}


void
gkStore::gkStore_delete(const char *path) {
  char   name[FILENAME_MAX];

  sprintf(name,"%s/inf", path);
  unlink(name);

  sprintf(name,"%s/bat", path);
  unlink(name);

  sprintf(name,"%s/fsm", path);
  unlink(name);

  sprintf(name,"%s/fmd", path);
  unlink(name);

  sprintf(name,"%s/flg", path);
  unlink(name);

  sprintf(name,"%s/lib", path);
  unlink(name);

  sprintf(name,"%s/seq", path);
  unlink(name);

  sprintf(name,"%s/qlt", path);
  unlink(name);

  sprintf(name,"%s/hps", path);
  unlink(name);

  sprintf(name,"%s/src", path);
  unlink(name);

  sprintf(name,"%s/uid", path);
  unlink(name);

  sprintf(name,"%s/u2i", path);
  unlink(name);

  rmdir(path);
}



////////////////////////////////////////////////////////////////////////////////


void
gkStore::gkStore_loadPartition(uint32 partition) {
  char       name[FILENAME_MAX];
  int        i, f, e;

  assert(partmap    == NULL);
  assert(isReadOnly == 1);

  partnum = partition;

  sprintf(name,"%s/fsm.%03d", storePath, partnum);
  if (AS_UTL_fileExists(name, 0, FALSE) == 0) {
    fprintf(stderr, "gkStore_loadPartition()--  Partition %d doesn't exist; normal store used instead.\n", partnum);
    return;
  }

  //  load all our data

  isReadOnly = 1;

  sprintf(name,"%s/fsm.%03d", storePath, partnum);
  partfsm = loadStorePartial(name, 0, 0);

  sprintf(name,"%s/fmd.%03d", storePath, partnum);
  partfmd = loadStorePartial(name, 0, 0);
  sprintf(name,"%s/qmd.%03d", storePath, partnum);
  partqmd = loadStorePartial(name, 0, 0);

  sprintf(name,"%s/flg.%03d", storePath, partnum);
  partflg = loadStorePartial(name, 0, 0);
  sprintf(name,"%s/qlg.%03d", storePath, partnum);
  partqlg = loadStorePartial(name, 0, 0);

  //  zip through the frags and build a map from iid to the frag record

  partmap = CreateScalarHashTable_AS();

  f = getFirstElemStore(partfsm);
  e = getLastElemStore(partfsm);
  for (i=f; i<=e; i++) {
    gkShortFragment *p = (gkShortFragment *)getIndexStorePtr(partfsm, i);
    if (InsertInHashTable_AS(partmap,
                             (uint64)p->readIID, 0,
                             (INTPTR)(p), 0) != HASH_SUCCESS)
      assert(0);
  }

  f = getFirstElemStore(partfmd);
  e = getLastElemStore(partfmd);
  for (i=f; i<=e; i++) {
    gkMediumFragment *p = (gkMediumFragment *)getIndexStorePtr(partfmd, i);
    if (InsertInHashTable_AS(partmap,
                             (uint64)p->readIID, 0,
                             (INTPTR)(p), 0) != HASH_SUCCESS)
      assert(0);
  }

  f = getFirstElemStore(partflg);
  e = getLastElemStore(partflg);
  for (i=f; i<=e; i++) {
    gkLongFragment *p = (gkLongFragment *)getIndexStorePtr(partflg, i);
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

  uint32  firstsm   = 0, lastsm = 0;
  uint32  firstmd   = 0, lastmd = 0;
  uint32  firstlg   = 0, lastlg = 0;

  //  Position the metadata stream -- this code is similar to
  //  gkStore_loadPartition.

  assert(partmap    == NULL);
  assert(isReadOnly == 1);

  gkStore_decodeTypeFromIID(beginIID, &stType, &stTiid);
  gkStore_decodeTypeFromIID(endIID,   &edType, &edTiid);

  if (stType == edType) {
    switch (stType) {
      case GKFRAGMENT_SHORT:
        firstsm = stTiid;
        lastsm  = edTiid;
        break;
      case GKFRAGMENT_MEDIUM:
        firstmd = stTiid;
        lastmd  = edTiid;
        break;
      case GKFRAGMENT_LONG:
        firstlg = stTiid;
        lastlg  = edTiid;
        break;
    }
  } else if ((stType == GKFRAGMENT_SHORT) && (edType == GKFRAGMENT_MEDIUM)) {
    firstsm = stTiid;
    lastsm  = STREAM_UNTILEND;
    firstmd = STREAM_FROMSTART;
    lastmd  = edTiid;
  } else if ((stType == GKFRAGMENT_SHORT) && (edType == GKFRAGMENT_LONG)) {
    firstsm = stTiid;
    lastsm  = STREAM_UNTILEND;
    firstmd = STREAM_FROMSTART;
    lastmd  = STREAM_UNTILEND;
    firstlg = STREAM_FROMSTART;
    lastlg  = edTiid;
  } else if ((stType == GKFRAGMENT_MEDIUM) && (edType == GKFRAGMENT_LONG)) {
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
    fsm = convertStoreToPartialMemoryStore(fsm, firstsm, lastsm);
  }

  if (firstmd != lastmd) {
    gkMediumFragment mdbeg;
    gkMediumFragment mdend;

    mdbeg.seqOffset = mdend.seqOffset = 0;
    mdbeg.qltOffset = mdend.qltOffset = 0;

    getIndexStore(fmd, firstmd, &mdbeg);

    if (lastmd + 1 <= getLastElemStore(fmd))
      getIndexStore(fmd, lastmd+1, &mdend);

    fmd = convertStoreToPartialMemoryStore(fmd, firstmd, lastmd);

    if (flags == GKFRAGMENT_SEQ)
      smd = convertStoreToPartialMemoryStore(smd, mdbeg.seqOffset, mdend.seqOffset);

    if (flags == GKFRAGMENT_QLT)
      qmd = convertStoreToPartialMemoryStore(qmd, mdbeg.qltOffset, mdend.qltOffset);
  }

  if (firstlg != lastlg) {
    gkLongFragment lgbeg;
    gkLongFragment lgend;

    lgbeg.seqOffset = lgend.seqOffset = 0;
    lgbeg.qltOffset = lgend.qltOffset = 0;

    getIndexStore(flg, firstlg, &lgbeg);

    if (lastlg + 1 <= getLastElemStore(flg))
      getIndexStore(flg, lastlg+1, &lgend);

    flg = convertStoreToPartialMemoryStore(flg, firstlg, lastlg);

    if (flags == GKFRAGMENT_SEQ)
      slg = convertStoreToPartialMemoryStore(slg, lgbeg.seqOffset, lgend.seqOffset);

    if (flags == GKFRAGMENT_QLT)
      qlg = convertStoreToPartialMemoryStore(qlg, lgbeg.qltOffset, lgend.qltOffset);
  }
}


////////////////////////////////////////////////////////////////////////////////

void
gkStore::gkStore_decodeTypeFromIID(AS_IID iid, uint32 *type, uint32 *tiid) {
  *type = 0;
  *tiid = 0;

  if (iid <= inf.numShort + inf.numMedium + inf.numLong) {
    *type = GKFRAGMENT_LONG;
    *tiid = iid - inf.numShort - inf.numMedium;
  }
  if (iid <= inf.numShort + inf.numMedium) {
    *type = GKFRAGMENT_MEDIUM;
    *tiid = iid - inf.numShort;
  }
  if (iid <= inf.numShort) {
    *type = GKFRAGMENT_SHORT;
    *tiid = iid;
  }
  if (*tiid == 0) {
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
  if (fr->type == GKFRAGMENT_SHORT) {
    fr->hasSEQ = 1;
    fr->hasQLT = 1;

    if (fr->enc == NULL) {
      fr->enc = (char *)safe_malloc(sizeof(char) * AS_READ_MAX_LONG_LEN + 1);
      fr->seq = (char *)safe_malloc(sizeof(char) * AS_READ_MAX_LONG_LEN + 1);
      fr->qlt = (char *)safe_malloc(sizeof(char) * AS_READ_MAX_LONG_LEN + 1);
    }

    decodeSequenceQuality(fr->fr.sm.enc, fr->seq, fr->qlt);

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
      fr->enc = (char *)safe_malloc(sizeof(char) * AS_READ_MAX_LONG_LEN + 1);
      fr->seq = (char *)safe_malloc(sizeof(char) * AS_READ_MAX_LONG_LEN + 1);
      fr->qlt = (char *)safe_malloc(sizeof(char) * AS_READ_MAX_LONG_LEN + 1);
    }

    assert(partmap == NULL);

    if (gst == NULL) {
      store = (fr->type == GKFRAGMENT_MEDIUM) ? smd : slg;
      getStringStore(store, seqOff, fr->enc, AS_READ_MAX_LONG_LEN, &actLen, &nxtOff);
    } else {
      stream = (fr->type == GKFRAGMENT_MEDIUM) ? gst->smd : gst->slg;
      nextStream(stream, fr->enc, AS_READ_MAX_LONG_LEN, &actLen);
    }

    decodeSequence(fr->enc, fr->seq, seqLen);

    assert(fr->seq[seqLen] == 0);
  }


  if (flags == GKFRAGMENT_QLT) {
    fr->hasSEQ = 1;
    fr->hasQLT = 1;

    if (fr->enc == NULL) {
      fr->enc = (char *)safe_malloc(sizeof(char) * AS_READ_MAX_LONG_LEN + 1);
      fr->seq = (char *)safe_malloc(sizeof(char) * AS_READ_MAX_LONG_LEN + 1);
      fr->qlt = (char *)safe_malloc(sizeof(char) * AS_READ_MAX_LONG_LEN + 1);
    }

    if (partmap) {
      store = (fr->type == GKFRAGMENT_MEDIUM) ? partqmd : partqlg;
      getStringStore(store, qltOff, fr->enc, AS_READ_MAX_LONG_LEN, &actLen, &nxtOff);
    } else if (gst == NULL) {
      store = (fr->type == GKFRAGMENT_MEDIUM) ? qmd : qlg;
      getStringStore(store, qltOff, fr->enc, AS_READ_MAX_LONG_LEN, &actLen, &nxtOff);
    } else {
      stream = (fr->type == GKFRAGMENT_MEDIUM) ? gst->qmd : gst->qlg;
      nextStream(stream, fr->enc, AS_READ_MAX_LONG_LEN, &actLen);
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

  gkStore_decodeTypeFromIID(iid, &fr->type, &fr->tiid);

  //fprintf(stderr, "gkStore_getFragment()--  Retrieving IID=%d from %d,%d\n", iid, fr->type, fr->tiid);

  if (partmap) {
    //  If partitioned, we have everything in memory.  This is keyed
    //  off of the global IID.
    void *p = (void *)(INTPTR)LookupValueInHashTable_AS(partmap, iid, 0);

    if (p == NULL)
      fprintf(stderr, "getFrag()-- ERROR!  IID "F_IID" not in partition!\n", iid);
    assert(p != NULL);

    switch (fr->type) {
      case GKFRAGMENT_SHORT:
        memcpy(&fr->fr.sm, p, sizeof(gkShortFragment));
        break;
      case GKFRAGMENT_MEDIUM:
        memcpy(&fr->fr.md, p, sizeof(gkMediumFragment));
        break;
      case GKFRAGMENT_LONG:
        memcpy(&fr->fr.lg, p, sizeof(gkLongFragment));
        break;
    }

  } else {
    //  Not paritioned, load from disk store.
    switch (fr->type) {
      case GKFRAGMENT_SHORT:
        getIndexStore(fsm, fr->tiid, &fr->fr.sm);
        break;
      case GKFRAGMENT_MEDIUM:
        getIndexStore(fmd, fr->tiid, &fr->fr.md);
        break;
      case GKFRAGMENT_LONG:
        getIndexStore(flg, fr->tiid, &fr->fr.lg);
        break;
    }
  }

  //fprintf(stderr, "gkStore_getFragment()--  Retrieved IID=%d (asked for %d) from %d,%d\n", fr->gkFragment_getReadIID(), iid, fr->type, fr->tiid);

  gkStore_getFragmentData(NULL, fr, flags);
}




void
gkStore::gkStore_setFragment(gkFragment *fr) {
  assert(partmap    == NULL);
  assert(isReadOnly == 0);

  uint32  type;
  uint32  tiid;

  gkStore_decodeTypeFromIID(fr->gkFragment_getReadIID(), &type, &tiid);

  assert(type == fr->type);
  assert(tiid == fr->tiid);

  //fprintf(stderr, "gkStore_setFragment()--  Setting IID=%d to %d,%d\n", fr->gkFragment_getReadIID(), fr->type, fr->tiid);

  switch (fr->type) {
    case GKFRAGMENT_SHORT:
      setIndexStore(fsm, fr->tiid, &fr->fr.sm);
      break;
    case GKFRAGMENT_MEDIUM:
      setIndexStore(fmd, fr->tiid, &fr->fr.md);
      break;
    case GKFRAGMENT_LONG:
      setIndexStore(flg, fr->tiid, &fr->fr.lg);
      break;
  }
}




void
gkStore::gkStore_delFragment(AS_IID iid) {
  gkFragment   fr;
  AS_IID       mateiid;

  assert(partmap    == NULL);
  assert(isReadOnly == 0);

  //  Delete fragment with iid from the store.  If the fragment has a
  //  mate, remove the mate relationship from both fragmentss.

  gkStore_getFragment(iid, &fr, GKFRAGMENT_INF);

  mateiid = fr.gkFragment_getMateIID();

  switch (fr.type) {
    case GKFRAGMENT_SHORT:
      fr.fr.sm.deleted = 1;
      fr.fr.sm.mateIID = 0;
      break;
    case GKFRAGMENT_MEDIUM:
      fr.fr.md.deleted = 1;
      fr.fr.md.mateIID = 0;
      break;
    case GKFRAGMENT_LONG:
      fr.fr.lg.deleted = 1;
      fr.fr.lg.mateIID = 0;
      break;
  }
  gkStore_setFragment(&fr);

  if (mateiid > 0) {
    gkStore_getFragment(mateiid, &fr, GKFRAGMENT_INF);
    switch (fr.type) {
      case GKFRAGMENT_SHORT:
        fr.fr.sm.mateIID = 0;
        break;
      case GKFRAGMENT_MEDIUM:
        fr.fr.md.mateIID = 0;
        break;
      case GKFRAGMENT_LONG:
        fr.fr.lg.mateIID = 0;
        break;
    }
    gkStore_setFragment(&fr);
  }
}



void
gkStore::gkStore_addFragment(gkFragment *fr,
                             uint32 clrBgn, uint32 clrEnd,
                             uint32 vecBgn, uint32 vecEnd,
                             uint32 maxBgn, uint32 maxEnd,
                             uint32 tntBgn, uint32 tntEnd) {
  int encLen;

  assert(partmap    == NULL);
  assert(isReadOnly == 0);

  assert(fr->type != GKFRAGMENT_ERROR);

  if (clrBgn < clrEnd)  clearRange[AS_READ_CLEAR_CLR]->gkClearRange_enableCreate();
  if (vecBgn < vecEnd)  clearRange[AS_READ_CLEAR_VEC]->gkClearRange_enableCreate();
  if (maxBgn < maxEnd)  clearRange[AS_READ_CLEAR_MAX]->gkClearRange_enableCreate();
  if (tntBgn < tntEnd)  clearRange[AS_READ_CLEAR_TNT]->gkClearRange_enableCreate();

  switch (fr->type) {
    case GKFRAGMENT_SHORT:
      fr->tiid          = ++inf.numShort;
      fr->fr.sm.readIID = gkStore_getNumFragments() + 1;

      assert(fr->tiid == getLastElemStore(fsm) + 1);

      if (clrBgn < clrEnd)  clearRange[AS_READ_CLEAR_CLR]->gkClearRange_makeSpaceShort(fr->tiid);
      if (vecBgn < vecEnd)  clearRange[AS_READ_CLEAR_VEC]->gkClearRange_makeSpaceShort(fr->tiid);
      if (maxBgn < maxEnd)  clearRange[AS_READ_CLEAR_MAX]->gkClearRange_makeSpaceShort(fr->tiid);
      if (tntBgn < tntEnd)  clearRange[AS_READ_CLEAR_TNT]->gkClearRange_makeSpaceShort(fr->tiid);

      fr->fr.sm.clearBeg = clrBgn;
      fr->fr.sm.clearEnd = clrEnd;
      break;

    case GKFRAGMENT_MEDIUM:
      fr->tiid          = ++inf.numMedium;
      fr->fr.md.readIID = gkStore_getNumFragments() + 1;

      assert(fr->tiid == getLastElemStore(fmd) + 1);

      if (clrBgn < clrEnd)  clearRange[AS_READ_CLEAR_CLR]->gkClearRange_makeSpaceMedium(fr->tiid);
      if (vecBgn < vecEnd)  clearRange[AS_READ_CLEAR_VEC]->gkClearRange_makeSpaceMedium(fr->tiid);
      if (maxBgn < maxEnd)  clearRange[AS_READ_CLEAR_MAX]->gkClearRange_makeSpaceMedium(fr->tiid);
      if (tntBgn < tntEnd)  clearRange[AS_READ_CLEAR_TNT]->gkClearRange_makeSpaceMedium(fr->tiid);

      fr->fr.md.clearBeg = clrBgn;
      fr->fr.md.clearEnd = clrEnd;
      break;

    case GKFRAGMENT_LONG:
      fr->tiid           = ++inf.numLong;
      fr->fr.lg.readIID = gkStore_getNumFragments() + 1;

      assert(fr->tiid == getLastElemStore(flg) + 1);

      if (clrBgn < clrEnd)  clearRange[AS_READ_CLEAR_CLR]->gkClearRange_makeSpaceLong(fr->tiid);
      if (vecBgn < vecEnd)  clearRange[AS_READ_CLEAR_VEC]->gkClearRange_makeSpaceLong(fr->tiid);
      if (maxBgn < maxEnd)  clearRange[AS_READ_CLEAR_MAX]->gkClearRange_makeSpaceLong(fr->tiid);
      if (tntBgn < tntEnd)  clearRange[AS_READ_CLEAR_TNT]->gkClearRange_makeSpaceLong(fr->tiid);

      fr->fr.lg.clearBeg = clrBgn;
      fr->fr.lg.clearEnd = clrEnd;
      break;
  }

  assert(fr->tiid != 0);

  //  Set clear ranges...if they're defined.  !!NOTE!!  The "CLR"
  //  range MUST be set last since it is "the" clear range, and
  //  setClearRange() always updates the latest clear range along with
  //  the named region.  (That is, setting VEC will set both VEC and
  //  LATEST.)
  //
  if (vecBgn < vecEnd)  clearRange[AS_READ_CLEAR_VEC]->gkClearRange_setClearRegion(fr, vecBgn, vecEnd);
  if (maxBgn < maxEnd)  clearRange[AS_READ_CLEAR_MAX]->gkClearRange_setClearRegion(fr, maxBgn, maxEnd);
  if (tntBgn < tntEnd)  clearRange[AS_READ_CLEAR_TNT]->gkClearRange_setClearRegion(fr, tntBgn, tntEnd);
  if (clrBgn < clrEnd)  clearRange[AS_READ_CLEAR_CLR]->gkClearRange_setClearRegion(fr, clrBgn, clrEnd);

  switch (fr->type) {
    case GKFRAGMENT_SHORT:
      encodeSequenceQuality(fr->fr.sm.enc, fr->seq, fr->qlt);

      gkStore_setUIDtoIID(fr->fr.sm.readUID, fr->fr.sm.readIID, AS_IID_FRG);
      appendIndexStore(fsm, &fr->fr.sm);
      break;

    case GKFRAGMENT_MEDIUM:
      fr->fr.md.seqOffset = getLastElemStore(smd) + 1;
      fr->fr.md.qltOffset = getLastElemStore(qmd) + 1;

      gkStore_setUIDtoIID(fr->fr.md.readUID, fr->fr.md.readIID, AS_IID_FRG);
      appendIndexStore(fmd, &fr->fr.md);

      encLen = encodeSequence(fr->enc, fr->seq);
      appendStringStore(smd, fr->enc, encLen);

      encodeSequenceQuality(fr->enc, fr->seq, fr->qlt);
      appendStringStore(qmd, fr->enc, fr->fr.md.seqLen);
      break;

    case GKFRAGMENT_LONG:
      fr->fr.lg.seqOffset = getLastElemStore(slg) + 1;
      fr->fr.lg.qltOffset = getLastElemStore(qlg) + 1;

      gkStore_setUIDtoIID(fr->fr.lg.readUID, fr->fr.lg.readIID, AS_IID_FRG);
      appendIndexStore(flg, &fr->fr.lg);

      encLen = encodeSequence(fr->enc, fr->seq);
      appendStringStore(slg, fr->enc, encLen);

      encodeSequenceQuality(fr->enc, fr->seq, fr->qlt);
      appendStringStore(qlg, fr->enc, fr->fr.lg.seqLen);
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
    if (p == -1)
      continue;


    if (fr.type == GKFRAGMENT_SHORT) {
      appendIndexStore(gkpart[p]->partfsm, &fr.fr.sm);
    }


    if (fr.type == GKFRAGMENT_MEDIUM) {
      fr.fr.md.seqOffset = -1;
      fr.fr.md.qltOffset = getLastElemStore(gkpart[p]->partqmd) + 1;

      appendIndexStore(gkpart[p]->partfmd, &fr.fr.md);

      encodeSequenceQuality(fr.enc, fr.seq, fr.qlt);
      appendStringStore(gkpart[p]->partqmd, fr.enc, fr.fr.md.seqLen);
    }


    if (fr.type == GKFRAGMENT_LONG) {
      fr.fr.lg.seqOffset = -1;
      fr.fr.lg.qltOffset = getLastElemStore(gkpart[p]->partqlg) + 1;

      appendIndexStore(gkpart[p]->partflg, &fr.fr.lg);

      encodeSequenceQuality(fr.enc, fr.seq, fr.qlt);
      appendStringStore(gkpart[p]->partqlg, fr.enc, fr.fr.lg.seqLen);
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
