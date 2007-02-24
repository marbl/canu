
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

static char CM_ID[] = "$Id: AS_PER_gkpStore.c,v 1.25 2007-02-24 15:42:33 brianwalenz Exp $";

//    A thin layer on top of the IndexStore supporing the storage and
// retrieval of records used by the gatekeeper records.
//
//    The idea is to provide easier to use shortcuts for the common
// operations, and let the other operations be accessed through the
// generic Index Store API.

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <unistd.h>

#include "AS_global.h"
#include "AS_PER_genericStore.h"
#include "AS_PER_gkpStore.h"
#include "AS_PER_encodeSequenceQuality.h"


static
int
fileExists(const char *path,
           int directory,
           int readwrite) {
  struct stat  s;
  int          r;

  errno = 0;
  r = stat(path, &s);
  if (errno) {
    //fprintf(stderr, "Failed to stat() file '%s'; assumed to not exist.\n", path);
    return(0);
  }

  if (directory == 1) {
    if ((readwrite == 0) &&
        (s.st_mode & S_IFDIR) &&
        (s.st_mode & (S_IRUSR | S_IRGRP | S_IROTH)) &&
        (s.st_mode & (S_IXUSR | S_IXGRP | S_IXOTH)))
      return(1);
    if ((readwrite == 1) &&
        (s.st_mode & S_IFDIR) &&
        (s.st_mode & (S_IRUSR | S_IRGRP | S_IROTH)) &&
        (s.st_mode & (S_IWUSR | S_IWGRP | S_IWOTH)) &&
        (s.st_mode & (S_IXUSR | S_IXGRP | S_IXOTH)))
      return(1);
  }

  if (directory == 0) {
    if ((readwrite == 0) &&
        (s.st_mode & (S_IRUSR | S_IRGRP | S_IROTH)))
      return(1);
    if ((readwrite == 1) &&
        (s.st_mode & (S_IRUSR | S_IRGRP | S_IROTH)) &&
        (s.st_mode & (S_IWUSR | S_IWGRP | S_IWOTH)))
      return(1);
  }

  return(0);
}


int
testOpenGateKeeperStore(const char *path,
                        int   writable) {

  char name[FILENAME_MAX];
  int  fileCount = 0;

  if (fileExists(path, 1, writable)) {
    sprintf(name,"%s/gkp", path);
    fileCount += fileExists(name, 0, writable);

    sprintf(name,"%s/bat", path);
    fileCount += fileExists(name, 0, writable);

    sprintf(name,"%s/frg", path);
    fileCount += fileExists(name, 0, writable);

    sprintf(name,"%s/lib", path);
    fileCount += fileExists(name, 0, writable);

    sprintf(name,"%s/lis", path);
    fileCount += fileExists(name, 0, writable);

    sprintf(name,"%s/seq", path);
    fileCount += fileExists(name, 0, writable);

    sprintf(name,"%s/qlt", path);
    fileCount += fileExists(name, 0, writable);

    sprintf(name,"%s/hps", path);
    fileCount += fileExists(name, 0, writable);

    sprintf(name,"%s/src", path);
    fileCount += fileExists(name, 0, writable);

    sprintf(name,"%s/phs", path);
    fileCount += fileExists(name, 0, writable);

  }

  return(fileCount == NUM_GKP_FILES);
}


GateKeeperStore *
openGateKeeperStore(const char *path,
                    int   writable) {

  char              name[FILENAME_MAX];
  FILE             *gkpinfo;

  GateKeeperStore  *gkpStore = (GateKeeperStore *)safe_calloc(1, sizeof(GateKeeperStore));

  strcpy(gkpStore->storePath, path);

  gkpStore->phs = 0L;
  gkpStore->bat = 0;
  gkpStore->frg = 0;
  gkpStore->lib = 0;
  gkpStore->lis = 0;

  gkpStore->seq = 0;
  gkpStore->qlt = 0;
  gkpStore->hps = 0;
  gkpStore->src = 0;
  

  sprintf(name,"%s/gkp", gkpStore->storePath);
  errno = 0;
  gkpinfo = fopen(name, "r");
  if (errno) {
    fprintf(stderr, "failed to open gatekeeper store '%s': %s\n", name, strerror(errno));
    exit(1);
  }

  fread(&gkpStore->gkp, sizeof(GateKeeperStoreInfo), 1, gkpinfo);
  fclose(gkpinfo);

  if (gkpStore->gkp.gkpMagic != 1) {
    fprintf(stderr, "invalid magic!\n");
  }
  if (gkpStore->gkp.gkpVersion != 1) {
    fprintf(stderr, "invalid version!\n");
  }
  if (gkpStore->gkp.gkpBatchRecordSize != sizeof(GateKeeperBatchRecord)) {
    fprintf(stderr, "ERROR!  Store built unsing GateKeeperBatchRecord of size %d bytes.\n", gkpStore->gkp.gkpBatchRecordSize);
    fprintf(stderr, "        Code compiled with GateKeeperBatchRecord of size %d bytes.\n", (int)sizeof(GateKeeperBatchRecord));
    assert(gkpStore->gkp.gkpBatchRecordSize == sizeof(GateKeeperBatchRecord));
  }
  if (gkpStore->gkp.gkpLibraryRecordSize != sizeof(GateKeeperLibraryRecord)) {
    fprintf(stderr, "ERROR!  Store built unsing GateKeeperLibraryRecord of size %d bytes.\n", gkpStore->gkp.gkpLibraryRecordSize);
    fprintf(stderr, "        Code compiled with GateKeeperLibraryRecord of size %d bytes.\n", (int)sizeof(GateKeeperLibraryRecord));
    assert(gkpStore->gkp.gkpLibraryRecordSize == sizeof(GateKeeperLibraryRecord));
  }
  if (gkpStore->gkp.gkpFragmentRecordSize != sizeof(GateKeeperFragmentRecord)) {
    fprintf(stderr, "ERROR!  Store built unsing GateKeeperFragmentRecord of size %d bytes.\n", gkpStore->gkp.gkpFragmentRecordSize);
    fprintf(stderr, "        Code compiled with GateKeeperFragmentRecord of size %d bytes.\n", (int)sizeof(GateKeeperFragmentRecord));
    assert(gkpStore->gkp.gkpFragmentRecordSize == sizeof(GateKeeperFragmentRecord));
  }


  char  mode[4];
  if (writable)
    strcpy(mode, "r+");
  else
    strcpy(mode, "r");

  sprintf(name,"%s/bat", gkpStore->storePath);
  gkpStore->bat   = openGateKeeperBatchStore(name, mode);

  sprintf(name,"%s/frg", gkpStore->storePath);
  gkpStore->frg   = openGateKeeperFragmentStore(name, mode);

  sprintf(name,"%s/lib", gkpStore->storePath);
  gkpStore->lib   = openGateKeeperLibraryStore(name, mode);

  sprintf(name,"%s/lis", gkpStore->storePath);
  gkpStore->lis = openGateKeeperLibraryStore(name, mode);

  sprintf(name,"%s/seq", gkpStore->storePath);
  gkpStore->seq = openStore(name, mode);

  sprintf(name,"%s/qlt", gkpStore->storePath);
  gkpStore->qlt = openStore(name, mode);

  sprintf(name,"%s/hps", gkpStore->storePath);
  gkpStore->hps = openStore(name, mode);

  sprintf(name,"%s/src", gkpStore->storePath);
  gkpStore->src = openStore(name, mode);

  if ((NULLSTOREHANDLE == gkpStore->bat) ||
      (NULLSTOREHANDLE == gkpStore->frg) ||
      (NULLSTOREHANDLE == gkpStore->lib) ||
      (NULLSTOREHANDLE == gkpStore->lis) ||
      (NULLSTOREHANDLE == gkpStore->seq) ||
      (NULLSTOREHANDLE == gkpStore->qlt) ||
      (NULLSTOREHANDLE == gkpStore->hps) ||
      (NULLSTOREHANDLE == gkpStore->src)) {
    fprintf(stderr,"**** Failure to open Gatekeeper Store ...\n");
    assert(0);
  }

  sprintf(name,"%s/phs", gkpStore->storePath);

  if (mode && mode[0] == 'r' && mode[1] == '\0') {
    gkpStore->phs = OpenReadOnlyPHashTable_AS(name);
  } else {
    gkpStore->phs = OpenPHashTable_AS(name);
  }

  if (gkpStore->phs == NULL) {
    fprintf(stderr,"**** Failed to open GateKeeper Persistent HashTable...\n");
    assert(0);
  }

  return(gkpStore);
}


GateKeeperStore *
createGateKeeperStore(const char *path) {
  char   name[FILENAME_MAX];
  FILE  *gkpinfo;

  GateKeeperStore  *gkpStore = (GateKeeperStore *)safe_calloc(1, sizeof(GateKeeperStore));

  strcpy(gkpStore->storePath, path);

  gkpStore->phs = 0L;
  gkpStore->bat = 0;
  gkpStore->frg = 0;
  gkpStore->lib = 0;
  gkpStore->lis = 0;

  gkpStore->seq = 0;
  gkpStore->qlt = 0;
  gkpStore->hps = 0;
  gkpStore->src = 0;

  errno = 0;
  mkdir(path, S_IRWXU | S_IRWXG | S_IROTH);
  if (errno) {
    fprintf(stderr, "CreateGateKeeperStore(): failed to create directory '%s': %s\n", gkpStore->storePath, strerror(errno));
    exit(1);
  }

  gkpStore->gkp.gkpMagic              = 1;
  gkpStore->gkp.gkpVersion            = 1;
  gkpStore->gkp.gkpBatchRecordSize    = sizeof(GateKeeperBatchRecord);
  gkpStore->gkp.gkpLibraryRecordSize  = sizeof(GateKeeperLibraryRecord);
  gkpStore->gkp.gkpFragmentRecordSize = sizeof(GateKeeperFragmentRecord);

  sprintf(name,"%s/gkp", path);
  errno = 0;
  gkpinfo = fopen(name, "w");
  if (errno) {
    fprintf(stderr, "failed to create gatekeeper store '%s': %s\n", name, strerror(errno));
    exit(1);
  }

  fwrite(&gkpStore->gkp, sizeof(GateKeeperStoreInfo), 1, gkpinfo);
  fclose(gkpinfo);

  sprintf(name,"%s/bat", path);
  gkpStore->bat = createGateKeeperBatchStore(name, "bat", 1);

  sprintf(name,"%s/frg", path);
  gkpStore->frg = createGateKeeperFragmentStore(name, "frg", 1);

  sprintf(name,"%s/lib", path);
  gkpStore->lib = createGateKeeperLibraryStore(name, "lib", 1);

  sprintf(name,"%s/lis", path);
  gkpStore->lis = createGateKeeperLibraryStore(name, "lib", 1);

  sprintf(name,"%s/seq", path);
  gkpStore->seq = createVLRecordStore(name, "seq", MAX_SEQ_LENGTH, 1);

  sprintf(name,"%s/qlt", path);
  gkpStore->qlt = createVLRecordStore(name, "qlt", MAX_SEQ_LENGTH, 1);

  sprintf(name,"%s/hps", path);
  gkpStore->hps = createVLRecordStore(name, "hps", MAX_HPS_LENGTH, 1);

  sprintf(name,"%s/src", path);
  gkpStore->src = createVLRecordStore(name, "src", MAX_SRC_LENGTH, 1);

  sprintf(name,"%s/phs", path);
  gkpStore->phs = CreatePHashTable_AS(2048,name);

  return(gkpStore);
}


void
closeGateKeeperStore(GateKeeperStore *gkpStore) {

  if(gkpStore->bat != NULLSTOREHANDLE)
    closeStore(gkpStore->bat);

  if(gkpStore->frg != NULLSTOREHANDLE)
    closeStore(gkpStore->frg);

  if(gkpStore->lib != NULLSTOREHANDLE)
    closeStore(gkpStore->lib);

  if(gkpStore->lis != NULLSTOREHANDLE)
    closeStore(gkpStore->lis);

  if(gkpStore->seq != NULLSTOREHANDLE)
    closeStore(gkpStore->seq);

  if(gkpStore->qlt != NULLSTOREHANDLE)
    closeStore(gkpStore->qlt);

  if(gkpStore->hps != NULLSTOREHANDLE)
    closeStore(gkpStore->hps);

  if(gkpStore->src != NULLSTOREHANDLE)
    closeStore(gkpStore->src);

  if(gkpStore->phs != NULL)
    ClosePHashTable_AS(gkpStore->phs);

  safe_free(gkpStore);
}




void clearGateKeeperBatchRecord(GateKeeperBatchRecord *g) {
  memset(g, 0, sizeof(GateKeeperBatchRecord));
}

void clearGateKeeperLibraryRecord(GateKeeperLibraryRecord *g) {
  memset(g, 0, sizeof(GateKeeperLibraryRecord));
}

void clearGateKeeperFragmentRecord(GateKeeperFragmentRecord *g) {
  memset(g, 0, sizeof(GateKeeperFragmentRecord));
}




////////////////////////////////////////////////////////////////////////////////



fragRecord *new_fragRecord(void) {
  fragRecord *fr = (fragRecord *)safe_malloc(sizeof(fragRecord));
  clr_fragRecord(fr);
  return(fr);
}

void        del_fragRecord(fragRecord *fr) {
  safe_free(fr);
}

void        clr_fragRecord(fragRecord *fr) {
  clearGateKeeperFragmentRecord(&fr->gkfr);
  fr->seq[0] = 0;
  fr->qlt[0] = 0;
  fr->hps[0] = 0;
  fr->src[0] = 0;
}

//  Set clear region 'which' and all later clear regions.  You are
//  explicitly not allowed to set the original clear range.
//
void        setFragRecordClearRegion(fragRecord *fr,
                                     uint32 start,
                                     uint32 end,
                                     uint32 which) {
  assert(which > AS_READ_CLEAR_VEC);
  assert(which < AS_READ_CLEAR_NUM);
  for (; which < AS_READ_CLEAR_NUM; which++) {
    fr->gkfr.clearBeg[which] = start;
    fr->gkfr.clearEnd[which] = end;
  }
}


void        getFragRecordClearRegion(fragRecord *fr, uint32 *start, uint32 *end, uint32 which) {
  assert(which <  AS_READ_CLEAR_NUM);
  *start = fr->gkfr.clearBeg[which];
  *end   = fr->gkfr.clearEnd[which];
}


uint32      getFragRecordClearRegionBegin(fragRecord *fr, uint32 which) {
  assert(which < AS_READ_CLEAR_NUM);
  return(fr->gkfr.clearBeg[which]);
}


uint32      getFragRecordClearRegionEnd  (fragRecord *fr, uint32 which) {
  assert(which < AS_READ_CLEAR_NUM);
  return(fr->gkfr.clearEnd[which]);
}



////////////////////////////////////////////////////////////////////////////////


static
void
getFragData(GateKeeperStore *gkp, fragRecord *fr, int streamFlags) {
  VLSTRING_SIZE_T   actualLength = 0;

  fr->hasSEQ = 0;
  fr->hasQLT = 0;
  fr->hasHPS = 0;
  fr->hasSRC = 0;

  fr->seq[0] = 0;
  fr->qlt[0] = 0;
  fr->hps[0] = 0;
  fr->src[0] = 0;

  if ((streamFlags & FRAG_S_SEQ) &&
      !(streamFlags & FRAG_S_QLT)) {
    fr->hasSEQ = 1;
    if (fr->gkfr.seqLen > 0) {
      getVLRecordStore(gkp->seq,
                       fr->gkfr.seqOffset,
                       fr->seq,
                       VLSTRING_MAX_SIZE,
                       &actualLength);
      fr->seq[actualLength] = 0;
    }
  }
  if (streamFlags & FRAG_S_QLT) {
    assert(fr->hasSEQ == 0);
    fr->hasSEQ = 1;
    fr->hasQLT = 1;
    if (fr->gkfr.seqLen > 0) {
      getVLRecordStore(gkp->qlt,
                       fr->gkfr.qltOffset,
                       fr->qlt,
                       VLSTRING_MAX_SIZE,
                       &actualLength);
      fr->qlt[actualLength] = 0;
      decodeSequenceQuality(fr->qlt, fr->seq, fr->qlt);
    }
  }
  if (streamFlags & FRAG_S_HPS) {
    fr->hasHPS = 1;
    if (fr->gkfr.hpsLen > 0) {
      getVLRecordStore(gkp->hps,
                       fr->gkfr.hpsOffset,
                       fr->hps,
                       VLSTRING_MAX_SIZE,
                       &actualLength);
      fr->hps[actualLength] = 0;
    }
  }
  if (streamFlags & FRAG_S_SRC) {
    fr->hasSRC = 1;
    if (fr->gkfr.srcLen > 0) {
      getVLRecordStore(gkp->src,
                       fr->gkfr.srcOffset,
                       fr->src,
                       VLSTRING_MAX_SIZE,
                       &actualLength);
      fr->src[actualLength] = 0;
    }
  }
}


void    getFrag(GateKeeperStore *gkp, int64 iid, fragRecord *fr, int32 flags) {
  getGateKeeperFragmentStore(gkp->frg, iid, &fr->gkfr);
  getFragData(gkp, fr, flags);
}

void    setFrag(GateKeeperStore *gkp, int64 iid, fragRecord *fr) {
  setGateKeeperFragmentStore(gkp->frg, iid, &fr->gkfr);
}

void    delFrag(GateKeeperStore *gkp, int64 iid) {
  GateKeeperFragmentRecord   gkfr;
  CDS_IID_t                  miid;

  //  Delete fragment with iid from the store.  If the fragment has a
  //  mate, remove the mate relationship from both fragmentss.

  getGateKeeperFragmentStore(gkp->frg, iid, &gkfr);
  miid = gkfr.mateIID;
  gkfr.deleted = 1;
  gkfr.mateIID = 0;
  setGateKeeperFragmentStore(gkp->frg, iid, &gkfr);

  if (miid > 0) {
    UnRefPHashTable_AS(gkp->phs, UID_NAMESPACE_AS, gkfr.readUID);

    getGateKeeperFragmentStore(gkp->frg, miid, &gkfr);
    gkfr.mateIID = 0;
    setGateKeeperFragmentStore(gkp->frg, miid, &gkfr);

    UnRefPHashTable_AS(gkp->phs, UID_NAMESPACE_AS, gkfr.readUID);
  }
}



////////////////////////////////////////////////////////////////////////////////


FragStream      *openFragStream(GateKeeperStore *gkp, int flags) {
  FragStream  *fs = (FragStream *)safe_malloc(sizeof(FragStream));
  fs->gkp   = gkp;
  fs->frg   = NULLSTOREHANDLE;
  fs->seq   = NULLSTOREHANDLE;
  fs->qlt   = NULLSTOREHANDLE;
  fs->hps   = NULLSTOREHANDLE;
  fs->src   = NULLSTOREHANDLE;
  fs->flags = flags;

#if 0
  int bufferSize = 1048576;
  fs->frgBuffer = (char *)safe_malloc(sizeof(char) * bufferSize);
  fs->seqBuffer = (char *)safe_malloc(sizeof(char) * bufferSize);
  fs->qltBuffer = (char *)safe_malloc(sizeof(char) * bufferSize);
  fs->hpsBuffer = (char *)safe_malloc(sizeof(char) * bufferSize);
  fs->seqBuffer = (char *)safe_malloc(sizeof(char) * bufferSize);
#else
  int bufferSize = 0;
  fs->frgBuffer = NULL;
  fs->seqBuffer = NULL;
  fs->qltBuffer = NULL;
  fs->hpsBuffer = NULL;
  fs->seqBuffer = NULL;
#endif

  fs->frg   = openStream(fs->gkp->frg, fs->frgBuffer, bufferSize);

  if ((fs->flags & FRAG_S_SEQ) && !(fs->flags & FRAG_S_QLT))
    fs->seq   = openStream(fs->gkp->seq, fs->seqBuffer, bufferSize);
  if (fs->flags & FRAG_S_QLT)
    fs->qlt   = openStream(fs->gkp->qlt, fs->qltBuffer, bufferSize);
  if (fs->flags & FRAG_S_HPS)
    fs->hps   = openStream(fs->gkp->hps, fs->hpsBuffer, bufferSize);
  if (fs->flags & FRAG_S_SRC)
    fs->src   = openStream(fs->gkp->src, fs->srcBuffer, bufferSize);

  return(fs);
}


void             resetFragStream(FragStream *fs, int64 startIndex, int64 endIndex) {
  int64  seqOffset, qltOffset, hpsOffset, srcOffset;

  if (startIndex == STREAM_FROMSTART) {
    seqOffset = STREAM_FROMSTART;
    qltOffset = STREAM_FROMSTART;
    hpsOffset = STREAM_FROMSTART;
    srcOffset = STREAM_FROMSTART;
  } else {

    //  Read the frag the hard way so we can find the offsets the other
    //  stores need to use.

    GateKeeperFragmentRecord  gkpf;

    getGateKeeperFragmentStore(fs->gkp->frg, startIndex, &gkpf);

    seqOffset = gkpf.seqOffset;
    qltOffset = gkpf.qltOffset;
    hpsOffset = gkpf.hpsOffset;
    srcOffset = gkpf.srcOffset;
  }

  resetStream(fs->frg, startIndex, endIndex);

  if ((fs->flags & FRAG_S_SEQ) && !(fs->flags & FRAG_S_QLT))
    resetStream(fs->seq, seqOffset, STREAM_UNTILEND);
  if (fs->flags & FRAG_S_QLT)
    resetStream(fs->qlt, qltOffset, STREAM_UNTILEND);
  if (fs->flags & FRAG_S_HPS)
    resetStream(fs->hps, hpsOffset, STREAM_UNTILEND);
  if (fs->flags & FRAG_S_SRC)
    resetStream(fs->src, srcOffset, STREAM_UNTILEND);
}


void             closeFragStream(FragStream *fs) {
  closeStream(fs->frg);
  if ((fs->flags & FRAG_S_SEQ) && !(fs->flags & FRAG_S_QLT))
    closeStream(fs->seq);
  if (fs->flags & FRAG_S_QLT)
    closeStream(fs->qlt);
  if (fs->flags & FRAG_S_HPS)
    closeStream(fs->hps);
  if (fs->flags & FRAG_S_SRC)
    closeStream(fs->src);

  safe_free(fs->frgBuffer);
  safe_free(fs->seqBuffer);
  safe_free(fs->qltBuffer);
  safe_free(fs->hpsBuffer);
  safe_free(fs->seqBuffer);
}


int64            getStartIndexFragStream(FragStream *fs) {
  return(getStartIndexStream(fs->frg));
}


int              nextFragStream(FragStream *fs, fragRecord *fr) {
  VLSTRING_SIZE_T  actualLength = 0;

  if (nextStream(fs->frg, &fr->gkfr) == 0)
    return(0);

  //  So we can use the stream, we can't use getFragData() here.  We
  //  need to duplicate it.

  fr->seq[0] = 0;
  fr->qlt[0] = 0;
  fr->hps[0] = 0;
  fr->src[0] = 0;

  fr->hasSEQ = 0;
  fr->hasQLT = 0;
  fr->hasHPS = 0;
  fr->hasSRC = 0;

  if ((fs->flags & FRAG_S_SEQ) && !(fs->flags & FRAG_S_QLT)) {
    fr->hasSEQ = 1;
    nextVLRecordStream(fs->seq, fr->seq, MAX_SEQ_LENGTH, &actualLength);
    fr->seq[actualLength] = 0;
  }

  if (fs->flags & FRAG_S_QLT) {
    assert(fr->hasSEQ == 0);
    fr->hasSEQ = 1;
    fr->hasQLT = 1;
    nextVLRecordStream(fs->qlt, fr->qlt, MAX_SEQ_LENGTH, &actualLength);
    fr->qlt[actualLength] = 0;
    decodeSequenceQuality(fr->qlt, fr->seq, fr->qlt);
  }

  if (fs->flags & FRAG_S_HPS) {
    fr->hasHPS = 1;
    nextVLRecordStream(fs->hps, fr->hps, MAX_HPS_LENGTH, &actualLength);
    fr->hps[actualLength] = 0;
  }

  if (fs->flags & FRAG_S_SRC) {
    fr->hasSRC = 1;
    nextVLRecordStream(fs->src, fr->src, MAX_SRC_LENGTH, &actualLength);
    fr->src[actualLength] = 0;
  }

  return(1);
}





GateKeeperStore *
loadFragStorePartial(const char *path,
                     int64       firstElem,
                     int64       lastElem,
                     int         flags) {
  GateKeeperStore           *gkp = openGateKeeperStore(path, FALSE);
  GateKeeperFragmentRecord   gkf;

  int64  seqFirst = 0, seqLast = STREAM_UNTILEND;
  int64  qltFirst = 0, qltLast = STREAM_UNTILEND;
  int64  hpsFirst = 0, hpsLast = STREAM_UNTILEND;
  int64  srcFirst = 0, srcLast = STREAM_UNTILEND;

  int64  lastFrag = getLastElemFragStore(gkp);

  //  Making the seq, qlt, etc, stores into memory stores is more
  //  trouble.

  getGateKeeperFragmentStore(gkp->frg, firstElem, &gkf);

  seqFirst = gkf.seqOffset;
  qltFirst = gkf.qltOffset;
  hpsFirst = gkf.hpsOffset;
  srcFirst = gkf.srcOffset;

  if ((lastElem != STREAM_UNTILEND) &&
      (lastElem != lastFrag)) {
    getGateKeeperFragmentStore(gkp->frg, lastElem + 1, &gkf);

    seqLast = gkf.seqOffset;
    qltLast = gkf.qltOffset;
    hpsLast = gkf.hpsOffset;
    srcLast = gkf.srcOffset;
  }

  if ((flags & FRAG_S_SEQ) && !(flags & FRAG_S_QLT))
    gkp->seq = convertStoreToPartialMemoryStore(gkp->seq, seqFirst, seqLast);

  if (flags & FRAG_S_QLT)
    gkp->qlt = convertStoreToPartialMemoryStore(gkp->qlt, qltFirst, qltLast);

  if (flags & FRAG_S_HPS)
    gkp->hps = convertStoreToPartialMemoryStore(gkp->hps, hpsFirst, hpsLast);

  if (flags & FRAG_S_SRC)
    gkp->src = convertStoreToPartialMemoryStore(gkp->src, srcFirst, srcLast);

  //  Making the frag info a memory store is easy!  We do this last,
  //  since we need to get a frag that isn't included here -- being
  //  lastElem+1.  The cost of doing this last is two fragment loads
  //  that otherwise we could do from the memory store.
  //
  gkp->frg = convertStoreToPartialMemoryStore(gkp->frg, firstElem, lastElem);

  return(gkp);
}
