
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

static char CM_ID[] = "$Id: AS_PER_gkpStore.c,v 1.18 2007-02-14 07:20:13 brianwalenz Exp $";

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
    return(0);
  }

  if (directory == 0) {
    if ((readwrite == 0) &&
        (s.st_mode & (S_IRUSR | S_IRGRP | S_IROTH)))
      return(1);
    if ((readwrite == 1) &&
        (s.st_mode & (S_IRUSR | S_IRGRP | S_IROTH)) &&
        (s.st_mode & (S_IWUSR | S_IWGRP | S_IWOTH)))
      return(1);
    return(0);
  }
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

  gkpStore->gkp.gkpMagic   = 1;
  gkpStore->gkp.gkpVersion = 1;

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

  fprintf(stderr, "WARNING!  USING 2048 AS LENGTH FOR VLRecordStore!\n");

  sprintf(name,"%s/seq", path);
  gkpStore->seq = createVLRecordStore(name, "seq", 2048, 1);

  sprintf(name,"%s/qlt", path);
  gkpStore->qlt = createVLRecordStore(name, "qlt", 2048, 1);

  sprintf(name,"%s/hps", path);
  gkpStore->hps = createVLRecordStore(name, "hps", 2048, 1);

  sprintf(name,"%s/src", path);
  gkpStore->src = createVLRecordStore(name, "src", 2048, 1);

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

  g->UID            = 0;
  g->name[0]        = 0;
  g->comment[0]     = 0;
  g->created        = 0;
  g->deleted        = 0;
  g->spare          = 0;
  g->numFragments   = 0;
  g->numLibraries   = 0;
  g->numLibraries_s = 0;

  memset(g->name,    0, AS_PER_NAME_LEN);
  memset(g->comment, 0, AS_PER_COMMENT_LEN);
}

void clearGateKeeperLibraryRecord(GateKeeperLibraryRecord *g) {

  memset(g, 0, sizeof(GateKeeperLibraryRecord));

  g->UID             = 0;
  g->name[0]         = 0;
  g->comment[0]      = 0;
  g->created         = 0;
  g->deleted         = 0;
  g->redefined       = 0;
  g->orientation     = 0;
  g->spare           = 0;
  g->mean            = 0.0;
  g->stddev          = 0.0;
  g->numFeatures     = 0;
  g->prevInstanceID  = 0;
  g->prevID          = 0;
  g->birthBatch      = 0;
  g->deathBatch      = 0;

  memset(g->name,    0, AS_PER_NAME_LEN);
  memset(g->comment, 0, AS_PER_COMMENT_LEN);
}

void clearGateKeeperFragmentRecord(GateKeeperFragmentRecord *g) {

  memset(g, 0, sizeof(GateKeeperFragmentRecord));

  g->UID           = 0;
  g->readIID       = 0;
  g->mateIID       = 0;
  g->libraryIID    = 0;
  g->plateUID      = 0;
  g->plateLocation = 0;

  g->deleted     = 0;
  g->nonrandom   = 0;
  g->status      = 0;
  g->hasQLT      = 0;
  g->hasHPS      = 0;
  g->hasOVLclr   = 0;
  g->hasCNSclr   = 0;
  g->hasCGWclr   = 0;
  g->orientation = 0;
  g->spare       = 0;

  g->clrSta = 0;
  g->clrEnd = 0;
  g->ovlSta = 0;
  g->ovlEnd = 0;
  g->cnsSta = 0;
  g->cnsEnd = 0;
  g->cgwSta = 0;
  g->cgwEnd = 0;

  g->seqOffset = 0;
  g->qltOffset = 0;
  g->hpsOffset = 0;
  g->srcOffset = 0;

  g->birthBatch = 0;
  g->deathBatch = 0;
}




////////////////////////////////////////////////////////////////////////////////



static
void
getFragData(GateKeeperStore *gkp, ReadStruct *rs, int streamFlags) {
  VLSTRING_SIZE_T   actualLength = 0;

  rs->seq[0] = 0;
  rs->qlt[0] = 0;
  rs->hps[0] = 0;
  rs->src[0] = 0;

  if (streamFlags & FRAG_S_SEQ) {
    getVLRecordStore(gkp->seq,
                     rs->gkfr.seqOffset,
                     rs->seq,
                     VLSTRING_MAX_SIZE,
                     &actualLength);
    rs->seq[actualLength] = 0;
  }
  if (streamFlags & FRAG_S_QLT) {
    getVLRecordStore(gkp->qlt,
                     rs->gkfr.qltOffset,
                     rs->qlt,
                     VLSTRING_MAX_SIZE,
                     &actualLength);
    rs->qlt[actualLength] = 0;
  }
  if (streamFlags & FRAG_S_HPS) {
    getVLRecordStore(gkp->hps,
                     rs->gkfr.hpsOffset,
                     rs->hps,
                     VLSTRING_MAX_SIZE,
                     &actualLength);
    rs->hps[actualLength] = 0;
  }
  if (streamFlags & FRAG_S_SRC) {
    getVLRecordStore(gkp->src,
                     rs->gkfr.srcOffset,
                     rs->src,
                     VLSTRING_MAX_SIZE,
                     &actualLength);
    rs->src[actualLength] = 0;
  }
}


int     getFrag(GateKeeperStore *gkp, int64 iid, ReadStruct *rs, int32 flags) {
  getGateKeeperFragmentStore(gkp->frg, iid, &rs->gkfr);
  getFragData(gkp, rs, flags);
  return(0);
}

int     setFrag(GateKeeperStore *gkp, int64 iid, ReadStruct *rs) {
  setGateKeeperFragmentStore(gkp->frg, iid, &rs->gkfr);
  return(0);
}

int     delFrag(GateKeeperStore *gkp, int64 iid) {
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
    UnRefPHashTable_AS(gkp->phs, UID_NAMESPACE_AS, gkfr.UID);

    getGateKeeperFragmentStore(gkp->frg, miid, &gkfr);
    gkfr.mateIID = 0;
    setGateKeeperFragmentStore(gkp->frg, miid, &gkfr);

    UnRefPHashTable_AS(gkp->phs, UID_NAMESPACE_AS, gkfr.UID);
  }
}



////////////////////////////////////////////////////////////////////////////////


FragStream      *openFragStream(GateKeeperStore *gkp) {
  FragStream  *fs = (FragStream *)safe_malloc(sizeof(FragStream));
  fs->gkp = gkp;
  fs->frg = openStream(fs->gkp->frg, NULL, 0);
  fs->seq = openStream(fs->gkp->seq, NULL, 0);
  fs->qlt = openStream(fs->gkp->qlt, NULL, 0);
  fs->hps = openStream(fs->gkp->hps, NULL, 0);
  fs->src = openStream(fs->gkp->src, NULL, 0);
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
  resetStream(fs->seq, seqOffset, STREAM_UNTILEND);
  resetStream(fs->qlt, qltOffset, STREAM_UNTILEND);
  resetStream(fs->hps, hpsOffset, STREAM_UNTILEND);
  resetStream(fs->src, srcOffset, STREAM_UNTILEND);
}


void             closeFragStream(FragStream *fs) {
  closeStream(fs->frg);
  closeStream(fs->qlt);
  closeStream(fs->hps);
  closeStream(fs->src);
}


int64            getStartIndexFragStream(FragStream *fs) {
  return(getStartIndexStream(fs->frg));
}


int              nextFragStream(FragStream *fs, ReadStruct *rs, int streamFlags) {
  if (nextStream(fs->frg, &rs->gkfr) == 0)
    return(0);
  getFragData(fs->gkp, rs, streamFlags);
  return(1);
}
