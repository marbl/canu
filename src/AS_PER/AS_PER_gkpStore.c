
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

static char CM_ID[] = "$Id: AS_PER_gkpStore.c,v 1.44 2007-11-08 12:38:15 brianwalenz Exp $";

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

    sprintf(name,"%s/seq", path);
    fileCount += fileExists(name, 0, writable);

    sprintf(name,"%s/qlt", path);
    fileCount += fileExists(name, 0, writable);

    sprintf(name,"%s/hps", path);
    fileCount += fileExists(name, 0, writable);

    sprintf(name,"%s/src", path);
    fileCount += fileExists(name, 0, writable);

    sprintf(name,"%s/u2i", path);
    fileCount += fileExists(name, 0, writable);

    sprintf(name,"%s/uid", path);
    fileCount += fileExists(name, 0, writable);
  }

  return(fileCount == 10);
}


GateKeeperStore *
openGateKeeperStore(const char *path,
                    int   writable) {

  char              name[FILENAME_MAX];
  FILE             *gkpinfo;

  GateKeeperStore  *gkpStore = (GateKeeperStore *)safe_calloc(1, sizeof(GateKeeperStore));

  gkpStore->writable = 0;

  strcpy(gkpStore->storePath, path);

  gkpStore->bat = NULL;
  gkpStore->frg = NULL;
  gkpStore->lib = NULL;

  gkpStore->seq = NULL;
  gkpStore->qlt = NULL;
  gkpStore->hps = NULL;
  gkpStore->src = NULL;

  gkpStore->uid = NULL;

  gkpStore->UIDtoIID  = NULL;
  gkpStore->STRtoUID  = NULL;

  gkpStore->lib_cache = NULL;
  gkpStore->frgUID    = NULL;

  gkpStore->partnum = -1;
  gkpStore->partfrg = NULL;
  gkpStore->partqlt = NULL;
  gkpStore->parthps = NULL;
  gkpStore->partsrc = NULL;
  gkpStore->partmap = NULL;


  sprintf(name,"%s/gkp", gkpStore->storePath);
  errno = 0;
  gkpinfo = fopen(name, "r");
  if (errno) {
    fprintf(stderr, "failed to open gatekeeper store '%s': %s\n", name, strerror(errno));
    exit(1);
  }

  if (1 != AS_UTL_safeRead(gkpinfo, &gkpStore->gkp, "openGateKeeperStore:header", sizeof(GateKeeperStoreInfo), 1)) {
    fprintf(stderr, "failed to open gatekeeper store '%s': couldn't read the header (%s)\n", name, strerror(errno));
    exit(1);
  }
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

  //  writable is -1 if we are called from createGateKeeperPartition()

  if (writable >= 0) {
    char  mode[4];

    if (writable)
      strcpy(mode, "r+");
    else
      strcpy(mode, "r");

    sprintf(name,"%s/bat", gkpStore->storePath);
    gkpStore->bat   = openStore(name, mode);

    sprintf(name,"%s/frg", gkpStore->storePath);
    gkpStore->frg   = openStore(name, mode);

    sprintf(name,"%s/lib", gkpStore->storePath);
    gkpStore->lib   = openStore(name, mode);

    sprintf(name,"%s/seq", gkpStore->storePath);
    gkpStore->seq = openStore(name, mode);

    sprintf(name,"%s/qlt", gkpStore->storePath);
    gkpStore->qlt = openStore(name, mode);

    sprintf(name,"%s/hps", gkpStore->storePath);
    gkpStore->hps = openStore(name, mode);

    sprintf(name,"%s/src", gkpStore->storePath);
    gkpStore->src = openStore(name, mode);

    sprintf(name,"%s/uid", gkpStore->storePath);
    gkpStore->uid = openStore(name, mode);

    if ((NULL == gkpStore->bat) ||
        (NULL == gkpStore->frg) ||
        (NULL == gkpStore->lib) ||
        (NULL == gkpStore->seq) ||
        (NULL == gkpStore->qlt) ||
        (NULL == gkpStore->hps) ||
        (NULL == gkpStore->src) ||
        (NULL == gkpStore->uid)) {
      fprintf(stderr,"**** Failure to open Gatekeeper Store ...\n");
      assert(0);
    }
  }

  AS_UID_setGatekeeper(gkpStore);
  return(gkpStore);
}


GateKeeperStore *
createGateKeeperStore(const char *path) {
  char   name[FILENAME_MAX];
  FILE  *gkpinfo;

  GateKeeperStore  *gkpStore = (GateKeeperStore *)safe_calloc(1, sizeof(GateKeeperStore));

  strcpy(gkpStore->storePath, path);

  gkpStore->writable = 1;

  gkpStore->bat = NULL;
  gkpStore->frg = NULL;
  gkpStore->lib = NULL;

  gkpStore->seq = NULL;
  gkpStore->qlt = NULL;
  gkpStore->hps = NULL;
  gkpStore->src = NULL;
  gkpStore->uid = NULL;

  gkpStore->UIDtoIID  = NULL;
  gkpStore->STRtoUID  = NULL;

  gkpStore->lib_cache = NULL;
  gkpStore->frgUID    = NULL;

  gkpStore->partnum = -1;
  gkpStore->partfrg = NULL;
  gkpStore->partqlt = NULL;
  gkpStore->parthps = NULL;
  gkpStore->partsrc = NULL;
  gkpStore->partmap = NULL;

  AS_UTL_mkdir(path);

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

  AS_UTL_safeWrite(gkpinfo, &gkpStore->gkp, "createGateKeeperStore:header", sizeof(GateKeeperStoreInfo), 1);
  if (fclose(gkpinfo)) {
    fprintf(stderr, "failed to create gatekeeper store '%s': %s\n", name, strerror(errno));
    exit(1);
  }

  sprintf(name,"%s/bat", path);
  gkpStore->bat = createIndexStore(name, "bat", sizeof(GateKeeperBatchRecord), 1);

  sprintf(name,"%s/frg", path);
  gkpStore->frg = createIndexStore(name, "frg", sizeof(GateKeeperFragmentRecord), 1);

  sprintf(name,"%s/lib", path);
  gkpStore->lib = createIndexStore(name, "lib", sizeof(GateKeeperLibraryRecord), 1);

  sprintf(name,"%s/seq", path);
  gkpStore->seq = createStringStore(name, "seq");

  sprintf(name,"%s/qlt", path);
  gkpStore->qlt = createStringStore(name, "qlt");

  sprintf(name,"%s/hps", path);
  gkpStore->hps = createStringStore(name, "hps");

  sprintf(name,"%s/src", path);
  gkpStore->src = createStringStore(name, "src");

  sprintf(name,"%s/uid", path);
  gkpStore->uid = createStringStore(name, "uid");

  sprintf(name,"%s/u2i", path);
  gkpStore->UIDtoIID = CreateScalarHashTable_AS(32 * 1024);
  SaveHashTable_AS(name, gkpStore->UIDtoIID);

  AS_UID_setGatekeeper(gkpStore);
  return(gkpStore);
}


void
closeGateKeeperStore(GateKeeperStore *gkpStore) {
  char  name[FILENAME_MAX];
  FILE *gkpinfo;

  if (gkpStore == NULL)
    return;

  if (gkpStore->writable) {
    sprintf(name,"%s/gkp", gkpStore->storePath);
    errno = 0;
    gkpinfo = fopen(name, "w");
    if (errno) {
      fprintf(stderr, "failed to write gatekeeper store into to '%s': %s\n", name, strerror(errno));
      exit(1);
    }
    AS_UTL_safeWrite(gkpinfo, &gkpStore->gkp, "closeGateKeeperStore:header", sizeof(GateKeeperStoreInfo), 1);
    if (fclose(gkpinfo)) {
      fprintf(stderr, "failed to close gatekeeper store '%s': %s\n", name, strerror(errno));
      exit(1);
    }
  }

  if(gkpStore->bat != NULL)
    closeStore(gkpStore->bat);

  if(gkpStore->frg != NULL)
    closeStore(gkpStore->frg);

  if(gkpStore->lib != NULL)
    closeStore(gkpStore->lib);

  if(gkpStore->seq != NULL)
    closeStore(gkpStore->seq);

  if(gkpStore->qlt != NULL)
    closeStore(gkpStore->qlt);

  if(gkpStore->hps != NULL)
    closeStore(gkpStore->hps);

  if(gkpStore->src != NULL)
    closeStore(gkpStore->src);

  if(gkpStore->uid != NULL)
    closeStore(gkpStore->uid);

  if(gkpStore->UIDtoIID != NULL)
    DeleteHashTable_AS(gkpStore->UIDtoIID);

  if(gkpStore->STRtoUID != NULL)
    DeleteHashTable_AS(gkpStore->STRtoUID);

  safe_free(gkpStore->lib_cache);
  safe_free(gkpStore->frgUID);

  if(gkpStore->partfrg != NULL)
    closeStore(gkpStore->partfrg);

  if(gkpStore->partqlt != NULL)
    closeStore(gkpStore->partqlt);

  if(gkpStore->parthps != NULL)
    closeStore(gkpStore->parthps);

  if(gkpStore->partsrc != NULL)
    closeStore(gkpStore->partsrc);

  if (gkpStore->partmap != NULL)
    DeleteHashTable_AS(gkpStore->partmap);

  safe_free(gkpStore);
}



////////////////////////////////////////////////////////////////////////////////




GateKeeperStore *createGateKeeperPartition(const char *path, uint32 partnum) {
  char       name[FILENAME_MAX];

  GateKeeperStore *gkp = openGateKeeperStore(path, -1);

  gkp->partnum = partnum;

  sprintf(name,"%s/frg.%03d", gkp->storePath, partnum);
  gkp->partfrg = createIndexStore(name, "partfrg", sizeof(GateKeeperFragmentRecord), 1);

  sprintf(name,"%s/qlt.%03d", gkp->storePath, partnum);
  gkp->partqlt = createStringStore(name, "partqlt");

  sprintf(name,"%s/hps.%03d", gkp->storePath, partnum);
  gkp->parthps = createStringStore(name, "parthps");

  sprintf(name,"%s/src.%03d", gkp->storePath, partnum);
  gkp->partsrc = createStringStore(name, "partsrc");

  return(gkp);
}




void       loadGateKeeperPartition(GateKeeperStore *gkp, uint32 partnum) {
  char       name[FILENAME_MAX];
  int        i;

  if (gkp->partnum != -1) {
    fprintf(stderr, "WARNING:  Throwing out partition %d to load %d\n",
            gkp->partnum, partnum);

    if(gkp->partfrg != NULL)
      closeStore(gkp->partfrg);

    if(gkp->partqlt != NULL)
      closeStore(gkp->partqlt);

    if(gkp->parthps != NULL)
      closeStore(gkp->parthps);

    if(gkp->partsrc != NULL)
      closeStore(gkp->partsrc);

    if (gkp->partmap != NULL)
      DeleteHashTable_AS(gkp->partmap);
  }

  sprintf(name,"%s/frg.%03d", gkp->storePath, partnum);
  if (fileExists(name, 0, FALSE) == 0) {
    fprintf(stderr, "loadGateKeeperPartition()--  Partition %d doesn't exist; normal store used instead.\n", partnum);
    return;
  }

  gkp->partnum = partnum;

  //  load all our data

  sprintf(name,"%s/frg.%03d", gkp->storePath, partnum);
  gkp->partfrg = loadStorePartial(name, 0, 0);

  sprintf(name,"%s/qlt.%03d", gkp->storePath, partnum);
  gkp->partqlt = loadStorePartial(name, 0, 0);

  sprintf(name,"%s/hps.%03d", gkp->storePath, partnum);
  gkp->parthps = loadStorePartial(name, 0, 0);

  sprintf(name,"%s/src.%03d", gkp->storePath, partnum);
  gkp->partsrc = loadStorePartial(name, 0, 0);

  //  zip through the frg and build a map from iid to the frg record

  int64  firstElem = getFirstElemStore(gkp->partfrg);
  int64  lastElem  = getLastElemStore(gkp->partfrg);

  gkp->partmap = CreateScalarHashTable_AS(lastElem + 1);

  for(i = firstElem; i <= lastElem; i++) {
    GateKeeperFragmentRecord *p = getIndexStorePtr(gkp->partfrg, i);

    if (InsertInHashTable_AS(gkp->partmap,
                             (uint64)p->readIID, 0,
                             (INTPTR)(p), 0) != HASH_SUCCESS)
      assert(0);
  }
}


////////////////////////////////////////////////////////////////////////////////


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
//
//  UID lookups.
//
//  For numeric UIDs (those with both isString and UID valid), you can
//  do a UID to IID mapping.
//
//  String UIDs (those with only UIDstring defined) must be added to
//  the store first.  The UIDstring itself is stored in a String
//  store, just like sequence and source.  The index into this store
//  is saved in the (isString, UID) pair.
//
//  A hash must be maintained that maps the UIDstring to the location
//  in the store (which is also exactly the (isString, UID) pair.
//  This will be used to lookup a string UID.
//
//  The gkp files are:
//    'u2i' -> the UID to IID mapping
//    'uid' -> the string store itself


static
void
loadGatekeeperSTRtoUID(GateKeeperStore *gkp) {
  if (gkp->STRtoUID == NULL) {
    gkp->uid = convertStoreToMemoryStore(gkp->uid);

    gkp->STRtoUID = CreateStringHashTable_AS(32 * 1024);
    //SaveHashTable_AS(name, gkpStore->STRtoUID);

    char          *uidptr = NULL;
    int64          uidoff = 0;
    uint32         actlen = 0;

    while ((uidptr = getStringStorePtr(gkp->uid, uidoff, &actlen)) != NULL) {
      if (InsertInHashTable_AS(gkp->STRtoUID,
                               (INTPTR)uidptr, actlen,
                               uidoff, 0) == HASH_FAILURE) {
        fprintf(stderr, "loadGatekeeperSTRtoUID()-- failed to insert uid '%s' into store; already there?!\n", uidptr);
        assert(0);
      }

      uidoff += actlen + 1;
    }
  }
  assert(gkp->STRtoUID != NULL);
}


//  Takes a uid, with UIDstring set, populates the rest of the uid
//  (isString, UID).
//
AS_UID
AS_GKP_getUIDfromString(GateKeeperStore *gkp, AS_UID uid) {
  uint64  loc = 0;

  loadGatekeeperSTRtoUID(gkp);

  if (LookupInHashTable_AS(gkp->STRtoUID,
                           (INTPTR)uid.UIDstring, strlen(uid.UIDstring),
                           &loc, 0)) {
    uid.isString  = 1;
    uid.UID       = loc;
  } else {
    uid = AS_UID_undefined();
  }
  return(uid);
}

//  Takes a uid with (isString, UID) set (and without UIDstring set),
//  returns a uid with all three set.
//
AS_UID
AS_GKP_getUID(GateKeeperStore *gkp, AS_UID uid) {

  uid.UIDstring = NULL;

  if (uid.isString) {
    uint32  actlen = 0;
    loadGatekeeperSTRtoUID(gkp);
    uid.UIDstring = getStringStorePtr(gkp->uid, uid.UID, &actlen);
  }

  return(uid);
}


//  Adds a uid (UIDstring) to the store, populates the rest of the UID.
//
AS_UID
AS_GKP_addUID(GateKeeperStore *gkp, AS_UID uid) {

  //  Could probably just return...might want to verify that this
  //  UIDstring is the same as the UID.

  assert((uid.isString == 0) && (uid.UID == 0) && (uid.UIDstring != NULL));

  uint64     loc    = 0;
  uint64     len    = strlen(uid.UIDstring);

  loadGatekeeperSTRtoUID(gkp);

  //  If the UID is already in the store, just return as if it was
  //  new.  Otherwise, add it to the store.

  if (LookupInHashTable_AS(gkp->STRtoUID, (INTPTR)uid.UIDstring, len, &loc, 0) == FALSE) {
    char    *str = NULL;
    uint32   act = 0;

    loc = getLastElemStore(gkp->uid);

    //  Stash the UID on disk.
    appendStringStore(gkp->uid, uid.UIDstring, len);

    str = getStringStorePtr(gkp->uid, loc, &act);

    if (InsertInHashTable_AS(gkp->STRtoUID,
                             (INTPTR)str, len,
                             loc, 0) == HASH_FAILURE) {
      fprintf(stderr, "setGatekeeperUID()-- failed to insert uid '%s' into store; already there?!\n", uid.UIDstring);
      assert(0);
    }
  }

  uid.isString  = 1;
  uid.UID       = loc;
  uid.UIDstring = NULL;  //  its not valid until we get().

  return(uid);
}




static
void
loadGatekeeperUIDtoIID(GateKeeperStore *gkp) {
  if (gkp->UIDtoIID == NULL) {
    char  name[FILENAME_MAX];
    sprintf(name,"%s/u2i", gkp->storePath);
    gkp->UIDtoIID = LoadUIDtoIIDHashTable_AS(name);
  }
  assert(gkp->UIDtoIID != NULL);
}

//  The only public accessor for the persistent hash in the
//  gatekeeper.  Returns the IID, or 0 if the uid was not found.
//
AS_IID
getGatekeeperUIDtoIID(GateKeeperStore *gkp, AS_UID uid, uint32 *type) {
  uint64   iid = 0;
  loadGatekeeperUIDtoIID(gkp);
  if (AS_UID_isDefined(uid))
    LookupInHashTable_AS(gkp->UIDtoIID, AS_UID_toInteger(uid), 0, &iid, type);
  return((AS_IID)iid);
}

int
setGatekeeperUIDtoIID(GateKeeperStore *gkp, AS_UID uid, AS_IID iid, uint32 type) {
  loadGatekeeperUIDtoIID(gkp);
  assert(AS_UID_isDefined(uid) == TRUE);
  assert(AS_IID_isDefined(iid) == TRUE);
  return(InsertInHashTable_AS(gkp->UIDtoIID, AS_UID_toInteger(uid), 0, (uint64)iid, type));
}






static
void
loadGatekeeperIIDtoUID(GateKeeperStore *gkp) {

  if (gkp->frgUID)
    return;

  uint64   lastIID = getNumGateKeeperFragments(gkp);

  gkp->frgUID = (uint64 *)safe_calloc(lastIID, sizeof(uint64));

  HashTable_Iterator_AS   iterator  = {0};
  uint64                  key       = 0;
  uint64                  value     = 0;
  uint32                  valuetype = 0;
  uint32                  added     = 0;

  loadGatekeeperUIDtoIID(gkp);
  InitializeHashTable_Iterator_AS(gkp->UIDtoIID, &iterator);

  while (NextHashTable_Iterator_AS(&iterator, &key, &value, &valuetype)) {
    if (valuetype == AS_IID_FRG)
      gkp->frgUID[value] = key;
  }
}

AS_UID
getGatekeeperIIDtoUID(GateKeeperStore *gkp, AS_IID iid, uint32 type) {
  AS_UID  uid = AS_UID_undefined();

  loadGatekeeperIIDtoUID(gkp);

  switch (type) {
    case AS_IID_FRG:
      uid = AS_UID_fromInteger(gkp->frgUID[iid]);
      break;
    case AS_IID_LIB:
      uid = getGateKeeperLibrary(gkp, iid)->libraryUID;
      break;
    case AS_IID_BAT:
      assert(0);
      break;
    default:
      break;
  }

  return(uid);
}




////////////////////////////////////////////////////////////////////////////////


static
int
AS_PER_decodeLibraryFeaturesBoolean(char *feature, char *value) {
  int  ret = 0;

  //  Decodes a string with 0/1, false/true, no/yes into an integer flag.

  switch (value[0]) {
    case '0':
    case 'f':
    case 'F':
    case 'n':
    case 'N':
      ret = 0;
      break;
    case '1':
    case 't':
    case 'T':
    case 'y':
    case 'Y':
      ret = 1;
      break;
    default:
      fprintf(stderr, "AS_PER_decodeLibraryFeatures()-- Found feature '%s' but has unknown boolean value '%s'\n",
              feature, value);
      break;
  }

  return(ret);
}



void
AS_PER_decodeLibraryFeatures(GateKeeperLibraryRecord *gkpl,
                             LibraryMesg             *lmesg) {
  int f;
  for (f=0; f<lmesg->num_features; f++) {
    char *fea = lmesg->features[f];
    char *val = lmesg->values[f];


    //  isNotRandom --
    if        (strcasecmp(fea, "isNotRandom") == 0) {
      gkpl->isNotRandom = AS_PER_decodeLibraryFeaturesBoolean("isNotRandom", val);
    }

    //  doNotOverlapTrim -- 
    else if (strcasecmp(fea, "doNotOverlapTrim") == 0) {
      gkpl->doNotOverlapTrim = AS_PER_decodeLibraryFeaturesBoolean("doNotOverlapTrim", val);
    }

    //  doNotTrustHomopolymerRuns -- 
    else if (strcasecmp(fea, "doNotTrustHomopolymerRuns") == 0) {
      gkpl->doNotTrustHomopolymerRuns = AS_PER_decodeLibraryFeaturesBoolean("doNotTrustHomopolymerRuns", val);
    }

    //  hpsIsPeakSpacing -- 
    else if (strcasecmp(fea, "hpsIsPeakSpacing") == 0) {
      gkpl->hpsIsPeakSpacing = AS_PER_decodeLibraryFeaturesBoolean("hpsIsPeakSpacing", val);
    }

    //  hpsIsFlowGram -- 
    else if (strcasecmp(fea, "hpsIsFlowGram") == 0) {
      gkpl->hpsIsFlowGram = AS_PER_decodeLibraryFeaturesBoolean("hpsIsFlowGram", val);
    }

    else {
      fprintf(stderr, "AS_PER_decodeLibraryFeatures()-- Found feature '%s' but don't understand it.\n",
              fea);
    }
  }
}


void
AS_PER_encodeLibraryFeaturesCleanup(LibraryMesg *lmesg) {
  while (lmesg->num_features > 0) {
    lmesg->num_features--;
    safe_free(lmesg->features[lmesg->num_features]);
    safe_free(lmesg->values  [lmesg->num_features]);
  }
  safe_free(lmesg->features);
  safe_free(lmesg->values);
}


void
AS_PER_encodeLibraryFeatures(GateKeeperLibraryRecord *gkpl,
                             LibraryMesg             *lmesg) {

  //  Examine the gkpl, allocate space to encode the features into
  //  features/values, return the number of features encoded.
  //
  //  Be sure to call AS_PER_encodeLibraryFeaturesCleanup to properly
  //  cleanup the LibraryMesg after it is written!

  //  We can hardcode the maximum number of features we expect to be
  //  writing.  Otherwise, we should count the number of features we
  //  want to encode, allocate....but what a pain.
  //
  lmesg->num_features = 0;
  lmesg->features     = (char **)safe_malloc(5 * sizeof(char*));
  lmesg->values       = (char **)safe_malloc(5 * sizeof(char*));

  int    nf  = 0;
  char **fea = lmesg->features;
  char **val = lmesg->values;

  //  Mostly for debugging, but just might be generally a
  //  GoodThing(tm) to always specify optional features.
  int    alwaysEncode = 1;

  if (gkpl->isNotRandom || alwaysEncode) {
    fea[nf] = (char *)safe_malloc(32 * sizeof(char));
    val[nf] = (char *)safe_malloc(32 * sizeof(char));
    sprintf(fea[nf], "isNotRandom");
    sprintf(val[nf], "%d", gkpl->isNotRandom);
    nf++;
  }

  if (gkpl->doNotOverlapTrim || alwaysEncode) {
    fea[nf] = (char *)safe_malloc(32 * sizeof(char));
    val[nf] = (char *)safe_malloc(32 * sizeof(char));
    sprintf(fea[nf], "doNotOverlapTrim");
    sprintf(val[nf], "%d", gkpl->doNotOverlapTrim);
    nf++;
  }

  if (gkpl->doNotTrustHomopolymerRuns || alwaysEncode) {
    fea[nf] = (char *)safe_malloc(32 * sizeof(char));
    val[nf] = (char *)safe_malloc(32 * sizeof(char));
    sprintf(fea[nf], "doNotTrustHomopolymerRuns");
    sprintf(val[nf], "%d", gkpl->doNotTrustHomopolymerRuns);
    nf++;
  }

  if (gkpl->hpsIsPeakSpacing || alwaysEncode) {
    fea[nf] = (char *)safe_malloc(32 * sizeof(char));
    val[nf] = (char *)safe_malloc(32 * sizeof(char));
    sprintf(fea[nf], "hpsIsPeakSpacing");
    sprintf(val[nf], "%d", gkpl->hpsIsPeakSpacing);
    nf++;
  }

  if (gkpl->hpsIsFlowGram || alwaysEncode) {
    fea[nf] = (char *)safe_malloc(32 * sizeof(char));
    val[nf] = (char *)safe_malloc(32 * sizeof(char));
    sprintf(fea[nf], "hpsIsFlowGram");
    sprintf(val[nf], "%d", gkpl->hpsIsFlowGram);
    nf++;
  }

  lmesg->num_features = nf;
}


////////////////////////////////////////////////////////////////////////////////

//  Set clear region 'which' and all later clear regions.  You are
//  explicitly not allowed to set the original clear range.
//
void        setFragRecordClearRegion(fragRecord *fr,
                                     uint32 start,
                                     uint32 end,
                                     uint32 which) {
  assert(which >  AS_READ_CLEAR_VEC);
  assert(which <= AS_READ_CLEAR_LATEST);
  for (; which <= AS_READ_CLEAR_LATEST; which++) {
    fr->gkfr.clearBeg[which] = start;
    fr->gkfr.clearEnd[which] = end;
  }
}


void        getFragRecordClearRegion(fragRecord *fr, uint32 *start, uint32 *end, uint32 which) {
  if (which == AS_READ_CLEAR_UNTRIM) {
    *start = 0;
    *end   = fr->gkfr.seqLen;
  } else {
    assert(which <= AS_READ_CLEAR_LATEST);
    *start = fr->gkfr.clearBeg[which];
    *end   = fr->gkfr.clearEnd[which];
  }
}


uint32      getFragRecordClearRegionBegin(fragRecord *fr, uint32 which) {
  if (which == AS_READ_CLEAR_UNTRIM) {
    return(0);
  } else {
    assert(which <= AS_READ_CLEAR_LATEST);
    return(fr->gkfr.clearBeg[which]);
  }
}


uint32      getFragRecordClearRegionEnd  (fragRecord *fr, uint32 which) {
  if (which == AS_READ_CLEAR_UNTRIM) {
    return(fr->gkfr.seqLen);
  } else {
    assert(which <= AS_READ_CLEAR_LATEST);
    return(fr->gkfr.clearEnd[which]);
  }
}


////////////////////////////////////////////////////////////////////////////////


static
void
getFragData(GateKeeperStore *gkp, fragRecord *fr, int streamFlags) {
  uint32            actualLength = 0;
  StoreStruct      *store        = NULL;

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
      assert(gkp->partmap == NULL);
      getStringStore(gkp->seq,
                       fr->gkfr.seqOffset,
                       fr->seq,
                       MAX_SEQ_LENGTH,
                       &actualLength);
      fr->seq[actualLength] = 0;
    }
  }
  if (streamFlags & FRAG_S_QLT) {
    assert(fr->hasSEQ == 0);
    fr->hasSEQ = 1;
    fr->hasQLT = 1;
    if (fr->gkfr.seqLen > 0) {
      store = gkp->qlt;
      if (gkp->partmap)
        store = gkp->partqlt;
      getStringStore(store,
                       fr->gkfr.qltOffset,
                       fr->qlt,
                       MAX_SEQ_LENGTH,
                       &actualLength);
      fr->qlt[actualLength] = 0;
      decodeSequenceQuality(fr->qlt, fr->seq, fr->qlt);
    }
  }
  if (streamFlags & FRAG_S_HPS) {
    fr->hasHPS = 1;
    if (fr->gkfr.hpsLen > 0) {
      store = gkp->hps;
      if (gkp->partmap)
        store = gkp->parthps;
      getStringStore(store,
                       fr->gkfr.hpsOffset,
                       fr->hps,
                       MAX_HPS_LENGTH,
                       &actualLength);
      fr->hps[actualLength] = 0;
    }
  }
  if (streamFlags & FRAG_S_SRC) {
    fr->hasSRC = 1;
    if (fr->gkfr.srcLen > 0) {
      store = gkp->src;
      if (gkp->partmap)
        store = gkp->partsrc;
      getStringStore(store,
                       fr->gkfr.srcOffset,
                       fr->src,
                       MAX_SRC_LENGTH,
                       &actualLength);
      fr->src[actualLength] = 0;
    }
  }
}


void    getFrag(GateKeeperStore *gkp, AS_IID    iid, fragRecord *fr, int32 flags) {
  if (gkp->partmap == NULL) {
    getIndexStore(gkp->frg, iid, &fr->gkfr);
  } else {
    GateKeeperFragmentRecord *gkfr;

    gkfr = (GateKeeperFragmentRecord *)(INTPTR)LookupValueInHashTable_AS(gkp->partmap, iid, 0);
    if (gkfr == NULL) {
      fprintf(stderr, "getFrag()-- ERROR!  IID "F_IID" not in partition!\n", iid);
      assert(0);
    }

    memcpy(&fr->gkfr, gkfr, sizeof(GateKeeperFragmentRecord));
  }

  getFragData(gkp, fr, flags);
}


void    setFrag(GateKeeperStore *gkp, AS_IID    iid, fragRecord *fr) {
  assert(gkp->partmap == NULL);
  setIndexStore(gkp->frg, iid, &fr->gkfr);
}

void    delFrag(GateKeeperStore *gkp, AS_IID    iid) {
  GateKeeperFragmentRecord   gkfr;
  AS_IID                     miid;

  assert(gkp->partmap == NULL);

  //  Delete fragment with iid from the store.  If the fragment has a
  //  mate, remove the mate relationship from both fragmentss.

  getIndexStore(gkp->frg, iid, &gkfr);
  miid = gkfr.mateIID;
  gkfr.deleted = 1;
  gkfr.mateIID = 0;
  setIndexStore(gkp->frg, iid, &gkfr);

  if (miid > 0) {
    getIndexStore(gkp->frg, miid, &gkfr);
    gkfr.mateIID = 0;
    setIndexStore(gkp->frg, miid, &gkfr);
  }
}



////////////////////////////////////////////////////////////////////////////////


FragStream      *openFragStream(GateKeeperStore *gkp, int flags) {
  FragStream  *fs = (FragStream *)safe_malloc(sizeof(FragStream));
  fs->gkp   = gkp;
  fs->frg   = NULL;
  fs->seq   = NULL;
  fs->qlt   = NULL;
  fs->hps   = NULL;
  fs->src   = NULL;
  fs->flags = flags;

  fs->frg   = openStream(fs->gkp->frg);

  if ((fs->flags & FRAG_S_SEQ) && !(fs->flags & FRAG_S_QLT))
    fs->seq   = openStream(fs->gkp->seq);
  if (fs->flags & FRAG_S_QLT)
    fs->qlt   = openStream(fs->gkp->qlt);
  if (fs->flags & FRAG_S_HPS)
    fs->hps   = openStream(fs->gkp->hps);
  if (fs->flags & FRAG_S_SRC)
    fs->src   = openStream(fs->gkp->src);

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

    getIndexStore(fs->gkp->frg, startIndex, &gkpf);

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
}



int              nextFragStream(FragStream *fs, fragRecord *fr) {
  uint32    actualLength = 0;

  if (nextStream(fs->frg, &fr->gkfr, 0, NULL) == 0)
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
    nextStream(fs->seq, fr->seq, MAX_SEQ_LENGTH, &actualLength);
    fr->seq[actualLength] = 0;
  }

  if (fs->flags & FRAG_S_QLT) {
    assert(fr->hasSEQ == 0);
    fr->hasSEQ = 1;
    fr->hasQLT = 1;
    nextStream(fs->qlt, fr->qlt, MAX_SEQ_LENGTH, &actualLength);
    fr->qlt[actualLength] = 0;
    decodeSequenceQuality(fr->qlt, fr->seq, fr->qlt);
  }

  if (fs->flags & FRAG_S_HPS) {
    fr->hasHPS = 1;
    nextStream(fs->hps, fr->hps, MAX_HPS_LENGTH, &actualLength);
    fr->hps[actualLength] = 0;
  }

  if (fs->flags & FRAG_S_SRC) {
    fr->hasSRC = 1;
    nextStream(fs->src, fr->src, MAX_SRC_LENGTH, &actualLength);
    fr->src[actualLength] = 0;
  }

  return(1);
}



void
loadGateKeeperStorePartial(GateKeeperStore *gkp,
                           int64            firstElem,
                           int64            lastElem,
                           int              flags) {
  GateKeeperFragmentRecord   gkf;

  int64  seqFirst = 0, seqLast = 0;
  int64  qltFirst = 0, qltLast = 0;
  int64  hpsFirst = 0, hpsLast = 0;
  int64  srcFirst = 0, srcLast = 0;

  if (firstElem != 0) {
    getIndexStore(gkp->frg, firstElem, &gkf);

    seqFirst = gkf.seqOffset;
    qltFirst = gkf.qltOffset;
    hpsFirst = gkf.hpsOffset;
    srcFirst = gkf.srcOffset;
  }

  if ((lastElem != 0) &&
      (lastElem != getLastElemFragStore(gkp))) {
    getIndexStore(gkp->frg, lastElem + 1, &gkf);

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
}
