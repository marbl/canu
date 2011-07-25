
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

static char *rcsid = "$Id: AS_PER_gkStore.C,v 1.21 2011-07-25 20:00:47 mkotelbajcvi Exp $";

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


#define AS_GKP_CURRENT_VERSION    8


gkStore::gkStore() {

  gkStore_clear();

  uid      = createStringStore(NULL, "uid");
  STRtoUID = CreateStringHashTable_AS();
}



gkStore::gkStore(const char *path, int partition) {
  char       name[FILENAME_MAX];

  gkStore_clear();

  if ((path == NULL) || (path[0] == 0)) {
    fprintf(stderr, "gkStore::gkStore()--  ERROR!  Empty gkpStore path supplied.\n");
    exit(1);
  }

  strcpy(storePath, path);

  assert(partition > 0);
  partnum = partition;

  sprintf(name,"%s/fpk.%03d", storePath, partnum);
  partfpk = createIndexStore(name, "partfpk", sizeof(gkPackedFragment), 1);
  sprintf(name,"%s/qpk.%03d", storePath, partnum);
  partqpk = createIndexStore(name, "partqpk", sizeof(gkPackedSequence), 1);

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
gkStore::gkStore_open(int writable, int doNotUseUIDs) {
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
    fprintf(stderr, "failed to open gatekeeper store '%s': couldn't read the header (%s)\n",
            name, strerror(errno));
    exit(1);
  }

  if (!feof(gkpinfo)) {
    uint32 nr = inf.numPacked + inf.numNormal + inf.numStrobe + 1;
    uint32 na = 0;
    uint32 nb = 0;

    IIDtoTYPE = (uint8  *)safe_malloc(sizeof(uint8)  * nr);
    IIDtoTIID = (uint32 *)safe_malloc(sizeof(uint32) * nr);

    na = AS_UTL_safeRead(gkpinfo, IIDtoTYPE, "gkStore_open:header", sizeof(uint8), nr);
    nb = AS_UTL_safeRead(gkpinfo, IIDtoTIID, "gkStore_open:header", sizeof(uint32), nr);

    //  If EOF was hit, and nothing was read, there is no index saved.  Otherwise, something was
    //  read, and we fail if either was too short.

    if ((feof(gkpinfo)) && (na == 0) && (nb == 0)) {
      safe_free(IIDtoTYPE);
      safe_free(IIDtoTIID);
    } else if ((na != nr) || (nb != nr)) {
      fprintf(stderr, "failed to open gatekeeper store '%s': couldn't read the IID maps (%s)\n",
              name, strerror(errno));
      exit(1);
    }
  }

  fclose(gkpinfo);

  if (inf.gkMagic != 1) {
    fprintf(stderr, "gkStore_open()-- Invalid magic; corrupt %s/inf?\n", name);
    exit(1);
  }

  if (inf.gkVersion != AS_GKP_CURRENT_VERSION) {
    fprintf(stderr, "gkStore_open()-- Invalid version!  Found version %d, code supports version %d.\n",
            inf.gkVersion, AS_GKP_CURRENT_VERSION);
    exit(1);
  }

  if ((inf.gkLibrarySize        != sizeof(gkLibrary)) ||
      (inf.gkPackedSequenceSize != sizeof(gkPackedSequence)) ||
      (inf.gkPackedFragmentSize != sizeof(gkPackedFragment)) ||
      (inf.gkNormalFragmentSize != sizeof(gkNormalFragment)) ||
      (inf.gkStrobeFragmentSize != sizeof(gkStrobeFragment)) ||
      (inf.readMaxLenBits       != AS_READ_MAX_NORMAL_LEN_BITS)) {
    fprintf(stderr, "gkStore_open()--  ERROR!  Incorrect element sizes; code and store are incompatible.\n");
    fprintf(stderr, "  gkLibrary:                    store %5d   code %5d bytes\n", inf.gkLibrarySize, (int)sizeof(gkLibrary));
    fprintf(stderr, "  gkPackedSequence:             store %5d   code %5d bytes\n", inf.gkPackedSequenceSize, (int)sizeof(gkPackedSequence));
    fprintf(stderr, "  gkPackedFragment:             store %5d   code %5d bytes\n", inf.gkPackedFragmentSize, (int)sizeof(gkPackedFragment));
    fprintf(stderr, "  gkNormalFragment:             store %5d   code %5d bytes\n", inf.gkNormalFragmentSize, (int)sizeof(gkNormalFragment));
    fprintf(stderr, "  gkStrobeFragment:             store %5d   code %5d bytes\n", inf.gkStrobeFragmentSize, (int)sizeof(gkStrobeFragment));
    fprintf(stderr, "  AS_READ_MAX_NORMAL_LEN_BITS:  store %5d   code %5d\n", inf.readMaxLenBits, AS_READ_MAX_NORMAL_LEN_BITS);
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
  sprintf(name,"%s/qpk", storePath);
  qpk   = openStore(name, mode);

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

  //  UIDtoIID and STRtoUID are loaded on demand.
  UIDtoIID = NULL;
  STRtoUID = NULL;
  doNotLoadUIDs = doNotUseUIDs;

  if ((NULL == fpk) || (NULL == qpk) ||
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

  isReadOnly      = 0;
  isCreating      = 1;
  doNotLoadUIDs   = 0;

  sprintf(name,"%s/inf", storePath);
  if (AS_UTL_fileExists(name, FALSE, TRUE)) {
    fprintf(stderr, "GateKeeper Store '%s' exists; will not create a new one on top of it.\n", storePath);
    exit(1);
  }

  AS_UTL_mkdir(storePath);

  inf.gkMagic              = 1;
  inf.gkVersion            = AS_GKP_CURRENT_VERSION;
  inf.gkLibrarySize        = sizeof(gkLibrary);
  inf.gkPackedSequenceSize = sizeof(gkPackedSequence);
  inf.gkPackedFragmentSize = sizeof(gkPackedFragment);
  inf.gkNormalFragmentSize = sizeof(gkNormalFragment);
  inf.gkStrobeFragmentSize = sizeof(gkStrobeFragment);
  inf.gkPlacementSize      = sizeof(gkPlacement);
  inf.readMaxLenBits       = AS_READ_MAX_NORMAL_LEN_BITS;

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
  sprintf(name,"%s/qpk", storePath);
  qpk = createIndexStore(name, "qpk", sizeof(gkPackedSequence), 1);

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

  IIDmax    = 0;
  IIDtoTYPE = NULL;
  IIDtoTIID = NULL;
}


void
gkStore::gkStore_construct(const char * path, int creatable, int writable, int doNotUseUIDs) {
  char   name[FILENAME_MAX];

  gkStore_clear();

  if ((path == NULL) || (path[0] == 0)) {
    fprintf(stderr, "gkStore::gkStore()--  ERROR!  Empty gkpStore path supplied.\n");
    exit(1);
  }

  strcpy(storePath, path);

  sprintf(name,"%s/inf", storePath);

  if ((AS_UTL_fileExists(name, FALSE, writable)) && (creatable == 0)) {
    gkStore_open(writable, doNotUseUIDs);
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


gkStore::gkStore(const char *path, int creatable, int writable) {
   gkStore_construct(path, creatable, writable, FALSE);
}

gkStore::gkStore(const char *path, int creatable, int writable, int doNotUseUIDs) {
  gkStore_construct(path, creatable, writable, doNotUseUIDs);
}



gkStore::~gkStore() {
  char  name[FILENAME_MAX];
  FILE *gkpinfo;

  if (isCreating) {
    sprintf(name, "%s/inf", storePath);
    errno = 0;
    gkpinfo = fopen(name, "w");
    if (errno) {
      fprintf(stderr, "failed to write gatekeeper store into to '%s': %s\n", name, strerror(errno));
      exit(1);
    }

    AS_UTL_safeWrite(gkpinfo, &inf, "closegkStore:header", sizeof(gkStoreInfo), 1);

    if (IIDtoTYPE) {
      AS_UTL_safeWrite(gkpinfo, IIDtoTYPE, "closegkStore::IIDtoTYPE", sizeof(uint8),  inf.numPacked + inf.numNormal + inf.numStrobe + 1);
      AS_UTL_safeWrite(gkpinfo, IIDtoTIID, "closegkStore::IIDtoTIID", sizeof(uint32), inf.numPacked + inf.numNormal + inf.numStrobe + 1);
    }

    if (fclose(gkpinfo)) {
      fprintf(stderr, "failed to close gatekeeper store '%s': %s\n", name, strerror(errno));
      exit(1);
    }
  }

  closeStore(fpk);
  closeStore(qpk);

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
  closeStore(partqpk);
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
  isCreating = 0;

  memset(&inf, 0, sizeof(gkStoreInfo));

  fpk = NULL;
  qpk = NULL;

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
  IIDtoTYPE = NULL;
  IIDtoTIID = NULL;

  clearRange = NULL;

  partnum = 0;

  partfpk = partqpk = NULL;
  partfnm = partqnm = NULL;
  partfsb = partqsb = NULL;

  partmap = NULL;
  
  doNotLoadUIDs = FALSE;
}


void
gkStore::gkStore_delete(void) {
  char   name[FILENAME_MAX];

  //  This function does both a ~gkStore (needed to close open files, and etc) and then removes the
  //  files from disk.  It does not handle a partitioned store.

  closeStore(fpk);
  closeStore(qpk);

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
  closeStore(partqpk);
  closeStore(partfnm);
  closeStore(partfsb);
  closeStore(partqnm);
  closeStore(partqsb);

  DeleteHashTable_AS(partmap);

  //  Remove files (and close/purge clear ranges).

  sprintf(name,"%s/inf", storePath);  unlink(name);

  sprintf(name,"%s/fpk", storePath);  unlink(name);
  sprintf(name,"%s/qpk", storePath);  unlink(name);

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



uint64
gkStore::gkStore_metadataSize(void) {
  uint64   totalSize = 0;
  char     name[FILENAME_MAX];

  sprintf(name, "%s/fpk", gkStore_path());
  totalSize += AS_UTL_sizeOfFile(name);

  sprintf(name, "%s/fmd", gkStore_path());
  totalSize += AS_UTL_sizeOfFile(name);

  sprintf(name, "%s/fsb", gkStore_path());
  totalSize += AS_UTL_sizeOfFile(name);

  return(totalSize);
}


void
gkStore::gkStore_metadataCaching(bool enable) {

  assert(partnum == 0);

  if (enable == true) {
    fpk = convertStoreToMemoryStore(fpk);
    fnm = convertStoreToMemoryStore(fnm);
    fsb = convertStoreToMemoryStore(fsb);
  }
}
