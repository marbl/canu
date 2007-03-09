
/**************************************************************************
 * This file is part of Celera Assembler, a software program that 
 * assembles whole-genome shotgun reads into contigs and scaffolds.
 * Copyright (C) 2007, J. Craig Venter Institute. All rights reserved.
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

static char CM_ID[] = "$Id: overlapStore_merge.c,v 1.1 2007-03-09 07:29:23 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <fcntl.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "AS_global.h"
#include "AS_UTL_fileIO.h"
#include "AS_OVS_overlap.h"
#include "AS_OVS_overlapFile.h"
#include "AS_OVS_overlapStore.h"

#include "overlapStore.h"


void
renameToBackup(char const *storeName, char const *name) {
  char   orig[FILENAME_MAX];
  char   bkup[FILENAME_MAX];

  sprintf(orig, "%s/%s", storeName, name);
  sprintf(bkup, "%s/%s~", storeName, name);

  errno = 0;
  rename(orig, bkup);
  if (errno) {
    fprintf(stderr, "overlapStore: ERROR: failed to make backup of '%s' into '%s': %s\n", orig, bkup, strerror(errno));
    assert(0);
  }
}


void
nukeBackup(char const *storeName, char const *name) {
  char   bkup[FILENAME_MAX];

  sprintf(bkup, "%s/%s~", storeName, name);

  errno = 0;
  //unlink(bkup);
  if ((errno) && (errno != ENOENT))
    fprintf(stderr, "overlapStore: WARNING: failed to remove backup '%s': %s\n", bkup, strerror(errno));
}


void
mergeStore(char *storeName, char *mergeName) {
  int             i;
  char            name[FILENAME_MAX];

  OverlapStore   *merge = NULL;
  OverlapStore   *store = NULL;

  uint32          numberOfOrigStoreFiles = 0;

  OVSoverlap      movl;
  OVSoverlap      sovl;

  int             mValid = 1;
  int             sValid = 1;

  int             saveSpace = 1;

  //  Open the two stores.  We open the original only to get the number
  //  of files in it.  It'll be reopened later.
  //
  merge = AS_OVS_openOverlapStore(mergeName);
  store = AS_OVS_openOverlapStore(storeName);

  numberOfOrigStoreFiles = store->ovs.highestFileIndex + 1;

  AS_OVS_closeOverlapStore(store);


  //  Rename the original store to backup files.
  //
  renameToBackup(storeName, "ovs");
  renameToBackup(storeName, "idx");
  for (i=1; i<numberOfOrigStoreFiles; i++) {
    sprintf(name, "%04d", i);
    renameToBackup(storeName, name);
  }


  //  Recreate a store in the same place as the original store.
  //
  store = AS_OVS_createOverlapStore(storeName, FALSE);


  //  Now just add stuff to the new store.
  //
  //  Yes, it would have been trivial to merge stores by creating a
  //  completely new store, and reading the old ones.  We try to save
  //  a little space by removing the old store files as soon as we're
  //  done with them.
  //

  mValid = AS_OVS_readOverlapFromStore(merge, &movl);

  for (i=1; i<numberOfOrigStoreFiles; i++) {
    sprintf(name, "%s/%04d~", storeName, i);
    BinaryOverlapFile  *bof = AS_OVS_createBinaryOverlapFile(name, TRUE, FALSE);

    sValid = AS_OVS_readOverlap(bof, &movl);

    while (sValid) {
      if ((mValid) &&
          (OVSoverlap_sort(&sovl, &movl) > 0)) {
        AS_OVS_writeOverlapToStore(store, &movl);
        mValid = AS_OVS_readOverlapFromStore(merge, &movl);
      } else {
        AS_OVS_writeOverlapToStore(store, &sovl);
        sValid = AS_OVS_readOverlap(bof, &sovl);
      }
    }

    AS_OVS_closeBinaryOverlapFile(bof);

    if (saveSpace) {
      sprintf(name, "%04d", i);
      nukeBackup(storeName, name);
    }
  }

  //  All the original data is in, we now just need to copy any
  //  remaining new data into the new store.
  //
  while (mValid) {
    AS_OVS_writeOverlapToStore(store, &movl);
    mValid = AS_OVS_readOverlapFromStore(merge, &movl);
  }

  //  ALL DONE!  Close the stores, nuke the backups and get outta here.
  //
  AS_OVS_closeOverlapStore(merge);
  AS_OVS_closeOverlapStore(store);

  nukeBackup(storeName, "ovs");
  nukeBackup(storeName, "idx");
  for (i=1; i<numberOfOrigStoreFiles; i++) {
    sprintf(name, "%04d", i);
    nukeBackup(storeName, name);
  }

  exit(0);
}
