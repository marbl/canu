
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

static const char *rcsid = "$Id: AS_OVS_overlapStore.c,v 1.19 2008-10-08 22:02:58 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <fcntl.h>
#include <assert.h>

#include "AS_OVS_overlapStore.h"
#include "AS_OVS_overlapFile.h"
#include "AS_UTL_fileIO.h"



static
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

static
void
renameFromBackup(char const *storeName, char const *name) {
  char   orig[FILENAME_MAX];
  char   bkup[FILENAME_MAX];

  sprintf(orig, "%s/%s", storeName, name);
  sprintf(bkup, "%s/%s~", storeName, name);

  errno = 0;
  rename(bkup, orig);
  if (errno) {
    fprintf(stderr, "overlapStore: ERROR: failed to restore backup of '%s' into '%s': %s\n", bkup, orig, strerror(errno));
    assert(0);
  }
}


static
void
nukeBackup(char const *storeName, char const *name) {
  char   bkup[FILENAME_MAX];

  sprintf(bkup, "%s/%s~", storeName, name);

  errno = 0;
  unlink(bkup);
  if ((errno) && (errno != ENOENT))
    fprintf(stderr, "overlapStore: WARNING: failed to remove backup '%s': %s\n", bkup, strerror(errno));
}




OverlapStore *
AS_OVS_openOverlapStorePrivate(const char *path, int useBackup, int saveSpace) {
  char            name[FILENAME_MAX];
  FILE           *ovsinfo;

  OverlapStore   *ovs = (OverlapStore *)safe_calloc(1, sizeof(OverlapStore));

  //  Overlap store cannot be from stdin!
  assert((path != NULL) && (strcmp(path, "-") != 0));

  strcpy(ovs->storePath, path);

  ovs->isOutput  = FALSE;
  ovs->useBackup = (useBackup) ? '~' : 0;
  ovs->saveSpace = saveSpace;

  ovs->ovs.ovsMagic              = 1;
  ovs->ovs.ovsVersion            = 1;
  ovs->ovs.numOverlapsPerFile    = 0;  //  not used for reading
  ovs->ovs.smallestIID           = 1000000000;
  ovs->ovs.largestIID            = 0;
  ovs->ovs.numOverlapsTotal      = 0;
  ovs->ovs.highestFileIndex      = 0;

  sprintf(name, "%s/ovs", path);
  errno = 0;
  ovsinfo = fopen(name, "r");
  if (errno) {
    fprintf(stderr, "failed to open info file '%s': %s\n", name, strerror(errno));
    exit(1);
  }
  AS_UTL_safeRead(ovsinfo, &ovs->ovs, "AS_OVS_openOverlapStore info", sizeof(OverlapStoreInfo), 1);
  fclose(ovsinfo);


  //  If we're not supposed to be using the backup, load the stats.
  //
#if 0
  if (ovs->useBackup == 0) {
    FILE *ost;

    sprintf(name, "%s/ost", ovs->storePath);
    errno = 0;
    ost = fopen(name, "r");
    if (errno) {
      fprintf(stderr, "failed to open the stats file '%s': %s\n", name, strerror(errno));
      exit(1);
    }

    AS_UTL_safeRead(ost, &ovs->stats, "AS_OVS_openOverlapStore", sizeof(OverlapStoreStats), 1);
    fclose(ost);
  }
#endif

  //  If we're supposed to be using the backup, actually make the
  //  backup.  OK, it's not much of a "backup", it's really just a
  //  stashing the current store somewhere, so we can recreate this
  //  store with merged in data.
  //
  if (ovs->useBackup) {
    int i;
    renameToBackup(ovs->storePath, "ovs");
    renameToBackup(ovs->storePath, "idx");
    for (i=1; i<=ovs->ovs.highestFileIndex; i++) {
      sprintf(name, "%04d", i);
      renameToBackup(ovs->storePath, name);
    }
  }


  sprintf(name, "%s/idx%c", path, ovs->useBackup);
  errno = 0;
  ovs->offsetFile      = fopen(name, "r");
  ovs->offset.a_iid    = 0;
  ovs->offset.fileno   = 0;
  ovs->offset.offset   = 0;
  ovs->offset.numOlaps = 0;
  if (errno) {
    fprintf(stderr, "AS_OVS_openOverlapStore()-- failed to open offset file '%s': %s\n", name, strerror(errno));
    exit(1);
  }

  ovs->missing.a_iid    = 0;
  ovs->missing.fileno   = 0;
  ovs->missing.offset   = 0;
  ovs->missing.numOlaps = 0;


  ovs->firstIIDrequested = ovs->ovs.smallestIID;
  ovs->lastIIDrequested  = ovs->ovs.largestIID;

  ovs->overlapsThisFile = 0;
  ovs->currentFileIndex = 1;

  sprintf(name, "%s/%04d%c", ovs->storePath, ovs->currentFileIndex, ovs->useBackup);
  ovs->bof = AS_OVS_openBinaryOverlapFile(name, TRUE);

  return(ovs);
}

void
AS_OVS_restoreBackup(OverlapStore *ovs) {
  char            name[FILENAME_MAX];

  //  Restore the backup for an overlap store.
  //
  if (ovs->useBackup) {
    int i;
    renameFromBackup(ovs->storePath, "ovs");
    renameFromBackup(ovs->storePath, "idx");
    for (i=1; i<=ovs->ovs.highestFileIndex; i++) {
      sprintf(name, "%04d", i);
      renameFromBackup(ovs->storePath, name);
    }
  }
}

int
AS_OVS_readOverlapFromStore(OverlapStore *ovs, OVSoverlap *overlap, uint32 type) {

  assert(ovs->isOutput == FALSE);

  //  If we've finished reading overlaps for the current a_iid, get
  //  another a_iid.  If we hit EOF here, we're all done, no more
  //  overlaps.
  //
 again:
  while (ovs->offset.numOlaps == 0)
    if (0 == AS_UTL_safeRead(ovs->offsetFile, &ovs->offset, "AS_OVS_readOverlap offset",
                             sizeof(OverlapStoreOffsetRecord), 1))
      return(FALSE);

  //  And if we've exited the range of overlaps requested, return.
  //
  if (ovs->offset.a_iid > ovs->lastIIDrequested)
    return(FALSE);

  while (AS_OVS_readOverlap(ovs->bof, overlap) == FALSE) {
    char name[FILENAME_MAX];

    //  We read no overlap, open the next file and try again.

    AS_OVS_closeBinaryOverlapFile(ovs->bof);

    if (ovs->saveSpace) {
      sprintf(name, "%04d", ovs->currentFileIndex);
      nukeBackup(ovs->storePath, name);
    }

    ovs->currentFileIndex++;

    sprintf(name, "%s/%04d%c", ovs->storePath, ovs->currentFileIndex, ovs->useBackup);
    ovs->bof = AS_OVS_openBinaryOverlapFile(name, TRUE);

    //  AS_OVS_openBinaryOverlapFile() actually bombs if it can't open; this test is useless....
    if (ovs->bof == NULL) {
      fprintf(stderr, "AS_OVS_readOverlapFromStore()-- failed to open overlap file '%s': %s\n", name, strerror(errno));
      exit(1);
    }
  }

  overlap->a_iid   = ovs->offset.a_iid;

  ovs->offset.numOlaps--;

  if (type == AS_OVS_TYPE_ANY)
    return(TRUE);

  if (type != overlap->dat.ovl.type)
    goto again;

  return(TRUE);
}




void
AS_OVS_setRangeOverlapStore(OverlapStore *ovs, uint32 firstIID, uint32 lastIID) {
  char            name[FILENAME_MAX];

  //  make the index be one record per read iid, regardless, then we
  //  can quickly grab the correct record, and seek to the start of
  //  those overlaps

  if (firstIID >= ovs->ovs.largestIID)
    firstIID = ovs->ovs.largestIID;
  if (lastIID >= ovs->ovs.largestIID)
    lastIID = ovs->ovs.largestIID;

  //  If our range is invalid (firstIID > lastIID) we keep going, and
  //  let AS_OVS_readOverlapFromStore() deal with it.

  AS_UTL_fseek(ovs->offsetFile, (size_t)firstIID * sizeof(OverlapStoreOffsetRecord), SEEK_SET);

  //  Unfortunately, we need to actually read the record to figure out
  //  where to position the overlap stream.  If the read fails, we
  //  silently return, letting AS_OVS_readOverlapFromStore() deal with
  //  the problem.

  ovs->offset.a_iid    = 0;
  ovs->offset.fileno   = 0;
  ovs->offset.offset   = 0;
  ovs->offset.numOlaps = 0;

  if (0 == AS_UTL_safeRead(ovs->offsetFile, &ovs->offset, "AS_OVS_readOverlap offset",
                           sizeof(OverlapStoreOffsetRecord), 1))
    return;

  ovs->overlapsThisFile = 0;
  ovs->currentFileIndex = ovs->offset.fileno;

  AS_OVS_closeBinaryOverlapFile(ovs->bof);

  sprintf(name, "%s/%04d%c", ovs->storePath, ovs->currentFileIndex, ovs->useBackup);
  ovs->bof = AS_OVS_openBinaryOverlapFile(name, TRUE);

  AS_OVS_seekOverlap(ovs->bof, ovs->offset.offset);

  ovs->firstIIDrequested = firstIID;
  ovs->lastIIDrequested  = lastIID;
}



void
AS_OVS_resetRangeOverlapStore(OverlapStore *ovs) {
  char            name[FILENAME_MAX];

  rewind(ovs->offsetFile);

  ovs->offset.a_iid    = 0;
  ovs->offset.fileno   = 0;
  ovs->offset.offset   = 0;
  ovs->offset.numOlaps = 0;

  ovs->overlapsThisFile = 0;
  ovs->currentFileIndex = 1;

  AS_OVS_closeBinaryOverlapFile(ovs->bof);

  sprintf(name, "%s/%04d%c", ovs->storePath, ovs->currentFileIndex, ovs->useBackup);
  ovs->bof = AS_OVS_openBinaryOverlapFile(name, TRUE);

  ovs->firstIIDrequested = ovs->ovs.smallestIID;
  ovs->lastIIDrequested  = ovs->ovs.largestIID;
}





////////////////////////////////////////////////////////////////////////////////


void
AS_OVS_closeOverlapStore(OverlapStore *ovs) {
  char name[FILENAME_MAX];

  if (ovs == NULL)
    return;

  if (ovs->useBackup) {
    int i;
    nukeBackup(ovs->storePath, "ovs");
    nukeBackup(ovs->storePath, "idx");
    for (i=1; i<=ovs->ovs.highestFileIndex; i++) {
      sprintf(name, "%04d", i);
      nukeBackup(ovs->storePath, name);
    }
  }

  if (ovs->isOutput) {
    FILE *ovsinfo = NULL;

    //  Write the last index element, maybe, and don't forget to fill
    //  in gaps!
    //
    if (ovs->offset.numOlaps > 0) {
      while (ovs->missing.a_iid < ovs->offset.a_iid) {
        ovs->missing.fileno    = ovs->offset.fileno;
        ovs->missing.offset    = ovs->offset.offset;
        ovs->missing.numOlaps  = 0;
        AS_UTL_safeWrite(ovs->offsetFile,
                         &ovs->missing,
                         "AS_OVS_closeOverlapStore offset",
                         sizeof(OverlapStoreOffsetRecord),
                         1);
        ovs->missing.a_iid++;
      }

      AS_UTL_safeWrite(ovs->offsetFile, &ovs->offset, "AS_OVS_closeOverlapStore offset",
                       sizeof(OverlapStoreOffsetRecord), 1);
    }

    //  Update the info
    //

    sprintf(name, "%s/ovs", ovs->storePath);
    errno = 0;
    ovsinfo = fopen(name, "w");
    if (errno) {
      fprintf(stderr, "failed to create overlap store '%s': %s\n", ovs->storePath, strerror(errno));
      exit(1);
    }
    ovs->ovs.highestFileIndex = ovs->currentFileIndex;
    AS_UTL_safeWrite(ovsinfo, &ovs->ovs, "AS_OVS_closeOverlapStore", sizeof(OverlapStoreInfo), 1);
    fclose(ovsinfo);

    fprintf(stderr, "Closing the new store:\n");
    fprintf(stderr, "ovs->ovs.ovsMagic           = "F_U64"\n", ovs->ovs.ovsMagic);
    fprintf(stderr, "ovs->ovs.ovsVersion         = "F_U64"\n", ovs->ovs.ovsVersion);
    fprintf(stderr, "ovs->ovs.numOverlapsPerFile = "F_U64"\n", ovs->ovs.numOverlapsPerFile);
    fprintf(stderr, "ovs->ovs.smallestIID        = "F_U64"\n", ovs->ovs.smallestIID);
    fprintf(stderr, "ovs->ovs.largestIID         = "F_U64"\n", ovs->ovs.largestIID);
    fprintf(stderr, "ovs->ovs.numOverlapsTotal   = "F_U64"\n", ovs->ovs.numOverlapsTotal);
    fprintf(stderr, "ovs->ovs.highestFileIndex   = "F_U64"\n", ovs->ovs.highestFileIndex);
  }

#if 0
  if (ovs->statsUpdated) {
    FILE *ost;

    fprintf(stderr, "Writing new stats.\n");

    sprintf(name, "%s/ost", ovs->storePath);
    errno = 0;
    ost = fopen(name, "w");
    if (errno) {
      fprintf(stderr, "failed to write overlap stats '%s': %s\n", name, strerror(errno));
      exit(1);
    }

    AS_UTL_safeWrite(ost, &ovs->stats, "AS_OVS_closeOverlapStore", sizeof(OverlapStoreStats), 1);
    fclose(ost);
  }
#endif

  AS_OVS_closeBinaryOverlapFile(ovs->bof);

  fclose(ovs->offsetFile);
  safe_free(ovs);
}


////////////////////////////////////////////////////////////////////////////////


//  Create a new overlap store.  By default, the new
//  store is write-only.
//
OverlapStore *
AS_OVS_createOverlapStore(const char *path, int failOnExist) {
  char            name[FILENAME_MAX];
  FILE           *ovsinfo;

  assert((path != NULL) && (strcmp(path, "-") != 0));

  OverlapStore   *ovs = (OverlapStore *)safe_calloc(1, sizeof(OverlapStore));
  strcpy(ovs->storePath, path);

  ovs->isOutput  = TRUE;
  ovs->useBackup = 0;
  ovs->saveSpace = 0;

  AS_UTL_mkdir(path);

  sprintf(name, "%s/ovs", path);
  errno = 0;
  ovsinfo = fopen(name, "w");
  if (errno) {
    fprintf(stderr, "failed to create overlap store '%s': %s\n", ovs->storePath, strerror(errno));
    exit(1);
  }
  ovs->ovs.ovsMagic              = 1;
  ovs->ovs.ovsVersion            = 1;
  ovs->ovs.numOverlapsPerFile    = 1024 * 1024 * 1024 / sizeof(OVSoverlapINT);
  ovs->ovs.smallestIID           = 1000000000;
  ovs->ovs.largestIID            = 0;
  ovs->ovs.numOverlapsTotal      = 0;
  ovs->ovs.highestFileIndex      = 0;
  AS_UTL_safeWrite(ovsinfo, &ovs->ovs, "AS_OVS_createOverlapStore", sizeof(OverlapStoreInfo), 1);
  fclose(ovsinfo);


  sprintf(name, "%s/idx", path);
  errno = 0;
  ovs->offsetFile      = fopen(name, "w");
  ovs->offset.a_iid    = 0;
  ovs->offset.fileno   = 0;
  ovs->offset.offset   = 0;
  ovs->offset.numOlaps = 0;
  if (errno) {
    fprintf(stderr, "AS_OVS_createOverlapStore()-- failed to open offset file '%s': %s\n", name, strerror(errno));
    exit(1);
  }


  ovs->missing.a_iid    = 0;
  ovs->missing.fileno   = 0;
  ovs->missing.offset   = 0;
  ovs->missing.numOlaps = 0;


  ovs->overlapsThisFile = 0;
  ovs->currentFileIndex = 0;
  ovs->bof              = NULL;


  return(ovs);
}





void
AS_OVS_writeOverlapToStore(OverlapStore *ovs, OVSoverlap *overlap) {
  char            name[FILENAME_MAX];

  assert(ovs->isOutput == TRUE);

  if (ovs->offset.a_iid > overlap->a_iid) {
    //  Woah!  The last overlap we saw is bigger than the one we have now?!
    fprintf(stderr, "LAST:  a:"F_U32"\n", ovs->offset.a_iid);
    fprintf(stderr, "THIS:  a:"F_U32" b:"F_U32"\n", overlap->a_iid, overlap->b_iid);
  }
  assert(ovs->offset.a_iid <= overlap->a_iid);

  if (ovs->ovs.smallestIID > overlap->a_iid)
    ovs->ovs.smallestIID = overlap->a_iid;
  if (ovs->ovs.largestIID < overlap->a_iid)
     ovs->ovs.largestIID = overlap->a_iid;

  //  If we don't have an output file yet, or the current file is
  //  too big, open a new file.
  //
  if (ovs->overlapsThisFile >= ovs->ovs.numOverlapsPerFile) {
    AS_OVS_closeBinaryOverlapFile(ovs->bof);

    ovs->bof              = NULL;
    ovs->overlapsThisFile = 0;
  }
  if (ovs->bof == NULL) {
    char  name[FILENAME_MAX];

    ovs->currentFileIndex++;

    sprintf(name, "%s/%04d", ovs->storePath, ovs->currentFileIndex);
    ovs->bof = AS_OVS_createBinaryOverlapFile(name, TRUE);
  }


  //  Put the index to disk, filling any gaps
  //
  if ((ovs->offset.numOlaps != 0) &&
      (ovs->offset.a_iid != overlap->a_iid)) {

    while (ovs->missing.a_iid < ovs->offset.a_iid) {
      ovs->missing.fileno    = ovs->offset.fileno;
      ovs->missing.offset    = ovs->offset.offset;
      ovs->missing.numOlaps  = 0;
      AS_UTL_safeWrite(ovs->offsetFile,
                       &ovs->missing,
                       "AS_OVS_writeOverlapToStore offset",
                       sizeof(OverlapStoreOffsetRecord),
                       1);
      ovs->missing.a_iid++;
    }

    //  One more, since this iid is not missing -- we write it next!
    ovs->missing.a_iid++;

    AS_UTL_safeWrite(ovs->offsetFile,
                     &ovs->offset,
                     "AS_OVS_writeOverlapToStore offset",
                     sizeof(OverlapStoreOffsetRecord),
                     1);
    ovs->offset.numOlaps  = 0;
  }


  //  Update the index if this is the first time we've seen this a_iid.
  //
  if (ovs->offset.numOlaps == 0) {
    ovs->offset.a_iid     = overlap->a_iid;
    ovs->offset.fileno    = ovs->currentFileIndex;
    ovs->offset.offset    = ovs->overlapsThisFile;
  }

  //AS_OVS_accumulateStats(ovs, overlap);
  AS_OVS_writeOverlap(ovs->bof, overlap);
  ovs->offset.numOlaps++;
  ovs->ovs.numOverlapsTotal++;
  ovs->overlapsThisFile++;
}




uint64
AS_OVS_numOverlapsInRange(OverlapStore *ovs) {
  size_t                     originalposition = 0;
  uint64                     i = 0;
  uint64                     len = 0;
  OverlapStoreOffsetRecord  *offsets = NULL;
  uint64                     numolap = 0;

  originalposition = AS_UTL_ftell(ovs->offsetFile);

  AS_UTL_fseek(ovs->offsetFile, (size_t)ovs->firstIIDrequested * sizeof(OverlapStoreOffsetRecord), SEEK_SET);

  //  Even if we're doing a whole human-size store, this allocation is
  //  (a) temporary and (b) only 512MB.  The only current consumer of
  //  this code is FragCorrectOVL.c, which doesn't run on the whole
  //  human, it runs on ~24 pieces, which cuts this down to < 32MB.

  len = ovs->lastIIDrequested - ovs->firstIIDrequested + 1;
  offsets = (OverlapStoreOffsetRecord *)safe_malloc(sizeof(OverlapStoreOffsetRecord) * len);

  if (len != AS_UTL_safeRead(ovs->offsetFile, offsets, "AS_OVS_numOverlapsInRange",
                             sizeof(OverlapStoreOffsetRecord), len)) {
    fprintf(stderr, "AS_OVS_numOverlapsInRange()-- short read on offsets!\n");
    exit(1);
  }

  for (i=0; i<len; i++)
    numolap += offsets[i].numOlaps;

  safe_free(offsets);

  AS_UTL_fseek(ovs->offsetFile, originalposition, SEEK_SET);

  return(numolap);
}
