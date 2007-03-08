
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "AS_OVS_overlapStore.h"
#include "AS_OVS_overlapFile.h"
#include "AS_UTL_fileIO.h"

OverlapStore *
AS_OVS_openOverlapStore(const char *path) {
  char            name[FILENAME_MAX];
  FILE           *ovsinfo;

  OverlapStore   *ovs = (OverlapStore *)safe_malloc(sizeof(OverlapStore));

  //  Overlap store cannot be from stdin!
  assert((path != NULL) && (strcmp(path, "-") != 0));

  strcpy(ovs->storePath, path);

  ovs->ovs.ovsMagic              = 1;
  ovs->ovs.ovsVersion            = 1;
  ovs->ovs.numOverlapsPerFile    = 0;
  ovs->ovs.smallestIID           = 0;
  ovs->ovs.largestIID            = 0;
  ovs->ovs.numOverlapsTotal      = 0;

  sprintf(name, "%s/ovs", path);
  errno = 0;
  ovsinfo = fopen(name, "r");
  if (errno) {
    fprintf(stderr, "failed to open info file '%s': %s\n", name, strerror(errno));
    exit(1);
  }
  AS_UTL_safeRead(ovsinfo, &ovs->ovs, "AS_OVS_openOverlapStore info", sizeof(OverlapStoreInfo), 1);
  fclose(ovsinfo);


  sprintf(name, "%s/idx", path);
  errno = 0;
  ovs->offsetFile      = fopen(name, "r");
  ovs->offset.a_iid    = 0;
  ovs->offset.offset   = 0;
  ovs->offset.numOlaps = 0;
  if (errno) {
    fprintf(stderr, "AS_OVS_openOverlapStore()-- failed to open offset file '%s': %s\n", name, strerror(errno));
    exit(1);
  }

  ovs->bufferLen        = 0;
  ovs->bufferPos        = 1048576;
  ovs->bufferMax        = 1048576;
  ovs->buffer           = (OVSoverlapINT *)safe_malloc(sizeof(OVSoverlapINT) * ovs->bufferMax);
  ovs->isOutput         = FALSE;

  ovs->overlapsThisFile = 0;
  ovs->currentFileIndex = 1;

  sprintf(name, "%s/%04d", ovs->storePath, ovs->currentFileIndex);
  errno = 0;
  ovs->file = fopen(name, "r");
  if (errno) {
    fprintf(stderr, "AS_OVS_openOverlapStore()-- failed to open overlap file '%s': %s\n", name, strerror(errno));
    exit(1);
  }

  return(ovs);
}



int
AS_OVS_readOverlapFromStore(OverlapStore *ovs, OVSoverlap *overlap) {

  assert(ovs->isOutput == FALSE);

  //  If we've finished reading overlaps for the current a_iid, get
  //  another a_iid.

  if (ovs->offset.numOlaps == 0) {

    //  Read the next offset.  If we're at the end of file, there must
    //  be no more overlaps.  Might as well check that the overlap
    //  file is at the end too...

    if (0 == AS_UTL_safeRead(ovs->offsetFile,
                             &ovs->offset,
                             "AS_OVS_readOverlap offset",
                             sizeof(OverlapStoreOffsetRecord),
                             1)) {
      ovs->bufferLen = AS_UTL_safeRead(ovs->file,
                                       ovs->buffer,
                                       "AS_OVS_readOverlap",
                                       sizeof(OVSoverlapINT),
                                       ovs->bufferMax);
      if (ovs->bufferLen > 0) {
        fprintf(stderr, "AS_OVS_readOverlap()--  WARNING!  The index file is at the end,\n");
        fprintf(stderr, "AS_OVS_readOverlap()--  indicating no more overlaps, but the overlap\n");
        fprintf(stderr, "AS_OVS_readOverlap()--  file had more overlaps in it!\n");
        fprintf(stderr, "AS_OVS_readOverlap()--  Something is wrong, and I'll just stop here.\n");
        assert(ovs->bufferLen == 0);
      }

      return(FALSE);
    }

    //fprintf(stdout, "READ new offset: A:%d Num:%d\n", ovs->offset.a_iid, ovs->offset.numOlaps);
  }


  while (ovs->bufferPos >= ovs->bufferLen) {
    ovs->bufferPos = 0;
    ovs->bufferLen = AS_UTL_safeRead(ovs->file,
                                     ovs->buffer,
                                     "AS_OVS_readOverlap",
                                     sizeof(OVSoverlapINT),
                                     ovs->bufferMax);

    //  If we read no overlaps, open the next file and try again.
    //
    if (ovs->bufferLen == 0) {
      char name[FILENAME_MAX];

      fclose(ovs->file);

      ovs->currentFileIndex++;

      sprintf(name, "%s/%04d", ovs->storePath, ovs->currentFileIndex);
      errno = 0;
      ovs->file = fopen(name, "r");
      if (errno) {
        fprintf(stderr, "AS_OVS_openOverlapStore()-- failed to open overlap file '%s': %s\n", name, strerror(errno));
        exit(1);
      }
    }
  }


  overlap->a_iid   = ovs->offset.a_iid;
  overlap->b_iid   = ovs->buffer[ovs->bufferPos].b_iid;
  overlap->dat.dat = ovs->buffer[ovs->bufferPos].dat.dat;

  ovs->bufferPos++;
  ovs->offset.numOlaps--;

  return(TRUE);
}


////////////////////////////////////////////////////////////////////////////////


void
AS_OVS_closeOverlapStore(OverlapStore *ovs) {

  if (ovs == NULL)
    return;

  //  Flush the current buffer to disk.
  //
  if ((ovs->isOutput == TRUE) && (ovs->bufferLen > 0)) {
    AS_UTL_safeWrite(ovs->file, ovs->buffer, "AS_OVS_closeOverlapStore",
                     sizeof(OVSoverlapINT), ovs->bufferLen);
    AS_UTL_safeWrite(ovs->offsetFile, &ovs->offset, "AS_OVS_writeOverlapToStore offset",
                     sizeof(OverlapStoreOffsetRecord), 1);
  }

  //  Update the info
  //


  fclose(ovs->file);
  fclose(ovs->offsetFile);
  safe_free(ovs->buffer);
  safe_free(ovs);
}


////////////////////////////////////////////////////////////////////////////////


//  Create a new overlap store.  By default, the new
//  store is write-only.
//
OverlapStore *
AS_OVS_createOverlapStore(const char *path) {
  char            name[FILENAME_MAX];
  FILE           *ovsinfo;

  OverlapStore   *ovs = (OverlapStore *)safe_malloc(sizeof(OverlapStore));

  assert((path != NULL) && (strcmp(path, "-") != 0));

  strcpy(ovs->storePath, path);

  errno = 0;
  mkdir(path, S_IRWXU | S_IRWXG | S_IROTH);
  if (errno) {
    fprintf(stderr, "CreateOverlapStore(): failed to create directory '%s': %s\n", ovs->storePath, strerror(errno));
    exit(1);
  }


  sprintf(name, "%s/ovs", path);
  errno = 0;
  ovsinfo = fopen(name, "w");
  if (errno) {
    fprintf(stderr, "failed to create overlap store '%s': %s\n", ovs->storePath, strerror(errno));
    exit(1);
  }
  ovs->ovs.ovsMagic              = 1;
  ovs->ovs.ovsVersion            = 1;
  ovs->ovs.numOverlapsPerFile    = 0;
  ovs->ovs.smallestIID           = 0;
  ovs->ovs.largestIID            = 0;
  ovs->ovs.numOverlapsTotal      = 0;
  AS_UTL_safeWrite(ovsinfo, &ovs->ovs, "AS_OVS_createOverlapStore", sizeof(OverlapStoreInfo), 1);
  fclose(ovsinfo);


  sprintf(name, "%s/idx", path);
  errno = 0;
  ovs->offsetFile      = fopen(name, "w");
  ovs->offset.a_iid    = 0;
  ovs->offset.offset   = 0;
  ovs->offset.numOlaps = 0;
  if (errno) {
    fprintf(stderr, "AS_OVS_createOverlapStore()-- failed to open offset file '%s': %s\n", name, strerror(errno));
    exit(1);
  }


  ovs->bufferLen        = 0;
  ovs->bufferPos        = 1048576;
  ovs->bufferMax        = 1048576;
  ovs->buffer           = (OVSoverlapINT *)safe_malloc(sizeof(OVSoverlapINT) * ovs->bufferMax);
  ovs->isOutput         = TRUE;


  ovs->overlapsThisFile = 0;
  ovs->currentFileIndex = 0;
  ovs->file             = NULL;


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
    exit(1);
  }

  if (ovs->bufferLen >= ovs->bufferMax) {

    //  If we don't have an output file yet, or the current file it
    //  too big, open a new file.
    //
    if (ovs->overlapsThisFile > ovs->ovs.numOverlapsPerFile) {
      fclose(ovs->file);

      ovs->file             = NULL;
      ovs->overlapsThisFile = 0;
    }
    if (ovs->file == NULL) {
      char  name[FILENAME_MAX];

      ovs->currentFileIndex++;

      sprintf(name, "%s/%04d", ovs->storePath, ovs->currentFileIndex);
      ovs->file = fopen(name, "w");
      if (errno) {
        fprintf(stderr, "AS_OVS_writeOverlap()-- Failed to open '%s' for writing: %s\n",
                name, strerror(errno));
        exit(1);
      }
    }

    //
    //  OK, got a file.  Now put stuff in it.
    //

    AS_UTL_safeWrite(ovs->file, ovs->buffer, "AS_OVS_writeOverlap", sizeof(OVSoverlap), ovs->bufferLen);

    ovs->overlapsThisFile += ovs->bufferLen;

    ovs->bufferLen = 0;

    //  XXXX should we bother making sure the file has exactly
    //  overlapsThisFile in it?  Easy but kind of ugly.
  }

  //  Update the index, and maybe put it out to disk
  //
  if ((ovs->offset.numOlaps != 0) && (ovs->offset.a_iid != overlap->a_iid)) {
    AS_UTL_safeWrite(ovs->offsetFile,
                     &ovs->offset,
                     "AS_OVS_writeOverlapToStore offset",
                     sizeof(OverlapStoreOffsetRecord),
                     1);
    ovs->offset.numOlaps  = 0;
  }
  if (ovs->offset.numOlaps == 0) {
    ovs->offset.a_iid     = overlap->a_iid;
    ovs->offset.fileno    = ovs->currentFileIndex;
    ovs->offset.offset    = ovs->overlapsThisFile;
  }

  ovs->buffer[ovs->bufferLen].b_iid   = overlap->b_iid;
  ovs->buffer[ovs->bufferLen].dat.dat = overlap->dat.dat;

  ovs->bufferLen++;
  ovs->offset.numOlaps++;
}
