
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

#include "AS_OVS_overlapFile.h"
#include "AS_UTL_fileIO.h"


BinaryOverlapFile *
AS_OVS_openBinaryOverlapFile(const char *name, int isInternal) {
  char     cmd[1024 + FILENAME_MAX];

  BinaryOverlapFile   *bof = (BinaryOverlapFile *)safe_malloc(sizeof(BinaryOverlapFile));

  bof->bufferLen  = 0;
  bof->bufferPos  = 16384 * 12;
  bof->bufferMax  = 16384 * 12;
  bof->buffer     = (uint32 *)safe_malloc(sizeof(uint32) * bof->bufferMax);
  bof->isOutput   = FALSE;
  bof->isSeekable = FALSE;
  bof->isPopened  = FALSE;
  bof->isInternal = isInternal;
  bof->file       = NULL;

  //  The size of the buffer MUST be divisible by 3 and 4, otherwise
  //  our writer will lose data.  We carefully chose 16384*12 to be
  //  close to a 1MB buffer.
  //
  assert((bof->bufferMax % 3) == 0);
  assert((bof->bufferMax % 4) == 0);

  if ((name == NULL) || (strcmp(name, "-") == 0))
    name = NULL;

  errno = 0;
  if        (name == NULL) {
    bof->file = stdin;
  } else if (strcmp(name+strlen(name)-3, ".gz") == 0) {
    sprintf(cmd, "gzip -dc %s", name);
    bof->file = popen(cmd, "r");
    bof->isPopened = TRUE;
  } else if (strcmp(name+strlen(name)-4, ".bz2") == 0) {
    sprintf(cmd, "bzip2 -dc %s", name);
    bof->file = popen(cmd, "r");
    bof->isPopened = TRUE;
  } else {
    bof->file       = fopen(name, "r");
    bof->isSeekable = TRUE;
  }
  if (errno) {
    fprintf(stderr, "AS_OVS_createBinaryOverlapFile()-- Failed to open '%s' for reading: %s\n",
            name, strerror(errno));
    exit(1);
  }

  return(bof);
}




BinaryOverlapFile *
AS_OVS_createBinaryOverlapFile(const char *name, int isInternal) {
  char     cmd[1024 + FILENAME_MAX];

  BinaryOverlapFile   *bof = (BinaryOverlapFile *)safe_malloc(sizeof(BinaryOverlapFile));

  bof->bufferLen  = 0;
  bof->bufferPos  = 16384 * 12;
  bof->bufferMax  = 16384 * 12;
  bof->buffer     = (uint32 *)safe_malloc(sizeof(uint32) * bof->bufferMax);
  bof->isOutput   = TRUE;
  bof->isSeekable = FALSE;
  bof->isPopened  = FALSE;
  bof->isInternal = isInternal;
  bof->file       = NULL;

  //  The size of the buffer MUST be divisible by 3 and 4, otherwise
  //  our writer will lose data.  We carefully chose 16384*12 to be
  //  close to a 1MB buffer.
  //
  assert((bof->bufferMax % 3) == 0);
  assert((bof->bufferMax % 4) == 0);

  if ((name == NULL) || (strcmp(name, "-") == 0))
    name = NULL;

  errno = 0;
  if (name == NULL) {
    bof->file = stdout;
  } else if (strcmp(name+strlen(name)-3, ".gz") == 0) {
    sprintf(cmd, "gzip -9c %s", name);
    bof->file = popen(cmd, "w");
    bof->isPopened = TRUE;
  } else if (strcmp(name+strlen(name)-4, ".bz2") == 0) {
    sprintf(cmd, "bzip2 -9c %s", name);
    bof->file = popen(cmd, "w");
    bof->isPopened = TRUE;
  } else {
    bof->file = fopen(name, "w");
  }
  if (errno) {
    fprintf(stderr, "AS_OVS_createBinaryOverlapFile()-- Failed to open '%s' for writing: %s\n",
            name, strerror(errno));
    exit(1);
  }

  return(bof);
}



void
AS_OVS_flushBinaryOverlapFile(BinaryOverlapFile *bof) {
  if ((bof->isOutput) && (bof->bufferLen > 0)) {
    AS_UTL_safeWrite(bof->file, bof->buffer, "AS_OVS_flushBinaryOverlapFile", sizeof(uint32), bof->bufferLen);
    bof->bufferLen = 0;
  }
}


void
AS_OVS_closeBinaryOverlapFile(BinaryOverlapFile *bof) {

  if (bof == NULL)
    return;

  if (bof->isOutput)
    AS_OVS_flushBinaryOverlapFile(bof);

  if (bof->isPopened)
    pclose(bof->file);
  else
    fclose(bof->file);

  safe_free(bof->buffer);
  safe_free(bof);
}


void
AS_OVS_writeOverlap(BinaryOverlapFile *bof, OVSoverlap *overlap) {

  assert(bof->isOutput == TRUE);

  if (bof->bufferLen >= bof->bufferMax) {
    AS_UTL_safeWrite(bof->file, bof->buffer, "AS_OVS_outputOverlap", sizeof(uint32), bof->bufferLen);
    bof->bufferLen = 0;
  }

  if (bof->isInternal == FALSE)
    bof->buffer[bof->bufferLen++] = overlap->a_iid;

  bof->buffer[bof->bufferLen++] = overlap->b_iid;
  bof->buffer[bof->bufferLen++] = (overlap->dat.dat >> 32) & 0xffffffff;
  bof->buffer[bof->bufferLen++] = (overlap->dat.dat)       & 0xffffffff;
}


int
AS_OVS_readOverlap(BinaryOverlapFile *bof, OVSoverlap *overlap) {

  assert(bof->isOutput == FALSE);

  if (bof->bufferPos >= bof->bufferLen) {
    bof->bufferPos = 0;
    bof->bufferLen = AS_UTL_safeRead(bof->file,
                                     bof->buffer,
                                     "AS_OVS_readOverlap",
                                     sizeof(uint32),
                                     bof->bufferMax);
  }

  if (bof->bufferPos >= bof->bufferLen)
    return(FALSE);

  if (bof->isInternal == FALSE)
    overlap->a_iid = bof->buffer[bof->bufferPos++];

  overlap->b_iid     = bof->buffer[bof->bufferPos++];
  overlap->dat.dat   = bof->buffer[bof->bufferPos++];
  overlap->dat.dat <<= 32;
  overlap->dat.dat  |= bof->buffer[bof->bufferPos++];

  return(TRUE);
}


void
AS_OVS_seekOverlap(BinaryOverlapFile *bof, uint32 overlap) {

  assert(bof->isSeekable == TRUE);

  //  Move to the correct spot, and force a load on the next read

  CDS_FSEEK(bof->file,
            (off_t)overlap * sizeof(uint32) * ((bof->isInternal) ? 3 : 4),
            SEEK_SET);
  bof->bufferPos = bof->bufferLen;
}
