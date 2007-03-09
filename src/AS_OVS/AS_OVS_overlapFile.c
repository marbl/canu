
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

//  XXX  split this into create() and open()

BinaryOverlapFile *
AS_OVS_createBinaryOverlapFile(const char *name, int isInternal, int isOutput) {

  BinaryOverlapFile   *bof = (BinaryOverlapFile *)safe_malloc(sizeof(BinaryOverlapFile));

  bof->bufferLen  = 0;
  bof->bufferPos  = 16384 * 12;
  bof->bufferMax  = 16384 * 12;
  bof->buffer     = (uint32 *)safe_malloc(sizeof(uint32) * bof->bufferMax);
  bof->isOutput   = isOutput;
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

  if (isOutput) {
    errno = 0;
    if (name == NULL)
      bof->file = stdout;
    else
      bof->file = fopen(name, "w");
    if (errno) {
      fprintf(stderr, "AS_OVS_createBinaryOverlapFile()-- Failed to open '%s' for writing: %s\n",
              name, strerror(errno));
      exit(1);
    }
  } else {
    errno = 0;
    if (name == NULL)
      bof->file = stdin;
    else
      bof->file = fopen(name, "r");
    if (errno) {
      fprintf(stderr, "AS_OVS_createBinaryOverlapFile()-- Failed to open '%s' for reading: %s\n",
              name, strerror(errno));
      exit(1);
    }
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

  fclose(bof->file);

  safe_free(bof->buffer);
  safe_free(bof);
}


void
AS_OVS_writeOverlap(BinaryOverlapFile *bof, OVSoverlap *overlap) {

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
