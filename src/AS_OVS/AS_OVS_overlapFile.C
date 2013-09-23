
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

static const char *rcsid = "$Id$";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

#include "AS_OVS_overlapFile.H"
#include "AS_UTL_fileIO.H"

int AS_OVS_BOF_BufferMax = 1 * 1024 * 1024;


void
AS_OVS_setBinaryOverlapFileBufferSize(int size) {
  if (size < 16 * 1024)
    size = 16 * 1024;

  AS_OVS_BOF_BufferMax = size;
}


int
AS_OVS_getBinaryOverlapFileBufferSize(void) {
  return(AS_OVS_BOF_BufferMax);
}


static void AS_OVS_initializeBOF(BinaryOverlapFile *bof, int isInternal, int isOutput) {
  assert(bof != NULL);

  uint32  lcf = (AS_OVS_NWORDS + 1) * (AS_OVS_NWORDS + 2);

  bof->bufferLen  = 0;
  bof->bufferPos  = (AS_OVS_BOF_BufferMax / (lcf * sizeof(uint32))) * lcf;
  bof->bufferMax  = (AS_OVS_BOF_BufferMax / (lcf * sizeof(uint32))) * lcf;
  bof->buffer     = (uint32 *)safe_malloc(sizeof(uint32) * bof->bufferMax);
  bof->isOutput   = isOutput;
  bof->isSeekable = FALSE;
  bof->isInternal = isInternal;
  bof->reader     = NULL;
  bof->writer     = NULL;
  bof->file       = NULL;

  //  The size of the buffer MUST be divisible by our two overlap sizes, otherwise our writer will
  //  lose data.  We carefully chose a size that also gives us a 4MB buffer.

  if ((bof->bufferMax % (AS_OVS_NWORDS + 1)) != 0)
    fprintf(stderr, "ERROR: bufferMax="F_U32" AS_OVS_NWORDS=%d\n", bof->bufferMax, AS_OVS_NWORDS);
  assert((bof->bufferMax % (AS_OVS_NWORDS + 1)) == 0);

  if ((bof->bufferMax % (AS_OVS_NWORDS + 2)) != 0)
    fprintf(stderr, "ERROR: bufferMax="F_U32" AS_OVS_NWORDS=%d\n", bof->bufferMax, AS_OVS_NWORDS);
  assert((bof->bufferMax % (AS_OVS_NWORDS + 2)) == 0);
}

BinaryOverlapFile *
AS_OVS_openBinaryOverlapFile(const char *name, int isInternal) {
  BinaryOverlapFile   *bof = (BinaryOverlapFile *)safe_malloc(sizeof(BinaryOverlapFile));

  AS_OVS_initializeBOF(bof, isInternal, FALSE);

  bof->reader = new compressedFileReader(name);
  bof->file   = bof->reader->file();

  if (bof->reader->isCompressed() == false)
    bof->isSeekable = true;

  return(bof);
}




BinaryOverlapFile *
AS_OVS_createBinaryOverlapFile(const char *name, int isInternal) {
  BinaryOverlapFile   *bof = (BinaryOverlapFile *)safe_malloc(sizeof(BinaryOverlapFile));

  AS_OVS_initializeBOF(bof, isInternal, TRUE);

  bof->writer = new compressedFileWriter(name);
  bof->file   = bof->writer->file();

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

  delete bof->reader;
  delete bof->writer;

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
  bof->buffer[bof->bufferLen++] = overlap->dat.dat[0];
  bof->buffer[bof->bufferLen++] = overlap->dat.dat[1];
#if AS_OVS_NWORDS > 2
  bof->buffer[bof->bufferLen++] = overlap->dat.dat[2];
#endif
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

  overlap->b_iid      = bof->buffer[bof->bufferPos++];
  overlap->dat.dat[0] = bof->buffer[bof->bufferPos++];
  overlap->dat.dat[1] = bof->buffer[bof->bufferPos++];
#if AS_OVS_NWORDS > 2
  overlap->dat.dat[2] = bof->buffer[bof->bufferPos++];
#endif
  assert(bof->bufferPos <= bof->bufferLen);

  return(TRUE);
}


void
AS_OVS_seekOverlap(BinaryOverlapFile *bof, uint32 overlap) {

  assert(bof->isSeekable == TRUE);

  //  Move to the correct spot, and force a load on the next readOverlap by setting the position to
  //  the end of the buffer.

  AS_UTL_fseek(bof->file,
               (off_t)overlap * sizeof(uint32) * ((bof->isInternal) ? (AS_OVS_NWORDS + 1) : (AS_OVS_NWORDS + 2)),
               SEEK_SET);
  bof->bufferPos = bof->bufferLen;
}
