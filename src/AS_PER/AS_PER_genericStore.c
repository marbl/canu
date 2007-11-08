
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

static char CM_ID[] = "$Id: AS_PER_genericStore.c,v 1.22 2007-11-08 12:38:15 brianwalenz Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <string.h>
#include <assert.h>
#include <errno.h>

#include "AS_global.h"
#include "AS_PER_genericStore.h"
#include "AS_UTL_fileIO.h"

static int INITIAL_ALLOCATION = 4096;
static int WRITING_BUFFER     = 8 * 1024;

#define INDEX_STORE    1
#define STRING_STORE   2

void
AS_PER_setBufferSize(int wb) {
  WRITING_BUFFER = wb;
}


#define  computeOffset(S, I)  (((I) - (S)->firstElem) * (int64)((S)->elementSize) + sizeof(StoreStruct))



StoreStruct *
openStore(const char *path, const char *rw) {
  StoreStruct *s = (StoreStruct *)safe_calloc(sizeof(StoreStruct), 1);

  errno = 0;
  FILE *fp = fopen(path, rw);
  if (errno) {
    fprintf(stderr,"openStore()-- Failed to open store %s (%s): %s\n", path, rw, strerror(errno));
    exit(1);
  }

  if (1 != AS_UTL_safeRead(fp, s, "openStore", sizeof(StoreStruct), 1)) {
    fprintf(stderr, "openStore()-- Failed to read the header for store %s.\n");
    exit(1);
  }

  s->allocatedSize   = 0;
  s->memoryBuffer    = NULL;
  s->fp              = fp;
  s->diskBuffer      = NULL;
  s->isDirty         = 0;
  s->readOnly        = 0;
  s->lastWasWrite    = 0;

  if (strcmp(rw, "r") == 0)
    s->readOnly = 1;

  return(s);
}


void
closeStore(StoreStruct *s) {

  //  Flush the memory-resident store to our disk backing.
  //
  if ((s->memoryBuffer) && (s->isDirty)) {
    assert(s->readOnly == FALSE);
    rewind(s->fp);
    AS_UTL_safeWrite(s->fp,
                     s->memoryBuffer, "writeUpdate",
                     sizeof(char),
                     sizeof(StoreStruct) + (s->lastElem - s->firstElem) * s->elementSize);
  }

  //  If we're dirty write the header.  Yes, we overwrite the header
  //  stored in the memory buffer.
  //
  if (s->isDirty) {
    rewind(s->fp);
    AS_UTL_safeWrite(s->fp, s, "writeHeader", sizeof(StoreStruct), 1);
  }

  if ((s->fp) && (fclose(s->fp) != 0)) {
    fprintf(stderr, "Failed to close the store; this usually means your disk is full.\n");
    exit(1);
  }

  safe_free(s->memoryBuffer);
  safe_free(s->diskBuffer);
  memset(s, 0xfe, sizeof(StoreStruct));
  safe_free(s);
}


////////////////////////////////////////////////////////////////////////////////

StoreStruct *
createIndexStore(const char *path, const char *storeLabel, int32 elementSize, int64 firstID) {
  StoreStruct *s = (StoreStruct *)safe_calloc(sizeof(StoreStruct), 1);

  assert(path != NULL);
  assert(strlen(path) < FILENAME_MAX);

  strncpy(s->storeLabel, storeLabel, 7);

  s->storeType       = INDEX_STORE;
  s->elementSize     = elementSize;
  s->firstElem       = firstID;
  s->lastElem        = firstID - 1;

  s->allocatedSize   = 0;
  s->memoryBuffer    = NULL;
  s->fp              = NULL;
  s->diskBuffer      = NULL;
  s->isDirty         = 0;
  s->readOnly        = 0;
  s->lastWasWrite    = 0;

  errno = 0;
  s->fp = fopen(path, "w+");
  if (errno) {
    fprintf(stderr,"createIndexStore()-- Failure opening Store %s for w+: %s\n", path, strerror(errno));
    exit(1);
  }

  if (WRITING_BUFFER > 0)
    setbuffer(s->fp,
              s->diskBuffer = (char *)safe_calloc(sizeof(char), WRITING_BUFFER),
              WRITING_BUFFER);

  AS_UTL_safeWrite(s->fp, s, "writeHeader", sizeof(StoreStruct), 1);

  return(s);
}

StoreStruct *
createStringStore(const char *path, const char *storeLabel) {
  StoreStruct *s = (StoreStruct *)safe_calloc(sizeof(StoreStruct), 1);

  assert(path != NULL);
  assert(strlen(path) < FILENAME_MAX);

  strncpy(s->storeLabel, storeLabel, 7);

  s->storeType       = STRING_STORE;
  s->elementSize     = sizeof(char);
  s->firstElem       = 0;
  s->lastElem        = 0;

  s->allocatedSize   = 0;
  s->memoryBuffer    = NULL;
  s->fp              = NULL;
  s->diskBuffer      = NULL;
  s->isDirty         = 0;
  s->readOnly        = 0;
  s->lastWasWrite    = 0;

  errno = 0;
  s->fp = fopen(path, "w+");
  if (errno) {
    fprintf(stderr,"createStringStore()-- Failure opening Store %s for w+: %s\n", path, strerror(errno));
    exit(1);
  }

  if (WRITING_BUFFER > 0)
    setbuffer(s->fp,
              s->diskBuffer = (char *)safe_calloc(sizeof(char), WRITING_BUFFER),
              WRITING_BUFFER);

  AS_UTL_safeWrite(s->fp, s, "writeHeader", sizeof(StoreStruct), 1);

  return(s);
}

////////////////////////////////////////////////////////////////////////////////

void
getIndexStore(StoreStruct *s, int64 index, void *buffer) {
  int64 offset = computeOffset(s,index);

  if (s->memoryBuffer) {
    memcpy(buffer, s->memoryBuffer + offset, s->elementSize);
  } else {
    AS_UTL_fseek(s->fp, (off_t)offset, SEEK_SET);
    if (s->lastWasWrite)
      fflush(s->fp);

    if (1 != AS_UTL_safeRead(s->fp,buffer,"getIndexStore",s->elementSize, 1)) {
      fprintf(stderr, "getIndexStore()-- Failed to read the record.  Incomplete store?\n");
      exit(1);
    }
    s->lastWasWrite = 0;
  }
}


void *
getIndexStorePtr(StoreStruct *s, int64 index) {
  assert(s->memoryBuffer);
  if (index > s->lastElem)
    return(NULL);
  return((void *)(s->memoryBuffer + computeOffset(s,index)));
}


void
getStringStore(StoreStruct *s, int64 offset, void *buffer, uint32 maxLength, uint32 *actualLength) {
  uint32 length = 0;
  int64 actualOffset;

  assert(s->storeType == STRING_STORE);
  assert(offset <= s->lastElem);

  actualOffset = offset - s->firstElem;

  if (s->memoryBuffer) {
    memcpy(&length, s->memoryBuffer + actualOffset + sizeof(StoreStruct), sizeof(uint32));
    assert((length + actualOffset <= s->lastElem) && (length <= maxLength));

    if (length > 0)
      memcpy(buffer, s->memoryBuffer + actualOffset + sizeof(StoreStruct) + sizeof(uint32), length);
  } else {
    AS_UTL_fseek(s->fp, (off_t) (actualOffset + sizeof(StoreStruct)), SEEK_SET);
    if (s->lastWasWrite)
      fflush(s->fp);

    if (1 != AS_UTL_safeRead(s->fp,&length,"getStringStore",sizeof(uint32), 1)) {
      fprintf(stderr, "getStringStore()-- Failed to read the length of the record.  Incomplete store?\n");
      exit(1);
    }

    assert((length + offset + sizeof(uint32) <= s->lastElem) && (length <= maxLength));

    if (length != AS_UTL_safeRead(s->fp,buffer,"getStringStore",sizeof(char), length)) {
      fprintf(stderr, "getStringStore()-- Failed to read all "F_U32" bytes.  Incomplete store?\n", length);
      exit(1);
    }

    s->lastWasWrite = 0;
  }

  *actualLength = length;
}

char *
getStringStorePtr(StoreStruct *s, int64 offset, uint32 *actualLength) {
  assert(s->memoryBuffer);
  assert(s->storeType == STRING_STORE);

  if (offset > s->lastElem) {
    *actualLength = 0;
    return(NULL);
  }

  memcpy(actualLength, s->memoryBuffer + offset - s->firstElem + sizeof(StoreStruct), sizeof(uint32));
  return(s->memoryBuffer + offset - s->firstElem + sizeof(StoreStruct) + sizeof(uint32));
}


////////////////////////////////////////////////////////////////////////////////

void
setIndexStore(StoreStruct *s, int64 index, void *element) {

  assert(s->readOnly == FALSE);
  assert(s->firstElem <= index);
  assert(s->lastElem >= index);

  if (s->memoryBuffer) {
    memcpy(s->memoryBuffer + computeOffset(s, index), element, s->elementSize);
  } else {
    AS_UTL_fseek(s->fp, (off_t)computeOffset(s, index), SEEK_SET);
    if (s->lastWasWrite == 0)
      fflush(s->fp);

    AS_UTL_safeWrite(s->fp, element, "setIndexStore", s->elementSize, 1);
    s->lastWasWrite = 1;
  }

  s->isDirty = 1;
}

////////////////////////////////////////////////////////////////////////////////

void
appendIndexStore(StoreStruct *s, void *element) {

  if (s->memoryBuffer) {
    int64 desiredSize = computeOffset(s, s->lastElem + 1) + (s->lastElem - s->firstElem + 1) * s->elementSize;

    if (s->allocatedSize <= desiredSize) {
      s->allocatedSize = desiredSize + 32 * 1024 * 1024;
      s->memoryBuffer  = (char *)safe_realloc(s->memoryBuffer, s->allocatedSize);
    }
  }

  setIndexStore(s, ++s->lastElem, element);
}

void
appendStringStore(StoreStruct *s, void *str, uint32 len) {

  assert(s->readOnly == FALSE);
  assert(s->storeType == STRING_STORE);

  int64 offset = sizeof(StoreStruct) + s->lastElem;

  if (str)
    len++;  //  write the null byte

  if (s->memoryBuffer) {
    int64 desiredSize = offset + len + sizeof(uint32);
    if (s->allocatedSize <= desiredSize) {
      s->allocatedSize = desiredSize + 32 * 1024 * 1024;
      s->memoryBuffer  = (char *)safe_realloc(s->memoryBuffer, s->allocatedSize);
    }
    memcpy(s->memoryBuffer + offset, &len, sizeof(uint32));
    memcpy(s->memoryBuffer + offset + sizeof(uint32), str, len);
  } else {
    AS_UTL_fseek(s->fp, (off_t) offset, SEEK_SET);
    if (s->lastWasWrite == 0)
      fflush(s->fp);

    AS_UTL_safeWrite(s->fp, &len, "appendStringStore", sizeof(uint32), 1);
    AS_UTL_safeWrite(s->fp, str, "appendStringStore", sizeof(char), len);
    s->lastWasWrite = 1;
  }

  s->lastElem += len + sizeof(uint32);
  s->isDirty = 1;
}


////////////////////////////////////////////////////////////////////////////////

StreamStruct *
openStream(StoreStruct *fs) {
  StreamStruct *ss = (StreamStruct *)safe_calloc(sizeof(StreamStruct), 1);
  ss->store       = fs;
  ss->startIndex  = fs->firstElem;
  ss->endIndex    = fs->lastElem;
  return(ss);
}

void
resetStream(StreamStruct *ss, int64 startIndex, int64 endIndex) {
  ss->startIndex = (startIndex == STREAM_FROMSTART) ? ss->store->firstElem : startIndex;
  ss->endIndex   = (endIndex   == STREAM_UNTILEND)  ? ss->store->lastElem  : endIndex;
}

int
nextStream(StreamStruct *ss, void *buffer, uint32 maxLength, uint32 *actualLength) {
  if (ss->startIndex > ss->endIndex)
    return(0);
  if (ss->store->storeType == INDEX_STORE) {
    getIndexStore(ss->store, ss->startIndex, buffer);
    ss->startIndex++;
  } else {
    getStringStore(ss->store, ss->startIndex, buffer, maxLength, actualLength);
    ss->startIndex += *actualLength + sizeof(uint32);
  }
  return(1);
}

void
closeStream(StreamStruct *ss) {
  memset(ss, 0xfe, sizeof(StreamStruct));
  safe_free(ss);
}


////////////////////////////////////////////////////////////////////////////////

StoreStruct *
convertStoreToMemoryStore(StoreStruct *source) {

  assert(source->allocatedSize == 0);
  assert(source->memoryBuffer == NULL);

  source->allocatedSize   = sizeof(StoreStruct) + (source->lastElem - source->firstElem) * source->elementSize;
  source->memoryBuffer    = (char *)safe_calloc(sizeof(char), source->allocatedSize);

  fflush(source->fp);
  rewind(source->fp);

  //  Copy the data to the memory store.  NOTE that even though we
  //  copy the header, it is UNUSED -- we use the values stored in the
  //  object itself.

  if (AS_UTL_safeRead(source->fp,
                      source->memoryBuffer,
                      "convertStoreMemoryStore",
                      sizeof(char),
                      source->allocatedSize) != source->allocatedSize) {
    fprintf(stderr, "convertStoreToMemoryStore()-- failed to convert store (label '%s') to memory store.\n",
            source->storeLabel);
    exit(1);
  }

  return(source);
}


StoreStruct *
convertStoreToPartialMemoryStore(StoreStruct *source,
                                 int64         firstElem,
                                 int64         lastElem) {
  StoreStruct *target         = NULL;
  int64        sourceOffset    = 0;
  int64        sourceMaxOffset = 0;

  if (firstElem <= 0)
    firstElem = source->firstElem;

  if (lastElem <= 0)
    lastElem = source->lastElem;

  if (source->storeType == STRING_STORE) {
    sourceOffset    = sizeof(StoreStruct) + firstElem;
    sourceMaxOffset = sizeof(StoreStruct) + lastElem;
  } else {
    sourceOffset    = computeOffset(source, firstElem);
    sourceMaxOffset = computeOffset(source, lastElem + 1);
  }

  target  = (StoreStruct *)safe_calloc(sizeof(StoreStruct), 1);
  *target = *source;

  //  This partial memory store is READ ONLY.  We do not allow writes
  //  because the writing code doesn't support it, because nobody
  //  needs to use it yet, and because it gets complicated -- the
  //  partial store cannot allow writes past what is loaded, unless
  //  what is loaded is at the end of the store.  Consider a partial
  //  store of the first 100 elements out of 1000, vs the first 100
  //  elements out of 100.  Writing element 101 to the former will
  //  cause grief, but to the latter is ok.

  target->allocatedSize   = sizeof(StoreStruct) + sourceMaxOffset - sourceOffset;
  target->memoryBuffer    = (char *)safe_calloc(sizeof(char), target->allocatedSize);
  target->fp              = NULL;
  target->diskBuffer      = NULL;
  target->readOnly        = 1;
  target->isDirty         = 0;
  target->lastWasWrite    = 0;

  target->firstElem       = firstElem;
  target->lastElem        = lastElem;

  //  Copy the data to the memory store.  The first
  //  sizeof(StoreStruct) bytes are unused -- this would be the
  //  header.

  AS_UTL_fseek(source->fp, (off_t)sourceOffset, SEEK_SET);
  fflush(source->fp);
  if (AS_UTL_safeRead(source->fp,
                      target->memoryBuffer + sizeof(StoreStruct),
                      "convertStoreToPartialMemoryStore",
                      sizeof(char),
                      sourceMaxOffset - sourceOffset) != sourceMaxOffset - sourceOffset) {
    fprintf(stderr, "convertStoreToPartialMemoryStore()-- failed to convert store (label '%s') to memory store.\n",
            source->storeLabel);
    exit(1);
  }

  closeStore(source);

  return(target);
}
