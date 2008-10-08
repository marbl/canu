
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

static char *rcsid = "$Id: AS_PER_genericStore.c,v 1.31 2008-10-08 22:02:58 brianwalenz Exp $";

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
  if (s->fp) {
    if ((s->memoryBuffer) && (s->isDirty)) {
      assert(s->readOnly == FALSE);
      rewind(s->fp);
      AS_UTL_safeWrite(s->fp,
                       s->memoryBuffer, "writeUpdate",
                       sizeof(char),
                       sizeof(StoreStruct) + (s->lastElem - s->firstElem + 1) * s->elementSize);
    }

    //  If we're dirty write the header.  Yes, we overwrite the header
    //  stored in the memory buffer.
    //
    if (s->isDirty) {
      rewind(s->fp);
      AS_UTL_safeWrite(s->fp, s, "writeHeader", sizeof(StoreStruct), 1);
    }

    if (fclose(s->fp) != 0) {
      fprintf(stderr, "Failed to close the store; this usually means your disk is full.\n");
      exit(1);
    }
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

  if (path) {
    assert(strlen(path) < FILENAME_MAX);

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
  } else {
    s = convertStoreToMemoryStore(s);
  }

  return(s);
}

StoreStruct *
createStringStore(const char *path, const char *storeLabel) {
  StoreStruct *s = (StoreStruct *)safe_calloc(sizeof(StoreStruct), 1);

  //  If path is NULL, the newly created store is an in-core store, with no
  //  disk backing.

  strncpy(s->storeLabel, storeLabel, 7);

  s->storeType       = STRING_STORE;
  s->elementSize     = sizeof(char);
  s->firstElem       = 1;
  s->lastElem        = 0;

  s->allocatedSize   = 0;
  s->memoryBuffer    = NULL;
  s->fp              = NULL;
  s->diskBuffer      = NULL;
  s->isDirty         = 0;
  s->readOnly        = 0;
  s->lastWasWrite    = 0;

  if (path) {
    assert(strlen(path) < FILENAME_MAX);

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
  } else {
    s = convertStoreToMemoryStore(s);
  }

  return(s);
}

////////////////////////////////////////////////////////////////////////////////

void
getIndexStore(StoreStruct *s, int64 index, void *buffer) {
  int64 offset = computeOffset(s,index);

  assert(s->storeType == INDEX_STORE);
  assert(index > 0);

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
  assert(s->storeType == INDEX_STORE);
  assert(s->memoryBuffer);
  assert(index > 0);
  if (index > s->lastElem)
    return(NULL);
  return((void *)(s->memoryBuffer + computeOffset(s,index)));
}


void
getStringStore(StoreStruct *s, int64 offset, char *buffer, uint32 maxLength, uint32 *actualLength, int64 *nextOffset) {
  uint32 length = 0;
  int64 actualOffset;

  buffer[0] = 0;

  assert(s->storeType == STRING_STORE);
  assert(offset > 0);

  actualOffset = offset - s->firstElem;

  assert(actualOffset <= s->lastElem);

  if (s->memoryBuffer) {
    memcpy(&length, s->memoryBuffer + actualOffset + sizeof(StoreStruct), sizeof(uint32));

    assert(length <= maxLength);
    assert(length + actualOffset + sizeof(uint32) <= s->lastElem);

    if (length > 0)
      memcpy(buffer, s->memoryBuffer + actualOffset + sizeof(StoreStruct) + sizeof(uint32), length + 1);
  } else {
    AS_UTL_fseek(s->fp, (off_t) (actualOffset + sizeof(StoreStruct)), SEEK_SET);
    if (s->lastWasWrite)
      fflush(s->fp);

    if (1 != AS_UTL_safeRead(s->fp,&length,"getStringStore",sizeof(uint32), 1)) {
      fprintf(stderr, "getStringStore()-- Failed to read the length of the record.  Incomplete store?\n");
      exit(1);
    }

    assert(length <= maxLength);
    assert(length + actualOffset + sizeof(uint32) <= s->lastElem);

    if (length > 0) {
      if (length + 1 != AS_UTL_safeRead(s->fp,buffer,"getStringStore",sizeof(char), length + 1)) {
        fprintf(stderr, "getStringStore()-- Failed to read all "F_U32" bytes.  Incomplete store?\n", length);
        exit(1);
      }
    }

    s->lastWasWrite = 0;
  }

  //  Careful!  Precedence of ? sucks.
  *nextOffset   = offset + sizeof(uint32) + ((length > 0) ? length + 1 : 0);
  *actualLength = length;
}

char *
getStringStorePtr(StoreStruct *s, int64 offset, uint32 *actualLength, int64 *nextOffset) {
  assert(s->memoryBuffer);
  assert(s->storeType == STRING_STORE);
  assert(offset > 0);

  if (offset + sizeof(uint32) > s->lastElem) {
    *actualLength = 0;
    *nextOffset   = 0;
    return(NULL);
  }

  memcpy(actualLength, s->memoryBuffer + sizeof(StoreStruct) + offset - s->firstElem, sizeof(uint32));

  //  Not exactly ideal, but we don't store long stuff here anyway.
  //  There's no real reason for asserting this, just in the hope that
  //  if there is corruption of some goofy crud, we'll get a bogus
  //  length.
  //
  assert(*actualLength <= 1048576);
  assert(offset - s->firstElem + *actualLength + sizeof(uint32) <= s->lastElem);

  //  Careful!  Precedence of ? sucks.
  *nextOffset = offset + sizeof(uint32) + ((*actualLength > 0) ? *actualLength + 1 : 0);

  return(s->memoryBuffer + sizeof(StoreStruct) + offset - s->firstElem + sizeof(uint32));
}


////////////////////////////////////////////////////////////////////////////////

void
setIndexStore(StoreStruct *s, int64 index, void *element) {

  assert(s->readOnly == FALSE);
  assert(s->firstElem <= index);
  assert(s->lastElem >= index);

  assert(index > 0);

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

int64
appendStringStore(StoreStruct *s, char *str, uint32 len) {

  assert(s->readOnly == FALSE);
  assert(s->storeType == STRING_STORE);

  int64 retval = 0;
  int64 offset = sizeof(StoreStruct) + s->lastElem;

  //  Make sure the string is zero terminated.
  if ((str) && (len > 0))
    assert(str[len] == 0);

  if (s->memoryBuffer) {
    int64 desiredSize = offset + len + sizeof(uint32);
    if (s->allocatedSize <= desiredSize) {
      char *oldbuffer = s->memoryBuffer;

#if 0
      //
      //  Useful for playing with string UID reallocation bugs.
      //
      fprintf(stderr, "REALLOCATE from %d to %d\n",
              s->allocatedSize, desiredSize + 2 * 1024 * 1024);
      size_t oldsize = s->allocatedSize;
      s->allocatedSize = desiredSize + 2 * 1024 * 1024;
      char *n = safe_calloc(s->allocatedSize, 1);
      memcpy(n, s->memoryBuffer, oldsize);
      safe_free(s->memoryBuffer);
      s->memoryBuffer = n;
#else
      s->allocatedSize = desiredSize + 32 * 1024 * 1024;
      s->memoryBuffer  = (char *)safe_realloc(s->memoryBuffer, s->allocatedSize);
#endif

      retval = s->memoryBuffer - oldbuffer;
    }
    memcpy(s->memoryBuffer + offset, &len, sizeof(uint32));
    if (len > 0)
      memcpy(s->memoryBuffer + offset + sizeof(uint32), str, len + 1);
  } else {
    AS_UTL_fseek(s->fp, (off_t) offset, SEEK_SET);
    if (s->lastWasWrite == 0)
      fflush(s->fp);

    AS_UTL_safeWrite(s->fp, &len, "appendStringStore", sizeof(uint32), 1);
    if (len > 0)
      AS_UTL_safeWrite(s->fp, str, "appendStringStore", sizeof(char), len + 1);
    s->lastWasWrite = 1;
  }

  //  Just painful.  If the string was empty, we wrote nothing.
  //  Otherwise, we wrote len+1 bytes.
  //
  //  Careful!  Precedence of ? sucks.
  s->lastElem += sizeof(uint32) + ((len > 0) ? len + 1 : 0);
  s->isDirty   = 1;

  return(retval);
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
  assert(startIndex != 0);
  assert(endIndex   != 0);
  ss->startIndex = (startIndex == STREAM_FROMSTART) ? ss->store->firstElem : startIndex;
  ss->endIndex   = (endIndex   == STREAM_UNTILEND)  ? ss->store->lastElem  : endIndex;
}

int
nextStream(StreamStruct *ss, void *buffer, uint32 maxLength, uint32 *actualLength) {
  uint32 nextOffset = 0;

  if (ss->startIndex > ss->endIndex)
    return(0);
  if (ss->store->storeType == INDEX_STORE) {
    getIndexStore(ss->store, ss->startIndex, buffer);
    ss->startIndex++;
  } else {
    getStringStore(ss->store, ss->startIndex, buffer, maxLength, actualLength, &ss->startIndex);
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

  //  These should be true for a disk-based store, and false for a
  //  memory store.

  assert(source->firstElem    == 1);
  assert(source->memoryBuffer == NULL);

  source->allocatedSize   = sizeof(StoreStruct) + (source->lastElem - source->firstElem + 1) * source->elementSize;
  source->memoryBuffer    = (char *)safe_calloc(sizeof(char), source->allocatedSize);

  //  Copy the data to the memory store.  NOTE that even though we
  //  copy the header, it is UNUSED -- we use the values stored in the
  //  object itself.

  if (source->fp) {
    fflush(source->fp);
    rewind(source->fp);

    if (AS_UTL_safeRead(source->fp,
                        source->memoryBuffer,
                        "convertStoreMemoryStore",
                        sizeof(char),
                        source->allocatedSize) != source->allocatedSize) {
      fprintf(stderr, "convertStoreToMemoryStore()-- failed to convert store (label '%s') to memory store.\n",
              source->storeLabel);
      exit(1);
    }
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
  int64        bytesRead       = 0;

  if (firstElem <= 0)
    firstElem = source->firstElem;

  if (lastElem <= 0)
    lastElem = source->lastElem;

  //  These should be true for a disk-based store, and false for a
  //  memory store.

  assert(source->firstElem    == 1);
  assert(source->memoryBuffer == NULL);

  if (source->storeType == STRING_STORE) {
    sourceOffset    = sizeof(StoreStruct) + firstElem - source->firstElem;
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
  //
  if (source->fp) {
    fflush(source->fp);
    AS_UTL_fseek(source->fp, (off_t)sourceOffset, SEEK_SET);

    bytesRead = AS_UTL_safeRead(source->fp,
                                target->memoryBuffer + sizeof(StoreStruct),
                                "convertStoreToPartialMemoryStore",
                                sizeof(char),
                                sourceMaxOffset - sourceOffset);
    if (bytesRead != sourceMaxOffset - sourceOffset) {
      fprintf(stderr, "convertStoreToPartialMemoryStore()-- failed to convert store (label '%s') to memory store.\n", source->storeLabel);
      fprintf(stderr, "convertStoreToPartialMemoryStore()-- wanted to read %d bytes, actually read %d.\n", sourceMaxOffset - sourceOffset, bytesRead);
      exit(1);
    }
  }

  closeStore(source);

  return(target);
}
