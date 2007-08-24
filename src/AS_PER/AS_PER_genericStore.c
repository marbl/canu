
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

static char CM_ID[] = "$Id: AS_PER_genericStore.c,v 1.21 2007-08-24 15:29:48 brianwalenz Exp $";

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

#define  computeOffset(S, I)  (((I) - (S)->header.firstElem) * (int64)((S)->header.elementSize) + sizeof(StoreStat))

static
void
writeHeader(StoreStruct *s){
  s->header.lastUpdateTime = time(0);
  if (s->isMemoryStore) {
    memcpy(s->buffer,&s->header, sizeof(StoreStat));
  } else {
    rewind(s->fp);
    AS_UTL_safeWrite(s->fp, &s->header, "writeHeader", sizeof(StoreStat), 1);
    //fwrite(&s->header, sizeof(StoreStat), 1, s->fp);
  }
}


static
StoreStruct  *
allocateStore(void){
  StoreStruct *ss = (StoreStruct *)safe_calloc(sizeof(StoreStruct), 1);

  ss->fp                    = NULL;
  ss->buffer                = NULL;
  ss->allocatedSize         = 0;
  ss->status                = UnInitializedStore;
  ss->header.firstElem      = -1;
  ss->header.lastElem       = -1;
  ss->lastCommittedElem     = -1;

  ss->header.isDeleted      = 0;
  ss->header.type           = 0;
  ss->header.p1             = 0;
  ss->header.p2             = 0;
  ss->header.storeType[0]   = 0;
  ss->header.storeType[1]   = 0;
  ss->header.storeType[2]   = 0;
  ss->header.storeType[3]   = 0;
  ss->header.storeType[4]   = 0;
  ss->header.storeType[5]   = 0;
  ss->header.storeType[6]   = 0;
  ss->header.storeType[7]   = 0;
  ss->header.firstElem      = -1;
  ss->header.lastElem       = -1;
  ss->header.elementSize    = -1;
  ss->header.creationTime   = 0;
  ss->header.lastUpdateTime = 0;

  ss->lastCommittedElem     = -1;
  ss->isMemoryStore         = 0;
  ss->isDirty               = 0;
  ss->lastWasRead           = 1;

  return(ss);
}


static
StreamStruct  *
allocateStream(StoreStruct *s, void *buffer, int32 bufferSize  ){
  StreamStruct *ss = (StreamStruct *)safe_calloc(sizeof(StreamStruct), 1);

  ss->store       = s;
  ss->status      = ActiveStore;
  ss->buffer      = buffer;
  ss->bufferSize  = bufferSize;
  ss->startIndex  = s->header.firstElem;
  ss->endIndex    = s->header.lastElem;

  return(ss);
}


void
AS_PER_setBufferSize(int wb) {
  WRITING_BUFFER = wb;
}


StoreStruct *
openStore(const char *path, const char *rw) {
  StoreStruct *s = allocateStore();

  errno = 0;
  s->fp     = fopen(path,rw);
  s->status = ActiveStore;

  if(errno){
    fprintf(stderr,"openStore()-- Failed to open store %s (%s): %s\n",
            path, rw, strerror(errno));
    exit(1);
  }

  if (1 != AS_UTL_safeRead(s->fp, (void *)&s->header, "openStore", sizeof(StoreStat), 1)) {
    fprintf(stderr, "openStore()-- Failed to read the header for store %s.\n");
    exit(1);
  }

  return s;
}


void
statsStore(StoreStruct *s, StoreStat *stats){
  *(stats) = s->header;
}


void
closeStore(StoreStruct *s){
  assert(s->status == ActiveStore);
  if (s->isMemoryStore){
    s->allocatedSize = 0;
  } else {
    if (s->isDirty)
      writeHeader(s);
    if (fclose(s->fp) != 0) {
      fprintf(stderr, "Failed to close the store; this usually means your disk is full.\n");
      exit(1);
    }
  }
  safe_free(s->buffer);
  s->status = UnAllocatedStore;
}


int64
getLastElemStore(StoreStruct *s){
  return s->header.lastElem;
}


int64
getFirstElemStore(StoreStruct *s){
  return s->header.firstElem;
}


StoreStruct *
createIndexStore(const char *path, const char *storeType, int32 elementSize, int64 firstID) {
  StoreStruct *s = allocateStore();

  s->header.type = INDEX_STORE;
  strncpy(s->header.storeType, storeType,7);
  s->header.elementSize = elementSize;
  s->header.firstElem = firstID;
  s->header.lastElem =  firstID -1; /* Sentinel */
  s->lastCommittedElem = -1;
  s->status = UnInitializedStore;
  s->header.creationTime = time(0);
  s->status = ActiveStore;

  if(path){ /* File based store */
    assert(strlen(path) < FILENAME_MAX);
    errno = 0;
    s->fp = fopen(path, "w+");
    if(errno){
      fprintf(stderr,"createIndexStore()-- Failure opening Store %s for w+: %s\n", path, strerror(errno));
      assert(errno == 0);
    }
    s->buffer = NULL;
    if (WRITING_BUFFER > 0) {
      //fprintf(stderr, "createIndexStore()--  Allocating fp buffer of "F_SIZE_T" bytes for '%s'.\n", WRITING_BUFFER, path);
      s->buffer = (char *)safe_calloc(sizeof(char), WRITING_BUFFER);
      setbuffer(s->fp, s->buffer, WRITING_BUFFER);
    }
    s->isMemoryStore = 0;
  }else{
    s->allocatedSize = INITIAL_ALLOCATION * elementSize;
    s->buffer = (char *)safe_calloc(sizeof(char), s->allocatedSize);
    s->isMemoryStore = 1;
    s->fp = NULL;
  }
  writeHeader(s); 
  s->isDirty = 0;
  return s;
}


void
getIndexStore(StoreStruct *s, int64 index, void *buffer){
  int64 offset;

  assert(s->status == ActiveStore);
  offset = computeOffset(s,index);

  if(s->isMemoryStore){
    memcpy(buffer, s->buffer + offset, s->header.elementSize);
    //return(0);
  }else{
    AS_UTL_fseek(s->fp, (off_t) offset, SEEK_SET);

    if (s->lastWasRead == 0)
      fflush(s->fp);
    s->lastWasRead = 1;

    assert(1 == AS_UTL_safeRead(s->fp,buffer,"getIndexStore",s->header.elementSize, 1));
  }
}


void *
getIndexStorePtr(StoreStruct *s, int64 index) {

  assert(s->status == ActiveStore);
  assert(s->isMemoryStore);

  return((void *)(s->buffer + computeOffset(s,index)));
}


void
setIndexStore(StoreStruct *s, int64 index, void *element){
  int64 offset;

  assert(s->status == ActiveStore);
  assert(s->header.firstElem <= index);
  assert(s->header.lastElem >= index);

  s->isDirty = 1;

  offset = computeOffset(s, index);

  if (s->isMemoryStore) {
    memcpy(s->buffer + offset, element, s->header.elementSize);
  } else {
    if (s->lastWasRead)
      fflush(s->fp);
    s->lastWasRead = 0;

    AS_UTL_fseek(s->fp, (off_t) offset, SEEK_SET);
    AS_UTL_safeWrite(s->fp, element, "setIndexStore", s->header.elementSize, 1);
  }
}


void
appendIndexStore(StoreStruct *s, void *element){
  int64 offset;

  assert(s->status == ActiveStore);
  s->isDirty = 1;

  offset = computeOffset(s, s->header.lastElem + 1);

  if (s->isMemoryStore){
    int64 desiredSize = offset + (s->header.lastElem - s->header.firstElem + 1) * s->header.elementSize;

    if (s->allocatedSize <= desiredSize) {
      while(s->allocatedSize <= desiredSize)
	s->allocatedSize *= 2;
      s->buffer = (char *)safe_realloc(s->buffer, s->allocatedSize);
    }
  }

  s->header.lastElem++;

  setIndexStore(s, s->header.lastElem, element);
}


StoreStruct *
createVLRecordStore(const char *path, const char *storeType, 
                    int32 expectedVLRecordSize) {
  StoreStruct *s = allocateStore();

  s->header.type = VLRECORD_STORE;

  strncpy(s->header.storeType, storeType,7);
  s->header.elementSize = expectedVLRecordSize;
  s->header.firstElem = 0;
  s->header.lastElem =  0; /* We use this to track total length */
  s->lastCommittedElem = -1;
  s->status = UnInitializedStore;
  s->header.creationTime = time(0);

  if(path){ /* File based store */
    s->isMemoryStore = 0;
    assert(strlen(path) < FILENAME_MAX);
    s->isMemoryStore = 0;
    errno = 0;
    s->fp = fopen(path,"w+");
    if(errno){
      fprintf(stderr,"createVLRecordStore()-- Failure opening Store %s for w+: %s\n", path, strerror(errno));
      assert(errno == 0);
    }
    s->buffer = NULL;
    if (WRITING_BUFFER > 0) {
      //fprintf(stderr, "createVLRecordStore()--  Allocating fp buffer of "F_SIZE_T" bytes for '%s'.\n", WRITING_BUFFER, path);
      s->buffer = (char *)safe_calloc(sizeof(char), WRITING_BUFFER);
      setbuffer(s->fp, s->buffer, WRITING_BUFFER);
    }
  }else{
    s->allocatedSize = expectedVLRecordSize;
    s->buffer = (char *)safe_calloc(sizeof(char), s->allocatedSize);
    AssertPtr(s->buffer);
    s->isMemoryStore = 1;
    s->fp = NULL;
  }

  s->status = ActiveStore;

  writeHeader(s); 

  return s;
}


void
getVLRecordStore(StoreStruct *s, int64 offset, void *buffer, VLSTRING_SIZE_T maxLength, VLSTRING_SIZE_T *actualLength){
  VLSTRING_SIZE_T length = 0;
  int64 actualOffset;

  assert(s->header.type == VLRECORD_STORE);
  if (offset > s->header.lastElem)
    fprintf(stderr," ERROR: VLRecordStoreGet at offset "F_S64"...lastElem = "F_S64"\n",
	    offset, s->header.lastElem);
  assert(offset <= s->header.lastElem);
  
  actualOffset = offset - s->header.firstElem;

  if(s->isMemoryStore){
    memcpy(&length,(s->buffer + actualOffset + sizeof(StoreStat)),
           sizeof(VLSTRING_SIZE_T));
    if(!(length + actualOffset <= s->header.lastElem &&
	 length <= maxLength)){
      fprintf(stderr,
              "getVLRecordStore(memory) FAILURE at offset "F_S64" actualOffset "F_S64"\n\t"
	      "length = " F_VLS " implied max offset = "F_S64"  "
              "lastElem = "F_S64" maxLength = " F_VLS "\n",
	      offset, actualOffset,
	      length, length + actualOffset,
              s->header.lastElem, maxLength);
      assert(0);
    }
    if(length > 0)
      memcpy(buffer, s->buffer + actualOffset +
             sizeof(StoreStat) + sizeof(VLSTRING_SIZE_T), length);

    
  }else{
    AS_UTL_fseek(s->fp, (off_t) (actualOffset + sizeof(StoreStat)), SEEK_SET);

    if (s->lastWasRead == 0)
      fflush(s->fp);
    s->lastWasRead = 1;

    if (1 != AS_UTL_safeRead(s->fp,&length,"getVLRecordStore",sizeof(VLSTRING_SIZE_T), 1)) {
      fprintf(stderr, "getVLRecordStore()-- Failed to read the length of the record.  Incomplete store?\n");
      exit(1);
    }

    if(!(length + offset + sizeof(VLSTRING_SIZE_T) <= s->header.lastElem &&
	 length <= maxLength)){
      fprintf(stderr,"* getVLRecordStore(file) FAILURE: Length = " F_VLS " offset = "F_S64" "
	      "maxLength = " F_VLS " lastElem = "F_S64"\n",
	      length,offset, maxLength, s->header.lastElem);
      assert(0);
    }

    if (length != AS_UTL_safeRead(s->fp,buffer,"getVLRecordStore",sizeof(char), length)) {
      fprintf(stderr, "getVLRecordStore()-- Failed to read all "F_VLS" bytes.  Incomplete store?\n", length);
      exit(1);
    }
  }
  *actualLength = length;
}


void
appendVLRecordStore(StoreStruct *s, void *vlr, VLSTRING_SIZE_T length){
  int64 offset;

  assert(s->status == ActiveStore);
  assert(s->header.type == VLRECORD_STORE);
  

  offset = sizeof(StoreStat) + s->header.lastElem;

  if(s->isMemoryStore){
    if(s->allocatedSize <= offset + length + sizeof(VLSTRING_SIZE_T) ){
      while(s->allocatedSize <= offset + length + sizeof(VLSTRING_SIZE_T)){
	s->allocatedSize *= 2;
      }

      s->buffer = (char *)
	safe_realloc(s->buffer, s->allocatedSize);
    }
    memcpy(s->buffer + offset, &length, sizeof(VLSTRING_SIZE_T));
    if(length > 0)
      memcpy(s->buffer + offset + sizeof(VLSTRING_SIZE_T), vlr, length);
  }else{
    AS_UTL_fseek(s->fp, (off_t) offset, SEEK_SET);

    if (s->lastWasRead)
      fflush(s->fp);
    s->lastWasRead = 0;

    AS_UTL_safeWrite(s->fp, &length, "appendVLRecordStore", sizeof(VLSTRING_SIZE_T), 1);
    AS_UTL_safeWrite(s->fp, vlr, "appendVLRecordStore", sizeof(char), length);
  }

  s->header.lastElem += length + sizeof(VLSTRING_SIZE_T);
  s->isDirty = 1;
}


StreamStruct *
openStream(StoreStruct *fs,
           void *buffer,
           int32 bufferSize) {
  return(allocateStream(fs,buffer,bufferSize));
}


void
resetStream(StreamStruct *ss, int64 startIndex, int64 endIndex){

  if(startIndex == STREAM_FROMSTART)
    ss->startIndex = ss->store->header.firstElem;
  else
    ss->startIndex = startIndex;
  if(endIndex == STREAM_UNTILEND)
    ss->endIndex = ss->store->header.lastElem;
  else
    ss->endIndex = endIndex;
}


int
nextStream(StreamStruct *ss, void *buffer){

  assert(ss->status == ActiveStore);

  if (ss->startIndex > ss->endIndex)
    return(0);

  getIndexStore(ss->store, ss->startIndex++, buffer);
  return(1);
}


void
closeStream(StreamStruct *ss){
  safe_free(ss);
}


int
nextVLRecordStream(StreamStruct *ss, void *buffer, VLSTRING_SIZE_T maxLength, VLSTRING_SIZE_T *actualLength){
  int result;

  assert(ss->status == ActiveStore);

  if (ss->startIndex >  ss->endIndex)
    return(0);

  getVLRecordStore(ss->store, ss->startIndex, buffer, maxLength, actualLength);
  ss->startIndex += (*actualLength + sizeof(VLSTRING_SIZE_T));

  return(1);
}




StoreStruct *
convertStoreToPartialMemoryStore(StoreStruct *source,
                                 int64         firstElem,
                                 int64         lastElem) {
  StoreStruct *target         = NULL;
  int64        sourceOffset    = 0;
  int64        sourceMaxOffset = 0;

  if (firstElem <= 0)
    firstElem = source->header.firstElem;

  if (lastElem <= 0)
    lastElem = source->header.lastElem;

  if (source->header.type == VLRECORD_STORE) {
    sourceOffset    = sizeof(StoreStat) + firstElem;
    sourceMaxOffset = sizeof(StoreStat) + lastElem;
  }else{
    sourceOffset    = computeOffset(source, firstElem);
    sourceMaxOffset = computeOffset(source, lastElem + 1);
  }

  target  = allocateStore();
  *target = *source;

  target->allocatedSize    = sizeof(StoreStat) + sourceMaxOffset - sourceOffset;
  target->buffer           = (char *)safe_calloc(sizeof(char), target->allocatedSize);
  target->isMemoryStore    = 1;
  target->fp               = NULL;
  target->header.firstElem = firstElem;
  target->header.lastElem  = lastElem;

  //  Copy the store into the memory store.  The first
  //  sizeof(StoreStat) bytes are for a copy of target->header
  //  (written by writeHeader()), which is never used.
  //
  AS_UTL_fseek(source->fp, (off_t)sourceOffset, SEEK_SET);

  if (source->lastWasRead == 0)
    fflush(source->fp);
  source->lastWasRead = 1;

  fprintf(stderr, "reading "F_S64" bytes\n", sourceMaxOffset - sourceOffset);

  if (sourceMaxOffset - sourceOffset != AS_UTL_safeRead(source->fp, target->buffer + sizeof(StoreStat),
                                                        "convertStoreToPartialMemoryStore",
                                                        sizeof(char), sourceMaxOffset - sourceOffset))
    assert(0);

  closeStore(source);

  return(target);
}
